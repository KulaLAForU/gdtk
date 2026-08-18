// characteristic_outflow.d
//
// Supersonic characteristic outflow based on Christian Stemmer's char_bc
// routine (Paul Harris' University of Adelaide thesis formulation).

module lmr.bc.ghost_cell_effect.characteristic_outflow;

import std.math;

import gas;
import geom;

import lmr.bc;
import lmr.flowstate;
import lmr.fluidblock;
import lmr.fluidfvcell;
import lmr.fvinterface;
import lmr.globalconfig;
import lmr.sfluidblock;


class GhostCellCharacteristicOutflow : GhostCellEffect {
public:
    this(int id, int boundary)
    {
        super(id, boundary, "CharacteristicOutflow");
    }

    override string toString() const
    {
        return "CharacteristicOutflow()";
    }

    @nogc
    override void apply_for_interface_unstructured_grid(double t, int gtl, int ftl,
                                                         FVInterface f)
    {
        throw new Error("CharacteristicOutflow is only available for structured grids.");
    }

    @nogc
    override void apply_unstructured_grid(double t, int gtl, int ftl)
    {
        throw new Error("CharacteristicOutflow is only available for structured grids.");
    }

    @nogc
    override void apply_for_interface_structured_grid(double t, int gtl, int ftl,
                                                       FVInterface f)
    {
        applyFace(f.i_bndry, gtl);
    }

    @nogc
    override void apply_structured_grid(double t, int gtl, int ftl)
    {
        auto sblk = cast(SFluidBlock) blk;
        assert(sblk !is null);
        if (sblk.myConfig.dimensions != 2 ||
            (which_boundary != Face.north && which_boundary != Face.south)) {
            throw new Error("CharacteristicOutflow currently requires a north or south boundary of a 2-D structured grid.");
        }
        foreach (i; 0 .. sblk.bc[which_boundary].faces.length) applyFace(i, gtl);
    }

private:
    @nogc
    FluidFVCell interiorCell(FVInterface f, size_t layer)
    {
        auto bc = blk.bc[which_boundary];
        return bc.outsigns[f.i_bndry] == 1 ? f.left_cells[layer] : f.right_cells[layer];
    }

    @nogc
    FluidFVCell ghostCell(FVInterface f, size_t layer)
    {
        auto bc = blk.bc[which_boundary];
        return bc.outsigns[f.i_bndry] == 1 ? f.right_cells[layer] : f.left_cells[layer];
    }

    @nogc
    void applyFace(size_t ifaceIndex, int gtl)
    {
        auto sblk = cast(SFluidBlock) blk;
        assert(sblk !is null);
        auto bc = sblk.bc[which_boundary];
        auto face = bc.faces[ifaceIndex];
        auto nearest = interiorCell(face, 0);

        const haveInteriorStencil = bc.outsigns[ifaceIndex] == 1 ?
            face.left_cells.length >= 2 : face.right_cells.length >= 2;
        if (sblk.myConfig.dimensions != 2 ||
            (which_boundary != Face.north && which_boundary != Face.south) ||
            !haveInteriorStencil) {
            copyIntoGhosts(nearest.fs, face);
            return;
        }

        const speed = sqrt(nearest.fs.vel.x*nearest.fs.vel.x +
                           nearest.fs.vel.y*nearest.fs.vel.y +
                           nearest.fs.vel.z*nearest.fs.vel.z);
        const mach = speed / nearest.fs.gas.a;
        if (!isFinite(mach) || mach <= 1.0) {
            copyIntoGhosts(nearest.fs, face);
            return;
        }

        // The outgoing characteristic is theta+mu at a north boundary and
        // theta-mu at a south boundary.  Trace it from the boundary into the
        // first two cell layers, using the actual (possibly nonuniform) grid.
        const outsign = bc.outsigns[ifaceIndex];
        const theta = atan2(nearest.fs.vel.y, nearest.fs.vel.x);
        const mu = asin(1.0/mach);
        const tanCharacteristic = tan(theta + outsign*mu);
        if (!isFinite(tanCharacteristic) || abs(tanCharacteristic) < 1.0e-12) {
            copyIntoGhosts(nearest.fs, face);
            return;
        }

        auto layer1Here = interiorCell(face, 0);
        auto layer2Here = interiorCell(face, 1);
        const yBoundary = face.pos.y;
        const x1 = face.pos.x + (layer1Here.pos[gtl].y-yBoundary)/tanCharacteristic;
        const x2 = face.pos.x + (layer2Here.pos[gtl].y-yBoundary)/tanCharacteristic;

        size_t i10, i11, i20, i21;
        double w1, w2;
        if (!bracketLayer(x1, 0, gtl, i10, i11, w1) ||
            !bracketLayer(x2, 1, gtl, i20, i21, w2)) {
            copyIntoGhosts(nearest.fs, face);
            return;
        }

        auto a1 = interiorCell(bc.faces[i10], 0).fs;
        auto b1 = interiorCell(bc.faces[i11], 0).fs;
        auto a2 = interiorCell(bc.faces[i20], 1).fs;
        auto b2 = interiorCell(bc.faces[i21], 1).fs;
        auto dest = ghostCell(face, 0).fs;
        extrapolateState(a1, b1, w1, a2, b2, w2, dest);
        foreach (n; 1 .. sblk.n_ghost_cell_layers) ghostCell(face, n).fs.copy_values_from(dest);
    }

    @nogc
    bool bracketLayer(double x, size_t layer, int gtl,
                      ref size_t i0, ref size_t i1, ref double weight)
    {
        auto faces = blk.bc[which_boundary].faces;
        if (faces.length == 0) return false;
        const firstX = interiorCell(faces[0], layer).pos[gtl].x;
        const lastX = interiorCell(faces[$-1], layer).pos[gtl].x;
        const increasing = lastX >= firstX;
        if ((increasing && (x < firstX || x > lastX)) ||
            (!increasing && (x > firstX || x < lastX))) return false;
        foreach (i; 0 .. faces.length-1) {
            const xa = interiorCell(faces[i], layer).pos[gtl].x;
            const xb = interiorCell(faces[i+1], layer).pos[gtl].x;
            if ((x-xa)*(x-xb) <= 0.0) {
                i0 = i; i1 = i+1;
                weight = abs(xb-xa) > 1.0e-14 ? (x-xa)/(xb-xa) : 0.0;
                return true;
            }
        }
        // A one-cell-wide boundary can only supply an exactly coincident foot.
        if (faces.length == 1 && abs(x-firstX) < 1.0e-14) {
            i0 = i1 = 0; weight = 0.0; return true;
        }
        return false;
    }

    @nogc
    void extrapolateState(FlowState* a1, FlowState* b1, double w1,
                          FlowState* a2, FlowState* b2, double w2,
                          FlowState* dest)
    {
        auto gmodel = blk.myConfig.gmodel;
        // Interpolate on each characteristic helper line, then apply the
        // original (4*q1-q2)/3 extrapolation on primitive quantities.
        dest.gas.rho = (4.0*((1.0-w1)*a1.gas.rho+w1*b1.gas.rho) -
                            ((1.0-w2)*a2.gas.rho+w2*b2.gas.rho))/3.0;
        dest.gas.u = (4.0*((1.0-w1)*a1.gas.u+w1*b1.gas.u) -
                          ((1.0-w2)*a2.gas.u+w2*b2.gas.u))/3.0;
        version(multi_T_gas) {
            foreach (m; 0 .. blk.myConfig.n_modes) {
                dest.gas.u_modes[m] = (4.0*((1.0-w1)*a1.gas.u_modes[m]+w1*b1.gas.u_modes[m]) -
                                       ((1.0-w2)*a2.gas.u_modes[m]+w2*b2.gas.u_modes[m]))/3.0;
            }
        }
        version(multi_species_gas) {
            foreach (isp; 0 .. blk.myConfig.n_species) {
                dest.gas.massf[isp] = (4.0*((1.0-w1)*a1.gas.massf[isp]+w1*b1.gas.massf[isp]) -
                                       ((1.0-w2)*a2.gas.massf[isp]+w2*b2.gas.massf[isp]))/3.0;
            }
            scale_mass_fractions(dest.gas.massf);
        }
        gmodel.update_thermo_from_rhou(dest.gas);
        dest.vel.x = (4.0*((1.0-w1)*a1.vel.x+w1*b1.vel.x) - ((1.0-w2)*a2.vel.x+w2*b2.vel.x))/3.0;
        dest.vel.y = (4.0*((1.0-w1)*a1.vel.y+w1*b1.vel.y) - ((1.0-w2)*a2.vel.y+w2*b2.vel.y))/3.0;
        dest.vel.z = (4.0*((1.0-w1)*a1.vel.z+w1*b1.vel.z) - ((1.0-w2)*a2.vel.z+w2*b2.vel.z))/3.0;
        // Auxiliary quantities are less central to the characteristic model;
        // retain a bounded, consistent value from the nearest interior cell.
        dest.mu_t = a1.mu_t;
        dest.k_t = a1.k_t;
        version(turbulence) dest.turb = a1.turb;
        version(MHD) { dest.B = a1.B; dest.psi = a1.psi; dest.divB = a1.divB; }
    }

    @nogc
    void copyIntoGhosts(FlowState* src, FVInterface face)
    {
        auto sblk = cast(SFluidBlock) blk;
        foreach (n; 0 .. sblk.n_ghost_cell_layers) ghostCell(face, n).fs.copy_values_from(src);
    }
}
