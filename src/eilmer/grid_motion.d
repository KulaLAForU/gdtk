// Authors: RG, PJ, KD & IJ
// Date: 2015-11-20

module grid_motion;

import std.string;
import std.conv;

import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;
import nm.luabbla;
import lua_helper;
import fvcore;
import fvvertex;
import fvinterface;
import globalconfig;
import globaldata;
import geom;
import geom.luawrap;
import fluidblock;
import sfluidblock;
import std.stdio;


@nogc
int set_gcl_interface_properties(SFluidBlock blk, size_t gtl, double dt) {
    size_t i, j, k;
    FVInterface IFace;
    Vector3 pos1, pos2, temp;
    Vector3 averaged_ivel, vol;
    if (blk.myConfig.dimensions == 2) {
    FVVertex vtx1, vtx2;
    k = blk.kmin;
    // loop over i-interfaces and compute interface velocity wif'.
    for (j = blk.jmin; j <= blk.jmax; ++j) {
        for (i = blk.imin; i <= blk.imax+1; ++i) {//  i <= blk.imax+1
            vtx1 = blk.get_vtx(i,j,k);
            vtx2 = blk.get_vtx(i,j+1,k);
            IFace = blk.get_ifi(i,j,k);         
            pos1 = vtx1.pos[gtl];
            pos1 -= vtx2.pos[0];
            pos2 = vtx2.pos[gtl];
            pos2 -= vtx1.pos[0];
            averaged_ivel = vtx1.vel[0];
            averaged_ivel += vtx2.vel[0];
            averaged_ivel.scale(0.5);
            // Use effective edge velocity
            // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
            // Full Potential and Euler solutions for transonic unsteady flow
            // Aeronautical Journal November 1994 Eqn 25
            cross(vol, pos1, pos2);
            if (blk.myConfig.axisymmetric == false) {
                // vol = 0.5*cross(pos1, pos2);
                vol.scale(0.5);
            } else {
                // vol=0.5*cross(pos1, pos2)*((vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y)/4.0);
                vol.scale(0.125*(vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y));
            }
            temp = vol; temp /= dt*IFace.area[0];
            // temp is the interface velocity (W_if) from the GCL
            // interface area determined at gtl 0 since GCL formulation
            // recommends using initial interfacial area in calculation.
            IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            IFace.gvel.set(temp.z, averaged_ivel.y, averaged_ivel.z);
            averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
            IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);                           
        }
    }
    // loop over j-interfaces and compute interface velocity wif'.
    for (j = blk.jmin; j <= blk.jmax+1; ++j) { //j <= blk.jmax+
        for (i = blk.imin; i <= blk.imax; ++i) {
            vtx1 = blk.get_vtx(i,j,k);
            vtx2 = blk.get_vtx(i+1,j,k);
            IFace = blk.get_ifj(i,j,k);
            pos1 = vtx2.pos[gtl]; pos1 -= vtx1.pos[0];
            pos2 = vtx1.pos[gtl]; pos2 -= vtx2.pos[0];
            averaged_ivel = vtx1.vel[0]; averaged_ivel += vtx2.vel[0]; averaged_ivel.scale(0.5);
            cross(vol, pos1, pos2);
            if (blk.myConfig.axisymmetric == false) {
                // vol=0.5*cross( pos1, pos2 );
                vol.scale(0.5);
            } else {
                // vol=0.5*cross(pos1, pos2)*((vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y)/4.0);
                vol.scale(0.125*(vtx1.pos[gtl].y+vtx1.pos[0].y+vtx2.pos[gtl].y+vtx2.pos[0].y));
            }
            IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
            averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);  
            if (blk.myConfig.axisymmetric && j == blk.jmin && IFace.area[0] == 0.0) {
                // For axi-symmetric cases the cells along the axis of symmetry have 0 interface area,
                // this is a problem for determining Wif, so we have to catch the NaN from dividing by 0.
                // We choose to set the y and z directions to 0, but take an averaged value for the
                // x-direction so as to not force the grid to be stationary, defeating the moving grid's purpose.
                IFace.gvel.set(averaged_ivel.x, to!number(0.0), to!number(0.0));
            } else {
                temp = vol; temp /= dt*IFace.area[0];
                IFace.gvel.set(temp.z, averaged_ivel.y, averaged_ivel.z);
            }
            averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
            IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
        }
    } 
    }

    // Do 3-D cases, where faces move and create a new hexahedron
    if (blk.myConfig.dimensions == 3) {
        FVVertex vtx0, vtx1, vtx2, vtx3;
        Vector3 p0, p1, p2, p3, p4, p5, p6, p7;
        Vector3 centroid_hex, sub_centroid;
        number volume, sub_volume, temp2;
        // loop over i-interfaces and compute interface velocity wif'.
        for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (j = blk.jmin; j <= blk.jmax; ++j) {
                for (i = blk.imin; i <= blk.imax+1; ++i) {//  i <= blk.imax+1
                    // Calculate volume generated by sweeping face 0123 from pos[0] to pos[gtl]
                    vtx0 = blk.get_vtx(i,j  ,k  );
                    vtx1 = blk.get_vtx(i,j+1,k  );
                    vtx2 = blk.get_vtx(i,j+1,k+1);
                    vtx3 = blk.get_vtx(i,j  ,k+1  );
                    p0 = vtx0.pos[0]; p1 = vtx1.pos[0];
                    p2 = vtx2.pos[0]; p3 = vtx3.pos[0];
                    p4 = vtx0.pos[gtl]; p5 = vtx1.pos[gtl];
                    p6 = vtx2.pos[gtl]; p7 = vtx3.pos[gtl];
                    // use 6x pyramid approach as used to calculate internal volume of hex cells
                    centroid_hex.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    volume = 0.0; 
                    pyramid_properties(p6, p7, p3, p2, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p5, p6, p2, p1, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p4, p5, p1, p0, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p4, p0, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p6, p5, p4, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p0, p1, p2, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    //
                    IFace = blk.get_ifi(i,j,k);         
                    averaged_ivel = vtx0.vel[0];
                    averaged_ivel += vtx1.vel[0];
                    averaged_ivel += vtx2.vel[0];
                    averaged_ivel += vtx3.vel[0];
                    averaged_ivel.scale(0.25);
                    // Use effective face velocity, analoguous to edge velocity concept
                    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
                    // Full Potential and Euler solutions for transonic unsteady flow
                    // Aeronautical Journal November 1994 Eqn 25
                    temp2 = volume; temp /= dt*IFace.area[0];
                    // temp2 is the interface velocity (W_if) from the GCL
                    // interface area determined at gtl 0 since GCL formulation
                    // recommends using initial interfacial area in calculation.
                    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    IFace.gvel.set(temp2, averaged_ivel.y, averaged_ivel.z);
                    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
                    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);  
                }
            }
        }
        // loop over j-interfaces and compute interface velocity wif'.
        for (k = blk.kmin; k <= blk.kmax; ++k) {
            for (j = blk.jmin; j <= blk.jmax+1; ++j) { //j <= blk.jmax+
                for (i = blk.imin; i <= blk.imax; ++i) {
                    // Calculate volume generated by sweeping face 0123 from pos[0] to pos[gtl]
                    vtx0 = blk.get_vtx(i  ,j,k  );
                    vtx1 = blk.get_vtx(i  ,j,k+1);
                    vtx2 = blk.get_vtx(i+1,j,k+1);
                    vtx3 = blk.get_vtx(i+1,j,k  );
                    p0 = vtx0.pos[0]; p1 = vtx1.pos[0];
                    p2 = vtx2.pos[0]; p3 = vtx3.pos[0];
                    p4 = vtx0.pos[gtl]; p5 = vtx1.pos[gtl];
                    p6 = vtx2.pos[gtl]; p7 = vtx3.pos[gtl];
                    // use 6x pyramid approach as used to calculate internal volume of hex cells
                    centroid_hex.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    volume = 0.0; 
                    pyramid_properties(p6, p7, p3, p2, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p5, p6, p2, p1, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p4, p5, p1, p0, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p4, p0, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p6, p5, p4, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p0, p1, p2, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    // 
                    IFace = blk.get_ifj(i,j,k);         
                    averaged_ivel = vtx0.vel[0];
                    averaged_ivel += vtx1.vel[0];
                    averaged_ivel += vtx2.vel[0];
                    averaged_ivel += vtx3.vel[0];
                    averaged_ivel.scale(0.25);
                    // Use effective face velocity, analoguous to edge velocity concept
                    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
                    // Full Potential and Euler solutions for transonic unsteady flow
                    // Aeronautical Journal November 1994 Eqn 25
                    temp2 = volume; temp /= dt*IFace.area[0];
                    // temp2 is the interface velocity (W_if) from the GCL
                    // interface area determined at gtl 0 since GCL formulation
                    // recommends using initial interfacial area in calculation.
                    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    IFace.gvel.set(temp2, averaged_ivel.y, averaged_ivel.z);
                    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
                    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);  
                }
            } 
        }
        // loop over k-interfaces and compute interface velocity wif'.
        for (k = blk.kmin; k <= blk.kmax+1; ++k) {
            for (j = blk.jmin; j <= blk.jmax; ++j) { //j <= blk.jmax+
                for (i = blk.imin; i <= blk.imax; ++i) {
                    // Calculate volume generated by sweeping face 0123 from pos[0] to pos[gtl]
                    vtx0 = blk.get_vtx(i  ,j  ,k);
                    vtx1 = blk.get_vtx(i+1,j  ,k);
                    vtx2 = blk.get_vtx(i+1,j+1,k);
                    vtx3 = blk.get_vtx(i  ,j+1,k);
                    p0 = vtx0.pos[0]; p1 = vtx1.pos[0];
                    p2 = vtx2.pos[0]; p3 = vtx3.pos[0];
                    p4 = vtx0.pos[gtl]; p5 = vtx1.pos[gtl];
                    p6 = vtx2.pos[gtl]; p7 = vtx3.pos[gtl];
                    // use 6x pyramid approach as used to calculate internal volume of hex cells
                    centroid_hex.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
                    volume = 0.0; 
                    pyramid_properties(p6, p7, p3, p2, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p5, p6, p2, p1, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p4, p5, p1, p0, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p4, p0, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p7, p6, p5, p4, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    pyramid_properties(p0, p1, p2, p3, centroid_hex, sub_centroid, sub_volume);
                    volume += sub_volume; 
                    //
                    IFace = blk.get_ifk(i,j,k);         
                    averaged_ivel = vtx0.vel[0];
                    averaged_ivel += vtx1.vel[0];
                    averaged_ivel += vtx2.vel[0];
                    averaged_ivel += vtx3.vel[0];
                    averaged_ivel.scale(0.25);
                    // Use effective face velocity, analoguous to edge velocity concept
                    // Reference: D. Ambrosi, L. Gasparini and L. Vigenano
                    // Full Potential and Euler solutions for transonic unsteady flow
                    // Aeronautical Journal November 1994 Eqn 25
                    temp2 = volume; temp2 /= dt*IFace.area[0];
                    // temp2 is the interface velocity (W_if) from the GCL
                    // interface area determined at gtl 0 since GCL formulation
                    // recommends using initial interfacial area in calculation.
                    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    averaged_ivel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
                    IFace.gvel.set(temp2, averaged_ivel.y, averaged_ivel.z);
                    averaged_ivel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);           
                    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);  
                }
            } 
        }
    } 
    return 0;
} // end set_gcl_interface_properties()

@nogc
void predict_vertex_positions(SFluidBlock blk, double dt, int gtl)
{
    // [TODO] PJ 2018-09-10: Could generalize by just looping over all vertices in FluidBlock list.
    size_t krangemax = ( blk.myConfig.dimensions == 2 ) ? blk.kmax : blk.kmax+1;
    for (size_t k = blk.kmin; k <= krangemax; ++k) {
        for (size_t j = blk.jmin; j <= blk.jmax+1; ++j) {
            for (size_t i = blk.imin; i <= blk.imax+1; ++i) {
                FVVertex vtx = blk.get_vtx(i,j,k);
                if (gtl == 0) {
                    // predictor/sole step; update grid
                    vtx.pos[1] = vtx.pos[0] + dt * vtx.vel[0];
                } else {
                    // corrector step; keep grid fixed
                    vtx.pos[2] = vtx.pos[1];
                }
            }
        }
    }
    return;
} // end predict_vertex_positions()

// ------------------ Lua interface functions below this point. ----------------------

void setGridMotionHelperFunctions(lua_State *L)
{
    lua_pushcfunction(L, &luafn_getVtxPosition);
    lua_setglobal(L, "getVtxPosition");
    lua_pushcfunction(L, &luafn_getVtxPositionXYZ);
    lua_setglobal(L, "getVtxPositionXYZ");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForDomain);
    lua_setglobal(L, "setVtxVelocitiesForDomain");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForDomainXYZ);
    lua_setglobal(L, "setVtxVelocitiesForDomainXYZ");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForBlock);
    lua_setglobal(L, "setVtxVelocitiesForBlock");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForBlockXYZ);
    lua_setglobal(L, "setVtxVelocitiesForBlockXYZ");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesForRotatingBlock);
    lua_setglobal(L, "setVtxVelocitiesForRotatingBlock");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesByCorners);
    lua_setglobal(L, "setVtxVelocitiesByCorners");
    lua_pushcfunction(L, &luafn_setVtxVelocitiesByCornersReg);
    lua_setglobal(L, "setVtxVelocitiesByCornersReg");
    lua_pushcfunction(L, &luafn_setVtxVelocity);
    lua_setglobal(L, "setVtxVelocity");
    lua_pushcfunction(L, &luafn_setVtxVelocityXYZ);
    lua_setglobal(L, "setVtxVelocityXYZ");
}

extern(C) int luafn_getVtxPosition(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);

    // Grab the appropriate vtx
    auto vtx = globalFluidBlocks[blkId].get_vtx(i, j, k);
    
    // Return the interesting bits as a table with entries x, y, z.
    lua_newtable(L);
    lua_pushnumber(L, vtx.pos[0].x.re); lua_setfield(L, -2, "x");
    lua_pushnumber(L, vtx.pos[0].y.re); lua_setfield(L, -2, "y");
    lua_pushnumber(L, vtx.pos[0].z.re); lua_setfield(L, -2, "z");
    return 1;
}

extern(C) int luafn_getVtxPositionXYZ(lua_State *L)
{
    // Get arguments from lua_stack
    auto blkId = lua_tointeger(L, 1);
    auto i = lua_tointeger(L, 2);
    auto j = lua_tointeger(L, 3);
    auto k = lua_tointeger(L, 4);

    // Grab the appropriate vtx
    auto vtx = globalFluidBlocks[blkId].get_vtx(i, j, k);
    
    // Return the components x, y, z on the stack.
    lua_pushnumber(L, vtx.pos[0].x.re);
    lua_pushnumber(L, vtx.pos[0].y.re);
    lua_pushnumber(L, vtx.pos[0].z.re);
    return 3;
}

extern(C) int luafn_setVtxVelocitiesForDomain(lua_State* L)
{
    // Expect a single argument: a Vector3 object
    auto vel = checkVector3(L, 1);

    foreach ( blk; localFluidBlocks ) {
        foreach ( vtx; blk.vertices ) {
            /* We assume that we'll only update grid positions
               at the start of the increment. This should work
               well except in the most critical cases of time
               accuracy.
            */
            vtx.vel[0].set(vel);
        }
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

extern(C) int luafn_setVtxVelocitiesForDomainXYZ(lua_State* L)
{
    // Expect three velocity components.
    double velx = lua_tonumber(L, 1);
    double vely = lua_tonumber(L, 2);
    double velz = lua_tonumber(L, 3);

    foreach ( blk; localFluidBlocks ) {
        foreach ( vtx; blk.vertices ) {
            /* We assume that we'll only update grid positions
               at the start of the increment. This should work
               well except in the most critical cases of time
               accuracy.
            */
            vtx.vel[0].set(velx, vely, velz);
        }
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

extern(C) int luafn_setVtxVelocitiesForBlock(lua_State* L)
{
    // Expect two arguments: 1. a block id
    //                       2. a Vector3 object
    auto blkId = lua_tointeger(L, 1);
    auto vel = checkVector3(L, 2);

    foreach ( vtx; globalFluidBlocks[blkId].vertices ) {
        /* We assume that we'll only update grid positions
           at the start of the increment. This should work
           well except in the most critical cases of time
           accuracy.
        */
        vtx.vel[0].set(vel);
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

extern(C) int luafn_setVtxVelocitiesForBlockXYZ(lua_State* L)
{
    // Expect two arguments: 1. a block id
    //                       2. a Vector3 object
    auto blkId = lua_tointeger(L, 1);
    double velx = lua_tonumber(L, 2);
    double vely = lua_tonumber(L, 3);
    double velz = lua_tonumber(L, 4);

    foreach ( vtx; globalFluidBlocks[blkId].vertices ) {
        /* We assume that we'll only update grid positions
           at the start of the increment. This should work
           well except in the most critical cases of time
           accuracy.
        */
        vtx.vel[0].set(velx, vely, velz);
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

/**
 * Sets the velocity of an entire block, for the case
 * that the block is rotating about an axis with direction 
 * (0 0 1) located at a point defined by vector (x y 0).
 *
 * setVtxVelocitiesRotatingBlock(blkId, omega, vector3)
 *      Sets rotational speed omega (rad/s) for rotation about 
 *      (0 0 1) axis defined by Vector3.
 *
 * setVtxVelocitiesRotatingBlock(blkId, omega)
 *      Sets rotational speed omega (rad/s) for rotation about 
 *      Z-axis.
 */
extern(C) int luafn_setVtxVelocitiesForRotatingBlock(lua_State* L)
{
    // Expect two/three arguments: 1. a block id
    //                             2. a float object
    //                             3. a vector (optional) 
    int narg = lua_gettop(L);
    auto blkId = lua_tointeger(L, 1);
    double omega = lua_tonumber(L, 2);
    double velx, vely;

    if ( narg == 2 ) {
        // assume rotation about Z-axis 
        foreach ( vtx; globalFluidBlocks[blkId].vertices ) {
            velx = - omega * vtx.pos[0].y.re;
            vely =   omega * vtx.pos[0].x.re;
            vtx.vel[0].set(velx, vely, 0.0);
        }
    }
    else if ( narg == 3 ) {
        auto axis = checkVector3(L, 3);
        foreach ( vtx; globalFluidBlocks[blkId].vertices ) {
            velx = - omega * (vtx.pos[0].y.re - axis.y.re);
            vely =   omega * (vtx.pos[0].x.re - axis.x.re);
            vtx.vel[0].set(velx, vely, 0.0);
        }
    }
    else {
        string errMsg = "ERROR: Too few arguments passed to luafn: setVtxVelocitiesRotatingBlock()\n";
        luaL_error(L, errMsg.toStringz);
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

/**
 * Sets the velocity of vertices in a block based on
 * specified corner velocities. The velocity of any cell 
 * is estimated using decomposition of the quad into 4x 
 * triangles and then calculates velocities using barycentric 
 * interpolation within the respective triangles. 
 * Works for clustered grids.
 *
 *   p3-----p2
 *   |      |
 *   |      |
 *   p0-----p1
 *
 * This function can be called for 2-D structured and
 * extruded 3-D grids only. The following calls are allowed:
 *
 * setVtxVelocitiesByCorners(blkId, p0vel, p1vel, p2vel, p3vel)
 *   Sets the velocity vectors for block with corner velocities 
 *   specified by four corner velocities. This works for both 
 *   2-D and 3- meshes. In 3-D a uniform velocity is applied 
 *   in k direction.
 */

extern(C) int luafn_setVtxVelocitiesByCorners(lua_State* L)
{
    // Expect five/nine arguments: 1. a block id
    //                             2-5. corner velocities for 2-D motion
    //                             6-9. corner velocities for full 3-D motion 
    //                                  (optional)
    int narg = lua_gettop(L);
    auto blkId = lua_tointeger(L, 1);
    auto blk = cast(SFluidBlock) globalFluidBlocks[blkId];
    // get corner velocities
    auto p00vel = checkVector3(L, 2);
    auto p10vel = checkVector3(L, 3);
    auto p11vel = checkVector3(L, 4);
    auto p01vel = checkVector3(L, 5);  
    // get coordinates for corner points
    size_t  k = blk.kmin;
    Vector3 p00 = blk.get_vtx(blk.imin,blk.jmin,k).pos[0];
    Vector3 p10 = blk.get_vtx(blk.imax+1,blk.jmin,k).pos[0];
    Vector3 p11 = blk.get_vtx(blk.imax+1,blk.jmax+1,k).pos[0];
    Vector3 p01 = blk.get_vtx(blk.imin,blk.jmax+1,k).pos[0];

    size_t i, j;
    Vector3 centroidVel = 0.25 * ( *p00vel + *p10vel + *p01vel + *p11vel);
    Vector3 centroid, n, pos, Coords;
    number area;
    // get centroid location
    quad_properties(p00, p10, p11, p01, centroid,  n, n, n, area);

    @nogc
    void setAsWeightedSum(ref Vector3 result,
                          double w0, Vector3* v0,
                          double w1, Vector3* v1,
                          double w2, Vector3* v2)
    {
        result.set(v0); result.scale(w0);
        result.add(v1, w1);
        result.add(v2, w2);
    }
    
    if ( narg == 5 ) {
        if (blk.myConfig.dimensions == 2) {
            // deal with 2-D meshes
            k = blk.kmin;
            for (j = blk.jmin; j <= blk.jmax+1; ++j) {
                for (i = blk.imin; i <= blk.imax+1; ++i) {
                    // get position of current point                    
                    pos = blk.get_vtx(i,j,k).pos[0];
                    //writeln("pos",pos);

                    // find baricentric ccordinates, always keep p0 at centroid
                    // sequentially try the different triangles.
                    // try south triangle
                    P_barycentricCoords(pos, centroid, p00, p10, Coords);
                    if ((Coords.x <= 0 ) || (Coords.y >= 0 && Coords.z >= 0)) {
                        setAsWeightedSum(blk.get_vtx(i,j,k).vel[0], Coords.x.re, &centroidVel, 
                                         Coords.y.re, p00vel, Coords.z.re, p10vel);
                        //writeln("Vel-S",  Coords.x * centroidVel 
                        //    + Coords.y * *p00vel + Coords.z * *p10vel,i,j);
                        continue;
                    }
                    // try east triangle
                    P_barycentricCoords(pos, centroid, p10, p11, Coords);
                    if ((Coords.x <= 0 ) || (Coords.y >= 0 && Coords.z >= 0)) {
                        setAsWeightedSum(blk.get_vtx(i,j,k).vel[0], Coords.x.re, &centroidVel, 
                                         Coords.y.re, p10vel, Coords.z.re, p11vel);
                        //writeln("Vel-E",  Coords.x * centroidVel 
                        //    + Coords.y * *p10vel + Coords.z * *p11vel,i,j);
                        continue;
                    }
                    // try north triangle
                    P_barycentricCoords(pos, centroid, p11, p01, Coords);
                    if ((Coords.x <= 0 ) || (Coords.y >= 0 && Coords.z >= 0)) {
                        setAsWeightedSum(blk.get_vtx(i,j,k).vel[0], Coords.x.re, &centroidVel, 
                                         Coords.y.re, p11vel, Coords.z.re, p01vel);
                        //writeln("Vel-N",  Coords.x * centroidVel 
                        //    + Coords.y * *p11vel + Coords.z * *p01vel,i,j);
                        continue;
                    }
                    // try west triangle
                    P_barycentricCoords(pos, centroid, p01, p00, Coords);
                    if ((Coords.x <= 0 ) || (Coords.y >= 0 && Coords.z >= 0)) {
                        setAsWeightedSum(blk.get_vtx(i,j,k).vel[0], Coords.x.re, &centroidVel, 
                                         Coords.y.re, p01vel, Coords.z.re, p00vel);
                        //writeln("Vel-W",  Coords.x * centroidVel 
                        //    + Coords.y * *p01vel + Coords.z * *p00vel,i,j);
                        continue;
                    }
                    // One of the 4 continue statements should have acted by now. 
                    writeln("Pos", pos, p00, p10, p11, p01);
                    writeln("Cell-indices",i,j,k);
                    string errMsg = "ERROR: Barycentric Calculation failed 
                                     in luafn: setVtxVelocitiesByCorners()\n";
                    luaL_error(L, errMsg.toStringz);
                }
            }
        }
        else { // deal with 3-D meshesv (assume constant properties wrt k index)
            writeln("Aaah How did I get here?");

        }
    }
    else {
        string errMsg = "ERROR: Wrong number of arguments passed to luafn: setVtxVelocitiesByCorners()\n";
        luaL_error(L, errMsg.toStringz);
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}



/**
 * Sets the velocity of vertices in a block based on
 * specified corner velocities. The velocity of any cell 
 * is estimated using interpolation based on cell indices. 
 * This should only be used for blocks with regular spacing.
 * On clustered grids deformation of the mesh occurs.
 *
 *   p3-----p2
 *   |      |
 *   |      |
 *   p0-----p1
 *
 * This function can be called for structured 
 * grids only. The following calls are allowed:
 *
 * setVtxVelocitiesByCornersReg(blkId, p0vel, p1vel, p2vel, p3vel)
 *   Sets the velocity vectors for block with corner velocities 
 *   specified by four corner velocities. This works for both 
 *   2-D and 3- meshes. In 3-D a uniform velocity is applied 
 *   in k direction.
 * 
 * setVtxVelocitiesByCornersReg(blkId, p0vel, p1vel, p2vel, p3vel, 
 *                                  p4vel, p5vel, p6vel, p7vel)
 *   As above but suitable for 3-D meshes with eight specified 
 *   velocities.
 */
extern(C) int luafn_setVtxVelocitiesByCornersReg(lua_State* L)
{
    // Expect five/nine arguments: 1. a block id
    //                             2-5. corner velocities for 2-D motion
    //                             6-9. corner velocities for full 3-D motion 
    //                                  (optional)
    int narg = lua_gettop(L);
    double u, v;
    auto blkId = lua_tointeger(L, 1);
    size_t i, j, k;
    Vector3 velw, vele, veln, vels, vel;
    auto blk = cast(SFluidBlock) globalFluidBlocks[blkId];
    // get corner velocities
    auto p00vel = checkVector3(L, 2);
    auto p10vel = checkVector3(L, 3);
    auto p11vel = checkVector3(L, 4);
    auto p01vel = checkVector3(L, 5);  

    @nogc
    void setAsWeightedSum(ref Vector3 result,
                          double w0, Vector3* v0,
                          double w1, Vector3* v1)
    {
        result.set(v0); result.scale(w0);
        result.add(v1, w1);
    }

    if ( narg == 5 ) {
        if (blk.myConfig.dimensions == 2) {
            // deal with 2-D meshes
            k = blk.kmin;
            for (j = blk.jmin; j <= blk.jmax+1; ++j) {
                // find velocity along west and east edge
                v = to!double(j-blk.jmin) / to!double(blk.jmax+1-blk.jmin); 
                setAsWeightedSum(velw, v, p01vel, 1-v, p00vel);
                setAsWeightedSum(vele, v, p11vel, 1-v, p10vel);

                for (i = blk.imin; i <= blk.imax+1; ++i) {
                    //// interpolate in i direction
                    u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin); 
                    setAsWeightedSum(vel, u, &vele, 1-u, &velw);
                    //// set velocity
                    blk.get_vtx(i,j,k).vel[0].set(vel);

                    // transfinite interpolation is yielding same result, but omitted as more expensive.
                    //u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin);
                    //vels = u * *p10vel + (1 - u) * *p00vel;
                    //veln = u * *p11vel + (1 - u) * *p01vel;
                    //// do transfinite interpolation
                    //vel = (1-v)*vels + v*veln + (1-u)*velw + u*vele
                    //    - ( (1-u)*(1-v)* *p00vel + u*v* *p11vel + u*(1-v)* *p10vel + (1-u)*v* *p01vel );
                    //// set velocity
                    //blk.get_vtx(i,j,k).vel[0] = vel;
                }
            }
        }
        else { // deal with 3-D meshesv (assume constant properties wrt k index)
            for (j = blk.jmin; j <= blk.jmax+1; ++j) {
                // find velocity along west and east edge
                v = to!double(j-blk.jmin) / to!double(blk.jmax+1-blk.jmin); 
                setAsWeightedSum(velw, v, p01vel, 1-v, p00vel);
                setAsWeightedSum(vele, v, p11vel, 1-v, p10vel);


                for (i = blk.imin; i <= blk.imax+1; ++i) {
                    // interpolate in i direction
                    u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin); 
                    setAsWeightedSum(vel, u, &vele, 1-u, &velw);

                    // set velocity for all k
                    for (k = blk.kmin; k <= blk.kmax+1; ++i) {
                        blk.get_vtx(i,j,k).vel[0].set(vel);
                    }
                }
            }
        }
    }
    else if ( narg == 9 ) {
        // assume all blocks are 3-D
        writeln("setVtxVelocitiesByCorners not verified for 3-D. 
                Proceed at own peril. See /src/eilmer/grid_motion.d");
        double w;
        Vector3 velwt, velet, velt;
        // get corner velocities
        auto p001vel = checkVector3(L, 6);
        auto p101vel = checkVector3(L, 7);
        auto p111vel = checkVector3(L, 8);
        auto p011vel = checkVector3(L, 9);  
        for (j = blk.jmin; j <= blk.jmax+1; ++j) {
            // find velocity along west and east edge (four times)
            v = to!double(j-blk.jmin) / to!double(blk.jmax+1-blk.jmin); 
            setAsWeightedSum(velw, v, p01vel, 1-v, p00vel);
            setAsWeightedSum(vele, v, p11vel, 1-v, p10vel);
            setAsWeightedSum(velwt, v, p011vel, 1-v, p001vel);
            setAsWeightedSum(velet, v, p111vel, 1-v, p101vel);

            for (i = blk.imin; i <= blk.imax+1; ++i) {
                // interpolate in i direction (twice)
                u = to!double(i-blk.imin) / to!double(blk.imax+1-blk.imin); 
                setAsWeightedSum(vel, u, &vele, 1-u, &velw);
                setAsWeightedSum(velt, u, &velet, 1-u, &velwt);

                // set velocity by interpolating in k.
                for (k = blk.kmin; k <= blk.kmax+1; ++i) {
                    w = to!double(k-blk.kmin) / to!double(blk.kmax+1-blk.kmin); 
                    setAsWeightedSum(blk.get_vtx(i,j,k).vel[0], w, &velt, 1-w, &vel);
                }
            }
        }
    }
    else {
        string errMsg = "ERROR: Wrong number of arguments passed to luafn: setVtxVelocitiesByCornersReg()\n";
        luaL_error(L, errMsg.toStringz);
    }
    // In case, the user gave use more return values than
    // we used, just set the lua stack to empty and let
    // the lua garbage collector do its thing.
    lua_settop(L, 0);
    return 0;
}

/**
 * Sets the velocity of an individual vertex.
 *
 * This function can be called for structured 
 * or unstructured grids. We'll determine what
 * type grid is meant by the number of arguments
 * supplied. The following calls are allowed:
 *
 * setVtxVelocity(vel, blkId, vtxId)
 *   Sets the velocity vector for vertex vtxId in
 *   block blkId. This works for both structured
 *   and unstructured grids.
 *
 * setVtxVelocity(vel, blkId, i, j)
 *   Sets the velocity vector for vertex "i,j" in
 *   block blkId in a two-dimensional structured grid.
 *
 * setVtxVelocity(vel, blkId, i, j, k)
 *   Set the velocity vector for vertex "i,j,k" in
 *   block blkId in a three-dimensional structured grid.
 */
extern(C) int luafn_setVtxVelocity(lua_State* L)
{
    int narg = lua_gettop(L);
    auto vel = checkVector3(L, 1);
    auto blkId = lua_tointeger(L, 2);

    if ( narg == 3 ) {
        auto vtxId = lua_tointeger(L, 3);
        globalFluidBlocks[blkId].vertices[vtxId].vel[0].set(vel);
    }
    else if ( narg == 4 ) {
        auto i = lua_tointeger(L, 3);
        auto j = lua_tointeger(L, 4);
        globalFluidBlocks[blkId].get_vtx(i,j).vel[0].set(vel);
    }
    else if ( narg >= 5 ) {
        auto i = lua_tointeger(L, 3);
        auto j = lua_tointeger(L, 4);
        auto k = lua_tointeger(L, 5);
        globalFluidBlocks[blkId].get_vtx(i,j,k).vel[0].set(vel);
    }
    else {
        string errMsg = "ERROR: Too few arguments passed to luafn: setVtxVelocity()\n";
        luaL_error(L, errMsg.toStringz);
    }
    lua_settop(L, 0);
    return 0;
}

/**
 * Sets the velocity components of an individual vertex.
 *
 * This function can be called for structured 
 * or unstructured grids. We'll determine what
 * type grid is meant by the number of arguments
 * supplied. The following calls are allowed:
 *
 * setVtxVelocityXYZ(velx, vely, velz, blkId, vtxId)
 *   Sets the velocity vector for vertex vtxId in
 *   block blkId. This works for both structured
 *   and unstructured grids.
 *
 * setVtxVelocityXYZ(velx, vely, velz, blkId, i, j)
 *   Sets the velocity vector for vertex "i,j" in
 *   block blkId in a two-dimensional structured grid.
 *
 * setVtxVelocityXYZ(velx, vely, velz, blkId, i, j, k)
 *   Set the velocity vector for vertex "i,j,k" in
 *   block blkId in a three-dimensional structured grid.
 */
extern(C) int luafn_setVtxVelocityXYZ(lua_State* L)
{
    int narg = lua_gettop(L);
    double velx = lua_tonumber(L, 1);
    double vely = lua_tonumber(L, 2);
    double velz = lua_tonumber(L, 3);
    auto blkId = lua_tointeger(L, 4);

    if ( narg == 5 ) {
        auto vtxId = lua_tointeger(L, 5);
        globalFluidBlocks[blkId].vertices[vtxId].vel[0].set(velx, vely, velz);
    }
    else if ( narg == 6 ) {
        auto i = lua_tointeger(L, 5);
        auto j = lua_tointeger(L, 6);
        globalFluidBlocks[blkId].get_vtx(i,j).vel[0].set(velx, vely, velz);
    }
    else if ( narg >= 7 ) {
        auto i = lua_tointeger(L, 5);
        auto j = lua_tointeger(L, 6);
        auto k = lua_tointeger(L, 7);
        globalFluidBlocks[blkId].get_vtx(i,j,k).vel[0].set(velx, vely, velz);
    }
    else {
        string errMsg = "ERROR: Too few arguments passed to luafn: setVtxVelocityXYZ()\n";
        luaL_error(L, errMsg.toStringz);
    }
    lua_settop(L, 0);
    return 0;
}

void assign_vertex_velocities_via_udf(double sim_time, double dt)
{
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "assignVtxVelocities");
    lua_pushnumber(L, sim_time);
    lua_pushnumber(L, dt);
    int number_args = 2;
    int number_results = 0;

    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
        string errMsg = "ERROR: while running user-defined function assignVtxVelocities()\n";
        errMsg ~= to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
}
