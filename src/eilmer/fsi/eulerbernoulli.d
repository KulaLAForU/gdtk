module fsi.eulerbernoulli;

import std.math;
import std.conv;
import std.algorithm;
import std.stdio;
import core.time;
import std.format;
import std.array;

import nm.number;
import nm.bbla;
import nm.smla;
import fsi;
import geom;

class EulerBernoulliBeam : FEMModel {
public:

    this(string jobName, int id) { super(jobName, id); }

    // Use snake case for formatting this function since it appears in the main code
    override void model_setup() {

        // Set some constants used in the master initialiser to allocate the correct amount 
        // of memory to the vectors/matrices- 2 DoFs per node
        nNodes = myConfig.Nx + 1;
        nDoF = nNodes * 2;
        nQuadPoints = myConfig.Nx * 2;

        super.model_setup();
    } // end model_setup

    // begin generateMassStiffnessMatrices
    override void GenerateMassStiffnessMatrices() {
        number l = myConfig.length / myConfig.Nx;

        // Generate the local stiffness and mass matrices- for uniform meshes, these are the same for every element
        // See function signatures for descriptions on where to find the specifications of the matrices.
        Matrix!number KL = LocalStiffnessMatrix(l);
        KL.scale(myConfig.youngsModulus * (pow(myConfig.thickness, 3) / 12) / pow(l, 3));
        Matrix!number ML = LocalMassMatrix(l);
        ML.scale(myConfig.thickness * l * myConfig.density / 420);

        // Iterate through the elements, then the nodes per element, then DoFs etc.
        size_t GlobalNodeIndx, GlobalRowIndx, GlobalColIndx, LocalRowIndx, LocalColIndx;
        foreach (i; 0 .. myConfig.Nx) {
            foreach (node; 0 .. 2) {
                GlobalNodeIndx = i + node;
                foreach (DoF; 0 .. 2) {
                    LocalRowIndx = node * 2 + DoF;
                    GlobalRowIndx = GlobalNodeIndx * 2 + DoF;
                    foreach (n_node; 0 .. 2) {
                        foreach (D_DoF; 0 .. 2) {
                            LocalColIndx = n_node * 2 + D_DoF;
                            GlobalColIndx = (i + n_node) * 2 + D_DoF;
                            M[GlobalRowIndx, GlobalColIndx] += ML[LocalRowIndx, LocalColIndx];
                            K[GlobalRowIndx, GlobalColIndx] += KL[LocalRowIndx, LocalColIndx];
                        } // end foreach D_DoF
                    } // end foreach n_node
                } // end foreach DoF
            } // end foreach node
        } // end foreach i

        // Set the boundary conditions, by setting the rows/cols corresponding
        // to fixed DoFs to diagonal ones
        foreach (ZeroedIndx; zeroedIndices) {
            foreach (DoF; 0 .. (myConfig.Nx + 1) * 2) {
                if (ZeroedIndx == DoF) {
                    K[ZeroedIndx, DoF] = 1.0;
                    M[ZeroedIndx, DoF] = 1.0;
                } else {
                    K[ZeroedIndx, DoF] = 0.0;
                    M[ZeroedIndx, DoF] = 0.0;
                    K[DoF, ZeroedIndx] = 0.0;
                    M[DoF, ZeroedIndx] = 0.0;
                }
            }
        }
    } // end GenerateMassStiffnessMatrices

    override void ConvertToNodeVel() {
        // Convert the rates of change of the degrees of freedom to velocities usable by the mesh
        foreach (node; 0 .. myConfig.Nx + 1) {
            FEMNodeVel[node].x = V[2 * node];
            FEMNodeVel[node].y = 0.0;
            FEMNodeVel[node].z = 0.0;
            FEMNodeVel[node].transform_to_global_frame(plateNormal, plateTangent1, plateTangent2);
        }
    } // end ConvertToNodeVel

    number[4] ShapeFunctionEval(number L, number x) {
        // Evaluate the shape functions, Eq 2.22 in
        // "Programming the Finite Element Method"
        // at the quadrature points
        number[4] N;
        N[0] = (1 / pow(L, 3)) * (pow(L, 3) - 3 * L * pow(x, 2) + 2 * pow(x, 3));
        N[1] = (1 / pow(L, 2)) * (pow(L, 2) * x - 2 * L * pow(x, 2) + pow(x, 3));
        N[2] = (1 / pow(L, 3)) * (3 * L * pow(x, 2) - 2 * pow(x, 3));
        N[3] = (1 / pow(L, 2)) * (pow(x, 3) - L * pow(x, 2));
        return N;
    } // end ShapeFunctionEval

    override void UpdateForceVector() {
        // Update the force vector based on the external pressures using equation 2.25 in
        // "Programming the Finite Element Method" by Smith et al.
        // with the modification that q be a function x and moced inside the integral.
        
        number l = myConfig.length / myConfig.Nx;

        // Evaluate the shape functions
        number[4] Nq1, Nq2;
        Nq1 = ShapeFunctionEval(l, (l / 2) * (-1 / sqrt(3.) + 1));
        Nq2 = ShapeFunctionEval(l, (l / 2) * (1 / sqrt(3.) + 1));

        // Perform two point gauss quadrature to compute the integral
        foreach (i; 0 .. myConfig.Nx) {
            number q1 = southPressureAtQuads[2*i] - northPressureAtQuads[2*i];
            number q2 = southPressureAtQuads[2*i+1] - northPressureAtQuads[2*i+1];

            F._data[i*2 .. (i+2)*2] += (l / 2) * (q1 * Nq1[] + q2 * Nq2[]);
        }

        // Apply the boundary conditions
        foreach (ZeroedIndx; zeroedIndices) {
            F._data[ZeroedIndx] = 0.0;
        }
    } // end updateForceVector

    // begin determineBoundaryConditions
    override void DetermineBoundaryConditions(string BCs) {
        // Determine the boundary conditions. The BC string should be 2 characters long,
        // each character denoting a boundary. The order of the BCs are "(-x)(+x)".
        // The boundary may be:
        //      F: Free, no constraints on the boundary
        //      C: Clamped, all 3 degrees of freedom are fixed to 0
        //      P: Pinned, the displacement and slope along the boundary are fixed to 0

        // Negative x
        switch (BCs[0]) {
            case 'C':
                foreach (DoF; 0 .. 2) {
                    zeroedIndices ~= DoF;
                }
                break;
            case 'P':
                zeroedIndices ~= 0;
                break;
            case 'F':
                break;
            default:
                throw new Error("Unrecognised BC specification in FSI; should be 'F', 'C' or 'P'");
        }

        // Positive x
        switch (BCs[1]) {
            case 'C':
                foreach (DoF; 0 .. 2) {
                    zeroedIndices ~= (myConfig.Nx + 1) * 2 + DoF;
                }
                break;
            case 'P':
                zeroedIndices ~= myConfig.Nx * 2;
                break;
            case 'F':
                break;
            default:
                throw new Error("Unrecognised BC specification in FSI; should be 'F', 'C' or 'P'");
        }

    } // end determineBoundaryConditions

    // begin LocalStiffnessMatrix
    Matrix!number LocalStiffnessMatrix(number l) {
        // Generate the element stiffness matrix from equation 2.26 in
        // "Programming the Finite Element Method" by Smith et al.
        // (with scaling performed in the higher function)
        Matrix!number KL = new Matrix!number(4);
        KL[0, 0] = 12; KL[0, 1] = 6 * l; KL[0, 2] = -12; KL[0, 3] = 6 * l;
        KL[1, 1] = 4 * pow(l, 2); KL[1, 2] = -6 * l; KL[1, 3] = 2 * pow(l, 2);
        KL[2, 2] = 12; KL[2, 3] = -6 * l;
        KL[3, 3] = 4 * pow(l, 2);

        foreach (i; 1 .. 4) {
            foreach (j; 0 .. i) {
                KL[i, j] = KL[j, i];
            }
        }

        return KL;
    } // end LocalStiffnessMatrix

    // begin LocalMassMatrix
    Matrix!number LocalMassMatrix(number l) {
        // Generate the element mass matrix from equation 2.30 in
        // "Programming the Finite Element Method" by Smith et al.
        // (with scaling performed in the higher level)
        Matrix!number ML = new Matrix!number(4);
        ML[0, 0] = 156; ML[0, 1] = 22 * l; ML[0, 2] = 54; ML[0, 3] = -13 * l;
        ML[1, 1] = 4 * pow(l, 2); ML[1, 2] = 13 * l; ML[1, 3] = -3 * pow(l, 2);
        ML[2, 2] = 156; ML[2, 3] = -22 * l;
        ML[3, 3] = 4 * pow(l, 2);

        foreach (i; 1 .. 4) {
            foreach (j; 0 .. i) {
                ML[i, j] = ML[j, i];
            }
        }

        return ML;
    } // end LocalMassMatrix

    // begin write
    override void WriteToFile(size_t tindx) {
        auto writeFile = File(format("FSI/t%04d.dat", tindx), "w+");
        writeFile.write("# x\ttheta_x\tdxdt\tdtheta_xdt\n");
        foreach (i; 0 .. (myConfig.Nx + 1)) {
            writeFile.write(format("%1.8e %1.8e %1.8e %1.8e\n", X[i * 2].re, X[i * 2 + 1].re, V[i * 2].re, V[i * 2 + 1].re));
        }
    }

    override void ReadFromFile(size_t tindx) {
        auto readFile = File(format("FSI/t%04d.dat", tindx), "r").byLine();
        // Pop the header line
        readFile.popFront();
        double[4] line;
        foreach (i; 0 .. (myConfig.Nx + 1)) {
            line = map!(to!double)(splitter(readFile.front())).array; readFile.popFront();
            X[i * 2 .. (i + 1) * 2] = line[0 .. 2];
            V[i * 2 .. (i + 1) * 2] = line[2 .. 4];
        }
    }
}
