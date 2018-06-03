/** 
 * adaptive_lut_CEA.d
 *
 * Gas model class for reading adaptive look-up tables generated
 * by the program, build_cea_adaptive_lut.py in directory src/gas/LUT/
 *
 * 
 * Author: James M. Burgess 06-Mar-2016
 */
module gas.adaptive_lut_CEA;

import std.stdio;
import std.math;
import std.algorithm;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;


class AdaptiveLUT : GasModel {
public:
    // Default Implementation
    this()
    {
        _n_species = 1;
        _n_modes = 0; 
        _species_names ~= "AdaptiveLUT";
        assert(_species_names.length == 1);
        create_species_reverse_lookup();
        // _mol_masses is defined at the end of the constructor  
        version(complex_numbers) {
            throw new Error("Do not use with complex numbers.");
        }
    }

    // Passing a lua-state which contains the look-up table generated by
    // program: build_cea_adaptive_lut.py
    this(lua_State *L) 
    {
        this();    // Call the default constructor
        string interp_method;
        _interp_method = getString(L, LUA_GLOBALSINDEX, "interpolation_method");
        _emin = getDouble(L, LUA_GLOBALSINDEX, "emin");
        _emax = getDouble(L, LUA_GLOBALSINDEX, "emax");     
        _log_rho_min = getDouble(L, LUA_GLOBALSINDEX, "log_rho_min");   
        _log_rho_max = getDouble(L, LUA_GLOBALSINDEX, "log_rho_max");   

        lua_getglobal(L, "data_tree");
        if ( !lua_istable(L, -1) ) {
            string msg;
            msg ~= format("AdaptiveLUT constructor\n");
            msg ~= format("   Error in look-up table input file.\n"); 
            msg ~= format("   A table of 'data' is expected, but not found.\n");
            throw new Exception(msg);
        }

        size_t np = lua_objlen(L, -1); // Determine number of patches
        for (int i = 1; i != np+1 ; ++i) {
            // For each patch in the data_tree
            lua_rawgeti(L, -1, i); // Push patch i in data_tree on to stack
                
            lua_rawgeti(L, -1, 1); // Push table 1 containing ID's/split data
            lua_rawgeti(L, -1, 1); // Push ID
            size_t ID = luaL_checkinteger(L, -1);
            lua_pop(L, 1);
            lua_rawgeti(L, -1, 2); // Push splitID
            string splitID = to!string(lua_tostring(L, -1));
            lua_pop(L, 1);
            lua_rawgeti(L, -1, 3); // Push child_left_ID
            size_t child_left_ID = luaL_checkinteger(L, -1);
            lua_pop(L, 1);
            lua_rawgeti(L, -1, 4); // Push child_right_ID
            size_t child_right_ID = luaL_checkinteger(L, -1);
            lua_pop(L, 1);
            lua_pop(L, 1); // Pop off table 1 of patch i

            lua_rawgeti(L, -1, 2); // Push table 2 of patch i containing patch coords
            lua_rawgeti(L, -1, 1); // Push lr_lo
            double lr_lo = luaL_checknumber(L, -1);
            lua_pop(L, 1);
            lua_rawgeti(L, -1, 2); // Push lr_hi
            double lr_hi = luaL_checknumber(L,-1);
            lua_pop(L, 1);
            lua_rawgeti(L, -1, 3); // Push e_lo
            double e_lo = luaL_checknumber(L, -1);
            lua_pop(L, 1);
            lua_rawgeti(L, -1, 4);
            double e_hi = luaL_checknumber(L, -1);
            lua_pop(L, 1);

            lua_pop(L, 1); // Pop off table 2 of patch i
                
            if (splitID == "n") {
                /* Patches with splitID=="n" are leaf nodes. Only leaf nodes have 
                   bezier control point data. This block reads the data differently
                   depending on whether interpolation method is bezier or linear */

                switch ( _interp_method ) {
                case "bezier":
                    lua_rawgeti(L, -1, 3); // push data table of all control points
                    // create arrays for 16 control points for all 7 properties
                    double[16] Cv_hat, Cv, R_hat, Cp_hat, gamma_hat, mu, k;
                        
                    /* Loop through all 16 control points.The '+1' in the for statement and
                       '-1' in the arary index is because lua indexing starts at 1 */
                        
                    for(int j= 1; j<16+1; j++) { 
                        lua_rawgeti(L, -1, j); // Push lua_table of control point j 
                            
                        lua_rawgeti(L, -1, 1); // Push control point j for Cv_hat
                        Cv_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1); //Pop cp j, Cv_hat
                        lua_rawgeti(L, -1, 2); // Push CP j for Cv
                        Cv[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 3); // Push CP j for R_hat
                        R_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 4); // Push CP j for Cp_hat
                        Cp_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 5); // Push CP j for gamma_hat
                        gamma_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 6); // Push CP j for mu
                        mu[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 7); // Push Cp j for k
                        k[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                            
                        lua_pop(L, 1); // Pop lua_table of control point j
                    } // End for loop (through all control points)
                    lua_pop(L, 1); // Pop 2D lua_table of all control points
                        
                    _tree[ID] = new BezierPatch( lr_lo,  lr_hi,  e_lo, e_hi, ID, 
                                                 child_left_ID, child_right_ID, splitID,
                                                 Cv_hat, Cv, R_hat, Cp_hat, gamma_hat, mu, k);
                    break; // End case 'bezier' for leaf nodes of binary tree
                    
                case "linear":

                    lua_rawgeti(L, -1, 3); // push data table of corner properties
                        
                    double[4] Cv_hat, Cv, R_hat, Cp_hat, gamma_hat, mu, k;
                
                    for(int j= 1; j<4+1; j++) { 
                        lua_rawgeti(L, -1, j); // Push lua_table of data at point j
                            
                        lua_rawgeti(L, -1, 1); // Push  j for Cv_hat
                        Cv_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1); //Pop cp j, Cv_hat
                        lua_rawgeti(L, -1, 2); // Push point j for Cv
                        Cv[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 3); // Push point j for R_hat
                        R_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 4); // Push point j for Cp_hat
                        Cp_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 5); // Push point j for gamma_hat
                        gamma_hat[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 6); // Push point j for mu
                        mu[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);
                        lua_rawgeti(L, -1, 7); // Push point j for k
                        k[j-1] = luaL_checknumber(L, -1);
                        lua_pop(L, 1);

                        lua_pop(L, 1); // Pop lua_table of control point j
                    }
                        
                    lua_pop(L, 1); // Pop 2D lua_table of all properties
                    _tree[ID] = new LinearPatch( lr_lo,  lr_hi,  e_lo, e_hi, ID, 
                                                 child_left_ID, child_right_ID, splitID,
                                                 Cv_hat, Cv, R_hat, Cp_hat, gamma_hat, mu, k);
                    break; // end case 'linear' for leaf nodes of binary tree

                default:
                    string msg;
                    msg ~= "No valid value of interpolation_method could be ";
                    msg ~= "found in adaptive lut file. Valid options are ";
                    msg ~= "'bezier' or 'linear'.";
                    throw new Exception(msg);
                } // End switch/case block
                
            } // End if statement (writing data to leaf nodes - patches with data)
                    
            else if (splitID == "e" || splitID == "lr") {
                /* Non-leaf nodes still need to be in the tree, but they don't need cntrl
                   point data - so they are Patch (non BezierPatch objects) */
                        
                _tree[ID] = new Patch( lr_lo,  lr_hi,  e_lo,  e_hi, ID, child_left_ID,
                                       child_right_ID, splitID);
            }
            else {
                string msg;
                msg ~= "Constructing the AdaptiveLUT gas model failed in the constructor ";
                msg ~= "the patches because splitID of patch ID " ~ to!string(ID);
                msg ~= " was not defined as 'n', 'lr' or 'e'.";
                throw new Exception(msg);                      
            }

            lua_pop(L, 1); // Pop patch i in data tree 
        }
            
        // The patches have all been read - assign pointers to 'child nodes'
        foreach(ID, ref node; _tree) {
            if (node._splitID != "n") {
                node._child_left = &_tree[node._child_left_ID];
                node._child_right = &_tree[node._child_right_ID];
            }
            // The interger address of the child nodes are no longer needed
            // since we are using pointers instead
            destroy(node._child_right_ID);
            destroy(node._child_left_ID);
        }

        // Save the first patch as the root node as a starting point for tree searches
        this._root_node = &_tree[0];

        // Save the bounds of the first patch (which is at _tree[1]) as global table bounds
        // as a static variable in the patch class
        Patch._lr_lo_global = _tree[0]._lr_lo;
        Patch._lr_hi_global = _tree[0]._lr_hi;
        Patch._e_lo_global = _tree[0]._e_lo;
        Patch._e_hi_global = _tree[0]._e_hi;

        Patch._de_global = Patch._e_hi_global - Patch._e_lo_global;
        Patch._dlr_global = Patch._lr_hi_global - Patch._lr_lo_global;
              
        // Read entropy reference conditions
        with_entropy = getInt(L, LUA_GLOBALSINDEX, "with_entropy");
        assert(with_entropy == 1, "Error in adaptive lut file");   
        _p1 = getDouble(L, LUA_GLOBALSINDEX, "p1");
        _T1 = getDouble(L, LUA_GLOBALSINDEX, "T1");
        _s1 = getDouble(L, LUA_GLOBALSINDEX, "s1");

        lua_pop(L, 1); // Finish reading the table
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "AdaptiveLUT =(";
        repr ~= "Interpolation_method = " ~ _interp_method;
        repr ~= ", with_entropy = 1";
        repr ~= ", entropy reference properties ={";
        repr ~= "p1 = " ~ to!string(_p1);
        repr ~= ", T1 = " ~ to!string(_T1);
        repr ~= ", s1 = " ~ to!string(_s1);
        repr ~= "}, valid interpolation range ={";
        repr ~= "e_min = " ~ to!string(_emin); 
        repr ~= ", e_max = " ~ to!string(_emax);
        repr ~= ", log_rho_min = " ~ to!string(_log_rho_min);
        repr ~= ", log_rho_max = " ~ to!string(_log_rho_max);
        repr ~= "}}";
        return to!string(repr);
    }
        
    const(Patch)* search_tree(in double lr, in double e) const
    {
        // Start search at patch with ID = 1 that covers the whole tree
        // The patches are considered nodes of the binary search tree
        const(Patch) *node = _root_node;
        while ((*node)._splitID != "n") {
            if ((*node)._splitID == "e") {
                if (e < (*(*node)._child_left)._e_hi ) {
                    node = (*node)._child_left;
                } else {
                    node = (*node)._child_right;
                }
            }
            else if ((*node)._splitID == "lr") {
                // SplitID must be "lr" because contructor asserts is "n", "e" or "lr"
                if (lr < (*(*node)._child_left)._lr_hi ) {
                    node = (*node)._child_left;
                } else {
                    node = (*node)._child_right;
                }
            }
        }
        return node;
    }
 
    override void update_thermo_from_rhou(GasState Q) {
        double Cv_eff, R_eff, g_eff;

        double lr = log10(Q.rho.re);
        const(Patch)* node = this.search_tree(lr, Q.u.re);

        Cv_eff = node.interpolate(lr, Q.u.re, "Cv_hat");
        R_eff =  node.interpolate(lr, Q.u.re, "R_hat"); 
        g_eff = node.interpolate(lr, Q.u.re, "gamma_hat");

        // Reconstruct the thermodynamic properties.
        Q.T = Q.u / Cv_eff;
        Q.p = Q.rho*R_eff*Q.T;
        Q.a = sqrt(g_eff*R_eff*Q.T);

        // Fix meaningless values if they arise
        if ( Q.p < 0.0 ) Q.p = 0.0;
        if ( Q.T < 0.0 ) Q.T = 0.0;
        if ( Q.a < 0.0 ) Q.a = 0.0;

    }

    override void update_sound_speed(GasState Q){ 
        if (isNaN(Q.rho) || isNaN(Q.u)) {
            string err;
            err ~= format("update_sound_speed() method for LUT requires e and rho of ");
            err ~= format("GasState Q to be defined: rho = %.5s, e = .8s", Q.rho, Q.u);
            throw new Exception(err);
        }

        double lr = log10(Q.rho.re);
        const(Patch)* node = this.search_tree(lr, Q.u.re);

        double R_eff = node.interpolate(lr, Q.u.re, "R_hat");
        double g_eff = node.interpolate(lr, Q.u.re, "gamma_hat");

        // Reconstruct the thermodynamic properties.
        Q.a = sqrt(g_eff*R_eff*Q.T);
    }

    override void update_trans_coeffs(GasState Q) const
    {
        double lr = log10(Q.rho.re);
        const(Patch)* node = this.search_tree(lr, Q.u.re);
        
        Q.mu = node.interpolate(lr, Q.u.re, "mu");
        Q.k = node.interpolate(lr, Q.u.re, "k"); 
    }
    
    override void update_thermo_from_pT(GasState Q) {
        update_thermo_state_pT(this, Q);
    }
    
    override void update_thermo_from_rhoT(GasState Q) {
        update_thermo_state_rhoT(this, Q);
    }

    override void update_thermo_from_rhop(GasState Q) {
        update_thermo_state_rhop(this, Q);
    }

    override void update_thermo_from_ps(GasState Q, number s) {
        update_thermo_state_ps(this, Q, s);
    }
    
    override void update_thermo_from_hs(GasState Q, number h, number s) {
        update_thermo_state_hs(this, Q, h, s); 
    }
    

    // Factor for computing forward difference derivatives
    immutable double Epsilon = 1.0e-6;

    override number dpdrho_const_T(in GasState Q)
    {
        // Cannot be directly computed from look up table
        // Apply forward finite difference in calling update from p, rho
        // Make a copy of the GasState object, to avoid modifying Q
        auto Q_temp = new GasState(this);
        Q_temp.copy_values_from(Q);

        // Save values for derivative calculation
        double p = Q_temp.p.re;
        double rho = Q_temp.rho.re;
        double drho = Epsilon * rho; // Step size for varying rho
        Q_temp.rho = rho + drho; // Modify density slightly

        // Evaluate gas state due to perturbed rho, holding constant T
        this.update_thermo_from_rhoT(Q_temp);
        double p_step = Q_temp.p.re;
        return to!number((p_step - p) / drho);
    }



    // const void update_diff_coeffs(ref GasState Q) {}

    override number dudT_const_v(in GasState Q) 
    {
        double lr = log10(Q.rho.re);
        const(Patch)* node = this.search_tree(lr, Q.u.re);
        
        double Cv_actual = node.interpolate(lr, Q.u.re, "Cv");
        return to!number(Cv_actual);
    }
    override number dhdT_const_p(in GasState Q) 
    { 
        double lr = log10(Q.rho.re);
        const(Patch)* node = this.search_tree(lr, Q.u.re);
        
        double Cv_actual = node.interpolate(lr, Q.u.re, "Cv"); 
        double R_eff = node.interpolate(lr, Q.u.re, "R_hat"); 

        return to!number(Cv_actual + R_eff);
    }
   
    override number gas_constant(in GasState Q)  
    {
        double lr = log10(Q.rho.re);
        const(Patch)* node = this.search_tree(lr, Q.u.re);
        
        double R_eff = node.interpolate(lr, Q.u.re, "R_hat"); 
        return to!number(R_eff);
    }
    
    override number internal_energy(in GasState Q)  
    {
        // This method should never be called expecting quality data
        // because the LUT gas doesn not keep the species information.
        // This implementation is here to keep D happy that
        // all of the methods are implemented as required,
        // and there may be times when having this function
        // return something reasonable may make other code
        // simpler because it doesn't have to treat the
        // LUT gas specially.
        // This is the one-species implementation

        return Q.u;
    }
    override number enthalpy(in GasState Q)  
    {
        // This method assumes that the internal energy,
        // pressure and density are up-to-date in the
        // GasData object. Then enthalpy is computed
        // from definition.

        number h = Q.u + Q.p/Q.rho;
        return h;
    }

    override number entropy(in GasState Q)  {
        assert(with_entropy == 1);
        // Eilmer3 implementation uses ideal gas as a backup for calculating 
        // entropy if the table is without entropy data. 
        // Since the table-generating script always produces entropy info,
        // this function does not do this

        double lr = log10(Q.rho.re);
        const(Patch)* node = this.search_tree(lr.re, Q.u.re);
        
        double Cv_eff = node.interpolate(lr, Q.u.re, "Cv_hat");
        double R_eff = node.interpolate(lr, Q.u.re, "R_hat"); 
        double Cp_eff = node.interpolate(lr, Q.u.re, "Cp_hat"); 

        double T = Q.u.re / Cv_eff;
        double p = Q.rho.re * R_eff * T;
        double s = _s1 + Cp_eff*log(T/_T1) - R_eff*log(p/_p1); 

        return to!number(s);
    }

    double s_molecular_weight(int isp) const
    { 
        // This method is not very meaningful for an equilibrium
        // gas.  The molecular weight is best obtained from
        // the mixture molecular weight methods which IS a function
        // of gas composition and thermodynamic state, however,
        // there are times when a value from this function makes
        // other code simpler, in that the doesn't have to treat
        // the look-up gas specially.
        if ( isp != 0 ) {
            throw new Exception("LUT gas: should not be looking up isp != 0");
        }

        // Find some cold value in the table - these are the low limits of patch 1
        double lr = _tree[1]._lr_lo;
        double e = _tree[1]._e_lo;

        const(Patch)* node = this.search_tree(lr, e);
        
        double Rgas = node.interpolate(lr, e, "R_hat");// J/kg/deg-K
        double M = R_universal / Rgas;    
        return M;
    }
       
private:
    int with_entropy;
    double _p1, _T1, _s1;
    Patch[size_t] _tree;
    Patch *_root_node;
    double _emin, _emax, _log_rho_min, _log_rho_max;
    string _interp_method; // For toString() method

} // End of AdaptiveLUT class



// The Patch class is used by the AdaptiveLUT model for patches in the binary
// search tree that are not leaf nodes.
// They contain infomration about their bounds and references to their child
// patches. Patches that are leaf nodes, and actually have property data
// inherit from Patch: they are BezierPatch and LinearPatch, depending on
// the interpolation method
class Patch {
public:

    this( double lr_lo,  double lr_hi,  double e_lo,  double e_hi,  size_t ID,
          size_t child_left_ID,  size_t child_right_ID, string splitID)
    {
        this._lr_lo = lr_lo;
        this._lr_hi = lr_hi;
        this._e_lo = e_lo;
        this._e_hi = e_hi;
        this._ID = ID;
        this._child_left_ID = child_left_ID;
        this._child_right_ID = child_right_ID;
        this._splitID = splitID;
    }

    void check_interpolation_bounds(ref double lr, ref double  e) const
    {
        // Deal with the possibility of data being outside the table
        // If the requested value is outside the bounds of the table, the search
        // function return the nearest patch, if we move directly toward the table
        // We will assign any value outside the table as the projected value
        // onto the  edge
        if (lr < _lr_lo ) {
            // Data is to the left of the table
            lr = _lr_lo;
        } else if  (lr > _lr_hi) {
            // Data is to the right of the table
            lr = _lr_hi;
        } 
        
        if (e < _e_lo) {
            // Data is below the table
            e = _e_lo;
        } else if (e > _e_hi) {
            // Data is above table
            e = _e_hi;
        }
    }

    double interpolate( double lr,  double e, string prop) const
    {
        // Interpolate functions to be overridden
        // The error message is the same for all these functions
        // The message explains what the problem is
        string msg;
        msg ~= "The interpolate function is been called on a Patch object. ";
        msg ~= "This should only be called on leaf nodes of the search tree, ";
        msg ~= "which must be of the inherited class BezierPatch or LinearPatch.";
        throw new Exception(msg);
    }
      
private:
    double _lr_lo, _lr_hi;
    double _e_lo, _e_hi;
    Patch *_child_left;
    Patch *_child_right;
    size_t _ID, _child_left_ID, _child_right_ID;
    string _splitID;
    static double _e_lo_global, _e_hi_global, _de_global;
    static double _lr_lo_global, _lr_hi_global, _dlr_global;

}


class BezierPatch : Patch
{
public:
    this( double lr_lo,  double lr_hi,  double e_lo,  double e_hi,  size_t ID,
          size_t child_left_ID,  size_t child_right_ID, string splitID,  double[16] Cv_hat,
          double[16] Cv,  double[16] R_hat,  double[16] Cp_hat, 
          double[16] gamma_hat,  double[16] mu,  double[16] k)
    {   
        super(lr_lo, lr_hi, e_lo, e_hi, ID, child_left_ID, 
              child_right_ID, splitID); // Call parent class contructor (Patch) 
        this._Cv_hat = this.map_control_points(Cv_hat);
        this._Cv =  this.map_control_points(Cv);
        this._R_hat =  this.map_control_points(R_hat);
        this._Cp_hat =  this.map_control_points(Cp_hat);
        this._gamma_hat =  this.map_control_points(gamma_hat);
        this._mu =  this.map_control_points(mu);
        this._k =  this.map_control_points(k);
    }


    double[4][4][4] map_control_points(in double[16] bezier_cps)
    {
        /*  Maps the control points to a form ready for the interpolation algorithm
            The 16 points are on a 4 by 4 grid nested inside a 4-elemnent array.
            The grid is in the last element of that array, and the remaining points are 0.
            The outer array is used for the iteration procedure */
        double[4][4][4] bs = 0;
        size_t i;
        foreach(ref row; bs[3]){
            foreach(ref element; row){
                element = bezier_cps[i];
                i++;
            }
        }
        return bs;
    }

    
    override double interpolate(double lr, double e, string prop) const
    {
        // Check if the values requrested are oustide the table
        check_interpolation_bounds(lr.re,  e.re);

        // Now the points will be bounded by the patch
        double u = (lr - _lr_lo) / (_lr_hi - _lr_lo);
        double v = (e - _e_lo) / (_e_hi - _e_lo);

        double[4][4][4] data;
        final switch (prop) {
        case "Cv_hat":
            data = _Cv_hat;
            break;
        case "Cv":
            data = _Cv;
            break;
        case "R_hat":
            data = _R_hat;
            break;
        case "Cp_hat":
            data = _Cp_hat;
            break;
        case "gamma_hat":
            data = _gamma_hat;
            break;
        case "mu":
            data = _mu;
            break;
        case "k":
            data = _k;
            break;
        }
        for(int k = 2; k != -1; k--){
            for(int i = 0; i != k+1; i++){
                for(int j = 0; j != k+1; j++){
                    data[k][j][i] = (1-u)*(1-v)*data[k+1][j][i] + u*(1-v)*data[k+1][j][i+1] +
                        (1-u)*v*data[k+1][j+1][i] + u*v*data[k+1][j+1][i+1];
                }
            }
        }
        return data[0][0][0];           
    }

private:
    /* 16 bezier control points per property, stored in the array form that is 
       required for the interpolation algorithm to read */
    double[4][4][4] _Cv_hat, _Cv, _R_hat;
    double[4][4][4] _Cp_hat, _gamma_hat, _mu, _k;
}



class LinearPatch : Patch 
{
public:
    this(double lr_lo, double lr_hi, double e_lo, double e_hi, size_t ID,
         size_t child_left_ID, size_t child_right_ID, string splitID,  double[4] Cv_hat,
         double[4] Cv,  double[4] R_hat,  double[4] Cp_hat, 
         double[4] gamma_hat,  double[4] mu,  double[4] k )
    {   
        super(lr_lo, lr_hi, e_lo, e_hi, ID, child_left_ID, 
              child_right_ID, splitID);
        this._Cv_hat = Cv_hat;
        this._Cv = Cv;
        this._R_hat = R_hat;
        this._Cp_hat = Cp_hat;
        this._gamma_hat = gamma_hat;
        this._mu = mu;
        this._k = k;
        
        
    }

    override double interpolate( double lr,  double e, string prop) const
    {
        // Check if the values requrested are oustide the table
        check_interpolation_bounds(lr.re,  e.re);
 
        // Now the points are bounded by the patch
        double[4] data;
        final switch (prop) {
        case "Cv_hat":
            data = _Cv_hat;
            break;
        case "Cv":
            data = _Cv;
            break;
        case "R_hat":
            data = _R_hat;
            break;
        case "Cp_hat":
            data = _Cp_hat;
            break;
        case "gamma_hat":
            data = _gamma_hat;
            break;
        case "mu":
            data = _mu;
            break;
        case "k":
            data = _k;
            break;
        }
        
        // Properties are in this order:
        // (lr_lo, e_lo), (lr_hi, e_lo),(lr_lo, e_hi),(lr_hi, e_hi) 

        double lrfrac = (lr - _lr_lo) / (_lr_hi - _lr_lo);
        double efrac = (e - _e_lo) / (_e_hi - _e_lo);
         
        double res = (1.0 - efrac) * (1.0 - lrfrac) * data[0] +
            efrac         * (1.0 - lrfrac)          * data[2] +
            efrac         * lrfrac                  * data[3] +
            (1.0 - efrac) * lrfrac                  * data[1];

        return res;
    }

private:
    double[4] _Cv_hat, _Cv, _R_hat; // 7 properties at patch corners
    double[4] _Cp_hat, _gamma_hat, _mu, _k;
}




// These test functions are the same as in uniform_lut.d
version(adaptive_lut_CEA_test) 
{
    import util.msg_service;
    
    int main() {
        GasModel gm;
        
        try {
            lua_State* L = init_lua_State();
            doLuaFile(L, "sample-data/cea-adaptive-lut-air-bezier.lua");
            gm = new AdaptiveLUT(L);
        }
        catch (Exception e) {
            writeln(e.msg);
            string msg;
            msg ~= "Test of look up table in adaptive_lut_CEA.d require file:";
            msg ~= "cea-adaptive-lut-air-bezier.lua ";
            msg ~= " in directory: gas/sample_data";
            throw new Exception(msg);
        }
   
        // An arbitrary state was defined for 'Air', massf=1, in CEA2
        // using the utility cea2_gas.py
        double p_given = 1.0e6; // Pa
        double T_given = 1.0e3; // K
        double rho_given = 3.4837; // kg/m^^3
        // CEA uses a reference temperature of 298K (Eilmer uses 0K) so the
        // temperature was offset by amount e_offset 
        double e_CEA =  456600; // J/kg
        double e_offset = 303949.904; // J/kg
        double e_given = e_CEA + e_offset; // J/kg
        double h_CEA = 743650; // J/kg
        double h_given = h_CEA + e_offset; // J/kg 
        double a_given = 619.2; // m/s
        double s_given = 7475.7; // J(kg.K)
        double R_given = 287.036; // J/(kg.K)
        double gamma_given = 1.3866;
        double Cp_given = 1141; // J/(kg.K)
        double mu_given = 4.3688e-05; // Pa.s
        double k_given = 0.0662; // W/(m.K)
        double Cv_given = e_given / T_given; // J/(kg.K)
        
        
        auto Q = new GasState(gm, p_given, T_given);
        // Return values not stored in the GasState
        // The constructor of the gas state will call update_thermo_from_pT
        // which itself calls update_thermo_from_rhou, so we are testing both
        double Cv = gm.dudT_const_v(Q);
        double Cp = gm.dhdT_const_p(Q);
        double R = gm.gas_constant(Q);
        double h = gm.enthalpy(Q);
        double s = gm.entropy(Q);
        
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(approxEqual(e_given, Q.u, 1.0e-4), failedUnitTest());
        assert(approxEqual(rho_given, Q.rho, 1.0e-4), failedUnitTest());
        assert(approxEqual(a_given, Q.a, 1.0e-4), failedUnitTest());
        assert(approxEqual(Cp_given, Cp, 1.0e-3), failedUnitTest());
        assert(approxEqual(h_given, h, 1.0e-4), failedUnitTest());
        assert(approxEqual(mu_given, Q.mu, 1.0e-4), failedUnitTest());
        assert(approxEqual(k_given, Q.k, 1.0e-4), failedUnitTest());
        assert(approxEqual(s_given, s, 1.0e-4), failedUnitTest());
        assert(approxEqual(R_given, R, 1.0e-4), failedUnitTest());
        

        
        // Now do the same test for the linear implementation
        try {
            lua_State* L = init_lua_State();
            doLuaFile(L, "sample-data/cea-adaptive-lut-air-linear.lua");
            gm = new AdaptiveLUT(L);
        }
        catch (Exception e) {
            writeln(e.msg);
            string msg;
            msg ~= "Test of look up table in adaptive_lut_CEA.d require file:";
            msg ~= "cea-adaptive-lut-air-linear.lua";
            msg ~= " in directory: gas/sample_data";
            throw new Exception(msg);
        }
        
        Q = new GasState(gm, p_given, T_given);

        // Note that the unit tests are to a higher error tolerance
        // because the look-up table was generated higher error allowed
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(approxEqual(e_given, Q.u, 1.0e-3), failedUnitTest());
        assert(approxEqual(rho_given, Q.rho, 1.0e-3), failedUnitTest());
        assert(approxEqual(a_given, Q.a, 1.0e-3), failedUnitTest());
        assert(approxEqual(Cp_given, Cp, 1.0e-3), failedUnitTest());
        assert(approxEqual(h_given, h, 1.0e-3), failedUnitTest());
        assert(approxEqual(mu_given, Q.mu, 1.0e-3), failedUnitTest());
        assert(approxEqual(k_given, Q.k, 1.0e-3), failedUnitTest());
        assert(approxEqual(s_given, s, 1.0e-3), failedUnitTest());
        assert(approxEqual(R_given, R, 1.0e-3), failedUnitTest());
        
        return 0;
    }
}
  
