/*
    Prototyping for vibrationally specific CO gas model, intended for use in simulating 
    gas dynamic lasers.

    Notes:
    "Kinetic Modelling of the High-Power Carbon Monoxide Laser*"
    Joseph W. Rich, Journal of Applied Physics, Volume 42, Number 7, June 1971

    @author: Nick Gibbons (24/03/22)
*/

module gas.vib_specific_co;

import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;

import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;

immutable uint _n_vibe_states = 10;

class VibSpecificCO: GasModel {
    this(lua_State *L)
    {
        type_str = "VibSpecificCO";
        _n_modes = 0;
        _n_species = _n_vibe_states;
        _species_names.length = _n_vibe_states;
        foreach(i; 0 .. _n_vibe_states){
            _species_names[i] = format("CO-%02d", i);
        }
        create_species_reverse_lookup();

        _mol_masses.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _mol_masses[isp] = _M;
        }

        //gs = new GasState(1, 0); // Fake gas state for interfacing with transport properties.
        //mTransProps = new GasMixtureTransProps(L, ["CO"]);

        // The energy levels of the individual vibration levels are constant, once computed.
        _vib_energy.length = _n_vibe_states;
        _vib_energy_per_kg.length = _n_vibe_states;
        foreach (i; 0 .. _n_vibe_states) {
            _vib_energy[i] = vib_energy(i);
            _vib_energy_per_kg[i]= (Avogadro_number/_M) * _vib_energy[i];
        }

        // Let's not bother with diffusion stuff at the moment.
    }

    this(string gas_file_name){
        lua_State* L = init_lua_State();
        doLuaFile(L, gas_file_name);
        this(L);
        lua_close(L);
    }

    override string toString() const
    {
        char[] repr;
        repr ~= format("VibSpecificCO=(numVibLevels=%d)", _n_vibe_states);
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const
    {
        Q.rho = Q.p/(_R*Q.T);
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_rhou(GasState Q) const
    {
        // From internal energy, remove vibrational energy before computing trans-rotational temperature.
        number u = Q.u;
        foreach (i; 0 .. _n_vibe_states) { u -= Q.massf[i]*_vib_energy_per_kg[i]; }
        Q.T = (0.4/_R)*u;
        Q.p = Q.rho*_R*Q.T;
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
        Q.p = Q.rho*_R*Q.T;
        // Start with trans-rotational component of internal energy and add vibrational energy.
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
        Q.T = Q.p/(Q.rho*_R);
        // Start with trans-rotational component of internal energy and add vibrational energy.
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_ps(GasState Q, number s) const
    {
        number sum=0.0;
        foreach(v; 0 .. _n_vibe_states) sum += Q.massf[v]*log(Q.massf[v]);
        Q.T = Tref*exp((s+_R*log(Q.p/pref) + _R*sum)/_cp);
        Q.rho = Q.p/(_R*Q.T);
        Q.u = internal_energy(Q);
    }
    override void update_thermo_from_hs(GasState Q, number h, number s) const
    {
        throw new Error("VibSpecificCO:update_thermo_from_hs NOT IMPLEMENTED.");
    }
    override void update_sound_speed(GasState Q) const
    {
        Q.a = sqrt(_gamma*_R*Q.T);
    }
    override void update_trans_coeffs(GasState Q)
    {
        // FIXME: Put the fake gas state stuff here
        //mTransProps.updateTransProps(gs);
        //throw new Error("ERROR: VibSpecificCO:update_trans_coeffs NOT IMPLEMENTED.");
        Q.mu = 0.0; Q.k = 0.0;
    }
    override number dudT_const_v(in GasState Q) const
    {
        return to!number(_cv);
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return to!number(_cp);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        return _R*Q.T;
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(_R);
    }
    override number internal_energy(in GasState Q) const
    {
        // For full internal energy, start with trans-rotational mode
        // and then add vibrational modes.
        number u = _cv*Q.T;
        foreach (i; 0 .. _n_vibe_states) { u += Q.massf[i]*_vib_energy_per_kg[i]; }
        return u;
    }
    override number enthalpy(in GasState Q) const
    {
        number h = _cp*Q.T;
        foreach (i; 0 .. _n_vibe_states) { h += Q.massf[i]*_vib_energy_per_kg[i]; }
        return h;
    }
    override number entropy(in GasState Q) const
    {
    /*
        Not completely sure about this. See 06/04/22 derivation. I think I can get away with
        setting s @ 300K and 1 atm to zero, which eliminates the reference entropy.
    */
        number s=0.0;
        foreach(i; 0 .. n_vibe_states){
            number Xi = Q.massf[i]; // Molefraction = mass fraction in CO only mode
            number si = _cp*log(Q.T/Tref) - _R*log(Xi*Q.p/pref);
            s += Q.massf[i]*si;
        }
        return s;
    }

    @nogc uint n_vibe_states() const { return _n_vibe_states; }

    @nogc const
    number boltzmann_eq_population_fraction(size_t v, number T){
    /*
        Returns the equilibrium population fraction for quantum-level v, given temperature.
        v==0 is the ground state.
    */
        if (v >= _n_vibe_states) return to!number(0.0);

        number summ = 0.0;
        foreach(ej; 0 .. _n_vibe_states) { summ += exp(-_vib_energy[ej] / (Boltzmann_constant*T)); }
        number nvf = exp(-_vib_energy[v]/(Boltzmann_constant*T)) / summ;
        return nvf;
    }
private:
    // General Gas Stuff
    //GasState gs;
    //TransportPropertiesModel mTransProps;
    immutable double Tref = 300.0;
    immutable double pref = 101.325e3;
    immutable double Ru = 8.31446261815324; // Exact as of 2019 SI redefinition
    immutable double _M = 0.0280101;
    immutable double _R = Ru/_M;
    immutable double _Cv = 5.0/2.0*Ru;
    immutable double _Cp = _Cv + 1.0*Ru;
    immutable double _cv = _Cv/_M;
    immutable double _cp = _Cp/_M;
    immutable double _gamma = _cp/_cv;

    // State Specific Constants 
    immutable double E01 = 4.31e-13*1e-7; // erg -> J Ground state energy (FIXME: Is this E1? A typo???)
    immutable double d   = 0.00598;       //          CO Anharmonicity

    // Internal storage
    double[] _vib_energy; // quantum level energies
    double[] _vib_energy_per_kg;

    @nogc const double vib_energy(uint v){
    // Returns the quantum level energy for mode v, subtracting off the zero-point so that E(v=0)=0
        return E01*(v - d*v*(v+1.0));
    }

    // There used to be a Tvib function here but I don't think we need it
} // end class VibSpecificCO

//class EquilibriumVibSpecificCO: VibSpecificCO {
//    this(lua_State *L)
//    {
//        type_str = "EquilibriumVibSpecificCO";
//        _n_modes = 0;
//        _n_species = 1;
//        _species_names.length = 1;
//        _species_names[0] = "CO";
//        create_species_reverse_lookup();
//
//        _mol_masses.length = 1;
//        _mol_masses[0] = _M;
//
//        // Fake gas state for  with transport properties.
//        gs = new GasState(_n_vibe_states, 0); 
//        vsco = new VibSpecificCO(L);
//
//        //gs = new GasState(1, 0); // Fake gas state for interfacing with transport properties.
//        //mTransProps = new GasMixtureTransProps(L, ["CO"]);
//
//        // The energy levels of the individual vibration levels are constant, once computed.
//        _vib_energy.length = _n_vibe_states;
//        _vib_energy_per_kg.length = _n_vibe_states;
//        foreach (i; 0 .. _n_vibe_states) {
//            _vib_energy[i] = vib_energy(i);
//            _vib_energy_per_kg[i]= (Avogadro_number/_M) * _vib_energy[i];
//        }
//        // Let's not bother with diffusion stuff at the moment.
//    }
//
//    this(string gas_file_name){
//        lua_State* L = init_lua_State();
//        doLuaFile(L, gas_file_name);
//        this(L);
//        lua_close(L);
//    }
//
//    override string toString() const
//    {
//        char[] repr;
//        repr ~= format("EquilibriumVibSpecificCO=(numVibLevels=%d)", _n_vibe_states);
//        return to!string(repr);
//    }
//
//    override void update_thermo_from_pT(GasState Q)
//    {
//        gs.p = Q.p;
//        gs.T = Q.T;
//        foreach(v; 0 .. _n_vibe_states) gs.massf[v] = boltzmann_eq_population_fraction(v, Q.T);
//        vsco.update_thermo_from_pT(gs);
//        Q.rho = gs.rho;
//        Q.u = gs.u;
//    }
//    override void update_thermo_from_rhou(GasState Q)
//    {
//        gs.rho = Q.rho;
//        gs.u = Q.u;
//        foreach(v; 0 .. _n_vibe_states) gs.massf[v] = boltzmann_eq_population_fraction(v, Q.T);
//
//        // From internal energy, remove vibrational energy before computing trans-rotational temperature.
//        number u = Q.u;
//        foreach (i; 0 .. _n_vibe_states) { u -= Q.massf[i]*_vib_energy_per_kg[i]; }
//        Q.T = (0.4/_R)*u;
//        Q.p = Q.rho*_R*Q.T;
//    }
//    override void update_thermo_from_rhoT(GasState Q)
//    {
//        gs.rho = Q.rho;
//        gs.T = Q.T;
//        foreach(v; 0 .. _n_vibe_states) gs.massf[v] = boltzmann_eq_population_fraction(v, Q.T);
//        vsco.update_thermo_from_rhoT(gs);
//        Q.p = gs.p;
//        Q.u = gs.u;
//    }
//    override void update_thermo_from_rhop(GasState Q) const
//    {
//        Q.T = Q.p/(Q.rho*_R);
//        // Start with trans-rotational component of internal energy and add vibrational energy.
//        Q.u = 2.5*_R*Q.T;
//        foreach (i; 0 .. _n_vibe_states) { Q.u += Q.massf[i]*_vib_energy_per_kg[i]; }
//    }
//    override void update_thermo_from_ps(GasState Q, number s) const
//    {
//        throw new Error("VibSpecificCO:update_thermo_from_ps NOT IMPLEMENTED.");
//    }
//    override void update_thermo_from_hs(GasState Q, number h, number s) const
//    {
//        throw new Error("VibSpecificCO:update_thermo_from_hs NOT IMPLEMENTED.");
//    }
//    override void update_sound_speed(GasState Q) const
//    {
//        Q.a = sqrt(_gamma*_R*Q.T);
//    }
//    override void update_trans_coeffs(GasState Q)
//    {
//        // FIXME: Put the fake gas state stuff here
//        //mTransProps.updateTransProps(gs);
//        throw new Error("ERROR: VibSpecificCO:update_trans_coeffs NOT IMPLEMENTED.");
//    }
//    override number dudT_const_v(in GasState Q) const
//    {
//        return to!number(_R/(_gamma - 1.0));
//    }
//    override number dhdT_const_p(in GasState Q) const
//    {
//        throw new Error("ERROR: VibSpecificNitrogen:dhdT_const_p NOT IMPLEMENTED.");
//    }
//    override number dpdrho_const_T(in GasState Q) const
//    {
//        throw new Error("ERROR: VibSpecificNitrogen:dpdrho_const_T NOT IMPLEMENTED.");
//    }
//    override number gas_constant(in GasState Q) const
//    {
//        return to!number(_R);
//    }
//    override number internal_energy(in GasState Q) const
//    {
//        // For full internal energy, start with trans-rotational mode
//        // and then add vibrational modes.
//        gs.T = Q.T;
//        boltzmann_eq_population_fractions(gs.massf, Q.T);
//        return vsco.internal_energy(gs);
//    }
//    override number enthalpy(in GasState Q) const
//    {
//        gs.T = Q.T;
//        boltzmann_eq_population_fractions(gs.massf, Q.T);
//        return vsco.enthalpy(gs);
//    }
//    override number entropy(in GasState Q) const
//    {
//        /* This is different and kind of weird compared to the nonequilibrium case */
//        return 0.0;
//    }
//
//private:
//    GasState gs;
//    VibSpecificCO vsco;
//
//    @nogc const
//    void boltzmann_eq_population_fractions(ref number[] nvf, number T){
//    /*
//        Returns the equilibrium population fraction for all quantum-levels v, given temperature.
//    */
//        number summ = 0.0;
//        foreach(ej; 0 .. _n_vibe_states) { summ   += exp(-_vib_energy[ej]/(Boltzmann_constant*T)); }
//        foreach(ej; 0 .. _n_vibe_states) { nvf[ej] = exp(-_vib_energy[ej]/(Boltzmann_constant*T)) / summ }
//        return;
//    }
//} // end class EquilibriumVibSpecificCO

version(vib_specific_co_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        //lua_State* L = init_lua_State();
        //doLuaFile(L, );
        //auto gm = new VibSpecificCO(L);
        //lua_close(L);
        auto gm = new VibSpecificCO("sample-data/vib-specific-CO-gas.lua");
        auto Q = new GasState(gm.n_species, 0);

        // Practice problem to match the data in table II
        Q.p = 26.7;  // Pa
        Q.T = 175.0; // K
        // Set up the species mass fractions assuming equilibrium.
        foreach (v; 0 .. gm.n_vibe_states) { Q.massf[v] = gm.boltzmann_eq_population_fraction(v, Q.T); }

        double R_CO = 296.8379191791532; // gas constant for CO
        double M_CO = 0.0280101; // kg/mole
        double gamma = 7.0/5.0; // ratio of specific heats.

        gm.update_thermo_from_pT(Q);
        double my_rho = 26.7 / (R_CO * 175.0);
        assert(isClose(Q.rho, my_rho, 1.0e-6), failedUnitTest());

        double my_u = 2.5 * R_CO * 175.0;
        foreach (i; 0 .. gm.n_vibe_states) {
            my_u += (Avogadro_number/M_CO) * gm.vib_energy(i) * Q.massf[i];
        }
        assert(isClose(Q.u, my_u, 1.0e-6), failedUnitTest());

        //gm.update_trans_coeffs(Q);
        //assert(isClose(Q.mu, 0.0, 1.0e-6), failedUnitTest());
        //assert(isClose(Q.k, 0.0, 1.0e-6), failedUnitTest());

        gm.update_sound_speed(Q);
        double my_a = sqrt(gamma * R_CO * 175.0);
        assert(isClose(Q.a, my_a, 1.0e-6), failedUnitTest());

        return 0;
    }
}

