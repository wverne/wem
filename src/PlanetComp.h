/*
 * Represents a planet composition, or the class of planets 
 * with that composition
 * Wesley Verne
 */

#include "stdafx.h"
#include "EOS.h"
#include "Planet.h"

// inside PlanetComp.h
#ifndef PLANET_COMP_H
#define PLANET_COMP_H

using namespace cons;

// note: all values in SI
class PlanetComp
{
private:
    // private constants
    static const double ITERATE_PRECISION; // precision when iterating
    static const int    ITERATE_CUTOFF;    // max number of iterations
    static const double DIFF_PRECISION;    // precision when differentiating
    static const double INITIAL_P;         // initial pressure for fixMass

    // private variables
    double h;     // radius step size
    EOS *eosc;    // central EOS
    bool verbose; // verbose output
    double T;     // temperature

    // storage for eos boundaries
    // the give ndouble is the mass fraction at which the EOS changes
    typedef pair<EOS*, double> EOSBoundaryFrac;
    vector<EOSBoundaryFrac> boundaries;

public:
    // --- public methods ---
    // --- setup methods ---
    PlanetComp(double setH, EOS* EOSc);
    ~PlanetComp();
    // reset central EOS
    void setEOSc(EOS* newEOSc);
    // reset step size
    void setH(double newH);
    // add a mass fraction boundary and eos
    void addEOS(double newM, EOS* newEOS);  
    // set verbose output
    void setVerbose(bool verboseSet);
    // set temperature
    void setT(double newT);

    // --- main methods ---
    // return a planet of the given mass with this composition
    Planet fixMass(double mass); 
    // return a planet of the given radius with this composition
    // NOT IMPLEMENTED: Difficult because EOS boundaries are defined
    // on mass
    Planet fixRadius(double radius); 

private:
    // --- private methods ---
    // initializes and integrates a planet with the given central pressure
    // with this composition, using mGuess to create boundaries
    Planet createPlanet(double Pc, double mGuess = M_EARTH);

    // determines dM/dP for a planet with the given central pressure
    // with this composition, using mGuess to create boundaries
    double dMdP(double Pc, double mGuess = M_EARTH);

    // determines dR/dP for a planet with the given central pressure
    // with this composition, using mGuess to create boundaries
    double dRdP(double Pc, double mGuess = M_EARTH);
};

#endif
