/*
 * Represents a Planet
 * Wesley Verne
 */

#include "stdafx.h"
#include "EOS.h"

// inside Planet.h
#ifndef PLANET_H
#define PLANET_H

// note: all values in SI
class Planet
{
private:
    // private variables
    double h;     // radius step size
    EOS *eos;     // current EOS
    bool verbose; // verbose output
    double T;     // temperature
	
    // storage for values
    vector<double> rVals;   // radius
    vector<double> rhoVals; // density
    vector<double> mVals;   // mass
    vector<double> pVals;   // pressure

    // storage for eos boundaries
    // the given double is the mass at which the EOS changes
    typedef pair<EOS*, double> EOSBoundary;
    struct EOSBoundaryCompare {
	bool operator()(const EOSBoundary& p1, const EOSBoundary& p2) const;
    };
    priority_queue<EOSBoundary, vector<EOSBoundary>, 
	EOSBoundaryCompare> boundaries; 

public:

    // --- public methods ---
    // --- setup methods ---
    Planet(double setH, double Pc, EOS* EOSc);
    // add a mass boundary and eos
    void addEOS(double newM, EOS* newEOS);  
    // set verbose output
    void setVerbose(bool verboseSet);
    // set temperature
    void setT(double newT);
    // integrate the planet (this may only be done if surface pressure > 0)
    void integrate();  

    // --- main methods ---
    // print the planet's profile to outFile
    void printRecord(string outFile, int interval = 1);
    int getNumLayers();
    double getRTotal();
    double getRhoTotal();
    double getMTotal();
    double getPSurface();
    double getPc();
    double getT();
    double getUg(); // gravitational potential
    double getW();
    double getUt();
    double getE();
    double getI();
    double getk2();

    double getR(  int layer);
    double getRho(int layer);
    double getM(  int layer);
    double getP(  int layer);

    // --- debug methods ---
    // prints the EOS boundaries in order (destroys PQ)
    void printBoundaries();

private:
    // --- private methods ---
    // --- main methods ---
    // step the iteration forward once
    void step();
    // check if EOS boundary reached, change EOS (conserving P) if so
    void checkBoundary(); 

    // --- runge-kutta methods ---
    // dRho/dR = g(r) * rho / (dP/dRho)
    double f(double rn, double rhon, double mn); 
    // dM/dR = area of shell of radius r times density
    double g(double rn, double rhon);
    double kRho1();
    double kM1();
    double kRho2(double kRho, double kM);
    double kM2(double kRho, double kM);
    double kRho3(double kRho, double kM);
    double kM3(double kRho, double kM);
    double kRho4(double kRho, double kM);
    double kM4(double kRho, double kM);

    // Love number runge-kutta methods
    void stepT();
    // dT/ds = X
    double fT(double Xn);
    // dX/dx = -(2/s) * X - A(s) * T
    double gX(double sn, double Tn, double Xn); 
    double A(double sn);
    double kT1();
    double kX1();
    double kT2(double kX);
    double kX2(double kT, double kX);
    double kT3(double kX);
    double kX3(double kT, double kX);
    double kT4(double kX);
    double kX4(double kT, double kX);
};


#endif // PLANET_H
