#pragma once

#include "stdafx.h"

// inside Planet.h
#ifndef PLANET_H
#define PLANET_H

// note: all values in SI
class Planet
{
private:
	// private constants
	static const double ITERATE_PRECISION; // precision when iterating
	static const double DIFF_PRECISION;    // precision when differentiating

	// private variables
	bool isRecording;   // is this recording its intermediate values?
	double h;           // radius step size
	double PCentral;    // central pressure

	// current (intermediate) values
	double r;     // radius
	double rho;   // density
	double rhoLast; // density of previous layer (for Love number)
	double rhoLast2;// density of 2nd previous layer (for Love number)
	double dRhods;  // derivative of density wrt radius (for Love number)
	double dVds;    // derivative of gravitational potential wrt radius (for Love number)
	double m;     // mass
	double mLast; // mass of previous layer (mass at P<0 is undefined)
	double T;     // temperature
	double Ug;    // gravitational potential energy
	double W;     // internal energy (P dV)
	double Ut;    // thermal energy
	double I;     // moment of inertia
	double Tn;    // function for Love number
	double X;     // dT/ds
	EOS* eos;     // equation of state

	// storage for values
	vector<double> rVals;   // radius
	vector<double> rhoVals; // density
	vector<double> mVals;   // mass
	vector<double> pVals;   // pressure

	// storage for eos boundaries (must be in order of mass)
	queue<double> mBoundaries;  // mass at which EOS changes
	queue<EOS*> eosBoundaries;   // new EOS

public:
	// public methods
	// setup methods
	Planet(double setH, double Pc, EOS* EOSc, double Tc);
	void addEOS(double newM, EOS* newEOS);  // add a mass boundary and eos
	void setRecord();  // set the planet to record intermediate values
	void integrate();  // integrate the planet (this may only be done if surface pressure > 0)

	void printRecord(string outFile); // print the recorded intermediate values to outFile

	// meta-methods
	void fixMass(double mass); // find and integrate from the central pressure which yields argument mass

	// main methods
	double getR();
	double getRho();
	double getM();
	double getT();
	double getPc();
	double getUg();
	double getW();
	double getUt();
	double getE();
	double getI();
	double getk2();

	// static methods
        // integrate a planet of given core mass fraction and mass
	static Planet* PlanetCMF(double setH, double CMF, double M);
        // PlanetCMF at given temperature
	static Planet* PlanetCMF(double setH, double CMF, double M, 
                                 double T); 
	static Planet* PlanetCMF(double setH, double CMF, double M, 
                                 double T, double guessP);
	/* integrate planets from mass Mi to Mf, with given logarithmic
           step, and output results to outFile */
	static void PrintCMF(double setH, double CMF, double Mi, 
                             double Mf, double step, double T, 
                             string outFile);
	static void Print(   double setH, int EOSnum, double Mi, 
                             double Mf, double step, string outFile);
	static void PrintTab(double setH, string inFile, double Mi, 
                             double Mf, double step, double T, 
                             string outFile);
	static void PrintTabLin(double setH, string inFile, double Mi, 
                                double Mf, double step, double T,
                                string outFile);
	static void Print2TabLin(double setH, string inFileC, 
                                 string inFileM, double CMF, double Mi,
                                 double Mf, double step, double T, 
                                 string outFile);
	static void PrintTabLinCSV(double setH, string inFile, 
                                   double Mi, double Mf, double step, 
                                   double T, string outFile);
	static void Print2TabLinCSV(double setH, string inFileC, 
                                    string inFileM, double CMF, 
                                    double Mi, double Mf, double step, 
                                    double T, string outFile);
private:
	// private methods
	// main methods
	void step();          // step the iteration forward once
	void record();        // record the current intermediate values
	void checkBoundary(); // check if the planet has reached an EOS boundary, change EOS (conserving P) if so
	void updateE();       // update energy values

	// static print methods
	static void printHeader(std::ofstream outputFile); // print header to output file
	static void printPlanet(std::ofstream outputFile, Planet* planet, EOS eosC, double T); // print current planet's stats to output file

	// meta-methods
	double dMdP(EOS* initialEOS, queue<double>* mBoundsFixed, queue<EOS*>* eosBoundsFixed);
	void clearIntegration();

	// runge-kutta methods
	double f(double rn, double rhon, double mn); // dRho/dR = g(r) * rho / (dP/dRho)
	double g(double rn, double rhon);            // dM/dR = area of shell of radius r times density
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
	double fT(double Xn);                       // dT/ds = X
	double gX(double sn, double Tn, double Xn); // dX/dx = -(2/s) * X - A(s) * T
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
