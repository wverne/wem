/*
 * Represents a planet composition, or the class of planets
 * with that composition
 * Wesley Verne
 */

#include "stdafx.h"
#include "files.h"

using namespace cons;

// constants
const double PlanetComp::ITERATE_PRECISION = 0.0001;
const int    PlanetComp::ITERATE_CUTOFF    = 10;
const double PlanetComp::DIFF_PRECISION    = 0.001;
const double PlanetComp::INITIAL_P         = 1e10;

// --- public functions ---
// setup functions
PlanetComp::PlanetComp(double setH, EOS* EOSc)
{
    h = setH;
    eosc = EOSc;
    verbose = false;
    T = 0.0;
}

PlanetComp::~PlanetComp()
{
    for (int i = 0; i < boundaries.size(); i++)
	delete boundaries[i].first;
}

void PlanetComp::setEOSc(EOS* newEOSc) { eosc = newEOSc; }

void PlanetComp::setH(double newH) { h = newH; }

void PlanetComp::addEOS(double newM, EOS* newEOS)
{ boundaries.push_back(EOSBoundaryFrac(newEOS, newM)); }

void PlanetComp::setVerbose(bool verboseSet) { verbose = verboseSet; }

void PlanetComp::setT(double newT) { T = newT; }

// --- main methods ---
Planet PlanetComp::fixMass(double mass)
{
    if (mass <= 0)
	throw "Invalid mass";

    double pCentral = INITIAL_P;
    int i = 0;
    
    Planet cPlanet = createPlanet(pCentral);

    double m = cPlanet.getMTotal();

    while ((abs(m - mass) / mass) > ITERATE_PRECISION)
    {
	if (verbose)
	    cout << m/M_EARTH << "\n";

	double pCentralLast = pCentral;
	pCentral -= (m - mass) / dMdP(pCentral, mass);
	if (pCentral <= 0) // protect against negative pressures
	    pCentral = pCentralLast / 2.0;

	m = createPlanet(pCentral, mass).getMTotal();

	if (i > ITERATE_CUTOFF)
	{
	    cerr << "Iteration cutoff exceeded\n";
	    break;
	}
	i++;
    }
    if (verbose)
	cout << m/M_EARTH << "\n";

    return createPlanet(pCentral, m);
}

Planet PlanetComp::fixRadius(double radius)
{
    throw "Method not implemented";
    // difficult because EOS boundaries are defined on mass
}

void PlanetComp::printMR(double startMass, double endMass, double step,
			 string outFile)
{
    ofstream outputFile (outFile.c_str());
    outputFile << "Mass (Earth Masses) | Radius (Earth Radii) | ";
    outputFile << "Central Pressure (Mbar) | Central Density (g cm^-3) | ";
    outputFile << "Gravitational Potential (J) | Mechanical Energy (J) | ";
    outputFile << "Moment of Inertia Coefficient\n";

    for (double mass = startMass, mass <= endMass, mass += step)
    {
	Planet planet = fixMass(mass);
	outputFile << (planet.getMTotal() / M_EARTH) << " ";
	outputFile << (planet.getRTotal() / R_EARTH) << " ";
	outputFile << (planet.getPc() / 1e11) << " ";
	outputFile << (planet.getRho(0) / 1e3) << " ";
	outputFile << planet.getUg() << " ";
	outputFile << planet.getW() << " ";
	outputFile << planet.getICoeff() << "\n";
    }
}

// --- private methods ---
Planet PlanetComp::createPlanet(double Pc, double mGuess)
{
    if ((Pc <= 0.0) || (mGuess <= 0.0))
	throw "Bad parameters";

    Planet pOut = Planet(h, Pc, eosc);
    pOut.setT(T);
    
    for (int i = 0; i < boundaries.size(); i++)
	pOut.addEOS(boundaries[i].second * mGuess, boundaries[i].first);

    pOut.integrate();

    return pOut;
}

double PlanetComp::dMdP(double Pc, double mGuess)
{
    Planet planetLess = createPlanet(Pc * (1 - DIFF_PRECISION), mGuess);
    Planet planetMore = createPlanet(Pc * (1 + DIFF_PRECISION), mGuess);
    return (planetMore.getMTotal() - planetLess.getMTotal()) / 
	   (2 * Pc * DIFF_PRECISION);
}

double PlanetComp::dRdP(double Pc, double mGuess)
{
    Planet planetLess = createPlanet(Pc * (1 - DIFF_PRECISION), mGuess);
    Planet planetMore = createPlanet(Pc * (1 + DIFF_PRECISION), mGuess);
    return (planetMore.getRTotal() - planetLess.getRTotal()) / 
	   (2 * Pc * DIFF_PRECISION);
}
