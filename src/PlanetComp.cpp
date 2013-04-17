/*
 * Represents a planet composition, or the class of planets
 * with that composition
 * Wesley Verne
 */

#include "stdafx.h"
#include "files.h"

using namespace cons;

// constants
const double PlanetComp::ITERATE_PRECISION = 0.00001;
const int    PlanetComp::ITERATE_CUTOFF    = 20;
const double PlanetComp::DIFF_PRECISION    = 0.001;
const double PlanetComp::INITIAL_P         = 1e10;

// --- public functions ---
// setup functionsa
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
// return a planet of the given mass with this composition
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
	    cerr << "Integration cutoff exceeded\n";
	    break;
	}
	i++;
    }
    if (verbose)
	cout << "\n";

    return createPlanet(pCentral, m);
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
    Planet planetM = createPlanet(Pc * (1 - DIFF_PRECISION), mGuess);
    Planet planetP = createPlanet(Pc * (1 + DIFF_PRECISION), mGuess);
    return (planetP.getMTotal() - planetM.getMTotal()) / 
	   (2 * Pc * DIFF_PRECISION);
}
