/*
 * Represents a Planet
 * Wesley Verne
 */

#include "stdafx.h"
#include "files.h"

using namespace cons;

// private constants
const double Planet::ITERATE_PRECISION = 0.00001;
const double Planet::DIFF_PRECISION    = 0.001;

// boundary comparator setup
bool Planet::EOSBoundaryCompare::operator()(const EOSBoundary& p1, 
					    const EOSBoundary& p2) const
{ return p1.second > p2.second; }

// public functions
// setup functions
Planet::Planet(double setH, double Pc, EOS* EOSc)
{
    h   = setH;
    eos = EOSc;
    rVals.push_back(0.0);
    rhoVals.push_back(eos->getRho(Pc));
    mVals.push_back(0.0);
    pVals.push_back(Pc);
    verbose = false;
    T = 0.0;
}

void Planet::addEOS(double newM, EOS* newEOS)
{ boundaries.push(EOSBoundary(newEOS, newM)); }

void Planet::setVerbose(bool verboseSet) { verbose = verboseSet; }

void Planet::setT(double newT) { T = newT; }

void Planet::integrate()
{
    if (pVals.back() <= 0.0)
	throw "Cannot integrate complete planet";

    int i = 0;
    while (pVals.back() > 0.0)
    {
	checkBoundary();
	step();
	stepT();
	if (verbose && (i == 10000))
	{
	    cout << rhoVals.back() << endl;
	    i = 0;
	}
	i++;
    }
}

// meta-functions
/*void Planet::fixMass(double mass)
{
	if (eos->getP(rho, T) < 0.0)
		throw "Cannot fix mass of complete planet";
	EOS* initialEOS = eos;
	queue<double> mBoundsFixed = mBoundaries;
	queue<EOS*> eosBoundsFixed = eosBoundaries;

	integrate();

	while ((abs(m - mass) / mass) > ITERATE_PRECISION)
	{
		cout << m/M_EARTH << endl;
		PCentral -= ((m - mass) / dMdP(initialEOS, &mBoundsFixed, 
					       &eosBoundsFixed));
		clearIntegration();
		eos = initialEOS;
		rho = eos->getRho(PCentral, T);
		mBoundaries   = mBoundsFixed;
		eosBoundaries = eosBoundsFixed;
		setRecord();
		integrate();
	}
	cout << endl;
}//*/

// main functions
void Planet::printRecord(string outFile, int interval)
{
    ofstream outputFile (outFile.c_str());

    outputFile << "Radius (Earth Radii) | Density (g cm^-3) | ";
    outputFile << "Pressure (Mbar) | Mass (Earth Masses)\n";

    for (int i = 0; i < (getNumLayers() - 1); i += interval)
    {
	outputFile << (rVals[i] / R_EARTH) << " ";
	outputFile << (rhoVals[i] / 1e3) << " ";
	outputFile << (pVals[i] / 1e11) << " ";
	outputFile << (mVals[i] / M_EARTH) << "\n";
    }
    outputFile << (rVals.back() / R_EARTH) << " ";
    outputFile << (rhoVals.back() / 1e3) << " ";
    outputFile << (pVals.back() / 1e11) << " ";
    outputFile << (mVals.back() / M_EARTH) << "\n";
}

int Planet::getNumLayers() { return rVals.size(); }

double Planet::getRTotal() { return rVals.back(); }

double Planet::getRhoTotal()
{
    return (3.0 * mVals.back()) / (4.0 * PI * pow(rVals.back(), 3));
}

double Planet::getMTotal() { return mVals.back(); }

double Planet::getPSurface() { return pVals.back(); }

double Planet::getPc() { return pVals.front(); }

double Planet::getT() { return T; }

double Planet::getR(int layer)
{
    if ((layer >= 0) && (layer < getNumLayers()))
	return rVals[layer];
    else
	throw "Not a vaild layer";
}

double Planet::getRho(int layer)
{
    if ((layer >= 0) && (layer < getNumLayers()))
	return rhoVals[layer];
    else
	throw "Not a vaild layer";
}

double Planet::getM(int layer)
{
    if ((layer >= 0) && (layer < getNumLayers()))
	return mVals[layer];
    else
	throw "Not a vaild layer";
}

double Planet::getP(int layer)
{
    if ((layer >= 0) && (layer < getNumLayers()))
	return pVals[layer];
    else
	throw "Not a vaild layer";
}

double Planet::getUg()
{
    throw "Method not implemented";
}

double Planet::getW()
{
    throw "Method not implemented";
}

double Planet::getUt()
{
    throw "Method not implemented";
}

double Planet::getE()
{
    throw "Method not implemented";
}

double Planet::getI()
{
    throw "Method not implemented";
}

double Planet::getk2()
{
    return ((getT() * getRTotal()) / (G * getMTotal())) - 1.0;
}

// debug functions
void Planet::printBoundaries()
{
    while (!boundaries.empty())
    {
	printf("boundary: %d, %f\n", boundaries.top().first->getEosNum(),
	       boundaries.top().second);
	boundaries.pop();
    }
}

// private functions
// main functions
void Planet::step()
{
	double rho1 = kRho1();
	double m1   =   kM1();
	double rho2 = kRho2(rho1, m1);
	double m2   =   kM2(rho1, m1);
	double rho3 = kRho3(rho2, m2);
	double m3   =   kM3(rho2, m2);
	double rho4 = kRho4(rho3, m3);
	double m4   =   kM4(rho3, m3);

	rhoVals.push_back(rhoVals.back() + 
			  ((h / 6) * 
			   (rho1 + (2 * rho2) + (2 * rho3) + rho4)));
	mVals.push_back(mVals.back() + 
			((h / 6) * 
			 (m1   + (2 * m2) + (2 * m3) + m4)));
	rVals.push_back(rVals.back() + h);
	pVals.push_back(eos->getP(rhoVals.back(), getT()));
}

void Planet::checkBoundary()
{
    if (!boundaries.empty() && (mVals.back() > boundaries.top().second))
    {
	eos = boundaries.top().first;
	rhoVals[rhoVals.size() - 1] = eos->getRho(pVals.back(), T);
	boundaries.pop();
    }
}

// runge-kutta functions
double Planet::f(double rn, double rhon, double mn)
{
    if (rn > 0)
	return -((G * mn) / (rn * rn)) * rhon / eos->getdPdRho(rhon, T);
    else
	return 0.0;
}

double Planet::g(double rn, double rhon)
{
    return (4 * PI * rn * rn * rhon);
}

double Planet::kRho1()
{
    return f(rVals.back(), rhoVals.back(), mVals.back());
}

double Planet::kM1()
{
    return g(rVals.back(), rhoVals.back());
}

double Planet::kRho2(double kRho, double kM)
{
    return f(rVals.back()   + (0.5 * h), 
	     rhoVals.back() + (0.5 * h * kRho), 
	     mVals.back()   + (0.5 * h * kM));
}

double Planet::kM2(double kRho, double kM)
{
    return g(rVals.back()   + (0.5 * h), 
	     rhoVals.back() + (0.5 * h * kRho));
}

double Planet::kRho3(double kRho, double kM)
{
	return f(rVals.back()   + (0.5 * h), 
		 rhoVals.back() + (0.5 * h * kRho), 
		 mVals.back()   + (0.5 * h * kM));
}

double Planet::kM3(double kRho, double kM)
{
	return g(rVals.back()   + (0.5 * h), 
		 rhoVals.back() + (0.5 * h * kRho));
}

double Planet::kRho4(double kRho, double kM)
{
	return f(rVals.back()   + h, 
		 rhoVals.back() + (h * kRho), 
		 mVals.back()   + (h * kM));
}

double Planet::kM4(double kRho, double kM)
{
	return g(rVals.back()   + h, 
		 rhoVals.back() + (h * kRho));
}

// Love number runge-kutta methods
void Planet::stepT()
{/*
    if (r < 2*h)
	return;
    double T1 = kT1();
    double X1 = kX1();
    double T2 = kT2(X1);
    double X2 = kX2(T1, X1);
    double T3 = kT3(X2);
    double X3 = kX3(T2, X2);
    double T4 = kT4(X3);
    double X4 = kX4(T3, X3);

    Tn += (h / 6) * (T1 + (2 * T2) + (2 * T3) + T4);
    X  += (h / 6) * (X1 + (2 * X2) + (2 * X3) + X4);
//*/
}

/*double Planet::fT(double Xn)
{
	return Xn;
}

double Planet::gX(double sn, double Tn, double Xn) // dX/dx = -(2/s) * X - A(s) * T
{
	return -((2 / sn) * X) - (A(sn) * Tn);
}

double Planet::A(double sn)
{
	return ((4 * PI * G * dRhods) / dVds) - (6 / (sn * sn));
}

double Planet::kT1()
{
	return fT(X);
}

double Planet::kX1()
{
	return gX(r - h, Tn, X);
}

double Planet::kT2(double kX)
{
	return fT(X + (0.5 * h * kX));
}

double Planet::kX2(double kT, double kX)
{
	return gX(r - (0.5 * h), Tn + (0.5 * h * kT), X + (0.5 * h * kX));
}

double Planet::kT3(double kX)
{
	return fT(X + (0.5 * h * kX));
}

double Planet::kX3(double kT, double kX)
{
	return gX(r - (0.5 * h), Tn + (0.5 * h * kT), X + (0.5 * h * kX));
}

double Planet::kT4(double kX)
{
	return fT(X + (h * kX));
}

double Planet::kX4(double kT, double kX)
{
	return gX(r, Tn + (h * kT), X + (h * kX));
}
//*/
