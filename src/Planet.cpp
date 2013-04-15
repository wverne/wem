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
}

void Planet::addEOS(double newM, EOS* newEOS)
{ boundaries.push(EOSBoundary(newEOS, newM)); }

void Planet::setVerbose(bool verboseSet) { verbose = verboseSet; }

void Planet::setT(double newT) { T = newT; }

/*void Planet::integrate()
{
    if (pVals.back() <= 0.0)
	throw "Cannot integrate complete planet";

    int i = 0;
    while (pVals.back() > 0.0)
    {
	if (boundaries.size() > 0)
	{
	    checkBoundary();
	}
	step();
	stepT();
	if (verbose && (i == 10000))
	{
	    cout << rhoVals.back() << endl;
	    i = 0;
	}
	i++;
    }
}//*/

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
void Planet::printRecord(string outFile)
{
    ofstream outputFile (outFile.c_str());

    outputFile << "Radius (Earth Radii) | Density (g cm^-3) | Pressure (Mbar) | Mass (Earth Masses)\n";

    for (int i = 0; i < getNumLayers(); i++)
	printf("%d %d %d %d", rVals[i] / R_EARTH, rhoVals[i] / 1e3,
	       pVals[i] / 1e11, mVals[i] / M_EARTH); 
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

double Planet::getUg(int layer)
{
    throw "Method not implemented";
}

double Planet::getW(int layer)
{
    throw "Method not implemented";
}

double Planet::getUt(int layer)
{
    throw "Method not implemented";
}

double Planet::getE(int layer)
{
    throw "Method not implemented";
}

double Planet::getI(int layer)
{
    throw "Method not implemented";
}

/*double Planet::getk2()
{
	return ((Tn * r) / (G * mLast)) - 1.0;
}//*/

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


/*
void Planet::Print(double setH, int EOSnum, double Mi, double Mf, double step, 
                   string outFile)
{
   std::ofstream outputFile (outFile.c_str());
	
	outputFile << "Mass (Earth Masses) | Radius (Earth Radii) | Central Pressure (Mbar) | Central Density (g cm^-3) | ";
	outputFile << "Gravitational Potential (J * 10^32) | Mechanical Energy (J * 10^32) | Thermal Energy (J * 10^32) | ";
	outputFile << "Total Energy (J * 10^32) | Moment of Inertia Coefficient | Love Number k2" << endl;

	EOS eos = EOS();
	eos.setNum(EOSnum);

	double M = Mi;
	double logM = log10(Mi);
	double guessP = eos.getP(10000, 0.0);
	
	Planet* planet = new Planet(setH, 1e9, &eos, 0.0);

	while (M < Mf)
	{
		M = pow(10, logM);

		planet->fixMass(M);
		
		outputFile << planet->getM()/M_EARTH << " " << planet->getR()/R_EARTH << " " << planet->getPc()/1e11 << " ";
		outputFile << eos.getRho(planet->getPc(), 0.0)/1e3 << " " << planet->getUg()/1e32 << " " << planet->getW()/1e32;
		outputFile << " " << planet->getUt()/1e32 << " " << planet->getE()/1e32 << " ";
		outputFile << planet->getI()/(planet->getM()*planet->getR()*planet->getR()) << " " << planet->getk2() << endl;

		guessP = planet->getPc();
		planet = new Planet(setH, guessP, &eos, 0.0);
		logM += step;
	}
}

void Planet::PrintTab(double setH, const string inFile, double Mi, 
		      double Mf, double step, double T, 
		      const string outFile)
{
   std::ofstream outputFile (outFile.c_str());
	
	outputFile << "Mass (Earth Masses) | Radius (Earth Radii) | Central Pressure (Mbar) | Central Density (g cm^-3) | ";
	outputFile << "Gravitational Potential (J * 10^32) | Mechanical Energy (J * 10^32) | Thermal Energy (J * 10^32) | ";
	outputFile << "Total Energy (J * 10^32) | Moment of Inertia Coefficient | Love Number k2" << endl;

	EOS* eosC = new EOS();
	eosC->setMixTab(inFile);

	double M = Mi;
	double logM = log10(Mi);
	double guessP = 1.1e11;

	while (M < Mf)
	{
		M = pow(10, logM);

		Planet* planet = new Planet(setH, guessP, eosC, T);
		planet->fixMass(M);
		
		outputFile << planet->getM()/M_EARTH << " " << planet->getR()/R_EARTH << " " << planet->getPc()/1e11 << " ";
		outputFile << eosC->getRho(planet->getPc(), T)/1e3 << " " << planet->getUg()/1e32 << " " << planet->getW()/1e32;
		outputFile << " " << planet->getUt()/1e32 << " " << planet->getE()/1e32 << " ";
		outputFile << planet->getI()/(planet->getM()*planet->getR()*planet->getR()) << " " << planet->getk2() << endl;

		guessP = planet->getPc();
		logM += step;
	}
}

void Planet::PrintTabLin(double setH, const string inFile, double Mi, 
			 double Mf, double step, double T, 
			 const string outFile)
{
   std::ofstream outputFile (outFile.c_str());
	
	outputFile << "Mass (Earth Masses) | Radius (Earth Radii) | Central Pressure (Mbar) | Central Density (g cm^-3) | ";
	outputFile << "Gravitational Potential (J * 10^32) | Mechanical Energy (J * 10^32) | Thermal Energy (J * 10^32) | ";
	outputFile << "Total Energy (J * 10^32) | Moment of Inertia Coefficient | Love Number k2" << endl;

	EOS* eosC = new EOS();
	eosC->setMixTab(inFile);

	double M = Mi;
	double guessP = 1.1e11;

	while ((M - step) < Mf)
	{
		Planet* planet = new Planet(setH, guessP, eosC, T);
		planet->fixMass(M);
		
		outputFile << planet->getM()/M_EARTH << " " << planet->getR()/R_EARTH << " " << planet->getPc()/1e11 << " ";
		outputFile << eosC->getRho(planet->getPc(), T)/1e3 << " " << planet->getUg()/1e32 << " " << planet->getW()/1e32;
		outputFile << " " << planet->getUt()/1e32 << " " << planet->getE()/1e32 << " ";
		outputFile << planet->getI()/(planet->getM()*planet->getR()*planet->getR()) << " " << planet->getk2() << endl;

		guessP = planet->getPc();
		M += step;
	}
}

void Planet::Print2TabLin(double setH, const string inFileC, 
                          const string inFileM, double CMF, double Mi, 
                          double Mf, double step, double T, 
			  const string outFile)
{
   std::ofstream outputFile (outFile.c_str());
	
	outputFile << "Mass (Earth Masses) | Radius (Earth Radii) | Central Pressure (Mbar) | Central Density (g cm^-3) | ";
	outputFile << "Gravitational Potential (J * 10^32) | Mechanical Energy (J * 10^32) | Thermal Energy (J * 10^32) | ";
	outputFile << "Total Energy (J * 10^32) | Moment of Inertia Coefficient | Love Number k2" << endl;

	EOS* eosC = new EOS();
	eosC->setMixTab(inFileC);

	EOS* eosM = new EOS();
	eosM->setMixTab(inFileM);

	double M = Mi;
	double guessP = 1e11;

	while ((M - step) < Mf)
	{
		Planet* planet = new Planet(setH, guessP, eosC, T);
		planet->addEOS(M * CMF, eosM);
		planet->fixMass(M);
		
		outputFile << planet->getM()/M_EARTH << " " << planet->getR()/R_EARTH << " " << planet->getPc()/1e11 << " ";
		outputFile << eosC->getRho(planet->getPc(), T)/1e3 << " " << planet->getUg()/1e32 << " " << planet->getW()/1e32;
		outputFile << " " << planet->getUt()/1e32 << " " << planet->getE()/1e32 << " ";
		outputFile << planet->getI()/(planet->getM()*planet->getR()*planet->getR()) << " " << planet->getk2() << endl;

		guessP = planet->getPc();
		M += step;
	}
}

void Planet::PrintTabLinCSV(double setH, string inFile, double Mi, double Mf, double step, double T, string outFile)
{
   std::ofstream outputFile (outFile.c_str());
	
	outputFile << "Mass (Earth Masses),Radius (Earth Radii),Central Pressure (Mbar),Central Density (g cm^-3),";
	outputFile << "Moment of Inertia Coefficient,Love Number" << endl;

	EOS* eosC = new EOS();
	eosC->setMixTab(inFile);

	double M = Mi;
	double guessP = 1.1e11;

	while ((M - step) < Mf)
	{
		Planet* planet = new Planet(setH, guessP, eosC, T);
		planet->fixMass(M);
		
		double Cnd = planet->getI()/(planet->getM()*planet->getR()*planet->getR());
		
		outputFile << planet->getM()/M_EARTH << "," << planet->getR()/R_EARTH << "," << planet->getPc()/1e11 << ",";
		outputFile << eosC->getRho(planet->getPc(), T)/1e3 << "," << Cnd << ",";
		outputFile << ((5/(((2.5*(1-(1.5*Cnd)))*(2.5*(1-(1.5*Cnd))))+1)) - 1) << endl;

		guessP = planet->getPc();
		M += step;
	}
}

void Planet::Print2TabLinCSV(double setH, string inFileC, string inFileM, double CMF, double Mi, double Mf, double step, double T, string outFile)
{
   std::ofstream outputFile (outFile.c_str());
	
	outputFile << "Mass (Earth Masses),Radius (Earth Radii),Central Pressure (Mbar),Central Density (g cm^-3),";
	outputFile << "Moment of Inertia Coefficient,Love Number" << endl;

	EOS* eosC = new EOS();
	eosC->setMixTab(inFileC);

	EOS* eosM = new EOS();
	eosM->setMixTab(inFileM);

	double M = Mi;
	double guessP = 1e11;

	while ((M - step) < Mf)
	{
		Planet* planet = new Planet(setH, guessP, eosC, T);
		planet->addEOS(M * CMF, eosM);
		planet->fixMass(M);

		double Cnd = planet->getI()/(planet->getM()*planet->getR()*planet->getR());
		
		outputFile << planet->getM()/M_EARTH << "," << planet->getR()/R_EARTH << "," << planet->getPc()/1e11 << ",";
		outputFile << eosC->getRho(planet->getPc(), T)/1e3 << "," << Cnd << ",";
		outputFile << ((5/(((2.5*(1-(1.5*Cnd)))*(2.5*(1-(1.5*Cnd))))+1)) - 1) << endl;

		guessP = planet->getPc();
		M += step;
	}
}

// private functions
// main functions
void Planet::step()
{
	mLast   = m;
	rhoLast = rho;
	double rho1 = kRho1();
	double m1   =   kM1();
	double rho2 = kRho2(rho1, m1);
	double m2   =   kM2(rho1, m1);
	double rho3 = kRho3(rho2, m2);
	double m3   =   kM3(rho2, m2);
	double rho4 = kRho4(rho3, m3);
	double m4   =   kM4(rho3, m3);

	rho += (h / 6) * (rho1 + (2 * rho2) + (2 * rho3) + rho4);
	m   += (h / 6) * (m1   + (2 * m2)   + (2 * m3)   + m4);
	r += h;

	dRhods = (rho - rhoLast) / h;
}

void Planet::record()
{
	rVals.push_back(r);
	rhoVals.push_back(rho);
	mVals.push_back(m);
	pVals.push_back(eos->getP(rho, T));
}

void Planet::checkBoundary()
{
	if (m > mBoundaries.front())
	{
		double boundaryP = eos->getP(rho, T);
		eos = eosBoundaries.front();
		rhoLast = rho;
		rho = eos->getRho(boundaryP, T);

		mBoundaries.pop();
		eosBoundaries.pop();

		X += ((4 * PI * G) / dVds) * (rhoLast - rho) * Tn;
	}
}

void Planet::updateE()
{
	dVds = (G / r) * ((m / r) - ((m - mLast) / h));
	Ug -= G * (m - mLast) * mLast / (r - (h / 2));
	W  += eos->getP(rho, T) * (4.0 / 3.0) * PI * ((r * r * r) - ((r - h) * (r - h) * (r - h)));
	if (eos->isTherm())
	{
		double molarMass;
		if (eos->getEosNum() == 1)
			molarMass = M_FE;
		else if (eos->getEosNum() == 2)
			molarMass = M_PV;
		else
			throw "Other materials not yet supported";
		Ut += eos->getMolETherm(T) * (m - mLast) / molarMass;
	}
	I += (2.0 / 3.0) * (m - mLast) * (r - (h / 2)) * (r - (h / 2));
}


// meta-functions
double Planet::dMdP(EOS* initialEOS, queue<double>* mBoundsFixed, queue<EOS*>* eosBoundsFixed)
{
	double PcCopy = PCentral;
	double mCopy  = m;
	clearIntegration();
	eos = initialEOS;
	mBoundaries   = *mBoundsFixed;
	eosBoundaries = *eosBoundsFixed;
	PCentral = PcCopy * (1 + DIFF_PRECISION);
	rho = eos->getRho(PCentral, T);
	integrate();
	PCentral = PcCopy;
	return (m - mCopy) / (PCentral * DIFF_PRECISION);
}

void Planet::clearIntegration()
{
	r     = 0.0;
	m     = 0.0;
	mLast = 0.0;
	Ug    = 0.0;
	W     = 0.0;
	Ut    = 0.0;
	I     = 0.0;

	if (isRecording)
	{
		isRecording = false;
		rVals.erase(rVals.begin(), rVals.end());
		rhoVals.erase(rhoVals.begin(), rhoVals.end());
		mVals.erase(mVals.begin(), mVals.end());
	}
}

// runge-kutta functions
double Planet::f(double rn, double rhon, double mn)
{
	if (rn > 0)
		return -((G * mn) / (rn * rn)) * rhon / eos->getdPdRho(rhon, T);
	else
		return 0;
}

double Planet::g(double rn, double rhon)
{
	return (4 * PI * rn * rn * rhon);
}

double Planet::kRho1()
{
	return f(r, rho, m);
}

double Planet::kM1()
{
	return g(r, rho);
}

double Planet::kRho2(double kRho, double kM)
{
	return f(r + (0.5 * h), rho + (0.5 * h * kRho), m + (0.5 * h * kM));
}

double Planet::kM2(double kRho, double kM)
{
	return g(r + (0.5 * h), rho + (0.5 * h * kRho));
}

double Planet::kRho3(double kRho, double kM)
{
	return f(r + (0.5 * h), rho + (0.5 * h * kRho), m + (0.5 * h * kM));
}

double Planet::kM3(double kRho, double kM)
{
	return g(r + (0.5 * h), rho + (0.5 * h * kRho));
}

double Planet::kRho4(double kRho, double kM)
{
	return f(r + h, rho + (h * kRho), m + (h * kM));
}

double Planet::kM4(double kRho, double kM)
{
	return g(r + h, rho + (h * kRho));
}

// Love number runge-kutta methods
void Planet::stepT()
{
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
}

double Planet::fT(double Xn)
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
