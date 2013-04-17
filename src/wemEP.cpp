/* 
 * wemEP: Exoplanet Profiler
 * Wesley Verne
 */

#include "stdafx.h"
#include "files.h"

using namespace cons;
using namespace std;

void printHelp(const char* progName);
bool isNum(const char* arg);

int main(int argc, const char* argv[])
{
    // default values
    EOS* eosc = new EOS();
    eosc->setNum(1);
    int interval = 100;
    double mass = M_EARTH;
    string oFilename = "wemEP_out";
    bool verbose = false;

    PlanetComp planetComp = PlanetComp(100, eosc);

    // accept user input
    if (argc == 1)
    {
	printHelp(argv[0]);
	return 0;
    }

    int i = 1;
    while (i < argc)
    {
	if (!strcmp("-h", argv[i]) || !strcmp("-help", argv[i]))
	{
	    printHelp(argv[0]);
	    return 0;
	}
	else if (!strcmp("-c", argv[i]) || !strcmp("-central", argv[i]))
	{
	    if (i == argc-1)
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    eosc->setMixTab(argv[i+1]);
	    i += 2;
	}
	else if (!strcmp("-cn", argv[i]) || !strcmp("-centralnum", argv[i]))
	{
	    if ((i == argc-1) || (!isNum(argv[i+1])))
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    eosc->setNum(atoi(argv[i+1]));
	    i += 2;
	}
	else if (!strcmp("-e", argv[i]) || !strcmp("-eos", argv[i]))
	{
	    if ((i >= argc-2) || (!isNum(argv[i+2])))
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    EOS* eosn = new EOS();
	    eosn->setMixTab(argv[i+1]);
	    planetComp.addEOS(atof(argv[i+2]), eosn);
	    i += 3;
	}
	else if (!strcmp("-en", argv[i]) || !strcmp("-eosnum", argv[i]))
	{
	    if ((i >= argc-2) || (!isNum(argv[i+1])) || (!isNum(argv[i+2])))
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    EOS* eosn = new EOS();
	    eosn->setNum(atoi(argv[i+1]));
	    planetComp.addEOS(atof(argv[i+2]), eosn);
	    i += 3;
	}
	else if (!strcmp("-i", argv[i]) || !strcmp("-interval", argv[i]))
	{
	    if ((i == argc-1) || (!isNum(argv[i+1])))
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    interval = atoi(argv[i+1]);
	    i += 2;
	}
	else if (!strcmp("-m", argv[i]) || !strcmp("-mass", argv[i]))
	{
	    if ((i == argc-1) || (!isNum(argv[i+1])))
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    mass = atof(argv[i+1]) * M_EARTH;
	    i += 2;
	}
	else if (!strcmp("-s", argv[i]) || !strcmp("-step", argv[i]))
	{
	    if ((i == argc-1) || (!isNum(argv[i+1])))
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    planetComp.setH(atof(argv[i+1]));
	    i += 2;
	}
	else if (!strcmp("-o", argv[i]) || !strcmp("-outfile", argv[i]))
	{
	    if (i == argc-1)
	    {
		printHelp(argv[0]);
		return 0;
	    }
	    oFilename = argv[i+1];
	    i += 2;
	}
	else if (!strcmp("-v", argv[i]) || !strcmp("-verbose", argv[i]))
	{
	    planetComp.setVerbose(true);
	    verbose = true;
	    i++;
	}
	else
	{
	    fprintf(stderr, "Unrecognized argument: %s\n", argv[i]);
	    i++;
	}
    }

    if (verbose)
	printf("Calculating planet profile...\n");

    // calculate and output planet profile
    Planet planet = planetComp.fixMass(mass);
    planet.printRecord(oFilename, interval);

    delete eosc;
    return 0;
}

void printHelp(const char* progName)
{
    printf("Usage: %s arguments\n", progName);
    printf("Arguments: -help/-h\n");
    printf("           -central/-c eos_file\n");
    printf("           -centralnum/-cn eos_num\n");
    printf("           -eos/-e eos_file mass_fraction\n");
    printf("           -eosnum/-en eos_num mass_fraction\n");
    printf("           -interval/-i output_interval\n");
    printf("           -mass/-m total_mass (Earth masses)\n");
    printf("           -step/-s step_size (meters)\n");
    printf("           -outfile/-o filename\n");
    printf("           -verbose/-v\n");
    printf("See readme for details and default values\n");
}

bool isNum(const char* arg)
{
    if (*arg == '-')
    {
	if (!isdigit(arg[1]) && !(arg[1] == '.'))
	    return false;
    }
    return true;
}
