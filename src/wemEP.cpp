/* 
 * wemEG: Equation of State Generator
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
   bool useP = true;
   bool verbose = false;
   int num = 1;
   double step = 0.01;
   double startP = 1e6;
   double endP = 1e15;
   double startRho = 6220;
   double endRho = 1e5;
   double T = 0.0;
   string oFilename = "wemEP_out";
   EOS* eos = new EOS();
   eos->setNum(1);

   // accept user input
   //if (argc == 1)
   {
      printHelp(argv[0]);
      return 0;
   }
   /*int i = 1;
   while (i < argc)
   {
      if (!strcmp("-h", argv[i]) || !strcmp("-help", argv[i]))
      {
	 printHelp(argv[0]);
	 return 0;
      }
      else if (!strcmp("-n", argv[i]) || !strcmp("-num", argv[i]))
      {
	 if ((i == argc-1) || (!isNum(argv[i+1])))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 eos->setNum(atoi(argv[i+1]));
	 i += 2;
      }
      else if (!strcmp("-d", argv[i]) || !strcmp("-density", argv[i]))
      {
	 useP = false;
	 i++;
      }
      else if (!strcmp("-p", argv[i]) || !strcmp("-pressure", argv[i]))
      {
	 useP = true;
	 i++;
      }
      else if (!strcmp("-s", argv[i]) || !strcmp("-step", argv[i]))
      {
	 if ((i == argc-1) || (!isNum(argv[i+1])))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 step = atof(argv[i+1]);
	 i += 2;
      }
      else if (!strcmp("-b", argv[i]) || !strcmp("-start", argv[i]))
      {
	 if ((i == argc-1) || (!isNum(argv[i+1])))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 startP = startRho = atof(argv[i+1]);
	 i += 2;
      }
      else if (!strcmp("-e", argv[i]) || !strcmp("-end", argv[i]))
      {
	 if ((i == argc-1) || (!isNum(argv[i+1])))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 endP = endRho = atof(argv[i+1]);
	 i += 2;
      }
      else if (!strcmp("-tn", argv[i]) || !strcmp("-tempnum", argv[i]))
      {
	 if ((i == argc-1) || (!isNum(argv[i+1])))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 eos->setThermal(atoi(argv[i+1]));
	 i += 2;
      }
      else if (!strcmp("-t", argv[i]) || !strcmp("-temp", argv[i]))
      {
	 if ((i == argc-1) || (!isNum(argv[i+1])))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 T = atof(argv[i+1]);
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
      else if (!strcmp("-m", argv[i]) || !strcmp("-mix", argv[i]))
      {
	 if ((i == argc-1) || (!isNum(argv[i+1])))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 double frac = atof(argv[i+1]);
	 if ((frac < 0.0) || (frac > 1.0))
	 {
	    printHelp(argv[0]);
	    return 0;
	 }
	 eos->setMix(1, 2, frac);
	 i += 2;
      }
      else if (!strcmp("-v", argv[i]) || !strcmp("-verbose", argv[i]))
      {
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

   // calculate and print EOS
   if (useP)
      eos->printEOS(step, startP, endP, T, oFilename);
   else
      eos->printEOSRho(step, startRho, endRho, T, oFilename);//*/

   return 0;
}

void printHelp(const char* progName)
{
   printf("Usage: %s arguments\n", progName);
   printf("Arguments: -help/-h\n");
   printf("           -central/-c eos_file\n");
   printf("           -eos/-e eos_file mass_fraction\n");
   printf("           -mass/-m total_mass (Earth masses)\n");
   printf("           -radius/-r radius (Earth radii)\n");
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
