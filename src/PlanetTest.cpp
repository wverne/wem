/* 
 * Tests Planet class
 * Wesley Verne
 */

#include "stdafx.h"
#include "files.h"

int main(int argc, const char* argv[])
{
  EOS *eosC = new EOS();
  eosC->setNum(0);
  EOS *eos1 = new EOS();
  eos1->setNum(1);
  EOS *eos2 = new EOS();
  eos2->setNum(2);
  EOS *eos3 = new EOS();
  eos3->setNum(3);
  
  Planet planet1(1000.0, 1e30, eosC);
  planet1.addEOS(100.0, eos1);
  planet1.addEOS(300.0, eos3);
  planet1.addEOS(200.0, eos2);

  planet1.printBoundaries();

  Planet planet2(1000.0, 1e11, eos1);
  planet2.integrate();
  planet2.printRecord("testplanet1", 100);

  Planet planet3(10.0, 1e11, eos1);
  planet3.integrate();
  planet3.printRecord("testplanet2", 10000);
}
