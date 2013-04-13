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
  
  Planet planet(1000.0, 1e30, eosC);
  planet.addEOS(100.0, eos1);
  planet.addEOS(300.0, eos3);
  planet.addEOS(200.0, eos2);

  planet.printBoundaries();
}
