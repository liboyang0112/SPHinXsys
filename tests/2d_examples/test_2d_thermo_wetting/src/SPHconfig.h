/**
 * @file 	SPHconfig.h
 * @brief A config file header.
 * @author Boyang Li
 */
/**
 * @brief 	SPHinXsys Library.
 */

#ifndef __SPHCONFIG
#define __SPHCONFIG
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>

// This example reads the configuration file 'example.cfg' and displays
// some of its contents.


class SPHconfig
{
public:
  libconfig::Config cfg;
  libconfig::Setting* vdWFluids;
  libconfig::Setting* Solids;
  libconfig::Setting* Job;
  libconfig::Setting* ExternalForce;
  SPHconfig();
  ~SPHconfig(){};
};
#endif