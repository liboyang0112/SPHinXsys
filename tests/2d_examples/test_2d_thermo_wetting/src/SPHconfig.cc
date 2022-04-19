/* ----------------------------------------------------------------------------
   libconfig - A library for processing structured configuration files
   Copyright (C) 2005-2010  Mark A Lindner

   This file is part of libconfig.

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   as published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, see
   <http://www.gnu.org/licenses/>.
   ----------------------------------------------------------------------------
*/

#include "SPHconfig.h"
#include  <memory>
using namespace std;
using namespace libconfig;

// This example reads the configuration file 'example.cfg' and displays
// some of its contents.

int readConfigFile(const char * filename, Config &cfg)
{

  // Read the file. If there is an error, report it and exit.
  try
  {
    cfg.readFile(filename);
  }
  catch(const FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}

SPHconfig::SPHconfig(){
  int ret = readConfigFile("fluid.cfg", cfg);
  if(ret==EXIT_FAILURE) exit(ret);

  // Get the store name.
  try
  {
    string name = cfg.lookup("name");
    cout << "config file name: " << name << endl << endl;
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No 'name' setting in configuration file." << endl;
  }

  const Setting& root = cfg.getRoot();

  // Output a list of all vdWFluids in the inventory.
  try
  {
    vdWFluids = &(root["Materials"]["vdWFluid"]);
    int count = vdWFluids->getLength();

    cout << "vdWFluids:"<<endl
         << setw(10) << left << "name" << "  "
         << setw(10) << left << "visc" << "   "
         << setw(10) << left << "Kappa" << "  "
         << setw(10) << left << "T0" << "  "
         << "rho_0"
         << endl;

    for(int i = 0; i < count; ++i)
    {
      const Setting &vdWFluid = (*vdWFluids)[i];

      // Only output the record if all of the expected fields are present.
      string name;
      double viscousity, heatconduction, temperature, rho_0;            //initial density

      if(!(vdWFluid.lookupValue("name", name)
           && vdWFluid.lookupValue("viscousity", viscousity)
           && vdWFluid.lookupValue("heatconduction", heatconduction)
           && vdWFluid.lookupValue("temperature", temperature)
           && vdWFluid.lookupValue("rho_0", rho_0)))
        continue;

    cout << setw(10) << left << name << "  "
         << setw(10) << left << viscousity << "   "
         << setw(10) << left << heatconduction << "  "
         << setw(10) << left << temperature << "  "
         << rho_0
         << endl;
    }
    cout << endl;
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No vdWFluids setting in configuration file." << endl;
  }

  // Output a list of all vdWFluids in the inventory.
  try
  {
    Solids = &(root["Materials"]["Solid"]);
    int count = Solids->getLength();

    cout << "Solids:"<<endl
         << setw(10) << left << "name" << "  "
         << setw(10) << left << "heat cond" << "  "
         << setw(10) << left << "Temperature"
         << endl;

    for(int i = 0; i < count; ++i)
    {
      const Setting &Solid = (*Solids)[i];

      // Only output the record if all of the expected fields are present.
      string name;
      double heatconduction, temperature;            //initial density

      if(!(Solid.lookupValue("name", name)
           && Solid.lookupValue("heatconduction", heatconduction)
           && Solid.lookupValue("temperature", temperature)))
        continue;

      cout << setw(10) << left << name << "  "
         << setw(10) << left << heatconduction << "  "
         << setw(10) << left << temperature
         << endl;
    }
    cout << endl;
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No Solids setting in configuration file." << endl;
  }
  try{
    Job = &(root["Job"]);
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No Job setting in configuration file." << endl;
  }
  try{
    ExternalForce = &(root["ExternalForce"]);
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No ExternalForce setting in configuration file." << endl;
  }
}

// eof