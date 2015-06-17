#ifndef readinput_header
#define readinput_header

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <time.h>
#include "../core/memhandler_class.h"
#include "tools.h"

using namespace std;

/*----------------------------------------
       ****************************
       |  Store Parameters in MH  |
       ****************************
------------------------------------------*/
void SaveData(string flag, string value, MemoryHandler &MH)
{
  if (flag.compare("TargetGauParams") == 0)
  {
    MH.IP.TargetGauParams = value;
  }
  else if (flag.compare("WorkingDir") == 0)
  {
    MH.IP.WkDir = value;
  }
  else if (flag.compare("AM1GauFilename") == 0)
  {
    MH.IP.AM1PFile = value;
  }
  else if (flag.compare("J") == 0)
  {
    MH.IP.J = atoi(value.c_str());
  }
  else if (flag.compare("N") == 0)
  {
    MH.IP.N = atoi(value.c_str());
  }
  else if (flag.compare("nproc") == 0)
  {
    MH.IP.nproc = atoi(value.c_str());
  }
  else if (flag.compare("cew") == 0)
  {
    MH.IP.cew = atof(value.c_str());
  }
  else if (flag.compare("ew") == 0)
  {
    MH.IP.ew = atof(value.c_str());
  }
  else if (flag.compare("few") == 0)
  {
    MH.IP.few = atof(value.c_str());
  }
  else if (flag.compare("a") == 0)
  {
    MH.IP.a = atof(value.c_str());
  }
  else if (flag.compare("b") == 0)
  {
    MH.IP.b = atof(value.c_str());
  }
  else if (flag.compare("forcebest") == 0)
  {
    MH.IP.forcebest = atoi(value.c_str());
  }
  else if (flag.compare("InitialMutations") == 0)
  {
    MH.IP.InitialMutations = atof(value.c_str());
  }
  else if (flag.compare("InitialPrecision") == 0)
  {
    MH.IP.InitialPrecision = atof(value.c_str());
  }
  else if (flag.compare("MainloopPrecision") == 0)
  {
    MH.IP.MainloopPrecision = atof(value.c_str());
  }
  else if (flag.compare("G") == 0)
  {
    MH.IP.G = atoi(value.c_str());
  }
  else if (flag.compare("Gstd") == 0)
  {
    MH.IP.Gstd = atof(value.c_str());
  }
  else if (flag.compare("conv") == 0)
  {
    MH.IP.conv = atof(value.c_str());
  }
  else
  {
    cout << "Unknown flag found input file. Flag: " << flag
         << " Value: " << value << endl;
  }
};

/*----------------------------------------
       ****************************
       |        Parse Line        |
       ****************************
Parse and store the parameters in the memory
handler.
------------------------------------------*/
void Parser(string line, MemoryHandler &MH)
{
  if (line.find("#") == string::npos)
  {
    size_t d1 = line.find(":");
    size_t d2 = line.find("!");

    string flag = trim(line.substr(0, d1));
    string value = trim(line.substr(d1 + 1, d2 - d1 - 1));

    SaveData(flag, value, MH);
  }
};

/*----------------------------------------
       ***************************
       |   Read Input Function   |
       ***************************
Function reads the inputs and obtains
program parameters.
------------------------------------------*/
extern void ReadInput(string filename, MemoryHandler &MH, outputFile &out)
{
  string line;
  ifstream inputfile(filename.c_str());
  if (inputfile.is_open())
  {
    while (getline(inputfile, line))
    {
      Parser(line, MH);
    }
    inputfile.close();

    out.PrintIptParms(MH);
  }
  else
    cout << "Unable to open file: " << filename.c_str() << endl;
}

#endif
