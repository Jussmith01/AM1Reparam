#ifndef am1paramhandler_header
#define am1paramhandler_header

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
#include "../structs/AM1params_struct.h"
#include "../structs/coords_struct.h"
#include "../utils/tools.h"

using namespace std;

//      *************************************************************     //
//                       AM1 Parameter Handler Class			  //
//      Handles all AM1 Parameters loading, storage, and manipulation	  //
//      *************************************************************     //
class AM1ParameterHandler
{
	//-----------------------------
        //	    Private
        //-----------------------------
	/*Struct used to store the AM1 Parameters and ipt coords for each atom type*/
	vector<string> AM1GauFilePostfix;
	vector<AM1ParamsStore> APS;
	
	public:
        //-----------------------------
        //	    Public
        //-----------------------------
	vector<string> AM1GauFilePrefix;
	vector<coordstore> CS;

	/*Function loads AM1 parameters from a g09 input file into APS*/
	void LoadParametersFromG09Input(string infile);

	/*Function produces the initial parameters vector*/
	vector<float> ProduceInitialParamsVector(void);

	/*Function Creates a Gaussian input file with a set of modified params*/
	void ProduceGaussianInput(vector<float> &DNA,string fileprefix,int geom);
};
#endif

