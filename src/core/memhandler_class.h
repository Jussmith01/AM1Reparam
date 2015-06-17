#ifndef memhandler_header
#define memhandler_header

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <limits.h>
#include <bitset>
#include <time.h>
#include "../structs/iptparams_struct.h"
#include "../utils/flaghandler_class.h"
#include "AM1Parameters_class.h"

using namespace std;

//      *************************************************************     //
//                           Memory Handler Class
//               Handles all data required for the computation
//      *************************************************************     //
class MemoryHandler
{
	//-----------------------------
        //	    Private
        //-----------------------------
	/*Construct AM1 parameters class: Stores and Manipulate AM1 Params*/
	AM1ParameterHandler APH;

	public:
        //-----------------------------
        //	    Public
        //-----------------------------
	int NParms;

	vector<float> IPV; //Initial AM1 Parameters Vector (floats)

	/*Construct input parameters structure: Stores input parms*/
	InputParameters IP;
	flaghandler *flags;

    MemoryHandler (flaghandler *flags)
    {
        this->flags=flags;
    };

	//-----------------------------
        //    Public Class Functions
        //-----------------------------
	void LoadAPH()
	{
		APH.LoadParametersFromG09Input(IP.AM1PFile);
		IPV=APH.ProduceInitialParamsVector();
		NParms = (int)IPV.size();
	};

	void ProduceRandomGeometries()
	{
		int G = IP.G;
		float Gstd = IP.Gstd;

		int i=0;
		while (i<G-1)
		{
			//cout << "LOOPTEST-" << i+1 << " of " << G-1 << endl;
			APH.CS.push_back(APH.CS[0].normalperturbation(Gstd));
			++i;
		}
	};

	void ProduceGauInput(vector<float> &fltvec,string file,int geom)
		{APH.ProduceGaussianInput(fltvec,file,geom);}

	/*Produce a binary vector given a vector of floats*/
	vector<float> ProduceBinaryVector(vector<float> &fltvec);

	/*Produce a float vector from the bit vector*/
	vector<float> ProduceFloatVector(vector<float>  &bitvec);

	/*Function Creates a Gaussian input file for the high level of theory*/
        string ProduceControlGaussianOutputs(int i,int npg);
};
#endif

