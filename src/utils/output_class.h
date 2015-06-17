#ifndef output_header
#define output_header

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include "../core/memhandler_class.h"

using namespace std;

//      *************************************************************     //
//                            Program Output File
//               Holds program output variables and functions
//      *************************************************************     //
class outputFile
{
	public:
        //--------------------------
        //Private Class Declarations
        //--------------------------
	ofstream ofile;

        //------------------------------
        //Private Member Class Functions
        //------------------------------
	int verbosity;

        //-----------------------------
        //      Class Constructor
        //-----------------------------
	outputFile(string filename,int verbosity)
	{
		this->verbosity = verbosity; 

        	//Define output name string
        	string actname;

                //Check if name is not included as the argument, then set output name.
                if (filename.compare("")==0) 
			{actname = "DEFAULT_OUTPUT.opt";} 
		else 
			{actname = filename;}

                //Open the output file by given filename.
                ofile.open(actname.c_str());
	};
        //-----------------------------
        //Public Member Class Functions
        //-----------------------------

	/*____________________________________________________
		----Close Output File       ----
		----Author: Justin Smith           ----
		----Date Modified: 5/11/2015       ----
		----Modified By:                   ----
	This saves the data and closes the output file. Do this 
	before ending the program. When input int close_type = 0 
	closes normally, when 1, abnormal termination.
	*/
	void close_output (int close_type)
	{
	        switch(close_type)
	        {
	                case 0: {ofile << "Normal Termination\n";break;}
	                case 1: {ofile << "Abnormal Termination\n";break;}
	        }

	        ofile.close();
	};

	/*____________________________________________________
		----Error Print and Exit Prog    ----
		----Author: Justin Smith         ----
		----Date Modified: 5/11/2015     ----
		----Modified By:                 ----
	This prints a closing error message and ends the program. 
	*/
	void Prog_err_quit(string Message)
	{
	        ofile << "\n" << Message << "\n";
	        ofile << "Freeing Memory Allocations...\n";
	        close_output(1);
	        exit(1);
	};

        /*____________________________________________________
                ----Print Input Parameters       ----
                ----Author: Justin Smith         ----
                ----Date Modified: 5/12/2015     ----
                ----Modified By:                 ----
        This prints the input parameters to the output file. The
	parameters are stored in the memory handler class.
        */
        void PrintIptParms(MemoryHandler &MH)
        {
		ofile << "________________________________________\n\n";
                ofile << "Program Input Parameters: \n";
		ofile << "Working Directory: " << MH.IP.WkDir << endl;
		ofile << "Target Gaussian Parameters: " << MH.IP.TargetGauParams << endl;
		ofile << "AM1 Gaussian Infile: " << MH.IP.AM1PFile << endl;
		ofile << "J (Maximum Generations): " << MH.IP.J << endl;
		ofile << "N (Population Size): " << MH.IP.N << endl;
		ofile << "Force Inclusion of Best: " << MH.IP.forcebest << endl;
		ofile << "G (Geometries to Fit): " << MH.IP.G << endl;
		ofile << "Gstd (std. dev. of geoms.): " << MH.IP.Gstd << endl;
		ofile << "Charge Energy Weight: " << MH.IP.cew << endl;
		//ofile << "K (Number of States): " << MH.IP.K << endl;
		ofile << "a (Mutation Percent): " << MH.IP.a << endl;
		ofile << "b (Keep Best %): " << MH.IP.b << endl;
		ofile << "Initial Mutations: " << MH.IP.InitialMutations << endl;
		ofile << "Initial Precision: " << MH.IP.InitialPrecision << endl;
		ofile << "Mainloop Precision: " << MH.IP.MainloopPrecision << endl;
		ofile << "RMSD Convergence: " << MH.IP.conv << endl;
		ofile << "________________________________________\n\n";
        };
};
#endif

