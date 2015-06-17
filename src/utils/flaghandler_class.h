#ifndef flaghandler_header
#define flaghandler_header

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

using namespace std;

//      *************************************************************     //
//                            Flag Handler Class
//                  	    Handles all input flags
//      *************************************************************     //
class flaghandler
{
    //--------------------------
    //   Vector to hold flags
    //--------------------------
    vector<string> flags;

    //------------------------------
    //Private Member Class Functions
    //------------------------------
    void FlagParser()
    {
        //Set Defaults
        inputfilename = "DEFAULT.in";
        outputfilename = "DEFAULT.out";
        verbosity = 0;
        existinghl = false;

        for (int i=1; i<(int)flags.size(); ++i)
        {
            if (flags[i].compare("-i")==0)
            {
                inputfilename=flags[i+1];
            }
            else if (flags[i].compare("-o")==0)
            {
                outputfilename=flags[i+1];
            }
            else if (flags[i].compare("-v")==0)
            {
                verbosity=atoi(flags[i+1].c_str());/*cout << "VERB: " << verbosity << endl;*/
            }
            else if (flags[i].compare("--existinghl")==0)
            {
                existinghl=true;
            };
        }
    };

public:
    //--------------------------
    //     Public Variable
    //--------------------------
    string inputfilename;
    string outputfilename;
    int verbosity;
    bool existinghl;

    //-----------------------------
    //	Class Constructor
    //-----------------------------
    flaghandler(int argc, char *argv[])
    {
        for (int i=0; i<argc; ++i)
        {
            flags.push_back(argv[i]);
        }

        FlagParser();
    };
};
#endif

