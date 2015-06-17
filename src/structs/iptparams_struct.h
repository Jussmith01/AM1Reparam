#ifndef iptstruct_header
#define iptstruct_header

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

/*----------------------------------------
       ***************************
       | Input Parameters Struct |
       ***************************

Holds all of the input parameters.
------------------------------------------*/
struct InputParameters
{
    //*******************//
    //     Filenames     //
    //*******************//

    /*Output of the gaussian target
    calculated at a high level of
    theory - obtain target energies
    from here.*/
    string TargetGauParams;

    /*Initial AM1 Gaussian 09
    input file with added parms.*/
    string AM1PFile;

    /*Working Directory blank is current*/
    string WkDir;

    InputParameters ()
    {
        nproc=1;
        ew=1.0f;
        cew=0.0f;
        few=0.0f;
        forcebest=0;
    };

    //*******************//
    //     Set Sizes     //
    //*******************//
    /*Number of states to calculate*/
    int nproc;


    //*******************//
    //     Set Sizes     //
    //*******************//

    /*Size of the population to use*/
    int N;

    //*******************//
    //    Parameters     //
    //*******************//
    /*Number of Generations*/
    float J;

    /*Energy Weight*/
    float ew;

    /*Charge Energy Weight*/
    float cew;

    /*Force to Energy Weight*/
    float few;

    /*Percentage of Evolution*/
    float a;

    /*Percent of best subset*/
    float b;

    /*Forces inclusion of new best in every population*/
    int forcebest;

    /*This is the inital a value*/
    float InitialMutations;

    /*Precision of the inital float evolution*/
    float InitialPrecision;

    /*Precision of the mainloop float evolution*/
    float MainloopPrecision;

    /*Precision of the mainloop float evolution*/
    int G;

    /*Precision of the mainloop float evolution*/
    float Gstd;

    /*Convergence of RMSD*/
    float conv;
};
#endif

