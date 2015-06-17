#ifndef genomestruct_header
#define genomestruct_header

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
#include "../utils/randlib.h"

using namespace std;

/*----------------------------------------
       **********************************
       |Binary Genome Population Struct |
       **********************************

------------------------------------------*/
struct Genome
{
    vector<float> strand;

    //Class Constructor
    Genome() {};

    //Class Assignment
    Genome& operator=(const Genome& instance)
    {
        this->strand = instance.strand;
        return *this;
    };

    //Function for Cross Breeding
    void CrossBreed(vector<float> &Parent1,vector<float> &Child1,vector<float> &Child2,int RandPos)
    {
        int SS = strand.size();

        for(int i=0; i<RandPos; ++i)
        {
            Child1[i] = strand[i];
        }

        for(int i=0; i<RandPos; ++i)
        {
            Child2[i] = Parent1[i];
        }

        for(int i=RandPos; i<SS; ++i)
        {
            Child1[i] = Parent1[i];
        }

        for(int i=RandPos; i<SS; ++i)
        {
            Child2[i] = strand[i];
        }
    };

    Genome GetMutant(float precision,float a,int NParms)
    {
        Genome Mutant;
        Mutant = *this;
        Mutant.Mutate(precision,a,NParms);
        return Mutant;
    };

    int Mutate(float precision,float a,int NParms)
    {   int Ngene = NParms;

        int Nmut=0;
        UniformRandomReal UR(2*Ngene,clock());
        for(int i=0; i<Ngene; ++i)
        {
            float RV = UR.UniformRandReal(0.0,1.0);
            //cout << "a: " << a << " rv: " << RV << endl;
            if (RV < a)
            {
                float RM = UR.UniformRandReal(-precision,precision);
                strand[i]*=(1.0f+RM);
                ++Nmut;
            }
        }
        return Nmut;
    };
};
#endif

