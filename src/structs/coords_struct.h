#ifndef coordstruct_header
#define coordstruct_header

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <time.h>
#include <omp.h>
#include "../utils/tools.h"
#include "../utils/randlib.h"

using namespace std;

/*----------------------------------------
       ***************************
       |  Coords Storage Struct  |
       ***************************
	
------------------------------------------*/
struct coordstore
{
	//Constructor 1
	coordstore() {};

	//Constructor 2
	coordstore(vector<string> Z,vector<float> x,vector<float> y,vector<float> z) 
	{
		int S = (int)Z.size();

		this->Z.resize(S);
		this->x.resize(S);
		this->y.resize(S);
		this->z.resize(S);

		for (int i=0;i<S;++i)
		{
			this->Z[i]=Z[i];
			this->x[i]=x[i];
			this->y[i]=y[i];
			this->z[i]=z[i];
		}
	};	

	//Assignment Operator
        coordstore& operator=(const coordstore& instance)
        {
                int S = (int)instance.Z.size();

                this->Z.resize(S);
                this->x.resize(S);
                this->y.resize(S);
                this->z.resize(S);

                for (int i=0;i<S;++i)
                {
                        this->Z[i]=instance.Z[i];
                        this->x[i]=instance.x[i];
                        this->y[i]=instance.y[i];
                        this->z[i]=instance.z[i];
                }

                return *this;
        };

	/*Atomic Letter*/
	vector<string> Z;

	/*Coordinates*/
	vector<float> x;
	vector<float> y; 
	vector<float> z;

	void StoreCoordsfromString (string line)
	{
		string ltmp = trim(line);

		//Save atomic number
		Z.push_back(trim(ltmp).substr(0,2));

		//Save x coord
		ltmp = trim(ltmp.substr(2));
		size_t pos = ltmp.find_first_of(" \t\r\n\x0b");
		string xtmp = trim(ltmp.substr(0,pos));
		x.push_back(atof(xtmp.c_str()));

		//Save y coord
		ltmp = trim(ltmp.substr(pos));
                pos = ltmp.find_first_of(" \t\r\n\x0b");
                string ytmp = trim(ltmp.substr(0,pos));
                y.push_back(atof(ytmp.c_str()));

		//Save z coord
                ltmp = trim(ltmp.substr(pos));
                pos = ltmp.find_first_of(" \t\r\n\x0b");
                string ztmp = trim(ltmp.substr(0,pos));
                z.push_back(atof(ztmp.c_str()));
	};

	//Returns number of coords in the class
	int size() {return (int)Z.size();};

	//Returns a string
	string getline(int idx)
	{
		stringstream ss;
		
		ss.setf( std::ios::fixed, std::ios::floatfield );
		ss.precision(10);
		ss <<  Z[idx] << "     " << x[idx] << "     " << y[idx] << "     " << z[idx] << "  ";

		return ss.str();
	}

	coordstore normalperturbation(float std)
	{
		int C = Z.size();
		NormRandomReal NR(3*C,clock());

		vector<string> Zt;
		vector<float> xt;
		vector<float> yt;
		vector<float> zt;

		for (int i=0;i<C;++i)
		{
			Zt.push_back(this->Z[i]);

			float xR = NR.GenRandReal(x[i],std);
			float yR = NR.GenRandReal(y[i],std);
			float zR = NR.GenRandReal(z[i],std);

			xt.push_back(this->x[i] + xR);
			yt.push_back(this->y[i] + yR);
			zt.push_back(this->z[i] + zR);
		}

		coordstore CDtmp(Zt,xt,yt,zt);
		
		return CDtmp;
	};
};
#endif
