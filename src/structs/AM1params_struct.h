#ifndef am1paramstruct_header
#define am1paramstruct_header

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
#include "../utils/tools.h"

using namespace std;

/*----------------------------------------
       ***************************
       |  AM1 Parameters Struct  |
       ***************************
	
Holds all of the AM1 parameters for a given
atom type. Used by the AM1 Paramters class.
------------------------------------------*/
struct AM1ParamsStore
{
	AM1ParamsStore(string Z)
	{
		/*Construct w. Atomic Number*/
		this->Z = Z;

		/*Construct w. defaults*/
		for(int i=0;i<17;++i)
			{SCP[i]=false;}
		NZO = false; NU = false; NB = false;
		NDDN = false; NCKO = false; NKON = false;
		NEH = false; NS = false; NDH = false;
		NA = false; NGC = false;
	};

	/*Atomic Letter*/
	string Z;

	/*Principal Quantum Numbers for each shell (s,p,d)*/
	vector<int> PQN;
	
	/*Number of valence electrons*/
        int NValence;

	/*Slater-Condon parameters for one-center two-electron integrals*/
	/*SCP stores true for param if found else false*/
	bool SCP[17];
	float F0ss;float F0sp;float F0pp;
	float F0sd;float F0pd;float F0dd;
	float F2pp;float F2pd;float F2dd;float F4dd;
	float G1sp;float G1pd;float G2sd;float G3pd;
	float Rsppd;float Rsdpp;float Rsddd;

	/*Slater exponents for basis functions*/
	/*NZO stores true if ZO calue found*/
	bool NZO; vector<float> ZetaOverlap;

	/*Diagonal core Hamiltonian matrix elements*/
	/*NU is true if param found. # of U vals (One per Ang. Momen.)*/
	bool NU; vector<float> U;

	/*Off-diagonal core Hamiltonian parameters*/
	/*NB is true if param found. # of Beta vals (One per Ang. Momen.)*/
	bool NB; vector<float> Beta;
	
        /*Point-charge distance parameters for multipole-approximated two-center two-electron integrals*/
        /*NDDN is true if param found. Form L1,L2,Value*/
        bool NDDN; vector<int> DDNL1;vector<int> DDNL2;vector<float> DDNval;

        /*Klopman-Ohno parameter used in nuclear attraction terms*/
        /*NDDN is true if param found.*/
        bool NCKO; float CoreKO;

        /*Klopman-Ohno parameters for two-center two-electron integrals*/
        /*NKON is true if param found. Form LT,L1,L2,Value*/
        bool NKON; vector<int> KONLT; vector<int> KONL1; vector<int> KONL2; vector<float> KONval;

        /*Heat of formation of the isolated atom*/
        /*NCKO hold true if value is found*/
        bool NEH; float EHeat;

        /*Energy of the isolated atom*/
        /*NS hold true if value is found*/
        bool NS; float EISol;

        /*Dipole moment hybridization parameter*/
        /*NS holds true if found*/
        bool NDH; float DipHyp;

        /*Alpha ???????*/
        /*NA holds true if found*/
        bool NA; float Alpha;

	/*GCore ???????*/
        /*NGC holds true if found*/
	bool NGC; vector<float> GCorex; vector<float> GCorey; vector<float> GCorez;

	//A function that adds a parameter value to the struct with a flag and value
	void AddParameter(string flag,string value)
	{
		if (flag.compare("PQN")==0)
		{
			vector<string> vals = commaparse(value);
			for(int i=0;i<(int)vals.size();++i)
				{PQN.push_back(atoi(vals[i].c_str()));}	
		} else if (flag.compare("NValence")==0) {NValence=atof(value.c_str());
		} else if (flag.compare("F0ss")==0) {SCP[0]=true;F0ss=atof(value.c_str());
		} else if (flag.compare("F0sp")==0) {SCP[1]=true;F0sp=atof(value.c_str());
		} else if (flag.compare("F0pp")==0) {SCP[2]=true;F0pp=atof(value.c_str());
		} else if (flag.compare("F0sd")==0) {SCP[3]=true;F0sd=atof(value.c_str());
		} else if (flag.compare("F0pd")==0) {SCP[4]=true;F0pd=atof(value.c_str());
		} else if (flag.compare("F0dd")==0) {SCP[5]=true;F0dd=atof(value.c_str());
		} else if (flag.compare("F2pp")==0) {SCP[6]=true;F2pp=atof(value.c_str());
		} else if (flag.compare("F2pd")==0) {SCP[7]=true;F2pd=atof(value.c_str());
		} else if (flag.compare("F2dd")==0) {SCP[8]=true;F2dd=atof(value.c_str());
		} else if (flag.compare("F4dd")==0) {SCP[9]=true;F4dd=atof(value.c_str());
		} else if (flag.compare("G1sp")==0) {SCP[10]=true;G1sp=atof(value.c_str());
		} else if (flag.compare("G1pd")==0) {SCP[11]=true;G1pd=atof(value.c_str());
		} else if (flag.compare("G2sd")==0) {SCP[12]=true;G2sd=atof(value.c_str());
		} else if (flag.compare("G3pd")==0) {SCP[13]=true;G3pd=atof(value.c_str());
		} else if (flag.compare("Rsppd")==0) {SCP[14]=true;Rsppd=atof(value.c_str());
		} else if (flag.compare("Rsdpp")==0) {SCP[15]=true;Rsdpp=atof(value.c_str());
		} else if (flag.compare("Rsddd")==0) {SCP[16]=true;Rsddd=atof(value.c_str());
		} else if (flag.compare("ZetaOverlap")==0) {
			NZO=true; vector<string> vals = commaparse(value);
                        for(int i=0;i<(int)vals.size();++i)
                                {ZetaOverlap.push_back(atof(vals[i].c_str()));}
                } else if (flag.compare("U")==0) {
                        NU=true; vector<string> vals = commaparse(value);
                        for(int i=0;i<(int)vals.size();++i)
                                {U.push_back(atof(vals[i].c_str()));}
                } else if (flag.compare("Beta")==0) {
                        NB=true; vector<string> vals = commaparse(value);
                        for(int i=0;i<(int)vals.size();++i)
                                {Beta.push_back(atof(vals[i].c_str()));}
                } else if (flag.compare("DDN")==0) {
                        NDDN=true; vector<string> vals = commaparse(value);
			DDNL1.push_back(atoi(vals[0].c_str()));
			DDNL2.push_back(atoi(vals[1].c_str()));
			DDNval.push_back(atof(vals[2].c_str()));
                } else if (flag.compare("CoreKO")==0) {NCKO=true;CoreKO=atof(value.c_str());
		} else if (flag.compare("KON")==0) {
                        NKON=true; vector<string> vals = commaparse(value);
                        KONLT.push_back(atoi(vals[0].c_str()));
                        KONL1.push_back(atoi(vals[1].c_str()));
                        KONL2.push_back(atoi(vals[2].c_str()));
                        KONval.push_back(atof(vals[3].c_str()));
                } else if (flag.compare("EHeat")==0) {NEH=true;EHeat=atof(value.c_str());
                } else if (flag.compare("EISol")==0) {NS=true;EISol=atof(value.c_str());
                } else if (flag.compare("DipHyp")==0) {NDH=true;DipHyp=atof(value.c_str());
                } else if (flag.compare("Alpha")==0) {NA=true;Alpha=atof(value.c_str());
		} else if (flag.compare("GCore")==0) {
                        NGC=true; vector<string> vals = commaparse(value);
                        GCorex.push_back(atof(vals[0].c_str()));
                        GCorey.push_back(atof(vals[1].c_str()));
                        GCorez.push_back(atof(vals[2].c_str()));
		}
	};
};
#endif
