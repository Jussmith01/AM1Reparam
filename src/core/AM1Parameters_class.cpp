#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include "AM1Parameters_class.h"
#include "../utils/variables.h"

using namespace std;

bool testline (string line)
{
	bool test=false;

	string tline = trim(line);

	if (tline.compare("0    1")==0)
		{test=true;}
	else if (tline.compare("0   1")==0)
		{test=true;}
        else if (tline.compare("0  1")==0)
                {test=true;}
        else if (tline.compare("0 1")==0)
                {test=true;}
        else if (tline.compare("1    1")==0)
                {test=true;}
        else if (tline.compare("1   1")==0)
                {test=true;}
        else if (tline.compare("1  1")==0)
                {test=true;}
        else if (tline.compare("1 1")==0)
                {test=true;}

	return test;
};

//      *************************************************************     //
//                   AM1 Parameter Handler Class Functions		  //
//      *************************************************************     //
void AM1ParameterHandler::LoadParametersFromG09Input(string infile)
{
	cout << "APH loader running. Loading " << infile << ".\n";

	string line;
        ifstream inputfile(infile.c_str());
	int i=0; 
	bool sec1=true; 
	bool sec2=false; 
	bool sec3=false; 
	bool sec4=false; 
	bool delim=false; int delimP=-1;

	vector<string> Params;
	
	coordstore CStmp;

	//Read file into memory
        if (inputfile.is_open())
        {
                while ( getline (inputfile,line) )
                {
			if (sec1)
				{AM1GauFilePrefix.push_back(line);}
	
                        if (sec2)
                        {		
					if(trim(line).size()!=0)
                                        	{CStmp.StoreCoordsfromString(line);}
					else
						{sec3=true;sec2=false;}
			}

                        if (sec3)
                        	{AM1GauFilePostfix.push_back(line);}

			if (sec4)
                        {
				if(delim) {Params.push_back(trim(line));delim=false;}
				else {Params[delimP] += strvar::spc + trim(line);}
                        }

			if (testline(line)) {sec2=true;sec1=false;}
			if (line.compare("****")==0) {sec4=true;sec3=false;delim=true;++delimP;}
			//Parser(line,MH);
                        ++i;
		}
                inputfile.close();
        }
        else cout << "Unable to open file: " << infile.c_str() << endl;

	CS.push_back(CStmp);

	//Parse up data in memory into
	for (int k=0;k<(int)Params.size()-1;++k)
	{
		//cout << k << ") " << Params[k] << endl;
		//cout<< " " << k <<") READ TERM!!\n";
		APS.push_back(ObtainStringLetter(Params[k]));
		
		bool chkfin = true;
		while(chkfin)
		{
			string flag; string value;
			ParseString(flag,value,Params[k]);

			//cout << "LINE: " << Params[k] << " **END**\n";
			//cout << "flag: " << flag << " value: " << value << endl;
			
			APS[k].AddParameter(flag,value);

			if (trim(Params[k]).compare("****")==0) {chkfin=false;}
		}
	}
};

//      *************************************************************     //
//        Produce Parameters Vector for use in the Genetic Algorithm      //
//      *************************************************************     //
vector<float> AM1ParameterHandler::ProduceInitialParamsVector()
{
	int J = (int)APS.size();//Number of parameter sets
	vector<float> DNA;

	for (int i=0;i<J;++i)
	{
		if(APS[i].SCP[0]){DNA.push_back(APS[i].F0ss);}
		if(APS[i].SCP[1]){DNA.push_back(APS[i].F0sp);}
		if(APS[i].SCP[2]){DNA.push_back(APS[i].F0pp);}
		if(APS[i].SCP[3]){DNA.push_back(APS[i].F0sd);}
		if(APS[i].SCP[4]){DNA.push_back(APS[i].F0pd);}
		if(APS[i].SCP[5]){DNA.push_back(APS[i].F0dd);}
		if(APS[i].SCP[6]){DNA.push_back(APS[i].F2pp);}
		if(APS[i].SCP[7]){DNA.push_back(APS[i].F2pd);}
		if(APS[i].SCP[8]){DNA.push_back(APS[i].F2dd);}
		if(APS[i].SCP[9]){DNA.push_back(APS[i].F4dd);}
		if(APS[i].SCP[10]){DNA.push_back(APS[i].G1sp);}
		if(APS[i].SCP[11]){DNA.push_back(APS[i].G1pd);}
		if(APS[i].SCP[12]){DNA.push_back(APS[i].G2sd);}
		if(APS[i].SCP[13]){DNA.push_back(APS[i].G3pd);}
		if(APS[i].SCP[14]){DNA.push_back(APS[i].Rsppd);}
		if(APS[i].SCP[15]){DNA.push_back(APS[i].Rsdpp);}
		if(APS[i].SCP[16]){DNA.push_back(APS[i].Rsddd);}
		if(APS[i].NZO)
		{
			for (int k=0;k<(int)APS[i].ZetaOverlap.size();++k)
				{DNA.push_back(APS[i].ZetaOverlap[k]);}
		}
                if(APS[i].NU)
                {
                        for (int k=0;k<(int)APS[i].U.size();++k)
                                {DNA.push_back(APS[i].U[k]);}
                }
                if(APS[i].NB)
                {
                        for (int k=0;k<(int)APS[i].Beta.size();++k)
                                {DNA.push_back(APS[i].Beta[k]);}
                }
                if(APS[i].NDDN)
                {
                        for (int k=0;k<(int)APS[i].DDNval.size();++k)
                                {DNA.push_back(APS[i].DDNval[k]);}
                }
		if(APS[i].NCKO){DNA.push_back(APS[i].CoreKO);}
                if(APS[i].NKON)
                {
                        for (int k=0;k<(int)APS[i].KONval.size();++k)
                                {DNA.push_back(APS[i].KONval[k]);}
                }
		//if(APS[i].NEH){DNA.push_back(APS[i].EHeat);}
		//if(APS[i].NS){DNA.push_back(APS[i].EISol);}
		if(APS[i].NA){DNA.push_back(APS[i].Alpha);}
		if(APS[i].NDH){DNA.push_back(APS[i].DipHyp);}
                if(APS[i].NGC)
                {
                        for (int k=0;k<(int)APS[i].GCorex.size();++k)
                        {
				DNA.push_back(APS[i].GCorex[k]);
				DNA.push_back(APS[i].GCorey[k]);
				DNA.push_back(APS[i].GCorez[k]);
			}
                }
	}

	return DNA;
};

//      *************************************************************     //
//             Produce Gaussian input from an parameter vector            //
//      *************************************************************     //
void AM1ParameterHandler::ProduceGaussianInput(vector<float> &DNA,string filename,int geom)
{
	//cout << "Producing Gaussian Inputfile: " << filename << endl;

	ofstream GauOut;
	GauOut.open(filename.c_str());
	GauOut.setf( std::ios::fixed, std::ios::floatfield );

        int K = (int)AM1GauFilePrefix.size();
        for (int i=0;i<K;++i)
        	{GauOut << AM1GauFilePrefix[i] << endl;}

        K = CS[geom].size();
        for (int i=0;i<K;++i)
                {GauOut << CS[geom].getline(i) << endl;}

        K = (int)AM1GauFilePostfix.size();
        for (int i=0;i<K;++i)
                {GauOut << AM1GauFilePostfix[i] << endl;}


        int J = (int)APS.size();//Number of parameter sets
	int DNAidx=0;
        for (int i=0;i<J;++i)
        {
		GauOut << APS[i].Z << endl;

		GauOut << "PQN=";
		for(int k=0;k<(int)APS[i].PQN.size();++k)
		{GauOut << setprecision(8) << APS[i].PQN[k]; if (k!=(int)APS[i].PQN.size()-1) {GauOut << ",";}}

		GauOut << " NValence=" << APS[i].NValence;

                if(APS[i].SCP[0]){GauOut << " F0ss=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[1]){GauOut << " F0sp=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[2]){GauOut << " F0pp=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[3]){GauOut << " F0sd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[4]){GauOut << " F0pd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[5]){GauOut << " F0dd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[6]){GauOut << " F2pp=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[7]){GauOut << " F2pd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[8]){GauOut << " F2dd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[9]){GauOut << " F4dd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[10]){GauOut << " G1sp=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[11]){GauOut << " G1pd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[12]){GauOut << " G2sd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[13]){GauOut << " G3pd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[14]){GauOut << " Rsppd=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[15]){GauOut << " Rsdpp=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].SCP[16]){GauOut << " Rsddd=" << DNA[DNAidx];++DNAidx;}
		GauOut << endl;
                if(APS[i].NZO)
                {
                	GauOut << "ZetaOverlap=";
                	for(int k=0;k<(int)APS[i].ZetaOverlap.size();++k)
                	{GauOut << DNA[DNAidx];++DNAidx;if (k!=(int)APS[i].ZetaOverlap.size()-1) {GauOut << ",";}}
                }
                if(APS[i].NU)
                {
                        GauOut << " U=";
                        for(int k=0;k<(int)APS[i].U.size();++k)
                        {GauOut << DNA[DNAidx];++DNAidx;if (k!=(int)APS[i].U.size()-1) {GauOut << ",";}}
                }
		GauOut << endl;
                if(APS[i].NB)
                {
                        GauOut << "Beta=";
                        for(int k=0;k<(int)APS[i].Beta.size();++k)
                        {GauOut << DNA[DNAidx];++DNAidx;if (k!=(int)APS[i].Beta.size()-1) {GauOut << ",";}}
                }
                if(APS[i].NDDN)
                {
                        for(int k=0;k<(int)APS[i].DDNval.size();++k)
                        {GauOut << " DDN=" << APS[i].DDNL1[k] << "," << APS[i].DDNL2[k] << "," << DNA[DNAidx];++DNAidx;}
                }
		GauOut << endl;
                if(APS[i].NCKO){GauOut << "CoreKO=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].NKON)
                {
                        for(int k=0;k<(int)APS[i].KONval.size();++k)
                        {GauOut << " KON=" << APS[i].KONLT[k] << "," << APS[i].KONL1[k] << "," << APS[i].KONL2[k] << "," << DNA[DNAidx];++DNAidx;}

                }
		GauOut << endl;
                if(APS[i].NEH){GauOut << "EHeat=" << APS[i].EHeat;}
                if(APS[i].NS){GauOut << " EISol=" << APS[i].EISol;}

		GauOut << endl;
                if(APS[i].NA){GauOut << "Alpha=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].NDH){GauOut << " DipHyp=" << DNA[DNAidx];++DNAidx;}
                if(APS[i].NGC)
                {
			GauOut << endl;
                        for (int k=0;k<(int)APS[i].GCorex.size();++k)
                        {
                                GauOut << "GCore=" << DNA[DNAidx] << "," << DNA[DNAidx+1] << "," << DNA[DNAidx+2] << endl;
                                DNAidx+=3;
                        }
                }
		GauOut << "****" << endl;
        }

	//cout << "INDEX: " << DNAidx << " of " << DNA.size() << endl;

	GauOut << endl;
	GauOut.close();
};

