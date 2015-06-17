#include "memhandler_class.h"
#include <iomanip>

using namespace std;

//      *************************************************************     //
//                        Memory Handler Class Functions
//      *************************************************************     //
/*----------------------------------------
  ***************************************
  | Produce an input file for gaussian  |
  ***************************************
------------------------------------------*/
string MemoryHandler::ProduceControlGaussianOutputs(int j,int npg)
{
	stringstream ss1,ss2;

	ss1 << "originalhl_gauin_geom_" << j << ".com";
	ss2 << "originalhl_gauout_geom_" << j << ".log";

	ofstream GauOut;
        GauOut.open(ss1.str().c_str());

	GauOut << "\%nproc=" << npg << endl;
	GauOut << "#P " << IP.TargetGauParams << endl;

	for (int i=1;i<(int)APH.AM1GauFilePrefix.size();++i)
		{GauOut << APH.AM1GauFilePrefix[i] << endl;}

	//Store Coords
	int NA = APH.CS[j].size();
        for (int i=0;i<NA;++i)
	{
		GauOut << APH.CS[j].getline(i) << endl;
	}

	GauOut <<  endl;
	GauOut.close();

    if (!flags->existinghl)
    {
        stringstream command;
        command << "g09 " << ss1.str() << " " << ss2.str();
        string gaurtn = exec(command.str());//Run gaussian command
    }

	return ss2.str();
};

/*----------------------------------------
  ***************************************
  |  Produce a Binary Vector of Ints    |
  ***************************************
Input a vector of floating points and return
vector of their bit representation.
------------------------------------------*/
vector<float> MemoryHandler::ProduceBinaryVector(vector<float> &fltvec)
{
        /*vector<bool> bitvec;

        union
        {
                float input;   // assumes sizeof(float) == sizeof(int)
                int   output;
        } data;

        const int S = sizeof(float) * CHAR_BIT;

        for (int i=0;i<(int)fltvec.size();++i)
        {
                data.input = fltvec[i];
                std::bitset<S> bits(data.output);
                for (int j=0;j<S;++j)
                        {bitvec.push_back(bits[j]);}
        }*/

        return fltvec;
};

/*----------------------------------------
  ***************************************
  |  Produce a Binary Vector of Ints    |
  ***************************************
Input a vector of floating points and return
vector of their bit representation.
------------------------------------------*/
vector<float> MemoryHandler::ProduceFloatVector(vector<float>  &bitvec)
{
        /*vector<float> fltvec;

        union
        {
                float input;   // assumes sizeof(float) == sizeof(int)
                int   output;
        } data;

        const int S = sizeof(float) * CHAR_BIT;
        int NumFlt = bitvec.size()/S;

        for (int i=0;i<NumFlt;++i)
        {
                std::bitset<S> bits;
                for (int j=0;j<S;++j)
                        {bits[j]=bitvec[j+i*S];}

                data.output = bits.to_ulong();

		if (data.input!=data.input || abs(data.input) > 1.0E4)
		{fltvec.push_back(IPV[i]);cout << "NAN DETECTED!! PULLING INITIAL PARAMETER VALUE!!" << endl;}
		else {fltvec.push_back(data.input);}
        }*/

        return bitvec;
};


