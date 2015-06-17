#ifndef RandHeader
#define RandHeader

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <random>
#include <omp.h>
#include <math.h>

//      *************************************************************     //
//                   Functions Random Number Generation
//      *************************************************************     //

/*----------------------------------------
  ***************************************
  |  Class for generating random Ints   |
  ***************************************
Example Use:
RandomInt RI(w,i); //int w=num seeds
		   //int i=thread seed

int someInt = GenRandInt(high,low)
		   //Parameters give the
		   //range of the int rtn'd
------------------------------------------*/
class RandomInt
{
        std::default_random_engine generator;
        std::vector<int> array;
        int index;

        public:
        RandomInt(){};
        RandomInt(int w,int i){Setup(w,i);};

        void Setup(int w,int i)
        {
                time_t Time;
                time(&Time);
                int seedOffset=(int)Time;

                array.resize(w);

                int t = (int)omp_get_wtime()+i;
                std::seed_seq seed = {seedOffset,t,i+100};
                seed.generate(array.begin(),array.end());//Seed the generator
                index = 0;
        };

        int GenRandInt(int high,int low)
        {
                generator.seed(array[index]);//Seed the generator
                std::uniform_int_distribution<int> distribution(low,high);//Setup the distribution
                int RN = (int)distribution(generator);//Denerate the random number
                //std::cout << "RandomNumber: " << RN << std::endl;
                ++index;//Increase seed offset
                return RN;
        };
};

/*-----------------------------------------------
  ***********************************************
  |Class for Generating a Flt from a Normal Dist|
  ***********************************************
Example Use:
NormRandomReal NR(w,i); //int w=num seeds
                    	//int i=thread seed

float someflt = GenRandReal(float mean,float std)
                   //mean is mean, duh
                   //std is standard deviation
-------------------------------------------------*/
class NormRandomReal
{
        std::default_random_engine generator;
        std::vector<int> array;
        int index;

        public:
        NormRandomReal(){};
        NormRandomReal(int w,int seed){Setup(w,seed);};

        void Setup(int w,int i)
        {
                time_t Time;
                time(&Time);
                int seedOffset=(int)Time;

                array.resize(w);

                int t = (int)omp_get_wtime()+i;
                std::seed_seq seed = {seedOffset,t,i+100};
                seed.generate(array.begin(),array.end());//Seed the generator
                index = 0;
        };

        float GenRandReal(float mean,float std)
        {
                generator.seed(array[index]);//Seed the generator
                std::normal_distribution<float> distribution(mean,std);//Setup the distribution
                float RN = (float)distribution(generator);//Denerate the random number
                //std::cout << "RandomNumber: " << RN << std::endl;
                ++index;//Increase seed offset
                return RN;
        };
};

/*-----------------------------------------------
  ***********************************************
  |Class for Generating a Flt from a Normal Dist|
  ***********************************************
Example Use:
UniformRandomReal UR(w,i); //int w=num seeds
                           //int i=thread seed

float someflt = UniformRandReal(float low,float high)
                   //Parameters give the
                   //range of the real rtn'd
-------------------------------------------------*/
class UniformRandomReal
{
        std::default_random_engine generator;
        std::vector<int> array;
        int index;

        public:
        UniformRandomReal(){};
        UniformRandomReal(int w,int seed){Setup(w,seed);};

        void Setup(int w,int i)
        {
                time_t Time;
                time(&Time);
                int seedOffset=(int)Time;

                array.resize(w);

                int t = (int)omp_get_wtime()+i;
                std::seed_seq seed = {seedOffset,t,i+100};
                seed.generate(array.begin(),array.end());//Seed the generator
                index = 0;
        };

        float UniformRandReal(float low,float high)
        {
                generator.seed(array[index]);//Seed the generator
                std::uniform_real_distribution<float> distribution(low,high);//Setup the distribution
                float RN = (float)distribution(generator);//Denerate the random number
                //std::cout << "RandomNumber: " << RN << std::endl;
                ++index;//Increase seed offset
                return RN;
        };
};
#endif
