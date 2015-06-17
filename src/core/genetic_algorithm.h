#ifndef GA_header
#define GA_header

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
#include <iomanip>
#include <limits.h>
#include <bitset>
#include <time.h>
#include "memhandler_class.h"
#include "../utils/tools.h"
#include "../utils/randlib.h"
#include "../utils/output_class.h"
#include "../structs/genomes_struct.h"

using namespace std;

//      *************************************************************     //
//                           Memory Handler Class
//               Handles all data required for the computation
//      *************************************************************     //
class GeneticAlgorithm
{
  //*****************************//
  //	     Private	       //
  //*****************************//

  //-----------------------------
  //     Registered Classes
  //-----------------------------
  /*Register Memory Handler Class*/
  MemoryHandler *MH;
  /*Register Output Handler Class*/
  outputFile *out;

  //-----------------------------
  //         Ancestry
  //-----------------------------
  Genome OA;          // Oldest Ancestor
  vector<Genome> CP;  // Current Generation
  vector<Genome> NP;  // Next Generation

  vector<bool> BL;  // Genome Black List
  int NBL;

  float mu;
  float std;

  float bestfitness;
  float besteVFitness;
  float bestForceFitness;
  float bestChargeFitness;

  float orgam1fit;
  float orgam1ev;
  float orgam1charge;
  float orgam1force;

  //-----------------------------
  //      Storage Vectors
  //-----------------------------
  vector<vector<string>> inputfiles;  // Input filenames
  vector<vector<string>> outputfiles;  // Output filenames

  vector<vector<float>> eVControl;  // Stores the control set values
  vector<vector<float>> eVAM1Original;  //

  vector<vector<float>> fControl;  // Stores original force values
  vector<vector<float>> fAM1Original;  //

  vector<vector<float>> parameters;  // Store the parameters in float

  //-----------------------------
  //      Class Functions
  //-----------------------------
  // Obtain excited state energies and ground state charges
  /*These values are used to determine fitness*/
  vector<float> excitedStateEnergies(string something);

  // Fitness Function
  /*Function determines the fitness values of a set of
    values obtained from the excitedStateEnergies function*/
  float eVFit(std::vector<float> eVTest, std::vector<float> eVControl);

  // creates an index of survivors
  /*Function returns the index of survivors. Also fills fitvals
    with the fitness values of each survivor. b=number of surv-
    ivors to find.*/
  vector<int> survivors(vector<float> &fitvals, vector<float> &evfitvals, vector<float> &cfitvals,
                        vector<float> &ffitvals, int b);

  // Creates a vector reprenting the forces
  std::vector<float> fReader(string inFile);

  // Calculate the charge fitness
  float cFit(std::vector<float> eVTest, std::vector<float> eVControl);

  // Calculates the force fitness
  float fFit(std::vector<float> fTest, std::vector<float> forceControl);

  public:
  //-----------------------------
  //	  Constructor
  //-----------------------------
  /*The class constructor registers the memory handler and the output handler*/
  GeneticAlgorithm(MemoryHandler *MH, outputFile *out)
  {
    this->MH = MH;
    this->out = out;
  };

  //-----------------------------
  //    Public Class Functions
  //-----------------------------
  /*Prepares the genetic algorithm*/
  void GeneticAlgorithmPreparation(void);

  /*Begin main phase of the genetic algorithm*/
  void GeneticAlgorithmMainloop(void);
};
#endif
