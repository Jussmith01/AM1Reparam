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
#include <limits.h>
#include <bitset>
#include <unistd.h>
#include <time.h>
#include "../utils/timer_class.h"
#include "genetic_algorithm.h"

using namespace std;

//      *************************************************************     //
//                 Functions for the Genetic Algorithm Class
//      *************************************************************     //

/*-----------------------------------------
  *****************************************
  |Function prepares the genetic algorithm|
  *****************************************
Written by: Justin Smith
Last modified by: Justin Smith

1)Creates Oldest Ancestor Gene (OAG)
2)Produces initial population N from (OAG)
3)Obtain target excited energies and ground
  state charges. TESTING
-------------------------------------------*/
void GeneticAlgorithm::GeneticAlgorithmPreparation()
{
    out->ofile << "\n________________________________________\n";
    out->ofile << "       Preparing Genetic Algorithm       \n\n";
    out->ofile.setf(std::ios::fixed, std::ios::floatfield);
    out->ofile << setprecision(7);

    //***************
    // Initial Defines
    //***************
    int N = MH->IP.N;
    int G = MH->IP.G;
    int P = MH->NParms;
    string dir = MH->IP.WkDir;

    out->ofile << "Total Parameters being Optimized: " << P << endl;

    //***************
    // Class Allocation
    //****************
    inputfiles.resize(N);
    outputfiles.resize(N);
    CP.resize(N);
    NP.resize(N);
    BL.resize(N);
    eVControl.resize(G);
    eVAM1Original.resize(G);
    fControl.resize(G);
    fAM1Original.resize(G);
    parameters.resize(N);
    for (int i = 0; i < N; ++i)
    {
        NP[i].strand.resize(P * 32);
        parameters[i].resize(P);
    }

    for (int parameterSet = 0; parameterSet < N; ++parameterSet)
    {
        inputfiles[parameterSet].resize(G);
        outputfiles[parameterSet].resize(G);
    }

    //************************************
    // Produce Original Ancestor(OA) Genome
    //************************************
    OA.strand = MH->ProduceBinaryVector(MH->IPV);

    //*************************
    // Produce Working Filenames
    //*************************
    for (int parameterSet = 0; parameterSet < N; ++parameterSet)
    {
        for (int geometry = 0; geometry < G; ++geometry)
        {
            stringstream ssi, sso;
            ssi << "gauin_mbr" << parameterSet << "_geom" << geometry;
            sso << "gauout_mbr" << parameterSet << "_geom" << geometry;

            stringstream a1, a2;
            a1 << dir << ssi.str() << ".com";
            a2 << dir << sso.str() << ".log";

            inputfiles[parameterSet][geometry] = a1.str();
            outputfiles[parameterSet][geometry] = a2.str();
        }
    }

    //*************************
    // Produce Random Geometries
    //*************************
    MH->ProduceRandomGeometries();

    //**************************
    // Setup target fitness array
    //**************************
    vector<float> begfits;
    begfits.resize(G);
    vector<float> begevfits;
    begevfits.resize(G);
    vector<float> begffits;
    begffits.resize(G);
    vector<float> begcfits;
    begcfits.resize(G);
    vector<float> begevfitsam1;
    begevfitsam1.resize(G);
    vector<float> begcfitsam1;
    begcfitsam1.resize(G);
    vector<float> begffitsam1;
    begffitsam1.resize(G);

    // Setup Threads
    int np = MH->IP.nproc;
    omp_set_num_threads(np);
    int npg =
        (int)floor(np / (float)G);  // Number of threads to allow per gaussian run
    npg = (npg != 0) ? npg : 1;
    out->ofile << "Number of Threads per High-level Gaussian Calculation: " << npg
               << "\n";

    out->ofile << "Calculating High-level Control Set and Original AM1 Fitness "
               "Values...\n";
    #pragma omp parallel for shared(G)
    for (int i = 0; i < G; ++i)  // Setup to be parallelized
    {
        //-------------------
        // Run initial AM1 G09
        //-------------------
        stringstream infile, outfile;
        infile << dir << "orgam1_gauin_" << i << ".com";
        outfile << dir << "orgam1_gauout_" << i << ".log";

        MH->ProduceGauInput(MH->IPV, infile.str(), i);

        stringstream command2;
        command2 << "g09 " << infile.str() << " "
                 << outfile.str();  // Build command
        // cout << " Running Original Gau AM1..." << endl;
        if (execg09(command2.str(), dir, i))
        {
            out->ofile << "ERROR **Original AM1 Parameters failed to converge for a "
                       "geometry**\n";
            out->ofile << "This will give an incorrect original fitness\n";
            out->ofile << "Try lowering the 'Gstd' value in the input\n";
            out->ofile.close();
            exit(1);
        };  // Run gaussian command

        eVAM1Original[i] = excitedStateEnergies(outfile.str());
        fAM1Original[i] = fReader(outfile.str());

        //------------------
        // Run HL Control G09
        //------------------
        string ctrfile = dir +
                         MH->ProduceControlGaussianOutputs(
                             i, npg);  // Turn TargetGFile into vector
        eVControl[i] = excitedStateEnergies(ctrfile);
        fControl[i] = fReader(ctrfile);

        //------------------------------------
        // Obtain Beginning Fits per Geometry
        //------------------------------------
        begevfits[i] = eVFit(excitedStateEnergies(outfile.str()), eVControl[i]);
        begffits[i] = fFit(fReader(outfile.str()), fControl[i]);
        begcfits[i] = cFit(excitedStateEnergies(outfile.str()), eVControl[i]);

        begfits[i] = begevfits[i] + begffits[i] + begcfits[i];

        begevfitsam1[i] =
            eVFit(excitedStateEnergies(outfile.str()), eVAM1Original[i]);
        begffitsam1[i] = fFit(fReader(outfile.str()), fAM1Original[i]);
        begcfitsam1[i] =
            cFit(excitedStateEnergies(outfile.str()), eVAM1Original[i]);

        stringstream r1, r2;

        r1 << "rm " << dir << infile.str();
        r2 << "rm " << dir << outfile.str();

        exec(r1.str());  // Delete Input Files
        exec(r2.str());  // Delete Output Files
    }

    cout << "Starting fitness testing\n";
    //*********************************************
    // Initialize the best fitness as the original
    //*********************************************
    besteVFitness = 0;
    bestForceFitness = 0;
    bestChargeFitness = 0;

    for (int i = 0; i < G; ++i)
    {
        besteVFitness += begevfits[i];
        bestForceFitness += begffits[i];
        bestChargeFitness += begcfits[i];
    }

    bestfitness = besteVFitness + bestForceFitness + bestChargeFitness;

    if (bestfitness != bestfitness)
    {
        out->ofile << "NaN Detected in original bestfitness value\n";
        exit(1);
    }

    orgam1ev = 0;
    orgam1force = 0;
    orgam1charge = 0;

    for (int i = 0; i < G; ++i)
    {
        orgam1ev += begevfitsam1[i];
        orgam1force += begffitsam1[i];
        orgam1charge += begcfitsam1[i];
    }

    orgam1fit = orgam1ev + orgam1force + orgam1charge;

    //*********************************************
    //     Print the best fitness original Data
    //*********************************************
    // Print total new best Data for AM1->HL
    out->ofile << "\n-------------------------------------\n";
    out->ofile << "Original (AM1 -> HL): " << bestfitness << endl;
    ;
    out->ofile << "    Geom #: {   Fit   ,  eVFit  ,  fFit   ,  cFit   }\n";

    for (int geometry = 0; geometry < G; geometry++)
    {
        out->ofile << "    Geom " << geometry << ": {";
        out->ofile << begfits[geometry] << ",";
        out->ofile << begevfits[geometry] << ",";
        out->ofile << begffits[geometry] << ",";
        out->ofile << begcfits[geometry] << "}" << endl;
    }

    out->ofile << "     Total: {";
    out->ofile << bestfitness << ",";
    out->ofile << besteVFitness << ",";
    out->ofile << bestForceFitness << ",";
    out->ofile << bestChargeFitness << "}" << endl;

    out->ofile << "\nOriginal AM1->AM1 Fitness(Should be 0): " << orgam1fit;
    out->ofile << " Excited Stated Fitness: " << orgam1ev << "\n";
    out->ofile << " Force Fitness: " << orgam1force << "\n";
    out->ofile << " Charge Fitness: " << orgam1charge << "\n\n";

    //***********************************************************
    // Produce initial mutant population and run gaussian for each
    //***********************************************************
    // Setup openmp
    out->ofile << "Producing Initial Population...\n";

    float IM = MH->IP.InitialMutations;
    float IP = MH->IP.InitialPrecision;

    #pragma omp parallel for shared(P, G)
    for (int i = 0; i < N; ++i)
    {
        CP[i] = OA.GetMutant(IP, IM, P);  // Get mutant
        vector<float> IPV2 =
            MH->ProduceFloatVector(CP[i].strand);  // Produce vector of floats

        BL[i] = false;

        for (int j = 0; j < G; ++j)
        {
            MH->ProduceGauInput(IPV2, inputfiles[i][j],
                                j);  // Produce the gaussian inputfile

            stringstream command;
            command << "g09 " << inputfiles[i][j].c_str() << " "
                    << outputfiles[i][j].c_str();  // Build command
            // cout << "EXEC: " << command.str() << endl;

            if (execg09(command.str(), dir, i))
            {
                BL[i] = true;  // Run gaussian command
            }
        }
    }

    NBL = CountFails(BL);
    out->ofile << "Number of Parameter sets that Failed to Converge: " << NBL
               << endl;
};

/*----------------------------------------
  ***************************************
  |	Genetic Algorithm Main Loop	|
  ***************************************
Written by: Justin Smith
Last modified by: Justin Smith
------------------------------------------*/
void GeneticAlgorithm::GeneticAlgorithmMainloop()
{
    int J = MH->IP.J;
    int N = MH->IP.N;
    int G = MH->IP.G;
    float a = MH->IP.a;
    int b = N * MH->IP.b;
    int P = MH->NParms;
    string dir = MH->IP.WkDir;

    /*This is variable dependent upon the std. dev.*/
    float precbits =
        MH->IP.MainloopPrecision;  // This value should shrink per generation

    out->ofile << "\n________________________________________\n";
    out->ofile << "      Beginning Genetic Algorithm       \n\n";

    //---------------//
    // Begin Main Loop//
    //---------------//
    int cntr = 0;  // Generation counter
    while (true)
    {
        out->ofile << "****************\n";
        out->ofile << "Generation: " << cntr << endl;
        out->ofile << "****************\n\n";

        // Cycle Timer
        timer LoopTimer;
        LoopTimer.set_timer();

        cout << "Generation: " << cntr + 1 << " of " << J << "\n";

        //*********************************
        // Determine probability of survival
        //*********************************
        double t1 = omp_get_wtime();
        out->ofile << "Choosing Survivors... \n";
        vector<float> fitvals;
        vector<float> evfitvals;
        vector<float> cfitvals;
        vector<float> ffitvals;
        vector<int> Survivors =
            survivors(fitvals, evfitvals, cfitvals, ffitvals, b);

        if (MH->flags->verbosity == 1)
        {
            for (int i = 0; i < b; ++i)
            {
                out->ofile << " Survivors: " << Survivors[i]
                           << " Fitness: " << fitvals[i] << endl;
            }
        }

        out->ofile << " Survivors Average: " << mean(fitvals, 0) << endl;
        out->ofile << "  eV Fit Average: " << mean(evfitvals, 0) << endl;
        out->ofile << "  Charge Fit Average: " << mean(cfitvals, 0) << endl;
        out->ofile << "  Force Fit Average: " << mean(ffitvals, 0) << endl;
        out->ofile << " Survivors Time: " << omp_get_wtime() - t1 << "s" << endl;

        if (J - 1 == cntr)
        {
            break;
        };

        //**************************************
        // Build new population through crossover
        //**************************************
        t1 = omp_get_wtime();
        RandomInt RI(5 * N, clock());  // Class for generating random parents
        out->ofile << "\nBeginning Crossover... \n";
        for (int i = 0; i < N / 2; ++i)
        {
            int P1 = RI.GenRandInt(b - 1, 0);
            int P2 = RI.GenRandInt(b - 1, 0);

            while (P1 == P2)
            {
                P2 = RI.GenRandInt(b - 1, 0);
            }

            int RP = RI.GenRandInt(P - 1, 0);

            int idxP1 = Survivors[P1];
            int idxP2 = Survivors[P2];

            // out->ofile << " Splitting Parent1: " << idxP1 << " Parent2: " << idxP2
            // << " at " << RP << endl;

            // Cross breed Parent1 (P1) with Parent2 (P2)
            CP[idxP1].CrossBreed(CP[idxP2].strand, NP[2 * i].strand,
                                 NP[2 * i + 1].strand, RP);
        }
        out->ofile << " Crossover Time: " << omp_get_wtime() - t1 << "s" << endl;

        //*********************
        // Mutate New Population
        //*********************
        out->ofile << "\nBeginning Population Mutation... \n";
        t1 = omp_get_wtime();
        vector<int> Nmutations;
        for (int i = 0; i < N; ++i)
        {
            Nmutations.push_back((float)NP[i].Mutate(precbits, a, P));
            // out->ofile << " Genome(" << i << ") Mutations: " << Nmutation << endl;
        }

        out->ofile << " Average Mutations: " << mean(Nmutations, 0) << endl;
        out->ofile << " Pop. Mut. Time: " << omp_get_wtime() - t1 << "s" << endl;
        //*****************
        // Run Gaussian Jobs
        //*****************
        out->ofile << "\nRunning Gaussian Jobs... \n";
        t1 = omp_get_wtime();
        #pragma omp parallel for
        for (int i = 0; i < N; ++i)
        {
            BL[i] = false;

            vector<float> IPVtmp = MH->ProduceFloatVector(NP[i].strand);

            for (int p = 0; p < P; ++p)
            {
                parameters[i][p] = IPVtmp[p];
            }

            for (int j = 0; j < G; ++j)
            {
                MH->ProduceGauInput(IPVtmp, inputfiles[i][j], j);
                stringstream command;
                command << "g09 " << inputfiles[i][j].c_str() << " "
                        << outputfiles[i][j].c_str();
                if (execg09(command.str(), dir, i))
                {
                    BL[i] = true;  // Run gaussian command
                }
            }

            CP[i] = NP[i];  // Transfer New population to the Current Population
        }

        NBL = CountFails(BL);
        out->ofile << " Number of Parameter sets that Failed to Converge: " << NBL
                   << endl;

        out->ofile << "\nGaussian Jobs Complete... \n";
        out->ofile << " Gaussian Jobs Time: " << omp_get_wtime() - t1 << "s"
                   << "\n\n";

        //******************
        // Clean up the cycle
        //******************
        LoopTimer.end_timer();  // End the program wall timer
        LoopTimer.print_time("Cycle Time: ",
                             out->ofile);  // Print the wall time to output

        out->ofile << "\n|---------End Cycle---------|\n\n";
        // Run Gaussian
        ++cntr;
    }

    out->ofile << "\n|---------End Cycle---------|\n\n";

    out->ofile << "!!!!!!!GENETIC ALGORITHM DONE!!!!!!!\n";

    out->ofile << "\nBest Fitness Obtained: " << bestfitness << "\n";
    out->ofile << "   with Force Fitness : " << bestForceFitness << "\n";
    out->ofile << "\nDifference from Original AM1: " << orgam1fit << "\n";
    out->ofile << "   with Force Fitness : " << orgam1force << "\n";
    out->ofile << "Output of best stored as: BESTOUTPUT.log"
               << "\n";

    out->ofile << "\nCleaning up...\n";
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < G; ++j)
        {
            stringstream command;
            command << "rm " << outputfiles[i][j].c_str();
            string gaurtn1 = exec(command.str());
            // system(command.str().c_str());

            stringstream command2;
            command2 << "rm " << inputfiles[i][j].c_str();
            string gaurtn2 = exec(command2.str());
            // system(command2.str().c_str());
        }
    }
};

//******************//
// EREADER FUNCTION //
//******************//

// The excitedStateEnergies function reads the gaussian output file
// and return a vector containing the excited state energies

vector<float> GeneticAlgorithm::excitedStateEnergies(string inputFile)
{
    // open file gaussian input file
    std::ifstream inFile;
    inFile.open(inputFile.c_str());
    float energyWeight = MH->IP.ew;

    std::vector<float> excitedStates;
    if (energyWeight > 1.0E-5)
    {
        // Check for errors in file opening
        if (inFile.fail())
        {
            std::cerr << "bad read for excitedStateEnergies -- File not found: "
                      << inputFile.c_str() << std::endl;
            std::exit(1);
        }

        // Search file for excited states and add to vector eV
        std::string search;
        float eVi;  // the number of excited state energies
        while (!inFile.eof())
        {
            std::getline(inFile, search);

            // The commented section is used to read ground
            // state energies
            /*************************************************
            if (search.find("SCF Done") < 100) {
              std::string scfs;
              float scfn;
              scfs = search.substr(search.find("=") + 2, search.find("A.U.")
                                   - search.find("="));
              scfn = atof(scfs.c_str());
              excitedStates.push_back(scfn);
            }
            **************************************************/
            if (search.find("Excited State") < 100)
            {
                std::string eVs;
                float eVn;
                eVs = search.substr(search.find("eV") - 7, 6);
                eVn = atof(eVs.c_str());
                excitedStates.push_back(eVn);
                eVi = excitedStates.size();  // count the nuber of
                // of excited states
            }

            if (search.find("Mulliken charges:") < 100)
            {
                std::getline(inFile, search);
                while (search.find("Sum") > 100)
                {
                    std::getline(inFile, search);
                    if (search.find("Sum") > 100)
                    {
                        std::string charges;
                        float chargen;
                        charges = search.substr(12, 9);
                        chargen = atof(charges.c_str());
                        excitedStates.push_back(chargen);
                    }
                }
            }
        }

        excitedStates.push_back(eVi);  // add the number of excited
        // state as the last term of the vector
    }

    return excitedStates;
}

//*******************//
// FRREADER FUNCTION //
//*******************//

// Reads the forces from an input and returns a vector of forces <x,y,z>

std::vector<float> GeneticAlgorithm::fReader(string inputFile)
{
    std::vector<float> xvector;
    std::vector<float> yvector;
    std::vector<float> zvector;
    float forceEnergyWeight = MH->IP.few;
    if (forceEnergyWeight > 1.0E-5)
    {
        {
            std::ifstream inFile(inputFile.c_str());
            if (!inFile)
            {
                std::cerr << "bad read for fReader -- File not found: "
                          << inputFile.c_str() << std::endl;
            }

            std::string search;

            while (inFile)
            {
                std::getline(inFile, search);

                if (search.find("Forces") != string::npos)
                {
                    std::getline(inFile, search);
                    std::getline(inFile, search);
                    std::getline(inFile, search);
                    while (search.find("Cartesian") == string::npos)
                    {
                        std::string xs;
                        std::string ys;
                        std::string zs;
                        float x;
                        float y;
                        float z;
                        xs = search.substr(26, 12);
                        x = atof(xs.c_str());
                        xvector.push_back(x);
                        ys = search.substr(41, 12);
                        y = atof(ys.c_str());
                        yvector.push_back(y);
                        zs = search.substr(56, 12);
                        z = atof(zs.c_str());
                        zvector.push_back(z);
                        std::getline(inFile, search);
                    }
                }
            }
        }

        xvector.erase(xvector.end() - 1);
        yvector.erase(yvector.end() - 1);
        zvector.erase(zvector.end() - 1);
        xvector.insert(xvector.end(), yvector.begin(), yvector.end());
        xvector.insert(xvector.end(), zvector.begin(), zvector.end());
    }
    return xvector;
}

//****************//
// EVFIT FUNCTION //
//****************//

// The eVFit function produces a value that will be used
// as a measure of fittness of a parameter set

float GeneticAlgorithm::eVFit(std::vector<float> eVTest,
                              std::vector<float> eVControl)
{
    float energyWeight = MH->IP.ew;
    float fitValue = 0.0;
    if (energyWeight > 1.0E-5)
    {
        {
            int numberOfEnergies = eVTest[eVTest.size() - 1];

            for (int i = 0; i < numberOfEnergies; i++)
            {
                fitValue += (pow(eVTest[i] - eVControl[i], 2)) / (numberOfEnergies);
            }
        }
    }
    return fitValue;
}

//***************//
// CFIT Function //
//***************//

// The cFit function produces a value that will be used
// as a measue of charge fitness of a parameter set

float GeneticAlgorithm::cFit(std::vector<float> eVTest,
                             std::vector<float> eVControl)
{
    float fitValue = 0.0;
    float chargeEnergyWeight = MH->IP.cew;
    if (chargeEnergyWeight > 1.0E-5)
    {
        int numberOfEnergies = eVTest[eVTest.size() - 1];
        for (int i = numberOfEnergies; i < (int)eVTest.size(); i++)
        {
            fitValue += (chargeEnergyWeight * pow(eVTest[i] - eVControl[i], 2)) /
                        (float)(eVTest.size() - numberOfEnergies - 1);
        }
    }
    return fitValue;
}

//***************//
// FFIT FUNCTION //
//***************//

// Returns a fitness measurment of the forces
float GeneticAlgorithm::fFit(std::vector<float> fTest,
                             std::vector<float> forceControl)
{
    float value = 0;
    float forceEnergyWeight = MH->IP.few;
    if (forceEnergyWeight > 1.0E-5)
    {
        for (unsigned int i = 0; i < fTest.size(); ++i)
        {
            value += (forceEnergyWeight * pow((fTest[i] - forceControl[i]), 2)) /
                     (float)fTest.size();
        }
    }
    return value;
}

//********************//
// SURVIVORS FUNCTION //
//********************//

// The survivors function chooses the survivors from an input list
// of the guassian outputs
std::vector<int> GeneticAlgorithm::survivors(vector<float>& fitvals,
        vector<float>& evfitvals,
        vector<float>& cfitvals,
        vector<float>& ffitvals, int b)
{
    int frcbst = MH->IP.forcebest;
    int NS = b;

    std::vector<int> survivorsn;

    std::vector<float> fitness;
    int G = MH->IP.G;
    float forceEnergyWeight = MH->IP.few;
    float chargeEnergyWeight = MH->IP.cew;

    // This for statment loops over all of the parameters sets (members of the
    // population)
    for (unsigned int parameterSet = 0; parameterSet < outputfiles.size();
            ++parameterSet)
    {
        float fitValue = 0;
        float eVFitValue = 0;
        float fFitValue = 0;
        float cFitValue = 0;

        // This if statment checks if the parameter set is blacklisted
        /*
        Blacklisting occurs when a parameters set causes gaussian to fail
        to converge and therefore fitness values are unobtainable.
        */
        if (!BL[parameterSet])
        {
            // Allocate per geometry data storage
            vector<float> GeomeVFit;
            GeomeVFit.resize(G);
            vector<float> GeomfeFit;
            GeomfeFit.resize(G);
            vector<float> GeomceFit;
            GeomceFit.resize(G);
            vector<float> GeomFit;
            GeomFit.resize(G);

            //*************************************************************************
            // This statement loops over all the geometries for the given parameter
            // set
            //*************************************************************************
            /*
            This happens for all non-blacklisted parameters sets for each geometry
            */
            for (int geometry = 0; geometry < G; geometry++)
            {
                // Calculate the eVFitValues
                GeomeVFit[geometry] =
                    eVFit(excitedStateEnergies(outputfiles[parameterSet][geometry]),
                          eVControl[geometry]);

                // If the force energy weight is set, calculate the forceFitValue
                if (forceEnergyWeight > 1.0E-5)
                {
                    GeomfeFit[geometry] = fFit(
                                              fReader(outputfiles[parameterSet][geometry]), fControl[geometry]);
                }
                else
                {
                    GeomfeFit[geometry] = 0;
                }

                // If the charge energy weight is set, calculate the chargeFitValue
                if (chargeEnergyWeight > 1.0E-5)
                {
                    GeomceFit[geometry] =
                        cFit(excitedStateEnergies(outputfiles[parameterSet][geometry]),
                             eVControl[geometry]);
                }
                else
                {
                    GeomceFit[geometry] = 0;
                }

                // Sum all values together
                eVFitValue += GeomeVFit[geometry];
                fFitValue += GeomfeFit[geometry];
                cFitValue += GeomceFit[geometry];

                GeomFit[geometry] =
                    GeomeVFit[geometry] + GeomfeFit[geometry] + GeomceFit[geometry];
            }

            // Combine into total fitness value
            fitValue = eVFitValue + fFitValue + cFitValue;

            //***********************************************
            // Determine if calculated fitness is a new Best
            //***********************************************
            if (fitValue < bestfitness && fitValue > 1.0E-7)
            {
                // Store Best Values
                bestfitness = fitValue;
                besteVFitness = eVFitValue;
                bestForceFitness = fFitValue;
                bestChargeFitness = cFitValue;

                // Print total new best Data for AM1->HL
                out->ofile << "-------------------------------------\n";
                out->ofile << "New Best (AM1-OPT -> HL): " << fitValue << endl;
                ;
                out->ofile << "    Geom #: {   Fit   ,  eVFit  ,  fFit   ,  cFit   }\n";

                for (int geometry = 0; geometry < G; geometry++)
                {
                    out->ofile << "    Geom " << geometry << ": {";
                    out->ofile << GeomFit[geometry] << ",";
                    out->ofile << GeomeVFit[geometry] << ",";
                    out->ofile << GeomfeFit[geometry] << ",";
                    out->ofile << GeomceFit[geometry] << "}" << endl;
                }

                out->ofile << "     Total: {";
                out->ofile << fitValue << ",";
                out->ofile << eVFitValue << ",";
                out->ofile << fFitValue << ",";
                out->ofile << cFitValue << "}" << endl;

                stringstream command;
                command << "cp " << outputfiles[parameterSet][0] << " BESTOUTPUT.log";

                string chk=exec(command.str());
                //system(command.str().c_str());

                //**************************************
                // Calculate the fits to the Original AM1
                //**************************************
                // Calculate the excited state fit between this and the original
                // Moved it here, so it would always be called
                vector<float> GeomAM1eVFit;
                GeomAM1eVFit.resize(G);
                vector<float> GeomAM1feFit;
                GeomAM1feFit.resize(G);
                vector<float> GeomAM1ceFit;
                GeomAM1ceFit.resize(G);
                vector<float> GeomAM1Fit;
                GeomAM1Fit.resize(G);

                orgam1fit = 0;
                orgam1ev = 0;
                orgam1force = 0;
                orgam1charge = 0;

                // Loop over geometries
                for (int geometry = 0; geometry < G; geometry++)
                {
                    // Calculate the eVFitValues
                    GeomAM1eVFit[geometry] =
                        eVFit(excitedStateEnergies(outputfiles[parameterSet][geometry]),
                              eVAM1Original[geometry]);

                    // If the force energy weight is set, calculate the forceFitValue
                    if (forceEnergyWeight > 1.0E-5)
                    {
                        GeomAM1feFit[geometry] =
                            fFit(fReader(outputfiles[parameterSet][geometry]),
                                 fAM1Original[geometry]);
                    }
                    else
                    {
                        GeomAM1feFit[geometry] = 0;
                    }

                    // If the charge energy weight is set, calculate the chargeFitValue
                    // for AM1
                    if (chargeEnergyWeight > 1.0E-5)
                    {
                        GeomAM1ceFit[geometry] =
                            cFit(excitedStateEnergies(outputfiles[parameterSet][geometry]),
                                 eVAM1Original[geometry]);
                    }
                    else
                    {
                        GeomAM1ceFit[geometry] = 0;
                    }

                    // Sum all values together
                    orgam1ev += GeomAM1eVFit[geometry];
                    orgam1force += GeomAM1feFit[geometry];
                    orgam1charge += GeomAM1ceFit[geometry];

                    GeomAM1Fit[geometry] = GeomAM1eVFit[geometry] +
                                           GeomAM1feFit[geometry] +
                                           GeomAM1ceFit[geometry];
                }

                // Calculate total fitness
                orgam1fit = orgam1ev + orgam1force + orgam1charge;

                // Print all of the optimized to original AM1 data
                out->ofile << "         (AM1-OPT -> AM1-ORG): " << orgam1fit << endl;
                ;
                out->ofile << "    Geom #: {   Fit   ,  eVFit  ,  fFit   ,  cFit   }\n";

                for (int geometry = 0; geometry < G; geometry++)
                {
                    out->ofile << "    Geom " << geometry << ": {";
                    out->ofile << GeomAM1Fit[geometry] << ",";
                    out->ofile << GeomAM1eVFit[geometry] << ",";
                    out->ofile << GeomAM1feFit[geometry] << ",";
                    out->ofile << GeomAM1ceFit[geometry] << "}" << endl;
                }

                out->ofile << "     Total: {";
                out->ofile << orgam1fit << ",";
                out->ofile << orgam1ev << ",";
                out->ofile << orgam1force << ",";
                out->ofile << orgam1charge << "}" << endl;
                out->ofile << endl;

                // Forces inclusion of new bests in reproduction
                if (frcbst == 1)
                {
                    survivorsn.push_back(parameterSet);
                    fitvals.push_back(fitValue);
                    evfitvals.push_back(eVFitValue);
                    cfitvals.push_back(cFitValue);
                    ffitvals.push_back(fFitValue);
                    --NS;
                }
            }

            if (fitValue < 1.0E-5)
            {
                // If a value ~0 is found, this is a gaussian failure.
                // However, this should not happen anymore because of the
                // blacklisting
                fitValue = 1.0E6;
            }
        }
        else
        {
            // If blacklisted, the parameter set will not be chosen
            fitValue = 0;
        }

        fitness.push_back(fitValue);
        evfitvals.push_back(eVFitValue);
        cfitvals.push_back(cFitValue);
        ffitvals.push_back(fFitValue);
    }

    // Some statistics calculations
    mu = mean(fitness, NBL);  // Total population mean
    std = stddev(fitness, mu, NBL);

    out->ofile << "\nWhole Population Fitness Mean(Std. Dev.): " << mu << "("
               << std << ")\n\n";

    // Begin finding survivors
    for (int i = 0; i < NS; i++)
    {
        // Normalize
        float normalizer = 0;
        for (unsigned int j = 0; j < fitness.size(); j++)
        {
            if (!BL[j])  // Blacklist failed parameter sets
            {
                /* I modified this to fix the indexing problem stated below.
                There is probably a better way to do this, but it works,
                so think about it. -Justin*/
                bool found = false;
                for (unsigned int m = 0; m < survivorsn.size(); ++m)
                {
                    if ((unsigned int)survivorsn[m] == j)
                    {
                        found = true;
                    }
                }

                if (!found)
                {
                    normalizer += pow(1 / (fitness[j]), 2);
                }
            }
        }

        // Create Rank Vector
        float rankCount = 1.0;
        std::vector<float> rankv;
        for (unsigned int k = 0; k < fitness.size(); k++)
        {
            bool found = false;
            for (unsigned int m = 0; m < survivorsn.size(); ++m)
            {
                if ((unsigned int)survivorsn[m] == k)
                {
                    found = true;
                }
            }

            // Checks 1) Has genome already been selected and 2) Has genome failed
            // If either is true then set chance of choosing to zero
            if (!found && !BL[k])
            {
                float chance = pow(1 / fitness[k], 2) / normalizer;
                rankCount -= chance;
                rankv.push_back(rankCount);
            }
            else
            {
                float chance = 0.0;
                rankCount -= chance;
                rankv.push_back(rankCount);
            }
        }

        // Create random number and choose which gaussian input survives
        UniformRandomReal UR(1, clock());
        float selector = UR.UniformRandReal(0.0, 1.0);
        int l = 0;

        while (selector < rankv[l])
        {
            l++;
        }

        survivorsn.push_back(l);

        fitvals.push_back(fitness[l]);
    }
    return survivorsn;
};
