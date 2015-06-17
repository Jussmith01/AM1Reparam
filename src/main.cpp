//**C/C++ included libraries**
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include "utils/timer_class.h"
#include "utils/output_class.h"
#include "utils/flaghandler_class.h"
#include "utils/readinput_function.h"
#include "core/memhandler_class.h"
#include "core/genetic_algorithm.h"

int main (int argc, char *argv[])
{
        //********************************************//
        //     Initial Class/Struct Construction      //
        //********************************************//

	/*1) Set the flags from program execution and store in class*/
		flaghandler flags(argc,argv);
	/*2) Construct class to handle all output*/
		outputFile out(flags.outputfilename,flags.verbosity);
	/*3) Construct and set the program timer class*/
		timer ProgTimer;ProgTimer.set_timer();
    /*4) Construct memory handler class: Allocation occurs later*/
        MemoryHandler MH(&flags);
	/*5) Construct the Genetic Algorithm Class*/
		GeneticAlgorithm GA(&MH,&out);

	//********************************************//
        //          Pre-Execution Programs            //
        //********************************************//

	/*1) Read the input file, save program parameters*/
		ReadInput(flags.inputfilename,MH,out);
	/*2) Prepare the memory handler class*/
		MH.LoadAPH();//Loads the AM1 Params
	/*3) Build initial population and setup target fitness*/
		GA.GeneticAlgorithmPreparation();

        //********************************************//
        //             Execute Main Loop              //
        //********************************************//

	/*1) Launch the Genetic Algorithm Main Loop*/
		GA.GeneticAlgorithmMainloop();

        //********************************************//
        //          Post Execution Programs           //
        //********************************************//
        ProgTimer.end_timer();//End the program wall timer
        //ProgTimer.print_clock_time("Computation Clock Time: ",out.ofile);//Print the wall time to output
        ProgTimer.print_time("Computation Wall Time: ",out.ofile);//Print the wall time to output
	out.close_output(0);

	return 0;
}

