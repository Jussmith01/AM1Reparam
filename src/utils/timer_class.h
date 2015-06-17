#ifndef timer_header
#define timer_header

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>

using namespace std;

//________________________________________________________________________//
//      *************************************************************     //
//                               Timer Class
//                  Holds timer variables and class functions
//      *************************************************************     //
class timer
{
        //--------------------------
        //Private Class Declarations
        //--------------------------

        time_t start_time; //Holds start time
        time_t end_time; //Holds end time
        double run_time; //Holds run time = end_time - start_time
	clock_t t;
	double CLOCK_TIME;

        //------------------------------
        //Private Member Class Functions
        //------------------------------

	//Intakes a value of time in seconds and returns a string formmatted as:
	//Computation Wall Time: d days h hours m minutes s seconds
        string mk_time_string(double time_val); 

	public:
        //-----------------------------
        //Public Member Class Functions
        //-----------------------------

		//Sets the start time for the timer
        	void set_timer(void);

		//Sets the end time for the timer
        	void end_timer(void);

		//Prints run_timer = end_time - start_time 
        	void print_time(string message,ofstream &file1);
		
		//Prints the clock time
		void print_clock_time (string message,ofstream &file1);

                //Prints run_timer = end_time - start_time 
                void print_time_s(string message,ofstream &file1);

		//Resets the timer if needed
        	void reset(void);
};
#endif

