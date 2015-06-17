#ifndef variables_header
#define variables_header

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

using namespace std;

/*----------------------------------------
       ***************************
       |    String Variables     |
       ***************************
Possibly needed variables for string
manipulations.
------------------------------------------*/
namespace strvar {
	string spc=" ";
	string eq="=";

	/*AM1 Hamiltonian Names*/
	string AM1PVars[] =
	{
		/*0 */"PQN",
		/*1 */"NValence",
		/*2 */"F0ss",
		/*3 */"F0sp",
		/*4 */"F0pp",
		/*5 */"F0sd",
		/*6 */"F0pd",
		/*7 */"F0dd",
		/*8 */"F2pp",
		/*9 */"F2pd",
		/*10*/"F2dd",
		/*11*/"F4dd",
		/*12*/"G1sp",
		/*13*/"G1pd",
		/*14*/"G2sd",
		/*15*/"G3pd",
		/*16*/"Rsppd",
		/*17*/"Rsdpp",
		/*18*/"Rsddd",
		/*  */"****"
	};
};	
#endif

