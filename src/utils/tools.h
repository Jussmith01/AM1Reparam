#ifndef tools_header
#define tools_header

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

using namespace std;

/*----------------------------------------
       ***************************
       |  Trim String Function   |
       ***************************
Trim whitespace from a string
------------------------------------------*/
inline extern string trim (string line)
{
        string trimmed_line;
        unsigned int first = line.find_first_not_of(" \t\r\n\x0b");
        unsigned int last = line.find_last_not_of(" \t\r\n\x0b") + 1;


        last = last - first;

        //Check for string of blanks and set default to blank
        int length = (int)line.length();
        int cntr = 0;
        for (int i=0;i<length;++i)
        {
                string test = line.substr(i,1);
                if(test.find_first_of(" \t\r\n\x0b")==0)
                {++cntr;}
        }

        //Trim spaces from front and end of line
        if (line.length() == 0 || (length - cntr) == 0)
        {
                trimmed_line = "";
        } else {
                trimmed_line = line.substr(first,last);
        }

        return trimmed_line;
};

/*----------------------------------------
       ***************************
       |    Command Execution    |
       ***************************
------------------------------------------*/
inline extern int CountFails(vector<bool> &arr)
{
	int cnt=0;

	for(int i=0;i<(int)arr.size();++i)
	{
		if (arr[i]) {++cnt;}
	}

	return cnt;
};
/*----------------------------------------
       ***************************
       |    Command Execution    |
       ***************************
------------------------------------------*/
inline extern string exec(string cmd) {
    //cout << "EXEC: " << cmd << endl;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "ERROR";
    char buffer[4000];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 4000, pipe) != NULL)
                result += buffer;
    }
    pclose(pipe);
    return result;
};

/*----------------------------------------
       ***************************
       | G09 Command Execution   |
       ***************************
------------------------------------------*/
inline extern bool execg09(string cmd,string wkdir,int tid) {
	//cout << "EXEC: " << cmd << endl;
	stringstream sscmd,errout;
	
	errout << wkdir << "gauerr_tid" << tid << ".e";
	sscmd << cmd << " 2> " << errout.str();
	
	// Open a pipe and run command
	FILE* pipe = popen(sscmd.str().c_str(), "r");
	if (!pipe) return "ERROR";
	char buffer[100];
	std::string result = "";
	while(!feof(pipe)) {
	    if(fgets(buffer, 100, pipe) != NULL)
	            result += buffer;
	}
	pclose(pipe);

	//Check for gaussian error
	string line;
	bool gerr = false;
	ifstream errfile (errout.str().c_str());
	if (errfile.is_open())
	{
		getline (errfile,line);
		//cout << "FIRSTLINE: " << trim(line) << endl;
		if (trim(line).compare("Error: segmentation violation")==0) 
			{gerr=true;}
		errfile.close();
	}
	else cout << "**Unable to open gaussian error file: " << errout.str() << " **"<<endl; 

	stringstream rmcmd;
	rmcmd << "rm " << errout.str();

	// Open a pipe and run command
        FILE* pipe2 = popen(rmcmd.str().c_str(), "r");
        if (!pipe2) return "ERROR";
        char buffer2[100];
        std::string result2 = "";
        while(!feof(pipe2)) {
            if(fgets(buffer2, 100, pipe2) != NULL)
                    result2 += buffer2;
        }
        pclose(pipe2);

    	return gerr;
};

/*----------------------------------------
       ***************************
       |      Mean Function      |
       ***************************
Calculate the arithmetic mean of a set of
values.
------------------------------------------*/
template <class T>
inline extern T mean(vector<T> &values,int err)
{
	T SUM=0;

	unsigned int cnt = values.size();
	unsigned int actcnt = cnt;

	for (unsigned int i=0;i<cnt;++i)
	{
		float val = values[i];

		if (!isinf(val))
			{SUM+=val;} 
		else 
			{--actcnt;}
	}

	return SUM/(T)(actcnt-err);
};

/*----------------------------------------
       ***************************
       |   Std. Dev. Function    |
       ***************************
Calculate the standard deviation of a set of
values given a mean, mu.
------------------------------------------*/
template <class T>
inline extern T stddev(vector<T> &values, T mu,int err)
{
        T SUM2=0;
	unsigned int cnt = values.size();
        unsigned int actcnt = cnt;

        for (unsigned int i=0;i<cnt;++i)
	{
		
                float val = values[i];

                if (!isinf(val))
		{
			T diff = val - mu;
                	SUM2+=diff*diff;
		} else 
			{--actcnt;}
	}

        return sqrt(SUM2/(T)(actcnt-err));
};

/*----------------------------------------
       ***************************
       |  Comma Parsing Funtion  |
       ***************************
Separate a list of things delimited by commas.
------------------------------------------*/
inline extern vector<string> commaparse (string line)
{
        vector<string> vals;
        size_t pos = line.find_first_of(",");

	if (pos==string::npos)
		{vals.push_back(line);}
	else 
	{
		while(true)
		{
        		pos = line.find_first_of(",");
			
			if (pos!=string::npos)
			{
				vals.push_back(line.substr(0,pos));
				line = line.substr(pos+1);
			} else {
				vals.push_back(line);
				break;
			}
		}
	}
        return vals;
};

/*----------------------------------------
  ***************************************
  |  Obtain first two Chars of string   |
  ***************************************
Used with the AM1 parameter loader function
in the AM1 parameter class.
------------------------------------------*/
inline extern string ObtainStringLetter(string &line)
{
        string rtnstr = trim(line.substr(0,2));
        line = trim(line.substr(2));
        return rtnstr;
};

/*----------------------------------------
     ***************************
     |     Parse a String      |
     ***************************
Parse a string 'line' with the format below:
line: 'SomeFlag=xx.xxx'

Returns in the following referened variable: 
flag=SomeFlag and value=xx.xxx
------------------------------------------*/
inline extern void ParseString(string &flag,string &value,string &line)
{
        size_t d1 = line.find_first_of("=");
        size_t d2 = line.find_first_of(" ");

        flag = trim(line.substr(0,d1));
        value = trim(line.substr(d1+1,d2-d1-1));
        line = trim(line.substr(d2));
};

#endif

