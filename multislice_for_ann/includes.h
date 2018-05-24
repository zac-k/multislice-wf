	
#define _CRT_RAND_S

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <string>
using std::string;
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <complex>
using std::complex;
#include <iomanip>
using std::setprecision;
using std::setw;
#include <FLOAT.H>
#include <fstream>
using std::ofstream;
using std::ifstream;


#include "stdafx.h"
#include "consts.h"						/* Definition of physical and mathematical constants */ 
#include "structs.h"
//#include "./prvt.h"					/* Headers for pvrt.cpp (and others) */


#include <io.h>   // For access().
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().
#include <iostream>
#include "fftw3.h"      /* FFT routines from FFTW 3*/






#include <windows.h>
HANDLE atomicLocationsMutex;
HANDLE TSmutex;
HANDLE hiresMutex;
HANDLE shapeFunctionMutex;
HANDLE freeMutex;


//#include <msclr\marshal_cppstd.h>
//#include "Form1.h"						/* Windows form containing GUI controls */
//#include "./prvt.cpp"					/* Phase retrieval and Vector Tomography functions */
//#include "../phaseError.cpp"

#include <direct.h>

#include <atlstr.h>
#include <sstream>
#include <random>


#include <conio.h>




//#using <System.dll>
//#using <System.Drawing.dll>
//#using <System.Windows.Forms.dll>

using std::setprecision;
using std::setw;
using namespace std;
//using namespace My1TransferProjectedWF;
#include <cstdio>  /* ANSI C libraries */



#include "functions.h"

#include "slicelib.hpp"   /* misc. routines for multislice */
#include "autoslic.cpp"   /* Multislice exit-wave calculation*/