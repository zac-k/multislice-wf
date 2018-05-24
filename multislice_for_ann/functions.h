#pragma once



using std::setprecision;
using std::setw;



/*
Set colour of console text
*/
void setcolor(int fore, int back);

/*
Create and initialize a 2D real vector array
*/
vector<vector<double> > CreateWF(const size_t nx, const size_t ny, const double init_val);

void flipImage(vector<vector<double> > &Image, string axis);

/*
Create a 2D complex vector array
*/
vector<vector<complex<double> > > CreateWF(const size_t nx, const size_t ny, const complex<double> init_val);

/*
Create and initialize a 3D boolean vector array
*/
vector<vector<vector<bool> > > CreateWF(const size_t n1, const size_t n2, const size_t n3, const bool init_val);

/*
Create a 3D complex vector array
*/
vector<vector<vector<complex<double> > > > CreateWF(const size_t n1, const size_t n2, const size_t n3,
	const complex<double> init_val);


/*create a string from a number with a set number of decimal places (precision). This is
used for creating folders and files with the names of numeric variables*/
string stringFromNumber(const double &numberInput, const int &precision);

bool fexists(const char *filename);

/* Rotate vector 'x' around the 'rotDir' axis by an angle 'theta' */
void  rotation(const double x[3], const double rotDir[3], const double theta, double output[]);

void apCalculation(vector<vector<float> > &loc1, double a, const vector<vector<vector<bool> > > &shapeFunction1);

/*
1D Fourier transform, non-centred data
*/
void fft(vector<complex<double> > &vec, int sign);

/*
1D Fourier transform, centred data
*/
void fftEx(vector<complex<double> > &vec, int sign);

/*
2D Fourier transform, centred data
*/
void fftEx(vector<vector<complex<double> > > &vec2d, int sign);

/*
3D Fourier transform, centred data
*/
void fftEx(vector<vector<vector<complex<double> > > > &vec3d, int sign);


/*Low pass filter a square 2D array with a gaussian with SD sigma*/
void lowpassFilter(vector<vector<complex<double> > > &vec, const double &sigma);

