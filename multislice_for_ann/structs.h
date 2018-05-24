#include <string>
#include <vector>
#include <complex>

using namespace std;

struct bools{
	bool isSC;
	bool isPPR;
	bool isVP;
	bool isEP;
	bool isMS;
	bool isTR;
	bool isTO;
	bool isOM;
	bool isRB;
	bool isTV;
	bool isPE;
	bool isNoiseLoop;
	bool isDefocLoop;
	bool isAt;
	bool isGA;
	bool isCF;
	bool isGFBP;
	bool isAIFS;

};

struct boolsVis{
	bool isMicro;
	bool isProjPhases;
	bool isRetPhases;
	bool isMagPhases;
	bool isVPrec;
	bool isVPex;
	bool isImagesOut;

};


struct imagingArgs{
	int mode;
	float k2max;
	float objlx;
	float objly;
	float Cs3;
	float Cs5;
};

struct VEC{
	double x;
	double y;
	double z;
};

struct rockingBeamArgs{
	double maxIll;
	int sampling;
};

struct ARGS{
	string material;
	double magnetisationFactor;
	double absorptionFactor;
	string TOgeo;
	rockingBeamArgs rockingBeam;
	double accelPot;
	double domain;
	double deltaz;
	double noise;
	double window;
	VEC defocus;
	string samplePath_std;
	string outputPath_std;
	string inFocusMethod;
	string momentDistributionMethod;
	bools guibools;
	boolsVis visualisation;
	imagingArgs imaging;
	size_t defocusSteps;
	double temperature;
	double nwobble;
	size_t nt;
	size_t naa;
	size_t M;
	size_t Mhr;

	VEC pyr;
	VEC m;
	VEC t;
	VEC l;
	VEC alpha;
	VEC noiseLoop;
	VEC defocLoop;

	string method;
	int gammaTot;
	int gammaPrimeTot;
};

struct VEC3{
	vector<vector<vector<complex<double> > > > x;
	vector<vector<vector<complex<double> > > > y;
	vector<vector<vector<complex<double> > > > z;
};




struct MATERIAL{
	double massMagnetisation;//<------emu/g
	double density;//<--------------g/cm^3             Permalloy(~8.74), Magnetite(~5.18)
	complex<double> mip;//<------------Mean inner potential in Volts.         Permalloy(~26), Magnetite(~17)
	double unitCellSize[3];//<--------Angstrom
	double m0;//<------A/m
};



struct ORIENT{
	string dir;
	string ts;
	size_t tilt;
	double angle;
};

struct TFARGS{
	ARGS *pArgs;
	VEC3 Vp;
	const vector<vector<float> > &atomicLocations;
	ORIENT orientation;
	vector<vector<vector<double> > > hires;
	const vector<vector<vector<double> > > &shapeFunction;
	const MATERIAL material;
	string &outputPath;
	bool isArbMag;
	vector<vector<vector<vector<double> > > > momentDistribution;






	//Initialisation Constructor
	TFARGS(
		vector<vector<float> > &atomicLocations_,
		vector<vector<vector<double> > > &shapeFunction_,
		string &outputPath_) :
		atomicLocations(atomicLocations_),
		shapeFunction(shapeFunction_),
		outputPath(outputPath_) {
	}
			




};




