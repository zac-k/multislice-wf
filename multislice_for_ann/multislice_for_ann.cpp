// multislice_for_ann.cpp : Defines the entry point for the console application.
//



#include "other.cpp"





int main()
{
	cout << "Start number: ";
	int first_example;
	cin >> first_example;
	cout << "End number: ";
	int final_example;
	cin >> final_example;
	if (final_example < first_example) {
		cout << "End number must not be less than Start number. Press any key to exit.";
		system("pause");
		return 0;
	}
	string sample_path = ".\\training4\\";//"C:\\Users\\zac\\PycharmProjects\\phaseRetDL\\data\\specimens\\training4\\";
	string output_path = ".\\images\\";// C:\\Users\\zac\\PycharmProjects\\phaseRetDL\\data\\images\\multislice\\";
	int Mhr = 1024;
	int M = 1024;
	int spec_size = 64; // Size of specimen in pixels
	double a = 150e-9;
	double defoc = 10e-6;
	double deltaz = 5; // slice thickness in angstroms

	ARGS pArgs;
	pArgs.domain = a;
	VEC defocus;
	defocus.x = defoc;
	defocus.y = 0.0;
	defocus.z = 0.0;
	pArgs.defocus = defocus;
	pArgs.accelPot = 300e3;
	pArgs.deltaz = deltaz;
	
	pArgs.imaging.Cs3 = 0.5; // spherical aberration in mm
	pArgs.imaging.Cs5 = 0.0;
	pArgs.imaging.k2max = 1000; // in mrad

	pArgs.guibools.isTV = false;
	pArgs.guibools.isVP = false;

	/* random number generator */
	uniform_real_distribution<double> dist(0.0, 2 * PI);
	mt19937_64 rng;
	rng.seed(random_device{}());


	

	atomicLocationsMutex = CreateMutex(NULL, FALSE, NULL);
	if (atomicLocationsMutex == NULL)  cout << "CreateMutex error";
	TSmutex = CreateMutex(NULL, FALSE, NULL);
	if (TSmutex == NULL)  cout << "CreateMutex error";
	hiresMutex = CreateMutex(NULL, FALSE, NULL);
	if (hiresMutex == NULL)  cout << "CreateMutex error";
	shapeFunctionMutex = CreateMutex(NULL, FALSE, NULL);
	if (shapeFunctionMutex == NULL)  cout << "CreateMutex error";
	freeMutex = CreateMutex(NULL, FALSE, NULL);
	if (freeMutex == NULL)  cout << "CreateMutex error";

	for (size_t example = first_example; example < final_example; example++) {
		vector<vector<vector<bool> > > shapeFunctionHR = CreateWF(spec_size, spec_size, spec_size, false);
		string shapeFunctionFile = sample_path + "particle(" + to_string(example) + ")";
		cout << "Read specimen mask file..." << endl;
		ReadWF(shapeFunctionHR, shapeFunctionFile);
		pArgs.pyr.x = dist(rng);
		pArgs.pyr.y = dist(rng);
		pArgs.pyr.z = dist(rng);

		vector<vector<float> > atomicLocations(1, vector<float>(4));
		apCalculation(atomicLocations, a, shapeFunctionHR);

		{
			cout << "Generating sample " << to_string(example) << endl;
			
			vector<vector<vector<complex<double> > > > Vpz = CreateWF(M, M, M, complex<double>(0.0, 0.0));

			pArgs.Mhr = M;
			pArgs.deltaz = 0;

			vector<vector<complex<double>>> wavefield = CreateWF(M, M, complex<double>(0.0, 0.0));
			autoslic(wavefield, Vpz, atomicLocations, &pArgs);

			WriteWF(wavefield, output_path + "image(" + to_string(example) + ")");
		}

		{
			pArgs.Mhr = Mhr;
			pArgs.deltaz = 5;
			vector<vector<vector<complex<double> > > > Vpz = CreateWF(M, M, M, complex<double>(0.0, 0.0));
			
			vector<vector<complex<double>>> wavefield_hr = CreateWF(Mhr, Mhr, complex<double>(0.0, 0.0));
			autoslic(wavefield_hr, Vpz, atomicLocations, &pArgs);

			WriteWF(wavefield_hr, output_path + "image_hr(" + to_string(example) + ")");
		}
		

	}
    return 0;
}









