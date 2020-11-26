#include <vector>
#include <complex>
#include <fftw3.h>

using namespace std;

#define PI 3.14159265
#define NUM_SPECKLE 512
#define NA 0.9
#define LAMBDA 520

enum PSF_TYPE{
  PSF, PSFN
};

class InputGenerator
{
private:
  vector<vector<double> > psf;
  vector<vector<double> > psfn;
  fftw_complex *OTF;
  fftw_complex *OTFn;
  double NA_spec;
  int pattern_num; 
  int p_size;
  void GeneratePSFandOTF(double effect_NA, PSF_TYPE type);
  vector<vector<double> > GenerateObjective();
  vector<vector<vector<double> > > GeneratePatterns();

public:
  InputGenerator(double NA_spec, int pattern_num, int p_size);
  vector<vector<vector<double> > > GenerateInputs();
  vector<vector<double> > getPSF();
  vector<vector<double> > getPSFn();
  fftw_complex *getOTF();
  fftw_complex *getOTFn();
};

