#include <vector>
#include <limits>
#include <random>
#include <algorithm>

using namespace std;

#define PI 3.14159265
#define NUM_SPECKLE 512
#define NA 0.9
#define LAMBDA 520
const double eps=numeric_limits<double>::epsilon();

enum PSF_TYPE{
  PSF, PSFN
};

class InputGenerator
{
private:
  vector<vector<double> > psf;
  vector<vector<double> > psfn;
  double NA_spec;
  int pattern_num; 
  int p_size;
  void GeneratePSF(double effect_NA, PSF_TYPE type);
  vector<vector<double> > GenerateObjective();
  vector<vector<vector<double> > > GeneratePatterns();

public:
  InputGenerator(double NA_spec, int pattern_num, int p_size);
  vector<vector<vector<double> > > GenerateInputs();
  vector<vector<double> > getPSF();
  vector<vector<double> > getPSFn();
};

