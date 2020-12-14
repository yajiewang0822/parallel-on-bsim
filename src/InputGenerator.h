/**
 * This is the header file of our inpute generator. 
 * 
 * @author: Peicheng Tang 
 * @author: Yajie Wang
 */
#include <vector>
#include <limits>
#include <random>
#include <algorithm>

using namespace std;

#define PI 3.14159265
#define NUM_SPECKLE 512
#define NA 0.9
#define LAMBDA 520
#define REPEAT 4
const double eps=numeric_limits<double>::epsilon();

enum PSF_TYPE{
  PSF, PSFN
};

class InputGenerator
{
private:
  vector<double> psf;
  vector<double> psfn;
  double NA_spec;
  int pattern_num; 
  int p_size;
  void GeneratePSF(double effect_NA, PSF_TYPE type);
  vector<double> GenerateObjective();
  vector<vector<double> > GeneratePatterns();

public:
  InputGenerator(double NA_spec, int pattern_num, int p_size);
  vector<vector<double> > GenerateInputs();
  vector<double> getPSF();
  vector<double> getPSFn();
};

