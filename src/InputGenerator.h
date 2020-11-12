#include <vector>
using namespace std;

#define PI 3.14159265
#define NUM_SPECKLE 512
#define NA 0.9
#define LAMBDA 520
class InputGenerator
{
private:
  vector<vector<double> > psf;
  vector<vector<double> > psfn;
  double NA_spec;
  int pattern_num; 
  int img_size;
  int p_size;
  void GeneratePSFandOTF(double effect_NA);
  vector<vector<double> > GenerateObjective();
  vector<vector<vector<double> > > GeneratePatterns();

public:
  InputGenerator(double NA_spec, int pattern_num, int img_size, int p_size);
  vector<vector<vector<double> > > GenerateInputs();
  vector<vector<double> > getPSF();
  vector<vector<double> > getPSFn();
  double sumImage(vector<vector<double> > input);
  vector<vector<double> > fastConvolution(vector<vector<double> > obj, vector<vector<double> >filter);
};

