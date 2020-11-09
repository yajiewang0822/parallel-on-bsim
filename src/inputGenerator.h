#include<vector>
using namespace std;

class InputGenerator
{
private:
  vector<vector<double>> psf;
  vector<vector<double>> psfn;
  double NA;
  double lambda;
  int pattern_num; 
  int img_size;
  void GeneratePSF();
  void GenerateOTF();
  vector<vector<double>> GenerateObjective();
  vector<vector<vector<double>>> GeneratePatterns();

public:
  InputGenerator(double NA, double lambda, int pattern_num, int img_size);
  vector<vector<double>> getPSF();
  vector<vector<double>> getPSFn();
  void GenerateInputs();
};