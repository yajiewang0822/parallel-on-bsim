/**
 * This is header file of BSIM 
 * 
 * @author: Peicheng Tang 
 * @author: Yajie Wang
 */

#include <vector>
using namespace std;

#define ITER_NUM 50

class BSIM
{
private:
  int pat_num;
public:
  BSIM(int pat_num);
  //TODO
  vector<vector<double> > patternEstimation(vector<vector<double> > inputs,  vector<double> psf, vector<double> obj, double coefficient, int procID);
  void Reconstruction(vector<vector<double> > inputs, vector<double> psf, int procID);
};