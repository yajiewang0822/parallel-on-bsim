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
  vector<double> Reconstruction(vector<vector<double> > inputs, vector<double> psfn, vector<double> psf);
};