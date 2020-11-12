#include <vector>
using namespace std;

#define ITER_NUM 150

class BSIM
{
private:
  int pat_num;
  int img_size;
  vector<vector<double> > covariance(vector<vector<double> > input, vector<vector<double> > pattern);
public:
  BSIM(int pat_num, int img_size);
  vector<vector<double> > Reconstruction(vector<vector<vector<double> > > inputs, vector<vector<double> > psfn, vector<vector<double> > psf);
};