#include <vector>
using namespace std;

class BSIM
{
private:
  int pat_num;
public:
  BSIM(int pat_num);
  vector<vector<double> >Reconstruction(vector<vector<vector<double> > > inputs, vector<vector<double> > psfn, vector<vector<double> > psf);
};