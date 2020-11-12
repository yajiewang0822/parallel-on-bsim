#include <vector>
using namespace std;

vector<vector<double> > fastConvolution(vector<vector<double> > obj, vector<vector<double> > filter);
vector<vector<double> > matrixDotMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixAdd(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixSub(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixAbs(vector<vector<double> > matrix);
vector<vector<double> > matrixScalarMul(vector<vector<double> > matrix, double multiplier);
double sumImage(vector<vector<double> > input, int img_size);