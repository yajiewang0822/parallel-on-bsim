#include <vector>
#include <complex>
using namespace std;

vector<vector<double> > fastConvolution(vector<vector<double> > obj, vector<vector<double> > filter);
vector<vector<double> > deconvolution(vector<vector<double> > obj, vector<vector<double> > filter);
vector<vector<complex<double> > > fft2(vector<vector<double> > input);
vector<vector<double> > ifft2(vector<vector<complex<double> > > input);
vector<vector<complex<double> > > complexValue(vector<vector<complex<double> > > input, int size);
vector<vector<double> > matrixDotMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixAdd(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixSub(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixAbs(vector<vector<double> > matrix);
vector<vector<double> > matrixScalarMul(vector<vector<double> > matrix, double multiplier);
vector<vector<double> > matrixColMean(vector<vector<double> > matrix);
double sumImage(vector<vector<double> > input, int img_size);