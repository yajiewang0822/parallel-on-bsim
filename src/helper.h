#include <vector>
#include <complex>
#include <fftw3.h>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <string>

using namespace std;
#define NA_SPEC 0.9
#define PATTERN_NUM 384
#define IMG_SIZE 256
#define PIXEL_SIZE 30

vector<vector<double> > fconv2(vector<vector<double> > obj, vector<vector<double> > filter);
vector<vector<double> > deconvolution(vector<vector<double> > obj, vector<vector<double> > filter);
void fft2(vector<vector<double> > input, fftw_complex *output);
vector<vector<double > > ifft2(fftw_complex *input);
vector<vector<complex<double> > > getValueOfComplex(vector<vector<complex<double> > > input);
vector<vector<double> > matrixDotMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixAdd(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixSub(vector<vector<double> > matrix1, vector<vector<double> > matrix2);
vector<vector<double> > matrixAbs(vector<vector<double> > matrix);
vector<vector<double> > matrixScalarMul(vector<vector<double> > matrix, double multiplier);
vector<vector<double> > matrixColMean(vector<vector<double> > matrix);
double sumImage(vector<vector<double> > input);
void saveImage(vector<vector<double> > input, string filename);