/**
 * This is the header file containing all the arithematic functions declarations needed to run the simulation. 
 * Those funsitons include convlolution, FFT, and matrix operations. 
 * It also involves file IO to store the image. 
 * 
 * @author: Peicheng Tang 
 * @author: Yajie Wang
 */

#include <vector>
#include <complex>
#include <fftw3.h>
#include <opencv2/core/core.hpp>
//#include <opencv2/imgcodes/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <string>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
#define NA_SPEC 0.9
#define PATTERN_NUM 2048
#define IMG_SIZE 512
#define PIXEL_SIZE 30

/**
 * FFTSHIFT and IFFTSHIFT operation
 */ 
//fftshift the fftw_complex array based on circular shift
#define fftshift(out, in, x, y) circshift(out, in, x, y, (x/2), (y/2))
//ifftshift the fftw_complex array based on circular shift
#define ifftshift(out, in, x, y) circshift(out, in, x, y, ((x+1)/2), ((y+1)/2))
//fftshift the double array based on circular shift
#define fftshift_double(out, in, x, y) circshift_double(out, in, x, y, (x/2), (y/2))
//ifftshift the double array based on circular shift
#define ifftshift_double(out, in, x, y) circshift_double(out, in, x, y, ((x+1)/2), ((y+1)/2))

vector<double> fconv2(vector<double>obj, vector<double> filter);
void fft2(vector<double> input, fftw_complex *output);
vector<double> ifft2(fftw_complex *input);
void circshift(fftw_complex *out, fftw_complex *in, int xdim, int ydim, int xshift, int yshift);
void circshift_double(double *out, double *in, int xdim, int ydim, int xshift, int yshift);

double sumImage(vector<double> input);
void saveImage(vector<double> input, string filename);
void saveData(vector<double> input, string filename);
vector<double> readData(string filename);

/**
 * matrix operation
 */ 
vector<double> matrixEleMul(vector<double> matrix1, vector<double> matrix2);
vector<double> matrixAdd(vector<double> matrix1, vector<double> matrix2);
vector<double> matrixSub(vector<double> matrix1, vector<double> matrix2);
vector<double> matrixAbs(vector<double> matrix);
vector<double> matrixScalarMul(vector<double> matrix, double multiplier);
vector<double> matrixColMean(vector<double> matrix);
