#include "helper.h"
#include "cmath"

using namespace std;
//TODO: with Halide/FFTW
vector<vector<double> > fconv2(vector<vector<double> > obj, vector<vector<double> > filter)
{
  fftw_complex *obj_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  fftw_complex *filter_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  fft2(obj, obj_fft);
  fft2(filter, filter_fft);
  fftw_complex *result_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  for(int i = 0; i < IMG_SIZE * IMG_SIZE; i++){
    result_fft[i][0] = obj_fft[i][0] * filter_fft[i][0] - obj_fft[i][1] * filter_fft[i][1];
    result_fft[i][1] = obj_fft[i][0] * filter_fft[i][1] + obj_fft[i][1] * filter_fft[i][0];
  }
  double sum = sumImage(filter);
  vector<vector<double> > result(IMG_SIZE, vector<double>(IMG_SIZE,0));
  result = ifft2(result_fft);
  printf("%f\n", sum);
  result = matrixScalarMul(result, (1.0/sum));
  return result;
}

//TODO: with Halide/FFTW, different from inversion
vector<vector<double> > deconvolution(vector<vector<double> > obj, vector<vector<double> > filter)
{
  vector<vector<double> > result;
  return result;
}

void fft2(vector<vector<double> > input, fftw_complex *output){
  fftw_complex *input_temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  for (int i=0;i<IMG_SIZE;i++){
    for (int j=0;j<IMG_SIZE;j++){
      input_temp[(i*IMG_SIZE)+j][0]=input[i][j];
      input_temp[(i*IMG_SIZE)+j][1]=0;
    }
  }
  fftw_plan p=fftw_plan_dft_2d(IMG_SIZE, IMG_SIZE, input_temp, output, FFTW_FORWARD, FFTW_MEASURE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_free(input_temp);
}

vector<vector<double > > ifft2(fftw_complex *input){
  fftw_complex *output_temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  vector<vector<double> > output(IMG_SIZE,vector<double>(IMG_SIZE,0));
  fftw_plan p=fftw_plan_dft_2d(IMG_SIZE, IMG_SIZE, input, output_temp, FFTW_BACKWARD, FFTW_MEASURE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  for (int i=0;i<IMG_SIZE;i++){
    for (int j=0;j<IMG_SIZE;j++){
      output[i][j] = output_temp[(i*IMG_SIZE)+j][0];
    }
  }
  fftw_free(output_temp);
  return output;
}

vector<vector<complex<double> > > getValueOfComplex(vector<vector<complex<double> > > input){
  for (int i = 0; i < IMG_SIZE; i++){
    for (int j = 0; j < IMG_SIZE; j++){
      input[i][j]= input[i][j] * conj(input[i][j]);
    }
  }
  return input;
}

//TODO: with Halide
vector<vector<double> > matrixDotMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
  int size = matrix1.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix1[i][j] * matrix2[i][j];
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixSub(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
  int size = matrix1.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix1[i][j] - matrix2[i][j];
    }
  }
  return result;
}

//TODO: with Halide
vector<vector<double> > matrixAdd(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
  int size = matrix1.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix1[i][j] + matrix2[i][j];
    }
  }
  return result;
}

//TODO: with Halide
vector<vector<double> > matrixScalarMul(vector<vector<double> > matrix, double multiplier)
{
  int size = matrix.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix[i][j] * multiplier;
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixAbs(vector<vector<double> > matrix)
{
  int size = matrix.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = abs(matrix[i][j]);
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixColMean(vector<vector<double> > matrix){
  int size = matrix.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int j = 0; j < size; j++){
    for(int i = 0; i < size; i++){
      result[0][j] += matrix[i][j];
    }
    result[0][j] /= size;
    for (int i=1;i<size;i++){
      result[i][j]=result[0][j];
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2){
  int size = matrix1.size();  
  vector<vector<double> > result(size, vector<double>(size,0));
  //TODO multiply
  return result; 
}

double sumImage(vector<vector<double> > input)
{
  double sum = 0;
  for (int i = 0; i < IMG_SIZE; i++)
  {
    for (int j = 0; j < IMG_SIZE; j++)
    {
      sum += input[i][j];
    }
  }
  return sum;
}

void saveImage(vector<vector<double> > input, string filename){
  double img[IMG_SIZE][IMG_SIZE];
  for (int i = 0;i < IMG_SIZE; i++){
    for (int j = 0;j < IMG_SIZE; j++){
      img[i][j] = input[i][j]*100;
    }
  }
  cv::Mat img_mat(IMG_SIZE,IMG_SIZE,CV_64F);
  memcpy(img_mat.data, img, IMG_SIZE*IMG_SIZE*sizeof(double));
  printf("show img!\n");
  cv::imwrite(filename,img_mat);
}