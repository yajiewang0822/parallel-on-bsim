#include "helper.h"
#include "cmath"

using namespace std;
//TODO: with Halide/FFTW
vector<vector<double> > fastConvolution(vector<vector<double> > obj, vector<vector<double> > filter)
{
  vector<vector<double> > result;
  return result;
}

//TODO: with Halide/FFTW
vector<vector<double> > inverseFastConvolution(vector<vector<double> > obj, vector<vector<double> > filter)
{
  vector<vector<double> > result;
  return result;
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


double sumImage(vector<vector<double> > input, int img_size)
{
  double sum = 0;
  for (int i = 0; i < img_size; i++)
  {
    for (int j = 0; j < img_size;)
    {
      sum += input[i][j];
    }
  }
  return sum;
}