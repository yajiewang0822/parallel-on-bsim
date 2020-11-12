#include "helper.h"

//TODO: with Halide
vector<vector<double> > fastConvolution(vector<vector<double> > obj, vector<vector<double> > filter)
{
  vector<vector<double> > result;
  return result;
}

//TODO: with Halide
vector<vector<double> > inverseFastConvolution(vector<vector<double> > obj, vector<vector<double> > filter)
{
  vector<vector<double> > result;
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