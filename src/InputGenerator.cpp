#include "InputGenerator.h"
#include <cmath>
#include <stdint.h>
#include <random>
#include <algorithm>
#include "helper.h"
using namespace std;

InputGenerator::InputGenerator(double NA_spec, int pattern_num, int img_size, int p_size)
{
  this->NA_spec = NA_spec;
  this->pattern_num = pattern_num;
  this->img_size = img_size;
  this->p_size = p_size;
}


vector<vector<vector<double> > > InputGenerator::GenerateInputs()
{
  int pat_num = this->pattern_num;
  int size = this->img_size;
  vector<vector<double> > objective = this->GenerateObjective();
  // saveImage(objective, size, "objective.jpg");
  double obj_sum = sumImage(objective,size);
  GeneratePSFandOTF(NA, PSF);
  // GeneratePSFandOTF(NA_spec);
  // vector<vector<double> > widefield = fastConvolution(objective, this->psf);
  // vector<vector<vector<double> > > patterns = this->GeneratePatterns();
  vector<vector<double> > pat_mean(pat_num, vector<double>(pat_num, 0));
  vector<vector<vector<double> > > inputs(pat_num, vector<vector<double> >(size, vector<double>(size, 0)));
  // for (int i = 0; i < pat_num; i++)
  // {
  //   for (int j = 0; j < size; j++)
  //   {
  //     for (int k = 0; k < size; k++)
  //     {
  //       inputs[i][j][k] = objective[j][k] * patterns[i][j][k];
  //     }
  //   }
  //   inputs[i] = fastConvolution(inputs[i], this->psf);
  // }
  return inputs;
}

void InputGenerator::GeneratePSFandOTF(double effect_NA, PSF_TYPE type)
{
  int size = this->img_size;
  int xc = round(size / 2);
  int yc = round(size / 2);
  vector<vector<int> > X;
  vector<int> X_row;
  for (int i = 0; i < size; i++)
  {
    X_row.push_back(i + 1 - xc);
  }
  for (int i = 0; i < size; i++)
  {
    X.push_back(X_row);
  }
  vector<vector<int> > Y;
  for (int i = 0; i < size; i++)
  {

    vector<int> Y_row;
    for (int j = 0; j < size; j++)
    {
      Y_row.push_back(j + 1 - yc);
    }
    Y.push_back(Y_row);
  }
  double scale=2*PI/LAMBDA*NA*p_size;
  
  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      double temp = sqrt(X[i][j] * X[i][j] + Y[i][j] * Y[i][j]);
      switch (type)
      {
      case PSF:
        this->psf[i][j]=cyl_bessel_i(1.0, scale*temp);
        break;
      case PSFN:
        break;
      default:
        break;
      }
    }
  }
//TODO:  Use bessel function this->psf and FFT
// psf=abs(2*besselj(1,2*pi./lambda*NA*R*psize+eps,1)...
//     ./(2*pi./lambda*NA*R*psize+eps)).^2;
// psf0=psf/max(max(psf));

// %Generate OTF
// OTF2d=fftshift(fft2(psf));
// OTF2dmax = max(max(abs(OTF2d)));
// OTF2d = OTF2d./OTF2dmax;
// OTF2dc = abs(OTF2d);

}

vector<vector<double> > InputGenerator::GenerateObjective()
{
  double pixel_resolution = 0.5 * LAMBDA / this->NA_spec / this->p_size;
  int size = this->img_size;
  vector<vector<double> > result(size, vector<double>(size, 0));
  int width = round(pixel_resolution + 7);
  
  int gap = 5;
  int y0 = 1;
  int width_y = floor((size - gap * 3) / 4);
  int num_bar = 3;
  int x0;
  while (y0 + width_y <= size)
  {
    x0 = gap;
    num_bar++;
    while ((x0 + num_bar * width) < size)
    {
      for (int i = 0; i < num_bar * width; i++)
      {
        for (int j = 0; j < width_y; j++)
        {
          result[j + y0][i + x0] = sin(2 * PI / width * i) + 1;
        }
      }
      
      x0 += num_bar * width + gap;
      width--;
      if (width <= 0)
      {
        break;
      }
    }
    y0 += width_y + gap;
    if (width <= 0)
    {
      break;
    }
  }
  // for (int i=0;i<size;i++){
  //   printf("%f ",result[10][i]);
  // }
  return result;
}

vector<vector<double> > InputGenerator::getPSF()
{
  return this->psf;
}

vector<vector<double> > InputGenerator::getPSFn()
{
  return this->psfn;
}

vector<vector<vector<double> > > InputGenerator::GeneratePatterns()
{
  int pat_num = this->pattern_num;
  int size = this->img_size;
  vector<vector<vector<double> > > pattern(pat_num, vector<vector<double> >(size, vector<double>(size, 0)));
  vector<int> idx(size * size);
  for (int k = 0; k < 3; k++)
  {
    for (int i = 0; i < size * size; i++)
    {
      idx[i] = i;
    }
    random_shuffle(idx.begin(), idx.end());
    int offset = (k - 1) * pat_num / 3;
    for (int i = 0; i < pat_num / 3; i++)
    {
      int count = 0;
      while (count < NUM_SPECKLE)
      {
        int index = idx[(i - 1) * NUM_SPECKLE + count];
        int pos_y = index / size;
        int pos_x = index - pos_y * size;
        pattern[offset + i][pos_x][pos_y] = 1;
        count += 1;
      }
      pattern[offset + i]=fastConvolution(pattern[offset + i],this->psfn);
    }
  }

  return pattern;
}