#include "InputGenerator.h"
#include <cmath>
#include <stdint.h>
#include <random>
#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include "helper.h"

using namespace std;

InputGenerator::InputGenerator(double NA_spec, int pattern_num, int p_size)
{
  this->OTF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  this->OTFn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  for (int i=0;i<IMG_SIZE;i++){
    vector<double> v1(IMG_SIZE,0);
    this->psf.push_back(v1);
    vector<double> v2(IMG_SIZE,0);
    this->psfn.push_back(v2);
  }
  this->NA_spec = NA_spec;
  this->pattern_num = pattern_num;
  this->p_size = p_size;
}


vector<vector<vector<double> > > InputGenerator::GenerateInputs()
{
  int pat_num = this->pattern_num;
  
  vector<vector<double> > objective = this->GenerateObjective();
  saveImage(objective, "objective.jpg");
  double obj_sum = sumImage(objective);
  GeneratePSFandOTF(NA, PSF);
  saveImage(this->psf, "psf.jpg");
  this->psf[IMG_SIZE/2-1][IMG_SIZE/2-1]=1;
  GeneratePSFandOTF(NA_spec, PSFN);
  this->psfn[IMG_SIZE/2-1][IMG_SIZE/2-1]=1;
  vector<vector<double> > widefield = fconv2(objective, this->psf);
  
  saveImage(widefield, "widefield.jpg");
  // vector<vector<vector<double> > > patterns = this->GeneratePatterns();
  vector<vector<double> > pat_mean(pat_num, vector<double>(pat_num, 0));
  vector<vector<vector<double> > > inputs(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
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
  //Use bessel function this->psf and FFT
  // psf=abs(2*besselj(1,2*pi./lambda*NA*R*psize+eps,1)...
  //     ./(2*pi./lambda*NA*R*psize+eps)).^2;
  // psf0=psf/max(max(psf));
  int xc = round(IMG_SIZE / 2);
  int yc = round(IMG_SIZE / 2);
  double scale=2*PI/LAMBDA*NA*p_size;
  for (int i = 0; i < IMG_SIZE; i++)
  {
    double x=i + 1 - xc;
    for (int j = 0; j < IMG_SIZE; j++)
    {
      double y=j + 1 - yc;
      double temp = sqrt(x * x + y * y);
      switch (type)
      {
      case PSF:
        this->psf[i][j]=pow((2*boost::math::cyl_bessel_j(1.0, scale*temp+eps)/((scale*temp+eps))),2);
        break;
      case PSFN:
        this->psfn[i][j]=pow((2*boost::math::cyl_bessel_j(1.0, scale*temp+eps)/((scale*temp+eps))),2);
        break;
      default:
        break;
      }
    }
  }


  // %Generate OTF
  // OTF2d=fftshift(fft2(psf));
  // OTF2dmax = max(max(abs(OTF2d)));
  // OTF2d = OTF2d./OTF2dmax;
  // OTF2dc = abs(OTF2d);
  switch (type)
      {
      case PSF:
        fft2(this->psf, this->OTF);
        break;
      case PSFN:
        fft2(this->psfn, this->OTFn);
        break;
      default:
        break;
      }
}

vector<vector<double> > InputGenerator::GenerateObjective()
{
  double pixel_resolution = 0.5 * LAMBDA / this->NA_spec / this->p_size;
  vector<vector<double> > result(IMG_SIZE, vector<double>(IMG_SIZE, 0));
  int width = round(pixel_resolution + 7);
  
  int gap = 5;
  int y0 = 1;
  int width_y = floor((IMG_SIZE - gap * 3) / 4);
  int num_bar = 3;
  int x0;
  while (y0 + width_y <= IMG_SIZE)
  {
    x0 = gap;
    num_bar++;
    while ((x0 + num_bar * width) < IMG_SIZE)
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
  vector<vector<vector<double> > > pattern(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
  vector<int> idx(IMG_SIZE * IMG_SIZE);
  for (int k = 0; k < 3; k++)
  {
    for (int i = 0; i < IMG_SIZE * IMG_SIZE; i++)
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
        int pos_y = index / IMG_SIZE;
        int pos_x = index - pos_y * IMG_SIZE;
        pattern[offset + i][pos_x][pos_y] = 1;
        count += 1;
      }
      pattern[offset + i]=fconv2(pattern[offset + i],this->psfn);
    }
  }

  return pattern;
}