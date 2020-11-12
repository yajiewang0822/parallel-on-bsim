#include "InputGenerator.h"
#include "helper.h"
#include "BSIM.h"
#include <stdio.h>
#include <cmath>

#define NA_SPEC 0.9
#define PATTERN_NUM 384
#define IMG_SIZE 512
#define PIXEL_SIZE 30

BSIM::BSIM(int pat_num, int img_size){
  this->pat_num = pat_num;
  this->img_size = img_size;
  
}

vector<vector<double> > BSIM::Reconstruction(vector<vector<vector<double> > > inputs, vector<vector<double> > psfn, vector<vector<double> > psf){
  //initial guess
  int size = this->img_size;
  int pat_num = this->pat_num;
  vector<vector<double> > obj(size, vector<double>(size,0));
  double sum = 0;
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      double temp = 0;
      for(int k = 0; k < pat_num; k++){
        temp += inputs[k][i][j];
        sum += inputs[k][i][j];
      }
      obj[i][j] = sqrt(temp/pat_num);
    }
  }
  double avg = sum/size/size;
  vector<vector<double> > I0(size, vector<double>(size, avg));
  vector<vector<vector<double> > > patterns(pat_num, vector<vector<double> >(size, vector<double>(size, avg/pat_num)));

  double t_prev = 1.0, t_curr;
  double step = 1.0;
  vector<vector<vector<double> > > patterns_prev = patterns;
  vector<vector<vector<double> > > patterns_next(pat_num, vector<vector<double> >(size, vector<double>(size, 0)));
  //fast proximal gradient descent
  for (int i=0; i < ITER_NUM; i++){
    t_curr=0.5*(1+sqrt(1+4*(t_prev*t_prev)));
    double alpha=t_prev/t_curr;
    vector<vector<vector<double> > > residual(pat_num, vector<vector<double> >(size, vector<double>(size, 0)));
    vector<vector<double> > sum_pat(size, vector<double>(size, 0));
    double cost_value, cost_value_next;
    // claculate residual
    for (int j=0; j < pat_num; j++){
      residual[j] = matrixSub(inputs[j],fastConvolution(matrixDotMul(patterns[i],obj),psf));
      cost_value += sumImage(matrixAbs(residual[j]), size);
    }

    // calculate gradient
    vector<vector<double> > gradient(size,vector<double>(size, 0));
    for (int j=0; j < pat_num - 1; j++){
      gradient = matrixScalarMul(matrixDotMul(obj, fastConvolution(residual[j],psf)),-2.0);
      patterns[j] = matrixSub(patterns_prev[j], matrixScalarMul(gradient, step));
      patterns_next[j] = matrixAdd(patterns[j], matrixScalarMul(matrixSub(patterns[j], patterns_prev[j]), alpha));
      sum_pat = matrixAdd(sum_pat,patterns_next[j]);
    }
    patterns_next[pat_num - 1]=matrixSub(I0,sum_pat);

    // update step
    for (int j=0; j < pat_num; j++){
      residual[j] = matrixSub(inputs[j],fastConvolution(matrixDotMul(patterns[i],obj),psf));
      cost_value_next += sumImage(matrixAbs(residual[j]), size);
    }
    if (cost_value_next>cost_value){
      step = step / 2;
    }
  }
  //TODO covariance
  vector<vector<double> > covar(size, vector<double>(size, 0));
  return covar;
}

int main()
{
  InputGenerator *inputGenerator = new InputGenerator(NA_SPEC, PATTERN_NUM, IMG_SIZE, PIXEL_SIZE);
  inputGenerator->GenerateInputs();
  return 0;
}