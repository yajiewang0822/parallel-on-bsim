#include "InputGenerator.h"

#include "helper.h"
#include "BSIM.h"
#include <stdio.h>
#include <cmath>




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
  obj = deconvolution(obj,psf);
  double avg = sum/size/size;
  vector<vector<double> > I0(size, vector<double>(size, avg));
  vector<vector<vector<double> > > patterns(pat_num, vector<vector<double> >(size, vector<double>(size, avg/pat_num)));

  double t_prev = 1.0;
  double step = 1.0;
  
  vector<vector<vector<double> > > patterns_next(pat_num, vector<vector<double> >(size, vector<double>(size, 0)));
  //fast proximal gradient descent
  for (int i=0; i < ITER_NUM; i++){
    //update co-effective
    double t_curr=0.5*(1+sqrt(1+4*(t_prev*t_prev)));
    double alpha = t_prev/t_curr;
    t_prev = t_curr;

    // claculate residual
    vector<vector<vector<double> > > residual(pat_num, vector<vector<double> >(size, vector<double>(size, 0)));
    double cost_value = 0.0, cost_value_next = 0.0;
    for (int j=0; j < pat_num; j++){
      residual[j] = matrixSub(inputs[j],fconv2(matrixDotMul(patterns[i],obj),psf));
      cost_value += sumImage(matrixAbs(residual[j]));
    }

    // calculate gradient
    vector<vector<double> > gradient(size,vector<double>(size, 0));
    vector<vector<double> > sum_pat(size, vector<double>(size, 0));
    for (int j=0; j < pat_num - 1; j++){
      // g = -2 * obj * fconv2(res * psf)
      gradient = matrixScalarMul(matrixDotMul(obj, fconv2(residual[j],psf)),-2.0);
      // pat_next = pat - g * step
      patterns_next[j] = matrixSub(patterns[j], matrixScalarMul(gradient, step));
      // pat_next = pat ;
      patterns_next[j] = matrixAdd(patterns[j], matrixScalarMul(matrixSub(patterns_next[j], patterns[j]), alpha));
      sum_pat = matrixAdd(sum_pat,patterns_next[j]);
    }
    patterns_next[pat_num - 1]=matrixSub(I0,sum_pat);

    // if the cost value increases, descard the update and decrease the step
    for (int j=0; j < pat_num; j++){
      residual[j] = matrixSub(inputs[j],fconv2(matrixDotMul(patterns[i],obj),psf));
      cost_value_next += sumImage(matrixAbs(residual[j]));
    }
    if (cost_value_next>cost_value){
      patterns_next = patterns;
      step = step / 2;
    }
  }

  
  vector<vector<double> > covar(size, vector<double>(size, 0));
  patterns = patterns_next;
  //TODO covariance
  //step1: covariance of inputs&patterns
  for (int i = 0;i < pat_num; i++){
    covar = matrixAdd(covar, covariance(patterns[i],inputs[i])); 
  }
  covar = matrixScalarMul(covar, 1.0/(pat_num-1));
  //step2: deconvolution
  // return deconvolution(covar,matrixDotMul(ifft2(getValueOfComplex(fft2(psfn), size)), psf));
  return deconvolution(covar, psf);
}

vector<vector<double> > BSIM::covariance(vector<vector<double> > input, vector<vector<double> > pattern){
  vector<vector<double> > mean_input = matrixColMean(input);
  vector<vector<double> > mean_pattern = matrixColMean(pattern);
  input =  matrixSub(input, mean_input);
  pattern = matrixSub(pattern, mean_input);
  return matrixMul(input, pattern);
}

int main()
{
  InputGenerator *inputGenerator = new InputGenerator(NA_SPEC, PATTERN_NUM, PIXEL_SIZE);
  inputGenerator->GenerateInputs();
 
  printf("Success!\n");
  return 0;
}