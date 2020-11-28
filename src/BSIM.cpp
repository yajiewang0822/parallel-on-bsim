/**
 * This file contains the main function to run the simulation. 
 * It includes the reconstruction step involveing proximal gradient descent
 * 
 * @author: Peicheng Tang 
 * @author: Yajie Wang
 */
#include "InputGenerator.h"

#include "helper.h"
#include "BSIM.h"
#include <stdio.h>
#include <ctime>

BSIM::BSIM(int pat_num){
  this->pat_num = pat_num;
}

/**
 * Reconstruct the high-resolution result based on low-resolution inputs
 * Using fast proximal gradient descending algorithm to get the final result
 * 
 * @param inputs, low-resolution inputs
 * @param psfn, point spread function of illumination
 * @param psf, point spread function of the microscope
 * 
 * @return high-resolution result
 */
vector<vector<double> > BSIM::Reconstruction(vector<vector<vector<double> > > inputs, vector<vector<double> > psfn, vector<vector<double> > psf){
  
  //initial guess on objective
  int pat_num = this->pat_num;
  vector<vector<double> > obj(IMG_SIZE, vector<double>(IMG_SIZE,0));
  double sum = 0;
  for(int i = 0; i < IMG_SIZE; i++){
    for(int j = 0; j < IMG_SIZE; j++){
      double temp = 0;
      for(int k = 0; k < pat_num; k++){
        temp += inputs[k][i][j];
        sum += inputs[k][i][j];
      }
      obj[i][j] = sqrt(temp/pat_num);
    }
  }
  saveImage(inputs[0], "data/outputs/input.jpg");

  //initial guess on the illumination patterns
  double avg = sum/IMG_SIZE/IMG_SIZE;
  double t_prev = 1.0;
  double step = 1.0;
  vector<vector<double> > I_total(IMG_SIZE, vector<double>(IMG_SIZE, sum));
  vector<vector<vector<double> > > patterns(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, avg/pat_num)));
  vector<vector<vector<double> > > patterns_next(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
  
  
  //fast proximal gradient descent
  int i=0;
  do{

    //update co-effective
    double t_curr=0.5*(1+sqrt(1+4*(t_prev*t_prev)));
    double alpha = t_prev/t_curr;
    t_prev = t_curr;


    // claculate residual and the cost
    vector<vector<vector<double> > > residual(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
    double cost_value = 0.0, cost_value_next = 0.0;
    for (int j=0; j < pat_num; j++){
      // res = input - fconv(pat*obj, psf);
      residual[j] = matrixSub(inputs[j],fconv2(matrixEleMul(patterns[i],obj),psf));      
      // f = sum(abs(res),3);
      cost_value += sumImage(matrixAbs(residual[j]));
    }


    // calculate gradient
    vector<vector<double> > gradient(IMG_SIZE,vector<double>(IMG_SIZE, 0));
    vector<vector<double> > sum_pat(IMG_SIZE, vector<double>(IMG_SIZE, 0));
    for (int j=0; j < pat_num - 1; j++){
      // g = -2 * obj * fconv2(res * psf);
      gradient = matrixScalarMul(matrixEleMul(obj, fconv2(residual[j],psf)),-2.0);
      // pat_next = pat - g * step;
      patterns_next[j] = matrixSub(patterns[j], matrixScalarMul(gradient, step));
      // pat_next = pat + alpha * (pat_next - pat);
      patterns_next[j] = matrixAdd(patterns[j], matrixScalarMul(matrixSub(patterns_next[j], patterns[j]), alpha));
      // sum_pat = sum(pat_next, 3);
      sum_pat = matrixAdd(sum_pat,patterns_next[j]);
    }
    // pat_next[last]= I_total - sum_pat;
    patterns_next[pat_num - 1] = matrixSub(I_total,sum_pat);


    // recalculate the residual and cost
    for (int j=0; j < pat_num; j++){
      // res = input - fconv(pat * obj, psf);
      residual[j] = matrixSub(inputs[j], fconv2(matrixEleMul(patterns_next[i],obj),psf));
      // f = sum(abs(res),3);
      cost_value_next += sumImage(matrixAbs(residual[j]));
    }
    saveImage(patterns_next[0], "data/outputs/pat_next.jpg");


    // if the cost value increases, descard the update and decrease the step
    if (cost_value_next>cost_value){
      printf("decard update\n");
      patterns_next = patterns;
      step = step / 2;
    } else {
      patterns = patterns_next;
      i++;
    }
  } while (i < ITER_NUM);
  saveImage(patterns_next[0], "data/outputs/pat.jpg");


  //covariance of inputs&patterns
  printf("begin covariance\n");
  vector<vector<double> > covar(IMG_SIZE, vector<double>(IMG_SIZE, 0));
  vector<double> input(pat_num, 0);
  vector<double> pattern(pat_num, 0);
  for (int i = 0; i < IMG_SIZE; i++){
    for (int j = 0; j < IMG_SIZE; j++){
      double mean_input = 0.0;
      double mean_pat = 0.0;
      for (int k = 0; k < this->pat_num; k++){
        input[k] = inputs[k][i][j];
        pattern[k] = patterns[k][i][j];
        mean_input += input[k];
        mean_pat += pattern[k];
      }
      mean_input /= this->pat_num;
      mean_pat /= this->pat_num;
      for (int k = 0; k < this->pat_num; k++){
        input[k] -= mean_input;
        pattern[k] -= mean_pat;
        covar[i][j] += input[k] * pattern[k];
      }
      covar[i][j] /= (pat_num - 1);
    }
  }

  return covar;
}

int main()
{
  time_t now = time(0);
   
  // convert now to string form
  char* dt = ctime(&now);

  printf("The local date and time is: ");
  printf("%s", dt);
  
  // InputGenerator *inputGenerator = new InputGenerator(NA_SPEC, PATTERN_NUM, PIXEL_SIZE);
  // vector<vector<vector<double> > > inputs = inputGenerator->GenerateInputs();
  
  BSIM *bsim = new BSIM(PATTERN_NUM);
  vector<vector<vector<double> > > inputs(PATTERN_NUM, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
  vector<vector<double> > psf;
  vector<vector<double> > psfn;
  for (int i=0;i<PATTERN_NUM;i++){
    inputs[i]=readData("data/inputs/input_"+to_string(i)+".txt");
  }
  psf=readData("data/psf.txt");
  psfn=readData("data/psfn.txt");

  vector<vector<double> > result = bsim->Reconstruction(inputs, psfn, psf);
  saveImage(result, "data/outputs/result.jpg");

  time_t fin = time(0);
   
  // convert now to string form
  char* dt_fin = ctime(&fin);
  printf("The local date and time is: ");
  printf("%s", dt_fin);
  printf("Success!\n");
  return 0;
}