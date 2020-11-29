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
  double coeffect=0.0;
  vector<vector<double> > obj(IMG_SIZE, vector<double>(IMG_SIZE,0));
  for(int i = 0; i < IMG_SIZE; i++){
    for(int j = 0; j < IMG_SIZE; j++){
      double temp = 0;
      for(int k = 0; k < pat_num; k++){
        temp += inputs[k][i][j];
      }
      obj[i][j] = temp/pat_num;
      coeffect = max(coeffect, obj[i][j]);
    }
  }
  obj=matrixScalarMul(obj, (1.0/coeffect));


  //initial guess on the illumination patterns
  double t_prev = 1.0;
  double step = 1.0;
  vector<vector<vector<double> > > patterns(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, coeffect*0.1)));
  vector<vector<vector<double> > > patterns_next(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
  
  
  //fast proximal gradient descent
  int i=0;
  do{

    //update co-effective
    double t_curr=0.5*(1+sqrt(1+4*(t_prev*t_prev)));
    double alpha = (t_prev-1)/t_curr;
    t_prev = t_curr;


    // claculate residual and the cost
    vector<vector<vector<double> > > residual(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
    double cost_value = 0.0, cost_value_next = 0.0;
    for (int j=0; j < pat_num; j++){
      // res = input - fconv(pat*obj, psf);
      residual[j] = matrixSub(inputs[j],fconv2(matrixEleMul(patterns[j],obj),psf));      
      // f = sum(abs(res),3);
      cost_value += sumImage(matrixAbs(residual[j]));
    }


    // calculate gradient
    vector<vector<double> > gradient(IMG_SIZE,vector<double>(IMG_SIZE, 0));
    vector<vector<double> > sum_pat(IMG_SIZE, vector<double>(IMG_SIZE, 0));
    for (int j=0; j < pat_num; j++){
      // g = -2 * obj * fconv2(res * psf);
      gradient = matrixScalarMul(matrixEleMul(obj, fconv2(residual[j],psf)),-2.0);
      // pat_next = pat - g * step;
      patterns_next[j] = matrixSub(patterns[j], matrixScalarMul(gradient, step));
      // pat_next = pat + alpha * (pat_next - pat);
      patterns_next[j] = matrixAdd(patterns[j], matrixScalarMul(matrixSub(patterns_next[j], patterns[j]), alpha));
    }


    // recalculate the residual and cost
    for (int j=0; j < pat_num; j++){
      // res = input - fconv(pat * obj, psf);
      residual[j] = matrixSub(inputs[j], fconv2(matrixEleMul(patterns_next[j],obj),psf));
      // f = sum(abs(res),3);
      cost_value_next += sumImage(matrixAbs(residual[j]));
    }


    // if the cost value increases, discard the update and decrease the step
    if (cost_value_next>cost_value){
      printf("i: %d; discard update      ", i);
      fflush(stdout);
      patterns_next = patterns;
      step = step / 2;
    } else {
      printf("i: %d;     ", i);
      fflush(stdout);
      patterns = patterns_next;
      i++;
    }
  } while (i < ITER_NUM);


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

  BSIM *bsim = new BSIM(PATTERN_NUM);
  
  // Uncomment if want to get new generated data
  InputGenerator *inputGenerator = new InputGenerator(NA_SPEC, PATTERN_NUM, PIXEL_SIZE);
  vector<vector<vector<double> > > inputs = inputGenerator->GenerateInputs();
  
  // Uncomment if want to read data from the file
  // vector<vector<vector<double> > > inputs(PATTERN_NUM, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
  // vector<vector<double> > psf;
  // vector<vector<double> > psfn;
  // for (int i=0;i<PATTERN_NUM;i++){
  //   inputs[i]=readData("data/inputs/input_"+to_string(i)+".txt");
  // }
  // psf=readData("data/psf.txt");
  // psfn=readData("data/psfn.txt");
  // vector<vector<double> > result = bsim->Reconstruction(inputs, psfn, psf);

  vector<vector<double> > result = bsim->Reconstruction(inputs, inputGenerator->getPSFn(), inputGenerator->getPSF());
  saveImage(result, "imgs/outputs/result.jpg");

  time_t fin = time(0);
   
  // convert now to string form
  char* dt_fin = ctime(&fin);
  printf("The local date and time is: ");
  printf("%s", dt_fin);
  printf("Success!\n");
  return 0;
}