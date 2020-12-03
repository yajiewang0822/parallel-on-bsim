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
#include <mpi.h>
#include <stdio.h>
#include <ctime>


BSIM::BSIM(int pat_num{
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
vector<vector<double> > BSIM::Reconstruction(){

  if (procID == ROOT){
    
    // receive cost_value

    // TODO: save image in the end in ROOT
    //saveImage(result, "imgs/outputs/result.jpg");
  } else {
    vector<vector<vector<double> > > inputs(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
    vector<vector<double> > psf;
    vector<vector<double> > psfn;
    int startIndex = procID * pat_num;
    for (int i = startIndex; i < startIndex + pat_num; i++){
      inputs[i - startIndex]=readData("data/inputs/input_"+to_string(i)+".txt");
    }
    psf=readData("data/psf.txt");
    psfn=readData("data/psfn.txt");

    // patterns estimation
  vector<vector<vector<double> > > patterns = patternEstimation(inputs, psfn, psf);
  }

  

  //covariance of inputs&patterns
  printf("begin covariance\n"); 
  vector<vector<double> > covar(IMG_SIZE, vector<double>(IMG_SIZE, 0));
  vector<double> input(pat_num, 0);
  vector<double> pattern(pat_num, 0);
  for (int i = 0; i < IMG_SIZE; i++){
    for (int j = 0; j < IMG_SIZE; j++){
      double mean_input = 0.0;
      double mean_pat = 0.0;
      // MPI TODO
      // 1. Each processor added all local patterns based on pixels and return the sum array(256*256)
      // 2. Master added the returned sum array from each processor based on pixels
      for (int k = 0; k < this->pat_num; k++){
        input[k] = inputs[k][i][j];
        pattern[k] = patterns[k][i][j];
        mean_input += input[k];
        mean_pat += pattern[k];
      }
      // Same with sequential: calculate the covariance of input and pattern based on pixels
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

vector<vector<vector<double> > > BSIM::patternEstimation(vector<vector<vector<double> > > inputs, vector<vector<double> > psfn, vector<vector<double> > psf){
  //initial guess on objective
  int pat_num = this->pat_num;
  double coefficient=0.0;
  vector<vector<double> > obj(IMG_SIZE, vector<double>(IMG_SIZE,0));
  for(int i = 0; i < IMG_SIZE; i++){
    for(int j = 0; j < IMG_SIZE; j++){
      double temp = 0;
      for(int k = 0; k < pat_num; k++){
        temp += inputs[k][i][j];
      }
      obj[i][j] = temp/pat_num;
      coefficient = max(coefficient, obj[i][j]);
    }
  }
  obj=matrixScalarMul(obj, (1.0/coefficient));

//initial guess on the illumination patterns
  double t_prev = 1.0;
  double step = 1.0;
  vector<vector<vector<double> > > patterns(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, coefficient*0.1)));
  vector<vector<vector<double> > > patterns_next(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
  
//fast proximal gradient descent
// Move Part A out of while loop, make it a variable cost_prev
//
// Concerns: parallel work may prevent gradient descending
//
// MPI process:
// 1. master does work/ master does not do work
// 2. Divide num_patterns based on num_processes, assign patterns to each, locally generate 
// 3. Each processor runs Part B, C, D, Send the cost value in D if cost value is small enough.

  int i=0;
  do{

    //update coefficient
    double t_curr=0.5*(1+sqrt(1+4*(t_prev*t_prev)));
    double alpha = (t_prev-1)/t_curr;
    t_prev = t_curr;


    // Part A claculate residual and the cost
    vector<vector<vector<double> > > residual(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
    double cost_value = 0.0, cost_value_next = 0.0;
    for (int j=0; j < pat_num; j++){
      // res = input - fconv(pat*obj, psf);
      residual[j] = matrixSub(inputs[j],fconv2(matrixEleMul(patterns[j],obj),psf));      
      // f = sum(abs(res),3);
      cost_value += sumImage(matrixAbs(residual[j]));
    }
    

    // Part B calculate gradient
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


    // Part C recalculate the residual and cost
    for (int j=0; j < pat_num; j++){
      // res = input - fconv(pat * obj, psf);
      residual[j] = matrixSub(inputs[j], fconv2(matrixEleMul(patterns_next[j],obj),psf));
      // f = sum(abs(res),3);
      cost_value_next += sumImage(matrixAbs(residual[j]));
    }


    // Part D if the cost value increases, discard the update and decrease the step
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
  return patterns;
}

int main(int argc, char **argv)
{
  // TODO: change to cycleTime.h
  time_t now = time(0);
  char* dt = ctime(&now);
  printf("The local date and time is: ");
  printf("%s", dt);

  // MPI init
  int num_processors = 0;
  int procID;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
  int pat_num = PATTERN_NUM / num_processors;

  // TODO
  // master: generate all the inputs, broadcast to all
  // other processor: read data after receive broadcast
  if(procID == ROOT){
    InputGenerator *inputGenerator = new InputGenerator(NA_SPEC, PATTERN_NUM, PIXEL_SIZE);
    inputGenerator->GenerateInputs();
    BSIM *bsim = new BSIM(PATTERN_NUM);
  } else {
    BSIM *bsim = new BSIM(pat_num);
  }
  bsim->Reconstruction();

  //MPI Finish
  MPI_Finalize();

  // TODO: change to cycleTime.h
  time_t fin = time(0);
  char* dt_fin = ctime(&fin);
  printf("The local date and time is: ");
  printf("%s", dt_fin);

  
  printf("Success!\n");
  return 0;
}