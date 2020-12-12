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

#define ROOT 0
#define WORK_TAG 0
#define STOP_TAG 1
#define COVAR_TAG 2

int APPROVE = 1;
int DISCARD = 0;
int num_processors;

BSIM::BSIM(int pat_num){
  this->pat_num = pat_num;
}


/**
 * Reconstruct the high-resolution result based on low-resolution inputs
 * Using fast proximal gradient descending algorithm to get the final result
 * 
 * @param inputs, low-resolution inputs
 * @param psf, point spread function of the microscope
 * 
 * @return high-resolution result
 */
void BSIM::Reconstruction(vector<vector<double> > inputs, vector<double> psf, int procID){
  vector<double> obj(IMG_SIZE*IMG_SIZE, 0);
  double coefficient;
    
  if (procID == ROOT){

    vector<double> single_cost_value_array(num_processors-1,0);
    vector<int> iter(num_processors-1,0);
    int done_count = 0;
    double cost_value = 0.0;

    //initial guess on objective
    int pat_num = this->pat_num;
    coefficient=0.0;

// #pragma omp parallel for   
    for(int k = 0; k < pat_num; k++){
      obj=matrixAdd(obj, inputs[k]);
    }
    obj=matrixScalarMul(obj, 1.0/pat_num);
    
#pragma omp parallel for  reduction(max:coefficient)
    for(int i = 0; i < IMG_SIZE * IMG_SIZE; i++){
      coefficient = max(coefficient, obj[i]);
    }

    obj=matrixScalarMul(obj, (1.0/coefficient));
    coefficient*=0.1;
    //send initial obj guess and coefficient
    for (int i=1; i < num_processors; i++){
      MPI_Send(&obj[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, i, WORK_TAG, MPI_COMM_WORLD);
      MPI_Send(&coefficient, 1, MPI_DOUBLE, i, WORK_TAG, MPI_COMM_WORLD);
    }
    
    // Receive initial cost value from workers
    for (int i=1;i<num_processors;i++){
      double receive;
      MPI_Status status;
      MPI_Recv(&receive, 1, MPI_DOUBLE, MPI_ANY_SOURCE, WORK_TAG, MPI_COMM_WORLD, &status);
      single_cost_value_array[status.MPI_SOURCE-1] = receive;
      cost_value += receive;
    }
    
    vector<double> cost_value_next_array(num_processors-1,cost_value);
    
    do{
      double cost_value_next = 0.0;

      // receive iteration cost_value
      double receive;
      MPI_Status status;
      MPI_Recv(&receive, 1, MPI_DOUBLE, MPI_ANY_SOURCE, WORK_TAG, MPI_COMM_WORLD, &status);
      single_cost_value_array[status.MPI_SOURCE-1] = receive;

      for (int i=1;i<num_processors;i++){    
        cost_value_next += single_cost_value_array[i];
      }

      // Part D if the cost value increases, discard the update and decrease the step
      if (cost_value_next>cost_value_next_array[status.MPI_SOURCE-1]){

        //send discard the update
        printf("proc: %d; Iter: %d; discard update      \n", status.MPI_SOURCE, iter[status.MPI_SOURCE-1]);
        fflush(stdout);     
        MPI_Send(&DISCARD, 1, MPI_INT, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD);

      } else {

        //send approve the update and ask workers to keep working
        printf("proc: %d; Iter: %d;     cost value next = %f; \n", status.MPI_SOURCE, iter[status.MPI_SOURCE-1], cost_value_next);
        fflush(stdout);
        iter[status.MPI_SOURCE-1]++;
        cost_value_next_array[status.MPI_SOURCE-1] = cost_value_next;
        if (iter[status.MPI_SOURCE-1] == ITER_NUM){
          done_count++;
          MPI_Send(&APPROVE, 1, MPI_INT, status.MPI_SOURCE, STOP_TAG, MPI_COMM_WORLD);
        } else {
          MPI_Send(&APPROVE, 1, MPI_INT, status.MPI_SOURCE, WORK_TAG, MPI_COMM_WORLD);
        } 
        
      }
    } while (done_count < (num_processors-1));

    //covariance of inputs&patterns
    printf("begin covariance\n"); 
    {
      vector<double> covar(IMG_SIZE * IMG_SIZE, 0);
      vector<double> mean_input(IMG_SIZE*IMG_SIZE, 0.0);
      vector<double> mean_pat(IMG_SIZE*IMG_SIZE, 0.0);

      // TODO: need synchronization
// #pragma omp parallel for
      for (int k = 0; k < this->pat_num; k++){
        mean_input = matrixAdd(mean_input, inputs[k]);
      }
      
      // receive sum_pat from each processor and add them
      for (int k=1;k<num_processors;k++){
        vector<double> sum_pat(IMG_SIZE*IMG_SIZE, 0);
        MPI_Recv(&sum_pat[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, COVAR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        mean_pat = matrixAdd(mean_pat, sum_pat);
      }
      
      // calculate the mean_pat and mean_input and send to each processor
      mean_input = matrixScalarMul(mean_input, 1.0/this->pat_num);
      mean_pat = matrixScalarMul(mean_pat,1.0/this->pat_num);
      for (int i=1; i < num_processors; i++){
        MPI_Send(&mean_pat[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, i, COVAR_TAG, MPI_COMM_WORLD);
        MPI_Send(&mean_input[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, i, COVAR_TAG, MPI_COMM_WORLD);
      }
    
      // receive covar(img-size*img-size) and add them together
      for (int k=1;k<num_processors;k++){
        vector<double> covar_part(IMG_SIZE*IMG_SIZE, 0);
        MPI_Recv(&covar_part[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, COVAR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        covar = matrixAdd(covar, covar_part);
      }
      covar = matrixScalarMul(covar,1.0/(this->pat_num - 1));

      // save image in the end in ROOT
      saveImage(covar, "imgs/outputs/result.jpg");
    }
  } else {

    int shiftIndex = (PATTERN_NUM / (num_processors-1)) * (procID - 1);

    // Receive guess of objective and initial coefficient from ROOT
    MPI_Recv(&obj[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, ROOT, WORK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&coefficient, 1, MPI_DOUBLE, ROOT, WORK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   

    // patterns estimation
    vector<vector<double> > estimated_patterns = patternEstimation(inputs, psf, obj, coefficient, procID);
    
    //covariance of inputs&patterns
    printf("begin covariance\n"); 
    {
      vector<double> covar(IMG_SIZE * IMG_SIZE, 0);
      vector<double> sum_pat(IMG_SIZE*IMG_SIZE, 0.0);
      vector<double> mean_pat(IMG_SIZE*IMG_SIZE, 0.0);
      vector<double> mean_input(IMG_SIZE*IMG_SIZE, 0.0);

      // Each processor added all local patterns based on pixels and return the sum array(256*256)
      // TODO: need synchronization
// #pragma omp parallel for
      for (int k = 0; k < this->pat_num; k++){
        sum_pat = matrixAdd(sum_pat, estimated_patterns[k]);
      }
      
      MPI_Send(&sum_pat[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, ROOT, COVAR_TAG, MPI_COMM_WORLD);
      
      MPI_Recv(&mean_pat[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, ROOT, COVAR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&mean_input[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, ROOT, COVAR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // TODO: need synchronization
// #pragma omp parallel for
      for (int k = 0; k < this->pat_num; k++){
        inputs[k + shiftIndex] = matrixSub(inputs[k + shiftIndex], mean_input);
        estimated_patterns[k] = matrixSub(estimated_patterns[k], mean_pat);
        
        covar = matrixAdd(covar, matrixEleMul(inputs[k + shiftIndex],estimated_patterns[k]));
      }

      MPI_Send(&covar[0], IMG_SIZE*IMG_SIZE, MPI_DOUBLE, ROOT, COVAR_TAG, MPI_COMM_WORLD);
    }
  
  }
}

vector<vector<double> > BSIM::patternEstimation(vector<vector<double> > inputs, vector<double> psf, vector<double> obj, double coefficient, int procID){
  
  //initial guess on the illumination patterns
  vector<vector<double> > patterns(pat_num, vector<double>(IMG_SIZE*IMG_SIZE, coefficient));
  vector<vector<double> > patterns_next(pat_num, vector<double>(IMG_SIZE*IMG_SIZE, 0));
    
  // calculate initial cost_value and send
  vector<vector<double> > residual(pat_num, vector<double>(IMG_SIZE * IMG_SIZE, 0));
  double cost_value = 0.0;
  int shiftIndex = (PATTERN_NUM / (num_processors-1)) * (procID - 1);
  // TODO: need synchronization
#pragma omp parallel for reduction(+: cost_value)
  for (int j=0; j < pat_num; j++){
    // res = input - fconv(pat*obj, psf);
    residual[j] = matrixSub(inputs[j + shiftIndex],fconv2(matrixEleMul(patterns[j],obj),psf));      
    // f = sum(abs(res),3);
    double temp = sumImage(matrixAbs(residual[j]));
    cost_value += temp;
  }
  MPI_Send(&cost_value, 1, MPI_DOUBLE, ROOT, WORK_TAG, MPI_COMM_WORLD);
  
  //fast proximal gradient descent
  // Move Part A out of while loop, make it a variable cost_prev
  //
  // Concerns: parallel work may prevent gradient descending
  //
  // MPI process:
  // 1. master does work/ master does not do work
  // 2. Divide num_patterns based on num_processes, assign patterns to each, locally generate 
  // 3. Each processor runs Part B, C, D, Send the cost value in D if cost value is small enough.

  //initial parameter for iteration
  double t_prev = 1.0;
  double step = 1.0;
  MPI_Status status;
  do{

    double cost_value_next = 0.0;

    //update coefficient
    double t_curr=0.5*(1+sqrt(1+4*(t_prev*t_prev)));
    double alpha = (t_prev-1)/t_curr;
    t_prev = t_curr;

    // Part B calculate gradient
    vector<double> gradient(IMG_SIZE*IMG_SIZE, 0);
    vector<double> sum_pat(IMG_SIZE*IMG_SIZE, 0);
    
#pragma omp parallel for
    for (int j=0; j < pat_num; j++){
      // g = -2 * obj * fconv2(res * psf);
      gradient = matrixScalarMul(matrixEleMul(obj, fconv2(residual[j],psf)),-2.0);
      // pat_next = pat - g * step;
      patterns_next[j] = matrixSub(patterns[j], matrixScalarMul(gradient, step));
      // pat_next = pat + alpha * (pat_next - pat);
      patterns_next[j] = matrixAdd(patterns[j], matrixScalarMul(matrixSub(patterns_next[j], patterns[j]), alpha));
    }

#pragma omp parallel for reduction(+: cost_value_next)
    // Part C recalculate the residual and cost
    for (int j=0; j < pat_num; j++){
      // res = input - fconv(pat * obj, psf);
      residual[j] = matrixSub(inputs[j + shiftIndex], fconv2(matrixEleMul(patterns_next[j],obj),psf));
      // f = sum(abs(res),3);
      double temp = sumImage(matrixAbs(residual[j]));
      cost_value_next += temp;
    }

    MPI_Send(&cost_value_next, 1, MPI_DOUBLE, ROOT, WORK_TAG, MPI_COMM_WORLD);
    
    int update;
    
    MPI_Recv(&update, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    // if the cost value increases, discard the update and decrease the step
    if (update==DISCARD){
      patterns_next = patterns;
      step = step / 2;
    } else {
      patterns = patterns_next;
    }   
  } while (status.MPI_TAG != STOP_TAG);


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
  
  int procID;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
  int pat_num = PATTERN_NUM / (num_processors - 1);
  //Generate inputs
  InputGenerator *inputGenerator = new InputGenerator(NA_SPEC, PATTERN_NUM, PIXEL_SIZE);
  vector<vector<double> > inputs = inputGenerator->GenerateInputs();

  //start reconstruction
  BSIM *bsim;
  if(procID == ROOT){
    bsim = new BSIM(PATTERN_NUM);
  } else {
    if (procID == num_processors -1){
      pat_num += PATTERN_NUM % (num_processors - 1);
    }
    bsim = new BSIM(pat_num);
  }
  bsim->Reconstruction(inputs, inputGenerator->getPSF(),procID);

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