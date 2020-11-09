#include "InputGenerator.h"
#include <cmath>
#include <cstdint>
using namespace std;

InputGenerator::InputGenerator(double NA_spec, int pattern_num, int img_size,int p_size){
  this->NA_spec = NA_spec;
  this->lambda = lambda;
  this->pattern_num = pattern_num;
  this->img_size = img_size;
  this->p_size = p_size;
}

double InputGenerator::sumImage(vector<vector<double>> input){
  double sum=0;
  for (int i=0;i<img_size;i++){
    for (int j=0;j<img_size;){
      sum+=input[i][j];
    }
  }
  return sum;
}

vector<vector<vector<double>>> InputGenerator::GenerateInputs(){
  vector<vector<double>> objective=this.GenerateObjective();
  obj_sum=sumImage(objective);
  GeneratePSFandOTF(NA);
  GeneratePSFandOTF(NA_spec);
  vector<vector<double>> widefield=fastConvolution(objective,this->psf);
  vector<vector<vector<double>>> patterns=this.GeneratePatterns();
  for (int i=0;i<this->pattern_num;i++){
    
  }
  this->pattern_sum;
}

void GeneratePSFandOTF(double effect_NA){
  int size = this->img_size;
  int xc = round(size/2);
  int yc = round(size/2);
  vector<vector<int>> X;
  vector<int> X_row; 
  for(int i = 0; i < size; i++){
    X_row.push_back(i+1-xc);
  }
  for(int i = 0; i < size; i++){
    X.push_back(X_row);
  }
  vector<vector<int>> Y;
  for(int i = 0; i < size; i++){

    vector<int> Y_row; 
    for(int j = 0; j < size; j++){
      Y_row.push_back(j+1-yc);
    }
    Y.push_back(Y_row);
  }
  vector<vector<double>> R;
  for(int i = 0; i < size; i++){
    vector<int> R_row; 
    for(int j = 0; j < size; j++){
      double temp = sqrt( X[i][j] * X[i][j] + Y[i][j] * Y[i][j]);
      R_row.push_back(temp);
    }
    R.push_back(R_row);
  }
  // Use bessel function this->psf

  // Use FFT 
}

vector<vector<double>> InputGenerator::GenerateObjective(){
  double pixel_resolution = 0.5 * this->lambda / this->NA_spec / this->p_size;
  int size = this->img_size;
  vector<vector<double>> result(size, vector<double>(size, 0));
  uint32_t width = round(pixel_resolution + 7);
  uint32_t gap = 5;
  uint32_t y0 = 1;
  uint32_t width_y = floor((size - gap * 3) / 4);
  uint32_t num_bar = 3;
  uint32_t x0;
  while(y0 + width_y <= size){
    x0 = gap;
    num_bar++;
    while((x0 + num_bar*width) < size){
      for(int i = 0; i < num_bar*width; i++){
        for(int j = 0; j < width_y; j++){
          result[j+y0][i+x0] = sin(2*PI/width*i)+1;
        }
      }
      x0 += num_bar*width+gap;
      width--;
      if(width<=0){
        break;
      }
    }
    y0 += width_y+gap;
    if(width<=0){
      break;
    }
  }
}