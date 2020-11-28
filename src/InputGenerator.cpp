/**
 * This file contains functions to generate input images we need for the simulaiton
 * including input image, pattern, and point spread function and etc.
 * 
 * @author: Peicheng Tang 
 * @author: Yajie Wang
 */
#include "InputGenerator.h"
#include <boost/math/special_functions/bessel.hpp>
#include "helper.h"

using namespace std;

InputGenerator::InputGenerator(double NA_spec, int pattern_num, int p_size)
{
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

/**
 * This function generates inputs needed for the simulation
 * 
 * @return a 3-D matrix containing all the inputs. Each input is the result of convolution 
 * between images under specific illumination pattern and the microscopy system.
 */
vector<vector<vector<double> > > InputGenerator::GenerateInputs()
{
  
  // generate objective
  vector<vector<double> > objective = this->GenerateObjective();

  // generate psf and psfn
  GeneratePSF(NA, PSF);
  this->psf[IMG_SIZE/2-1][IMG_SIZE/2-1]=1;
  GeneratePSF(NA_spec, PSFN);
  this->psfn[IMG_SIZE/2-1][IMG_SIZE/2-1]=1;

  // generate widefield
  vector<vector<double> > widefield = fconv2(objective, this->psf);

  // NOTE: uncomment the following if you want to view the image or view the raw data 
  // saveImage(this->psf, "imgs/psf.jpg");
  // saveData(this->psf, "data/psf.txt");
  // saveImage(this->psfn, "imgs/psfn.jpg");
  // saveData(this->psfn, "data/psfn.txt");
  // saveImage(widefield, "imgs/widefield.jpg");
  // saveData(widefield, "data/widefield.txt");

  //genearte patterns
  int pat_num = this->pattern_num;
  vector<vector<vector<double> > > patterns = this->GeneratePatterns();

  vector<vector<double> > pat_mean(pat_num, vector<double>(pat_num, 0));
  vector<vector<vector<double> > > inputs(pat_num, vector<vector<double> >(IMG_SIZE, vector<double>(IMG_SIZE, 0)));
  for (int i = 0; i < pat_num; i++)
  {
    for (int j = 0; j < IMG_SIZE; j++)
    {
      for (int k = 0; k < IMG_SIZE; k++)
      {
        inputs[i][j][k] = objective[j][k] * patterns[i][j][k];
      }
    }
    inputs[i] = fconv2(inputs[i], this->psf);

    // NOTE: uncomment the following if you want to view the image or view the raw data 
    // saveImage(inputs[i], "imgs/inputs/input_" + to_string(i)+".jpg");
    // saveData(inputs[i], "data/inputs/input_" + to_string(i)+".txt");
  }
  return inputs;
}


/**
 * This function generates point spread function of specific numerical apecture based on bessel function
 * 
 * @param effect_NA
 * @param type
 */
void InputGenerator::GeneratePSF(double effect_NA, PSF_TYPE type)
{
  //Use bessel function to calculate psf
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
}
/**
 * Generate the objective used for simulation.
 * The objective is designed to quickly get the resolution result, which
 * contains serveral-level resolution blocks.
 */
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
  // NOTE: uncomment the following if you want to view the image or view the raw data 
  // saveImage(result, "imgs/objective.jpg");
  // saveData(result, "data/obj.txt");
  return result;
}

/**
 * This function obtains point spread function
 * 
 * @return the psf of the simulation
 */
vector<vector<double> > InputGenerator::getPSF()
{
  return this->psf;
}

/**
 * This function obtains point spread function of illumination
 * 
 * @return the psfn of the simulation
 */
vector<vector<double> > InputGenerator::getPSFn()
{
  return this->psfn;
}

/**
 * This function generates the illumination patterns required for the simulation
 * 
 * @return a 3-D matrix containing all the patterns
 */
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
    int offset = k * pat_num / 3;
    random_shuffle(idx.begin(), idx.end());
    
    for (int i = 0; i < pat_num / 3; i++)
    {
      int count = 0;
      while (count < NUM_SPECKLE)
      {

        int index = idx[i * NUM_SPECKLE + count];
        
        int pos_y = index / IMG_SIZE;
        int pos_x = index - pos_y * IMG_SIZE;
        pattern[offset + i][pos_x][pos_y] = 1;
        count += 1;
      }
      pattern[offset + i]=fconv2(pattern[offset + i],this->psfn);
      // NOTE: uncomment the following if you want to view the image
      // saveImage(pattern[offset + i], "imgs/pats/pat_" + to_string(offset + i)+".jpg");
    }
  }

  return pattern;
}