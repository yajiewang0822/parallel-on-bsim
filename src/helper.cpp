#include "helper.h"
#include "cmath"

using namespace std;
//TODO: with Halide/FFTW
vector<vector<double> > fconv2(vector<vector<double> > obj, vector<vector<double> > filter)
{
  // vector<vector<double> > result(IMG_SIZE, vector<double>(IMG_SIZE,0));
  // fftw_complex *obj_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  // fftw_complex *filter_fft= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  
  // fft2(obj, obj_fft);
  // fft2(filter, filter_fft);
  // cout<<"obj: "<<obj_fft[1][0]<<endl;
  // cout<<"filter: "<<filter_fft[1][0]<<endl;

  // fftw_complex *result_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  // for(int i = 0; i < IMG_SIZE * IMG_SIZE; i++){
  //   result_fft[i][0] = obj_fft[i][0] * filter_fft[i][0] - obj_fft[i][1] * filter_fft[i][1];
  //   result_fft[i][1] = obj_fft[i][0] * filter_fft[i][1] + obj_fft[i][1] * filter_fft[i][0];
  // }

  // result = ifft2(result_fft);

  // double *result_in = (double*) malloc(sizeof(double) * IMG_SIZE * IMG_SIZE);
  // double *result_out = (double*) malloc(sizeof(double) * IMG_SIZE * IMG_SIZE);
  // for (int i=0;i<IMG_SIZE;i++){
  //   for (int j=0;j<IMG_SIZE;j++){
  //     result_in[(i*IMG_SIZE)+j] = result[i][j];
  //   }
  // }
  // fftshift_double(result_out,result_in, IMG_SIZE, IMG_SIZE);
  // for (int i=0;i<IMG_SIZE;i++){
  //   for (int j=0;j<IMG_SIZE;j++){
  //     result[i][j] = result_out[(i*IMG_SIZE)+j];
  //   }
  // }
  // double sum = sumImage(filter);
  // result = matrixScalarMul(result, (1.0/sum));
  // fftw_free(obj_fft);
  // fftw_free(filter_fft);
  // fftw_free(result_fft);
  // return result;

  vector<vector<double> > result(IMG_SIZE, vector<double>(IMG_SIZE,0));
  // fftw_complex *obj_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  // fft2(obj, obj_fft);
  // result = ifft2(obj_fft);

  vector<vector<double> > in(3, vector<double>(3, 0));
  double *try_f=(double*) malloc(sizeof(double)*3*3);
  double *try_ff=(double*) malloc(sizeof(double)*3*3);
  // fftw_complex *try_f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 256*256);
  // fftw_complex *try_ff= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 256*256);
  fftw_complex *try_fft= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 3*(3/2+1));
  for (int i=0;i<3;i++){
    for (int j=0;j<3; j++){
      try_f[i*3+j]=i*3+j;
      in[i][j]=i*3+j;
      // try_f[i*256+j][1]=0.0;
    }
  }
  fftw_plan p1=fftw_plan_dft_r2c_1d(9, try_f, try_fft, FFTW_MEASURE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  for (int i=0;i<6;i++){
    cout<<try_fft[i][0]<<"+"<<try_fft[i][1]<<"i ";
  }
  cout<<endl;
  fftw_plan p2=fftw_plan_dft_c2r_1d(9, try_fft, try_ff, FFTW_MEASURE);
  fftw_execute(p2);
  fftw_destroy_plan(p2);
  for (int i=0;i<3;i++){
    for (int j=0;j<3; j++){
      cout<<try_ff[i*3+j]/(3*3)<<" ";
    }
  }
  // double *result_in = (double*) malloc(sizeof(double) * IMG_SIZE * IMG_SIZE);
  // double *result_out = (double*) malloc(sizeof(double) * IMG_SIZE * IMG_SIZE);
  // for (int i=0;i<IMG_SIZE;i++){
  //   for (int j=0;j<IMG_SIZE;j++){
  //     result_in[(i*IMG_SIZE)+j] = result[i][j];
  //   }
  // }
  // fftshift_double(result_out,result_in, IMG_SIZE, IMG_SIZE);
  // for (int i=0;i<IMG_SIZE;i++){
  //   for (int j=0;j<IMG_SIZE;j++){
  //     result[i][j] = result_out[(i*IMG_SIZE)+j];
  //   }
  // }
  return result;
}

void fft2(vector<vector<double> > input, fftw_complex *output){
  fftw_complex *input_temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  for (int i=0;i<IMG_SIZE;i++){
    for (int j=0;j<IMG_SIZE;j++){
      input_temp[(i*IMG_SIZE)+j][0]=input[i][j];
      input_temp[(i*IMG_SIZE)+j][1]=0.0;
    }
  }
  fftw_plan p=fftw_plan_dft_2d(IMG_SIZE, IMG_SIZE, input_temp, output, FFTW_FORWARD, FFTW_MEASURE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_free(input_temp);
}

vector<vector<double > > ifft2(fftw_complex *input){
  fftw_complex *output_temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * IMG_SIZE * IMG_SIZE);
  vector<vector<double> > output(IMG_SIZE,vector<double>(IMG_SIZE,0));
  fftw_plan p=fftw_plan_dft_2d(IMG_SIZE, IMG_SIZE, input, output_temp, FFTW_BACKWARD, FFTW_MEASURE);
  fftw_execute(p);
  for (int i=0;i<IMG_SIZE;i++){
    for (int j=0;j<IMG_SIZE;j++){
      output[i][j] = output_temp[(i*IMG_SIZE)+j][0]/(IMG_SIZE*IMG_SIZE);
      if ( output_temp[(i*IMG_SIZE)+j][1]>0.0001){
        cout<<output_temp[(i*IMG_SIZE)+j][1]<<" ";
      }
    }
  }
  fftw_destroy_plan(p);
  fftw_free(output_temp);
  return output;
}

//TODO: with Halide
vector<vector<double> > matrixEleMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
  int size = matrix1.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix1[i][j] * matrix2[i][j];
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixSub(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
  int size = matrix1.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix1[i][j] - matrix2[i][j];
    }
  }
  return result;
}

//TODO: with Halide
vector<vector<double> > matrixAdd(vector<vector<double> > matrix1, vector<vector<double> > matrix2)
{
  int size = matrix1.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix1[i][j] + matrix2[i][j];
    }
  }
  return result;
}

//TODO: with Halide
vector<vector<double> > matrixScalarMul(vector<vector<double> > matrix, double multiplier)
{
  int size = matrix.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = matrix[i][j] * multiplier;
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixAbs(vector<vector<double> > matrix)
{
  int size = matrix.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      result[i][j] = abs(matrix[i][j]);
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixColMean(vector<vector<double> > matrix){
  int size = matrix.size();
  vector<vector<double> > result(size, vector<double>(size,0));
  for(int j = 0; j < size; j++){
    for(int i = 0; i < size; i++){
      result[0][j] += matrix[i][j];
    }
    result[0][j] /= size;
    for (int i=1;i<size;i++){
      result[i][j]=result[0][j];
    }
  }
  return result; 
}

//TODO: with Halide
vector<vector<double> > matrixMul(vector<vector<double> > matrix1, vector<vector<double> > matrix2){
  int size = matrix1.size();  
  vector<vector<double> > result(size, vector<double>(size,0));
  
  Eigen::MatrixXd mat1(IMG_SIZE,IMG_SIZE);
  Eigen::MatrixXd mat2(IMG_SIZE,IMG_SIZE);
  for (int i=0;i<IMG_SIZE;i++){
    mat1.row(i)=Eigen::VectorXd::Map(&matrix1[i][0],IMG_SIZE);
    mat2.row(i)=Eigen::VectorXd::Map(&matrix2[i][0],IMG_SIZE);
  }
  Eigen::MatrixXd mat3=mat1*mat2;
  for (int i=0;i<IMG_SIZE;i++){
    double* begin = &mat3.row(i).data()[0];
    result.push_back(vector<double>(begin, begin + mat3.cols()));
  }
  return result; 
}

double sumImage(vector<vector<double> > input)
{
  double sum = 0;
  for (int i = 0; i < IMG_SIZE; i++)
  {
    for (int j = 0; j < IMG_SIZE; j++)
    {
      sum += input[i][j];
    }
  }
  return sum;
}

void saveImage(vector<vector<double> > input, string filename){
  double img[IMG_SIZE][IMG_SIZE];
  for (int i = 0;i < IMG_SIZE; i++){
    for (int j = 0;j < IMG_SIZE; j++){
      img[i][j] = input[i][j];
    }
  }
  cv::Mat img_mat(IMG_SIZE,IMG_SIZE,CV_64F);
  cv::Mat img_mat_tmp(IMG_SIZE,IMG_SIZE,CV_64F);
  memcpy(img_mat_tmp.data, img, IMG_SIZE*IMG_SIZE*sizeof(double));
  cv::normalize(img_mat_tmp,img_mat, 255.0, 0.0, cv::NORM_MINMAX,-1, cv::noArray());
  cv::imwrite(filename,img_mat);
}

void saveData(vector<vector<double> > input, string filename){
  ofstream myfile;
  myfile.open (filename);
  for (int i = 0;i < IMG_SIZE; i++){
    for (int j = 0;j < IMG_SIZE; j++){
      myfile << input[i][j] << " ";
    }
  }
  myfile.close();
}

vector<vector<double> > readData(string filename){
  vector<vector<double> > output(IMG_SIZE, vector<double>(IMG_SIZE,0));
  ifstream myfile;
  myfile.open (filename);
  for (int i = 0;i < IMG_SIZE; i++){
    for (int j = 0;j < IMG_SIZE; j++){
      myfile >> output[i][j];
    }
  }
  myfile.close();
  return output;
}

void circshift(fftw_complex *out, fftw_complex *in, int xdim, int ydim, int xshift, int yshift){
  for (int i = 0; i < xdim; i++) {
    int ii = (i + xshift) % xdim;
    for (int j = 0; j < ydim; j++) {
      int jj = (j + yshift) % ydim;
      out[ii * ydim + jj][0] = in[i * ydim + j][0];
      out[ii * ydim + jj][1] = in[i * ydim + j][1];
    }
  }
}

void circshift_double(double *out, double *in, int xdim, int ydim, int xshift, int yshift){
  for (int i = 0; i < xdim; i++) {
    int ii = (i + xshift) % xdim;
    for (int j = 0; j < ydim; j++) {
      int jj = (j + yshift) % ydim;
      out[ii * ydim + jj] = in[i * ydim + j];
    }
  }
}