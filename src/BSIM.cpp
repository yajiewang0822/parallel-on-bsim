#include "InputGenerator.h"
#include "helper.h"
#include <stdio.h>

#define NA_SPEC 0.9
#define PATTERN_NUM 384
#define IMG_SIZE 512
#define PIXEL_SIZE 30

BSIM::BSIM(int pat_num){
  this->pat_num=pat_num;
}

int main()
{
  InputGenerator *inputGenerator = new InputGenerator(NA_SPEC, PATTERN_NUM, IMG_SIZE, PIXEL_SIZE);
  inputGenerator->GenerateInputs();
  return 0;
}