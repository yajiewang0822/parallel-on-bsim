#include "InputGenerator.h"

int main(){
  InputGenerator* inputGenerator=new InputGenerator(0.9,384,512,30);
  inputGenerator->GenerateInputs();
  return 0;
}