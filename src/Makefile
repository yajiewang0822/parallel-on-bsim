all: BSIM 

BSIM: InputGenerator.o helper.o BSIM.o
	mpicxx -Wall $^ -o $@ -lopencv_core -lopencv_highgui -lfftw3 -MMD -MP

BSIM.o: BSIM_sequential.cpp 
	mpicxx -Wall -std=c++11 -c $< -o $@ -I /usr/include/ -MMD -MP

InputGenerator.o: InputGenerator.cpp InputGenerator.h helper.cpp helper.h
	mpicxx -Wall -std=c++11 -c InputGenerator.cpp -o $@ -I /usr/include/ -MMD -MP

helper.o: helper.cpp helper.h
	mpicxx -Wall -std=c++11 -c helper.cpp -o $@ -I /usr/include/ -MMD -MP

clean:
	rm BSIM *.o
