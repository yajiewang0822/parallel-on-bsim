all: BSIM BSIM.o InputGenerator.o

BSIM: InputGenerator.o BSIM.o
	g++ -Wall $^ -o $@

BSIM.o: BSIM.cpp 
	g++ -Wall -std=c++11 -c $< -o $@

InputGenerator.o: InputGenerator.cpp InputGenerator.h
	g++ -Wall -std=c++11 -c $< -o $@

clean:
	rm BSIM *.o
