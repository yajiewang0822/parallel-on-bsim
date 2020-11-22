all: BSIM 

BSIM: InputGenerator.o helper.o BSIM.o
	g++ -Wall $^ -o $@

BSIM.o: BSIM.cpp 
	g++ -Wall -std=c++11 -c $< -o $@ -I ../fft2-3.3.8/api/

helper.o: helper.cpp helper.h
	g++ -Wall -std=c++11 -c helper.cpp -o $@

InputGenerator.o: InputGenerator.cpp InputGenerator.h helper.cpp helper.h 
	g++ -Wall -std=c++11 -c InputGenerator.cpp -o $@

clean:
	rm BSIM *.o
