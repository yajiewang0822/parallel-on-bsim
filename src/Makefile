all: BSIM 

BSIM: InputGenerator.o helper.o BSIM.o
	g++ -Wall $^ -o $@ -lopencv_core -lopencv_imgcodecs -lfftw3

BSIM.o: BSIM.cpp 
	g++ -Wall -std=c++11 -c $< -o $@ -I /usr/local/include/opencv4 -I ../fftw-3.3.8/api/ -I ../eigen-3.3.8/

InputGenerator.o: InputGenerator.cpp InputGenerator.h helper.cpp helper.h 
	g++ -Wall -std=c++11 -c InputGenerator.cpp -o $@ -I /usr/local/include/opencv4 -I ../fftw-3.3.8/api/  -I ../eigen-3.3.8/

helper.o: helper.cpp helper.h
	g++ -Wall -std=c++11 -c helper.cpp -o $@ -I ../fftw-3.3.8/api/ -I /usr/local/include/opencv4 -I ../eigen-3.3.8/

clean:
	rm BSIM *.o
