all: BSIM 

BSIM: InputGenerator.o helper.o BSIM.o
	g++ -Wall $^ -o $@ -lopencv_core -lopencv_imgcodecs

BSIM.o: BSIM.cpp 
	g++ -Wall -std=c++11 -c $< -o $@ -I /usr/local/include/opencv4

InputGenerator.o: InputGenerator.cpp InputGenerator.h helper.cpp helper.h 
	g++ -Wall -std=c++11 -c InputGenerator.cpp -o $@ -I /usr/local/include/opencv4

helper.o: helper.cpp helper.h
	g++ -Wall -std=c++11 -c helper.cpp -o $@ -I ../fft2-3.3.8/api/ -I /usr/local/include/opencv4

clean:
	rm BSIM *.o
