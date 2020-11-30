ll: BSIM 

BSIM: InputGenerator.o helper.o BSIM.o
	g++ -Wall $^ -o $@ -lopencv_core -lopencv_highgui -lfftw3

BSIM.o: BSIM_sequential.cpp 
	g++ -Wall -std=c++11 -c $< -o $@ -I /usr/include/

InputGenerator.o: InputGenerator.cpp InputGenerator.h helper.cpp helper.h 
	g++ -Wall -std=c++11 -c InputGenerator.cpp -o $@ -I /usr/include/

helper.o: helper.cpp helper.h
	g++ -Wall -std=c++11 -c helper.cpp -o $@ -I /usr/include/

clean:
	rm BSIM *.o
