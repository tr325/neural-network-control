                                                                     
                                                                     
                                                                     
                                             
test:
	g++ test.cpp -lblas -llapack -o test

util_test.o: util_test.cpp utilities.h
	g++ -lblas -llapack -c util_test.cpp

utilities.o: utilities.cpp utilities.h 
	g++ -lblas -llapack -c utilities.cpp

util_test: util_test.o utilities.o
	g++ util_test.o utilities.o -lblas -llapack -o utest
