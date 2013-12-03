                                                                     
                                                                     
                                                                     
                                             
test:
	g++ test.cpp -lblas -llapack -o test

util_test.o: util_test.cpp utilities.h
	g++ -lblas -llapack -c util_test.cpp

utilities.o: utilities.cpp utilities.h 
	g++ -lblas -llapack -c utilities.cpp

util_test: util_test.o utilities.o
	g++ util_test.o utilities.o -lblas -llapack -o utest

optimise.o: optimise.cpp optimise.h utilities.h
	g++ -lblas -llapack -c optimise.cpp

opt_test.o: opt_test.cpp optimise.h
	g++ -lblas -llapack -c opt_test.cpp

opt_test: utilities.o optimise.o opt_test.o 
	g++ opt_test.o optimise.o utilities.o -lblas -llapack -o optest
	
gradMat_num_test.o: gradMat_num_test.cpp utilities.h
	g++ -lblas -llapack -c gradMat_num_test.cpp
	
grad_test: gradMat_num_test.o utilities.o optimise.o
	g++ gradMat_num_test.o optimise.o utilities.o -lblas -llapack -o gradtest

MatMult_test.o: MatMult_test.cpp utilities.h
	g++ -lblas -llapack -c MatMult_test.cpp

MatMult_test: MatMult_test.o utilities.o
	g++ MatMult_test.o utilities.o -lblas -llapack -o mmtest
	
tracetest.o: tracetest.cpp utilities.h
	g++ -lblas -llapack -c tracetest.cpp
	
trace_test: trace_test.o utilities.o
	g++ tracetest.o utilities.o -lblas -llapack -o ttest








