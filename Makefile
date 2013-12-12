                                                                     
                                                                     
                                                                     

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

dale_test.o: dale_test.cpp utilities.h optimise.h
	g++ -lblas -llapack -c dale_test.cpp
	
dale_test: dale_test.o utilities.o optimise.o
	g++ dale_test.o utilities.o optimise.o -lblas -llapack -o dtest

repar_test.o: repar_test.cpp utilities.h optimise.h
	g++ -lblas -llapack -c repar_test.cpp

repar_test: repar_test.o utilities.o optimise.o
	g++ repar_test.o utilities.o optimise.o -lblas -llapack -o reptest
	
reform_Syn_test.o: reform_Syn_test.cpp utilities.h optimise.h
	g++ -lblas -llapack -c reform_Syn_test.cpp

refsyn_test: reform_Syn_test.o utilities.o optimise.o
	g++ reform_Syn_test.o utilities.o optimise.o -lblas -llapack -o refsyntest	
	
generateW.o: generateW.cpp generateW.h
	g++ -c generateW.cpp
	
gen_test.o: gen_test.cpp utilities.h generateW.h
	g++ -c gen_test.cpp
	
gen_test: generateW.o utilities.o
	g++ gen_test.o generateW.o utilities.o -lblas -llapack -o gentest


