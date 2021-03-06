
## Why doesn't the brain explode?

This was my Master's project at university -- I apologise for the state of the codebase.

The write-up of the project can be found in the `Report` directory.

### Project outline

The project simulated complex neuronal dynamics, observing the stabilisation of large networks 
of neurons using different models of inhibitory action. I then analysed connectivity patterns
in the network and compared these to experimental data.

The files containing the core of the stabilisation procedure are:
`utilities.cpp`
`optimise.cpp`
`generateW.cpp`

These together generate an initial synaptic weight matrix, and then stabilise the matrix using an iterative technique,
varying the inhibitory connection locations and strengths to match the excitatory network. 

`*.m` files contain scripts for analysing the output of trials involving several runs of the stabilisation procedure. 
`*_test.cpp` files were used for testing of the functions in the core `*.cpp` files.

Please read the final report for the project (in `Report/FinalReport.pdf`) if you want more detail on the 
methods, results and conclusions of the investigation.
