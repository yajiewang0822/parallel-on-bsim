
# Parallel BSIM Reconstruction Simulation
## URL
https://yajiewang0822.github.io/parallel-on-bsim/
## Summary
We would like to speed up the Reconstruction simulation by parallelizing the Blind Structured Illumination Microscopy(BSIM) algorithm which involves proximal griadient descent.
## Background
Due to the Abbe diffraction limitation, imaging resolution for microscopy is limited to twice of the Numerical aperture. Therefore, methods for improving resolution has been widely researched. 
BSIM method, due to its high resolution result, low excitation density and wide-field imaging, has been widely used in medical, biology and material fields. Research people have to wait for long time, probably several hours, to obtain the result but no efforts, as far as we know, have been made on speeding up the simulation because research people tend to care more about the result instead of the execution time. Therefore, parallel on BSIM method is an interesting and valuable idea to do some researches. 
The steps for BSIM reconstruction simulation is:
1. Simulate the microscopic imaging process based on BSIM, including generating objective, illumination patterns and imagings. 
2. Based on simulation results, utilizing proximal gradient method to generate esimated illumination patterns.
3. Reconstrct the final result based on the covariance relationship between estimated illumination patterns and real imagings.

## The Challenge
1. How to parallelize gradient descent, an inherently sequential operation. 
2. With multicore setting, how should we reduce the communication cost/resource contention? 
3. The locality may be a problem since there are a lot of computation on large 2D array. Therefore, how to partition the work to reduce the cache miss is also a challenge we will face.

## Resources 
We have the MATLAB version of the entire simulation process. We will build our C++ code based on the MATLAB code.
### Paper Refernce (Change to citation form later)
#### Information about BSIM 
1. https://www.osapublishing.org/DirectPDFAccess/93B2A204-E697-E8E4-C2CE227D772C14A8_344732/optica-3-6-667.pdf?da=1&id=344732&seq=0&mobile=no
2. https://www.osapublishing.org/DirectPDFAccess/93B8A8BD-EE08-D223-AB9B0DE179A910A0_357203/boe-8-2-695.pdf?da=1&id=357203&seq=0&mobile=no

#### Parallelze Stochastic Gradient Descent
1. http://martin.zinkevich.org/publications/nips2010.pdf
2. https://arxiv.org/pdf/0911.0491.pdf
# Goals and Deliverables 
We plan to parallelize BSIM part of the simulation. If time allows, we may try to parallelize the pattern generation part of the simulation using GPU as it can be a highly data-independent part. 
During the poster session, depending on the result, we may be able to have a live comparison between the sequential and the parallel version of the code by executing them. However, so far, we do not have C++ code yet, and thus we do not know how much time the sequentail version would require. Our best guess for poster session would be that we will simply show the results as a table/graph. 

## Platform Choice 
We would like to use C++ and Halide to code our program running in Latedays. 
We are still deciding whether we use OpenMP or MPI to parallelize our code. We would like to choose wither OpenMP or MPI as we have a better understanding of how the parallel process should look like as we read the related paper.

## Schedule 
- Week1(11.5  - 11.11): Learn to use Halide and read about paper to parallelize SGD and start port code from MATLAB to C++
- Week2(11.12 - 11.18): Finish code in C++
- Week3(11.19 - 11.25): Work on parallelizing BSIM
- Week4(11.26 - 12.2):  Continue work and debug the parallel BSIM
- Week5(12.3  - 12.9):  Collect result, write reports, and prepare the poster sessoon. 
