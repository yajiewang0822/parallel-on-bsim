
# Parallel BSIM Reconstruction Simulation

---

## URL
[https://yajiewang0822.github.io/parallel-on-bsim/](https://yajiewang0822.github.io/parallel-on-bsim/)
[link to the final report](https://github.com/yajiewang0822/parallel-on-bsim/blob/main/report.pdf)
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
### Paper Refernce 
#### Information about BSIM 
1. Ströhl, Florian and C. Kaminski. “Frontiers in structured illumination microscopy.” (2016).
2. Yeh, Li-Hao et al. “Structured illumination microscopy with unknown patterns and a statistical prior.” Biomedical optics express 8 2 (2017): 695-711 .

#### Parallelze Stochastic Gradient Descent
1. Zinkevich, Martin et al. “Parallelized Stochastic Gradient Descent.” NIPS (2010).
2. Zinkevich, Martin et al. “Slow Learners are Fast.” NIPS (2009).

## Goals and Deliverables 
We plan to parallelize BSIM part of the simulation. If time allows, we may try to parallelize the pattern generation part of the simulation using GPU as it can be a highly data-independent part. 
During the poster session, depending on the result, we may be able to have a live comparison between the sequential and the parallel version of the code by executing them. However, so far, we do not have C++ code yet, and thus we do not know how much time the sequentail version would require. Our best guess for poster session would be that we will simply show the results as a table/graph. 

## Platform Choice 
We would like to use C++ and Halide to code our program running in Latedays. 
We are still deciding whether we use OpenMP or MPI to parallelize our code. We would like to choose wither OpenMP or MPI as we have a better understanding of how the parallel process should look like as we read the related paper.

 ## <s>Schedule</s> 
- <s>Week1(11.5  - 11.11): Learn to use Halide and read about paper to parallelize SGD and start port code from MATLAB to C++</s>
- <s>Week2(11.12 - 11.18): Finish code in C++</s>
- <s>Week3(11.19 - 11.25): Work on parallelizing BSIM </s>
- <s>Week4(11.26 - 12.2):  Continue work and debug the parallel BSIM </s>
- <s>Week5(12.3  - 12.9):  Collect result, write reports, and prepare the poster sessoon. </s>

# Parallel on BSIM Checkpoint

### Timeline
- 11.7-11.9 Write sequential code for the first part of BSIM: Simulate microscope imaging process to get all the input images(Yajie) and related helper functions (Peicheng) 
- 11.10-11.12 Try applying Halide for convolution (Both). 
- 11.13-11.15 Move on writing sequential code for the second part of BSIM: Reconstruct high-resolution image from all the input images (Both).
- 11.16-11.21 Change our plan to use FFTW(Yajie) and Eigen(Peicheng) for convolution and matrix operations respectively. 
- 11.22-11.24 Apply opencv to the project(Peicheng). Debug the input images based on opencv(Yajie).
- 11.25-11.27 Debug the reconstruction result(Both). 
- 11.27-11.30 Starting parallelizing BSIM
- (expected)12.1- 12.7 Finish parallel version of BSIM 
- (expected) 12.7-12.11 Debug and improve the parallel version 
- (expected) 12.12-12.14 Collect results, finish report and prepare the poster session 

### Completed work
We finished porting the Matlab code to the sequential C++ code of BSIM as our baseline. There were some changes of plan we originally had. Instead of using Halide, we decided to use both FFTW and Eigen libraries in our implementation, which involved major changes of data types in our code.  The reason we initially planned to use Halide was that it has FFT related functions that we could take advantage of. However, when we did some research about Halide, we found that Halide may not be so good with large image convolution in CPU[1]. In addition, the source code for Matlab doing FFT used the FFTW library. We believed that using the same library can better compare our C++ implementation to the Matlab version of the code. Using Halide in GPU might be a “nice-to have” topic as we found that Halide in GPU had better performance[2]. However, we probably would not have sufficient time to do this. During the period, we also thought about how to parallelize the proximal gradient descent portion of our code. 

### Goals and Deliverables 
Comparing our current work status to the proposal, we are roughly one week behind the schedule due to change of plan from Halide to FFTW and Eigen. We spent about half a week learning about Halide but ended up not using it. We now have a more clear goal to achieve in the end. Our sequential C++ version needs around 30 minutes to run the entire BSIM while the Matlab version only takes about 7 minutes. We would try to make our parallel version have the speedup comparable to the Matlab version with 8 cores. In other words, around 4x speedup with 8 cores for parallel version.  

### Plan for poster session
In the poster session, we would show how the resolution is improved based on current paralleled code. We would highlight the improvement on performance and execution time in terms of bar charts and/or tables. We plan not to have a live demo due to the time limitation. As for now, our sequential version needs around 30 minutes to generate the result image. We anticipate less time with our parallel version, but we believe that the time required to run both versions and compare would be too much for the poster session. 

### Reference
[1] https://cacm.acm.org/magazines/2018/1/223877-halide/fulltext?mobile=false
[2] https://dvicini.github.io/projects/bscthesis/vicini_bsc_thesis.pdf

### Concerns
We read about the Bridges machine system configuration and found out that they have multiple CPU nodes available. We need to read the Bridges user guide to see how we know which nodes we can use and how to connect to the specific nodes if possible. 
We probably also need to figure out how MPI and OpenMP would assign jobs in the Bridges machine. 
In our message passing process, we may need to communicate the image, which might be a huge data transfer and slow down the program. We need to test it when we finish the parallel version. 

