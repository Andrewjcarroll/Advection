This is an RK4 method to solve the advection equation. It will solve the differential equation and then output the reults to a gnuplot and save every 25 steps to a plots folder. 
The Filtering does not currently work. Currently two types of derivatives offered an explict finite schgeme (2nd order) for testing purposes. These can be set in the dif.h 
There are also three initial_data.h types. The current types offered are SIN a simple sin wave.  KIM which is the intial data from Kim's paper and GAUSSIAN which is a simple Guassian pulse. 
To build simply type make. To rebuild make clean.  