This is an RK4 method to solve the advection equation. It will solve the differential equation and then output the results to a gnuplot and save every 25 steps to a plots folder.
The filtering does not currently work. Currently, two types of derivatives are offered: an explicit finite scheme (2nd order) for testing purposes. These can be set in the dif.h file.
There are also three initial_data.h types. The current types offered are: SIN (a simple sine wave), KIM (which is the initial data from Kim's paper), and GAUSSIAN (which is a simple Gaussian pulse).
To build, simply type make. To rebuild, type make clean.




