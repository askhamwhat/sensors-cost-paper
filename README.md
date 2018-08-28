# sensors-cost-paper

This repository contains the software companion
to the paper "Greedy Sensor Placement With Cost Constraints"
[preprint on arXiv](https://arxiv.org/abs/1805.03717).

## How to use

The code takes in samples of data as the rows of
a matrix (so that columns correspond to sensor
locations) and a vector representing the costs
at each location. It is then possible to trace
out a cost-error curve by scaling the input cost
vector for the algorithm (this changes its relative
strength in the algorithm).

There is a simple demonstration of the code and
algorithm in simple_qrpc_example.m

## Figure generation

For reproducibility, we have included here the
codes used to generate most of the figures in the
preprint. These codes take a very long time to run
(especially using the pure matlab implementation
of the code, see below).

### Obtaining the data

The data for the various experiments can be
obtained as follows.

#### Yale Faces

You can download the relevant file for the
Yale faces data with
[this dropbox link](https://www.dropbox.com/s/vp1pl8jriy5twzf/YaleB_32x32.mat?dl=0)

#### Sea Surface Temperature

This data is hosted by NOAA. ADD DETAILS

#### Cylinder flow

This data is included with the companion to the
Dynamic Mode Decomposition book [here](http://dmdbook.com/).
ADD DETAILS

## Faster Option With MEX

We've included a MEX wrapper to a FORTRAN routine.
The codes seem to compile and run well on UNIX systems
with standard gnu compilers. These faster codes are
basically unsupported for now but worth a shot.

These can be compiled by changing directory to
the src directory and running build_xqrmc_m

Be sure to then test the code with test_xqrmc_m.
If a zero is output, things should be good.

If you find that this has broken everything, simply
delete the MEX binary to restore the slower operation.