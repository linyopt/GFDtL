# Group Fused Graphical Lasso (GFGL)

Accompanying MATLAB code for *The Group Fused Graphical Lasso for Estimation of Piecewise Stationary Gaussian Graphical Models*. This package contains the basic Alternating Directed Method of Multipliers algorithm, alongside files to simulate data and demonstrate recovery of graphical structures.

### Requirements

This package depends on a couple of libraries for efficiently computing projection operators:

1) Within GFGL, for solving the group lasso projection we utilise a block-coordinate descent routine developed by Bleakley et. al. *"The group fused Lasso for multiple change-point detection"*. 

http://cbio.ensmp.fr/~jvert/svn/GFLseg/html/

2) For the solving the independent Fused Graphical Lasso (IFGL) we utilise the efficient routines described in the paper by Liu et. al. *"An Efficient Algorithm for a Class of Fused Lasso Problems", 2010*. Such routines for sparse learning are conviniently packaged in the *Sparse Learning with Efficient Projections (SLEP)* package.

http://www.yelab.net/software/SLEP/

*Files are downloaded in the absence of required libraries being specified on the MATLAB path


### Acknowledgements (2015)

Primary author and code
Alex Gibberd - UCL Department of Statistical Science

Co-author and supervisor
Dr James Nelson - UCL Department of Statistical Science

Funding provided by the Defence Science Technology Laboratory (DSTL)


## INSTALATION

1) Extract the constents of this folder into a directory included on the MATLAB path

2) Run *install.m* to attempt to automatically download dependencies

NOTE: The instalation procedure has been tested in Mac OSX and Linux (Redhat) but relies 
on the use of wget, tar and unzip commands. This may not work automatically in Windows.
The required packages (SLEP and GFLseg) may be required to be installed by hand and 
can be found at the addresses given in "Requirements".

## EXAMPLES

We give two examples of GFGL and IFGL applied to simulated data. The first "normalDemo()" demonstrates 
recovery of graphical structure in the standard T>P setting; the second "hdDemo()" looks at the 
high-dimensional case. Hyper-parameters can be specified through setting lambda1 
(for sparsity), or lambda2 (for smoothness), note I/G refer to parameters for IFGL/GFGL respectively. 

### Standard setting

> \>>normalDemo(lambda1G,lambda2G,lambda1I,lambda2I)

Estimation of a graph with size P=5, T=30, with n=5 true edges. Data is from zero-mean gaussian with dynamic correlation structure specfied by the graph and a single changepoint located at cp=T/2.

*Plot 1* - presents the estimated graph within each segment alongside the recovered estimates for GFGL and IFGL.

*Plot 2* - plots the estimates for the active ground-truth edges (highlighting changepoints) as a function of t=1,..,T. The grouping propety of GFGL in this setting.

### High-dimensional setting

> \>>hdDemo(lambda1G,lambda2G,lambda1I,lambda2I)

Same as above but in *high dimensional* (P>T) setting, with P=20,T=10, and n=5 true edges



