# expRadonUSFLT
Evaluation of the exponential Radon transform so as its inversion by using Unequally-spaced fast Laplace transforms.

## Compilation
Execute file compile.m with specifying paths to MKL libraries, e.g.
setenv('MKL_PATH_INCLUDE','/home/beams/VNIKITIN/anaconda/include') 
setenv('MKL_PATH_LIB','/home/beams/VNIKITIN/anaconda/lib') 

## Run 
Example scripts example_expRadonUSFLT.m and example_expRadonIters.m  should be run in Matlab with having MKL libraries in LD_LIBRARY_PATH, e.g.

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/beams0/VNIKITIN/anaconda/lib matlab &

## Reference
Andersson, Fredrik, Marcus Carlsson, and Viktor V. Nikitin. "Fast Laplace transforms for the exponential Radon transform." Journal of Fourier Analysis and Applications 24.2 (2018): 431-450.

