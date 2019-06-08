% compilation of src files
% for unlocking mex run 'class_interface_mex' without arguments

% Specify path to mkl:
setenv('MKL_PATH_INCLUDE','/home/beams/VNIKITIN/anaconda/include') 
setenv('MKL_PATH_LIB','/home/beams/VNIKITIN/anaconda/lib') 
mex -v -output mfiles/class_interface_mex -I$MKL_PATH_INCLUDE -L$MKL_PATH_LIB -liomp5 -lmkl_rt src/class_interface_mex.cpp src/expradon.cpp src/convs.cpp src/gridproc.cpp 