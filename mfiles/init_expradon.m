function cid=init_expradon(N,mu,Nthreads)

epsilon=1e-5;A=exp(mu*N);nu=2;
cid=class_expRadon_matlab.getInstance(int32([N Nthreads nu]),single([A epsilon]));cid.mu=mu;cid.N=N;


