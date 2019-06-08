function [frec]=expradon_inv(cid,g)
g=single(g);
%zero-padding
g=[zeros(3*size(g,1)/2,size(g,2));g;zeros(3*size(g,1)/2,size(g,2));];
[Ns,Ntheta]=size(g);

theta=(-Ntheta/2:Ntheta/2-1)/Ntheta*2*pi;
rho=(-Ns/2:Ns/2-1)'/Ns;
%take weights for composite quadrature rule corrections
w=take_weights(rho,cid.mu);w=single(w);
id=find(w);rho=rho(id);w=w(id);%Reduce the number of sampling points in the Laplace domain.
%polar coordinates
eta1= cid.mu*rho.^0*sin(theta)+2*pi*1i*rho*cos(theta);eta1=single(eta1);x=-imag(eta1)/(2*pi);ax=real(eta1);
eta2=-cid.mu*rho.^0*cos(theta)+2*pi*1i*rho*sin(theta);eta2=single(eta2);y=-imag(eta2)/(2*pi);ay=real(eta2);
set_grids(cid,x,y,ax,ay,int32([numel(rho) Ntheta Ns]));

gt=expfft1d(cid,g);%gt=gt(1:2:end,:)+1i*gt(2:2:end,:);
frec=us2eq(cid,gt,w,int32(id-1));
[x1,x2]=meshgrid(linspace(-1,1,cid.N),linspace(-1,1,cid.N));circ0=(sqrt(x1.^2+x2.^2)<1-4/cid.N)*1.0;
frec=frec'*(theta(2)-theta(1)).*circ0;