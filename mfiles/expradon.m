function [R]=expradon(cid,f,Ns,Ntheta)
f=single(f);
theta=(-Ntheta/2:Ntheta/2-1)/Ntheta*2*pi;
rho=(-Ns/2:Ns/2-1)'/Ns;
eta1=-cid.mu*rho.^0*sin(theta)-2*pi*1i*rho*cos(theta);eta1=single(eta1);x=-imag(eta1)/(2*pi);ax=real(eta1);
eta2= cid.mu*rho.^0*cos(theta)-2*pi*1i*rho*sin(theta);eta2=single(eta2);y=-imag(eta2)/(2*pi);ay=real(eta2);
set_grids(cid,x,y,ax,ay,int32([Ns Ntheta Ns]));

g=eq2us(cid,f);%g=g(1:2:end,:)+1i*g(2:2:end,:);R=real(fftshift(ifft(ifftshift(g))));
R=expifft1d(cid,g);
end