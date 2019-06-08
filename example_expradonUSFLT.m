clear all
addpath('mfiles');

N=2^9;Ntheta=ceil(pi)*N+3;Ns=N;mu=log(1)/N;Nthreads=8;
%create filtered phantom
[x1,x2]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));circ0=(sqrt(x1.^2+x2.^2)<1-4/N)*1.0;
sigma=1.2;
gauss=1/sigma^2/(2*pi)*exp(-1/2*((x1*N/2/sigma).^2+(x2*N/2/sigma).^2));
f=phantom(N);
ff=fftshift(ifft2(ifftshift( fftshift(fft2(ifftshift(f))).*fftshift(fft2(ifftshift(gauss)))        )));
ff=ff.*circ0;ff=ff/max(abs(ff(:)));
imagesc(ff);

%create class
cid=init_expradon(N,mu,Nthreads);
%expRadon
R=expradon(cid,ff',Ns,Ntheta);
%inversion
frec=expradon_inv(cid,R);
delete(cid);
%plots
subplot(2,2,1);imagesc(ff);colormap(gray);title('f');
subplot(2,2,2);imagesc(R');ylabel('\theta');xlabel('s');colormap(gray);title('R\mu f');
subplot(2,2,3);imagesc(frec);colormap(gray);title('rec');
subplot(2,2,4);imagesc(abs(ff-frec));colormap(gray);title('error');

