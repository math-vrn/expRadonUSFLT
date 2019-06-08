clear all
addpath('mfiles');
N=2^9;Ntheta=ceil(pi*N)+3;Ns=N;mu=log(100)/N;Nthreads=8;

%filtered phantom
f=phantom(N);
sigma=1.2;
[x1,x2]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));circ0=(sqrt(x1.^2+x2.^2)<1-8/N)*1.0;
gauss=1/sigma^2/(2*pi)*exp(-1/2*((x1*N/2/sigma).^2+(x2*N/2/sigma).^2));
ff=fftshift(ifft2(ifftshift( fftshift(fft2(ifftshift(f))).*fftshift(fft2(ifftshift(gauss)))        )));
ff=ff.*circ0;ff=ff/max(abs(ff(:)));


%create class
cid=init_expradon(N,mu,Nthreads);

%take 180 degrees of data 
v=ones(Ns,1)*[zeros(1,Ntheta/4) ones(1,Ntheta/2) zeros(1,Ntheta/4)];v=v/sum(v(1,:))*Ntheta;
R=expradon(cid,ff',Ns,Ntheta).*v;
fr=expradon_inv(cid,R);

% Compute point spread function
psf=zeros(N);psf(N/2+1,N/2+1)=1;
Rp=expradon(cid,psf',Ns,Ntheta).*v;
delete(cid);
cid=init_expradon(2*N,mu,Nthreads);
psf=expradon_inv(cid,Rp);
delete(cid);

P=single(fftshift(fft2(fftshift(psf))));
%Solve deconvolution problem
frc=Toeplitz_apply(fr,conj(P),circ0);
K=@(f)reshape( Toeplitz_apply(Toeplitz_apply(reshape(f,N,N),P,circ0),conj(P),circ0),N^2,1);
frec=reshape(pcg(K,frc(:),1e-5,50),N,N);



% plots
%plots
subplot(2,2,1);imagesc(ff);colormap(gray);title('f');
subplot(2,2,2);imagesc(R');ylabel('\theta');xlabel('s');colormap(gray);title('R\mu f');
subplot(2,2,3);imagesc(frec);colormap(gray);title('rec');
subplot(2,2,4);imagesc(abs(ff-frec));colormap(gray);title('error');





