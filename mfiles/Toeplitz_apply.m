function h=Toeplitz_apply(f,P,circ)
N=size(f,1);
fe=zeros(2*N,2*N);fe(N+1+(-N/2:N/2-1),N+1+(-N/2:N/2-1))=f.*circ;
he=fftshift(ifft2(fftshift(  fftshift(fft2(fftshift(fe))).*P)));
h=he(N+1+(-N/2:N/2-1),N+1+(-N/2:N/2-1)).*circ;


