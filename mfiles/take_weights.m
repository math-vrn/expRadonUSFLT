function w = take_weights(rho,mu)
Ns=numel(rho);
%% Composite quadrature rule computation
K=14;
mup=mu/(2*pi);
l_min=min(find(rho>(mup)))-K+1;
l_max=max(find(rho<(mup)))+K-1;
w=zeros(Ns,1);v=zeros(Ns,1);
S=zeros(K);for k=0:K-1,S(:,k+1)=linspace(-1,1,K).^k;end;
for l=l_min:l_max,
  t=linspace(rho(l),rho(l+K-1),K);
  q=(max([t(1) mup])*2-t(1)-t(end))/(t(end)-t(1));
  ws=( (t(end)-t(1))/2*((1.^[2:K+1]-q.^[2:K+1])./(2:K+1))+(t(end)+t(1))/2*((1.^[1:K]-q.^[1:K])./(1:K)))*inv(S)*Ns/(K-1)*(t(end)-t(1))/2;
  for k=1:K, w(l+k-1)=w(l+k-1)+ws(k);end;
end;
w(l_max+1:end)=abs(rho(l_max+1:end));
w(1)=abs(rho(1));w(2:Ns)=(w(2:Ns)+flipud(w(2:Ns)));
w=w/2;

%% First undersampling
P=2*ceil((l_max-Ns/2-1)/2);
M=10;
v=linspace(0,1,2*M)';v=1-(126*v.^5-420*v.^6+540*v.^7-315*v.^8+70*v.^9);
w(Ns/2+1+P+(1:2:2*M))=w(Ns/2+1+P+(1:2:2*M)).*v(1:2:end);
w(Ns/2+1+P+(2:2:2*M))=w(Ns/2+1+P+(2:2:2*M)).*(2-v(2:2:end));
w(Ns/2+1+P+1+2*M:2:end)=0;
w(Ns/2+1+P+2+2*M:2:end)=2*w(Ns/2+1+P+2+2*M:2:end);
w(Ns/2+1-(P+(1:2:2*M)))=w(Ns/2+1-(P+(1:2:2*M))).*v(1:2:end);
w(Ns/2+1-(P+(2:2:2*M)))=w(Ns/2+1-(P+(2:2:2*M))).*(2-v(2:2:end));
w(Ns/2+1-(P+1+2*M):-2:1)=0;
w(Ns/2+1-(P+2+2*M):-2:1)=2*w(Ns/2+1-(P+2+2*M):-2:1);

%% Second undersampling
P=P+2*M;
M=60;%need larger value second undersampling
v=linspace(0,1,2*M)';v=1-(126*v.^5-420*v.^6+540*v.^7-315*v.^8+70*v.^9);
w(Ns/2+1+P+(2:4:4*M))=w(Ns/2+1+P+(2:4:4*M)).*v(1:2:end);
w(Ns/2+1+P+(4:4:4*M))=w(Ns/2+1+P+(4:4:4*M)).*(2-v(2:2:end));
w(Ns/2+1+P+2+4*M:4:end)=0;
w(Ns/2+1+P+4+4*M:4:end)=2*w(Ns/2+1+P+4+4*M:4:end);
w(Ns/2+1-(P+(2:4:4*M)))=w(Ns/2+1-(P+(2:4:4*M))).*v(1:2:end);
w(Ns/2+1-(P+(4:4:4*M)))=w(Ns/2+1-(P+(4:4:4*M))).*(2-v(2:2:end));
w(Ns/2+1-(P+2+4*M):-4:1)=0;
w(Ns/2+1-(P+4+4*M):-4:1)=2*w(Ns/2+1-(P+4+4*M):-4:1);
