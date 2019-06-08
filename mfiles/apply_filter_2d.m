function [g,filter,t] = apply_filter_simple2d(f,os,filter,d,ellipse)
N=size(f,1);
Ne=os*(N-1);
fe=zeros(Ne);
fe(Ne/2+1+(-floor(N/2):ceil(N/2)-1),Ne/2+1+(-floor(N/2):ceil(N/2)-1))=f;

fhate=fftshift(fft2(ifftshift(fe)));

% fhate2=zeros(size(fhate));
% for k=1:size(ellipse,1)
% fhate2=fhate2+fphantom(os,Ne,ellipse(k,:))';
% end
% keyboard
%fhate=fhate2;

[t1,t2]=ndgrid((0:(Ne/2))/Ne,(0:(Ne/2))/Ne);
t=sqrt(t1.^2+t2.^2);
switch filter
    case 'ram-lak'
        wfa=Ne*0.5*t/Ne;%.*(t/(2*d)<=1);%compute the weigths
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        wfa = sinc(t/(2*d)).*(t/d<=2*sqrt(2));
    case 'cosine'
        wfa = cos(pi*t./(2*d)).*(t/d<=1);
    case 'cosine2'
        wfa = (cos(pi*t./(2*d))).^2.*(t/d<=1);
    case 'hamming'
        wfa = (.54 + .46 * cos(pi*t./d)).*(t/d<=1);
    case 'hann'
        wfa=(1+cos(pi*t./d)) / 2.*(t/d<=1);
    otherwise
        error(message('images:iradon:invalidFilter'))
end
%   case 'shepp-logan'
% %         be careful not to divide by 0:
%         filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
%     case 'cosine'
%         filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
%     case 'hamming'
%         filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
%     case 'hann'
%         filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
%     otherwise
%         error(message('images:iradon:invalidFilter'))
wfa=wfa.*(wfa>=0);%.*(t/d<=1);
% wfamid=2*wfa(:,1);
% keyboard
wfa=[flipud(fliplr(wfa(2:end,:))) flipud(wfa(2:end,2:end));fliplr(wfa(:,:)) wfa(:,2:end)];
% wfa=[flipud(fliplr(wfa(2:end,:))) flipud(wfa(2:end,2:end));fliplr(wfa(:,:)) wfa(:,2:end)];
wfa=wfa(1:end-1,1:end-1);

g=fftshift(ifft2(ifftshift(wfa.*fhate)));
g=g(Ne/2+1+(-floor(N/2):ceil(N/2)-1),Ne/2+1+(-floor(N/2):ceil(N/2)-1));
g=g;
filter=wfa;
t=[fliplr(-t(2:end)) 0 t(2:end-1)]';
end
