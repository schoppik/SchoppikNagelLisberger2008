function r = myconv(A,B);

% prepare for the Fourier transforms
nfft = 2^(1+nextpow2(size(A,1)));

% to get rid of the edges from the xcorr/fft
ftind = nfft/2-floor(size(A,1)/2)+1:nfft/2+round(size(A,1)/2);

% transform and reshape
r = fftshift(real(ifft(fft(A,nfft).*conj(fft(B,nfft)))),1);

% resize
r = r(ftind,:,:);
