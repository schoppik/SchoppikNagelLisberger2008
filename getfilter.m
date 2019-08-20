%
%
%
%
%              by: david schoppik
%            date: 8/18/06
%         purpose: to return a decorrelated filter (by diving through by
%                  the appropriate power spectrum in the Fourier domain)
% 
%           usage: for the "reverse" filter, in units of deg/s/spike, call:
%                  getfilter(spikes,eye)
%                  for the "forward" filter, in units of spike/deg/s, call:
%                  getfilter(spikes,eye,cutoff) where cutoff is the index
%                  corresponding to the appropriate cutoff frequency.  The
%                  filter is an exponential with a fixed decay of 5 indices.
%
%            N.B.: the units on these filters should be correct, but they assume
%                  that the global mean of both spikes and eye is zero.
%                  Similarly, any effects on the scale of the filtering of
%                  the "forward" filter aren't accounted for.
%                  Also, the filters returned are of meaningful length (tau
%                  goes from -signal_length/2:signal_length/2 so that it
%                  returns a filter as long as the original signal.

function r = getfilter(A,B,cutoff)

if nargin <= 1, 
  help getfilter
  return;
end
  
if nargin == 2 % divide through by the power spectrum of the spike train
  
  % prepare for the Fourier transforms
  nfft = 2^(1+nextpow2(size(A,1)));
  % to get rid of the edges from the xcorr/fft
  ftind = nfft/2-floor(size(A,1)/2)+1:nfft/2+round(size(A,1)/2);
  
  % cross correlate
  % the conj is around A to flip the filter correctly
  sta = real(ifft(fft(A,nfft).*conj(fft(B,nfft))));
  
  % get the power spectrum
  ps = mean(fft(A,nfft).*conj(fft(A,nfft)),2)+eps;
  
  % decorrelate by division
  dc = fftshift(real(ifft(fft(mean(sta,2),nfft)./ps)),1);
  
  % return the piece we want
  r = dc(ftind,:,:);
end

if nargin == 3 % divide through by the power spectrum of the eye movement
  
  % prepare for the Fourier transforms
  nfft = 2^(1+nextpow2(size(A,1)));
  % to get rid of the edges from the xcorr/fft
  ftind = nfft/2-floor(size(A,1)/2)+1:nfft/2+round(size(A,1)/2);
  
  % cross correlate
  
  % the conj is around B to flip the filter correctly
  sta = mean(real(ifft(fft(B,nfft).*conj(fft(A,nfft)))),2);
  ps = mean(fft(B,nfft).*conj(fft(B,nfft)),2)+eps;
  
  % make a filter to smooth out the data
  tau = 5;
  expfilt = ones(1,nfft/2);
  expfilt(cutoff+1:end) = exp(-[1:nfft/2-cutoff]/tau);
  expfilt = [expfilt fliplr(expfilt)]';
  expfilt = repmat(expfilt,[1 1 size(A,3)]);
    
  dc = fftshift(real(ifft(fft(sta,nfft)./ps.*expfilt)));
   
  r = dc(ftind,:,:);
  
end