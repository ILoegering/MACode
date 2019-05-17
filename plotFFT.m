function [h]= plotFFT(x,y)
%plotFFT - plot fast Fourier transform of raw signal data
%If data are real, only the positive half of the FFT is plotted since the
%FFT is symmetric about the vertical axis.
%
% Syntax:  [h] = plotFFT(x,y,real)
%
% Inputs:
%    x (required) - double array (numSamples x 1)
%           x data
%    y (required) - double array (numSamples x 1)
%           y data
%
% Outputs:
%    h - plot handle (1 x 1)
%           Handle to FFT plot
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: isaacloegering@gmail.com
% November 2018; Last revision: 14-May-2019
%------------- BEGIN CODE --------------
% Ensure data are arranged in column vectors
if(size(x,2)>1)
    x = x';
end
if(size(y,2)>1)
    y = y';
end
% Perform FFT on signal
y_hat = fftshift(fft(y-mean(y)));
numData = size(y_hat,1);
% Calculate range of x for FFT plotting purposes
xRange = max(x)-min(x);
% Create properly scaled x_hat vector
if(mod(numData,2)==1)
    x_hat = ((ceil(-numData/2):floor(numData/2))./xRange)';
else
    x_hat = ((ceil(-numData/2):floor(numData/2)-1)./xRange)';
end
% Plot positive half of FT (since it's mirrored across vertical access for
% real signals).
if(isreal(y))
    ind = ceil(numData/2):numData;
    h = plot(x_hat(ind),abs(y_hat(ind)));
else
    h = plot(x_hat,abs(y_hat));
end
title('FFT')
xlabel('Frequency (Hz)')
axis tight
end
%------------- END OF CODE --------------