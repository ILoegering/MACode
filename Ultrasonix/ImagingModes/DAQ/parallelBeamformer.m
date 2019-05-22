function out = parallelBeamformer(in, maxAprSz, Fs)
% The following script take the pre-beamform data and creates beamfromed
% image by applying parallel beamformers to the same data. 
%
% input:
%   in:             the first dimension should be scan lines and the second
%                   dimension should be the channel samples
%   maxAprSz:       size of the aperture beamformer
%   Fs              sampling frequency 40*1e6 or 80*1e6
% output:
%   out: beamformed data
%
% Copyright Ultrasonix Medical corporation - Analogic Ultrasound Group
% Author: Reza Zahiri Azar - Nov 2013

% extract the parameters
fnumber = 2.0;      % rx fnumber used for beamforming
pitch = 300*1e-6;       % spacing between lines/channels in the lateral direction [m]
c = 1540;               % speed of sound 1540 [m/s]
sampleSpacing = c/Fs/2;         % spacing between samples in the axial direction in milimeter for 40MHz
numOfLine = size(in,2);         % number of lines
numOfSample = size(in,1);       % number of samples

% beamforming (the following steps can be speed up by precalculating the delays)
out = zeros(numOfSample, numOfLine);
h = waitbar(0,'Running parallel beamforming, please wait...');
for i = 1:numOfLine    % for each line
    waitbar( i / (numOfLine) );
    for j=1:numOfSample  % find the value for each sample
        % find the sample location
        Z = (j) * sampleSpacing ; 
        X = (i) * pitch;
        
        % calculate the aperture based on the F number
        a = Z/(2*fnumber);
        hlfAprSz = floor(a / pitch);
        
        % chk to make sure we do not go beyound max apr
        if (hlfAprSz > maxAprSz/2)
            hlfAprSz = floor(maxAprSz / 2);
        end
        if (hlfAprSz < 1)
            hlfAprSz = 1;
        end
        
        % aperture indices
        x = -hlfAprSz: hlfAprSz;    
        fulAprSz = 2*hlfAprSz;

        % find the aperture
        cntr = i;       % center of apreture
        apr = cntr + x; % full aperture indices

        % calc delays based on sample depth and receive aperture
        X1 = (cntr + x) * pitch;
        delays = ( Z + sqrt( Z^2 + (X-X1).^2 ) ) / c * Fs; 
        delays = round(delays); % no interpolation used in this version
        chnumOfLines = zeros(1,fulAprSz);  % will be filled with proper signal from each channel
        
        % find the corresponding delayed values from each channel
        for k = 1:fulAprSz  
            % find channel index
            % ignore the values if the fall outside
            chlIndx = apr(k);   
            if chlIndx<1, continue, end;
            if chlIndx>numOfLine, continue, end;
            % find sample index
            % ignore the values if the fall outside
            chlSmpl = delays(k);    
            if chlSmpl<1, continue, end;
            if chlSmpl>numOfSample, continue, end;
            % add the delays value to the list
            chnumOfLines(k) = in(chlSmpl, chlIndx); 
        end;
        
        % apodization by windowing
        win = hanning(fulAprSz)';
        chnumOfLines = win .* chnumOfLines;
        
        % standard summation
        out(j,i) = sum( chnumOfLines );

    end;
end
close(h);
