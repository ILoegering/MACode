function [original_data,header] = load_b32(varargin)
%load_b32 - open B-mode analog data files from ultrasound
%Modified version of convert2DICOM.m written by Barry Belmont
%Source: https://github.com/barrybelmont/ultrasound
%
% Syntax:  [original_data,header] = load_b32(infile,inpath)
%
% Inputs:
%    infile (optional) - char array (varies)
%           Filename of input .b32 file with extension
%    inpath (optional) - char array (varies)
%           Path where .b32 file is located
%
% Outputs:
%    original_data - double matrix (numHeightPix x numWidthPix x numFrames)
%           Matrix of B-mode ultrasound image data
%    header - struct (1 x 1)
%           Header information
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Barry Belmont
% Adapted by: Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: isaacloegering@gmail.com
% May 2018; Last revision: 21-May-2019
%------------- BEGIN CODE --------------
if (nargin == 1)
    filename = varargin{1};
elseif (nargin == 2)
    filename = [varargin{2} '\' varargin{1} '.b32'];
else
    [filename, pathname] = uigetfile('*.b32', ...
        'Choose files to work with', pwd, ...
        'MultiSelect', 'off');
    
    addpath(genpath(pathname));
    cd(pathname)
end

fid = fopen(filename, 'r');
fileExt = filename(end-3:end);

% Read header info,atopm
hinfo = fread(fid, 19, 'int32');

% Load  header information into a structure and save as a separate file
header = struct('filetype', 0, 'nframes', 0, 'w', 0, 'h', 0, 'ss', 0,...
    'ul', [0,0], 'ur', [0,0], 'br', [0,0], 'bl', [0,0], 'probe', 0,...
    'txf', 0, 'sf', 0, 'dr', 0, 'ld', 0, 'extra', 0);
header.filetype = hinfo(1);
header.nframes = hinfo(2);
header.w = hinfo(3);
header.h = hinfo(4);
header.ss = hinfo(5);
header.ul = [hinfo(6), hinfo(7)];
header.ur = [hinfo(8), hinfo(9)];
header.br = [hinfo(10), hinfo(11)];
header.bl = [hinfo(12), hinfo(13)];
header.probe = hinfo(14);
header.txf = hinfo(15);
header.sf = hinfo(16);
header.dr = hinfo(17);
header.ld = hinfo(18);
header.extra = hinfo(19);

% Determine which kind of file type is being used
switch header.filetype
    case {2, 4, 8, 16, 64, 256, 1024, 4096, 8192, 16384, 32768}
        original_data = zeros(header.h, header.w, header.nframes);
    case {32}
        original_data = zeros(header.h, header.nframes);
    case {128}
        original_data = zeros(header.h, 1, header.nframes);
    case {512}
        original_data = zeros(header.h, header.w * header.extra, header.nframes);
    case {2048}
        original_data = zeros(header.h, header.w * 2, header.nframes);
    case {524288}
        original_data = zeros(header.w, header.nframes);
    otherwise
        original_data = [];
end

% Load the data and display progress to the user
h = waitbar(0, 'Please wait while the data are being loaded.');
for frame_count = 1:header.nframes
    
    if(header.filetype == 2) %.bpr
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end % Each frame has 4 byte header for frame number in older versions
        [v,count] = fread(fid,header.w*header.h,'uchar=>uchar');
        original_data(:,:,frame_count) = reshape(v,header.h,header.w);
        
    elseif(header.filetype == 4) %postscan B .b8
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end
        [v,count] = fread(fid,header.w*header.h,'uint8');
        temp = int16(reshape(v,header.w,header.h));
        original_data(:,:,frame_count) = imrotate(temp, -90);
        
    elseif(header.filetype == 8) %postscan B .b32
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end
        [v,count] = fread(fid,header.w*header.h,'uint32');
        temp = reshape(v,header.w,header.h);
        original_data(:,:,frame_count) = imrotate(temp, -90);
        
    elseif(header.filetype == 16) %rf
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end
        [v,count] = fread(fid,header.w*header.h,'int16');
        original_data(:,:,frame_count) = int16(reshape(v,header.h,header.w));
        
    elseif(header.filetype == 32) %.mpr
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end
        [v,count] = fread(fid,header.h,'int16');
        original_data(:,frame_count) = v;
        
    elseif(header.filetype == 64) %.m
        [v,count] = fread(fid,'uint8');
        temp = reshape(v,header.w,header.h);
        original_data(:,:,frame_count) = imrotate(temp,-90);
        
    elseif(header.filetype == 128) %.drf
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end
        [v,count] = fread(fid,header.h,'int16');
        original_data(:,:,frame_count) = int16(reshape(v,header.h,1));
        
    elseif(header.filetype == 512) %crf
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end
        [v,count] = fread(fid,header.extra*header.w*header.h,'int16');
        original_data(:,:,frame_count) = reshape(v,header.h,header.w*header.extra);
        % to obtain data per packet size use:
        % Im(:,:,:,frame_count) = reshape(v,header.h,header.w,header.extra);
        
    elseif(header.filetype == 256) %.pw
        [v,count] = fread(fid,'uint8');
        temp = reshape(v,header.w,header.h);
        original_data(:,:,frame_count) = imrotate(temp,-90);
        
        
    elseif(header.filetype == 1024) %.col
        [v,count] = fread(fid,header.w*header.h,'int');
        temp = reshape(v,header.w,header.h);
        temp2 = imrotate(temp, -90);
        original_data(:,:,frame_count) = mirror(temp2,header.w);
        
    elseif((header.filetype == 2048)) %color .cvv (the new format as of SONIX version 3.1X)
        % velocity data
        [v,count] = fread(fid,header.w*header.h,'uint8');
        temp = reshape(v,header.w,header.h);
        temp2 = imrotate(temp, -90);
        tempIm1 = mirror(temp2,header.w);
        
        % sigma
        [v,count] =fread(fid, header.w*header.h,'uint8');
        temp = reshape(v,header.w, header.h);
        temp2 = imrotate(temp, -90);
        tempIm2 = mirror(temp2,header.w);
        
        original_data(:,:,frame_count) = [tempIm1 tempIm2];
        
    elseif(header.filetype == 4096) %color vel
        if (regexp(version, '5.*') == 1) tag = fread(fid,1,'int32'); end
        [v,count] = fread(fid,header.w*header.h,'uchar=>uchar');
        temp = reshape(v,header.w,header.h);
        temp2 = imrotate(temp, -90);
        original_data(:,:,frame_count) = mirror(temp2,header.w);
        
    elseif(header.filetype == 8192) %.el
        [v,count] = fread(fid,header.w*header.h,'int32');
        temp = reshape(v,header.w,header.h);
        temp2 = imrotate(temp, -90);
        original_data(:,:,frame_count) = mirror(temp2,header.w);
        
    elseif(header.filetype == 16384) %.elo
        [v,count] = fread(fid,header.w*header.h,'uchar=>uchar');
        temp = int16(reshape(v,header.w,header.h));
        temp2 = imrotate(temp, -90);
        original_data(:,:,frame_count) = mirror(temp2,header.w);
        
    elseif(header.filetype == 32768) %.epr
        [v,count] = fread(fid,header.w*header.h,'uchar=>uchar');
        original_data(:,:,frame_count) = int16(reshape(v,header.h,header.w));
        
    elseif(header.filetype == 65536) %.ecg
        [v,count] = fread(fid,header.w*header.h,'uchar=>uchar');
        original_data = v;
        
    elseif or(header.filetype == 131072, header.filetype == 262144) %.gps
        gps_posx(:, frame_count) =  fread(fid, 1, 'double');   %8 bytes for double
        gps_posy(:, frame_count) = fread(fid, 1, 'double');
        gps_posz(:, frame_count) = fread(fid, 1, 'double');
        if strcmp(version, '6.0.3')
            gps_a(:, frame_count) =  fread(fid, 1, 'double');
            gps_e(:, frame_count) = fread(fid, 1, 'double');
            gps_r(:, frame_count) = fread(fid, 1, 'double');
        end
        gps_s(:,:, frame_count) = fread(fid, 9, 'double');    %position matrix
        if strcmp(version, '6.0.3')
            gps_q(:, frame_count) =  fread(fid, 4, 'double');
        end
        gps_time(:, frame_count) = fread(fid, 1, 'double');
        gps_quality(:, frame_count) = fread(fid, 1, 'ushort'); % 2 bytes for unsigned short
        Zeros(:, frame_count) =  fread(fid, 3, 'uint16');     % 6 bytes of padded zeroes
        
    elseif (header.filetype == 524288) % time stamps
        LineStartInClkCycles = fread(fid, header.w, 'uint32');
        original_data(:, frame_count) = LineStartInClkCycles;
    else
        disp('Data not supported');
    end
    
    if (ishandle(h))
        waitbar(frame_count/header.nframes, h);
    else
        h = waitbar(frame_count/header.nframes, 'Please wait while the data is being loaded.');
    end
end

if or(header.filetype == 131072, header.filetype == 262144) %.gps
    original_data.gps_posx    = gps_posx;
    original_data.gps_posy    = gps_posy;
    original_data.gps_posz    = gps_posz;
    original_data.gps_s       = gps_s;
    original_data.gps_time    = gps_time;
    original_data.gps_quality = gps_quality;
    original_data.gps_Zeros   = Zeros;
    original_data.gps_a    = gps_a;
    original_data.gps_e    = gps_e;
    original_data.gps_r    = gps_r;
    original_data.gps_q    = gps_q;
    
end

pause(0.1);
if (ishandle(h))
    delete(h);
end

fclose(fid);

end