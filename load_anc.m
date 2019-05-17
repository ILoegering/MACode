function [anc,time,samp_rate,channel_names,range,voltdata]=load_anc(infile,inpath)
%load_anc - open analog data files from Motion Analysis Realtime output (*.anc)
%The number of channels, samp rate, and range are determined from the file
%header.
%
% Syntax:  [anc,time,samp_rate,channel_names,range,voltdata] = load_anc(infile,inpath)
%
% Inputs:
%    infile (optional) - char array (varies)
%           Filename of .anc file to load. If the file extension '.anc' is
%           not included in infile, then it is added by the function. If
%           unspecified, the file may be selected via a file browser.
%    inpath (optional) - char array (varies)
%           Path to folder in which file is located. If unspecified, the
%           current directory is used.
%
% Outputs:
%    anc - double matrix (numSamples x numChannels) -OR- struct (1 x 1)
%           Matrix of kinematic data if only one output argument.
%           Otherwise, a structure containing all data and header info from
%           .anc file
%    samp_rate - double array (1 x numChannels)
%           Sample rate for each channel from header info
%    channel_names - cell array (numChannels x 1)
%           Names of channels from header info
%    range - double array (1 x numChannels)
%           Range for each channel from header info
%    voltdata - double matrix (numSamples x numChannels)
%           Data converted to volts
%           
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Amy Silder
% Editors: Darryl Thelen, Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: dgthelen@wisc.edu
% October 2005; Last revision: 17-May-2019
%------------- BEGIN CODE --------------
% Load .anc file
% If filename unspecified, prompt user with file browser to select file
if (nargin==0)
    [infile, inpath]=uigetfile('*.anc','Select input file');
    if infile==0
        disp('No file selected');
        return;
    end
    fid = fopen([inpath '\' infile],'r');
% If filename specified, ensure it ends in '.anc' and load
elseif (nargin>0)
    if ~strcmp(infile((end-3):end),'.anc') 
        infile = [infile,'.anc'];
    end
    if (nargin==2)
        fid = fopen([inpath '\' infile],'r');
    else 
        inpath='';
        fid=fopen(infile,'r');
    end
end
% If file was not found, alert user
if (fid==-1)
    disp(['File not found: ',infile]);
    anc=[];
    return
end

% Disregard header info
for h=1:2
    fgetl(fid);
end
% Get Duration
fscanf(fid,'%s',5); % Skip characters
duration=fscanf(fid,'%f',1);
% Get number of channels
fscanf(fid,'%s',1); % Skip characters
num_channels=fscanf(fid,'%f',1);
% Disregard header info and blank lines
for h=1:5
    fgetl(fid);
end
% Get channel names
fscanf(fid,'%s',1); % Skip characters
channel_name=fgetl(fid);    % Line with channel names
j=1;
jl=length(channel_name);
channel_names = cell(num_channels,1);
for i=1:(num_channels)
    name=sscanf(channel_name(j:jl),'%s',1);
    ii=strfind(channel_name(j:jl),name);
    j=j+ii(1)+length(name);
    channel_names(i,1)=cellstr(name);
end
% Get sample rate for each channel
fscanf(fid,'%s',1); % Skip characters
for i=1:num_channels
    samp_rate(i)=fscanf(fid,'%f',1);
end
% Get range for each channel
fscanf(fid,'%s',1);
for i=1:num_channels
    range(i)=fscanf(fid,'%f',1);
end
% Get time and data
npts=round(duration*samp_rate(1)+1);
time=zeros(npts,1);
data=zeros(npts,num_channels);
i=1;
linetemp=(fscanf(fid,'%f',num_channels+1))';
while ((feof(fid)==0)&&(~isempty(linetemp)))
    time(i,1)=linetemp(1);
    data(i,1:num_channels)=linetemp(2:(num_channels+1));
    i=i+1;
    linetemp=(fscanf(fid,'%f',num_channels+1))';
end

fclose(fid);

% Convert to volts
% 12-bit AD board
% twelvebitConv=(2^16)/2;
% voltdata=(ones(size(data,1),1)*range).*(data/twelvebitConv)*0.001;
% 16-bit AD board
sixteenbitConv=(2^16)/2;
voltdata=(ones(size(data,1),1)*range).*(data/sixteenbitConv)*0.001;

if (nargout>1)
    % return trc as the matrix of kinematic data
    anc=data;
else
    % return all information in a single structure
    anc.data=voltdata;
    anc.addata=data;
    anc.time=time;
    anc.freq=samp_rate;
    anc.nframes=size(data,1);
    anc.nchan=num_channels;
    anc.channel_names=channel_names;
    anc.range=range;
    anc.file=infile;
    anc.inpath=inpath;
end
end
%------------- END OF CODE --------------