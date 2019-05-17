function [data,time,f,n,nmrk,mrk_names,file,inpath]=load_opto(infile,inpath)
%   [pos,time,f,n,nmrk,mrk_names]=load_opto(infile,inpath)
%   alternately
%   data=load_opto(infile)
%
%   LOAD_OPTO is used to open an Optotrak data file exported from First
%   Principles. The data file must be formatted as an ASCII
%   comma-separated values (CSV) file.
%
%   NOTE: When data is exported from First Principles, two data files are
%   created: collectionName_3d.csv and collectionName_6d.csv. Use the 3d
%   file.
%
%   Inputs:
%       infile - csv file to be loaded
%                If infile is unspecified, the user is prompted to select the input file
%       inpath - directory of location where data file is located
%               when no path is specified, it defaults to current directory
%
%   Outputs:
%       pos     contains - the meaured marker positions in order of the markers
%               that is columns 1-3 are the x,y,z components of marker 1
%                       columns 4-6 are the x,y,z components of marker 2
%                          ....
%       time - column vector of time (sec)
%       f - sample frequency (Hz)
%       n - number of data frames
%       nmrk - number of markers
%       mrk_names - marker names
%
%   Adapted from load_trc.m by Isaac Loegering (iloegering@wisc.edu)
%
%   MATLAB R2017a

narg = nargin;
if (narg==0)
    [infile, inpath]=uigetfile('*.csv','Select input file');
    if infile==0
        disp('No file selected');
        f = 0;
        n = 0;
        nmrk = 0;
        mrk_names = {};
        pos = [];
        time = [];
        return;
    end
    fid = fopen([inpath '\' infile],'r');
    file = infile(1:length(infile));
elseif (narg==1)
    file = [infile(1:length(infile)) '.csv'];
    fid=fopen(file,'r');
    inpath='';
else
    file = [infile(1:length(infile)) '.csv'];
    [inpath '\' infile];
    fid=fopen([inpath '\' file],'r');
end

if (fid==-1)
    disp('File not found');
    data=[];
    return;
end

% get number of frames from first line
line=fgetl(fid);
if contains(line, ',')
    n=strsplit(line,',');
    line=n{1};
end
n=strsplit(line,': ');
n=str2num(n{2});

% get frame rate from second line
line=fgetl(fid);
if contains(line, ',')
    f=strsplit(line,',');
    line=f{1};
end
f=strsplit(line,': ');
f=str2num(f{2});
% skip "Units" line and blank line
fgetl(fid);
fgetl(fid);
% get number of markers and marker names
line=fgetl(fid);
mrk_names=strsplit(line,',');
mrk_names=mrk_names(2:3:end);
mrk_names=mrk_names';
nmrk=size(mrk_names,1);
for i=1:nmrk
    mrk_names{i,1}=mrk_names{i,1}(1:end-2);
end

% initialize the data variables
pos = zeros(n,3*nmrk);
time = zeros(n,1);
i=1;

while feof(fid)==0
    line=fgetl(fid);
    tokens=strsplit(line,',','CollapseDelimiters',0);
    tokens=tokens(2:end);
    % calculate time point from frame rate
    time(i,1)=(i-1)*(1/f);
    % if a line misses some markers, the values of these markers are set to NaN:
    for k=1:3*nmrk
        if isempty(tokens{1,k})
            pos(i,k) = NaN;
        else
            pos(i,k) = str2double(tokens{1,k});
        end
    end
    i=i+1;
end

% Return the position data in m
pos = pos/1000;

if (nargout>1)
    % return trc as the matrix of kinematic data
    data=pos;
else
    % return all information in a single structure
    data.pos=zeros(size(pos,1),3,nmrk);
    ii=1:3;
    for i=1:nmrk
        data.pos(:,:,i)=pos(:,ii);
        ii=ii+3;
    end
    data.time=time;
    data.freq=f;
    data.nframes=n;
    data.nmrk=nmrk;
    data.mrk_names=mrk_names;
    data.file=file;
    data.inpath=inpath;
end
