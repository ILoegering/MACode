function [trc,time,f,n,nmrk,mrk_names,file,inpath]=load_trc_wVirtual(infile,inpath)
%   [pos,time,f,n,nmrk,mrk_names]=load_trc_wVirtual(infile,inpath)
%   alternately
%   trc=load_trc_wVirtual(infile)
%
%   LOAD_TRC_WVIRTUAL is used to open a data file from Motion Analysis Realtime
%   output (*.trc).
%
%   Inputs:
%       infile - trc file to be loaded
%                If infile is unspecified, the user is prompted to select the input file
%       inpath - directory of location where data file is located
%               when no path is specified, it defaults to current directory
%
%   Outputs:
%       pos     contains - the meaured marker positions in order of the markers
%               that is columns 1-3 are the x,y,z components of marker 1
%                       columns 4-6 are the x,y,z components of marker 2
%                          ....
%       time - column vector of time
%       f - sample frequency
%       n - number of data frames
%       nmrk - number of markers
%       mrk_names - marker names
%
%   Updated: Feb. 15, 2005 (JWF)
%   Updated: Oct. 12, 2006 (David Remy, remy@wisc.edu)
%   -> Initialized data fields for better performance
%   -> Enable handling of missing markers (position values are set to [NaN,
%      NaN, NaN] if a marker is missing).
%   -> Cleaned up structure
%
%   MATLAB Version 7.1

narg = nargin;
if (narg==0);
    [infile, inpath]=uigetfile('*.trc','Select input file');
    if (infile==0)
        disp('No file selected');
        f = 0;
        n = 0;
        nmrk = 0;
        mrk_names = {};
        pos = [];
        time = [];
        return;
    end
    fid = fopen([inpath infile],'r');
    file = infile(1:length(infile));
elseif (narg==1)
    file = [infile(1:length(infile)) '.trc'];
    fid=fopen(file,'r');
    inpath='';
else
    file = [infile(1:length(infile)) '.trc'];
%     [inpath infile]
    fid=fopen([inpath file],'r');
end

if (fid==-1)
    disp('File not found');
    trc=[];
    return;
end

%disregard header info (first two lines)
fgetl(fid);
fgetl(fid);

% scan the next line for the file info:
line=fgetl(fid);
file_info=sscanf(line,'%f');
f   = file_info(1);     % data rate
n   = file_info(3);     % nr frames
nmrk= file_info(4);     % nr markers

% the next line contains the marker names:
line=fgetl(fid);
i = 1;
name =  sscanf(line,[repmat('%*s',1,2),' %s'],1);
while (~isempty(name))
    mrk_names(i,1) = cellstr(name);
    name = sscanf(line,[repmat('%*s',1,i+2),' %s'],1);
    i = i + 1;
end
vm_count = i-1-nmrk;
nmrk = i-1;

% disregard the next two lines
fgetl(fid);
fgetl(fid);

% initialize the data variables
pos = nan(n,3*nmrk);
time = zeros(n,1);
i=0;
numColConfirmed = 0;
while feof(fid)==0
    i=i+1;
    line=fgetl(fid);
    data=sscanf(line,'%f');
    numCol = size(find(line==9),2);
    % read time from the second column
    time(i,1)=data(2);
    if (~numColConfirmed)
        % check if there exists virtual marker data
        if (numCol==3*nmrk+2)
            disp(['Data for ' num2str(nmrk-vm_count) ' physical markers and ' num2str(vm_count) ' virtual markers were loaded.'])
        % check if just physical marker data
        elseif (numCol==3*(nmrk-vm_count)+2)
            nmrk = nmrk-vm_count;
            mrk_names = mrk_names(1:nmrk,1);
            disp([num2str(vm_count) ' virtual markers were assigned during collection, but no data were collected.'])
            pos = nan(n,3*nmrk);
        else
            error('There was a problem loading the data. The number of data columns does not match the number of markers. This may be because there is only a partial set of virtual marker data. This function will need to be modified to handle this case.')
        end
        numColConfirmed = 1;
    end
    % check if the full set of markers is given
    if size(data,1)==3*nmrk+2
        % in this case, just shift and copy the data from the string to the
        % marker positions:
        pos(i,:) = data(3:end)';
    else
        % if a line misses some markers, the values of these markers are set to NaN:
        % Find the tabs (we have to do this manually)
        tabIndex = find(line==9);
        % scan the line, tab by tab
        for k=1:numel(tabIndex)-2
            % The first pice of data is between the second and third tab,
            % the second is between the third and fourth...
            dataPice = sscanf(line(tabIndex(k+1):tabIndex(k+2)),'%f');
            if isempty(dataPice)
                % No data between these two tabs -> marker is missing.
                pos(i,k) = NaN;
            else
                pos(i,k) = dataPice;
            end
        end
    end
end

% Return the position data in m
pos = pos/1000;

if (nargout>1)
    % return trc as the matrix of kinematic data
    trc=pos;
else
    % return all information in a single structure
    trc.pos=zeros(size(pos,1),3,nmrk);
    ii=1:3;
    for i=1:nmrk
        trc.pos(:,:,i)=pos(:,ii);
        ii=ii+3;
    end
    trc.time=time;
    trc.freq=f;
    trc.nframes=n;
    trc.nmrk=nmrk;
    trc.mrk_names=mrk_names;
    trc.file=file;
    trc.inpath=inpath;
end
