clear
close all
clc

%% Collection info
isCortex = 1;  % 0 = Optotrak, 1 = Cortex
% Equipment: Ultrasound
tw = 0.038; % Transducer width in m
td = 0.020; % Transducer depth in m
% Subject/trial-specific info
subjectID = 'TMD_TD01';
kinemID = 'Ex'; % For ankle, use '20' to indicate 20 deg knee angle. For knee, use 'Ex' or 'Fl' for extension or flexion.
trial = 'K2'; % For ankle, use 'A#'. For knee, use 'K#'.
angSpacing = 5; % Angle spacing: spacing between angles over which to measure MA
testDate = ''; % Only used for pilot Optotrak data (ankle). Empty string if Cortex.
% Unique identifier used to save analysis files
save_name = [subjectID '_' kinemID '_' trial];
if (isCortex==0)
    save_name = [save_name '_' trialDate];
end
% Check td
if ((contains(trial,'A') && td == 0.020) || (contains(trial,'K') && td == 0.040))
    error('Check tendon depth variable td.  Achilles is usually 0.040, and knee is usually 0.020.')
end
% Equipment: MoCap
if (contains(trial,'A'))
    oldMrkNames = {'R.Knee','R.SH1','R.SH2','R.SH3','R.SH4','R.Ankle',...
        'R.Mkn','R.Man','R.Ft1','R.Ft2','R.Ft3','US1','US2','US3','LL',...
        'UL','LR','UR','VMUL','VMLL','VMLR','VMUR'}; % MoCap marker names given during collection
    newMrkNames = {'lat_con','prox1','prox2','prox3','prox4','lat_mal',...
        'med_con','med_mal','dist1','dist2','dist3','trans1','trans2',...
        'trans3','LL','UL','LR','UR','w1','w2','w3','w4'}; % MoCap marker names used in processing (see organizeMocap4MA.m)
elseif (contains(trial,'K'))
    oldMrkNames = {'R.TH1','R.TH2','R.TH3','R.TH4','R.Knee','R.Mkn',...
        'R.SH1','R.SH2','R.SH3','R.SH4','R.Ankle','R.Man','US1','US2',...
        'US3','LL','UL','LR','UR','VMUL','VMLL','VMLR','VMUR'}; % MoCap marker names given during collection
    newMrkNames = {'prox1','prox2','prox3','prox4','lat_epi','med_epi',...
        'dist1','dist2','dist3','dist4','lat_mal','med_mal','trans1',...
        'trans2','trans3','LL','UL','LR','UR','w1','w2','w3','w4'}; % MoCap marker names used in processing (see organizeMocap4MA.m)
else
    error('Unrecognized trial type')
end

% Data file paths and filenames
root = 'G:\My Drive\UW NMBL\Tendon\Tendon Mechanics Database\TMD_TDchildren\TMD Subject Data\';
static_path = [root subjectID '\Collected Data\Motion Capture Data\Processed_US\'];
static_filename = [subjectID '_PatCalib_1'];
mot_path = [root subjectID '\Collected Data\Motion Capture Data\Processed_US\'];
mot_filename = save_name;
us_path = [root subjectID '\Collected Data\US\'];
us_filename = save_name;
anc_path = [root subjectID '\Collected Data\Motion Capture Data\Processed_US\Generated_files\'];
anc_filename = save_name;
plots_path = [root subjectID '\Data Analysis\MA\Plots\'];
data_path = [root subjectID '\Data Analysis\MA\Data\'];
% Save flags
save_data = 1;
save_plots = 1;

%% Load data
% Load and restructure mocap data
static = load_trc_wVirtual(static_filename, static_path);
mot = load_trc_wVirtual(mot_filename, mot_path);
[static_pos, static_info] = organizeMocap4MA(static,oldMrkNames,newMrkNames);
[mot_pos, mot_info] = organizeMocap4MA(mot,oldMrkNames,newMrkNames);
% Load ultrasound data and trim border
[us_data, us_header] = RPread(strcat(us_path,[us_filename '.b32']), Inf);
xCollapse = squeeze(sum(us_data,[2,3])); % Find extrema along x direction
xmin = find(xCollapse>0,1,'first'); xmax = find(xCollapse>0,1,'last');
yCollapse = squeeze(sum(us_data,[1,3])); % Find extrema along y direction
ymin = find(yCollapse>0,1,'first'); ymax = find(yCollapse>0,1,'last');
wVox = xmax-xmin+1; % Calculate number of voxels along width dimension
dVox = ymax-ymin+1; % Calculate number of voxels along depth dimension
us_data = us_data(xmin:xmax,ymin:ymax,:);   % Trim border from data
us_data = flipud(us_data);    % Flip y-axis direction (reflect data across x-axis so origin is in lower left corner when YDir is normal)

%% Check data
% DATA CHECK: Window markers labeled correctly?
if (0)
    figure
    title('Check transducer window markers.')
    hold on
    pntrange = 1:190;
    % Left window - red/magenta
    w1 = getMrkPos(mot_pos,'w1');
    w2 = getMrkPos(mot_pos,'w2');
    plotMrk(w1,'r.',pntrange);
    plotMrk(w2,'m.',pntrange);
    % Right window - blue/cerulean
    w3 = getMrkPos(mot_pos,'w3');
    w4 = getMrkPos(mot_pos,'w4');
    plotMrk(w3,'b.',pntrange);
    plotMrk(w4,'c.',pntrange);
    % Proximal
    % Find indices of proximal and distal markers
    proxIdx = find(contains({static_pos.mrk_name},'prox')==1);
    distIdx = find(contains({static_pos.mrk_name},'dist')==1);
    for i = 1:numel(proxIdx)
        prox = getMrkPos(mot_pos,mot_pos(proxIdx(i)).mrk_name);
        h1 = plotMrk(prox,'g.',pntrange);
    end
    % Distal
    for i = 1:numel(distIdx)
        dist = getMrkPos(mot_pos,mot_pos(distIdx(i)).mrk_name);
        h2 = plotMrk(dist,'y.',pntrange);
    end
    hold off
    axis equal
    view(0,0)
    answer = questdlg('Are the markers labeled correctly?','Check Markers','Yes','No','No');
    switch(answer)
        case 'Yes'
        case 'No'
            error('Markers are incorrectly labeled')
        case ''
            error('Markers are incorrectly labeled')
    end
end

% DATA CHECK: Fourier Transform of mot data
if(1)
    figure
    hold on
    % For each marker
    for i = 1:numel(mot_pos)
        % For each coordinate
        for j = 1:size(mot_pos(i).pos,2)
            plotFFT(mot_info.time,mot_pos(i).pos(:,j));
        end
    end
    hold off
end

%% Filter data
% Define filter
fc = 6;         % Cut-off frequency
fs = mot_info.freq;  % Sampling frequency
filtOrder = 6;  % Order of filter
[b,a] = butter(filtOrder,fc/(fs/2),'low');
% Filter STATIC data
for i=1:numel(static_pos)
    % DATA CHECK: Plot unfiltered data
    if(0)
        f = figure;
        subplot(3,1,1)
        hold on
        plot(static_info.time,static_pos(i).pos(:,1),'b.')   % x-coord
        xlabel('time (sec)')
        ylabel('x (m)')
        subplot(3,1,2)
        hold on
        plot(static.time,static_pos(i).pos(:,2),'b.')   % y-coord
        xlabel('time (sec)')
        ylabel('y (m)')
        subplot(3,1,3)
        hold on
        plot(static.time,static_pos(i).pos(:,3),'b.')   % z-coord
        xlabel('time (sec)')
        ylabel('z (m)')
    end
    % Remove dropout frames
    [cf_static,idxNaN_static] = removeNaN(static_pos(i).pos);
    % If all frames blank, skip
    if (isempty(idxNaN_static))
        try
            close f
        catch
        end
        continue
    end
    % Filter complete data
    filteredData_static = filtfilt(b,a,cf_static);
    % Reinsert dropout frames
    restoredData_static = insertNaN(filteredData_static,idxNaN_static);
    static_pos(i).pos = restoredData_static;
    % DATA CHECK: Plot filtered data
    if(0)
        subplot(3,1,1)
        plot(static_info.time,static_pos(i).pos(:,1),'r.')   % x-coord
        legend('raw','filtered')
        hold off
        subplot(3,1,2)
        plot(static_info.time,static_pos(i).pos(:,2),'r.')   % y-coord
        legend('raw','filtered')
        hold off
        subplot(3,1,3)
        plot(static_info.time,static_pos(i).pos(:,3),'r.')   % z-coord
        legend('raw','filtered')
        hold off
        sgtitle(f,static_pos(i).mrk_name,'Interpreter','none')
    end
end
% Filter MOTION data
for i=1:numel(mot_pos)
    % DATA CHECK: Plot unfiltered data
    if(0)
        f = figure;
        subplot(3,1,1)
        hold on
        plot(mot_info.time,mot_pos(i).pos(:,1),'b.')   % x-coord
        xlabel('time (sec)')
        ylabel('x (m)')
        subplot(3,1,2)
        hold on
        plot(mot.time,mot_pos(i).pos(:,2),'b.')   % y-coord
        xlabel('time (sec)')
        ylabel('y (m)')
        subplot(3,1,3)
        hold on
        plot(mot.time,mot_pos(i).pos(:,3),'b.')   % z-coord
        xlabel('time (sec)')
        ylabel('z (m)')
    end
    % Remove dropout frames
    [cf_mot,idxNaN_mot] = removeNaN(mot_pos(i).pos);
    % If all frames blank, skip
    if (isempty(idxNaN_mot))
        try
            close f
        catch
        end
        continue
    end
    % Filter complete data
    filteredData_mot = filtfilt(b,a,cf_mot);
    % Reinsert dropout frames
    restoredData_mot = insertNaN(filteredData_mot,idxNaN_mot);
    mot_pos(i).pos = restoredData_mot;
    % DATA CHECK: Plot filtered data
    if(0)
        subplot(3,1,1)
        plot(mot_info.time,mot_pos(i).pos(:,1),'r.')   % x-coord
        legend('raw','filtered')
        hold off
        subplot(3,1,2)
        plot(mot_info.time,mot_pos(i).pos(:,2),'r.')   % y-coord
        legend('raw','filtered')
        hold off
        subplot(3,1,3)
        plot(mot_info.time,mot_pos(i).pos(:,3),'r.')   % z-coord
        legend('raw','filtered')
        hold off
        sgtitle(f,mot_pos(i).mrk_name,'Interpreter','none')
    end
end

%% Downsample and synchronize data
if (isCortex)
    % Downsample and synchronize mocap to US data
    dsampleRate = 10;
    % Calibration data - just downsample since no US data
    staticInd = 1:dsampleRate:1 + (static_info.nframes./dsampleRate - 1)*dsampleRate;
    static_pos = selectFrames(static_pos,staticInd);
    static_info.time = static_info.time(1:dsampleRate:end);
    static_info.freq = static_info.freq/dsampleRate;
    static_info.nframes = static_info.nframes/dsampleRate;
    % Motion data - downsample and synchronize
    if(~exist([us_path '\motIndStart_' trial '.mat'],'file'))
        anc = load_anc(anc_filename,anc_path);
        trig = anc.data(:,13);
        trig = trig./max(trig);
        trigThreshold = 0.1;
        motIndStart = find(diff(trig) > trigThreshold,1);
        % Check indStart
        f = figure;
        hold on
        plot(trig,'b-')
        plot([motIndStart motIndStart],[min(trig),max(trig)],'r-')
        hold off
        xlim([motIndStart-100, motIndStart+100])
        title('Check start index')
        answer = questdlg('Correct start index?','Check start index','Yes','No','Yes');
        close(f)
        if(strcmp(answer,'No'))
            error('Please adjust trigger threshold and try again.')
        end
        motIndStart = roundn(motIndStart,2)/dsampleRate;
        save([us_path '\motIndStart_' trial],'motIndStart')
    else
        load([us_path '\motIndStart_' trial '.mat'],'motIndStart')
    end
    % If numFramesUS hasn't been calculated, then load US data and calculate.
    % Otherwise, load numFramesUS from saved .mat file.
    if(~exist([us_path '\numFramesUS_' trial '.mat'],'file'))
        US = load_b32(us_filename,us_path);
        numFrameUS = size(US,3);
        save([us_path '\numFramesUS_' trial],'numFrameUS')
    else
        load([us_path '\numFramesUS_' trial '.mat'],'numFrameUS')
    end
    motIndEnd = motIndStart + (numFrameUS - 1)*dsampleRate;
    motInd = motIndStart:dsampleRate:motIndEnd;
    % Check if mocap ended before US
    if (mean(motInd < mot_info.nframes-(dsampleRate-1)) < 1)
        disp('The motion capture concluded before the ultrasound.')
    end
    motInd = motInd(motInd < mot_info.nframes-(dsampleRate-1));
    mot_pos = selectFrames(mot_pos,motInd);
    mot_info.time = mot_info.time(1:dsampleRate:dsampleRate*size(mot_pos(1).pos,1));
    mot_info.freq = mot_info.freq/dsampleRate;
    mot_info.nframes = size(mot_pos(1).pos,1);
end

%% Calculate moment arms
[angle, ma, tendonCoords] = getMomentArm(static_pos, static_info, mot_pos, mot_info, angSpacing, us_data, us_header, tw, wVox, td, dVox, plots_path, save_name);
if (save_data)
    if (7~=exist(data_path,'dir'))
        mkdir(data_path);
    end
    save([data_path '\' save_name],'angle','ma','tendonCoords')
end
% Plot moment arm vs. joint angle
f = figure;
plot(angle,ma,'o')
title('Moment Arm vs. Joint Angle')
xlabel('Joint Angle (deg)')
ylabel('Moment Arm (mm)')
if (save_plots)
    if (7~=exist([plots_path '\Moment Arm\tiff\'],'dir'))
        mkdir([plots_path '\Moment Arm\tiff\']);
    end
    if (7~=exist([plots_path '\Moment Arm\fig\'],'dir'))
        mkdir([plots_path '\Moment Arm\fig\']);
    end
    saveas(f, [plots_path '\Moment Arm\tiff\' save_name],'tiff')
    saveas(f, [plots_path '\Moment Arm\fig\' save_name],'fig')
end