function [angle, ma, tendonCoords] = getMomentArm(static_pos, static_info, mot_pos, mot_info, angSpacing, us_data, us_header, tw, wVox, td, dVox, plots_path, save_name)
%getMomentArm - calculate moment arm using combined ultrasound-mocap method
%Calculates moment arm as the shortest distance between the tendon line of
%action and a function joint axis. The superficial and deep edges of the
%tendon in the ultrasound images must be manually selected so that the
%tendon midline may be calculated as the average of the superficial and
%deep boundaries. The tendon line of action is then calculated as the
%linear best fit to the tendon midline. The functional joint axis is
%calculated as the best-fit screw axis that describes the motion of the
%distal segment with respect to the proximal segment. A static calibration
%trial defines the configuration in which the joint angle is zero; joint
%angles can then be calculated for each frame of the motion trial based on
%relative position of marker triads on the proximal and distal segments.
%
%
% Syntax:  [angle, ma, tendonCoords] = getMomentArm(static_pos, static_info, mot_pos, mot_info, angSpacing, us_data, us_header, tw, wVox, td, dVox, plots_path, save_name)
%
% Inputs:
%    static_pos (required) - struct array (numMrks x 1)
%           output of organizeMocap4MA, containing a struct for each marker
%           for the static calibration trial with fields 'pos' and
%           'mrk_name'
%    static_info (required) - struct (1 x 1)
%           output of organizeMocap4MA, containing non-position data for
%           the static calibration trial
%    mot_pos (required) - struct array (numMrks x 1)
%           output of organizeMocap4MA, containing a struct for each marker
%           for the motion trial with fields 'pos' and 'mrk_name'
%    mot_info (required) - struct (1 x 1)
%           output of organizeMocap4MA, containing non-position data for
%           the motion trial
%    angSpacing (required) - double (1 x 1)
%           Number of degrees of joint motion between samples selected for
%           moment arm calculations.
%    us_data (required) - 3D double matrix (numVoxWidth x numVoxDepth x numUSFrames)
%           ultrasound B32 data
%    us_header (required) - struct (1 x 1)
%           ultrasound header info
%    tw (required) - double
%           width of transducer in meters
%    wVox (required) - double
%           number of voxels along width dimension in ultrasound images
%    td (required) - double
%           depth of transducer in meters
%    dVox (required) - double
%           number of voxels along depth dimension in ultrasound images
%    plots_path (required) - char array
%           absolute path of folder in which plots will be saved
%    save_name (required) - char array
%           unique identifier for this trial (e.g.,
%           subjectID_condition_trial); used for filename of plots to be
%           saved
%                       
%
% Outputs:
%    angle - double array (numAnglesOfInterest x 1)
%           joint angles at which moment arm was calculated
%    ma - double array (numAnglesOfInterest x 1)
%           calculated moment arms
%    tendonCoords - struct ()
%           coordinates of tendon midline in ultrasound image associated
%           with each joint angle for which moment arm was calculated
%
% Other m-files required: functionalJointAxis.m, tendonDepth.m
% Subfunctions: valueSearch
% MAT-files required: none
%
% Author: Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: isaacloegering@gmail.com
% February 2019; Last revision: 14-May-2019
%------------- BEGIN CODE --------------
%% Compute functional joint axis
% Find indices of proximal and distal markers
proxIdx = find(contains({static_pos.mrk_name},'prox')==1);
distIdx = find(contains({static_pos.mrk_name},'dist')==1);
% Restructure for functionalJointAxis
static_prox = nan(static_info.nframes,size(static_pos(1).pos,2),numel(proxIdx));
static_dist = nan(static_info.nframes,size(static_pos(1).pos,2),numel(distIdx));
mot_prox = nan(mot_info.nframes,size(mot_pos(1).pos,2),numel(proxIdx));
mot_dist = nan(mot_info.nframes,size(mot_pos(1).pos,2),numel(distIdx));
for i = 1:numel(proxIdx)
    static_prox(:,:,i) = static_pos(proxIdx(i)).pos;
    mot_prox(:,:,i) = mot_pos(proxIdx(i)).pos;
end
for i = 1:numel(distIdx)
    static_dist(:,:,i) = static_pos(distIdx(i)).pos;
    mot_dist(:,:,i) = mot_pos(distIdx(i)).pos;
end
% Compute the joint center and functional joint axes
[jc,fja,ang,~,~] = functionalJointAxis(static_prox,static_dist,mot_prox,mot_dist);

% Free memory
clearvars static_pos static_info

% DATA CHECK: Plot MoCap data with jc and fja
if (1)
    % Motion data
    figure
    title('Motion')
    hold on
    pntrange = 1:mot_info.nframes;
    % Proximal
    for i = 1:numel(proxIdx)
        prox = getMrkPos(mot_pos,mot_pos(proxIdx(i)).mrk_name);
        h1 = plotMrk(prox,'b.',pntrange);
    end
    % Distal
    for i = 1:numel(distIdx)
        dist = getMrkPos(mot_pos,mot_pos(distIdx(i)).mrk_name);
        h2 = plotMrk(dist,'r.',pntrange);
    end
    % JC
    h3 = plot3(jc(pntrange,1),jc(pntrange,2),jc(pntrange,3),'go');
    % FJA
    f1=jc(pntrange,:)-0.2*fja(pntrange,:);
    f2=jc(pntrange,:)+0.2*fja(pntrange,:);
    for i=1:size(f1,1)
        h4 = plot3([f1(i,1) f2(i,1)],[f1(i,2) f2(i,2)],[f1(i,3) f2(i,3)],'k-');
    end
    hold off
    legend([h1 h2 h3 h4],{'prox','dist','JC','FJA'})
    axis equal
end

%% Select range of data to analyze
% Get angles in degrees
angDeg = ang(:,1)*(180/pi);

% Plot joint angles and reverse sign if needed
fig = figure;
plot(angDeg, 'b-');
title('Reverse sign of joint angles?')
xlabel('Frame')
ylabel('Joint Angle (deg)')
answer = questdlg(['Which angles are positive and negative is a '...
    'matter of convention. Are angles plotted as expected, or '...
    'should the sign of the angles be reversed?'],['Reverse sign '...
    'of joint angles?'],'Reverse','Continue','Quit','Continue');
switch(answer)
    case 'Reverse'
        angDeg = -angDeg;
    case 'Continue'
    otherwise
        error(['There was a problem with the range of joint angle '...
            'data. The moment arm calculation was terminated early.']);
end

% Trim data as needed
plot(angDeg, 'b-');
title('Please select start and end frame')
xlabel('Frame')
ylabel('Joint Angle (deg)')
[x,~] = ginput(2);
if (x(1) < x(2))
    first = floor(x(1));
    last = ceil(x(2));
else 
    first = floor(x(2));
    last = ceil(x(1));
end
if (first <= 0)
    first = 1;
end
if (last > size(angDeg,1))
    last = size(angDeg,1);
end
angDegTrim = angDeg(first:last,:);
[angFrame, angle] = valueSearch(angDegTrim,angSpacing);
angFrame = angFrame + first - 1;    % Shift so that indices match original angle array

% Check frames selected for angles of interest
hold on
plot(first:last,angDegTrim,'b-','LineWidth',2);
for i = 1:numel(angle)
    plot(angFrame(i),angDeg(angFrame(i)),'ko')
    plot(1:15:angFrame(i),angle(i),'k.')
    vertLine = angle(1)-5:2.5:angDeg(angFrame(i));
    plot(repmat(angFrame(i),1,numel(vertLine)),vertLine,'k.')
end
hold off
title('Selected frames')
legend('All data','Analyzed data','Points of interest','Location','Southeast')

% Save Selected Frames figure
if (1)
    if (7~=exist([plots_path 'Selected Frames\tiff\'],'dir'))
        mkdir([plots_path 'Selected Frames\tiff\']);
    end
    if (7~=exist([plots_path 'Selected Frames\fig\'],'dir'))
        mkdir([plots_path 'Selected Frames\fig\']);
    end
    saveas(fig, [plots_path 'Selected Frames\tiff\' save_name],'tiff')
    saveas(fig, [plots_path 'Selected Frames\fig\' save_name],'fig')
end

%% Get tendon line of action and transform into global reference frame
% Get coords of tendon midline from ultrasound images
offset = 0; % Offset applied to indexing; becomes nonzero when saved tendon coords are loaded from previous session
% If code crashed in a previous session and saved tendon coords...
if (exist(['AAA_temp_' save_name '_tendonCoords.mat']))
    % ...load those tendon coords...
    load(['AAA_temp_' save_name '_tendonCoords.mat'])
    % ...and check to see if they correspond to frames currently being
    % processed. If so...
    if (~ismember(ismember(tendonCoords.frame,angFrame),0))
        % ...ask the user if they want to load saved tendon coords.
        answer = questdlg(['It looks like you have tendon coords from '...
            'a previous session that can be loaded to save time. Would '...
            'you like to load saved tendonCoords?'],['Load saved '...
            'tendonCoords?'],'Yes','No','Yes');
        % If yes,...
        if (strcmp(answer,'Yes'))
            % ...change offset to number of frames already processed,...
            offset = numel(tendonCoords.frame);
            % ...remove already processed frames from angFrame so they
            % aren't processed again, and leave tendonCoords in Workspace.
            angFrame = angFrame(~ismember(angFrame,tendonCoords.frame));
        % Otherwise,...
        else
            % ...clear tendonCoords from Workspace and process all frames.
            clear tendonCoords
        end
    % If frames in saved tendon coords don't correspond with frames
    % currently being processed, alert user and continue to process all
    % frames.
    else
        msgbox(['There are tendon coords saved, but they do not '...
            'correspond with the frames you are currently processing.']);
    end
end
% Get tendon coords for each frame in angFrame
try
    for i = 1:numel(angFrame)
        % Get tendon coords for single frame
        output = tendonDepth(us_data, us_header, td/dVox, tw/wVox, tw, angFrame(i));
        % Store tendon coords in tendonCoords structure
        tendonCoords.depth{1,i + offset} = output.depth{1};
        tendonCoords.axialLocation{1,i + offset} = output.axialLocation{1};
        tendonCoords.frame(1,i + offset) = output.frame;
        % Progress indicator
        disp([num2str(numel(tendonCoords.frame)) ' of '...
            num2str(numel(angFrame) + offset) ' processed'])
    end
% If error occurs...
catch
    % ...save tendonCoords to current directory.
    save(['AAA_temp_' save_name '_tendonCoords.mat'],'tendonCoords')
    error(['The tendonDepth.m function crashed unexpectedly, but '...
        'don''t worry! Your progress was saved in '...
        'AAA_temp_' save_name '_tendonCoords.mat, which is located in '...
        'the current directory. Just don''t change your current '...
        'directory and re-run your code. You should be able to load '...
        'your saved progress.'])
end
% Delete temporary tendonCoords file
delete(['AAA_temp_' save_name '_tendonCoords.mat'])

% Extract ultrasound window marker data
w1 = getMrkPos(mot_pos,'w1');
w2 = getMrkPos(mot_pos,'w2');
w3 = getMrkPos(mot_pos,'w3');   % Not used (only three points needed to define plane/create reference frame)
w4 = getMrkPos(mot_pos,'w4');

% Calculate x,y,z directions of ultrasound reference frame
y = (w1 - w2)./vecnorm((w1 - w2)')';    % Along width of transducer
z = (w4 - w2)./vecnorm((w4 - w2)')';    % Along thickness of transducer
x = cross(y,z);                         % Along depth (into body)
z = cross(x,y); % Ensure reference frame is orthogonal

% Calculate ultrasound transformation matrix, Tu, for each frame
Tu = cell(mot_info.nframes,1);
for i=1:mot_info.nframes
    % Calculate 3x3 rotation matrix R
    R = [x(i,:)' y(i,:)' z(i,:)'];
    % Get translation vector v to origin of ultrasound reference frame
    w2w1 = w1(i,:) - w2(i,:);
    dw2w1 = norm(w2w1);
    correction = dw2w1/2 - tw/2; % All units in meters
    v = (w2(i,:) + w4(i,:))./2 + correction*(w2w1/dw2w1);
    % Input R and v into Tu
    Tu{i,1} = [R v'; 0 0 0 1];
end
% If Tu is invalid (contains NaN) for frames of interest, find closest
% Tu that is valid
for i = angFrame
    flag = 1;
    if (sum(isnan(Tu{i,1}),'all'))
        j = 1;
        while (flag)
            Tu_alt = Tu{i-j,1};
            if (~sum(isnan(Tu_alt),'all'))
                Tu{i,1} = Tu_alt;
                flag = 0;
                continue
            end
            Tu_alt = Tu{i+j,1};
            if (~sum(isnan(Tu_alt),'all'))
                Tu{i,1} = Tu_alt;
                flag = 0;
                continue
            end
            j = j + 1;
        end
    end
end

% Transform ultrasound data into global reference frame
numFOI = size(tendonCoords.frame,2);    % Number of frames of interest
gTendonCoords = struct('frame',[],'coords',[]); % Coordinates of tendon in global reference frame
for i=1:numFOI
    gTendonCoords.frame(1,i) = tendonCoords.frame(1,i);
    homoCoords = [tendonCoords.depth{1,i}; tendonCoords.axialLocation{1,i}; zeros(size(tendonCoords.depth{1,i})); ones(size(tendonCoords.depth{1,i}))];
    transHomoCoords = Tu{tendonCoords.frame(1,i)}*homoCoords;
    gTendonCoords.coords{i} = transHomoCoords(1:3,:)';
end

%% Calculate moment arm
Q = zeros(length(gTendonCoords.frame),3);
tendonDir = zeros(length(gTendonCoords.frame),3);
ma = Inf(size(gTendonCoords.frame));
for i=1:size(gTendonCoords.frame,2)
    % Find best fit line through tendon coordinates
    points = gTendonCoords.coords{1,i};
    Q(i,:) = mean(points,1);
    [~, ~, V] = svd(bsxfun(@minus, points, Q(i,:)));
    tendonDir(i,:) = V(:,1)';
    % Unit vector along the tendon
    u = tendonDir(i,:)/norm(tendonDir(i,:));
    % Unit vector along functional joint axis
    v = fja(gTendonCoords.frame(1,i),:);
    % Vector from joint center (which is on the fja) to point on tendon
    P = jc(gTendonCoords.frame(1,i),:);
    PQ = Q(i,:) - P;
    % Moment arm
    ma(1,i) = abs(dot(cross(PQ,u),v));
end

% Convert from m to mm
ma = ma*1000;

% DATA CHECK: 3D plot of tendon lines of action, fja's, etc.
if (1)
    figure
    hold on
    % Proximal
    for i = 1:numel(proxIdx)
        prox = getMrkPos(mot_pos,mot_pos(proxIdx(i)).mrk_name);
        h1 = plotMrk(prox,'m.',pntrange);
    end
    % Distal
    for i = 1:numel(distIdx)
        dist = getMrkPos(mot_pos,mot_pos(distIdx(i)).mrk_name);
        h2 = plotMrk(dist,'g.',pntrange);
    end
    % Plot jc's
    h3 = plot3(jc(gTendonCoords.frame,1),jc(gTendonCoords.frame,2),jc(gTendonCoords.frame,3),'co');
    % Plot fja's
    f1=jc-0.2*fja; f2=jc+0.2*fja;
    for i=gTendonCoords.frame
        h4 = plot3([f1(i,1) f2(i,1)],[f1(i,2) f2(i,2)],[f1(i,3) f2(i,3)],'r-');
    end
    % Plot tendon centers
    h5 = plot3(Q(:,1),Q(:,2),Q(:,3),'yo');
    % Plot tendon lines of action
    t1=Q-0.1*tendonDir; t2=Q+0.1*tendonDir;
    for i=1:size(Q,1)
        h6 = plot3([t1(i,1) t2(i,1)],[t1(i,2) t2(i,2)],[t1(i,3) t2(i,3)],'b-');
    end
    hold off
    axis equal;
    legend([h1 h2 h3 h4 h5 h6],{'prox','dist','JC','FJA','TC','TLOA'})
end
%------------- END OF CODE --------------
end

%% HELPER FUNCTIONS
function [sValuesIdx, sValues] = valueSearch(values, varargin)
%valueSearch - constructs search list and finds nearest entry in 'values'
%Constructs a search list 'sValues' of all multiples of 'interval' (or 5)
%contained within 'values'. For each value in 'sValues', the entry in
%'values' with the closest value is found and its corresponding index is
%returned in 'sValuesIdx'.
%
% Syntax:  [sValuesIdx, sValues] = valueSearch(values, varargin)
%
% Inputs:
%    values (required) - double array (numValues x 1)
%           values to search
%    interval (optional) - double
%           spacing interval for search list
%
% Outputs:
%    sValuesIdx - double array (size of search list)
%           indices of nearest values found
%    sValues - double array (size of search list)
%           nearest values found
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
% February 2019; Last revision: 22-Feb-2019
%------------- BEGIN CODE --------------
if (nargin > 1)
    interval = varargin{1};
else
    interval = 5;
end
minAng = min(values);
maxAng = max(values);
sValues = ceil(minAng/5)*5:interval:floor(maxAng/interval)*interval;
sValuesIdx = zeros(size(sValues));

for i=1:size(sValues,2)
    test = abs(values - sValues(i));
    best = min(test);
    sValuesIdx(1,i) = find(test==best);
end
%------------- END OF CODE --------------
end