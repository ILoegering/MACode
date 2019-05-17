function [f] = getMAFit(data_path,fits_path,subjectID,cond,kinemID,trials,db_path,db_name)
%getMAFit - fit a second degree polynomial to the data from select trials
%Moment arm data from a subject's good trials specified by 'trials' are
%loaded, and a second-degree polynomial is fit to the data. These data are
%then stored in database 'db'.
%
% This function assumes data files are named in the following manner:
%'<subjectID>_<kinemID>_<cond><trial>.mat'. Also, the name of the database
%variable contained within the database file should match the filename
%(i.e., db_name.mat should contain a struct array called db_name).
%
% Syntax:  [f] = getMAFit(data_path,plot_path,subjectID,kinemID,trials,db_path,db_name)
%
% Inputs:
%    data_path (required) - char array (varies)
%           Path of folder in which data files can be found
%    fits_path (required) - char array (varies)
%           Path of folder in which plot should be saved
%    subjectID (required) - char array (varies)
%           Unique subject ID
%    cond (required) - char array (1 x 1)
%           Condition ('A' for active, 'P' for passive)
%    kinemID (required) - char array (varies)
%           Kinematic identifier
%    trials (required) - char array (varies)
%           Trials selected for fitting
%    db_path (required) - char array (varies)
%           Path of folder in which database file is located
%    db_name (required) - char array (varies)
%           Database name (should be db filename without '.mat')
%
% Outputs:
%    f - Figure handle (1 x 1)
%           Handle of figure created with fit
%
%
% Other m-files required: Curve Fitting Toolbox
% Subfunctions: none
% MAT-files required: none
%
% Author: Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: isaacloegering@gmail.com
% July 2018; Last revision: 17-May-2019
%------------- BEGIN CODE --------------
%% Load data
% Angle range
minAng = Inf;
maxAng = -Inf;
% For each trial...
for i = 1:numel(trials)
    % Load data
    data(i) = load([data_path '\' subjectID '_' kinemID '_' cond num2str(trials(i)) '.mat']);
    % Fit data (regular and robust)
    p(:,i) = polyfit(data(i).angle,data(i).ma,2);
    r{1,i} = fit(data(i).angle',data(i).ma','poly2','Robust','Bisquare');
    % Determine minimum and maximum angles for which the moment arm was calculated for all trials selected
    if (min(data(i).angle) < minAng)
        minAng = min(data(i).angle);
    end
    if (max(data(i).angle) > maxAng)
        maxAng = max(data(i).angle);
    end
end

% Calculate average polyfit
p_net = sum(p,2)/size(trials,2)

% Calculate average robust fit
sum_r1 = 0;
sum_r2 = 0;
sum_r3 = 0;
for i = 1:size(r,2)
    sum_r1 = sum_r1 + r{1,i}.p1;
    sum_r2 = sum_r2 + r{1,i}.p2;
    sum_r3 = sum_r3 + r{1,i}.p3;
end
r_net(1,1) = sum_r1/size(r,2);
r_net(2,1) = sum_r2/size(r,2);
r_net(3,1) = sum_r3/size(r,2)

% Save fits to database
load([db_path '\' db_name '.mat'],db_name)
eval(['if (~sum(strcmp({' db_name '.id}, subjectID))) '...
    db_name '(numel(' db_name ')+1).id = subjectID; '...
    'end'])
fieldMAType = [cond '_' kinemID '_MA'];  % Used for naming field in database
try
    eval([db_name '(strcmp({' db_name '.id},subjectID)).' fieldMAType '_poly=p_net'';']);
    eval([db_name '(strcmp({' db_name '.id},subjectID)).' fieldMAType '_robust=r_net'';']);
    eval([db_name '(strcmp({' db_name '.id},subjectID)).' fieldMAType '_valAt0=p_net(3)'';'])
catch
    disp('This subject is not currently in the database.')
end
save([db_path '\' db_name '.mat'],db_name)

% Plot data
colors = ['r' 'g' 'm' 'c','b'];
f = figure;
hold on
for i = 1:size(data,2)
    plot(data(i).angle,data(i).ma,[colors(i) '*'])
end

% Plot fits
x = minAng:0.1:maxAng;
for i = 1:size(p,2)
    y = p(1,i)*x.^2 + p(2,i)*x + p(3,i);
    h(i) = plot(x,y,[colors(i) ':']);
end
yp_net = p_net(1,1)*x.^2 + p_net(2,1)*x + p_net(3,1);
yr_net = r_net(1,1)*x.^2 + r_net(2,1)*x + r_net(3,1);
hp_net = plot(x,yp_net,'k');
hr_net = plot(x,yr_net,'--b');
hold off
title('MA vs. Joint Angle')
xlabel('Joint Angle (deg)')
ylabel('Tendon MA (mm)')
switch size(trials,2)
    case 1
        legend([hp_net hr_net],{[cond num2str(trials(1)) ' - polyfit'],[cond num2str(trials(1)) ' - robust fit']},'Location','Southeast')
    case 2
        legend([h(1),h(2),hp_net,hr_net],{[cond num2str(trials(1))],[cond num2str(trials(2))],'polyfit avg','robust fit avg'},'Location','Southeast')
    case 3
        legend([h(1),h(2),h(3),hp_net, hr_net],{[cond num2str(trials(1))],[cond num2str(trials(2))],[cond num2str(trials(3))],'polyfit avg','robust fit avg'},'Location','Southeast')
    case 4
        legend([h(1),h(2),h(3),h(4),hp_net, hr_net],{[cond num2str(trials(1))],[cond num2str(trials(2))],[cond num2str(trials(3))],[cond num2str(trials(4))],'polyfit avg','robust fit avg'},'Location','Southeast')
    case 5
        legend([h(1),h(2),h(3),h(4),h(5),hp_net, hr_net],{[cond num2str(trials(1))],[cond num2str(trials(2))],[cond num2str(trials(3))],[cond num2str(trials(4))],[cond num2str(trials(5))],'polyfit avg','robust fit avg'},'Location','Southeast')
end
% Save plots
saveas(f,[fits_path  '\fig\' subjectID '_' kinemID '_' cond '.fig'])
saveas(f,[fits_path  '\fig\' subjectID '_' kinemID '_' cond '.tiff'])
end
%------------- END OF CODE --------------