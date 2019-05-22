function [f] = createCmbdMAvAngPlot(data_path,plot_path,subjectID,kinemID)
%createCmbdMAvAngPlot - create plot of all moment arm vs. joint angle data for a single subject
%This function assumes the filenames of all data files contain the
%following regular expression: 'subjectID_kinemID'.
%
% Syntax:  [f] = createCmbdMAvAngPlot(data_path,plot_path,subjectID,kinemID)
%
% Inputs:
%    data_path (required) - char array (varies)
%           Path of folder in which data files can be found
%    plot_path (required) - char array (varies)
%           Path of folder in which plot should be saved
%    subjectID (required) - char array (varies)
%           Unique subject ID
%    kinemID (required) - char array (varies)
%           Kinematic identifier
%
% Outputs:
%    f - Figure handle (1 x 1)
%           Handle of figure created
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
% July 2018; Last revision: 21-May-2019
%------------- BEGIN CODE --------------
% Get all files in the data folder for this subject
extension = '.mat';
files = dir([data_path '\*' subjectID '_' kinemID '*' extension]);

% Plot data
f = figure;
lineSpec = {'k*','b*','r*','g*','m*'};
hold on
% Loop through data files
for id = 1:length(files)
    % Get filename of data file
    filename = files(id).name;
    % Load data from data file
    data = load([data_path '\' filename]);
    % Plot moment arm vs. joint angle data
    plot(data.angle,data.ma,lineSpec{id})
    clearvars -except f files data_path plot_path subjectID kinemID lineSpec
end
hold off
title('MA vs. Joint Angle')
xlabel('Joint Angle (deg)')
ylabel('Tendon MA (mm)')
% Create cell array containing names of plots
for i = 1:numel(files)
    plotNames{i} = ['Trial ' num2str(i)];
end
legend(plotNames)
% Save figure in tiff and fig formats
saveas(f,[plot_path '\Moment Arm\tiff\' subjectID '_' kinemID '_combined'],'tiff')
saveas(f,[plot_path '\Moment Arm\fig\' subjectID '_' kinemID '_combined'],'fig')
end
%------------- END OF CODE --------------