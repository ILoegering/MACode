function [output_data] = selectFrames(input_data,idx)
%selectFrames - extract frames with indices given by 'idx' from
%'input_data'
%Extracts frames with indices given by 'idx' from the struct array of mocap
%data given by 'input_data', which is the output of organizeMocap4MA.
%
% Syntax:  [output_data] = selectFrames(input_data,idx)
%
% Inputs:
%    input_data (required) - struct array (numMrks x 1)
%           marker data in format output by organizeMocap4MA.m
%    idx (required) - double array (1 x numFrames)
%           indices for frames to be extracted from 'input_data'
%
%
% Outputs:
%    output_data - struct array (numMrks x 1)
%           position data with extracted frames
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
% February 2019; Last revision: 16-Feb-2019
%------------- BEGIN CODE --------------
output_data = struct('mrk_name', cell(1, numel(input_data)), 'pos', cell(1, numel(input_data)));
% For each marker...
for i = 1:numel(input_data)
    % DATA CHECK: Plot data before selecting frames
    if(0)
        figure
        hold on
        plot(input_data(i).pos(:,1),'b.')    % x-coord
    end
    sample = zeros(numel(idx),size(input_data(i).pos,2));
    % For each frame to be selected
    for j = 1:numel(idx)
        % Extract selected frame using indices in idx
        sample(j,:) = input_data(i).pos(idx(j),:);
    end
    % Store result
    output_data(i).mrk_name = input_data(i).mrk_name;
    output_data(i).pos = sample;
    % DATA CHECK: Plot data after selecting frames
    if(0)
        plot(idx,sample(:,1),'r.')
        hold off
    end
end
%------------- END OF CODE --------------
end