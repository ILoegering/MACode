function [posData, info] = organizeMocap4MA(input_data,varargin)
%organizeMocap4MA - helper function to organize mocap data for moment arm calculation
%Organizes mocap data struct output by load_trc_wVirtual.m or load_opto.m
%into two structures: a 'posData' struct array containing all position data
%with a separate struct for each marker and an 'info' struct containing all
%other information from the 'input_data' struct. This organization scheme
%allows markers to be more easily identified by a 'mrk_name' field and
%avoids indexing errors that can arise when extracting marker data using
%explicit indices from a general matrix containing data from all markers.
%Markers can also be renamed during restructuring if their assigned names
%don't follow the convention defined below. Marker names included in the
%'oldMrkNames' and 'newMrkNames' string cell arrays need not be in the same
%order as the markers in 'input_data,' but names should correspond between
%them (i.e., if oldMrkNames = {'foot1','shank2','foot3','LR',...}, then
%newMrkNames = {'dist1','prox2','dist3','w3',...}).
%
% Marker naming convention:
%   prox1, prox2, etc.,...      (proximal markers)
%	dist1, dist2, etc.,...      (distal markers)
%   trans1, trans2, etc.,...    (transducer markers)
%   w1, w2, w3, w4,...          (window markers) (order: LL, UL, LR, UR)
%	other, etc.,...             (other markers)
%
% Syntax:  [posData, info] = restructureMocap4MA(input_data, oldMrkNames, newMrkNames)
%
% Inputs:
%    input_data (required) - struct (1 x 1)
%           output of load_trc_wVirtual.m or load_opto.m
%    oldMrkNames (optional, paired) - string cell array (numMrks x 1)
%           marker names assigned during collection (order is arbitrary)
%    newMrkNames (optional, paired) - string cell array (numMrks x 1)
%           new marker names that follow convention (order must match
%           oldMrkNames)
%                       
%
% Outputs:
%    posData - struct array (numMrks x 1)
%           restructured position data
%    info - struct (contains same fields as input_data except for pos)
%           non-position data from input_data
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
if (nargin == 1)
    oldMrkNames = input_data.mrk_names;
    newMrkNames = input_data.mrk_names;
elseif (nargin == 3)
    % Check if all markers in input_data are included in oldMrkNames
    membership = ismember(input_data.mrk_names,varargin{1});
    if (membership)
        membership = ismember(varargin{1},input_data.mrk_names);
        % If there are marker names in oldMrkNames without corresponding 
        % data in input_data, prompt user
        if (sum(ismember(membership,0)))
            missing = varargin{1}(membership==0);
            answer = questdlg([missing{1} ' and possibly other markers were '...
                'included in oldMrkNames but don''t exist in '...
                'input_data. Would you like to exclude missing markers '...
                'and continue?'],'Missing Markers','Yes','No','Yes');
            if (~strcmp(answer,'Yes'))
                error(['There are marker names in oldMrkNames without '...
                    'corresponding position data in input_data.'])
            end
        end
        oldMrkNames = reshape(varargin{1}(membership==1),size(input_data.mrk_names));
        newMrkNames = reshape(varargin{2}(membership==1),size(input_data.mrk_names));
    else
        disp(input_data.mrk_names(membership==0))
        error('Markers listed above are missing from oldMrkNames.');
    end
else
    error('Must input either input_data OR input_data, oldMrkNames, and newMrkNames.');
end

% Save non-position data into 'info'
info = rmfield(input_data,{'pos','mrk_names'});
info.mrk_names = newMrkNames;
info = orderfields(info,[1:4,7,5,6]);

% Create variables using names from input_data containing position data for
% each marker
for i = 1:size(input_data.pos,3)
    varname = genvarname(input_data.mrk_names{i});
    eval([varname ' = input_data.pos(:,:,i);']);
end

% Create struct within posData for each marker and rename using names in
% newMrkNames
nmrk = numel(oldMrkNames);
for i = 1:nmrk
    posData(i).mrk_name = newMrkNames{i};
    eval(['posData(i).pos(:,:) = ' genvarname(oldMrkNames{i}) ';']);
end
%------------- END OF CODE --------------
end