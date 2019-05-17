function [s] = nanFill(s)
%nanFill - Replace empty values in struct array with NaN
%
% Syntax:  [s] = nanFill(s)
%
% Inputs:
%    s (required) - struct array (1 x numStructs)
%           Structure array to fill
%
% Outputs:
%    s - struct array (1 x numStructs)
%           Filled structure array
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
% May 2019; Last revision: 14-May-2019
%------------- BEGIN CODE --------------
% Get field names
fn = fieldnames(s);
% For each field name...
for i = 1:numel(fn)
    % Get indices of empty values
    idx = find(cellfun(@isempty,{s.(fn{i})})==1);
    % Fill with NaN
    if (~isempty(idx))
        for j = 1:numel(idx)
            s(idx(j)).(fn{i}) = nan;
        end
    end
end
end
%------------- END OF CODE --------------