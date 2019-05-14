function [pos] = getMrkPos(posData,mrk_name)
%getMrkPos - extract marker position data from posData
%Returns the (x,y,z)-coordinates of the marker named 'mrk_name' in the
%struct array 'posData'
%
% Syntax:  [pos] = getMrkPos(posData,mrk_name)
%
% Inputs:
%    posData (required) - struct array (numMrks x 1)
%           output of organizeMocap4MA, containing a struct for each marker
%           with fields 'pos' and 'mrk_name'
%    mrk_name (required) - string
%           name of marker for which data is to be extracted
%
%
% Outputs:
%    pos - double matrix (numFrames x 3)
%           position data for marker named 'mrk_name'
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
pos = posData(strcmp({posData.mrk_name},mrk_name)).pos;
%------------- END OF CODE --------------
end