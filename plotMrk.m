function [h] = plotMrk(mrk_pos,varargin)
%plotMrk - plot 3D marker data
%Uses plot3 to plot a marker's position data
%
% Syntax:  [h] = plotMrk(mrk_pos,lineSpec,frames)
%
% Inputs:
%    mrk_pos (required) - struct (numFrames x 3)
%           marker position data
%    lineSpec (optional) - string
%           lineSpec for plotted marker data; default is black dots, 'k.'
%    frames (optional) - double array (1 x numFrames)
%           indices of frames to plot; default if all frames
%
% Outputs:
%    h - handle
%           plot handle
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
if (nargin==1)
    lineSpec = 'k.';
    frames = 1:size(mrk_pos,1);
elseif (nargin==2)
    lineSpec = varargin{1};
    frames = 1:size(mrk_pos,1);
else
    lineSpec = varargin{1};
    frames = varargin{2};
end
h = plot3(mrk_pos(frames,1),mrk_pos(frames,2),mrk_pos(frames,3),lineSpec);
%------------- END OF CODE --------------
end