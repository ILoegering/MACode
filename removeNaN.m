function [B,idxNaN] = removeNaN(A)
%removeNaN - Remove rows (1st dim) from matrix containing NaN
%
% Syntax:  [B,idx] = removeNaN(A)
%
% Inputs:
%    A (required) - double matrix (numRows x numCol x ... x ...)
%           Matrix from which to remove rows that contain NaN entries in
%           any dimension
%
% Outputs:
%    B - double matrix (numRowsWithoutNaN x numCol x ... x ...)
%           Original matrix minus those rows that contain NaN entries in
%           any dimension
%    idxNaN - double array (numRowsWithoutNaN x 1)
%           Indices of rows that contain NaN entries in any dimension
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
% July 2018; Last revision: 14-May-2019
%------------- BEGIN CODE --------------
% Identify NaN entries
temp = isnan(A);
% Collapse matrix by summing elements along its last dimension repeatedly
% until it is a column vector
while(~isvector(temp))
    temp = squeeze(sum(temp,ndims(temp)));
end
% Determine indices of rows with and without NaN
idxNaN = find(temp~=0);
idx = setxor(idxNaN,1:size(A,1));
% Extract rows with no NaN from A
S.subs = repmat({':'},1,ndims(A));
S.subs{1} = idxNaN;
S.type = '()';
B = subsasgn(A,S,[]);
end