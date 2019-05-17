function [B] = insertNaN(A,idxNaN)
%insertNaN - insert NaN into matrix at specified indices
%Insert rows with all elements equal to NaN into matrix A at indices given
%by idxNaN.
%
% Syntax:  [B] = insertNaN(A,idxNaN)
%
% Inputs:
%    A (required) - double matrix (numRowsWithoutNaN x numCols x ... x ...)
%           Matrix into which rows of NaN are to be inserted
%    idxNaN (required) - double array (numIdxNaN x 1)
%           Indices at which rows of NaN are to be inserted
%
% Outputs:
%    B - double matrix (numRows x numCols x ... x ...)
%           Original matrix with the addition of the inserted rows of NaN
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
% July 2018; Last revision: 17-May-2019
%------------- BEGIN CODE --------------
% Copy A to B
B = A;
% Create NaN insert row
sz = size(B);
insert = nan([1,sz(2:end)]);
% Insert NaN row at each index in idxNaN
for i = 1:numel(idxNaN)
    % If index is 1, insert NaN row at beginning
    if (idxNaN(i) == 1)
        B = [insert; B];
    % Otherwise, if index is greater than the current number of rows in A,
    % insert row of NaN at end
    elseif (idxNaN(i) > size(B,1))
        B = [B; insert];
    % Otherwise, insert row of NaN at given index, which is between rows in
    % A
    else
        S.type = '()';
        S.subs = repmat({':'},1,ndims(B));
        eval(['S.subs{1} = 1:' num2str(idxNaN(i)-1) ';'])
        pt1 = subsref(B,S);
        eval(['S.subs{1} = ' num2str(idxNaN(i)) ':' num2str(size(B,1)) ';'])
        pt2 = subsref(B,S);
        B = [pt1; insert; pt2];
    end
end
end