function angles = body312ang(R)
%body312ang - % Solves for successive Body 3-1-2 rotation angles q1, q2, q3
%that would generate the rotation matrix R
%
% Syntax:  [angles] = body312ang(R)
%
% Inputs:
%    R (required) - double matrix (3 x 3)
%           Rotation matrix
%
% Outputs:
%    angles - double array (1 x 3)
%           Rotation angles q1, q2, and q3 given in the range from -pi to
%           pi in radians
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Darryl Thelen
% Editor: Isaac Loegering
% UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: dgthelen@wisc.edu
% A long time ago in a biomechancs lab far, far away; Last revision: 14-May-2019
%------------- BEGIN CODE --------------
% Check size of rotation matrix (must be 3 x 3)
if size(R)~=[3,3]
    disp('Error: rotation matrix has to be a 3x3 matrix')
    return;
end
% If any element is NaN, return array of NaN
if sum(isnan(R(:)))~=0, angles=[NaN,NaN,NaN]; return; end
% Calculate q2
q2 = asin(R(3,2));  %'assumption' that cos(alpha)>0
% Calculate q3
q3sin = asin(-R(3,1)/cos(q2));
q3cos = acos(R(3,3)/cos(q2));
if (q3cos>pi/2 && q3sin>0); q3=pi-q3sin; end
if (q3cos>pi/2 && q3sin<0); q3=-pi-q3sin; end
if (q3cos<=pi/2); q3=q3sin; end
% Calculate q1
q1sin= asin(-R(1,2)/cos(q2));
q1cos= acos(R(2,2)/cos(q2));
if (q1cos>pi/2 && q1sin>0); q1=pi-q1sin;  end
if (q1cos>pi/2 && q1sin<0); q1=-pi-q1sin; end
if (q1cos<=pi/2);  q1=q1sin;     end
% Combine angles into output array
angles=[q1 q2 q3];
end