function [jc,fja,ang,JC,FJA] = functionalJointAxis(Pa,Da,Pb,Db)
%functionalJointAxis - Calculate functional joint axis
%At each frame, singular value decomposition is used to compute the
%homogeneous transformation between the proximal and distal segments from
%the positions of the cluster markers on each segment. A functional axis is
%then computed as the best-fit screw axis that described the motion of the
%distal segment with respect to the proximal segment.
%
%This code was adapted from Peter Adamczyk's Functional Joint Axis code
%and is based on the algorithm described in:
%'Evaluation of a new algorithm to determine the hip joint center'
%Siston and Delp, J Biomech 39:125-130, 2006
%
% Syntax:  [jc,fja,ang,JC,FJA] = functionalJointAxis(Pa,Da,Pb,Db)
%
% Inputs:
%    Pa (required) - double matrix (numFrames x 3 x numProxMarkers)
%           Calibration (static) data from proximal marker shell
%    Da (required) - double matrix (numFrames x 3 x numDistMarkers)
%           Calibration (static) data from distal marker shell
%    Pb (required) - double matrix (numFrames x 3 x numProxMarkers)
%           Motion (dynamic) data from proximal marker shell
%    Db (required) - double matrix (numFrames x 3 x numDistMarkers)
%           Motion (dynamic) data from distal marker shell
%
% Outputs:
%    jc - varType (size)
%           joint center in global coordinates at each frame
%    fja - varType (size)
%           functional joint axis in global coordinates at each frame
%    ang - double array (numFrames x 3)
%           computed body 312 angles, column 1 corresponds to rotation
%           about functional joint axis
%    JC - double array (6 x 1)
%           joint center in proximal segment, expressed in global coordinates
%    FJA - double array (1 x 3)
%           function joint axis relative to proximal segment
%
%   ...all marker kinematics assumed expressed in global coordinates
%
%
% Other m-files required: body312ang.m, insertNaN.m, removeNaN.m, soderk.m
% Subfunctions: nulltol
% MAT-files required: none
%
% Author: Peter Adamczyk, Darryl Thelen, and Isaac Loegering
% BADGER Lab and UW Neuromuscular Biomechanics Lab
% University of Wisconsin-Madison
% 1513 University Ave, Rm 3046
% Madison, WI 53706
% email: peter.adamczyk@wisc.edu, dgthelen@wisc.edu, isaacloegering@gmail.com
% March 2018; Last revision: 17-May-2019
%------------- BEGIN CODE --------------
disp('Computing functional joint center');

%-----------------------------REMOVE DROPOUT------------------------------%
% Save original size of motion data structure
originalSize = size(Pb,1);
% Remove frames with dropout from calibration data
Pa = removeNaN(Pa);
Da = removeNaN(Da);
% Remove frames where dropout occurred in either segment for motion data
PDb = cat(3,Pb,Db);
[~,idxNaN_mot] = removeNaN(PDb);
idx_mot = setxor(idxNaN_mot,1:originalSize);
Pb = Pb(idx_mot,:,:);
Db = Db(idx_mot,:,:);
%-----------------------------REMOVE DROPOUT------------------------------%

% Compute average marker positions in the calibration trials
Pa_mrk=(squeeze(mean(Pa,1)))';
Da_mrk=(squeeze(mean(Da,1)))';
j=1;
% Assemble matrices needed to compute joint center, see Delp and Siston
for i=1:size(Pb,1)
    % Extract a single frame of marker data from motion file
    Pb_mrk=(squeeze(Pb(i,:,:)))';
    Db_mrk=(squeeze(Db(i,:,:)))';
    % Find the transformation matrices for both segments relative to
    % calibration reference frames
    [Tpb,paTpb,pbTpa]=soderk(Pa_mrk,Pb_mrk);
    [Tdb,daTdb,dbTda]=soderk(Da_mrk,Db_mrk);
    dbTpb=inv(Tdb)*Tpb;
    % Assemble the matrices to compute joint center in least squares sense
    LHS(j:j+2,1:3)=Tpb(1:3,1:3);
    LHS(j:j+2,4:6)=-Tdb(1:3,1:3);          % formerly -eye(3)
    RHS(j:j+2,1)=Tdb(1:3,4)-Tpb(1:3,4);
    j=j+3;
end
% compute the joint center, which will be expressed in proximal ref frame
% ...the first 3 XYZ values are joint center relative to proximal segment
% ...the next 3 XYZ values are joint center relative to distal segment
XYZ=LHS\RHS;
JC=XYZ;
% check the Singular Values to see how bad this FJC computation might
% be (null space describes a line (1 vector) or plane (2 vectors) that
% are indeterminate)
nullAxis = NaN*[0 0 0];
[~,s,~] = svd(LHS,'econ'); s = diag(s)';
n=1;
tol = mean(s(end-n:end-n+1))-eps;   % put it between the relevant vectors
nullSpace = nulltol(LHS,tol); % get the 1-axis null space of the computation
nullAxis = nullSpace(1:3,:);
nullAxis = nullAxis/sqrt(sum(nullAxis(:).^2));
% NOTE this is only valid for a single NullAxis. If the Null Space has two or more bases, the computations need to be redone.
% Now have the functional joint axis expressed in the proximal ref frame
FJA = nullAxis';
% Define unit vectors normal to FJA to create a functional reference frame
uz=FJA';
uy=mean(Pa_mrk)'; %-XYZ(1:3);
ux=cross(uy,uz);    ux=ux/sqrt(sum(ux.^2));
uy=cross(uz,ux);    uy=uy/sqrt(sum(uy.^2));
Rfja0=[ux uy uz];
% Find the joint center and functional joint axes at each frame
for i=1:size(Pb,1)
    Pb_mrk=(squeeze(Pb(i,:,:)))';
    % Transform markers to the fja ref frame
    [fjaTpb,paTpb,pbTpa]=soderk(Pa_mrk*Rfja0,Pb_mrk*Rfja0);
    Db_mrk=(squeeze(Db(i,:,:)))';
    [fjaTdb,daTdb,dbTDa]=soderk(Da_mrk*Rfja0,Db_mrk*Rfja0);
    % Convert the joint center location to the ground reference frame
    [Tpb,paTpb,pbTpa]=soderk(Pa_mrk,Pb_mrk);
    loc=Tpb*[XYZ(1:3); 1]; % this was XYZ(1:3)
    jc(i,:)=loc(1:3)';
    fja(i,:)=(Tpb(1:3,1:3)*FJA')';
    % Determine the rotation angle
    pbTdb=inv(fjaTpb)*fjaTdb;
    R=pbTdb(1:3,1:3);
    q=body312ang(R);
    ang(i,:)=q;
end
%----------------------------REINSERT DROPOUT-----------------------------%
jc = insertNaN(jc,idxNaN_mot);
fja = insertNaN(fja,idxNaN_mot);
ang = insertNaN(ang,idxNaN_mot);
%----------------------------REINSERT DROPOUT-----------------------------%
end


function Z = nulltol(A,arg2,arg3)
%NULL   Null space.
%   Z = NULL(A) is an orthonormal basis for the null space of A obtained
%   from the singular value decomposition.  That is,  A*Z has negligible
%   elements, size(Z,2) is the nullity of A, and Z'*Z = I.
%
%   Z = NULL(A,'r') is a "rational" basis for the null space obtained
%   from the reduced row echelon form.  A*Z is zero, size(Z,2) is an
%   estimate for the nullity of A, and, if A is a small matrix with
%   integer elements, the elements of R are ratios of small integers.
%
%   The orthonormal basis is preferable numerically, while the rational
%   basis may be preferable pedagogically.
%
%   Example:
%
%       A =
%
%           1     2     3
%           1     2     3
%           1     2     3
%
%       Z = null(A);
%
%       Computing the 1-norm of the matrix A*Z will be
%       within a small tolerance
%
%       norm(A*Z,1)< 1e-12
%       ans =
%
%          1
%
%       null(A,'r') =
%
%          -2    -3
%           1     0
%           0     1
%
%   Class support for input A:
%      float: double, single
%
%   See also SVD, ORTH, RANK, RREF.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.12.4.5 $  $Date: 2010/09/02 13:36:17 $
%
%  CUSTOMIZED BY PGA TO ALLOW A TOLERANCE ARGUMENT!!!

[m,n] = size(A);
if (nargin > 1) && isequal(arg2,'r')
    
    % Rational basis
    
    %    [R,pivcol] = rref(A); % original.
    
    % added by PGA to use a Tolerance:
    if nargin > 2
        tol = arg3;
        [R,pivcol] = rref(A,tol);
    else
        [R,pivcol] = rref(A);
    end
    
    r = length(pivcol);
    nopiv = 1:n;
    nopiv(pivcol) = [];
    Z = zeros(n,n-r,class(A));
    if n > r
        Z(nopiv,:) = eye(n-r,n-r,class(A));
        if r > 0
            Z(pivcol,:) = -R(1:r,nopiv);
        end
    end
else
    % Orthonormal basis
    [~,S,V] = svd(A,0);
    if m > 1, s = diag(S);
    elseif m == 1, s = S(1);
    else s = 0;
    end
    
    if (nargin > 1)% added by PGA to put in a special Tolerance value, similar to GNU octave.
        tol = arg2;
    else
        tol = max(m,n) * max(s) * eps(class(A));
    end
    r = sum(s > tol);
    Z = V(:,r+1:n);
end
end
