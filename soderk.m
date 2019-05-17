function [Tb,aTb,bTa,res]=soderk(A,B)
% function [Tb,aTb,bTa,res]=soderk(A,B)
%
% Description:	Program calculates the transformation matrix T containing
%		the rotation matrix (3x3) and the translation translation 
%		vector d (3x1) for a rigid body segment using a singular
%		value decomposition method (Soederkvist & Wedin 1993).
% 
% Input:    A:   marker positions (XYZ are the columns) in an initial calibration frame
%           B:   marker positions (XYZ are the columns) in a current frame
%
% Output:   Tb:     4x4 Transformation Matrix relating ref frame B to ground reference frame
%                   containing the rotation matrix aRb and the translation from ground origin
%                   to origin at B, expressed in ground's reference frame, T = [R,d; 0 0 0 1]
%           aTb:    4x4 Transformation Matrix relating ref frame B to ref frame A
%                   containing the rotation matrix aRb and the translation from origin at
%                   A to origin at B, expressed in A's reference frame, T = [R,d; 0 0 0 1]
%           bTa:    4x4 Transformation Matrix relating ref frame A to ref frame B
%                   containing the rotation matrix bRa and the translation from origin at
%                   B to origin at A, expressed in B's reference frame, T = [R,d; 0 0 0 1]
%           res:    norm of residuals (measure of fit; "rigidity" of body
%
% References:     Soderkvist I. and Wedin P. -A., (1993). Determining the
%                 movements of the skeleton using well-configured markers.
%                 Journal of Biomechanics, 26:1473-1477      
%
% Author:	Christoph Reinschmidt, HPL, The University of Calgary
%               (Matlab code adapted from Ron Jacobs, 1993)
% Date:		February, 1995
% Last Changes: December 09, 1996
% Version:      3.1
%
% Modified: Darryl Thelen, Feb 11, 2010
% to accept different format of marker kinematic inputs and to output
% transformation matrices going both ways


if ((size(A,1)==1)&&(size(A,2)==3)&&(size(A,3)>=3))
    % marker positions are stored in a 3-d matrix
    nmrk=size(A,3);
    A=[reshape(A,3,nmrk)]';
elseif ((size(A,1)>=3)&&(size(A,2)==3)&&(size(A,3)==1))
    % markers are stored in a 2-D matrix with each row representing a marker
    % don't need to reshape, in preferred format
elseif ((size(A,1)==1)&&(size(A,2)>=9)&&(size(A,3)==1))
    A=reshape(A,3,size(A,2)/3);
else
    disp('ERROR: A matrix format is not recognized'); return
    return;
end

if ((size(B,1)==1)&&(size(B,2)==3)&&(size(B,3)>=3))
    % marker positions are stored in a 3-d matrix
    nmrk=size(B,3);
    B=[reshape(B,3,nmrk)]';
elseif ((size(B,1)>=3)&&(size(B,2)==3)&&(size(B,3)==1))
    % markers are stored in a 2-D matrix with each row representing a marker
    % don't need to reshape, in preferred format
elseif ((size(B,1)==1)&&(size(B,2)>=9)&&(size(B,3)==1))
    B=reshape(B,3,size(B,2)/3);
else
    disp('ERROR: B matrix format is not recognized'); return
    return;
end
    

% Checking for NaNs and also checking if still 3 pts left and if not
% T=[NaN...];
cut=[0];
qA=isnan(A); qB=isnan(B); qAB=[qA,qB];
qsum=sum(qAB'); cut=find(qsum~=0);
A([cut],:)=[]; B([cut],:)=[];
if size(A,1)<3,
    Tb=[NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN;];
    aTb=Tb;
    bTa=Tb;
    res=NaN;
    return;
end


Amean=mean(A)'; Bmean=mean(B)';


for i=1:size(A,1)-size(cut,2),
        Ai(:,i)=[A(i,:)-Amean']';
        Bi(:,i)=[B(i,:)-Bmean']';
end


C=Bi*Ai';
[P,T,Q]=svd(C);
R=P*diag([1 1 det(P*Q')])*Q';
d=Bmean-R*(Amean);  % Translation from A to B
bTa=[R,d;0 0 0 1];  % Transformation matrix from A reference frame to B reference frame
Tb=[R,Bmean;0 0 0 1];  % Transformation matrix from B reference frame to global reference frame

aTb=inv(bTa);

% Calculating the norm of residuals
A=A'; A(4,:)=ones(1,size(A,2));
B=B';
Bcalc=bTa*A; Bcalc(4,:)=[]; Diff=B-Bcalc; Diffsquare=Diff.^2;
%DOF=3*(number of points)-6 unknowns (Hx,Hy,Hz,alpha,beta,gamma):
DOF=size(B,1)*size(B,2)-6;
res=[sum(Diffsquare(:))/DOF].^0.5;