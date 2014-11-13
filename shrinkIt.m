function [shrink lambda] = shrinkIt(X1, X2, t)
%
% This function performs shrinkage towards the group mean of subject-level 
% observations of any voxel-by-voxel similarity or distance matrix.  This 
% can be used as a pre-processing step before performing clustering in 
% order to improve reliability of subject-level parcellations.  If there 
% are M variables (e.g. the elements in the upper triangle of the voxel-
% by-voxel similarity/distance matrix) observed for each subject, then a 
% different shrinkage parameter will be computed for each variables. The 
% shrinkage parameter (lambda) ranges from 0 (no shrinkage) to 1 (complete 
% shrinkage), and is determined by the relationship between within-subject 
% (noise) variance and between-subject (signal) variance.  Both variance 
% components and the shrinkage parameter are estimated from the data.
%
%Usage:
%   [shrink lambda] = shrinkIt(X1, X2, t)
%Inputs:
%   X1 - An m-by-n matrix, where the ith column contains the m observed 
%        values of subject i.  Half of the total scan time available should 
%        be used to compute these values.  (see below for more details)
%
%   X2 - An m-by-n matrix, where the ith column contains a second set of 
%        the m observed values of subject i.  The remaining half of the
%        total scan time available should be used to compute these values
%        (see below for more details)
% 
%   t -  A scalar equal to the total amount of scan time (in minutes) 
%        collected for each subject.  This value will be used to adjust the 
%        noise variance estimate, since it will tend to be over-estimated 
%        by splitting the data to compute X1 and X2.
%
%Outputs:
%   shrink - m-by-n matrix of shrunken observations
%   lambda - vector (length m) of shrinkage parameter values between 0 (no
%   shrinkage) and 1 (complete shrinkage)

%% Perform Checks

%Check that:
% - X1 and X2 are provided, are numeric, and have the same dimensions
% - t is provided and is a numeric scalar
if(nargin ~= 3)
    error('Must specify three inputs: the 2 observations from each subject and the total scan time (in minutes) collected for each subject')
end
dims = size(X1); %Returns the dimensions m by n of the observation matrix
if isempty(X1) || isempty(X2)
    error('Both X1 and X2 should be provided in order to compute noise variance')
end
if ~isequal(dims, size(X2))
    error('Dimensions of X1 and X2 do not match')
end      
if ~isnumeric(X1) || ~isnumeric(X2)
    error('X1 and X2 must be numeric')
end
if length(t)~=1 || ~isnumeric(t)
    error('t must be a numeric scalar')
end

%% Initialize m, n and X
    
m = dims(:,1);  %number of variables
n = dims(:,2);  %number of subjects 

%Compute Raw Estimate
X = (X1 + X2)/2;

%% Compute Variance Components
   
%Compute Noise Variance (global scalar value):
Varu = mean(0.5*var(X2 - X1, 0, 2));

%Adjust Noise Variance (to account for using scans of length t/2)
theta = min(0.590 + 0.129*log(t), 1); 
Varu = Varu*theta;
 
%Compute Total Variance
Varw = var(X,0,2);  %for rows. vector, length m

%Compute Signal Variance (between subjects):
Varx = Varw-Varu; %vector of length m
Varx(Varx<0) = 0;  %set values less than 0 to 0

%% Compute Shrinkage Parameter & Perform Shrinkage
%The higher the noise variance, the more shrinkage towards the group mean 

lambda = (Varu)./(Varx+Varu); %vector length #observations
lambda = repmat(lambda, 1, n); %use same lambda vector for all subjects

%% Perform Shrinkage
groupmean = mean(X,2);  %vector of length m
groupmean = repmat(groupmean, 1, n); %matrix
shrink = (lambda.*(groupmean))+((1-lambda).*X); %mxn matrix

%% Make lambda into a vector again
lambda = lambda(:,1);

