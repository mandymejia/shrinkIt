function [X_shrink lambda var_within var_between] = shrinkIt(X1_grp, X2_grp)
%
% This function performs shrinkage towards the group mean of subject-level 
% observations of any summary statistic computed from time series data.
% For example, consider an fMRI time series for each subject stored as a
% TxV array, and the VxV sample Fisher-transformed correlation matrix as 
% the summary statistic we wish to shrink.  
%
% Shrinkage computes a weighted average between subject-level estimates and
% the group mean.  Lambda represents the degree of shrinkage (weighting of 
% the group mean), ranging from 0 (subject-level estimates perfectly 
% reliable, so no shrinkage) to 1 (no reliable subject-level information, 
% so complete shrinkage to the group mean).  The optimal degree of shrinkage
% minimizes the MSE of the estimate relative to the truth, and is equal to the 
% ratio of within-subject variance to total variance, which is the sum of
% within-subject variance and between-subject variance.
%
% This function estimates the variance components to determine the optimal
% degree of shrinkage for each estimated parameter.  Therefore lambda has
% the same dimensions as each subject's array of estimated parameters.  
%
% The inputs X1 and X2 can be computed by applying the split_ts function to
% each subject's data.  X1 and X2 have the same dimensions, (q1, q2, ..., qk, n).  
% Each subject  1,...,n has an array of parameters we want to estimate, of 
% dimension (q1, q2, ..., qk).
%
%Usage:
%   [X_shrink lambda var_within var_between] = shrinkIt(X1, X2)
%Inputs:
%   X1 - An array of dimensions (q1, q2, ..., qk, n) containing parameter
%       estimates for each subject computed using the first half of the 
%       time series for each subject (see split_ts.m and Example.m)
%
%   X2 - An array of dimensions (q1, q2, ..., qk, n) containing parameter
%       estimates for each subject computed using the second half of the 
%       time series for each subject (see split_ts.m and Example.m)
%
%Outputs:
%   X_shrink - array of dimensions (q1, q2, ..., qk, n) containing the
%              shrinkage estimates of each parameter for each subject
%   lambda - array of dimensions (q1, q2, ..., qk, n) containing the degree
%            of shrinkage for each estimated parameter
%	var_within - within-subject variance estimates, dimension (q1, q2, ..., qk)
%   var_within - between-subject variance estimates, dimension (q1, q2, ..., qk)

%% Perform Checks

if(nargin ~= 2)
    error('Must specify two inputs')
end

if isempty(X1_grp) || isempty(X2_grp)
    error('one or more inputs is empty')
end

dims = size(X1_grp); %Returns the dimensions m by n of the observation matrix
if ~isequal(dims, size(X2_grp))
    error('dimensions of inputs do not match')
end      

if ~isnumeric(X1_grp) || ~isnumeric(X2_grp)
    error('all inputs must be numeric arrays')
end

if size(X1_grp, ndims(X1_grp)) == 1
    if(max(dims)==1)
        error('last dimension of inputs must equal number of subjects > 1')
    else
        %check that have the correct number of subjects
        qStr = sprintf('Is %s the number of subjects? Y/N [Y]:', num2str(max(dims)));
        userCheck = input(qStr);
        if(isempty(userCheck))
            %largest dimension indexes subjects
            [~, nd] = max(dims);
        else
            error('last dimension of inputs must equal number of subjects > 1')
        end %userCheck
    end %max(dims)==1
else
    %last dimension of arrays (indexes subjects)
    nd = ndims(X1_grp);
end

%% SET-UP

%compute array of estimates for each subject
X_grp = (X1_grp + X2_grp)/2; %subject-level estimates

%number of subjects
n = size(X_grp,nd);


%% COMPUTE WITHIN-SUBJECT VARIANCE USING X1 and X2

D = X2_grp - X1_grp; %compute intrasession differences
var_within = (1/4)*var(D, 0, nd); %consists of within-subject intrasession(signal) and sampling (noise) variance

%% COMPUTE TOTAL VARIANCE

var_total = var(X_grp, 0, nd);
var_between = var_total - var_within; %between-subject variance
var_between(var_between < 0) = 0;

%% COMPUTE LAMBDA (DEGREE OF SHRINKAGE)

lambda = var_within./var_total;
lambda(lambda > 1) = 1; %occurs when within-subject var > between-subject variance of estimates 
lambda(lambda < 0) = 0;
lambda(var_total==0) = 0; %for parameters with no variance, e.g. diagonal of correlation matrix 

%% PERFORM SHRINKAGE

%compute mean across subjects
X_bar = mean(X_grp, nd);

%make X_bar and lambda of same size as input arrays
X_bar = reshape(repmat(X_bar, 1, n), size(X1_grp));
lambda2 = reshape(repmat(lambda, 1, n), size(X1_grp));

%compute shrinkage estimates
X_shrink = (lambda2.*(X_bar))+((1-lambda2).*X_grp); 


