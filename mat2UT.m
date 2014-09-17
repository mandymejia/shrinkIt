%% Transform a vector of upper triangular elements to a V-by-V matrix 
function [ UT ] = mat2UT(MAT)

% This function takes a V-by-V matrix and creates a vector of the V(V-1)/2 
% elements on the upper triangle.

% INPUTS
%
% MAT - V-by-V matrix

% RETURNS
%
% UT - vector of upper triangular elements of mat

%% Perform Checks and Initialize diag, M and V

%Check that a single argument was provided
if nargin ~= 1
    error('Function takes exactly one argument')
end


%Check that mat is a numeric matrix
if ~isnumeric(MAT) || ~ismatrix(MAT)
    error('MAT must be a numeric matrix')
end

%Check that mat is square
if size(MAT,1) ~= size(MAT,2)
    error('MAT must be square')
end

%Compute M and V.  Check that V is integer.
V = size(MAT,1);
M = V*(V-1)/2;

%% Create UT from mat

I_UT = triu(ones(V,V),1);
UT = MAT(I_UT==1);

end

