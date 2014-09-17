%% Transform a vector of upper triangular elements to a V-by-V matrix 
function [ MAT ] = UT2mat(UT, d)

% This function takes a vector of V(V-1)/2 elements and creates a matrix
% that contains those elements on the upper triangle, their transpose on 
% the lower triangle, and a constant diagonal.

% UT - vector of V(V-1)/2 elements to populate the upper and lower
% triangles of MAT
% d - scalar value to populate the diagonal of MAT.  defaults to 1.

%% Perform Checks and Initialize diag, M and V

%Check that UT is a numeric vector
if ~isnumeric(UT) || ~isvector(UT)
    error('UT must be a numeric vector')
end

%If d is not provided, set to 1
%If d is provided, check that it is a numeric scalar
if nargin == 1
    d = 1;
elseif ~isnumeric(d) || length(d)>1
    error('diag must be a numeric scalar')
end

%Compute M and V.  Check that V is integer.
M = max(size(UT));
V = (1+sqrt(8*M+1))/2;
if round(V)~=V
    error('Length of UT implies non-integer number of voxels.  Length should be equal to V(V-1)/2, where V is integer.')
end

%% Create MAT from UT and diag

%Create upper triangle
MAT_UT = triu(ones(V,V),1);
MAT_UT(MAT_UT==1) = UT;

%Create lower triangle
MAT_LT = MAT_UT';

%Create diagonal
MAT_diag = diag(ones(V,1));
MAT_diag(MAT_diag==1) = d;

%Add together to form MAT
MAT = MAT_UT + MAT_LT + MAT_diag;


end

