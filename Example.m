
addpath /Users/mejiaa/Documents/Software/shrinkIt/

%% Randomly generate connectivity matrices MAT1 and MAT2 for 20 subjects

% - Suppose there are 100 voxels in the ROI
% - Size of connectivity matrix = 100 x 100
% - Elements in upper triangle = 100(100-1)/2 = 4950

V = 100;
m = 4950;
n = 20;

MAT1 = zeros(V,V,n);
MAT2 = zeros(V,V,n);
for i = 1:20
    
    %Randomly generate values in upper triangle
    %Generate scalar mean for subject i with signal sd=0.1
    %Then generate m samples around mean mu_i with noise sd=0.1
    mu_i = normrnd(0.2,0.1);
    vals1 = normrnd(mu_i,0.1,m,1);
    vals2 = normrnd(mu_i,0.1,m,1);
    
    %Put values into 100-by-100 matrix
    MAT1(:,:,i) = UT2mat(vals1);
    MAT2(:,:,i) = UT2mat(vals2);

end

%% RESHAPE ARRAYS MAT1 AND MAT2 INTO m-by-n MATRICES

% Extract the upper triangle from each 100-by-100 matrix.  Gives a vector
% of length m.  Store these in X1 and X2.

X1 = zeros(m,n);
X2 = zeros(m,n);
for i = 1:20
    
    X1(:,i) = mat2UT(MAT1(:,:,i));
    X2(:,i) = mat2UT(MAT2(:,:,i));

end

%% APPLY INVERSE FISHER TRANSFORMATION TO OBTAIN CORRELATION COEFFICIENTS

X1 = (exp(2*X1)-1)./(exp(2*X1)+1);
X2 = (exp(2*X2)-1)./(exp(2*X2)+1);

%% RUN shrinkIt TO OBTAIN SHRINKAGE ESTIMATE X (m-by-n MATRIX) AND SHRINKAGE PARAMETER

%Assume each subject has 10 minutes of scan time 
%(5 minutes for X1 and 5 minutes for X2)

[X lambda] = shrinkIt(X1, X2, 10);

%Look at distribution of shrinkage parameter values
%Since signal and noise were both generated with sd=0.1, 
%and lambda = var(noise)/[var(noise)+var(signal)],
%lambda should be centered around 0.5.
hist(lambda)

