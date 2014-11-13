%In this example, we first simulate scan and rescan connectivity matrices 
%for 20 subjects, using random normals with signal (across-subject) 
%variance of 0.1^2 and noise (within-subject) variance of 0.1^2 (on the 
%Normal scale).  We apply the inverse Fisher transform to obtain 
%correlations.
%
%We then use the function mat2UT to extract the upper-triangle of each
%matrix.  This is done because the shrinkIt function expects a matrix of
%size m-by-n, where m is the number of variables observed for each subject
%and n is the number of subjects. 
%
%We use the shrinkIt function to apply shrinkage.  This function requires
%two measures for each subject (scan and rescan), as well as the total scan
%time (across scan and rescan) per subject.  It then estimates the optimal 
%degree of shrinkage, which depends on the estimated signal and noise
%variance.  The estimated noise variance depends on the total scan time per
%subject. 

%% ADD PATH TO SHRINKIT TOOLBOX

%% SIMULATE CONNECTIVITY DATA 

%In a real application, the rs-fMRI time series would be read in, and the 
%inter-voxel connectivity matrices would be computed.  If scan-rescan data 
%is available, one connectivity matrix would be computed for each scan (per 
%subject).  Otherwise, the time series of each subject would be split
%into two parts, and the connectivity matrix would be computed within each
%part.  In either case, two connectivity matrices per subject are required 
%by the shrinkIt function.

% RANDOMLY GENERATE CONNECTIVITY MATRICES MAT1 AND MAT2 FOR 20 SUBJECTS
% MAT1 AND MAT2 REPRESENT SCAN AND RESCAN (OR A SINGLE SCAN SPLIT IN TWO)

% - Suppose there are 100 voxels in the ROI
% - Size of connectivity matrix = 100 x 100
% - Elements in upper triangle = 100(100-1)/2 = 4950

V = 100;
m = V*(V-1)/2;
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
    
    %Apply inverse Fisher transformation to obtain correlation coefficients
    MAT1(:,:,i) = (exp(2*MAT1(:,:,i))-1)./(exp(2*MAT1(:,:,i))+1);
    MAT2(:,:,i) = (exp(2*MAT2(:,:,i))-1)./(exp(2*MAT2(:,:,i))+1);

end

%% RESHAPE DATA

% TRANSFORM ARRAYS MAT1 AND MAT2 INTO m-by-n MATRICES WITH mat2UT FUNCTION

% Extract the upper triangle from each V-by-V matrix.  Gives a vector
% of length m=100(100-1)/2.  Store these in X1 and X2.

X1 = zeros(m,n);
X2 = zeros(m,n);
for i = 1:20
    
    X1(:,i) = mat2UT(MAT1(:,:,i));
    X2(:,i) = mat2UT(MAT2(:,:,i));

end

%% PERFORM SHRINKAGE

% RUN shrinkIt TO OBTAIN SHRINKAGE ESTIMATE X (m-by-n MATRIX) AND SHRINKAGE PARAMETER

%Assume each subject has 10 minutes of scan time total
%(two 5-minute scans, or one 10-minute scan split in two)

[X lambda] = shrinkIt(X1, X2, 10);

size(X)
size(lambda)

%Look at distribution of shrinkage parameter values (lambda):
%
%lambda = var(noise)/[var(noise)+var(signal)]
%
%Since signal and noise were both generated with sd=0.1, lambda values 
%should be centered around 0.5.
hist(lambda)

%% PERFORM CLUSTERING

% Transform correlation to distance, reshape to V-by-V matrix with 0s on
% diagonal, and perform clustering.  

% kmedioids function allows user to input a custom distance matrix
% http://www.mathworks.com/matlabcentral/fileexchange/28860-kmedioids

clusters = zeros(V, n);

for i = 1:20
    dist = (1-X(:,i)/2);
    dist_mat = UT2mat(dist, 0);
    clusters(:,i) = kmedioids(dist_mat, 4);
end

