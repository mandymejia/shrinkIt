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

%% ADD PATH TO SHRINKIT TOOLBOX AND KMEDIOIDS TOOLBOX

addpath '~/matlab_toolboxes/shrinkIt/'
addpath '~/matlab_toolboxes/kmedioids/'


%% POINT TO DATA AND GET FILE NAMES

data_dir = '/dcs01/oasis/hpc/HCP500_Parcellation_Timeseries_Netmats/Results/node_timeseries/3T_Q1-Q6related468_MSMsulc_d300_ts2';
fnames = dir(data_dir); %get list of files
fnames = fnames(3:end); %remove . and ..
n = numel(fnames); %number of subjects



%% SET FUNCTION(S) TO APPLY TO DATA

%single function example: compute VxV correlation matrix
fun_single = 'corrcoef';

%multiple function example: compute upper triangle of VxV correlation matrix
fun_multiple = {'corrcoef', 'mat2UT'};

%length of blocks in even/odd splitting
b = 10; 



%% LOOP THROUGH SUBJECTS TO CREATE DATA MATRICES

cd(data_dir)
for ii = 1:n
    
    ii
    
    % READ IN TIME SERIES DATA
    
    fnamei = fnames(ii);
    fnamei = fnamei.name;
    
    %data stored as text files
    Yi = readtable(fnamei, 'Delimiter',' ','ReadVariableNames', false);
        
    
    % SPLIT TIME SERIES DATA USING SPLIT_TS()
    % COMPUTE UT OF CORRELATION MATRIX
    
    %single function example: 
    %split data and compute VxV correlation matrix for each split
    [X1i X2i Xoddi Xeveni] = split_ts(Yi, b, fun_single);
    
    %multiple function example:
    %split data and compute upper triangle of VxV correlation matrix for each split
    [X1i X2i Xoddi Xeveni] = split_ts(Yi, b, fun_multiple);
    
    
    
    % COMBINE SUBJECTS
    
    if(ii==1) 
        %initialize group arrays for each split
        [X1 X2 Xodd Xeven] = deal(zeros([size(X1i), n]));
    end

    %add current subject to group arrays:
    %X's have dimensions (p1,p2,...,pk,n), so use (:,:,...,:,ii) to index subject ii 
    index = repmat({':'},1,ndims(X1)); %(:,:,...,:,:)
    index{end} = ii;                   %(:,:,...,:,ii)

    X1(index{:}) = X1i;
    X2(index{:}) = X2i;
    Xodd(index{:}) = Xoddi;
    Xeven(index{:}) = Xeveni;
    
end

%remove superfluous dimensions (of size 1)
X1 = squeeze(X1);
X2 = squeeze(X2);
Xodd = squeeze(Xodd);
Xeven = squeeze(Xeven);



%% PERFORM SHRINKAGE

% RUN shrinkIt TO OBTAIN SHRINKAGE ESTIMATE X (m-by-n MATRIX) AND SHRINKAGE PARAMETER

%Assume each subject has 10 minutes of scan time total
%(two 5-minute scans, or one 10-minute scan split in two)

[X_shrink lambda] = shrinkIt(X1, X2, Xodd, Xeven);

size(X_shrink)
size(lambda)

%Look at distribution of shrinkage parameter values (lambda):
%
%lambda = var(noise)/[var(noise)+var(signal)]
%

hist(lambda)


%% PERFORM CLUSTERING OF VOXELS OF EACH SUBJECT USING SHRINKAGE ESTIMATES OF CORRELATION MATRICES

% Transform correlation to distance, reshape to V-by-V matrix with 0s on
% diagonal, and perform clustering.  

% kmedioids function allows user to input a custom distance matrix
% http://www.mathworks.com/matlabcentral/fileexchange/28860-kmedioids

% cluster labels may differ across subjects

V = size(Yi,2);
clusters = zeros(V, n);

for ii = 1:n
    dist = (1-X_shrink(:,ii)/2) - 0.5;
    dist_mat = UT2mat(dist, 0);
    clusters(:,ii) = kmedioids(dist_mat, 4); 
end

