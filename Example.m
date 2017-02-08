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

addpath '~/matlab_toolboxes/shrinkIt/'


%% POINT TO DATA AND GET FILE NAMES

data_dir = '/dcs01/oasis/hpc/HCP500_Parcellation_Timeseries_Netmats/Results/node_timeseries/3T_Q1-Q6related468_MSMsulc_d300_ts2';
fnames = dir(data_dir); %get list of files
fnames = fnames(3:end); %remove . and ..
n = numel(fnames); %number of subjects



%% SET FUNCTION(S) TO APPLY TO DATA

%single function example: compute VxV correlation matrix
fun_single = 'corrcoef';

%multiple function example: compute upper triangle of VxV correlation matrix, then Fisher-transform
fun_multiple = {'corrcoef', 'mat2UT', 'fish'};

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
    Yi1 = Yi(1:1200,:); %visit 1, LR acquisition
    Yi2 = Yi(2401:3600,:); %visit 2, LR acquisition
    
    % SPLIT TIME SERIES DATA USING SPLIT_TS()
    % COMPUTE UT OF CORRELATION MATRIX
    
    %single function example: 
    %split data and compute VxV correlation matrix for each split
    [X1i X2i Xoddi Xeveni] = split_ts(Yi1, b, fun_single);
    
    %multiple function example:
    %split data and compute upper triangle of VxV correlation matrix for each split
    [X1i X2i Xoddi Xeveni] = split_ts(Yi1, b, fun_multiple);


    % COMPUTE ESTIMATE FROM SECOND VISIT FOR RELIABILITY ANALYSIS
    X_visit2i = mat2UT(corrcoef(table2array(Yi2)));

    
    % COMBINE SUBJECTS
    
    if(ii==1) 
        %initialize group arrays for each split
        [X1 X2 Xodd Xeven X_visit2] = deal(zeros([size(X1i), n]));
    end

    %add current subject to group arrays:
    %X's have dimensions (p1,p2,...,pk,n), so use (:,:,...,:,ii) to index subject ii 
    index = repmat({':'},1,ndims(X1)); %(:,:,...,:,:)
    index{end} = ii;                   %(:,:,...,:,ii)

    X1(index{:}) = X1i;
    X2(index{:}) = X2i;
    Xodd(index{:}) = Xoddi;
    Xeven(index{:}) = Xeveni;
    X_visit2(index{:}) = X_visit2i;
    
end

%remove superfluous dimensions (of size 1)
X1 = squeeze(X1);
X2 = squeeze(X2);
Xodd = squeeze(Xodd);
Xeven = squeeze(Xeven);
X_visit2 = squeeze(X_visit2);



%% PERFORM SHRINKAGE

% RUN shrinkIt TO OBTAIN SHRINKAGE ESTIMATE X (m-by-n MATRIX) AND SHRINKAGE PARAMETER

%Assume each subject has 10 minutes of scan time total
%(two 5-minute scans, or one 10-minute scan split in two)

[X_shrink lambda varU varW varX] = shrinkIt(X1, X2, Xodd, Xeven);
X_shrink = unfish(X_shrink); %inverse Fisher-transform to obtain Pearson correlations


% VISUALIZE VARIANCE COMPONENTS & DEGREE OF SHRINKAGE

varWTHN = varU + varW;
varBTWN = varX;
maxv = max(max(varWTHN(:)),max(varBTWN(:)))*.5;

figure
h = subplot(1,3,1);
image(UT2mat(varWTHN, 0),'CDataMapping','scaled')
colorbar; caxis([0,maxv]); axis off;
title('Within-Subject Variance')
h = subplot(1,3,2);
image(UT2mat(varBTWN, 0),'CDataMapping','scaled')
colorbar; caxis([0,maxv]); axis off;
title('Between-Subject Variance')
h = subplot(1,3,3);
image(UT2mat(lambda, 0),'CDataMapping','scaled')
colorbar; caxis([0,1]); axis off;
title('Degree of Shrinkage')




% VISUALIZE RSFC ESTIMATES & RELIABILITY

X1 = unfish(X1);
X2 = unfish(X2);
X_avg = (X1 + X2)./2;

figure
h = subplot(3,3,1);
image(UT2mat(X_avg(:,1), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 1 Raw')
h = subplot(3,3,2);
image(UT2mat(X_avg(:,2), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 2 Raw')
h = subplot(3,3,3);
image(UT2mat(X_avg(:,3), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 3 Raw')
h = subplot(3,3,4);
image(UT2mat(X_shrink(:,1), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 1 Shrink')
h = subplot(3,3,5);
image(UT2mat(X_shrink(:,2), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 2 Shrink')
h = subplot(3,3,6);
image(UT2mat(X_shrink(:,3), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 3 Shrink')
h = subplot(3,3,7);
image(UT2mat(X_visit2(:,1), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 1 Visit 2')
h = subplot(3,3,8);
image(UT2mat(X_visit2(:,2), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 2 Visit 2')
h = subplot(3,3,9);
image(UT2mat(X_visit2(:,3), 1),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-.5,.5]); axis off;
title('Subject 3 Visit 2')



% COMPUTE VISIT 2 MSE OF RAW AND SHRINKAGE ESTIMATES

MSE_raw = mean((X_avg - X_visit2).^2)
MSE_shrink = mean((X_shrink - X_visit2).^2)

%Percent change in MSE due to shrinkage (positive = reduction in MSE = improved reliability)
(MSE_raw - MSE_shrink)./MSE_raw
mean((MSE_raw - MSE_shrink)./MSE_raw)

sqerr_raw = (X_avg - X_visit2).^2;
sqerr_shrink = (X_shrink - X_visit2).^2;
sqerr_change = 100*(sqerr_raw - sqerr_shrink)./sqerr_raw;
sqerr_change = median(sqerr_change, 2);

figure
image(UT2mat(sqerr_change, 0),'CDataMapping','scaled')
colormap parula; colorbar; caxis([-50,50]); axis off;
title('% Change in Reliability')

