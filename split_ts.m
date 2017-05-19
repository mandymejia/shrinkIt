function [X1 X2 Xodd Xeven] = split_ts(Y, b, d, fun)

%This function takes the time series data for a single subject and splits 
%it into four sub-series: 
%X1=fun(Y1), where Y1 consists of the first half of the time series
%X2=fun(Y2), where Y2 consists of the second half of the time series
%Xodd=fun(Yodd), where Yodd consists of the odd blocks of time points
%Xeven=fun(Yeven), where Yeven consists of the even blocks of time points
%multiplier is the factor by which the 

%Usage:
%   [X1 X2 Xeven Xodd] = split_ts(Y, b, d, fun)
%Inputs:
%   Y -  An t-by-p array or table, where t is the number of time points in the time
%   series and p is the number of observed variables.  
%
%   b -  A scalar indicating the block length for Yeven and Yodd.
%   For example, if b=3 (and d=0), Yodd will consist of the time points
%   (1:3,7:9,...) and Yeven will consist of (4:6,10:12,...)
%
%   d - A scalar indicating the spacing between blocks.  For example, if 
%   d=1 and b=3, Yodd will consist of the time points (1:3,9:11,...), 
%   and Yeven will consist of the time points (5:7,13:15,...)
%
%   fun - A character string or array of character strings with the name or
%   names of of the function(s) to be applied to Y1, Y2, Yodd, Yeven.  The
%   function or functions should be computed over time points, resulting in
%   a qx1 vector, where q is a function of p (e.g. q = 1, p, p^2).  If a
%   function requires additional arguments or a transformation of the data
%   Y, the user should create a user-defined function before calling split_ts.

%Outputs:
%   Y1 - (t1-by-p) matrix, t1 = floor(t/2).  Y1 = Y(1:t1,:)
%   Y2 - (t2-by-p) matrix, t2 = t - t1.  Y2 = Y((t1+1):end,:)
%   Yodd = (todd-by-p) matrix, consisting of the odd blocks of time points
%   Yeven = (teven-by-p) matrix, consisting of the even blocks of time points

%% Perform Checks

if(nargin ~= 4)
    error('Must specify four inputs')
end

if isempty(Y) || isempty(b) || isempty(fun)
    error('one or more inputs is empty')
end

if ~isnumeric(b) || max(size(b)) > 1 || b - round(b) ~= 0
    error('block length b must be an integer')
end      

if ~isnumeric(d) || max(size(d)) > 1 || d - round(d) ~= 0
    error('block gap length d must be an integer less than b')
end      


if ~ischar(fun) && ~iscellstr(fun)
    error('fun must be a string or cell array of strings')
end

%check that function names exist
if ischar(fun)
    fun_exist = exist(fun);
elseif iscellstr(fun)
    fun_exist = min(cellfun(@exist, fun));
end
if ~fun_exist
    error('one or more of the function names in "fun" input does not exist')
end

%% CREATE SUB-SERIES

if istable(Y) 
    Y = table2array(Y);
end

t = size(Y, 1);

%create Y1 and Y2
t1 = floor(t/2);
Y1 = Y(1:t1,:);
Y2 = Y((t1+1):(t1*2),:);

%create Yodd and Yeven
bd = b+d;
block1 = [ones(1,b), zeros(1,d), ones(1,b)*2, zeros(1,d)]; %[1,1,1,0,0,0] for b=3
inds = repmat(block1, [1,ceil(t/bd/2)]);
inds = inds(1,1:t); %truncate to length of time series
inds_odd = (inds==1); %indicator vector for odd blocks
inds_even = (inds==2);  %indicator vector for even blocks
Yodd = Y(find(inds_odd),:); %find() converts binary to indices
Yeven = Y(find(inds_even),:);

%% APPLY FUNCTION(S) TO SUB-SERIES

%apply first (or only) function to initialize X1, X2, Xodd, Xeven
fun = cellstr(fun); %convert to cell array of strings
fun0 = str2func(fun{1});

X1 = fun0(Y1);
X2 = fun0(Y2);
Xodd = fun0(Yodd);
Xeven = fun0(Yeven);

num_fun = numel(fun);
if(num_fun>1)
    
    for ii = 2:num_fun
        
        funi = str2func(fun{ii});
        X1 = funi(X1);
        X2 = funi(X2);
        Xodd = funi(Xodd);
        Xeven = funi(Xeven);
        
    end
    
end


