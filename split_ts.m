function [X1 X2] = split_ts(Y, fun)

%This function takes the time series data for a single subject and splits 
%it into two sub-series, then applies the function(s) in fun to each
%X1=fun(Y1), where Y1 consists of the first half of the time series
%X2=fun(Y2), where Y2 consists of the second half of the time series
%
%Usage:
%   [X1 X2] = split_ts(Y, fun)
%Inputs:
%   Y -  An t-by-p array or table, where t is the number of time points in the time
%   series and p is the number of observed variables.  
%
%   fun - A character string or array of character strings with the name or
%   names of of the function(s) to be applied to Y1 and Y2.  The
%   function or functions should be computed over time points, resulting in
%   a qx1 vector, where q is a function of p (e.g. q = 1, p, p^2).  If a
%   function requires additional arguments or a transformation of the data
%   Y, the user should create a user-defined function before calling split_ts.

%Outputs:
%   X1 = fun(Y1) (qx1), where Y1 is a (t1-by-p) matrix, t1 = floor(t/2).  Y1 = Y(1:t1,:)
%   X2 = fun(Y2) (qx1), where Y2 is a (t1-by-p) matrix, Y2 = Y((t-t1+1):t,:)

%% Perform Checks

if(nargin ~= 2)
    error('Must specify two inputs')
end

if isempty(Y) || isempty(fun)
    error('one or more inputs is empty')
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
Y2 = Y((t-t1+1):t,:);

%% APPLY FUNCTION(S) TO SUB-SERIES

%apply first (or only) function to initialize X1 and X2
fun = cellstr(fun); %convert to cell array of strings
fun0 = str2func(fun{1});

X1 = fun0(Y1);
X2 = fun0(Y2);

num_fun = numel(fun);
if(num_fun>1)
    
    for ii = 2:num_fun
        
        funi = str2func(fun{ii});
        X1 = funi(X1);
        X2 = funi(X2);
        
    end
    
end


