clear, home, close all;

load('xSamp.mat');

% For testing for now, just use firms with data 
%   available for the entire sample

idx             =   ~any(isnan(xest));
r               =   xest(:,idx);

