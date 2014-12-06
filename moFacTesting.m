clear, home, close all;

addpath('functions');
addpath('dhfmCode');
load('xSamp.mat');

% For testing for now, just use firms with data 
%   available for the entire sample

idx             =   ~any(isnan(xest));
r               =   xest(:,idx);

% Initialize factors as PC estimates plus random noise
T               =   size(r,1);
K               =   5;

[ehat1,f,lambda_G1,ss1] ...
                =   pc_T(r,K);
[coeff, score]  =   princomp(zscore(r));   
score           =   score(:,1:K);
            
f = f + noise_var*randn(T,K_F);
for k = 1:K_F,
    f(:,k) = sign(corr(f(:,k),allG(:,k)))*f(:,k);
end;