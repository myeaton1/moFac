clear, home, close all;

addpath('functions');
addpath('dhfmCode');
load('xSamp.mat');

% For testing for now, just use firms with data 
%   available for the entire sample

idx             =   ~any(isnan(xest));
r               =   xest(:,idx);
r               =   zscore(r);

n_burn=50000; n_keep=50000; n_skip=50; n_gibbs=n_burn+n_keep;

% Define model parameters
q               =   [1 1 1 1 1];
l               =   [0 0 0 0 0];            

sig_fix         =   .1;
noise_var       =   .2;

% Initialize factors as PC estimates plus random noise
T               =   size(r,1);
K               =   5;
noise_var       =   .2;

[ehat,xi,lambda,ss] ...
                =   pc_T(r,K);
            
xi              =   xi + noise_var*randn(T,K);
for k           =   1:K,
    xi(:,k)         =   sign(corr(xi(:,k),r(:,k)))*xi(:,k);
end;

% Now initialize parameters based on starting values for the factors

% First regress starting values of aggregate factors on their own lags
% to obtain starting values for the psi_F's and sig2_F's
startprior.Phi.mean = 0; startprior.Phi.var = 1e10;
startprior.Sigma.shape = 0; startprior.Sigma.dof = 0;
paramsXi.q = q; paramsXi.l = l;
[Phi,sig2_Phi,mu]=   draw_phi_sigma(xi,.5*ones(K,1),ones(K,1),startprior,paramsXi,sig_fix,true);








