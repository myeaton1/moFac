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

% Initialize lambda as PC estimates plus random noise
T               =   size(r,1);
K               =   5;
N               =   size(r,2);
noise_var       =   .2;

[ehat,lambda,loading,ss] ...
                =   pc_T(r,K);
            
lambda          =   lambda + noise_var*randn(T,K);
for k           =   1:K,
    lambda(:,k)     =   sign(corr(lambda(:,k),r(:,k)))*lambda(:,k);
end;

% Initialize f as random noise
f               =   noise_var*randn(T,K);

% Now initialize parameters based on starting values for the factors

% First regress starting values of aggregate factors on their own lags
% to obtain starting values for Phi, Sig_v, and mu, initialize Sig_f
startprior.Phi.mean = 0; startprior.Phi.var = 1e10;
startprior.Sigma.shape = 0; startprior.Sigma.dof = 0;
paramsLambda.q = q; paramsLambda.l = l;
[Phi,Sig_v,mu]  =   draw_phi_sigma(lambda,.5*ones(K,1),ones(K,1),startprior,paramsLambda,sig_fix,true);


Sig_f           =   repmat(sig_fix, K, 1);

z               =   lambda(1:end-1,:) + f(2:end,:);
rt              =   r(2:end,:);
[B,Sig_e]       =   draw_b_sigma(rt,z,.5*ones(K,N),ones(N,1),startprior,paramsLambda,sig_fix);

xi              =   [ lambda(2:end,:), lambda(1:end-1,:), f(2:end,:) ];
H               =   [ zeros(size(B)), B, B ];
F               =   [ [ diag(Phi), zeros(K), zeros(K) ] ; [ eye(K), zeros(K), zeros(K) ]; ...
                        [ zeros(K), zeros(K), zeros(K) ] ];
Sig_eMat        =   diag(Sig_e);
Sig             =   [ [ diag(Sig_v), zeros(K), zeros(K) ] ; [ zeros(K), zeros(K), zeros(K) ]; ...
                        [ zeros(K), zeros(K), diag(Sig_f) ] ];
alpha           =   [ (eye(K) - diag(Phi))*mu ; zeros(K,1); zeros(K,1) ];

jj=1;
% Priors for parameters related to F
prior.F.Lambda.mean = 0;
prior.F.Lambda.var = 10;

prior.F.Psi.mean = 0.5;
prior.F.Psi.var = 10;

prior.F.Sigma.shape = .01;
prior.F.Sigma.dof = 4;

randn('state',999);  
for gibbs=1:n_gibbs,

    % Update factors and corresponding parameters using last iteration's draws
    %   of params and factors 
%     xi              =   sample_facs(r,alpha,H,psi_F,sig2_F,psi_G,sig2_G,paramsF);    

%     [alphaF,lambda_F_tilde]=cal_alpha(f,cell2mat(psi_G'),Lambda_F);
%     [psi_F,sig2_F]=draw_psi_sigma(f,psi_F,sig2_F,prior.F,paramsF,sig_fix);
    
    [Phi,Sig_v,mu]  =   draw_phi_sigma(lambda,Phi,Sig_v,startprior,paramsLambda,[],true);
    
    if gibbs > n_burn && mod(gibbs,n_skip)==0,   
        PHI{jj}         =   Phi;
        SIG_V{jj}       =   Sig_v;
        jj              =   jj + 1;
    end;  % end storing results
    if mod(gibbs,250) == 0;
        disp(sprintf('Completed %d Draws',gibbs));
    end
    
end
