%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Replication code for the paper "Dynamic Hierarchical Factor Models"                                         %%%
%%%                                                                                                             %%%
%%% Copyright: Emanuel Moench (emanuel.moench@ny.frb.org), Serena Ng(serena.ng@columbia.edu, and                %%%
%%%            Simon Potter (simon.potter@ny.frb.org                                                            %%% 
%%%                                                                                                             %%%
%%% Implements a Gibbs sampler to estimate the following factor model:                                          %%%
%%%                                                                                                             %%%
%%%     Z_bsnt = Lambda_Hbsn(L)'*H_bst + e_Zbsnt                                                                %%%             
%%%      H_bst = Lambda_Gbs(L)'*G_bt + e_Hbst                                                                   %%%
%%%       G_bt = Lambda_Fb(L)'*F_t + e_Gbt                                                                      %%%
%%%    e_Zbsnt = psi_Zbsn(L)*e_Zbsn_{t-1} + epsilon_Zbsnt                                                       %%% 
%%%     e_Hbst = psi_Hbs(L)*e_Hbs_{t-1} + epsilon_Hbst                                                          %%% 
%%%      e_Gbt = psi_Gb(L)*e_Gb_{t-1} + epsilon_Gbt                                                             %%% 
%%%        F_t = psi_F*F_{t-1} + epsilon_Ft                                                                     %%%
%%%                                                                                                             %%% 
%%% B is the number of blocks, Bsub(b) is the number of subblocks in block b,                                   %%%  
%%% Nsub is the number of series in subblock s of block b,                                                      %%%
%%%                                                                                                             %%%
%%% K_F is the number of common factors                                                                         %%%
%%% K_G(b) is the number of block-specific factors in block b,                                                  %%%
%%% K_H{b}(s) is the number of subblock-specific factors in subblock s of block b                               %%%
%%%                                                                                                             %%%
%%% Note that if for some b, K_H{b} is zero, then the observation equation becomes                              %%%
%%% Z_bsnt = Lambda_Gbn(L)'*G_bt + e_Zbnt                                                                       %%%
%%%                                                                                                             %%%
%%% q_F is the lag order of the transition equation F_t = psi_F*F_{t-1} + epsilon_Ft                            %%%
%%% q_G(b) is the lag order of the AR model for e_Gbt,                                                          %%%
%%% q_H(b) is the lag order of the AR model for all subblock-specific components e_Hbst                         %%%
%%% q_Z(b) is the lag order of the AR model for all idiosyncratic components e_Zbsnt                            %%%
%%%                                                                                                             %%%
%%% l_F(b) is the lag order of the factor loading polynomial Lambda_Fb(L)                                       %%%
%%% l_G(b) is the lag order of the factor loading polynomials Lambda_Gbs(L)                                     %%%
%%% l_H(b) is the lag order of the factor loading polynomials Lambda_Hbsn(L)                                    %%%
%%%                                                                                                             %%%
%%% Note that we set these to zero and restrict loading matrices to be lower-triangular with ones on the        %%% 
%%% diagonal in order to identify the factors; In principle, the loadings can be dynamic but then other         %%%
%%% restrictions need to be imposed in order to identify the factors and loadings                               %%%
%%%                                                                                                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all;

n_burn=50000; n_keep=50000; n_skip=50; n_gibbs=n_burn+n_keep;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B=5;  % number of blocks 
K_F=1;
K_G=[1; 1; 1; 1; 1];
K_H{1}=[1 1 2];
K_H{2}=[2 1];
K_H{3}=[1 1];
K_H{4}=[0]; % K_H{4}=[1 1];
K_H{5}=[0]; % K_H{5}=[1 1 1];

q_F=1;                 
q_G=[1 1 1 1 1];      
q_H=[1 1 1 1 1];
q_Z=[1 1 1 1 1];

l_F=[0 0 0 0 0];            
l_G=[0 0 0 0 0];
l_H=[0 0 0 0 0];

sig_fix = .1;
noise_var = .2;

matname=['dhfm_results'];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'dhfm_data'

% Bsub(b) is the number of sub-blocks in block b
% Nsub(b,s) is the size of subblock-j of block b
% K_H{b}(s) is the number of factors in sub-block s of block b
NB=[];
for b = 1:B,
    Bsub(b) = length(K_H{b});
    if K_H{b} == 0, Bsub(b) = 0; end;        
    for s=1:Bsub(b);
        Nsub(b,s)=ZNB{b}(s);
    end;
    NB=[NB; sum(K_H{b})];
end;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize sampler
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dhfm_initialize_sampler;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run sampler
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dhfm_run_sampler;



