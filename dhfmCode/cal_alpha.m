% function cal_alpha computes time-varying intercept as function of factors at next higher level of the hierarchy

function [alpha,Lambda_tilde]=cal_alpha_new(F,Psi,Lambda);

K_F = size(F,2);
l_F = size(Lambda,2)-1;
q_G = size(Psi,2);
Lambda_lags = [];
for l=1:l_F+1;        
    Lambda_lags = [Lambda_lags Lambda{l}];
end;
l_star = l_F + q_G;

B = size(Lambda{1},1); 
Lambda_tilde = [];
alpha = [];
for b=1:B,    
    lambdab_tilde = [];
    for j=1:K_F;
        indx=j:K_F:K_F*(l_F+1);
        lambdabj_tilde = conv([1 -Psi(b,:)],Lambda_lags(b,indx));
        lambdab_tilde = [lambdab_tilde lambdabj_tilde'];
    end;
    lambdab_tilde = reshape(lambdab_tilde',1,K_F*(l_star+1));
    Lambda_tilde=[Lambda_tilde; lambdab_tilde];
    Flags = []; 
    for l = 0:l_star;
        Flags=[Flags lagn(F,l)];
    end;    
    alpha{b}=Flags*lambdab_tilde';    
end; % end b loop   
