function LAMBDA = draw_lambda_test(G,F,psi_G,SIG2_G,prior,params);

B = params.B;
q_G = params.q_G;
max_q_G = max(q_G);
K_G = params.K_blocks;
K_F = params.K_facs;
l_F = params.l_F;
max_l_F = max(l_F);
Lamda_F = [];
count = 0; 
for b=1:B,
    
    Gb=G{b};
    [T,k_G]=size(Gb);
    l_Fb = l_F(b);
        
    for i=1:k_G,
        
        count = count + 1; 
            
        F_filt=filter([1 -psi_G{b}(i,:)],1,F);
        ystar=filter([1 -psi_G{b}(i,:)],1,Gb(:,i));
        Fstar=F_filt;
        for j=1:l_Fb;
            Fstar=[Fstar lagn(F_filt,j)];
        end;
        Fstar=trimr(Fstar,max(max_q_G,l_Fb),0);
        ystar=trimr(ystar,max(max_q_G,l_Fb),0);

        % bring to right dimension
        prior_mean = prior.mean*ones(K_F*(l_Fb+1),1);
        prior_var = prior.var*eye(K_F*(l_Fb+1));
        sig_inv = 1/SIG2_G{b}(i);

        % posterior mean and variance
        post_var = inv(inv(prior_var) + sig_inv*Fstar'*Fstar);
        post_mean = post_var*(inv(prior_var)*prior_mean + sig_inv*Fstar'*ystar);

        % draw from multivariate normal with posterior mean and variance
        C=chol(post_var);
        lamda_Fi=post_mean + C'*randn(K_F*(l_Fb+1),1);
        Lamda_F=[Lamda_F; lamda_Fi' zeros(1,max_l_F-l_Fb)];
        
    end; % end i

end; % for b=1:B

% Impose lower triangularity with ones on the diagonal for identification
for l = 1:max_l_F+1,
    LAMBDA{l} = Lamda_F(:,(l-1)*K_F+1:l*K_F);        
    for j = 1:K_F,
        if l == 1,
            LAMBDA{l}(j,j) = 1;
        end;
        LAMBDA{l}(j,j+1:K_F) = zeros(1,K_F-j); % make loadings lower triangular at all lags  
    end; 
end;
    








