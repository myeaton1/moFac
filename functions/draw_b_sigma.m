function [PHI,SIG] = draw_b_sigma(r,xit,B,sig2,prior,params,var_fix)

N               =   size(r,2);
K               =   size(xit,2);

PHI             =   B;
SIG             =   sig2; 

for n = 1:N
    q_k             =   K;
    xPhi            =   zeros(K,1);
    if q_k > 0
        ystar           =   r(:,n);
        xstar           =   xit;
        TT              =   rows(ystar);

        % bring to right dimension
        prior_mean      =   prior.Phi.mean*ones(q_k,1);
        inv_prior_var   =   inv(prior.Phi.var)*eye(1);
        sig_inv         =   1/sig2(n);
            
        % posterior mean and variance
        post_var        =   inv(inv_prior_var + sig_inv*xstar'*xstar);
        post_mean       =   post_var*(inv_prior_var*prior_mean + sig_inv*xstar'*ystar);
            
        % draw from multivariate normal with posterior mean and variance
        C               =   chol(sig2(n)*post_var);
        psi_F1          =   post_mean+C'*randn(q_k,1);
        SSE             =   (ystar-xstar*psi_F1)'*(ystar-xstar*psi_F1);        
           
    else % if q_Fk >0
        
        TT              =   rows(xit);
        ystar           =   xit(:,n);
        SSE             =   ystar'*ystar;        
    
    end; % q_Fk >0
    
    if ~isempty(var_fix)
        sig2_F1         =   var_fix;     
    else
        d               =   prior.Sigma.shape + SSE;
        c               =   chi2rnd(TT + prior.Sigma.dof,1);
        sig2_F1         =   d/c;
    end % if DO_FIX_VARIANCE
    
    xPhi(1:q_k)     =   psi_F1(1:q_k);
    if n == 1
        SIG             =   sig2_F1;
        PHI             =   xPhi';
    else
        SIG             =   [SIG; sig2_F1];
        PHI             =   [PHI; xPhi'];
    end
        
end % for k = 1:K_F



