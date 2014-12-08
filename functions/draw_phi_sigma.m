function [PHI,SIG,MU] = draw_phi_sigma(xit,phi,sig2,prior,params,var_fix,const)

K               =   size(xit,2);
l               =   params.l;
q               =   params.q;
max_q           =   max(q);
max_l           =   max(l);

T               =   rows(xit);
PHI             =   phi;
SIG             =   sig2; 

for k = 1:K
    q_k             =   q(k);
    xPhi            =   zeros(max_q,1);
    if q_k > 0
        ystar           =   xit(:,k);
        xstar           =   [];
        for i = 1:q
            xstar           =   [xstar lagn(xit(:,k),i)];
        end
        if const
            xstar           =   [ones(T,1), xstar];
            q_k             =   q_k + 1;
        end
        ystar           =   trimr(ystar,max(max_q,max_l),0);
        xstar           =   trimr(xstar,max(max_q,max_l),0);
        TT              =   rows(ystar);

        % bring to right dimension
        prior_mean      =   prior.Phi.mean*ones(q_k,1);
        inv_prior_var   =   inv(prior.Phi.var)*eye(q_k);
        sig_inv         =   1/sig2(k);
            
        % posterior mean and variance
        post_var        =   inv(inv_prior_var + sig_inv*xstar'*xstar);
        post_mean       =   post_var*(inv_prior_var*prior_mean + sig_inv*xstar'*ystar);
            
        % draw from multivariate normal with posterior mean and variance
        C               =   chol(sig2(k)*post_var);
        accept          =   0;
        counter         =   0;
        while accept == 0 && counter < 1000
            psi_F1          =   post_mean+C'*randn(q_k,1);
            ceof            =   flipud([1; -psi_F1]);
            root            =   roots(ceof);
            rootmod         =   abs(root);
            if min(rootmod) > 1.0001
                accept          =   1; 
            else
                counter         =   counter + 1;
                accept          =   0;
            end
        end % while         
        if ~accept && counter >= 1000
            disp(counter);
            return;
        end
        SSE             =   (ystar-xstar*psi_F1)'*(ystar-xstar*psi_F1);        
           
    else % if q_Fk >0
        
        TT              =   rows(xit);
        ystar           =   xit(:,k);
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
    if k == 1
        SIG             =   sig2_F1;
        if const
            PHI             =   xPhi(2:end)';
            muPhi           =   xPhi(1)';
        else
            PHI             =   xPhi';
        end
    else
        SIG             =   [SIG; sig2_F1];
        if const
            PHI             =   [PHI; xPhi(2:end)'];
            muPhi           =   [muPhi; xPhi(1)'];
        else
            PHI             =   [PHI; xPhi'];
        end
    end
        
end % for k = 1:K_F

MU              =   (eye(size(PHI,1)) - diag(PHI))\muPhi;



