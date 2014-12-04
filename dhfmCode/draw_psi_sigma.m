function [PSI,SIG] = draw_psi_sigma2(Ft,psi_F,sig2_F,prior,params,var_fix);

K_F = size(Ft,2);
l_F = params.l_F;
q_F = params.q_F;
max_q_F = max(q_F);
max_l_F = max(l_F);

T=rows(Ft);
PSI=psi_F; SIG=sig2_F; 

for k=1:K_F;
    q_Fk=q_F(k);
    xpsi_F1=zeros(max_q_F,1);
    if q_Fk > 0;
        ystar=Ft(:,k);
        xstar=[];
        for i=1:q_Fk;
            xstar=[xstar lagn(Ft(:,k),i)];
        end;
        ystar=trimr(ystar,max(max_q_F,max_l_F),0);
        xstar=trimr(xstar,max(max_q_F,max_l_F),0);
        TT=rows(ystar);

        % bring to right dimension
        prior_mean = prior.Psi.mean*ones(q_Fk,1);
        inv_prior_var = inv(prior.Psi.var)*eye(q_Fk);
        sig_inv = 1/sig2_F(k);
            
        % posterior mean and variance
        post_var = inv(inv_prior_var + sig_inv*xstar'*xstar);
        post_mean = post_var*(inv_prior_var*prior_mean + sig_inv*xstar'*ystar);
            
        % draw from multivariate normal with posterior mean and variance
        C=chol(sig2_F(k)*post_var);
        accept=0; counter=0;
        while accept==0 & counter < 1000;
            psi_F1=post_mean+C'*randn(q_Fk,1);
            ceof= flipud([1; -psi_F1]);
            root=roots(ceof);
            rootmod=abs(root);
            if min(rootmod) > 1.0001;
                accept=1; 
            else; 
                counter=counter+1;
                accept=0;
            end;
        end; % while         
        if ~accept & counter >= 1000,
            disp(counter);
            return;
        end;
        SSE = (ystar-xstar*psi_F1)'*(ystar-xstar*psi_F1);        
           
    else, % if q_Fk >0
        
        TT=rows(Ft);
        ystar=Ft(:,k);
        SSE = ystar'*ystar;        
    
    end; % q_Fk >0
    
    if ~isempty(var_fix),
        sig2_F1= var_fix;     
    else,
        d=prior.Sigma.shape + SSE;
        c=chi2rnd(TT + prior.Sigma.dof,1);
        sig2_F1=d/c;
    end; % if DO_FIX_VARIANCE
    
    xpsi_F1(1:q_Fk)=psi_F1(1:q_Fk);
    if k == 1,
        SIG = sig2_F1;
        PSI = xpsi_F1';
    else,
        SIG=[SIG; sig2_F1];
        PSI=[PSI; xpsi_F1'];
    end;
        
end; % for k = 1:K_F



