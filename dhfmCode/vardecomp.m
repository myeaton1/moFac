function shares=cal_vcv_new3(Lambda_F,Lambda_Gsub,Lambda_G,Lambda_H,psi_F,psi_G,psi_H,psi_Z,...
			sig2_F,sig2_G,sig2_H,sig2_Z,params);

B=params.B;
Nsub=params.Nsub;
Bsub=params.Bsub;
K_F=params.K_F;
K_G=params.K_G;
K_H=params.K_H;
q_F=params.q_F;
q_G=params.q_G;
q_H=params.q_H;
q_Z=params.q_Z;
l_F=params.l_F;
l_G=params.l_G;
l_H=params.l_H;

shares=[];
% for blocks with subblocks
for b=1:B;
    
    if b == 1,
        Gb_indx = 1:K_G(b);
    else,
        Gb_indx = sum(K_G(1:b-1))+1:sum(K_G(1:b));        
    end;
    k_b = K_G(b); 
    
    Lambda_F_mat = cell2mat(Lambda_F);
    Lambda_Fb = Lambda_F_mat(Gb_indx,:);
    sum_lamdaFb = zeros(k_b^2,K_F^2);
    for l = 1:l_F(b)+1,
        sum_lamdaFb = sum_lamdaFb + kron(Lambda_Fb(:,(l-1)*K_F+1:l*K_F),Lambda_Fb(:,(l-1)*K_F+1:l*K_F));
    end;   
    
    sum_psiF = zeros(K_F^2);
    for q = 1:q_F;
       sum_psiF = sum_psiF + kron(diag(psi_F(:,q)),diag(psi_F(:,q)));
    end;
    Isum_psiF = eye(K_F^2) - sum_psiF;    
    
    % Putting together the Psis
    Psi_Gb = psi_G{b};
    sum_psiGb = zeros(k_b^2);
    for q = 1:q_G(b),
        sum_psiGb = sum_psiGb + kron(diag(Psi_Gb(:,q)),diag(Psi_Gb(:,q)));
    end;
    Isum_psiGb = eye(k_b^2) - sum_psiGb;

    
    if Bsub(b) == 0,
        
        Psi_Zb = psi_Z{b};
        N_b = length(Psi_Zb);
        Isum_psiZb = zeros(N_b,1);
        for i = 1:N_b,
            sum_psiZbi = 0;
            for q = 1:q_Z(b),
                sum_psiZbi = sum_psiZbi + Psi_Zb(i,q).^2;
            end;
            Isum_psiZb(i) = 1 - sum_psiZbi;
        end;
        
        Lambda_Gb = cell2mat(Lambda_G{b});
        sum_lamdaGb = zeros(N_b,k_b^2);
        for i = 1:N_b,
            for l = 1:l_G(b)+1;
                sum_lamdaGb(i,:) = sum_lamdaGb(i,:) + kron(Lambda_Gb(i,(l-1)*k_b+1:l*k_b),Lambda_Gb(i,(l-1)*k_b+1:l*k_b));
            end;
        end;
        
        % Loadings decomposition
        gamma_F{b} = sum_lamdaGb*sum_lamdaFb;
        gamma_G{b} = sum_lamdaGb;
        gamma_Z{b} = 1;

        % Variances of the different components
        v_F=Isum_psiF\vec(diag(sig2_F));
        v_Gb{b}=Isum_psiGb\vec(diag(sig2_G{b}));      
        v_Zb{b}=diag(1./Isum_psiZb)*sig2_Z{b};

        Var_F{b} = gamma_F{b}*v_F;
        Var_G{b} = gamma_G{b}*v_Gb{b};
        Var_Z{b} = gamma_Z{b}*v_Zb{b};

        total_Var{b} = Var_F{b} + Var_G{b} + Var_Z{b};

        share_F{b} = Var_F{b}./total_Var{b};
        share_G{b} = Var_G{b}./total_Var{b};
        share_Z{b} = Var_Z{b}./total_Var{b};    

        shares{b} = [total_Var{b} share_F{b} share_G{b} share_Z{b}];
                
    else, % if Bsub(b) == 0

        % for subblocks
        for s=1:Bsub(b);       
            N_bs = Nsub(b,s);
            k_bs = K_H{b}(s);

            if s == 1,
                Hbs_indx = 1:K_H{b}(s);
            else,
                Hbs_indx = sum(K_H{b}(1:s-1))+1:sum(K_H{b}(1:s));        
            end;

            % Putting together the lambdas
            Lambda_Gsub_mat = cell2mat(Lambda_Gsub{b});
            Lambda_Gbs = Lambda_Gsub_mat(Hbs_indx,:);
            sum_lamdaGbs = zeros(k_bs^2,k_b^2);
            for l = 1:l_G(b)+1,
                sum_lamdaGbs = sum_lamdaGbs + kron(Lambda_Gbs(:,(l-1)*k_b+1:l*k_b),Lambda_Gbs(:,(l-1)*k_b+1:l*k_b));
            end;

            Lambda_Hbs = cell2mat(Lambda_H{b}{s});
            sum_lamdaHbs = zeros(N_bs,k_bs^2);
            for i = 1:N_bs,
                for l = 1:l_H(b)+1;
                    sum_lamdaHbs(i,:) = sum_lamdaHbs(i,:) + kron(Lambda_Hbs(i,(l-1)*k_bs+1:l*k_bs),Lambda_Hbs(i,(l-1)*k_bs+1:l*k_bs));
                end;
            end;

            % Putting together the Psis
            Psi_Hbs = psi_H{b}{s};
            sum_psiHbs = zeros(k_bs^2);
            for q = 1:q_H(b),
                sum_psiHbs = sum_psiHbs + kron(diag(Psi_Hbs(:,q)),diag(Psi_Hbs(:,q)));
            end;
            Isum_psiHbs = eye(k_bs^2) - sum_psiHbs;

            Psi_Zbs = psi_Z{b}{s};
            Isum_psiZbs = zeros(N_bs,1);
            for i = 1:N_bs,
                sum_psiZbsi = 0;
                for q = 1:q_Z(b),
                    sum_psiZbsi = sum_psiZbsi + Psi_Zbs(i,q).^2;
                end;
                Isum_psiZbs(i) = 1 - sum_psiZbsi;
            end;
                        
            % Loadings decomposition
            gamma_F{b}{s} = sum_lamdaHbs*sum_lamdaGbs*sum_lamdaFb;
            gamma_G{b}{s} = sum_lamdaHbs*sum_lamdaGbs;
            gamma_H{b}{s} = sum_lamdaHbs;
            gamma_Z{b}{s} = 1;

            % Variances of the different components
            v_F=Isum_psiF\vec(diag(sig2_F));
            v_Gb{b}{s}=Isum_psiGb\vec(diag(sig2_G{b}));
            v_Hbs{b}{s}=Isum_psiHbs\vec(diag(sig2_H{b}{s}));
            v_Zbs{b}{s}=diag(1./Isum_psiZbs)*sig2_Z{b}{s};

            Var_F{b}{s} = gamma_F{b}{s}*v_F;
            Var_G{b}{s} = gamma_G{b}{s}*v_Gb{b}{s};
            Var_H{b}{s} = gamma_H{b}{s}*v_Hbs{b}{s};
            Var_Z{b}{s} = gamma_Z{b}{s}*v_Zbs{b}{s};
            
            total_Var{b}{s} = Var_F{b}{s} + Var_G{b}{s} + Var_H{b}{s} + Var_Z{b}{s};

            share_F{b}{s} = Var_F{b}{s}./total_Var{b}{s};
            share_G{b}{s} = Var_G{b}{s}./total_Var{b}{s};
            share_H{b}{s} = Var_H{b}{s}./total_Var{b}{s};
            share_Z{b}{s} = Var_Z{b}{s}./total_Var{b}{s};    

            shares{b}{s}  = [total_Var{b}{s} share_F{b}{s} share_G{b}{s} share_H{b}{s} share_Z{b}{s}];
            
        end; % for s = 1:Bsub(b);
    
    end; % if Bsub(b) == 0

end; % for b = 1:B




