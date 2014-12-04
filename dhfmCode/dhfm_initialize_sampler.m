% Preliminary PC analysis 
[F_pc,G_pc,H_pc,FH_pc,FX_pc,GX_pc,K_F_pc,K_G_pc,K_H_pc]=do_PC_level4(bigZ_ns,Bsub,8,K_F,K_G,K_H);

paramsF.B = B; paramsF.q_G = q_G; paramsF.K_blocks = K_G; paramsF.q_F = q_F*ones(K_F,1); paramsF.K_facs = K_F; paramsF.l_F = l_F;

% Initialize factors as PC estimates plus random noise
T=rows(bigZ_ns{1}{1});
bigZ=bigZ_ns;
for b=1:B;
    g{b}=G_pc{b}(:,1:K_G(b)) + noise_var*randn(T,K_G(b));
    if Bsub(b) == 0,
        bigZ{b}=standard_miss(bigZ_ns{b},999.99);
        for k = 1:K_G(b),
            g{b}(:,k) = sign(corr(g{b}(:,k),bigZ{b}(:,k)))*g{b}(:,k);
        end;
    else,
        for s=1:Bsub(b);
            bigZ{b}{s}=standard_miss(bigZ_ns{b}{s},999.99);
            h{b}{s}=H_pc{b}{s}(:,1:K_H{b}(s)) + noise_var*randn(T,K_H{b}(s));     
            for k = 1:K_H{b}(s),
                h{b}{s}(:,k)=sign(corr(h{b}{s}(:,k),bigZ{b}{s}(:,k)))*h{b}{s}(:,k);
            end;
        end;
        allH = cell2mat(h{b});
        for k = 1:K_G(b),
            g{b}(:,k) = sign(corr(g{b}(:,k),allH(:,k)))*g{b}(:,k);
        end;
    end;
end;

allG = cell2mat(g);
[ehat1,f,lambda_G1,ss1]=pc_T(allG,K_F);
f = f + noise_var*randn(T,K_F);
for k = 1:K_F,
    f(:,k) = sign(corr(f(:,k),allG(:,k)))*f(:,k);
end;

% Now initialize parameters based on starting values for the factors

% First regress starting values of aggregate factors on their own lags
% to obtain starting values for the psi_F's and sig2_F's
startprior.Psi.mean = 0; startprior.Psi.var = 1e10;
startprior.Sigma.shape = 0; startprior.Sigma.dof = 0;
[psi_F,sig2_F]=draw_psi_sigma(f,.5*ones(K_F,1),ones(K_F,1),startprior,paramsF,sig_fix);
% Now regress starting values of block factors on starting values of
% aggregate factors to obtain starting value for Lambda_F
F_lags = f;
for l = 1:max(l_F),
    F_lags = [F_lags lagn(f,l)];
end;
lambda_F = ((F_lags'*F_lags)\(F_lags'*allG))';
for l = 1:max(l_F)+1,
    Lambda_F{l} = lambda_F(:,(l-1)*K_F+1:l*K_F);        
end;
for i = 1:K_F,                
    signlam(i) = sign(Lambda_F{1}(i,i));
    f(:,i) = f(:,i)*signlam(i);
    for l = 1:l_F+1,
        Lambda_F{l}(:,i) = Lambda_F{l}(:,i)*signlam(i);
    end;
end;     
% Fix upper left block of Lambda_F{1} to be lower triangular with ones on diagonal
for j = 1:K_F,
    Lambda_F{1}(j,j) = 1;
    Lambda_F{1}(j,j+1:K_F) = zeros(1,K_F-j);
end;

for b=1:B,

    if b == 1,
        Gb_indx = 1:K_G(b);
    else,
        Gb_indx = sum(K_G(1:b-1))+1:sum(K_G(1:b));        
    end;

    paramsGb.q_F = q_G(b)*ones(K_G(b),1); paramsGb.K_facs = K_G(b); paramsGb.l_F = l_G(b);

    lambda_F_mat = cell2mat(Lambda_F);
    lambda_Fb = lambda_F_mat(Gb_indx,:);
    eG{b} = compute_resids(g{b},f,lambda_Fb,paramsGb);
    [psi_G{b},sig2_G{b}]=draw_psi_sigma(eG{b},.5*ones(K_G(b),1),ones(K_G(b),1),startprior,paramsGb,sig_fix);

    if Bsub(b) == 0,

            yy=bigZ{b};
            paramsG.q_G = q_Z(b); paramsG.K_blocks = ones(size(yy,2),1); paramsG.q_F = q_G(b)*ones(K_G(b),1); paramsG.q_G = q_Z(b)*ones(size(yy,2),1); paramsG.K_facs = K_G(b); paramsG.l_F = l_G(b)*ones(size(yy,2),1); 

            % Now regress all Z_bn on G_b in order to get starting values for Lambda_Gb  
            gb_lags = g{b};
            for l = 1:l_G(b),
                gb_lags = [gb_lags lagn(g{b},l)];
            end;
            % regress individual variables on subblock factor starting values
            lambda_Gb = ((gb_lags'*gb_lags)\(gb_lags'*yy))';
            for l = 1:l_G(b)+1,
                Lambda_G{b}{l} = lambda_Gb(:,(l-1)*K_G(b)+1:l*K_G(b));        
            end;
            % Find Zbi whose loading on Gb is closest to 1 in absolute value
            [mindist, neworder] = sort(abs(Lambda_G{b}{1})-ones(size(Lambda_G{b}{1})),'descend');
            for l = 1:l_G(b)+1,                
                Lambda_G{b}{l} = switch_rows(Lambda_G{b}{l},neworder(1,:));                        
            end;
            bigZ{b} = switch_cols(bigZ{b},neworder(1,:));
            bigZ_ns{b} = switch_cols(bigZ_ns{b},neworder(1,:));
            for i = 1:K_G(b),                
                signlam(i) = sign(Lambda_G{b}{1}(i,i));
                g{b}(:,i) = g{b}(:,i)*signlam(i);
                for l = 1:l_G(b)+1,
                    Lambda_G{b}{l}(:,i) = Lambda_G{b}{l}(:,i)*signlam(i);
                end;
            end;           
            for j = 1:K_G(b),
                Lambda_G{b}{1}(j,j) = 1;
                Lambda_G{b}{1}(j,j+1:K_G(b)) = zeros(1,K_G(b)-j);
            end;  

            eZ{b} = compute_resids(yy,g{b},cell2mat(Lambda_G{b}),paramsG);
            paramsZ.q_F = q_Z(b)*ones(size(yy,2),1); paramsZ.K_facs = size(yy,2); paramsZ.l_F = l_H(b)*ones(size(yy,2),1); 
            [psi_Z{b},sig2_Z{b}]=draw_psi_sigma(eZ{b},.5*ones(size(yy,2),1),ones(size(yy,2),1),startprior,paramsZ,.5);


    else, % if Bsub(b) == 0

        Lambda_G = [];

        % Now regress all H_bs on G_b in order to get starting values for Lambda_Gb  
        g{b} = g{b};
        gb_lags = g{b};
        for l = 1:l_G(b),
            gb_lags = [gb_lags lagn(g{b},l)];
        end;
        Hb = cell2mat(h{b});
        % regress individual variables on subblock factor starting values
        lambda_Gb = ((gb_lags'*gb_lags)\(gb_lags'*Hb))';
        for l = 1:l_G(b)+1,
            Lambda_Gsub{b}{l} = lambda_Gb(:,(l-1)*K_G(b)+1:l*K_G(b));        
        end;
        for i = 1:K_G(b),                
            signlam(i) = sign(Lambda_Gsub{b}{1}(i,i));
            g{b}(:,i) = g{b}(:,i)*signlam(i);
            for l = 1:l_G(b)+1,
                Lambda_Gsub{b}{l}(:,i) = Lambda_Gsub{b}{l}(:,i)*signlam(i);
            end;
        end;           
        % Fix upper left block of Lambda_Gsub{b}{1} to be lower triangular
        % with ones on diagonal
        for j = 1:K_G(b),
            Lambda_Gsub{b}{1}(j,j) = 1;
            Lambda_Gsub{b}{1}(j,j+1:K_G(b)) = zeros(1,K_G(b)-j);
        end;   


        paramsG.B = Bsub(b); paramsG.q_G = q_H; paramsG.K_blocks = K_H{b}; paramsG.q_F = q_G; paramsG.K_facs = K_G(b); paramsG.l_F = l_G;         
        % Now find starting values for level four params
        for s = 1:Bsub(b),

            yy=bigZ{b}{s};
            if s == 1,
                Hbs_indx = 1:K_H{b}(s);
            else,
                Hbs_indx = sum(K_H{b}(1:s-1))+1:sum(K_H{b}(1:s));        
            end;            
            paramsH.B = Nsub(b,s); paramsH.q_G = q_Z; paramsH.K_blocks = ones(size(yy,2),1); paramsH.q_F = q_H(b)*ones(K_H{b}(s),1); paramsH.q_G = q_Z(b)*ones(size(yy,2),1); paramsH.K_facs = K_H{b}(s); paramsH.l_F = l_H(b)*ones(size(yy,2),1); 

            lamda_Gsub_mat = cell2mat(Lambda_Gsub{b});
            lamda_Gbs = lamda_Gsub_mat(Hbs_indx,:);
            eH{b}{s} = compute_resids(h{b}{s},g{b},lamda_Gbs,paramsG);
            [psi_H{b}{s},sig2_H{b}{s}]=draw_psi_sigma(eH{b}{s},.5*ones(K_H{b}(s),1),ones(K_H{b}(s),1),startprior,paramsH,.5);

            % Now regress all Z_bsn on H_bs in order to get starting values for Lambda_Hbs  
            hbs_lags = h{b}{s};
            for l = 1:l_H(b),
                hbs_lags = [hbs_lags lagn(h{b}{s},l)];
            end;
            % regress individual variables on subblock factor starting values
            lambda_Hbs = ((hbs_lags'*hbs_lags)\(hbs_lags'*yy))';
            for l = 1:l_H(b)+1,
                Lambda_H{b}{s}{l} = lambda_Hbs(:,(l-1)*K_H{b}(s)+1:l*K_H{b}(s));        
            end;
            % Find Zbsi whose loading on Hbs is closest to 1 in absolute value
            [mindist, neworder] = sort(abs(Lambda_H{b}{s}{1})-ones(size(Lambda_H{b}{s}{1})),'descend');
            for l = 1:l_H(b)+1,                
                Lambda_H{b}{s}{l} = switch_rows(Lambda_H{b}{s}{l},neworder(1,:));                        
            end;
            bigZ{b}{s} = switch_cols(bigZ{b}{s},neworder(1,:));
            bigZ_ns{b}{s} = switch_cols(bigZ_ns{b}{s},neworder(1,:));
            for i = 1:K_H{b}(s),                
                signlam(i) = sign(Lambda_H{b}{s}{1}(i,i));
                h{b}{s}(:,i) = h{b}{s}(:,i)*signlam(i);
                for l = 1:l_H(b)+1,
                    Lambda_H{b}{s}{l}(:,i) = Lambda_H{b}{s}{l}(:,i)*signlam(i);
                end;
            end;           

            for j = 1:K_H{b}(s),
                Lambda_H{b}{s}{1}(j,j) = 1;
                Lambda_H{b}{s}{1}(j,j+1:K_H{b}(s)) = zeros(1,K_H{b}(s)-j);
            end;  

            eZ{b}{s} = compute_resids(yy,h{b}{s},cell2mat(Lambda_H{b}{s}),paramsH);
            paramsZ.q_F = q_Z(b)*ones(size(yy,2),1); paramsZ.K_facs = size(yy,2); paramsZ.l_F = l_H(b)*ones(size(yy,2),1); 
            [psi_Z{b}{s},sig2_Z{b}{s}]=draw_psi_sigma(eZ{b}{s},.5*ones(Nsub(b,s),1),ones(Nsub(b,s),1),startprior,paramsZ,.5);

        end; % for s = 1:Bsub(b)

    end; % if Bsub(b) == 0
  
end; % for b=1:B,

save([matname,'_start']);