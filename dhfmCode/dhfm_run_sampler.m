   
jj=1;
set_priors;
set_params;

randn('state',999);  
for gibbs=1:n_gibbs,
    
    % Update global factors and corresponding parameters using last iteration's draws of params and block factors 
    [alphaF,lambda_F_tilde]=cal_alpha(f,cell2mat(psi_G'),Lambda_F);
    f = sample_facs(g,zeros(size(f)),lambda_F_tilde,psi_F,sig2_F,psi_G,sig2_G,paramsF);    
    Lambda_F=draw_lambda(g,f,psi_G,sig2_G,prior.F.Lambda,paramsF);
    Lambda_F_mat = cell2mat(Lambda_F);

    [alphaF,lambda_F_tilde]=cal_alpha(f,cell2mat(psi_G'),Lambda_F);
    [psi_F,sig2_F]=draw_psi_sigma(f,psi_F,sig2_F,prior.F,paramsF,sig_fix);
    
    for b=1:B;
       
        if b == 1,
            Gb_indx = 1:K_G(b);
        else,
            Gb_indx = sum(K_G(1:b-1))+1:sum(K_G(1:b));        
        end;
        Lambda_Fb = Lambda_F_mat(Gb_indx,:);
        
        % update the blocks without sub-blocks  
        if Bsub(b) == 0,
            
            yy=bigZ{b};
            yy_cell = mat2cell(yy,size(yy,1),repmat(1,1,size(yy,2)));
            psi_Z_cell = mat2cell(psi_Z{b},repmat(1,size(psi_Z{b},1),1),size(psi_Z{b},2));
            sig2_Z_cell = mat2cell(sig2_Z{b},repmat(1,size(sig2_Z{b},1),1),size(sig2_Z{b},2));
            
            paramsG.q_G = q_Z(b); paramsG.K_blocks = ones(size(yy,2),1); paramsG.q_F = q_G(b)*ones(K_G(b),1); paramsG.q_G = q_Z(b)*ones(size(yy,2),1); paramsG.K_facs = K_G(b); paramsG.l_F = l_G(b)*ones(size(yy,2),1); 

            [alphaG{b},lambda_Gb_tilde]=cal_alpha(g{b},psi_Z{b},Lambda_G{b});  

            paramsG.B = size(yy,2); paramsG.q_G = q_Z(b)*ones(size(yy,2),1); paramsG.K_blocks = ones(size(yy,2),1); paramsG.q_F = q_G(b); paramsG.K_facs = K_G(b); paramsG.l_F = l_G(b)*ones(size(yy,2),1); 
            g{b} = sample_facs(yy_cell,alphaF{b},lambda_Gb_tilde,psi_G{b},sig2_G{b},psi_Z_cell,sig2_Z_cell,paramsG);
            Lambda_G{b}=draw_lambda(yy_cell,g{b},psi_Z_cell,sig2_Z_cell,prior.G.Lambda,paramsG);
            eG{b} = compute_resids(g{b},f,Lambda_Fb,paramsF);
            
            paramsGb.q_F = q_G(b)*ones(K_G(b),1); paramsGb.K_facs = K_G(b); paramsGb.l_F = l_G(b);
            [psi_G{b},sig2_G{b}]=draw_psi_sigma(eG{b},psi_G{b},sig2_G{b},prior.G,paramsGb,sig_fix);
            
            eZ{b} = compute_resids(yy,g{b},cell2mat(Lambda_G{b}),paramsG);
            paramsZ.q_F = q_Z(b)*ones(size(yy,2),1); paramsZ.K_facs = size(yy,2); paramsZ.l_F = l_H(b)*ones(size(yy,2),1); 
            
            [psi_Z{b},sig2_Z{b}]=draw_psi_sigma(eZ{b},psi_Z{b},sig2_Z{b},prior.Z,paramsZ,[]);
            
        % update the blocks with subblocks    
        else, % if Bsub(b) == 0
            
            [alphaG{b},lambda_Gb_tilde]=cal_alpha(g{b},cell2mat(psi_H{b}'),Lambda_Gsub{b});            
            paramsGb.B = Bsub(b); paramsGb.q_G = q_H(b)*ones(sum(K_H{b}),1); paramsGb.K_blocks = K_H{b}; paramsGb.q_F = q_G(b)*ones(K_G(b),1); paramsGb.K_facs = K_G(b); paramsGb.l_F = l_G(b)*ones(sum(K_H{b}),1); 
            g{b} = sample_facs(h{b},alphaF{b},lambda_Gb_tilde,psi_G{b},sig2_G{b},psi_H{b},sig2_H{b},paramsGb);
            Lambda_Gsub{b}=draw_lambda(h{b},g{b},psi_H{b},sig2_H{b},prior.G.Lambda,paramsGb);            

            eG{b} = compute_resids(g{b},f,Lambda_Fb,paramsF);
            [psi_G{b},sig2_G{b}]=draw_psi_sigma(eG{b},psi_G{b},sig2_G{b},prior.G,paramsGb,sig_fix);

            % Update sub-blocks
            for s=1:Bsub(b);           

                yy=bigZ{b}{s};        
                if s == 1,
                    Hbs_indx = 1:K_H{b}(s);
                else,
                    Hbs_indx = sum(K_H{b}(1:s-1))+1:sum(K_H{b}(1:s));        
                end;

                [alphaHb,lambda_Hbs_tilde]=cal_alpha(h{b}{s},psi_Z{b}{s},Lambda_H{b}{s});    

                paramsH.B = Nsub(b,s); paramsH.q_G = q_Z; paramsH.K_blocks = ones(size(yy,2),1); paramsH.q_F = q_H(b)*ones(K_H{b}(s),1); paramsH.q_G = q_Z(b)*ones(size(yy,2),1); paramsH.K_facs = K_H{b}(s); paramsH.l_F = l_H(b)*ones(size(yy,2),1); 
                yy_cell = mat2cell(yy,size(yy,1),repmat(1,1,size(yy,2)));
                psi_Z_cell = mat2cell(psi_Z{b}{s},repmat(1,size(psi_Z{b}{s},1),1),size(psi_Z{b}{s},2));
                sig2_Z_cell = mat2cell(sig2_Z{b}{s},repmat(1,size(sig2_Z{b}{s},1),1),size(sig2_Z{b}{s},2));

                h{b}{s}=sample_facs(yy_cell,alphaG{b}{s},lambda_Hbs_tilde,psi_H{b}{s},sig2_H{b}{s},psi_Z_cell,sig2_Z_cell,paramsH);
                Lambda_H{b}{s} = draw_lambda(yy_cell,h{b}{s},psi_Z_cell,sig2_Z_cell,prior.H.Lambda,paramsH);

                Lambda_Gsub_mat = cell2mat(Lambda_Gsub{b});
                Lambda_Gbs = Lambda_Gsub_mat(Hbs_indx,:);
                paramsHbs.K_facs = K_G(b); paramsHbs.l_F = l_G(b); 
                eH{b}{s} = compute_resids(h{b}{s},g{b},Lambda_Gbs,paramsHbs);
                [psi_H{b}{s},sig2_H{b}{s}]=draw_psi_sigma(eH{b}{s},psi_H{b}{s},sig2_H{b}{s},prior.H,paramsH,sig_fix);

                eZ{b}{s} = compute_resids(yy,h{b}{s},cell2mat(Lambda_H{b}{s}),paramsH);
                paramsZ.q_F = q_Z(b)*ones(size(yy,2),1); paramsZ.K_facs = size(yy,2); paramsZ.l_F = l_H(b)*ones(size(yy,2),1); 
                [psi_Z{b}{s},sig2_Z{b}{s}]=draw_psi_sigma(eZ{b}{s},psi_Z{b}{s},sig2_Z{b}{s},prior.Z,paramsZ,[]);
   
            end; % for s=1:Bsub(b)         
            
        end; % if Bsub(b) == 0

    end; % for b=1:B2
 
    if gibbs > n_burn & mod(gibbs,n_skip)==0,    
       shares = vardecomp(Lambda_F,Lambda_Gsub,Lambda_G,Lambda_H,psi_F,psi_G,psi_H,psi_Z,sig2_F,sig2_G,sig2_H,sig2_Z,params);
       F{jj}=f;
       G{jj}=g;
       H{jj}=h;
       EG{jj} = eG;
       EH{jj} = eH;
       LAMBDA_F{jj}=Lambda_F;
       LAMBDA_Gsub{jj}=Lambda_Gsub;
       LAMBDA_G{jj}=Lambda_G;
       LAMBDA_H{jj}=Lambda_H;
       PSI_F{jj}= psi_F;
       PSI_G{jj}= psi_G;
       PSI_H{jj}= psi_H;
       PSI_Z{jj}= psi_Z;
       SIG2_F{jj}= sig2_F;
       SIG2_G{jj}=sig2_G;       
       SIG2_H{jj}=sig2_H;
       SIG2_Z{jj}= sig2_Z;
       SHARES{jj} = shares;
       jj=jj+1;
   end;  % end storing results

   if mod(gibbs,50) ==0;
        disp(sprintf('Completed %d Draws',gibbs));
        print_vec = psi_F';
        if ~isempty(psi_G{1}),
            print_vec = [print_vec psi_G{1}'];
        end;
        if ~isempty(psi_H{1}{1}),
            print_vec = [print_vec psi_H{1}{1}'];
        end;
        mymprint(print_vec);
        save(matname);
   end;

end; % end gibbs

L4_alldraws.LAMBDA_F = LAMBDA_F;
L4_alldraws.LAMBDA_Gsub = LAMBDA_Gsub;
L4_alldraws.LAMBDA_G = LAMBDA_G;
L4_alldraws.LAMBDA_H = LAMBDA_H;
L4_alldraws.SIG2_F = SIG2_F;
L4_alldraws.SIG2_G = SIG2_G;
L4_alldraws.SIG2_H = SIG2_H;
L4_alldraws.SIG2_Z = SIG2_Z;
L4_alldraws.PSI_F = PSI_F;
L4_alldraws.PSI_G = PSI_G;
L4_alldraws.PSI_H = PSI_H;
L4_alldraws.PSI_Z = PSI_Z;
L4_alldraws.F = F;
L4_alldraws.G = G;
L4_alldraws.H = H;
L4_alldraws.SHARES = SHARES;
L4_alldraws.EG = EG;
L4_alldraws.EH = EH;
L4_alldraws.ndraws = jj-1;
L4_alldraws.params = params;
L4_alldraws.matname = matname;

save([matname,'_alldraws'],'L4_alldraws');
