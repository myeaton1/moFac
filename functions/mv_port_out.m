clear, home, close all;

fileList            =   {'kalman_oos_k_2_startY_19562014_08_17_LR.mat', ...
                            'kalman_oos_k_3_startY_19562014_08_18_LR.mat', ...
                            'kalman_oos_k_4_startY_19562014_08_18_LR.mat', ...
                            'kalman_oos_k_5_startY_19562014_08_20_LR.mat', ...
                            'kalman_oos_k_6_startY_19562014_08_20_LR.mat'};
                        
fileList            =   {'kalman_oos_k_2_startY_19562014_09_24_LR.mat'};
                        
kVarTrue            =   0;                        
                        
for fileN           =   1:length(fileList)
    
    file                =   fileList{fileN};
    
    if kVarTrue == 1
        file                =   'kalman_oos_k_varFactors_startY_19562014_08_15_LR.mat';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Import data
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd 'W:\data\kalman_sort_res';

rf                  =   dlmread('rf.csv',',');
dateRf              =   rf(:,1);
yRf                 =   floor(dateRf/100);
mRf                 =   dateRf - yRf*100;
dRf                 =   ones(size(yRf));
dateRf              =   datenum(yRf,mRf,dRf);
rf                  =   rf(:,2)/100;

% load('kalman_oos_k_ffFactors_startY_19562014_06_27_LR.mat');
load(file);
fprintf(['\nk = ',num2str(k),'\n']);

cd 'W:\data\kalman_sort_res\output';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create port. weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma               =   0.5;
lambdaRidge         =   0.5;

retMarkAll          =   [];
lambda              =   [];
ztEnd               =   [];
for ye              =   1:length(startYears),
    for mon             =   1:12,
        
        if kVarTrue == 1
            k = kList{ye,mon};
        end
        
        if ~(ye             ==  length(startYears) && mon == 12),
            
            lambdaTemp          =   outputs{ye,mon}.xi_TT_out;
            lambdaTemp          =   reshape(lambdaTemp,size(lambdaTemp,1),size(lambdaTemp,3));
            
            lambdaGrid{ye,mon}  =   lambdaTemp(1:k,end);
            
            zt{ye,mon}          =   lambdaTemp(k+1:2*k,:) + lambdaTemp(2*k+1:end,:);
            lambdatm1{ye,mon}   =   lambdaTemp(k+1:2*k,:);
            ft{ye,mon}          =   lambdaTemp(2*k+1:end,:);
            
            Sige{ye,mon}        =   diag(nanvar(retData{ye,mon} - (B{ye,mon}*zt{ye,mon})'));

            Sigt{ye,mon}        =   B{ye,mon}*outputs{ye,mon}.Sigf*B{ye,mon}' + Sige{ye,mon} + ...
                                        lambdaRidge*eye(size(Sige{ye,mon},1));
            
            
            condSigt{ye,mon}    =   cond(Sigt{ye,mon});
            
            mut{ye,mon}         =   B{ye,mon}*lambdaGrid{ye,mon};
            
            wt{ye,mon}          =   (1/gamma)*(Sigt{ye,mon}\mut{ye,mon});
            
            onesN               =   ones(size(Sigt{ye,mon},1),1);
            SigtC               =   Sigt{ye,mon};
            
            wt2{ye,mon}         = (SigtC\onesN)/(onesN'*(SigtC\onesN)) ...
                                    + (1/gamma) * (inv(SigtC) - ...
                                    (SigtC\onesN*onesN'/SigtC)/(onesN'*(SigtC\onesN))) * mut{ye,mon};
        end,
        
    end,
end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Create return s.d.s, E[r]
%   Sort on them
%%%%%%%%%%%%%%%%%%%%%%%%%%%
sortType = {'er','erSD2'}; % 'erSD'
for sortNum = 1:length(sortType),
    s = sortType{sortNum};
    retLS.(s) = [];
    retQall.(s) = [];
end,

retMark = [];
retEq = [];
avgWts = [];
eqWts = [];
retM = [];
retMq = [];
mOpt = 1;
labT = cell(0);

for ye = 1:length(startYears)
    for mon = 1:12
        
        if kVarTrue == 1
            k = kList{ye,mon};
        end

        %%%%%%%%%%%%%%%%%%%%%
        % Build data
        %%%%%%%%%%%%%%%%%%%%%
        rSD{ye,mon} = nan(size(B{ye,mon},1),1);
        for i = 1:size(B{ye,mon},1),
            bi = B{ye,mon}(i,:)';
%             rSD{ye,mon}(i) = sqrt(bi'*outputs{ye,mon}.varZ*bi);
        end,
        
%         varZ = [varZ ; diag(outputs{ye,mon}.varZ)'];
        
        Bt = B{ye,mon};
        ret = retData{ye,mon}';
        
        rSD2{ye,mon} = std(ret,0,2);

        lambda = outputs{ye,mon}.xi_TT_out;
        lambda = reshape(lambda,size(lambda,1),size(lambda,3));
        lambda = lambda(1:k,end)';
        
        if ye == 29 && mon == 2,
%             rSD{ye,mon} = rSD{ye,mon}(1:59);
            Bt = Bt(1:59,:);
            rSD2{ye,mon} = std(ret(1:59,:),0,2);
        end,
        
        er.er{ye,mon} = Bt*lambda';
        
%         er.erSD{ye,mon} = Bt*lambda'./rSD{ye,mon};
        
        er.erSD2{ye,mon} = Bt*lambda'./rSD2{ye,mon};
        
        %%%%%%%%%%%%%%%%%%%%%
        % Sorts
        %%%%%%%%%%%%%%%%%%%%%
        
        if ~(ye == length(startYears) && mon == 12)
            
            if mon == 12,
                yeF = ye+1;
                monF = 1;
            else
                yeF = ye;
                monF = mon+1;
            end,
            
            %%%%%%%%%%%%%%
            % Factors - Markowitz
            %%%%%%%%%%%%%%
            firmIdx1 = ismember(pList{ye,mon},pList{yeF,monF});
            firmIdx2 = ismember(pList{yeF,monF},pList{ye,mon});
            
            retF = retData{yeF,monF}(end,firmIdx2);
            retC = retData{ye,mon}(end,firmIdx1);
            wtC = wt2{ye,mon}(firmIdx1);
            
            retTemp = nansum(retF'.*wtC);
            
            retMark = [retMark ; retTemp];
            
            avgWts = [avgWts ; (nanmean(abs(wtC)))*100];
            eqWts = [eqWts ; (1/length(wtC))*100];
            
            %%%%%%%%%%%%%%
            % Equal Weighted
            %%%%%%%%%%%%%%
            retEq = [retEq; nanmean(retData{yeF,monF}(end,:))];
            
            %%%%%%%%%%%%%%
            % Factors - L/S
            %%%%%%%%%%%%%%
            for sortNum = 1:length(sortType),
                s = sortType{sortNum};

                retL = er.(s){ye,mon};
                retF = retData{yeF,monF}(end,:);

                moQuan = .3:.1:.7;
                quan = [ min(retL), quantile(retL,moQuan),max(retL) ];

                for q = 1:(size(quan,2)-1)
                    idx = (retL >= quan(q)) & (retL <= quan(q+1));
                    pListtM = pList{ye,mon}(idx);
                    pListtM = ismember(pList{yeF,monF},pListtM);
                    retTemp(q) = nanmean(retF(pListtM));
                end,
                retLS.(s) = [retLS.(s) ; retTemp(end) - retTemp(1)];
                retQall.(s) = [retQall.(s) ; retTemp];
            end,

            %%%%%%%%%%%%%%
            % Mo
            %%%%%%%%%%%%%%
            retL = prod(1+retData{ye,mon}(end-11:end-mOpt,:))-1;
            retF = retData{yeF,monF}(end,:);

            moQuan = .3:.1:.7;
            quan = [ min(retL), quantile(retL,moQuan),max(retL) ];

            for q = 1:(size(quan,2)-1)
                idx = (retL >= quan(q)) & (retL <= quan(q+1));
                pListtM = pList{ye,mon}(idx);
                pListtM = ismember(pList{yeF,monF},pListtM);
                retMtemp(q) = nanmean(retF(pListtM));
            end,
            retM = [retM ; retMtemp(end) - retMtemp(1)];
            retMq = [retMq ; retMtemp];   
            
            % Labels
            [~, tempLab] = strtok(lab{yeF,monF},'-');
            labT = [labT ; {tempLab(2:end)}];
        
        end
    end
end

% for sortNum = 1:length(sortType),
%     s = sortType{sortNum};
%     retLS.([s,'c']) = cumprod(1+retLS.(s))-1;
% end,
% cumRetM = cumprod(1+retM)-1;
% 
% cumRetMark = cumprod(1+retMark)-1;
% cumRetEq = cumprod(1+retEq)-1;
% 
% % Plot of M-V vs. E-W cret for each
% plot([cumRetMark,cumRetEq])
% dist = 4;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% legend(['M-V Efficient, \lambda = ',num2str(lambdaRidge)],...
%     'Equal Weighted','Location','NorthWest');
% title('Cum. Returns, Various Factor Ports.')
% % saveas(gcf,['cumretMV_k_',num2str(k),'_lambda_',num2str(lambdaRidge*100)], 'epsc');
% saveas(gcf,['cumretMV_k_Var_lambda_',num2str(lambdaRidge*100)], 'epsc');
% 
% plot([retLS.erc,cumRetM]) %retLS.erSDc
% dist = 4;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% legend('E[r_i] (L/S) factor', 'Momentum Portfolio (prior 2-12m)',...
%     'Location','NorthWest'); %'E[r_i]/std[r_i] factor'
% title('Cum. Returns, Various Factor Ports. (L/S EW T30%/B30%)')
% % saveas(gcf,['cumretLS_k_',num2str(k),'_lambda_',num2str(lambdaRidge*100)], 'epsc');
% saveas(gcf,['cumretLS_k_Var_lambda_',num2str(lambdaRidge*100)], 'epsc');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   Track # of PCs chosen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% kL      =   kList';
% kL      =   kL(:);
% kL      =   cell2mat(kL);
% 
% ICpL    =   ICp';
% ICpL    =   ICpL(:);
% ICpL    =   cell2mat(ICpL);
% 
% tstatL  =   tstat';
% tstatL  =   tstatL(:);
% tstatL  =   cell2mat(tstatL);
% 
% plot(kL)
% dist = 4;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% title('Optimal number of PCs as chosen by IC_p')
% saveas(gcf,'optimal_pcs', 'pdf');
% 
% plot(ICpL)
% dist = 4;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% legend('1 PC','2 PCs','3 PCs','4 PCs','5 PCs','6 PCs','Location','NorthWest'); %'E[r_i]/std[r_i] factor'
% title('IC_p for various numbers of PCs')
% saveas(gcf,'icp', 'pdf');
% 
% plot(tstatL)
% dist = 4;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% legend('1 PC','2 PCs','3 PCs','4 PCs','5 PCs','6 PCs','Location','NorthEast');
% title('t-stats. of \alpha of Momentum portfolio on z_t using Newey-West S.E.s')
% saveas(gcf,'alpha_tstats', 'pdf');

% append_pdfs('var_pcs_plots.pdf','optimal_pcs.pdf','icp.pdf','alpha_tstats.pdf',...
%     'cumret_ls_ew_VarFac_comp_lambda_150.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Compare these sorts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rf              =   rf(dateRf >= datenum('1/1/1966') & dateRf <= datenum('11/1/2012'));

sortType        =   {'mark','ew','er','mo'};
retAll.mark     =   retMark - rf;
retAll.ew       =   retEq - rf;
retAll.er       =   retLS.er;
retAll.mo       =   retM;

if kVarTrue == 0
    retAllK.(['k',num2str(k)]) ...
                    =   retAll;
end

for sortNum = 1:length(sortType),
    s = sortType{sortNum};
    retAll.([s,'sum']) = [mean(retAll.(s)) , std(retAll.(s)) , mean(retAll.(s))./std(retAll.(s)), ...
                            (mean(retAll.(s))./std(retAll.(s)))*sqrt(12)];
end,

outSum = [{'Sort','Mean','SD','Sharpe','Sharpe (ann)'};...
    [{'M-V Efficient'}, num2cell(retAll.marksum)];...
    [{'E-W'}, num2cell(retAll.ewsum)]; ...
    [{'E[r] (L/S) factor'}, num2cell(retAll.ersum)];...
    [{'Momentum'}, num2cell(retAll.mosum)]];

if kVarTrue == 1
    xlswrite(['ret_sum_varK_fac_comp',num2str(lambdaRidge*100),'.xlsx'],outSum)
else
    xlswrite(['ret_sum_k_',num2str(k),'_fac_comp',num2str(lambdaRidge*100),'.xlsx'],outSum);
end


end

if kVarTrue == 1
    load(['portRetAll',num2str(lambdaRidge*100),'.mat'])
    retAllK.kVar    =   retAll;
end

save(['portRetAll',num2str(lambdaRidge*100),'.mat'],'retAllK');











