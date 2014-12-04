clear, home, close all;

%% Set relevant paths, initialize key variables
if isunix,
    resultspath = '/home/rcemry01/data/kalman_sort_res/';
    outputpath = '/home/rcemry01/data/output';
    codepath = '/home/rcemry01/code/factor';
else
    resultspath = 'Y:\data\kalman_sort_res';
    outputpath = 'Y:\data\kalman_sort_res\output';
    codepath = 'S:\MattY\emanuel\Momentum\Code';
end,

addpath(codepath);

k = 5;				% # PCAs in Kalman filter
win = 10;			% length (in years) of EM PCA and Kalman filter windows
start = 1955;
startYears = start:(2012-win);

outfile = [resultspath, 'kalman_oos_k_', num2str(k), '_startY_',num2str(start),...
    datestr(today,'yyyy_mm_dd'), '_LR.mat'];

magRandPhi  = 10e-2;
magRandSigw = 10e-4;
magRandSigf = 10e-2;
magRand     = [repmat(magRandPhi,k,1);repmat(magRandSigw,k,1);...
    repmat(magRandSigf,k,1)];

nRuns = 5;



for ye = 1:length(startYears),
	for mon = 1:12,

		year = startYears(ye);

		[f{ye,mon}, B{ye,mon}, retData{ye,mon}, pList{ye,mon}, lab{ye,mon}] = ...
			em_pca_it(k, win ,year, mon);
        
        if ye == 1 && mon == 1,
            % Initialize parameters 
            Phi = .7*eye(k);
            Sigw = .001*ones(k,1);
            Sigf = var(f{1,1})';
            params0 = [diag(Phi); Sigw; Sigf];
            param_in = params0;
            nDim = length(param_in);
        end,

        for i = 1:nRuns,
            op.xi_TT_out = [];
            if i > 1,
                while isempty(op.xi_TT_out),
                    n = randn(nDim,1).*magRand;
                    test{i} = param_in+n;
                    [~, op] = loglik_kalman_filter_model4(test{i},f{ye,mon},1);

                end,
            else,
                test{i} = param_in;
            end,
            
            [testOutputs{i}, testParam_out{i}, testExitflag{i},testLIK{i}] = ...
                kalman_it(k, win ,year, mon, f{ye,mon}, B{ye,mon}, test{i}, lab{ye,mon});
            
        end,
        
        [~,argmin] = min(cell2mat(testLIK));
        
        outputs{ye,mon}     = testOutputs{argmin};
        param_out{ye,mon}   = testParam_out{argmin};
        exitflag{ye,mon}    = testExitflag{argmin};
        LIK{ye,mon}         = testLIK{argmin};

        param_in = testParam_out{argmin};

		[port{ye,mon}, portIdx{ye,mon}] = ...
			port_sort_it(k, win ,year, mon, B{ye,mon}, outputs{ye,mon}, ...
			retData{ye,mon}, pList{ye,mon}, lab{ye,mon});
        

	end,
end,


save(outfile);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   Create return s.d.s, E[r]
% %   Sort on them
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% sortType = {'er','erSD','erSD2'};
% for sortNum = 1:length(sortType),
%     s = sortType{sortNum};
%     retLS.(s) = [];
%     retQall.(s) = [];
% end,
% 
% retM = [];
% retMq = [];
% labT = cell(0);
% mOpt = 1;
% varZ = [];
% 
% for ye = 1:length(startYears),
% 	for mon = 1:12,
%         
%         %%%%%%%%%%%%%%%%%%%%%
%         % Build data
%         %%%%%%%%%%%%%%%%%%%%%
%         rSD{ye,mon} = nan(size(B{ye,mon},1),1);
%         for i = 1:size(B{ye,mon},1),
%             bi = B{ye,mon}(i,:)';
%             rSD{ye,mon}(i) = sqrt(bi'*outputs{ye,mon}.varZ*bi);
%         end,
%         
%         varZ = [varZ ; diag(outputs{ye,mon}.varZ)'];
%         
%         Bt = B{ye,mon};
%         ret = retData{ye,mon}';
%         
%         rSD2{ye,mon} = std(ret,0,2);
% 
%         lambda = outputs{ye,mon}.xi_TT_out;
%         lambda = reshape(lambda,size(lambda,1),size(lambda,3));
%         lambda = lambda(1:k,end)';
%         
%         if ye == 29 && mon == 2,
%             rSD{ye,mon} = rSD{ye,mon}(1:59);
%             Bt = Bt(1:59,:);
%         end,
%         
%         er.er{ye,mon} = Bt*lambda';
%         
%         er.erSD{ye,mon} = Bt*lambda'./rSD{ye,mon};
%         
%         er.erSD2{ye,mon} = Bt*lambda'./rSD2{ye,mon};
%         
%         %%%%%%%%%%%%%%%%%%%%%
%         % Sorts
%         %%%%%%%%%%%%%%%%%%%%%
%         
%         if ~(ye == length(startYears) && mon == 12),
%             
%             if mon == 12,
%                 yeF = ye+1;
%                 monF = 1;
%             else
%                 yeF = ye;
%                 monF = mon+1;
%             end,
% 
%             %%%%%%%%%%%%%%
%             % Factors
%             %%%%%%%%%%%%%%
%             for sortNum = 1:length(sortType),
%                 s = sortType{sortNum};
% 
%                 retL = er.(s){ye,mon};
%                 retF = retData{yeF,monF}(end,:);
% 
%                 moQuan = .1:.1:.9;
%                 quan = [ min(retL), quantile(retL,moQuan),max(retL) ];
% 
%                 for q = 1:(size(quan,2)-1)
%                     idx = (retL >= quan(q)) & (retL <= quan(q+1));
%                     pListtM = pList{ye,mon}(idx);
%                     pListtM = ismember(pList{yeF,monF},pListtM);
%                     retTemp(q) = nanmean(retF(pListtM));
%                 end,
%                 retLS.(s) = [retLS.(s) ; retTemp(end) - retTemp(1)];
%                 retQall.(s) = [retQall.(s) ; retTemp];
%             end,
% 
%             %%%%%%%%%%%%%%
%             % Mo
%             %%%%%%%%%%%%%%
%             retL = prod(1+retData{ye,mon}(end-11:end-mOpt,:))-1;
%             retF = retData{yeF,monF}(end,:);
% 
%             moQuan = .1:.1:.9;
%             quan = [ min(retL), quantile(retL,moQuan),max(retL) ];
% 
%             for q = 1:(size(quan,2)-1)
%                 idx = (retL >= quan(q)) & (retL <= quan(q+1));
%                 pListtM = pList{ye,mon}(idx);
%                 pListtM = ismember(pList{yeF,monF},pListtM);
%                 retMtemp(q) = nanmean(retF(pListtM));
%             end,
%             retM = [retM ; retMtemp(end) - retMtemp(1)];
%             retMq = [retMq ; retMtemp];
% 
%             [~, tempLab] = strtok(lab{yeF,monF},'-');
%             labT = [labT ; {tempLab(2:end)}];
%         
%         end,
%     end,
% end,
% 
% for sortNum = 1:length(sortType),
%     s = sortType{sortNum};
%     retLS.([s,'c']) = cumprod(1+retLS.(s))-1;
% end,
% cumRetM = cumprod(1+retM)-1;
% 
% % Plot of L/S cret for each
% plot([retLS.erc,retLS.erSDc,retLS.erSD2c])
% dist = 3;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% legend('E[r_i] factor', 'E[r_i]/std[r_i] factor','E[r_i]/std[r_i] (2) factor','Location','NorthWest');
% title('Cum. Returns, Various Factor Ports. (L/S EW T10%/B10%)')
% saveas(gcf,['cumret_ls_ew_', num2str(k),'fac_comp'], 'pdf');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% % Compare these sorts
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for sortNum = 1:length(sortType),
%     s = sortType{sortNum};
%     retQall.([s,'sum']) = [mean(retQall.(s)); std(retQall.(s)) ; mean(retQall.(s))./std(retQall.(s))];
% end,
% 
% outSum = [{'Sort','Var','Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10'};...
%     [[{'E[r] factor','Mean'}; {'','SD'}; {'','Sharpe'}], num2cell(retQall.ersum)];...
%     [[{'E[r]/sd[r] factor','Mean'}; {'','SD'}; {'','Sharpe'}], num2cell(retQall.erSDsum)];...
%     [[{'E[r]/sd[r] (2) factor','Mean'}; {'','SD'}; {'','Sharpe'}], num2cell(retQall.erSD2sum)]];
% 
% xlswrite(['ret_sum_k_',num2str(k),'_fac_comp.xlsx'],outSum)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% % Cross-sectional dispersion
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% csSD = [];
% for ye = 1:length(startYears),
% 	for mon = 1:12,
%         csSD = [csSD ; std(B{ye,mon})];
%     end,
% end,
% 
% plot(csSD)
% dist = 3;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% title('Cross-sectional Dispersion of Factors')
% saveas(gcf,['cross_sec_disp_', num2str(k)], 'pdf');
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   Form end-or-period series
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambda = [];
% Phi = [];
% sigw = [];
% sigf = [];
% portProb.h = [];
% portProb.l = [];
% for ye = 1:length(startYears),
% 	for mon = 1:12,
%         
%         if ~(ye == length(startYears) && mon == 12),
%             
%             lambdaTemp = outputs{ye,mon}.xi_TT_out;
%             lambdaTemp = reshape(lambdaTemp,size(lambdaTemp,1),size(lambdaTemp,3));
%             lambda = [lambda ; lambdaTemp(1:k,end)'];
%             
%             lambdaGrid{ye,mon} = lambdaTemp(1:k,end);
% 
%             Phi = [Phi ; diag(outputs{ye,mon}.Phi)'];
%             
%             sigw = [sigw ; diag(outputs{ye,mon}.Sigw)'];
%             sigf = [sigf ; diag(outputs{ye,mon}.Sigf)'];
% 
%             if mon == 12,
%                 yeF = ye+1;
%                 monF = 1;
%             else
%                 yeF = ye;
%                 monF = mon+1;
%             end,
% 
%             portProb.h = [portProb.h ; ...
%                 length(intersect(port{ye,mon}{1,end},port{yeF,monF}{1,end}))...
%                 /length(port{ye,mon}{1,end})];
% 
%             portProb.l = [portProb.l ; ...
%                 length(intersect(port{ye,mon}{1,1},port{yeF,monF}{1,1}))...
%                 /length(port{ye,mon}{1,1})];
% 
%             if ye == 1 && mon == 1,
%                 pInt = pList{ye,mon};
%             else
%                 pInt = intersect(pInt,pList{ye,mon});
%             end,
%         
%         end,
%         
%     end,
% end,
% 
% % % for i = 1:length(pInt),
% % %     BSer{i} = [];
% % %     for ye = 1:length(startYears),
% % %         for mon = 1:12,
% % %             
% % %             BSer{i} = [BSer{i}, B{ye,mon}(pInt(i) == pList{ye,mon},1:k)];
% % %                 
% % %         end,
% % %     end,
% % % end,
% % % 
% % % dist = 5;
% % % labT = lab';
% % % for i = 1:6,
% % %     subplot(3,2,i)
% % %     plot(BSer{1,i}(1:end-1)'.*lambda)
% % %     set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% % %     set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% % %     title(['B (factor loadings), permno = ', num2str(pInt(i))])
% % % end,
% % % cd(outputpath)
% % % saveas(gcf,'b_ser', 'pdf');
% % % 
% plot(lambda)
% title(['\lambda,k = ', num2str(k),])
% saveas(gcf,['lambda_', num2str(k)], 'pdf');
% 
% plot(Phi)
% title(['\Phi,k = ', num2str(k),])
% saveas(gcf,['Phi_', num2str(k)], 'pdf');
% % % hist([portProb.h, portProb.l])
% % % 
% % maxN = [];
% % minN = [];
% % for ye = 1:length(startYears),
% %     for mon = 1:12,
% %         
% %         if ~(ye == length(startYears) && mon == 12),
% %         er{ye,mon} = B{ye,mon}(:,1:k)*lambdaGrid{ye,mon};
% %         
% %         tb10 = sort(er{ye,mon});
% %         
% %         minN = [minN ; mean(tb10(1:10))];
% %         maxN = [maxN ; mean(tb10(end-9:end))];
% %         end,
% %     end,
% % end,
% % % 
% % % a = sort(lambda);
% % % 
% % plot([maxN, minN])
% % % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % for ye = 1:length(startYears),
% % %     for mon = 1:12,
% % %         zCheck(ye,mon) = abs(sum(outputs{ye,mon}.xi_TT_out(2:end,:,end)) - ...
% % %             f{ye,mon}(end,1)) < 1e-10;
% % %     end,
% % % end,
% % 
% 
% % Generate portfolio returns
% cd(outputpath)
% 
% ret.h = [];
% ret.l = [];
% retM = [];
% retMq = [];
% retq = [];
% retC = [];
% mo.h = [];
% mo.l = [];
% retMkt = [];
% labT = cell(0);
% mOpt = 1;
% for ye = 1:length(startYears),
%     for mon = 1:12,
%         
%         if ~(ye == length(startYears) && mon == 12),
%         if mon == 12,
%             yeF = ye+1;
%             monF = 1;
%         else
%             yeF = ye;
%             monF = mon+1;
%         end,
%         
%         % Momentum Sort
%         
%         retL = prod(1+retData{ye,mon}(end-11:end-mOpt,:))-1;
%         retF = retData{yeF,monF}(end,:);
%         
%         moQuan = .1:.1:.9;
%         quan = [ min(retL), quantile(retL,moQuan),max(retL) ];
% 
%         for q = 1:(size(quan,2)-1)
%             idx = (retL >= quan(q)) & (retL <= quan(q+1));
%             pListtM = pList{ye,mon}(idx);
%             pListtM = ismember(pList{yeF,monF},pListtM);
%             retMtemp(q) = nanmean(retF(pListtM));
%         end,
%         retM = [retM ; retMtemp(end) - retMtemp(1)];
%         retMq = [retMq ; retMtemp];
%         %           %
%         
%         for q = 1:10,
%             pListtQ = pList{ye,mon}(portIdx{ye,mon}{1,q});
%             pListtQ = ismember(pList{yeF,monF},pListtQ);
%             retQtemp(q) = nanmean(retF(pListtQ));
%         end,
%         retq = [retq ; retQtemp];
%         
%         pListt.h = pList{ye,mon}(portIdx{ye,mon}{1,end});
%         pListt.l = pList{ye,mon}(portIdx{ye,mon}{1,1});
%         
%         pListt.h = ismember(pList{yeF,monF},pListt.h);
%         pListt.l = ismember(pList{yeF,monF},pListt.l);
%         
%         ret.h = [ ret.h ; nanmean(retData{yeF,monF}(end,pListt.h)) ];
%         ret.l = [ ret.l ; nanmean(retData{yeF,monF}(end,pListt.l)) ];
%         
%         [~, tempLab] = strtok(lab{yeF,monF},'-');
%         labT = [labT ; {tempLab(2:end)}];
%         
%         retMkt = [retMkt; nanmean(retData{yeF,monF}(end,:))];
%         end,
%     end,
% end,
% 
% retF = ret.h - ret.l;
% cumRet = cumprod(1+retF)-1;
% cumRetM = cumprod(1+retM)-1;
% cumRetMkt = cumprod(1+retMkt)-1;
% 
% plot([cumRet,cumRetM])
% dist = 3;
% set(gca,'XTick',1:(12*(2^dist)):numel(labT)) 
% set(gca,'XTickLabel',labT(1:(12*(2^dist)):end));
% legend(['Factor Portfolio, k = ',num2str(k)],['Momentum Portfolio (prior ',num2str(mOpt+1),'-12m)'],'Location','NorthWest');
% title('Cum. Returns, Factor & Momentum Ports. (L/S EW T10%/B10%)')
% saveas(gcf,['cumret_ls_ew_', num2str(k),'moLag',num2str(mOpt+1)], 'pdf');
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% % Compare this mo. with FF mo.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd Y:\data
% 
% ffMo = dlmread('ffMoMon.csv',',');
% retMff = ffMo(1:end-1,2)/100;
% cumRetMff = cumprod(1+retMff)-1;
% plot([cumRetM,cumRetMff])
% 
% corr(retM,retMff)
% 
% ffMoPort = dlmread('ffMoPorts.csv',',');
% retMffq = ffMoPort(1:end-1,2:end)/100;
% 
% sumM = [mean(retMq); std(retMq) ; mean(retMq)./std(retMq)];
% sumQ = [mean(retq); std(retq) ; mean(retq)./std(retq)];
% sumMff = [mean(retMffq); std(retMffq) ; mean(retMffq)./std(retMffq)];
% 
% outSum = [{'Sort','Var','Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10'};...
%     [[{'Factor','Mean'}; {'','SD'}; {'','Sharpe'}], num2cell(sumQ)];...
%     [[{'Momentum','Mean'}; {'','SD'}; {'','Sharpe'}], num2cell(sumM)];...
%     [[{'FF Momentum','Mean'}; {'','SD'}; {'','Sharpe'}], num2cell(sumMff)]];
% 
% xlswrite(['ret_sum_k_',num2str(k),'.xlsx'],outSum)
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Momentum Issue Testing...
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % retC = retData{1,1};
% % % retM = [];
% % % for t = 12:size(retC,1)-1,
% % %     
% % %     retCC = prod(1+retC(t-11:t,:))-1;
% % %     
% % %     moQuan = [.1 .9];
% % %     quan = [ min(retCC), quantile(retCC,moQuan),...
% % %         max(retCC) ];
% % %     
% % %     for q = 1:(size(quan,2)-1)
% % %         idx = (retCC >= quan(q)) & (retCC <= quan(q+1));
% % %         retMtemp(q) = nanmean(retC(t+1,idx));
% % %     end,
% % %     retM = [retM ; retMtemp(end) - retMtemp(1)];
% % %     
% % % end,
% % % 
% % % cumRetM = cumprod(1+retM)-1;
% % % plot(cumRetM)
% 
% 
% 
% 
% 
