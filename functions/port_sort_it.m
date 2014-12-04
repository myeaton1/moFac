function [port, portIdx] = port_sort_it(k, win ,startY, startM, B, outputs, retData, plist, lab, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Portfolio Sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd(outputpath)

TB = 0;					% Top/bottom or quantile sorts
% for TB = 0:1,
tbNum = 10;				% Number of firms in each top/bottom bin
moQuan = .2:.2:.8;		% If TB == 0, the quantiles used for
						%	high/low portfolios
year = startY;
mon = startM;
                        
tbLab = {['T', num2str(moQuan(2)*100), '%/B', num2str(moQuan(1)*100), '%'],...
    ['T', num2str(tbNum), '/B', num2str(tbNum)]};
tbFileLab = {['TB_quan',num2str(moQuan(1)*100)], ['TB_num', num2str(tbNum)]};

% fac.portRet = [];
% fac.lambdaSer = [];
% fac.phiSer = [];
% fac.lab = cell(0);

fprintf(['\n',lab, '...\n']);

if mon == 12,
    yearF = year+1;
    monF = 1;
else
    yearF = year;
    monF = mon+1;
end,

B = B(:,1:k);
ret = retData';
% retF = retData{yeF,monF}(:,end);

% ret1 = retData{ye,mon}(:,2:end)';
% ret1(isnan(ret1)) = -1000;
% ret2 = retData{yeF,monF}(:,1:end-1)';
% ret2(isnan(ret2)) = -1000;

% if size(ret2,2) < size(ret1,2),
%     retCorr = corr(ret2,ret1);
%     idx = logical(sum(abs(retCorr-1) < 1e-10));
%     B = B(idx,:);
% else
%     retCorr = corr(ret1,ret2);
%     idx = logical(sum(abs(retCorr-1) < 1e-10));
%     retF = retF(idx,:);
% end,

lambda = outputs.xi_TT_out;
lambda = reshape(lambda,size(lambda,1),size(lambda,3));
lambda = lambda(1:k,end)';

er = B*lambda';

if TB == 1,
    erSort = sort(er);
    quan = [ min(er), erSort(tbNum), erSort(end-tbNum+1), max(er) ];
else
    quan = [ min(er), quantile(er,moQuan), max(er) ];
end,

port = cell(1,(size(quan,2)-1));
portIdx = cell(1,(size(quan,2)-1));
for q = 1:(size(quan,2)-1)
    idx = (er >= quan(q)) & (er <= quan(q+1));
    port{q} = plist(idx);
    portIdx{q} = idx;
    % port(q) = nanmean(retF(idx));
end,

% fac.portRet = [fac.portRet ; (port(end) - port(1)) ];
% fac.lambdaSer = [fac.lambdaSer; lambda];

% fac.lab = [fac.lab ; [num2str(year(ye+win)), 'm', num2str(mon)]];

