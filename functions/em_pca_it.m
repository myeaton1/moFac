function [f, B, retData, pList, lab] = em_pca_it(pcs, win ,startY, startM, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   em_pca_it.m
%	Matt Yeaton, 05/02/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% EM-PCA Iteration:
% This function runs EM-PCA for a single block of time

%%% Mandatory Arguments:
% pcs: 		# of principal components over which to run the 
% 				EM algorithm
% win: 		length of the time window, in years
% startY: 	first year in the time block
% startM: 	firmst month in the time block

%%% Optional Arguments:
% thres: 	firms are excluded if they do not have 
% 				observations on at least this percentage 
% 				of days over the sample
% 				(default = 0.9)
% delta: 	convergence criteria for the EM algorithm
% 				(default = 1e-5)
% monthly: 	run EM algorithm for monthly returns (instead
% 				of daily)
% 				(default = 1)

%%% Output
% f:			(pcs x 1) vector of latent factor observations
% 				(principal component coefficients)
% B:			(pcs x # of firms) matrix of factor loadings
% retData:	returns data for the time period
% lab:		a cell containing the date range label


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Deal with optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numvarargs = length(varargin);

% set defaults for optional inputs
optargs = {0.9 1e-5 1};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[thres, delta, monthly] = optargs{:};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   If running as a script rather than a function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
script = 0;

if script,
    
    clear, home, close all;
    
    pcs = 5;
    win = 10;
    startY = 1955;
    startM = 1;
    thres = .9;
    delta = 1e-5;
    monthly = 1;

end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EM-PCA Estimation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OS check, set relevant data paths
if isunix,
    datapath = '/home/rcemry01/data/year_grids/';
else
    datapath = 'S:\MattY\Emanuel_BigDataFile\data\year_grids\';
end,

% shvar = nan(length(year), pcs);
% fSer = cell(length(year)-win,12);
% fSerD = cell(length(year)-win,12);
% loadList = cell(length(year)-win,12);	
% retData = cell(length(year)-win,12);

fprintf('EM-PCA Estimation...');

% Set up month ranges
mon = startM
if mon == 1,
    mon2 = 12;				% ending month of date range
    extra = 0;				% Adjust for extra year for mon ~= 1
    						%	(10-year block have portions in both
    						%	the 1st and 11th year if mon ~= 1)
else
    mon2 = mon - 1;
    extra = 1;				
end,

% Generate date range label
lab = [num2str(startY), 'm', num2str(mon), '-'...
    ,num2str(startY+win-1+extra), 'm', num2str(mon2)];

fprintf(['\n',lab, '...\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Find firms that span the sample period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for roll = 1:(win+extra),

    % set current year of loop
    y = startY+roll-1;

    % load returns data in a permno x date matrix
    yg = dlmread([datapath, 'y_grid_',num2str(y),'.csv'],',');
    permnos = yg(1,2:end)';

    % extract permnos that span the period
    if roll == 1,
        % running list of permnos in each year of the period
        pList = permnos;
    elseif roll > 1,
        % intersect the running list with the current list for
        %   each year in the period
        pList = intersect(pList,permnos);
    end,

end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Extract preliminary data for the firms
%       who span the sample period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = [];
for roll = 1:(win+extra),

    % set current year of loop
    y = startY+roll-1;

    % load returns data in a permno x date matrix
    yg = dlmread([datapath, 'y_grid_',num2str(y),'.csv'],',');
    permnos = yg(1,2:end)';
    date = yg(2:end,1);

    % index of firms who span the sample period
    idx = ismember(permnos, pList);

    % returns data for the given year
    xr = yg(2:end,2:end)';

    % only keep returns for those who span the sample period
    xr = xr(idx,:);

    % append returns and dates data to get data for the
    %   entire period
    if roll == 1,
        x = xr;             % returns matrix 
        fD = date;
    elseif roll > 1,
        x = [x, xr];
        fD = [fD ; date];
    end,

end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Transforms into monthly data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if monthly == 1,
    fDM = round(fD/100);    % generate monthly dates (YYYYMM) from daily dates
    fDMu = unique(fDM);     % find the unique monthly dates for final series

    % Produce monthly returns by cumulating daily returns
    xM = nan(size(x,1),length(fDMu));
    for mi = 1:length(fDMu),
        m = fDMu(mi);                       % current month
        idxM = fDM == m;                    % index of daily dates for  the 
                                            %   current month
        xM(:,mi) = prod(1+x(:,idxM),2)-1;   % assign current month's returns
    end,

    x = xM;
    
    if mon ~= 1,
        x = x(:,mon:(end-(12-mon2)));
        fDMu = fDMu(mon:(end-(12-mon2)));
    end,
end,

% Check threshold requirements
idxMiss = ~(sum(isnan(x),2)/size(x,2) > thres);     % threshold checker
pList = pList(idxMiss);                             % updated permno list

% EM PCA
err = 10;
loop = 1;
xest = x(idxMiss,:)';
xest_1 = xest;
while (err > delta),

    if (loop > 500),
        break
    end,

    idxMiss = ~(sum(isnan(xest),1) > 0);
    xc = xest(:,idxMiss);

    [coeff, ~, latent] = princomp(zscore(xc));
    loadings = coeff(:,1:pcs);
    eigs = latent(1:pcs);
    f = xc*loadings;

    lambda = nan(size(xc,1),pcs);
    for i = 1:size(x,1),
        xi = x(i,:)';
        idx = ~isnan(xi);
        xi = xi(idx);
        ftemp = f(idx,:);
        lambda(i,:) = (ftemp'*ftemp)\ftemp'*xi;
    end,

    for j = 1:size(xest,2),
        for t = 1:size(xest,1),
            if isnan(xest_1(t,j)),
                xest(t,j) = lambda(j,:)*f(t,:)';
            end,
        end,
    end,
    

    if loop > 1,
        err = sum(sum((f - fold).^2));
    end,

    loop = loop+1;
    if loop == 10,
        fprintf('Iteration: ');
    end,
    if loop/10 == round(loop/10),
        fprintf([num2str(loop),'...']);
    end,

    fold = f;
end,

B = lambda;
retData = xest_1;