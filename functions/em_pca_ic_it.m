function [f, B, k, tstat, ICp, retData, pList, lab] = em_pca_ic_it(win ,startY, startM, retM, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   em_pca_it.m
%	Matt Yeaton, 06/19/2014
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
% retM:     momentum factor struct

%%% Optional Arguments:
% thres: 	firms are excluded if they do not have 
% 				observations on at least this percentage 
% 				of days over the sample
% 				(default = 0.5)
% delta: 	convergence criteria for the EM algorithm
% 				(default = 1e-5)
% monthly: 	run EM algorithm for monthly returns (instead
% 				of daily)
% 				(default = 1)

%%% Output
% f:			(t x pcs) vector of latent factor observations
% 				(principal component coefficients)
% B:			(n x pcs) matrix of factor loadings
% retData:	returns data for the time period
% lab:		a cell containing the date range label


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Deal with optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numvarargs              =   length(varargin);

% set defaults for optional inputs
optargs                 =   {0.5 1e-5 1};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs)   =   varargin;

% Place optional args in memorable variable names
[thres, delta, monthly] =   optargs{:};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EM-PCA Estimation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OS check, set relevant data paths
if isunix,
    datapath    =   '/home/rcemry01/data/year_grids/';
else
    datapath    =   'S:\MattY\Emanuel_BigDataFile\data\year_grids\';
end,

fprintf('EM-PCA Estimation...');

% Set up month ranges
mon     =   startM;
if  mon ==  1,
    mon2        =   12;			% ending month of date range
    extra       =   0;			% Adjust for extra year for mon ~= 1
                    			%	(10-year block have portions in both
                    			%	the 1st and 11th year if mon ~= 1)
else
    mon2        =   mon - 1;
    extra       =   1;				
end,

% Generate date range label
lab     =   [num2str(startY), 'm', num2str(mon), '-'...
                ,num2str(startY+win-1+extra), 'm', num2str(mon2)];

fprintf(['\n',lab, '...\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Find firms that span the sample period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for     roll    =   1:(win+extra),

    % set current year of loop
    y           =   startY+roll-1;

    % load returns data in a permno x date matrix
    yg          =   dlmread([datapath, 'y_grid_',num2str(y),'.csv'],',');
    permnos     =   yg(1,2:end)';

    % extract permnos that span the period
    if      roll    ==  1,
        % running list of permnos in each year of the period
        pList       =   permnos;
    elseif  roll    >   1,
        % intersect the running list with the current list for
        %   each year in the period
        pList       =   union(pList,permnos);
    end,

end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Extract preliminary data for the firms
%       who span the sample period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x   =   [];
for     roll    =   1:(win+extra),

    % set current year of loop
    y           =   startY+roll-1;

    % load returns data in a permno x date matrix
    yg          =   dlmread([datapath, 'y_grid_',num2str(y),'.csv'],',');
    permnos     =   yg(1,2:end)';
    date        =   yg(2:end,1);

    % index of firms who span the sample period
    idx         =   ismember(pList,permnos);

    % returns data for the given year
    xrRaw       =   yg(2:end,2:end)';

    % only keep returns for those who span the sample period
    xr          =   nan(length(pList),size(xrRaw,2));
    xr(idx,:)   =   xrRaw;

    % append returns and dates data to get data for the
    %   entire period
    if roll     ==  1,
        x           =   xr;             % returns matrix 
        fD          =   date;
        
    elseif roll > 1,
        x           =   [x, xr];
        fD          =   [fD ; date];
    end,

end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Transforms into monthly data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  monthly     ==  1,
    fDM             =   round(fD/100);    % generate monthly dates (YYYYMM) from daily dates
    fDMu            =   unique(fDM);      % find the unique monthly dates for final series
    xNoNan          =   x;
    xNoNan(isnan(xNoNan)) ...
                    =   0;

    % Produce monthly returns by cumulating daily returns
    xM              =   nan(size(x,1),length(fDMu));
    for mi          =   1:length(fDMu),
        m               =   fDMu(mi);                   % current month
        idxM            =   fDM == m;                   % index of daily dates for  the 
                                                        %   current month
        xM(:,mi)        =   prod(1+xNoNan(:,idxM),2)-1;      % assign current month's returns
    end,
    
    if mon          ~=  1,
        xM              =   xM(:,mon:(end-(12-mon2)));
        fDMu            =   fDMu(mon:(end-(12-mon2)));
    end,
end,

% Check threshold requirements
idxMissInit     =   ~(sum(~isnan(x),2)/size(x,2) < thres);     % threshold checker
pList           =   pList(idxMissInit);                       % updated permno list
x               =   xM;

% EM PCA
xest_1          =   x(idxMissInit,:)';
xest_1(xest_1 == 0) ...
                =   NaN;
    
kList           =   1:6;

parfor k        =   kList,
    
    fprintf(['\nk = ', num2str(k),'...\n']);
    
    err             =   10;
    loop            =   1;
    xest            =   x(idxMissInit,:)';
    xest(xest == 0) =   NaN;
    xestAll         =   nan(size(xest));
    while (err      >   delta),

        if (loop        >   500),
            break
        end,

        idxMiss         =   ~(sum(isnan(xest),1) > 0);
        xc              =   xest(:,idxMiss);

        coeff           =   princomp(zscore(xc));
        loadings        =   coeff(:,1:k);
        f               =   xc*loadings;

        lambda          =   nan(size(xest_1,2),k);
        for i           =   1:size(xest_1,2),
            xi              =   xest_1(:,i);
            idx             =   ~isnan(xi);
            xi              =   xi(idx);
            ftemp           =   f(idx,:);
            lambda(i,:)     =   (ftemp'*ftemp)\ftemp'*xi;
        end,

        for j           =   1:size(xest_1,2),
            for t           =   1:size(xest,1),
                if isnan(xest_1(t,j)),
                    xest(t,j)       =   lambda(j,:)*f(t,:)';
                end,

                xestAll(t,j)    =   lambda(j,:)*f(t,:)';     % every x estimated via factor model
            end,
        end,

        if loop         ==  1
            fOld            =   f + 10;
        end
        
        err             =   sum(sum((f - fOld).^2));
        loop            =   loop+1;
        if loop         ==  10,
            fprintf('Iteration: ');
        end,
        if loop/10 == round(loop/10),
            fprintf([num2str(loop),'...']);
        end,
        fOld = f;
    end,
    
    % Momentum Regs
    X               =   [ones(120,1),f];

    startIdx        =   12*((startY-9) - retM.startY) + mon;
    y               =   retM.retM19(startIdx:startIdx+120-1);

    b               =   (X'*X)\X'*y;
    e               =   y - X*b;

    nwsei           =   NeweyWest(e,X);
    tstat           =   b./nwsei;
    tstat           =   tstat(1);
    

    % Bai,Ng IC_p
    [T,N]           =   size(xest);
    NT              =   sum(sum(~isnan(xest_1)));

    e               =   xestAll - xest;
    Sigma           =   mean(sum(e.*e/T));

    ICp1            =   log(Sigma) + k*((N+T)/NT)*log(NT/(N+T));    
    
    % Keep relevant output
    fList{k}        =   f;
    lambdaList{k}   =   lambda;
    ICp1List(k)     =   ICp1;
    tstatList(k)    =   tstat;
    

end,

[~,argmin]      =   min(ICp1List);

f               =   fList{argmin};
B               =   lambdaList{argmin};
k               =   argmin;
tstat           =   tstatList;
ICp             =   ICp1List;
retData         =   xest_1;













