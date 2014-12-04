function [b, ffData, retData, pList, lab] = ff_regs_it(win ,startY, startM, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ff_regs_it.m
%	Matt Yeaton, 06/27/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Fama-French Regression Iteration:
% This function runs a regression of raw returns data on
% the Fama-French 3 factors for a single block of time

%%% Mandatory Arguments:
% win: 		length of the time window, in years
% startY: 	first year in the time block
% startM: 	firmst month in the time block

%%% Optional Arguments:
% thres: 	firms are excluded if they do not have 
% 				observations on at least this percentage 
% 				of days over the sample
% 				(default = 0.25)
% monthly: 	run EM algorithm for monthly returns (instead
% 				of daily)
% 				(default = 1)

%%% Output
% b:        (n x 3) beta matrix
% ffData:   (t x 3) raw Fama-French 3 factors for the 
%               time period
% retData:	(n x t) returns data for the time period
% pList:    a cell array containing firms in the sample
% lab:		a cell array containing the date range label


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Deal with optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numvarargs              =   length(varargin);

% set defaults for optional inputs
optargs                 =   {0.25 1};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs)   =   varargin;

% Place optional args in memorable variable names
[thres, monthly]        =   optargs{:};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EM-PCA Estimation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OS check, set relevant data paths
if isunix,
    datapath    =   '/home/rcemry01/data/year_grids/';
    ffPath      =   '/home/rcemry01/data/';
else
    datapath    =   'S:\MattY\Emanuel_BigDataFile\data\year_grids\';
    ffPath      =   'W:\data\';
end,

fprintf('Fama-French beta Estimation...');

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
%       in the sample period
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
% pList           =   pList(idxMissInit);                       % updated permno list
x               =   xM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Import Fama-French Factors
%   (order: mkt-rf, smb, hml, rf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ffData          =   dlmread([ffPath,'f_f_factors_monthly.csv'],',');
ffYear          =   floor(ffData(:,1)/100);
ffMon           =   ffData(:,1) - ffYear*100;
ffData          =   ffData(:,2:end-1)/100;

idxFF           =   ffYear >= startY & ffYear <= startY + 9 + extra;
ffData          =   ffData(idxFF,:);
ffMon           =   ffMon(idxFF);
mon1FF          =   find(ffMon == mon,1,'first');
mon2FF          =   find(ffMon == mon2,1,'last');
ffData          =   ffData(mon1FF:mon2FF,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Regress returns data on factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X               =   [ones(120,1),ffData];

y               =   x';

b               =   (X'*X)\X'*y;
b               =   b';
b               =   b(:,2:end);

retData         =   x;













