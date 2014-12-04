function [retM19,retM37] = mo_ser_generator(startY)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Momentum Portfolios
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OS check, set relevant data paths
if isunix,
    datapath = '/home/rcemry01/data/year_grids/';
else
    datapath = 'S:\MattY\Emanuel_BigDataFile\data\year_grids\';
end,

endY = 2012;
yList = (startY-1):endY;

xAll = [];
pAll = [];
for y = yList

    yg          =   dlmread([datapath, 'y_grid_',num2str(y),'.csv'],',');
    permnos     =   yg(1,2:end)';
    fD          =   yg(2:end,1);
    x           =   yg(2:end,2:end)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Transforms into monthly data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    pAll = [pAll ; repmat(permnos',12,1)];
    xAll = [xAll ; x'];

end

%%%%%%%%%%%%%%
% Mo
%%%%%%%%%%%%%%
retM19 = [];
retM37 = [];
retMq = [];
for d = 12:size(xAll,1)-1
    
    for dN = d-11:d+1
        if dN == d-11,
            % running list of permnos in each year of the period
            pList = pAll(d-11,:);
        elseif dN > d-11,
            % intersect the running list with the current list for
            %   each year in the period
            pList = intersect(pList,pAll(dN,:));
        end,
    end
    
    for dN = d-11:d+1
        % index of firms who span the sample period
        idx = ismember(pAll(dN,:), pList);

        % only keep returns for those who span the sample period
        xd = xAll(dN,idx);
        pd = pAll(dN,idx);

        % append returns and dates data to get data for the
        %   entire period
        if dN == d-11,
            xD = xd;             % returns matrix 
            pD = pd;
        elseif dN > d-11,
            xD = [xD; xd];
        end,
        
    end 
   
    retL    =   prod(1+xD(1:end-2,:))-1;
    retF    =   xD(end,:);
    
    moQuan = .1:.1:.9;
    quan = [ min(retL), quantile(retL,moQuan),max(retL) ];
    
    for q = 1:(size(quan,2)-1)
        idx = (retL >= quan(q)) & (retL <= quan(q+1));
        retMtemp(q) = nanmean(retF(idx));
    end
    
    retM19 = [retM19 ; retMtemp(end) - retMtemp(1)];
    retMq = [retMq ; retMtemp];   
    
    moQuan = .3:.1:.7;
    quan = [ min(retL), quantile(retL,moQuan),max(retL) ];
    
    for q = 1:(size(quan,2)-1)
        idx = (retL >= quan(q)) & (retL <= quan(q+1));
        retMtemp(q) = nanmean(retF(idx));
    end
    
    retM37 = [retM37 ; retMtemp(end) - retMtemp(1)];
    
end

% retM19      =   reshape(retM19,12,endY-startY+1)';
% retM37      =   reshape(retM37,12,endY-startY+1)';

y           =   (startY:endY)';
% cumRetM19 = cumprod(1+retM19)-1;
% cumRetM37 = cumprod(1+retM37)-1;
% plot([cumRetM19,cumRetM37])            


















