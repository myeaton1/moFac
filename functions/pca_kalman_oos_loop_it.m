function pca_kalman_oos_loop_it(ye, mon, k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Set relevant paths, initialize key variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% k = 5;			% # PCAs in Kalman filter
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run PCA, Kalman filter, portfolio sort for a
%       given month
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
year = startYears(ye);

[f{ye,mon}, B{ye,mon}, retData{ye,mon}, pList{ye,mon}, lab{ye,mon}] = ...
    em_pca_it(k, win ,year, mon);


% Initialize parameters 
Phi = .7*eye(k);
Sigw = .001*ones(k,1);
Sigf = var(f{ye,mon})';
params0 = [diag(Phi); Sigw; Sigf];
param_in = params0;
nDim = length(param_in);


for i = 1:nRuns,
    op.xi_TT_out = [];
    if i > 1,
        while isempty(op.xi_TT_out),
            n = randn(nDim,1).*magRand;
            test{i} = param_in+n;
            [~, op] = loglik_kalman_filter_model4(test{i},f{ye,mon},1);

        end,
    else
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
        


save(outfile);






