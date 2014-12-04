%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ff_regs_kalman_oos_loop_par.m
%   Matt Yeaton, 06/27/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% EM-PCA Kalman Filter for full sample:
% Parallelized EM-PCA generator and Kalman filter for
% the resultant factors

%%% Input:
% matrices of daily returns (date x firm) for a given
% year

%%% Output:
% dated .mat file with output organized by (year, month)
% pairs

clear, home, close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Set relevant paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isunix,
    resultspath         =   '/home/rcemry01/data/kalman_sort_res/';
    outputpath          =   '/home/rcemry01/data/output';
    codepath            =   '/home/rcemry01/code/factor';
else
    resultspath         =   'Y:\data\kalman_sort_res';
    outputpath          =   'Y:\data\kalman_sort_res\output';
    codepath            =   'S:\MattY\emanuel\Momentum\Code';
end,

addpath(codepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Claim workers for matlabpool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRuns                   =   15;

matlabpool(nRuns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Initialize key variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
win                     =   10;                 % length (in years) of EM PCA and Kalman filter windows
start                   =   1956;               % first starting year
startYears              =   start:(2012-win);   % list of starting years

outfile                 =   fullfile(resultspath, ['kalman_oos_k_ffFactors_startY_',num2str(start),...
                                datestr(today,'yyyy_mm_dd'), '_LR.mat']);
k                       =   3;                  % Fama-French 3 factors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Assign appropriate magnitudes for 
%       Kalman filter starting search-
%       space perturbations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
magRandPhi              =   10e-2;
magRandSigw             =   10e-4;
magRandSigf             =   10e-2;
magRand                 =   [repmat(magRandPhi,k,1);repmat(magRandSigw,k,1);...
                                repmat(magRandSigf,k,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Year-Month loop, runs main functions for each
%       (year, month) pair
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ye                  =   1:length(startYears),

	for mon                 =   1:12,

		year                =   startYears(ye);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Fama-French factor regressions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [b{ye,mon},ffData{ye,mon},retData{ye,mon}, pList{ye,mon}, lab{ye,mon}] ...
                            =   ff_regs_it(win ,year, mon);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Kalman Filter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Initialize parameters
        if ye       ==  1   &&  mon     ==  1
            Phi                 =   .7*eye(k);
            Sigw                =   .001*ones(k,1);
            Sigf                =   var(ffData{ye,mon})';
            params0             =   [diag(Phi); Sigw; Sigf];
            param_in            =   params0;
            nDim                =   length(param_in);
        end

        n                   =   cell(nRuns,1);
        test                =   cell(nRuns,1);
        testOutputs         =   cell(nRuns,1);
        testParam_out       =   cell(nRuns,1);
        testExitflag        =   cell(nRuns,1);
        testLIK             =   cell(nRuns,1);
        
        parfor i            =   1:nRuns,
            op.xi_TT_out        =   [];

            if i                >   1,
                while isempty(op.xi_TT_out),
                    n{i}                =   randn(nDim,1).*magRand;
                    test{i}             =   param_in+n{i};
                    [~, op]             =   loglik_kalman_filter_model4(test{i},ffData{ye,mon},1);

                end,
            else
                test{i}             =   param_in;
            end,
            
            [testOutputs{i}, testParam_out{i}, testExitflag{i},testLIK{i}] ...
                                =   kalman_it(k, win ,year, mon, ffData{ye,mon}, b{ye,mon}, ...
                                        test{i}, lab{ye,mon});
            
        end,
        
        [~,argmin]          =   min(cell2mat(testLIK));
        
        outputs{ye,mon}     =   testOutputs{argmin};
        param_out{ye,mon}   =   testParam_out{argmin};
        exitflag{ye,mon}    =   testExitflag{argmin};
        LIK{ye,mon}         =   testLIK{argmin};

        param_in            =   testParam_out{argmin};
        
	end,
end,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Save output, close pool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(outfile);

matlabpool close