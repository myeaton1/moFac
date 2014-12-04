function [outputs, param_out, residual,LIK] = kalman_it(k, win ,startY, startM, F, B, param_in, lab)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   kalman_it.m
%	Matt Yeaton, 06/19/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Kalman Filter Iteration:
% This function runs the Kalman filter optimization 
% for a single block of time

%%% Mandatory Arguments:
% k: 		# of principal components
% win: 		length of the time window, in years
% startY: 	first year in the time block
% startM: 	firmst month in the time block
% F: 		(t x k) vector of latent factor observations
% 			(principal component coefficients)
% B: 		(k x n) matrix of factor loadings
% param_in: starting location in search space
% lab: 		a cell array containing the date range label

%%% Output
% outputs:
% param_out:
% residual: 
% LIK:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Set key variables, optimization options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

year 					= 	startY;
mon 					= 	startM;

% Adjust for extra year for mon != 1
if mon 					== 	1
    mon2 					= 	12;
    extra 					= 	0;
else
    mon2 					= 	mon - 1;
    extra 					= 	1;
end

maxIterations 			= 	3000; % was 5e3
maxFunctionEvals 		= 	maxIterations*10; % was 1e4
options 				= 	optimset('Display','iter','TolX',1e-6,'TolFun',1e-6,...
								'MaxIter',maxIterations,'MaxFunEvals',...
								maxFunctionEvals,'Display','final');

% number of times the optimization is run
num_est_rounds 			= 	10;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for est_round 			= 	1:num_est_rounds

	fprintf(['\n',lab, '...\n']);
    fprintf(['est_round = ',num2str(est_round), '\n']);

    if mod(est_round,2),
        [param_out,obj,residual,exitflag] ...
        						= 	fminunc(@loglik_kalman_filter_model4,param_in,options,F(:,1:k),0);
    else
        [param_out,obj,residual,exitflag] ...
        						= 	fminsearch(@loglik_kalman_filter_model4,param_in,options,F(:,1:k),0);
    end;

	% Update input parameters as last round's parameters
	param_in 				= 	param_out;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Update Zi_t, etc. using optimal parameters 
%		found above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LIK, outputs] 			= 	loglik_kalman_filter_model4(param_out,F(:,1:k),1);