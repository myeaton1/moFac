% Preset parametes
function [LIK, outputs] = loglik_kalman_filter_model1new(params,y,DO_PRODUCE_P);

[T,k] = size(y);

H = [zeros(k) eye(k) eye(k)]';
R = zeros(k);
Phi = diag(params(1:k));
params(1:k) = [];
Sigw = diag(abs(params(1:k)));
params(1:k) = [];
Sigf = diag(abs(params(1:k)));
   
if max(abs(eig(Phi))) > 0.99 || any(diag(Phi < 0)) || any(diag(Sigw < 0)) || any(diag(Sigf < 0)),% || ~isempty(find((diag(Sigw) > diag(Sigf)))),
    LIK = 1e30;
    outputs.xi_TT_out = [];
    
    return;
else,
    a = [zeros(3*k,1)];
    F = [Phi zeros(k) zeros(k);eye(k) zeros(k) zeros(k); zeros(k) zeros(k) zeros(k)]; 
    Q = [Sigw zeros(k) zeros(k);zeros(k) zeros(k) zeros(k);zeros(k) zeros(k) Sigf];
    P00 = eye(3*k);
    xi0 = [mean(y)' ; mean(y)' ; zeros(k,1)];
    
    %
%     varXi = reshape((eye(9*k^2)-kron(F,F))\Q(:),3*k,3*k);
%     varZ = H'*varXi*H;
    
    % Compute sequence of conditional variance covariance matrices which are the same for each individual and each simulation
    PTL = P00;
    xiTL = xi0;
    LIK = 0;
    xi_TT_out = []; xi_TL_out = [];
    P_TT_out = []; P_TL_out = [];
    for t = 1:T,
       inv_vTL = inv(H'*PTL*H + R);
       K = PTL*H*inv_vTL;
       eps = y(t,:)' - H'*xiTL;
       xiTT = xiTL + K*eps;
       xi_TT_out = cat(3,xi_TT_out,xiTT);
       PTT = PTL-K*H'*PTL; 
       if DO_PRODUCE_P,
           xi_TL_out = cat(3,xi_TL_out,xiTL);
           P_TL_out = cat(3,P_TL_out,PTL);
           P_TT_out = cat(3,P_TT_out,PTT);
       end;
       xiTL = a + F*xiTT;     
       yfit(t,:) = (H'*xiTT)';
       PTL = F*PTT*F' + Q;   
       LIK = LIK + 0.5*log(1/det(inv_vTL)) + 0.5*(eps'*inv_vTL*eps);
    end;

end;

outputs.xi_TT_out = xi_TT_out; 
outputs.P_TT_out = P_TT_out; 
outputs.xi_TL_out = xi_TL_out; 
outputs.P_TL_out = P_TL_out; 
% outputs.varZ = varZ;

outputs.yfit = yfit;
outputs.R = R;
outputs.a = a;
% outputs.mu = mu;
outputs.Phi = Phi;
outputs.Sigw = Sigw;
outputs.Sigf = Sigf;
outputs.F = F;
outputs.H = H;


