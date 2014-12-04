function FT_mat = sample_facs(g,alpha,lambda_tilde,psi_F,sig2_F,psi_G,sig2_G,params);

B = params.B;
q_G = params.q_G;
max_q_G = max(q_G);
K_G = params.K_blocks;
K_F = params.K_facs;
q_F = params.q_F;
max_q_F = max(q_F);

Gmat = [];
sig2_G_vec = [];
psi_G_mat = [];

for b=1:B;
    Gmat=[Gmat g{b}];
    sig2_G_vec=[sig2_G_vec; sig2_G{b}];
    
    if max_q_G > 0;
        psi_G_mat=[psi_G_mat; psi_G{b} zeros(K_G(b),max_q_G-q_G(b))];
    else;
        psi_G_mat=[psi_G_mat; zeros(K_G(b),1)];    
    end;
end;
[T,N]=size(Gmat);
ystar=zeros(T,N);
for i=1:N;
    ystar(:,i)=filter([1 -psi_G_mat(i,:)],1,Gmat(:,i));
end;
ystar = ystar';


% transition matrix of the filtered data is r= p+1+q+1
ndim=size(lambda_tilde,2);
psi_F_mat=zeros(K_F,max_q_F*K_F);
for r=1:K_F;
    indx=r:K_F:K_F*max_q_F;
    psi_F_mat(r,indx)=psi_F(r,:);
end;

A = zeros(ndim);
A(1:size(psi_F_mat,1),1:size(psi_F_mat,2)) = psi_F_mat;
A(K_F+1:end,1:ndim-K_F) = eye(ndim-K_F);
H = lambda_tilde;
QQ=zeros(ndim,ndim);  
QQ(1:K_F,1:K_F)=diag(sig2_F);
R=diag(sig2_G_vec);
Alpha=[alpha zeros(rows(alpha),ndim-cols(alpha))];

% initialize Kalman filter
Gtt=zeros(ndim,1);
P00=inv(eye(ndim^2)-kron(A,A))*vec(QQ);
Ptt=reshape(P00,ndim,ndim);

% Kalman filter recursion
for t=1:T,
    Gtt1=Alpha(t,:)'+A*Gtt;         % x(t|t-1)
    Ptt1=A*Ptt*A'+QQ;               % P(t|t-1)
    ett1=ystar(:,t)-H*Gtt1;         % eta(t|t-1)=y(t)- y(t|t-1)= y(t)-H*x(t|t-1)
    v_ett1= H*Ptt1*H'+R;            % var(eta(t|t-1))
    k_gn= Ptt1*H'/v_ett1;           % K=P(t|t-1)H'/ f(t|t-1)
    Gtt = Gtt1+k_gn*ett1;           % x(t|t)=x(t|t-1)+ K eta(t|t-1)
    Ptt= (eye(ndim)-k_gn*H)*Ptt1;   % P(t|t)= P(t|t-1)-K*H*P(t|t-1)
    Fmat(:,t)=Gtt;
    Pmat{t}=Ptt;
end;

% Carter-Kohn backward sampling algorithm
FT_mat=zeros(T,K_F); 
G =chol(Pmat{T})';
FT = Fmat(:,T)+G*randn(ndim,1);
FT_mat(T,1:K_F)=FT(1:K_F)';
for t=T-1:-1:1;
    jj=1:K_F;
    etT = FT(jj)-Alpha(t+1,jj)'-A(jj,:)*Fmat(:,t);                     
    v_etT=A(jj,:)*Pmat{t}*A(jj,:)'+QQ(jj,jj); 
    k_gn0=Pmat{t}*A(jj,:)'/v_etT;                
    Fmat(:,t)=Fmat(:,t)+k_gn0*etT;                       
    Pmat{t}= (eye(ndim)-k_gn0*A(jj,:))*Pmat{t};  
    G=chol(Pmat{t})';
    FT=Fmat(:,t)+G*randn(ndim,1);
    FT_mat(t,1:K_F)=FT(1:K_F,1)';
end; % end t

done = 1;

