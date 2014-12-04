function e = compute_resids(obs,fac,Lambda,params);

K_F = params.K_facs;
l_F = params.l_F;
max_l_F = max(l_F);

Fmat=fac;
for j=1:max_l_F;           
    Fmat=[Fmat lagn(fac,j)];  
end;  
fit=zeros(size(obs));
for j=1:K_F;
  indx=j:K_F:K_F*(max_l_F+1);
  fit=fit+Fmat(:,indx)*Lambda(:,indx)';
end; 
e=obs-fit;