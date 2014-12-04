function X=switch_cols(X,id);
id1=id ; id2=flipud(id')';
temp=X;
X(:,id1)=temp(:,id2);
Y(:,id2)=temp(:,id1);
