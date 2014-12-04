function X=switch_rows(X,id);
id1=id ; id2=flipud(id);
temp=X;
X(id1,:)=temp(id2,:);
X(id2,:)=temp(id1,:);
