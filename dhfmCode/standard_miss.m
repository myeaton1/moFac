function x=standard_miss(yy,mask);
T=size(yy,1);
N=size(yy,2);
x=yy;
one=ones(T,1);
for i=1:N;
    if isnan(mask),
        bad = find(isnan(yy(:,i)));
    else,
        bad=find(yy(:,i)==mask);
    end;
  good=one;
  good(bad)=0;
  igood=find(good);
  y=yy(igood,i);
  z= (y-mean(y))/std(y); 
  x(igood,i)=z;
end;
