function y=lagn(x,n)
[nt,nc]=size(x);
if n> 0
x1=trimr(x,0,n);
 y=[zeros(n,nc); x1];
end
if n<0
x1=trimr(x,abs(n),0);
y=[x1;zeros(abs(n),nc)];
end
if n==0
  y=x;
end;
