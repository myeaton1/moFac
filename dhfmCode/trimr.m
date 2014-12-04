function xx=trimr(x,a,b);


[nt,nc]=size(x);
if a >= 0 ;
xx=x(a+1:nt,:);
end;
if b >=0;
if a > 0; x=xx; end;
[nt,nc]=size(x);
xx=x(1:nt-b,1:nc);
end;
