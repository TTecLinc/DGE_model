function y=lininter(xd,yd,x)
global tau r w b na pp a beta e i
global sigma vold
  j=sum(xd<=x);
  if j>200;
      j=200;
  end
  y=yd(j)+(yd(j+1)-yd(j)).*(x-xd(j))/(xd(j+1)-xd(j));
end