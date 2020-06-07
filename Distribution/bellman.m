function z=bellman(a0,a1,y)
global tau r w b na pp a beta e i 
global sigma vold
   if y==1;
      c=(1+(1-tau)*r)*a0+(1-tau)*w-a1;
   else
      c=(1+(1-tau)*r)*a0+b-a1;
   end
   if c<0;
      z=-1e5;
   end
   if a1>=a(na);
      z=a1.^2*(-1e10);
   end
   z=u(c)+beta*(pp(y,1)*value(a1,1)+pp(y,2)*value(a1,2));
end

