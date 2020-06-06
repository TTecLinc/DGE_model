function z=u(c,l,s)
global t r pen aopt w tau tr
global s gam psi beta
if s==1;
    z=ln(c+psi)+gam*ln(l);
else
    z=(((c+psi)*l^gam)^(1-s)-1)/(1-s);
end
end