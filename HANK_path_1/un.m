function y=un(c,n)
global t r pen aopt w tau tr
global s gam psi beta
y= gam*(c+psi)^(1-s).*n^(gam*(1-s)-1);
end