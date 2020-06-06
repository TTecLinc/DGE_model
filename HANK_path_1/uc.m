function y=uc(c,n)
global t r pen aopt w tau tr
global s gam psi beta
y= (c+psi).^(-s).*n^(gam.*(1-s));
end