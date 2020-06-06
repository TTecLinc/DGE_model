function y=rfyoung(x,i)
global t r pen aopt w tau tr
global s gam psi beta nopt
    k0=x(1);
    n0=x(2);
    k1=aopt(i+1);
    k2=aopt(i+2);
    if i==t;
        n1=0;
        c1=(1+r)*k1+pen-k2;
    else
        n1=nopt(i+1);
        c1=(1+r)*k1+(1-tau)*w*n1-k2;
    end
    c0=(1+r)*k0+(1-tau)*w*n0-k1;
    y(1)=uc(c0,1-n0)/beta-uc(c1,1-n1)*(1+r);
    y(2)=(1-tau)*w*uc(c0,1-n0)-un(c0,1-n0);
end
