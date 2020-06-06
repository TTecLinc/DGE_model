function y=rfold(x,i)
global t r pen aopt w tau tr
global s gam psi beta
    k0=x;
    k1=aopt(i+t+1);
    if i==tr-1;
        k2=0;
    else
        k2=aopt(i+t+2);
    end
    c0=(1+r)*k0+pen-k1;
    c1=(1+r)*k1+pen-k2;
    y=uc(c0,1)/beta-uc(c1,1)*(1+r);
end

