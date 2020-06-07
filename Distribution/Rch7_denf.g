@ ------------------------------ Rch7_denf.g ----------------------------

Date: March 19, 2007

Author: Burkhard Heer

computation of the policy functions and the density function
in the heterogenous-agent neoclassical growth model
with value function iteration

monotonocity of a' in a

DENSITY FUNCTION, Algorithm 7.2.3.

concavity of value function v(.)

Linear interpolation between grid points

Maximization: golden section search method

input procedures:
golden -- routine for golden section search
lininter -- linear interpolation

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
graphset;
Macheps=1e-20;

h0=hsec;


@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @

/* numerical computation */
tolw=1e-7;      /* stopping criterion for value function */
tol=0.001;      /* stopping criterion for capital stock */
tol1=1e-7;      /* stopping criterion for golden section search */
tolg=1e-6;      /* distribution function */
neg=-1e10;
eps=0.05;


/* calibration of parameters */
alpha=0.36;
beta=0.995;
delta=0.005;
sigma=2;
r=0.005;
w=5;
tau=0.02;
rep=0.25;

pp=(0.9565~0.0435|0.5~0.5);
/* stationary employment / unemployment */
pp1=equivec1(pp);
nn0=pp1[1];

amin1=-2;                 /* asset grid */
amax1=3000;



na=201;
astep=(amax1-amin1)/(na-1);
a=seqa(amin1,astep,na);

save aden=a;
nk=3*na;            /* asset grid for distribution */
agstep=(amax1-amin1)/(nk-1);
ag=seqa(amin1,agstep,nk);
save agden=ag;
gk=zeros(nk,2);

nit=100;            /* number of maximum iterations of the value function */
ngk=0;          /* density function */    
crit1=1;
kritg=1;
kritw=1;
nq=100;

psi=0.95;

kk0=(alpha/(1/beta-1+delta))^(1/(1-alpha))*nn0;
kk1=kk0-10;
w0=(1-alpha)*kk0^alpha*nn0^(-alpha);
b=w0*rep;

nk0=sumc(ag.<=kk0);

/* initialization of the value function: */
v=u((1-tau)*r*a+(1-tau)*w)~u((1-tau)*r*a+b);
v=v/(1-beta);                /* agents consume their income */
copt=zeros(na,2);           /* optimal consumption */
aopt=zeros(na,2);           /* optimal next-period assets */

kbarq=zeros(nq,1);      /* convergence of kk0 */
k1barq=zeros(nq,1);      /* convergence of kk1 */
kritwq=zeros(nq,2);     /* convergence of value function v */
crit=1;

q=0;
do until q==nq or (abs((kk1-kk0)/kk0)<tol and q>50); 
    q=q+1;

    if ngk<25000;
        ngk=ngk+500;
    endif;


    kt=zeros(ngk+1,1);
    if q==10; nit=nit*2; endif;
    if q==40; nit=nit*2; endif;
    w=(1-alpha)*kk0^alpha*nn0^(-alpha);
    r=alpha*kk0^(alpha-1)*nn0^(1-alpha)-delta;
    kbarq[q]=kk0;

    crit=1;
    j=0;
    do until (j>50 and crit<tol) or (j==nit);
        j=j+1;
        vold=v;

        
        "iteration j~q: " j~q;
        "time elapsed: " etstr(hsec-h0);
        "error value function: " crit;
        "error k: " crit1;
        "error distribution function: " kritg;
        "kk0: " kk0;
        "kk1: " kk1;

@ ----------------------------------------------------------------

Step 2: Iteration of the value function

----------------------------------------------------------------- @

        e=0;    /* iteration over the employment status */
        do until e==2;  /* e=1 employed, e=2 unemployed */
            e=e+1;
            i=0;        /* iteration over asset grid a in period t */
            l0=0;
            do until i==na;
                i=i+1;
                l=l0;
                v0=neg;
                ax=a[1]; bx=a[1]; cx=a[na];
                do until l==na; /* iteration over a' in period t*1 */
                    l=l+1;
                    if e==1;
                        c=(1+(1-tau)*r)*a[i]+(1-tau)*w-a[l];
                    else;
                        c=(1+(1-tau)*r)*a[i]+b-a[l];
                    endif;
                    if c>0;
                        v1=bellman(a[i],a[l],e);
                        if v1>v0;
                            v[i,e]=v1;
                            if l==1;
                                ax=a[1]; bx=a[1]; cx=a[2];
                            elseif l==na;
                                ax=a[na-1]; bx=a[na]; cx=a[na];
                            else;
                                ax=a[l-1]; bx=a[l]; cx=a[l+1];
                            endif;
                            v0=v1;
                            l0=l-1;
                        else;
                            l=na;   /* concavity of value function */
                        endif;
                    else;
                        l=na;
                    endif;
                endo;   /* l=1,..,na */
                if ax==bx;  /* boundary optimum, ax=bx=a[1]  */
                    bx=ax+eps*(a[2]-a[1]);
                    if value(bx,e)<value(ax,e);
                        aopt[i,e]=a[1];
                        v[i,e]=bellman(a[i],a[1],e);
                    else;
                        aopt[i,e]=golden(&value1,ax,bx,cx,tol1);
                        v[i,e]=bellman(a[i],aopt[i,e],e);
                   endif;
                elseif bx==cx;  /* boundary optimum, bx=cx=a[n] */
                    bx=cx-eps*(a[na]-a[na-1]);
                    if value(bx,e)<value(cx,e);
                        aopt[i,e]=a[na];
                    else;
                        aopt[i,e]=golden(&value1,ax,bx,cx,tol1);
                    endif;
                else;
                    aopt[i,e]=golden(&value1,ax,bx,cx,tol1);
                endif;
                v[i,e]=bellman(a[i],aopt[i,e],e);
            endo;   /* i=1,..na */
        endo;   /* e=1,2 */
        crit=meanc(abs(vold-v));
    endo;   /* j=1,..nit */

    kritwq[q,.]=crit';
    save vden=v,aoptden=aopt;

    copt[.,1]=(1+(1-tau)*r)*a+(1-tau)*w-aopt[.,1];
    copt[.,2]=(1+(1-tau)*r)*a+b-aopt[.,2];


  @-------------------------------------------------------------------------------------

 iteration to find invariant distribution 

 --------------------------------------------------------------------------------------- @
    
    q1=0;
    kritg=1;
    /* initialization of the distribution functions */
    /*
    gk=ones(nk,2)/nk;
    gk=gk.*pp1';
    */
    kconv=zeros(ngk+1,1);
    gk=zeros(nk,2);
    gk[nk0,1]=pp1[1];
    gk[nk0,2]=pp1[2];
    "computation of invariant distribution of wealth..";
    do until (q1>ngk);
        q1=q1+1; 
        gk0=gk;
        gk=zeros(nk,2);
        l=0;
        do until l==2;
            l=l+1;
            i=0;
            do until i==nk;
                i=i+1;
                k0=ag[i];
                if k0<=amin1;
                    k1=aopt[1,l];
                elseif k0>=amax1;
                    k1=aopt[na,l];
                else;
                    k1=lininter(a,aopt[.,l],k0);
                endif;
                if k1<=amin1;
                    gk[1,1]=gk[1,1]+gk0[i,l]*pp[l,1];
                    gk[1,2]=gk[1,2]+gk0[i,l]*pp[l,2];
                elseif k1>=amax1;
                    gk[nk,1]=gk[nk,1]+gk0[i,l]*pp[l,1];
                    gk[nk,2]=gk[nk,2]+gk0[i,l]*pp[l,2];
                elseif (k1>amin1) and (k1<amax1);
                    j=sumc(ag.<=k1)+1;
                    n0=(k1-ag[j-1])/(ag[j]-ag[j-1]);
                    gk[j,1]=gk[j,1]+n0*gk0[i,l]*pp[l,1];
                    gk[j,2]=gk[j,2]+n0*gk0[i,l]*pp[l,2];
                    gk[j-1,1]=gk[j-1,1]+(1-n0)*gk0[i,l]*pp[l,1];
                    gk[j-1,2]=gk[j-1,2]+(1-n0)*gk0[i,l]*pp[l,2];
                endif;
            endo;
        endo;



        gk=gk/sumc(sumc(gk));
        kk1=(gk[.,1]+gk[.,2])'*ag;
        kconv[q1]=kk1;
        kt[q1]=kk1;

        nround=ngk/100;
        if round(q1/nround)==q1/nround;
            kk1=(gk[.,1]+gk[.,2])'*ag;
            "iteration q: " q;
            "time elapsed: " etstr(hsec-h0);
            "iteration~capital stock: " q1~kk1;
            "kbarq~kbarq1: ";
            kbarq[1:q]~k1barq[1:q];
            qt=q1/nround;
        endif;

        kritg=sumc(abs(gk0-gk));
    endo;   /* q1=1,.., invariant distribution */

    kk1=(gk[.,1]+gk[.,2])'*ag;
    k1barq[q]=kk1;
    "kk0: " kk0;
    "kk1: " kk1;
    "wage: " w;
    "interest rate: " r;
    "kbarq~k1barq";
    kbarq[1:q]~k1barq[1:q];

    save aden=a;
    save aoptden=aopt;
    save vden=v;
    save agden=ag;
    save gkden=gk;
    save k1den=k1barq;

    save kt;
    "kk0: " kk0;
/*
    tau0=w*nn0*rep*pp1[2]/(w*nn0*(1+rep*pp1[2])+r*kk0);
    "tau0: " tau0;
    tau=psi*tau+(1-psi)*tau0;
    "tau: " tau;
*/
    crit1=abs((kk1-kk0)/kk0);
    kk0=psi*kk0+(1-psi)*kk1;
    tau=b*pp1[2]/(w*nn0+r*kk0);
    "tau: " tau;
    "iteration: " q;
    "time elapsed: " etstr(hsec-h0);
    "error value function: " crit;
    "error distribution: " kritg;
    "error capital stock: " crit1;


    if q==10;
        title("Mean of the distribution function");
        xlabel("Number of iterations over distribution function F(.)");
        xy(seqa(1,1,ngk+1),kt);
        wait;
    endif;

    
    if q==nq;
        wait;
        title("invariant capital distribution");
        xy(ag,gk);
        wait;
        title("optimal savings");
        xy(a,aopt-a);
        wait;
    endif;
    
    if q/2==round(q/2); 
        cls;
    endif;
endo;   /* q=1,..,nq */




"iteration: " q;
"time elapsed: " etstr(hsec-h0);
"error value function: " crit;
"error distribution: " kritg;
"error capital stock: " crit1;
wait;

/* plotting the solution */

title("invariant capital distribution");
xy(ag,gk);
wait;
xlabel("asset level");
_plegstr="employed\000unemployed";
title("value function");
xy(a,v);
wait;
xy(a,v[.,1]);
xy(a,v[.,2]);

title("consumption");
xy(a,copt);
wait;
xy(a,copt[.,1]);
xy(a,copt[.,2]);

title("change in asset level");
xy(a,aopt-a);
wait;
xy(a,aopt[.,1]-a);
xy(a,aopt[.,2]-a);

title("next-period assets");
_plegstr="employed\000a";
xy(a,aopt[.,1]~a);

_plegstr="unemployed\000a";
xy(a,aopt[.,2]~a);


proc u(x);
   retp(x^(1-sigma)/(1-sigma));
endp;

@  ----------------------------  procedures -----------


u(x) -- utility function

value(a,e) -- returns the value of the value function for asset
                a and employment status e

value1(x) -- given a=a[i] and epsilon=e, returns the
             value of the bellman equation for a'=x

bellman -- value for the right-hand side of the Bellman equation

------------------------------------------------------- @

proc u(x);
   retp(x^(1-sigma)/(1-sigma));
endp;

proc value(x,y);
    retp(lininter(a,vold[.,y],x));
endp;

proc value1(x);
    retp(bellman(a[i],x,e));
endp;


proc bellman(a0,a1,y);
   local c;
   if y==1;
      c=(1+(1-tau)*r)*a0+(1-tau)*w-a1;
   else;
      c=(1+(1-tau)*r)*a0+b-a1;
   endif;
   if c<0;
      retp(neg);
   endif;
   if a1>=a[na];
      retp(a1^2*neg);
   endif;
   retp(u(c)+beta*(pp[y,1]*value(a1,1)+pp[y,2]*value(a1,2)));
endp;

proc lininter(xd,yd,x);
  local j;
  j=sumc(xd.<=x');
  retp(yd[j]+(yd[j+1]-yd[j]).*(x-xd[j])./(xd[j+1]-xd[j]));
endp;

proc golden(&f,ay,by,cy,tol);
    local f:proc,x0,x1,x2,x3,xmin,r1,r2,f1,f2;
    r1=0.61803399; r2=1-r1;
    x0=ay;
    x3=cy;
    if abs(cy-by)<=abs(by-ay);
        x1=by; x2=by+r2*(cy-by);
    else;
        x2=by; x1=by-r2*(by-ay);
    endif;
    f1=-f(x1);
    f2=-f(x2);
    do until abs(x3-x0)<=tol*(abs(x1)+abs(x2));
        if f2<f1;
            x0=x1;
            x1=x2;
            x2=r1*x1+r2*x3;
            f1=f2;
            f2=-f(x2);
        else;
            x3=x2;
            x2=x1;
            x1=r1*x2+r2*x0;
            f2=f1;
            f1=-f(x1);
        endif;
    endo;
    if f1<=f2;
        xmin=x1;
        else;
        xmin=x2;
    endif;
    retp(xmin);
endp;

proc equivec1(p);
    local n,x;
    n=rows(p);
    p=diagrv(p,diag(p)-ones(n,1));
    p=p[.,1:n-1]~ones(n,1);
    x=zeros(n-1,1)|1;
    retp((x'*inv(p))');
endp;
