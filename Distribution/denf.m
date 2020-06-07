clear
global tau r w b na pp a beta e i
global sigma vold
tol=0.001;          % stopping criterion for final solution of K 
tol1=1e-3;          % stopping criterion for golden section search 
tolg=0.00001;       % distribution function 
neg=-1e10;
eps=0.05;

ni=10000;       % number of individuals 

% calibration of parameters 
alpha=0.36;
beta=0.995;
delta=0.005;
sigma=2;
r=0.005;
w=5;
tau=0.02;
rep=0.25;
b=rep*(1-tau)*w;
pp=[0.9565 0.0435;0.5 0.5];
% stationary employment/unemployment 
% P*P*P*...
pp1=[0.9200 0.0800];
amin1=-2;                 % asset grid 
amax1=3000;
na=201;
a=linspace(amin1,amax1,na);

nk=3*na;            % asset grid for distribution 
ag=linspace(amin1,amax1,nk);
gk=zeros(nk,2);     % distribution 

% initialization of the value function: 
v=zeros(na,2);
for i=1:na;
    v(i,1)=u((1-tau)*r*a(i)+(1-tau)*w);
end
for i=1:na;
    v(i,2)=u((1-tau)*r*a(i)+b);
end
v=v./(1-beta);               % agents consume their income 
copt=zeros(na,2);           % optimal consumption 
aopt=zeros(na,2);           % optimal next-period assets 

nit=1000;                % number of maximum iterations over value function 
ngk=0;               % initial number of iterations over distribution function 
crit1=1;
kritg=1;
nq=50;

psi=0.95;

kritw=1;
nn0=pp1(1);
kk0=(alpha/(1/beta-1+delta)).^(1/(1-alpha))*nn0;
w0=(1-alpha)*kk0^alpha*nn0.^(-alpha);
b=w0*rep;
kkbar=kk0;
kk1=kk0-1;
kbarq=zeros(nq,1);      % convergence of kk0 
k1barq=zeros(nq,1);      % convergence of kk1 
kritwq=zeros(nq,2);     % convergence of value function v 
crit=1;
q=0;
while q<nq && (abs((kk1-kk0)/kk0)>tol); 
    q=q+1;
    if ngk<25000;
        ngk=ngk+500;
    end

    kt=zeros(ngk+1,1);
    if q==10; nit=nit*2; end
    if q==40; nit=nit*2; end
    w=(1-alpha)*kk0.^alpha*nn0.^(-alpha);
    r=alpha*kk0.^(alpha-1)*nn0.^(1-alpha)-delta;
    kbarq(q)=kk0;
    crit=1;
    while j<50 && crit>tol;
        j=j+1;
        vold=v;
        for e=1:2  % e=1 employed, e=2 unemployed                  
            for i=1:na; % iteration over asset grid a in period t 
                v0=neg;
                ax=a(1); bx=a(1); cx=a(na);
                for lt=1:na; % iteration over a' in period t*1 
                    if e==1;
                        c=(1+(1-tau)*r)*a(i)+(1-tau)*w-a(lt);
                    else
                        c=(1+(1-tau)*r)*a(i)+b-a(lt);
                    end
                    if c>0;
                        v1=bellman(a(i),a(lt),e);
                        if v1>v0;
                            v(i,e)=v1;
                            if lt==1;
                                ax=a(1); bx=a(1); cx=a(2);
                            elseif lt==na;
                                ax=a(na-1); bx=a(na); cx=a(na);
                            else
                                ax=a(lt-1); bx=a(lt); cx=a(lt+1);
                            end
                            v0=v1;
                        else
                            break   % concavity of value function and stop
                        end
                    else
                        break
                    end
                end   % lt=1,..,na 
                if ax==bx;  % boundary optimum, ax=bx=a(1)  
                    bx=ax+eps*(a(2)-a(1));
                    if value(bx,e)<value(ax,e);
                        aopt(i,e)=a(1);
                        v(i,e)=bellman(a(i),a(1),e);
                    else
                        %aopt(i,e)=golden(@(x)value1(x),ax,bx,cx,tol1);
                        aopt(i,e)=fminbnd(@(x)value1(x),ax,cx);
                        v(i,e)=bellman(a(i),aopt(i,e),e);
                    end
                elseif bx==cx;  % boundary optimum, bx=cx=a(n) 
                    bx=cx-eps*(a(na)-a(na-1));
                    if value(bx,e)<value(cx,e);
                        aopt(i,e)=a(na);
                    else
                        %aopt(i,e)=golden(@(x)value1(x),ax,bx,cx,tol1);
                        aopt(i,e)=fminbnd(@(x)value1(x),ax,cx);
                    end
                else
                    %aopt(i,e)=golden(@(x)value1(x),ax,bx,cx,tol1);
                    aopt(i,e)=fminbnd(@(x)value1(x),ax,cx);
                end
                v(i,e)=bellman(a(i),aopt(i,e),e);
            end;   % i=1,..na 
        end;   % e=1,2 
        crit=mean(abs(vold-v));
    end;   % j=1,..nit 
    crit
    %break %test
%     kritwq(q,.)=crit';
%     save vden=v,aoptden=aopt;
% 
%     copt(.,1)=(1+(1-tau)*r)*a+(1-tau)*w-aopt(.,1);
%     copt(.,2)=(1+(1-tau)*r)*a+b-aopt(.,2);
%     
%     q1=0;
%     kritg=1;
%     % initialization of the distribution functions 
%     %
%     gk=ones(nk,2)/nk;
%     gk=gk.*pp1';
%     
%     kconv=zeros(ngk+1,1);
%     gk=zeros(nk,2);
%     gk(nk0,1)=pp1(1);
%     gk(nk0,2)=pp1(2);
%     "computation of invariant distribution of wealth..";
%     do until (q1>ngk);
%         q1=q1+1; 
%         gk0=gk;
%         gk=zeros(nk,2);
%         lt=0;
%         do until lt==2;
%             lt=lt+1;
%             i=0;
%             do until i==nk;
%                 i=i+1;
%                 k0=ag(i);
%                 if k0<=amin1;
%                     k1=aopt(1,lt);
%                 elseif k0>=amax1;
%                     k1=aopt(na,lt);
%                 else;
%                     k1=lininter(a,aopt(.,lt),k0);
%                 endif;
%                 if k1<=amin1;
%                     gk(1,1)=gk(1,1)+gk0(i,lt)*pp(lt,1);
%                     gk(1,2)=gk(1,2)+gk0(i,lt)*pp(lt,2);
%                 elseif k1>=amax1;
%                     gk(nk,1)=gk(nk,1)+gk0(i,lt)*pp(lt,1);
%                     gk(nk,2)=gk(nk,2)+gk0(i,lt)*pp(lt,2);
%                 elseif (k1>amin1) and (k1<amax1);
%                     j=sumc(ag.<=k1)+1;
%                     n0=(k1-ag(j-1))/(ag(j)-ag(j-1));
%                     gk(j,1)=gk(j,1)+n0*gk0(i,lt)*pp(lt,1);
%                     gk(j,2)=gk(j,2)+n0*gk0(i,lt)*pp(lt,2);
%                     gk(j-1,1)=gk(j-1,1)+(1-n0)*gk0(i,lt)*pp(lt,1);
%                     gk(j-1,2)=gk(j-1,2)+(1-n0)*gk0(i,lt)*pp(lt,2);
%                 endif;
%             endo;
%         endo;
% 
% 
% 
%         gk=gk/sumc(sumc(gk));
%         kk1=(gk(.,1)+gk(.,2))'*ag;
%         kconv(q1)=kk1;
%         kt(q1)=kk1;
% 
%         nround=ngk/100;
%         if round(q1/nround)==q1/nround;
%             kk1=(gk(.,1)+gk(.,2))'*ag;
%             "iteration q: " q;
%             "time elapsed: " etstr(hsec-h0);
%             "iteration~capital stock: " q1~kk1;
%             "kbarq~kbarq1: ";
%             kbarq(1:q)~k1barq(1:q);
%             qt=q1/nround;
%         endif;
% 
%         kritg=sumc(abs(gk0-gk));
%     
end % q1=1,.., invariant distribution 
