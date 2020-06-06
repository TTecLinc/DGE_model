%% Peilin Yang
%% As Auerbach and Kotlikoff (1987) Auerbach and Kotlikoff (1987)
global t r pen aopt w tau tr
global s gam psi beta nopt copt
%% Parameter
beta=0.99;         % discount factor 
r=0.045;         % initial value of the interest rate 
s=2;            % coefficient of relative risk aversion 
alp=0.3;        % production elasticity of capital 
rep=0.3;        % replacement ratio 
del=0.1;          % rate of depreciation 
tr=20;          % retired 
t=40;           % working time 
tau=rep/(2+rep);    % income tax rate 
gam=2;            % disutility from working 

kmax=10;        % upper limit of capital grid 
kinit=0;
na=101;          % number of grid points on assets 
%a=seqa(0,kmax/(na-1),na);   % asset grid 
psi=0.001;      % parameter of utility function 

tol=0.001;       % percentage deviation of final solution 
tolk=0.01;       % percentage deviation of final solution for k_1 
tol1=1e-6;        % tolerance for golden section search 
neg=-1e10;        % initial value for value function 
nq1=30;

% Initialization
nbar=0.2;
kbar=(alp/(r+del)).^(1/(1-alp))*nbar;
kold=100;
nold=2;
aopt=zeros(t+tr,1);
copt=zeros(t+tr,1);   
nopt=0.3*ones(t,1);
total_time=50;
kseq=zeros(total_time,1);
% Guess the k60
k60q=zeros(nq1,1);
k60q(1)=0.15;
k60q(2)=0.2;
k1q=zeros(nq1,1);

krit=1;
krit0=1;
Time=1;
phi=0.8;
%% Caculate
while krit>tol && krit0>tol && Time<total_time
    krit
    krit=abs((kbar-kold)/kbar);
    krit0=abs((nbar-nold)/nbar);
    w=(1-alp)*kbar^alp*nbar^(-alp);
    r=alp*kbar^(alp-1)*nbar^(1-alp)-del;
    pen=rep*(1-tau)*w*nbar*3/2;
    kold=kbar;
    nold=nbar;
    kseq(Time)=kbar;
    Time=Time+1;
    for i=1:nq1
        if i>2          
           %% Gauss method to iterate the k60
            k60q(i)=k60q(i-1)-(k60q(i-1)-k60q(i-2))/(k1q(i-1)-k1q(i-2))*k1q(i-1);    
        end
        aopt(t+tr)=k60q(i);
        opt=optimoptions('fsolve','Display','none');
        for t_re=tr-1:-1:1
            [k_for res]=fsolve(@(x)rfold(x,t_re),aopt(t_re+t+1),opt);
            aopt(t_re+t)=k_for;
        end
        for t_re=t:-1:1
            [k_for res]=fsolve(@(x)rfyoung(x,t_re),[aopt(t_re+1) nopt(t_re)],opt);
            aopt(t_re)=k_for(1);
            nopt(t_re)=k_for(2);
        end
        if i>2
            if abs((k60q(i)-k60q(i-1))/k60q(i-1))<tol && aopt(1)<tol1
                break
            end
        end
        k1q(i)=aopt(1);
    end
    %%  Use extrapolation to stabilize the sequence
    knew=mean(aopt);
    nnew=mean(nopt)*2/3;
    kbar=phi*kold+(1-phi)*knew;
    nbar=phi*nold+(1-phi)*nnew;
end
plot(nopt)
figure();
plot(aopt)
figure();
plot(kseq(1:Time-1));
figure();
plot(k60q(1:i))