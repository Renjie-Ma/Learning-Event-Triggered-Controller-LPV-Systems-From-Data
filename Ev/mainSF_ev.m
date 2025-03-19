%% 数据矩阵生成

clear;
Nd=23;
[X0,Xp,X1,U,Up,G,n_x,n_u,n_p,W,A0,A1,A2,B0,B1,B2]=getDataMatrix2(Nd);

%%验证是否能还原系统矩阵
AB=(X1-W)*pinv(G);

%% 控制器求解

%%调度变量顶点
p_1=[-1;1];
p_2=[-1;-1];
p_3=[1;-1];
p_4=[1;1];

n_1=n_p*(4*n_x+Nd);
[U1,U2,U3,U4]=getUMatirx(p_1,p_2,p_3,p_4,n_x,Nd,n_1);
[M11,M12,M21,M22]=LFT(n_x,n_p,Nd);

%%可调参数
ep_1=0.01;
si_1=4;
beta_1=0.2;


%%扰动上界
delta=0.1;
Delta=delta*sqrt(Nd)*eye(n_x);
XX1=blkdiag(X1,kron(eye(n_p),X1));
WW=ep_1*(Delta*Delta');
Ep_1=blkdiag(ep_1*eye(Nd),zeros(Nd*n_p));
M=[M11,M12;eye(n_1),zeros(n_1,n_1/2);M21,M22];
p_ch=p_2;
[I1,I2,I3]=getIMatirx(n_x,p_ch,Nd);

%%便于检查与计算

Ze_1=zeros((1+n_p)*n_x,(1+n_p)*n_x);
Ze_2=zeros((1+n_p)*n_x,(1+n_p)*Nd);
Ze_3=zeros(2*n_1,(1+n_p)*n_x);
Ze_4=zeros(2*n_1,(1+n_p)*Nd);


%%定义决策变量
Xi=sdpvar(2*n_1,2*n_1);
P=sdpvar(n_x,n_x);
FQ=sdpvar(Nd*(1+n_p),n_x*(1+n_p),'full');
Z_0=sdpvar(n_u,n_x);
Z_bar=sdpvar(n_u,n_x*n_p);

%%定义中间变量
 Y=blkdiag(P,zeros(n_x*n_p));
 Y_bar=blkdiag(P-WW,zeros(n_x*n_p));

%%LMI约束

LMI_1=[M'*[ Xi,     Ze_3,        Ze_3,      Ze_3,      Ze_3,          Ze_4;
 Ze_3',    -Y,         Ze_1,    -(XX1*FQ)',   -Y,           -FQ';
 Ze_3',   Ze_1',     -si_1*Y,     -Y,        Ze_1,          Ze_2;
 Ze_3',  -(XX1*FQ)     -Y,       -Y_bar,     Ze_1,          Ze_2;
 Ze_3',    -Y,         Ze_1',      Ze_1',   -(1/beta_1)*Y,  Ze_2;
 Ze_4',   -FQ,         Ze_2',      Ze_2',    Ze_2',        -Ep_1; ]*M<=0];

LMI_2=[U1'*Xi*U1>=0];
LMI_3=[U2'*Xi*U2>=0];
LMI_4=[U3'*Xi*U3>=0];
LMI_5=[U4'*Xi*U4>=0];

LMI_6=[Xi(n_1+1:2*n_1,n_1+1:2*n_1)<=0];

%%等式约束

EqualCon1=[[P,zeros(n_x,n_x*n_p),zeros(n_x,n_x*n_p*n_p)]*I3==X0*I2'*FQ*I1];

EqualCon2=[[zeros(n_x*n_p,n_x),kron(eye(n_p),P),zeros(n_x*n_p,n_x*n_p*n_p)]*I3==Xp*I2'*FQ*I1];

EqualCon3=[[Z_0,Z_bar,zeros(n_u,n_x*n_p*n_p)]*I3==U*I2'*FQ*I1];

EqualCon4=[[zeros(n_p*n_u,n_x),kron(eye(n_p),Z_0),kron(eye(n_p),Z_bar)]*I3==Up*I2'*FQ*I1];
%%SDP求解

Con=[LMI_1,LMI_2,LMI_3,LMI_4,LMI_5,LMI_6,EqualCon1,EqualCon2,EqualCon3,EqualCon4,trace(P)>=0.1,trace(P)<=10];
sol=optimize(Con,[]);


FQv =value(FQ);
Pv  =value(P);

if sol.problem == 0
    disp('求解成功');
else
    disp('运行出错');
    yalmiperror(sol.problem);
end


%%求解
K0=value(Z_0)*pinv(Pv);
K_bar=value(Z_bar)*pinv(kron(eye(n_p),Pv));
K1=K_bar(:,1:2);
K2=K_bar(:,3:4);

%% 事件触发增益求解
[M11_ev,M12_ev,M21_ev,M22_ev]=LFT_ev(n_x,n_p,Nd);
n_2=n_p*(2*n_x+Nd);
[U1_ev,U2_ev,U3_ev,U4_ev]=getUMatirx_ev(p_1,p_2,p_3,p_4,n_x,Nd,n_2);

Ze_1_ev=zeros((1+n_p)*n_x,(1+n_p)*n_x);
Ze_2_ev=zeros((1+n_p)*n_x,(1+n_p)*Nd);
Ze_3_ev=zeros(2*n_2,(1+n_p)*n_x);
Ze_4_ev=zeros(2*n_2,(1+n_p)*Nd);

%%可调参数
mu=40;
ep_2=0.001;


beta_2=beta_1/2;

%%便于检查与计算
WW_ev=ep_2*(Delta*Delta');
Ep_2=blkdiag(ep_2*eye(Nd),zeros(Nd*n_p));
M_ev=[M11_ev,M12_ev;eye(n_2),zeros(n_2,n_2/2);M21_ev,M22_ev];


Xi_ev=sdpvar(2*n_2,2*n_2);
Psi_1=sdpvar(n_x,n_x);
Psi_2=sdpvar(n_x,n_x);


Omega_11=-mu*beta_2*Pv+Pv*Psi_2*Pv;
Omega_22=-Pv*Psi_1*Pv+mu^2*WW_ev;

Y_ev=blkdiag(Omega_11,zeros(n_x*n_p));
Y_bar_ev=blkdiag(Omega_22,zeros(n_x*n_p));

LMI_1_ev=[M_ev'*[ Xi_ev,       Ze_3_ev,        Ze_3_ev,             Ze_4_ev;
                  Ze_3_ev',     Y_ev,         mu*(XX1*FQv)',         FQv';
                  Ze_3_ev',   mu*(XX1*FQv),    Y_bar_ev,           Ze_2_ev;
                  Ze_4_ev',     FQv,           Ze_2_ev',           -Ep_2; ]*M_ev<=0];


LMI_2_ev=[U1_ev'*Xi_ev*U1_ev>=0];
LMI_3_ev=[U2_ev'*Xi_ev*U2_ev>=0];
LMI_4_ev=[U3_ev'*Xi_ev*U3_ev>=0];
LMI_5_ev=[U4_ev'*Xi_ev*U4_ev>=0];


LMI_6_ev=[Xi_ev(n_2+1:2*n_2,n_2+1:2*n_2)<=0];

Con_ev=[LMI_1_ev,LMI_2_ev,LMI_3_ev,LMI_4_ev,LMI_5_ev,LMI_6_ev,Psi_1>=0,Psi_2>=0];

sol_ev=optimize(Con_ev,[]);
if sol_ev.problem == 0
    disp('事件触发——求解成功');
else
    disp('事件触发——运行出错');
    yalmiperror(sol_ev.problem);
end


Psi_1v=value(Psi_1);
Psi_2v=value(Psi_2);


clear tmax x_ev x1_ev x2_ev x1 x2 x  k time e v condition seqiE intiE 
%%可调参数
nu=0.01;


x{1}=[2;-2];
x_ev{1}=[2;-2];

dt=0.01;
ki=1;
s=1;
TT=2;
tmax=TT/dt;
for k=1:1:tmax
    p_v{k}=[-0.6+k/(tmax*4);-0.6+k/(tmax*4)];
    w=0.1*rand(1,tmax);
    time(k)=(k-1)*dt;
    x{k+1}=(A0+p_v{k}(1,1)*A1+p_v{k}(2,1)*A2)*x{k}+(B0+p_v{k}(1,1)*B1+p_v{k}(2,1)*B2)*(K0+p_v{k}(1,1)*K1+p_v{k}(2,1)*K2)*x{k}+w(k);
    
    e{k}=x_ev{ki}-x_ev{k};
    v{k}=(B0+p_v{k}(1,1)*B1+p_v{k}(2,1)*B2)*(K0+p_v{k}(1,1)*K1+p_v{k}(2,1)*K2)*e{k}+(B0+p_v{k}(1,1)*B1+p_v{k}(2,1)*B2)*((K0+p_v{ki}(1,1)*K1+p_v{ki}(2,1)*K2)-(K0+p_v{k}(1,1)*K1+p_v{k}(2,1)*K2))*x_ev{ki};
    condition{k}=v{k}'*value(Psi_1)*v{k}-x_ev{k}'*value(Psi_2)*x_ev{k}-nu;
    if condition{k}>=0
        seqiE(s)=ki;
        ki=k;
        intiE(s)=ki-seqiE(s); 
        s=s+1;
    end


x_ev{k+1}=(A0+p_v{k}(1,1)*A1+p_v{k}(2,1)*A2)*x_ev{k}+(B0+p_v{k}(1,1)*B1+p_v{k}(2,1)*B2)*(K0+p_v{ki}(1,1)*K1+p_v{ki}(2,1)*K2)*x_ev{k}+w(k);

x1(k)=x{k}(1);
x2(k)=x{k}(2);

x1_ev(k)=x_ev{k}(1);


end

figure(1)  

  plot(time,x1,'b','Linewidth',1);
  xlim([0,dt*tmax]);
  ylim([-5,5]);
  xlabel('Time (Sec.)'); 
[up1,lo1] = envelope(x1,1000,'analytic');
hold on

plot(time,up1,'-','Linewidth',1);
hold on
plot(time,lo1,'--','Linewidth',1);
hh=legend('State $ x_{k} $','Signal Upper Bound','Signal Lower Bound');

set(hh,'Interpreter','latex');
hold off

figure(2)
  subplot(2,1,1)
  plot(time,x1_ev,'b','Linewidth',1);

  xlim([0,dt*tmax]);
  ylim([-5,5]);
 xlabel('Time (Sec.)'); 

[up2,lo2] = envelope(x1_ev,3000,'analytic');
hold on
plot(time,up2,'-','Linewidth',0.5);

hold on
plot(time,lo2,'--','Linewidth',0.5);

hh=legend('State under ETM $ x_{k} $','Signal Upper Bound','Signal Lower Bound');

set(hh,'Interpreter','latex');
hold off

subplot(2,1,2)
stem(dt*seqiE,dt*intiE,'color','(0.0,0.45,0.74)','linewidth',0.5);
xlabel('Time (Sec.)'); 
xlim([0,dt*tmax]);
hh = legend('$\tau $');
set(hh,'Interpreter','latex','fontsize', 12);
