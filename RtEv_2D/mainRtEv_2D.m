%% 数据矩阵生成

clear;

%%数据采集长度：Nd>=(1+n_p)*n_x*[(1+n_p)*n_u+1]-1=29
Nd=29;

[X0_hat,Xp_hat,X1_hat,U,Up,G_hat,n_psi,n_u,n_p,W_hat,A0,A1,B0,B1,C0,C1,D0,D1,p]=getDataMatirxRtEv_2D(Nd);

A_hat_0=[A0,zeros(1,2);C0,eye(2)];
A_hat_1=[A1,zeros(1,2);C1,zeros(2,2)];
B_hat_0=[B0;D0];
B_hat_1=[B1;D1];

[M11,M12,M21,M22]=LFT(n_psi,n_p,Nd);
%%验证是否能还原系统矩阵
AB_hat=(X1_hat-W_hat)*pinv(G_hat);
%%
%%调度变量顶点
p_1=-1;
p_2=1;

n_1_hat=n_p*(4*n_psi+Nd);
[U1,U2]=getUMatirx(p_1,p_2,n_psi,Nd,n_1_hat);

%%可调参数
ep_1=0.0001;
si_1=4;
beta_1=0.5;


%%便于检查与计算
delta=2.51;
Delta=delta*sqrt(Nd)*eye(n_psi);
XX1_hat=blkdiag(X1_hat,kron(eye(n_p),X1_hat));
WW_hat=ep_1*(Delta*Delta');
Ep_1=blkdiag(ep_1*eye(Nd),zeros(Nd*n_p));
M=[M11,M12;eye(n_1_hat),zeros(n_1_hat,n_1_hat);M21,M22];
p_ch=p_1;
[I1,I2,I3]=getIMatirx(n_psi,p_ch,Nd);

%%便于检查与计算

Ze_1=zeros((1+n_p)*n_psi,(1+n_p)*n_psi);
Ze_2=zeros((1+n_p)*n_psi,(1+n_p)*Nd);
Ze_3=zeros(2*n_1_hat,(1+n_p)*n_psi);
Ze_4=zeros(2*n_1_hat,(1+n_p)*Nd);


%%定义决策变量
Xi=sdpvar(2*n_1_hat,2*n_1_hat);
P=sdpvar(n_psi,n_psi);
FQ=sdpvar(Nd*(1+n_p),n_psi*(1+n_p),'full');
Z_0=sdpvar(n_u,n_psi);
Z_bar=sdpvar(n_u,n_psi*n_p);

ss1=sdpvar(1,1);
ss2=sdpvar(1,1);
% ss3=sdpvar(1,1);
% ss4=sdpvar(1,1);

%%定义中间变量
 Y=blkdiag(P,zeros(n_psi*n_p));
 Y_bar=blkdiag(P-WW_hat,zeros(n_psi*n_p));
%%LMI约束

LMI_1=[M'*[ Xi,     Ze_3,        Ze_3,      Ze_3,      Ze_3,          Ze_4;
 Ze_3',    -Y,         Ze_1,    -(XX1_hat*FQ)',   -Y,           -FQ';
 Ze_3',   Ze_1',     -si_1*Y,     -Y,        Ze_1,          Ze_2;
 Ze_3',  -(XX1_hat*FQ)     -Y,       -Y_bar,     Ze_1,          Ze_2;
 Ze_3',    -Y,         Ze_1',      Ze_1',   -(1/beta_1)*Y,  Ze_2;
 Ze_4',   -FQ,         Ze_2',      Ze_2',    Ze_2',        -Ep_1; ]*M<=0];

LMI_2=[U1'*Xi*U1>=0];
LMI_3=[U2'*Xi*U2>=0];

LMI_6=[Xi(n_1_hat+1:2*n_1_hat,n_1_hat+1:2*n_1_hat)<=0];

%%等式约束

%EqualCon1=[[P,zeros(n_psi,n_psi*n_p),zeros(n_psi,n_psi*n_p*n_p)]*I3==X0_hat*(I2'*FQ*I1)];

EqualCon1=[(X0_hat*(I2'*FQ*I1)-ss1<=[P,zeros(n_psi,n_psi*n_p),zeros(n_psi,n_psi*n_p*n_p)]*I3)&([P,zeros(n_psi,n_psi*n_p),zeros(n_psi,n_psi*n_p*n_p)]*I3<=X0_hat*(I2'*FQ*I1)+ss1)];


%%EqualCon2=[[zeros(n_psi*n_p,n_psi),kron(eye(n_p),P),zeros(n_psi*n_p,n_psi*n_p*n_p)]*I3==Xp_hat*(I2'*FQ*I1)];



EqualCon2=[(Xp_hat*(I2'*FQ*I1)-ss2<=[zeros(n_psi*n_p,n_psi),kron(eye(n_p),P),...
  zeros(n_psi*n_p,n_psi*n_p*n_p)]*I3)&([zeros(n_psi*n_p,n_psi),kron(eye(n_p),P),...
  zeros(n_psi*n_p,n_psi*n_p*n_p)]*I3<=Xp_hat*(I2'*FQ*I1)+ss2)];

EqualCon3=[[Z_0,Z_bar,zeros(n_u,n_psi*n_p*n_p)]*I3==U*(I2'*FQ*I1)];
%%EqualCon3=[([Z_0,Z_bar,zeros(n_u,n_psi*n_p*n_p)]*I3-ss3<=U*(I2'*FQ*I1))&([Z_0,Z_bar,zeros(n_u,n_psi*n_p*n_p)]*I3<=U*(I2'*FQ*I1)+ss3)];

EqualCon4=[[zeros(n_p*n_u,n_psi),kron(eye(n_p),Z_0),kron(eye(n_p),Z_bar)]*I3==Up*(I2'*FQ*I1)];
%EqualCon4=[([zeros(n_p*n_u,n_psi),kron(eye(n_p),Z_0),kron(eye(n_p),Z_bar)]*I3-ss4<=Up*(I2'*FQ*I1))&([zeros(n_p*n_u,n_psi),kron(eye(n_p),Z_0),kron(eye(n_p),Z_bar)]*I3<=Up*(I2'*FQ*I1)+ss4)];

%%SDP求解

Con=[LMI_1,LMI_2,LMI_3,LMI_6,EqualCon1,EqualCon2,EqualCon3,EqualCon4,trace(P)>=0.1,trace(P)<=10];
% Con=[LMI_1,LMI_2,LMI_3,LMI_4,LMI_5,LMI_6,EqualCon1,EqualCon2,EqualCon3,EqualCon4];

Objective = 0.5*norm(ss1, 2)^2+0.1*norm(ss2, 2)^2;


sol=optimize(Con,Objective);
if sol.problem == 0
    disp('求解成功');
else
    disp('运行出错');
    yalmiperror(sol.problem);
end

FQv =value(FQ);
Pv  =value(P);


    
%%求解

K0=value(Z_0)*pinv(Pv);
K_bar=value(Z_bar)*pinv(kron(eye(n_p),Pv));
K1=K_bar;

Xiv=value(Xi);
Yv=value(Y);
Y_barv=value(Y_bar);



LMI1_J = M'*[ Xiv,     Ze_3,        Ze_3,      Ze_3,      Ze_3,          Ze_4;
 Ze_3',    -Yv,         Ze_1,    -(XX1_hat*FQv)',   -Yv,           -FQv';
 Ze_3',   Ze_1',     -si_1*Yv,     -Yv,        Ze_1,          Ze_2;
 Ze_3',  -(XX1_hat*FQv)     -Yv,       -Y_barv,     Ze_1,          Ze_2;
 Ze_3',    -Yv,         Ze_1',      Ze_1',   -(1/beta_1)*Yv,  Ze_2;
 Ze_4',   -FQv,         Ze_2',      Ze_2',    Ze_2',        -Ep_1; ]*M; 
eigenvalues = eig(LMI1_J);
if all(eigenvalues <= 0)
    disp('矩阵是负定矩阵');
else
    disp('矩阵不是负定矩阵');
end


[M11_ev,M12_ev,M21_ev,M22_ev]=LFT_ev(n_psi,n_p,Nd);
n_2_hat=n_p*(2*n_psi+Nd);
[U1_ev,U2_ev]=getUMatirx_ev(p_1,p_2,n_psi,Nd,n_2_hat);

Ze_1_ev=zeros((1+n_p)*n_psi,(1+n_p)*n_psi);
Ze_2_ev=zeros((1+n_p)*n_psi,(1+n_p)*Nd);
Ze_3_ev=zeros(2*n_2_hat,(1+n_p)*n_psi);
Ze_4_ev=zeros(2*n_2_hat,(1+n_p)*Nd);

%%可调参数
mu=15;
ep_2=0.001;


beta_2=beta_1/2;

%%便于检查与计算
WW_ev=ep_2*(Delta*Delta');
Ep_2=blkdiag(ep_2*eye(Nd),zeros(Nd*n_p));
M_ev=[M11_ev,M12_ev;eye(n_2_hat),zeros(n_2_hat,n_2_hat);M21_ev,M22_ev];


Xi_ev=sdpvar(2*n_2_hat,2*n_2_hat);
Psi_1=sdpvar(n_psi,n_psi);
Psi_2=sdpvar(n_psi,n_psi);


Omega_11=-mu*beta_2*Pv+Pv*Psi_2*Pv;
Omega_22=-Pv*Psi_1*Pv+mu^2*WW_ev;

Y_ev=blkdiag(Omega_11,zeros(n_psi*n_p));
Y_bar_ev=blkdiag(Omega_22,zeros(n_psi*n_p));

LMI_1_ev=[M_ev'*[ Xi_ev,       Ze_3_ev,        Ze_3_ev,             Ze_4_ev;
                  Ze_3_ev',     Y_ev,         mu*(XX1_hat*FQv)',         FQv';
                  Ze_3_ev',   mu*(XX1_hat*FQv),    Y_bar_ev,           Ze_2_ev;
                  Ze_4_ev',     FQv,           Ze_2_ev',           -Ep_2; ]*M_ev<=0];


LMI_2_ev=[U1_ev'*Xi_ev*U1_ev>=0];
LMI_3_ev=[U2_ev'*Xi_ev*U2_ev>=0];


LMI_6_ev=[Xi_ev(n_2_hat+1:2*n_2_hat,n_2_hat+1:2*n_2_hat)<=0];

% LMI_7_ev=[Psi_1>=0];
Con_ev=[LMI_1_ev,LMI_2_ev,LMI_3_ev,LMI_6_ev,Psi_1>=0;Psi_2>=0];

sol_ev=optimize(Con_ev,[]);
if sol_ev.problem == 0
    disp('事件触发——求解成功');
else
    disp('事件触发——运行出错');
    yalmiperror(sol_ev.problem);
end

clear tmax u_ev y_ev x_ev chi_ev psi_ev u y x chi psi r y1 r1 k time e v condition seqiE intiE y2 r2
%%可调参数
nu=500;
ki=1;
s=1;

x{1}=1;
chi{1}=[1;1];


x_ev{1}=1;
chi_ev{1}=[1;1];
psi_ev{1}=[x_ev{1};chi_ev{1}];


dt=0.01;
TT=30;
tmax=TT/dt;
R=2.5;
psi{1}=[x{1};chi{1}];
[r] = getRef_RtEv_circle(tmax,R,pi);

for k=1:1:(tmax-1)
    %%Time
    time(k)=(k-1)*dt;
    %%p_v{k}
    p_v{k}=-0.8+k/((tmax-1)*4);
    %%E*d{k}
    d=0.1*rand(1,tmax);
    %%u{k}
    u{k}= (K0+p_v{k}(1,1)*K1)*psi{k};
    %%y{k}
    y{k} = (C0+p_v{k}(1,1)*C1)*x{k}+(D0+p_v{k}(1,1)*D1)*u{k};
    %%x{k}
    x{k+1} = (A0+p_v{k}(1,1)*A1)*x{k}+(B0+p_v{k}(1,1)*B1)*u{k}+d(:,k);
    %%chi{k}
    chi{k+1}=chi{k}+(y{k}-r{k});

    psi{k+1}=[x{k+1};chi{k+1}];


e{k}=psi_ev{ki}-psi_ev{k};
v{k}=(B_hat_0+p_v{k}(1,1)*B_hat_1)*(K0+p_v{k}(1,1)*K1)*e{k}+(B_hat_0+p_v{k}(1,1)*B_hat_1)*((K0+p_v{ki}(1,1)*K1)-(K0+p_v{k}(1,1)*K1))*psi_ev{ki};

condition{k}=v{k}'*value(Psi_1)*v{k}-psi_ev{k}'*value(Psi_2)*psi_ev{k}-nu;
    if condition{k}>=0
        seqiE(s)=ki;
        ki=k;
        intiE(s)=ki-seqiE(s); 
        s=s+1;
    end

   %%事件触发状态变量
    %%u{k}
    u_ev{k} = (K0+p_v{ki}(1,1)*K1)*psi_ev{k};
    %%y{k}
    y_ev{k} = (C0+p_v{k}(1,1)*C1)*x_ev{k}+(D0+p_v{k}(1,1)*D1)*u_ev{k};
    %%x{k}
    x_ev{k+1} = (A0+p_v{k}(1,1)*A1)*x_ev{k}+(B0+p_v{k}(1,1)*B1)*u_ev{k}+d(:,k);
    %%chi{k}
    chi_ev{k+1}=chi_ev{k}+(y_ev{k}-r{k});
    % psi_ev{k}=[x_ev{k};chi_ev{k}];
    psi_ev{k+1}=[x_ev{k+1};chi_ev{k+1}];




psi1_ev(k)=psi_ev{k}(1);


x1(k)=psi{k}(1);
x2(k)=psi{k}(2);
chi1(k)=psi{k}(3);

y1(k)=y{k}(1);
y2(k)=y{k}(2);
y1_ev(k)=y_ev{k}(1);
y2_ev(k)=y_ev{k}(2);
r1(k)=r{k}(1);
r2(k)=r{k}(2);

end

figure(1)

subplot(3,1,1)
 plot(time,y1,'b','Linewidth',1);
  hold on 
  plot(time,r1,'r--','Linewidth',1);
  xlim([0,tmax*dt]);
  ylim([-3,3]);
 xlabel('Time (Sec.)'); 
hh = legend('Output Signal under ETM $ y_{k} $ ','Reference Signal $ r_{k} $');
set(hh,'Interpreter','latex');


subplot(3,1,2)
  plot(time,y2,'b','Linewidth',1);
  hold on 
  plot(time,r2,'r--','Linewidth',1);
  xlim([0,tmax*dt]);
  ylim([-3,3]);
 xlabel('Time (Sec.)'); 
hh = legend('Output Signal under ETM $ y_{k} $ ','Reference Signal $ r_{k} $');
set(hh,'Interpreter','latex');

subplot(3,1,3)
stem(dt*seqiE,dt*intiE,'color','(0.0,0.45,0.74)','linewidth',0.5);
xlabel('Time (Sec.)'); 
xlim([0,tmax*dt]);
hh = legend('$\tau $');
set(hh,'Interpreter','latex','fontsize', 12);

figure(2)
plot(r1,r2);
hold on 
plot(y1,y2);
xlim([-3,3]);
ylim([-3,3]);


figure(3)
plot(r1,r2);
hold on 
plot(y1_ev,y2_ev);
xlim([-3,3]);
ylim([-3,3]);

%%

Psi_1v=value(Psi_1);
psi_2v=value(Psi_2);