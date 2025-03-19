%% 数据矩阵生成
clear;
%%数据采集长度：Nd>=(1+n_p)*n_x*[(1+n_p)*n_u+1]-1=23
Nd=23;

%%生成数据矩阵
[X0,Xp,X1,U,Up,G,n_x,n_u,n_p,W,A0,A1,A2,B0,B1,B2]=getDataMatrix(Nd);

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
ep_1=0.1;
si_1=2;
beta_1=0.2;

%%扰动上界
delta=0.1;
Delta=delta*sqrt(Nd)*eye(n_x);

WW=ep_1*(Delta*Delta');
XX1=blkdiag(X1,kron(eye(n_p),X1));
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

EqualCon1=[[P,zeros(n_x,n_x*n_p),zeros(n_x,n_x*n_p*n_p)]*I3==X0*(I2'*FQ*I1)];

EqualCon2=[[zeros(n_x*n_p,n_x),kron(eye(n_p),P),zeros(n_x*n_p,n_x*n_p*n_p)]*I3==Xp*(I2'*FQ*I1)];

EqualCon3=[[Z_0,Z_bar,zeros(n_u,n_x*n_p*n_p)]*I3==U*(I2'*FQ*I1)];

EqualCon4=[[zeros(n_p*n_u,n_x),kron(eye(n_p),Z_0),kron(eye(n_p),Z_bar)]*I3==Up*(I2'*FQ*I1)];

%%SDP求解

Con=[LMI_1,LMI_2,LMI_3,LMI_4,LMI_5,LMI_6,EqualCon1,EqualCon2,EqualCon3,EqualCon4,trace(P)>=0.1,trace(P)<=10];
sol=optimize(Con,[]);

%%提取求解出的值
FQv =value(FQ);
Pv  =value(P);

if sol.problem == 0
    disp('求解成功');
else
    disp('运行出错');
    yalmiperror(sol.problem);
end


%%反解出控制增益
K0=value(Z_0)*pinv(Pv);
K_bar=value(Z_bar)*pinv(kron(eye(n_p),Pv));

K1=K_bar(:,1:2);
K2=K_bar(:,3:4);


%% 画图
x{1}=[3;3];
dt=0.01;
TT=3;
tmax=TT/dt;
for k=1:1:(tmax-1)
    p_v{k}=[-1+k/((tmax-1)*2);-1+k/((tmax-1)*2)];
    w=0.1*rand(1,tmax);
    time(k)=(k-1)*dt;
x{k+1}=(A0+p_v{k}(1,1)*A1+p_v{k}(2,1)*A2)*x{k}+(B0+p_v{k}(1,1)*B1+p_v{k}(2,1)*B2)*(K0+p_v{k}(1,1)*K1+p_v{k}(2,1)*K2)*x{k}+w(k);

x1(k)=x{k}(1);
x2(k)=x{k}(2);


end

figure(1)
subplot(2,1,1)
  plot(time,x1,'b','Linewidth',1);
  xlim([0,dt*tmax]);
  ylim([-5,5]);
 xlabel('Time (Sec.)'); 
hh = legend('$ x_1 $');
set(hh,'Interpreter','latex');

subplot(2,1,2)
  plot(time,x2,'b','Linewidth',1);
  xlim([0,dt*tmax]);
  ylim([-5,5]);
 xlabel('Time (Sec.)'); 
hh = legend('$ x_2 $');
set(hh,'Interpreter','latex');
