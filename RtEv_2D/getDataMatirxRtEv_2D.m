%%函数功能：获取数据矩阵
%%输入参数：Nd--数据长度
%%输出参数：X0--状态矩阵；Xp--状态——调度矩阵；X1--状态矩阵；%%          
%           U--输入矩阵；Up--输入-调度矩阵；  W--扰动矩阵； G--[X0;Xp;U;Up]      
function [X0_hat,Xp_hat,X1_hat,U,Up,G_hat,n_psi,n_u,n_p,W_hat,A0,A1,B0,B1,C0,C1,D0,D1,p]=getDataMatirxRtEv_2D(Nd)
A0=rand(1,1);
A1=rand(1,1);
B0=rand(1,2);
B1=rand(1,2);
C0=rand(2,1);
C1=rand(2,1);
D0=rand(2,2);
D1=rand(2,2);
A_hat_0=[A0,zeros(1,2);C0,eye(2)];
A_hat_1=[A1,zeros(1,2);C1,zeros(2,2)];
B_hat_0=[B0;D0];
B_hat_1=[B1;D1];
%%矩阵维度定义
n_psi=size(A_hat_0,1);
n_u=size(B0,2);
n_p=1;
p1=-1+2*rand(1,Nd);              
p=p1;                                
psi_0=10*rand(n_psi,1);            
psi=zeros(n_psi,Nd+1);
u=10*rand(n_u,Nd);
w_hat=0.1*rand(n_psi,Nd);
Xp_hat=zeros(n_p*n_psi,Nd);
Up=zeros(n_p*n_u,Nd);
%%数据采集
psi(:,1)=(A_hat_0+A_hat_1*p(1,1))*psi_0+(B_hat_0+B_hat_1*p(1,1))*u(:,1)+w_hat(:,1);
for dl=1:Nd-1

    psi(:,dl+1)=(A_hat_0+A_hat_1*p(1,dl+1))*psi(:,dl)+(B_hat_0+B_hat_1*p(1,dl+1))*u(:,dl+1)+w_hat(:,dl+1);
end
Xp_hat(:,1)=kron(p(:,1),psi_0);
Up(:,1)=kron(p(:,1),u(:,1));
for k=2:Nd
    Xp_hat(:,k)=kron(p(:,k),psi(:,k-1));
    Up(:,k)=kron(p(:,k),u(:,k));
end

X0_hat=[psi_0,psi(:,1:Nd-1)];
X1_hat=psi(:,1:Nd);
U=u;
W_hat=w_hat;
G_hat=[X0_hat;Xp_hat;U;Up];
end

