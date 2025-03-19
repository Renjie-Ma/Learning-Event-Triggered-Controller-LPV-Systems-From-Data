%%函数功能：获取数据矩阵
%%输入参数：Nd--数据长度
%%输出参数：X0--状态矩阵；Xp--状态——调度矩阵；X1--状态矩阵；%%          
%           U--输入矩阵；Up--输入-调度矩阵；  W--扰动矩阵； G--[X0;Xp;U;Up]
        
function [X0,Xp,X1,U,Up,G,n_x,n_u,n_p,W,A0,A1,A2,B0,B1,B2]=getDataMatrix2(Nd)
A0=[0.2485,-1.0355;0.8910,0.4065];
A1=[-0.0063,-0.0938;0,0.0188];
A2=[-0.0063,-0.0938;0,0.0188];
B0=[0.3190;-1.3080];
B1=[0.3;1.4];
B2=zeros(2,1);
%%矩阵维度定义
n_x=size(A0,1);
n_u=size(B0,2);
n_p=2;
p1=-1+2*rand(1,Nd);    
p2=-1+2*rand(1,Nd);                
p=[p1;p2];                                
x0=10*rand(n_x,1);            
x=zeros(n_x,Nd+1);
u=10*rand(n_u,Nd);
w=0.1*rand(1,Nd); 
Xp=zeros(n_p*n_x,Nd);
Up=zeros(n_p*n_u,Nd);
%%数据采集
x(:,1)=(A0+A1*p(1,1)+A2*p(2,1))*x0+(B0+B1*p(1,1)+B2*p(2,1))*u(1,1)+w(1,1);
for dl=1:Nd-1

    x(:,dl+1)=(A0+A1*p(1,dl+1)+A2*p(2,dl+1))*x(:,dl)+(B0+B1*p(1,dl+1)+B2*p(2,dl+1))*u(:,dl+1)+w(1,dl+1);
end
Xp(:,1)=kron(p(:,1),x0);
Up(:,1)=kron(p(:,1),u(:,1));
for k=2:Nd
    Xp(:,k)=kron(p(:,k),x(:,k-1));
    Up(:,k)=kron(p(:,k),u(:,k));
end

X0=[x0,x(:,1:Nd-1)];
X1=x(:,1:Nd);
U=u;
W=w;
G=[X0;Xp;U;Up];
end
