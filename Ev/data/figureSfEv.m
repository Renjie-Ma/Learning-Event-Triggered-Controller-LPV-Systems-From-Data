figure(1)  
subplot(3,1,1)
  plot(time,x1,'b','Linewidth',1);
  xlim([0,dt*tmax]);
  ylim([-3,3]);
  xlabel('Time (Sec.)'); 
[up1,lo1] = envelope(x1,1000,'analytic');
hold on

plot(time,up1,'-','Linewidth',0.5);
hold on
plot(time,lo1,'--','Linewidth',0.5);
hh=legend('State $ x_{k} $','Signal Upper Bound','Signal Lower Bound');

set(hh,'Interpreter','latex');
hold off

%%设置局部放大图

subAxesPosition1 = [0.8,0.5,1.2,2];
zoomAreaPosition1 = [0, -1, 0.25, 2];


zp1=BaseZoom(subAxesPosition1,zoomAreaPosition1);
zp1.run;


  subplot(3,1,2)
  plot(time,x1_ev,'b','Linewidth',1);

  xlim([0,dt*tmax]);
  ylim([-3,3]);
 xlabel('Time (Sec.)'); 

[up2,lo2] = envelope(x1_ev,5000,'analytic');
hold on
plot(time,up2,'-','Linewidth',0.5);

hold on
plot(time,lo2,'--','Linewidth',0.5);

hh=legend('State under ETM $ x_{k} $','Signal Upper Bound','Signal Lower Bound');

set(hh,'Interpreter','latex');
hold off

subAxesPosition2 = [0.8,0.5,1.2,2];
zoomAreaPosition2 = [0, -1, 0.25, 2];


zp2=BaseZoom(subAxesPosition2,zoomAreaPosition2);
zp2.run;

subplot(3,1,3)
stem(dt*seqiE,dt*intiE,'color','(0.0,0.45,0.74)','linewidth',0.5);
xlabel('Time (Sec.)'); 
xlim([0,dt*tmax]);
hh = legend('$\tau $');
set(hh,'Interpreter','latex','fontsize', 12);

%%
for k=1:200
x2_ev(k)=x_ev{k}(2);
end
figure(1)  
subplot(2,1,1)
plot(time,x1_ev,'b','Linewidth',1.5);
 hold on
plot(time,x2_ev,'color','#EDB120','Linewidth',1.5);
% plot(time,x2_ev,'color','#D95319','Linewidth',1.5);


  xlim([0,dt*tmax]);
  ylim([-3,3]);
 xlabel('Time (Sec.)'); 

hh=legend('State under ETM $ {x_1}_{k} $','State under ETM $ {x_2}_{k} $','Signal Upper Bound','Signal Lower Bound');

set(hh,'Interpreter','latex');
hold off

idx=find((time<=0.1) & (time>=0));
x1_ev1=x1_ev(idx);
x2_ev2=x2_ev(idx);

ax2=axes();
ax2.Position=[1,1,1.2,1.2];
p2=plot(time(idx),x1_ev1,'b','Linewidth',1);
hold on
p3=plot(time(idx),x2_ev2,'color','#EDB120','Linewidth',1);

p2.Parent=ax2;
p3.Parent=ax2;


subplot(2,1,2)
stem(dt*seqiE,dt*intiE,'color','(0.0,0.45,0.74)','linewidth',0.5);
xlabel('Time (Sec.)'); 
xlim([0,dt*tmax]);
hh = legend('$\tau $');
set(hh,'Interpreter','latex','fontsize', 12);
