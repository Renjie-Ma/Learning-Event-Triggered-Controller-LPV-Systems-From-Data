figure(1)

  subplot(2,1,1)
  plot(time,y1,'b','Linewidth',1)
  hold on
  plot(time,r1,'r--','Linewidth',1);
  xlim([0,tmax*dt]);
  ylim([-2.5,2.5]);
 xlabel('Time (Sec.)'); 
hh = legend('Output Signal $ y_{k} $ under ETM ','Reference Signal $ r_{k} $');
set(hh,'Interpreter','latex');
hold off

idx=find((time<=0.05) & (time>=0));
y_1=y1(idx);
r_1=r1(idx);

ax2=axes();
ax2.Position=[1,1,1.2,1.2];
p2=plot(time(idx),y_1,'b','Linewidth',1);
hold on
p3=plot(time(idx),r_1,'r--','Linewidth',1);

p2.Parent=ax2;
p3.Parent=ax2;

subplot(2,1,2)
stem(dt*seqiE,dt*intiE,'color','(0.0,0.45,0.74)','linewidth',0.5);
xlabel('Time (Sec.)'); 
xlim([0,tmax*dt]);
hh = legend('$\tau $');
set(hh,'Interpreter','latex','fontsize', 12);
%%

idx_1=find((dt*seqiE<=0.05) & (dt*seqiE>=0));
y_1=dt*intiE(idx);

ax2=axes();
ax2.Position=[1,1,1.2,1.2];
p2=stem(dt*seqiE(idx_1),y_1,'color','(0.0,0.45,0.74)','linewidth',0.5);
p2.Parent=ax2;
%%

idx_2=find((dt*seqiE<=1.55) & (dt*seqiE>=1.5));
y_2=dt*intiE(idx_2);

ax3=axes();
ax3.Position=[1,1,1.2,1.2];
p3=stem(dt*seqiE(idx_2),y_2,'color','(0.0,0.45,0.74)','linewidth',0.5);
p3.Parent=ax3;

%%

idx_3=find((dt*seqiE<=3.05) & (dt*seqiE>=3));
y_3=dt*intiE(idx_3);

ax4=axes();
ax4.Position=[1,1,1.2,1.2];
p4=stem(dt*seqiE(idx_3),y_3,'color','(0.0,0.45,0.74)','linewidth',0.5);
p4.Parent=ax4;

%%

idx_4=find((dt*seqiE<=4.55) & (dt*seqiE>=4.5));
y_4=dt*intiE(idx_4);

ax5=axes();
ax5.Position=[1,1,1.2,1.2];
p5=stem(dt*seqiE(idx_4),y_4,'color','(0.0,0.45,0.74)','linewidth',0.5);
p5.Parent=ax5;