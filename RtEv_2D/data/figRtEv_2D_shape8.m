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

%%
y1_1=y1(1,1:3333);
y1_2=y1(1,1:6666);
y1_3=y1(1,1:9999);

y2_1=y2(1,1:3333);
y2_2=y2(1,1:6666);
y2_3=y2(1,1:9999);

figure(2)
subplot(1,3,1)

plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_1,y2_1,'LineWidth',0.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal $ y_{k} $ ');
set(hh,'Interpreter','latex');
hh1=text(1.5, -3, 'Time: 0 - 10 (Sec.)', ...
    'Color', 'red', 'FontSize', 12, 'HorizontalAlignment', 'center');
set(hh1,'Interpreter','latex');



subplot(1,3,2)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_2,y2_2,'LineWidth',0.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal $ y_{k} $ ');
set(hh,'Interpreter','latex');
hh1=text(1.5, -3, 'Time: 10 - 20 (Sec.)', ...
    'Color', 'red', 'FontSize', 12, 'HorizontalAlignment', 'center');
set(hh1,'Interpreter','latex');



subplot(1,3,3)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_3,y2_3,'LineWidth',0.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal $ y_{k} $ ');
set(hh,'Interpreter','latex');
hh1=text(1.5, -3, 'Time: 20 - 30 (Sec.)', ...
    'Color', 'red', 'FontSize', 12, 'HorizontalAlignment', 'center');
set(hh1,'Interpreter','latex');
%%

y1_ev_1=y1_ev(1,1:3333);
y1_ev_2=y1_ev(1,1:6666);
y1_ev_3=y1_ev(1,1:9999);

y2_ev_1=y2_ev(1,1:3333);
y2_ev_2=y2_ev(1,1:6666);
y2_ev_3=y2_ev(1,1:9999);

figure(3)
plot(r1,r2);
hold on 
plot(y1_ev,y2_ev);
xlim([-3,3]);
ylim([-3,3]);


subplot(2,3,1)

plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_ev_1,y2_ev_1,'LineWidth',0.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
hh1=text(1.5, -3, 'Time: 0 - 10 (Sec.)', ...
    'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
set(hh1,'Interpreter','latex');


subplot(2,3,2)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_ev_2,y2_ev_2,'LineWidth',0.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
hh1=text(1.5, -3, 'Time: 10 - 20 (Sec.)', ...
    'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
set(hh1,'Interpreter','latex');


subplot(2,3,3)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_ev_3,y2_ev_3,'LineWidth',0.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
hh1=text(1.5, -3, 'Time: 20 - 30 (Sec.)', ...
    'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
set(hh1,'Interpreter','latex');


subplot(2,3,[4:6])
stem(dt*seqiE,dt*intiE,'color','(0.0,0.45,0.74)','linewidth',0.5);
xlabel('Time (Sec.)'); 
xlim([0,tmax*dt]);
hh = legend('$\tau $');
set(hh,'Interpreter','latex','fontsize', 12);

%%
y1_ev_1=y1_ev(1,1:tmax/3);
y1_ev_2=y1_ev(1,1:tmax*2/3);
y1_ev_3=y1_ev(1,1:tmax-1);

y2_ev_1=y2_ev(1,1:tmax/3);
y2_ev_2=y2_ev(1,1:tmax*2/3);
y2_ev_3=y2_ev(1,1:tmax-1);

figure(3)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_ev_1,y2_ev_1,'LineWidth',1.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
% hh1=text(1.5, -3, 'Time: 0 - 10 (Sec.)', ...
%     'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
% set(hh1,'Interpreter','latex');
axis equal

figure(4)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_ev_2,y2_ev_2,'LineWidth',1.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
% hh1=text(1.5, -3, 'Time: 10 - 20 (Sec.)', ...
%     'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
% set(hh1,'Interpreter','latex');
axis equal

figure(5)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_ev_3,y2_ev_3,'LineWidth',1.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
% hh1=text(1.5, -3, 'Time: 20 - 30 (Sec.)', ...
%     'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
% set(hh1,'Interpreter','latex');
axis equal
%%
%%
y1_1=y1(1,1:tmax/3);
y1_2=y1(1,1:tmax*2/3);
y1_3=y1(1,1:tmax-1);

y2_1=y2(1,1:tmax/3);
y2_2=y2(1,1:tmax*2/3);
y2_3=y2(1,1:tmax-1);

figure(3)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_1,y2_1,'LineWidth',1.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
% hh1=text(1.5, -3, 'Time: 0 - 10 (Sec.)', ...
%     'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
% set(hh1,'Interpreter','latex');
axis equal

figure(4)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_2,y2_2,'LineWidth',1.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
% hh1=text(1.5, -3, 'Time: 10 - 20 (Sec.)', ...
%     'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
% set(hh1,'Interpreter','latex');
axis equal

figure(5)
plot(r1,r2,'--d','MarkerIndices',1:500:length(r1),'LineWidth',1);
hold on 
plot(y1_3,y2_3,'LineWidth',1.5);
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
hh = legend('Reference Signal $ r_{k} $','Output Signal under ETM $ y_{k} $ ');
set(hh,'Interpreter','latex');
% hh1=text(1.5, -3, 'Time: 20 - 30 (Sec.)', ...
%     'Color', 'red', 'FontSize', 6, 'HorizontalAlignment', 'center');
% set(hh1,'Interpreter','latex');
axis equal


%%
stem(dt*seqiE,dt*intiE,'color','(0.0,0.45,0.74)','linewidth',0.5);
xlabel('Time (Sec.)'); 
xlim([0,tmax*dt]);
hh = legend('$\tau $');
set(hh,'Interpreter','latex','fontsize', 12);

%%
idx_1=find((dt*seqiE<=0.4) & (dt*seqiE>=0));
y_1=dt*intiE(idx_1);

ax2=axes();
ax2.Position=[1,1,1.2,1.2];
p2=stem(dt*seqiE(idx_1),y_1,'color','(0.0,0.45,0.74)','linewidth',0.5);
p2.Parent=ax2;

%%
idx_2=find((dt*seqiE<=15) & (dt*seqiE>=14.25));
y_2=dt*intiE(idx_2);

ax3=axes();
ax3.Position=[1,1,1.2,1.2];
p3=stem(dt*seqiE(idx_2),y_2,'color','(0.0,0.45,0.74)','linewidth',0.5);
p3.Parent=ax3;