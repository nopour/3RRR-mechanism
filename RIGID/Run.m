clc
clear
close all

addpath('GFunction');

global A_1 A_2 Rho1 Rho2 l1 l2 Iz1 Iz2 Ize me a kk
kk=1;

r0=0.1;
y0=0.25;
ome_t=2*pi;
ome_r=12*pi;
a=8e-2;
b=90e-2;
l1=40e-2;
l2=40e-2;
A_1=23e-3*6e-3;
A_2=23e-3*2e-3;
Iz1=23*6^3/12*10^-12;
Iz2=23*2^3/12*10^-12;

Rho1=2740;
Rho2=7801;

E1=71.705E9;
E2=207E9;
me=0.1898;
Ize=1.0124;





% dt=0.1;
% timeSpan=0:dt:0.5;
timeSpan=linspace(0,0.25,200);
tsize=length(timeSpan);


q_0=zeros(9,1);
dq_0=zeros(9,1);

% options=odeset('AbsTol',0.5);

%%
t0=clock;
[~,zGM]=ode45(@GM_dynamics,timeSpan,[q_0;dq_0]);
t1=clock;
timeGM=etime(t1,t0);
disp(['GM sim time: ' num2str(timeGM) '(s)'])
%%
GM.q1=zGM(:,1);
GM.q2=zGM(:,2);
GM.q3=zGM(:,3);
GM.q4=zGM(:,4);
GM.q5=zGM(:,5);
GM.q6=zGM(:,6);
GM.q7=zGM(:,7);
GM.q8=zGM(:,8);
GM.q9=zGM(:,9);


GM.dq1=zGM(:,10);
GM.dq2=zGM(:,11);
GM.dq3=zGM(:,12);
GM.dq4=zGM(:,13);
GM.dq5=zGM(:,14);
GM.dq6=zGM(:,15);
GM.dq7=zGM(:,16);
GM.dq8=zGM(:,17);
GM.dq9=zGM(:,18);
%%
% AE_GM=AEfunc(a,GM.dq1,GM.dq2,GM.dq3,GM.dq4,GM.dq5,GM.dq6,GM.dq7,GM.dq8,GM.dq9,l1,l2,GM.q3,GM.q4,GM.q5,GM.q6,GM.q7,GM.q8,GM.q9);
% E_GM=Efunc(A_1,A_2,Iz1,Iz2,Ize,Rho1,Rho2,GM.dq1,GM.dq2,GM.dq3,GM.dq4,GM.dq5,GM.dq6,GM.dq7,GM.dq8,GM.dq9,l1,l2,me,GM.q4,GM.q5,GM.q6,GM.q7,GM.q8,GM.q9);;
% ER_GM=(E_GM-E_GM(1))/E_GM(1);
% plot(timeSpan,E_GM,'k-')

% 
% %%
% figure
% plot(timeSpan,AE_GM,'k--','linewidth',1.5)
% legend('IM','AM','EM','GM')
% xlabel('time(s)')
% grid minor
% title('Constraint Error')
%%
figure
plot(GM.q1,GM.q2,'color','#A2142F','linewidth',2.5)
grid minor
xlabel('{X [m]} ','Interpreter','latex')
ylabel('{Y [m] }','Interpreter','latex')
legend('$Trajectory$','Interpreter','latex','Location','best')
%%
figure
plot(timeSpan,GM.q7*180/pi,'color','#77AC30','linewidth',2.5)
grid minor
xlabel('{time [s]} ','Interpreter','latex')
ylabel('$\theta_4 [Deg] $','Interpreter','latex')
legend('$\theta_4$','Interpreter','latex','Location','best')
%%

% [~,zGM1]=ode45(@GM_dynamics2,timeSpan,[q_0;dq_0;zeros(7,1)],options);
% GM.t1=zGM1(:,19);
% GM.t2=zGM1(:,20);
% GM.t3=zGM1(:,21);
% GM.t4=zGM1(:,22);
% GM.t5=zGM1(:,23);
% %%
% figure
% plot(timeSpan,GM.t2,'color','#A2142F','linewidth',2.5);
% % hold on
% % plot(timeSpan,GM.t2);
% % hold on
% % plot(timeSpan,GM.t3);
%%
% syms x y z
% r2_bar=[x y z].';
% for j=1:size(GM.q5)
% ro1_prime=[-b/2+l1*cos(GM.q4(j)) -sqrt(3)/6*b+l1*sin(GM.q4(j)) 0].';
% Qz2=[cos(GM.q5(j)) -sin(GM.q5(j)) 0;sin(GM.q5(j)) cos(GM.q5(j)) 0;0 0 1];
% Qze=[cos(GM.q3(j)) -sin(GM.q3(j)) 0;sin(GM.q3(j)) cos(GM.q3(j)) 0;0 0 1];
% r2(x,y,z)=ro1_prime+Qz2*r2_bar;
% re1_band=[a/2 sqrt(3)/6*a 0].';
% reodprime1=Qze*re1_band;
% ro1dprime=r2(l2,0,0);
% re1=ro1dprime+reodprime1;
% 
% rex(j)=re1(1);
% rey(j)=re1(2);
% end
% %
% figure
% plot(rex,rey)