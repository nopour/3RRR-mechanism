clc
clear
tic

close all


syms t q1 q2 q3 q4 q5 q6 q7 q8 q9
syms ddq1 ddq2 ddq3 dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 dq9
assume([ddq1 ddq2 ddq3 dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 dq9 t q1 q2 q3 q4 q5 q6 q7 q8 q9],'real')
r0=0.1;
y0=0.25;
ome_t=4*pi;
ome_r=12*pi;

%%

addpath('GFunction')


r0=0.1;
y0=0.25;
ome_t=4*pi;
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
q11(t)=r0*sin(ome_r*t)*cos(ome_t*t);
q22(t)=y0+r0*sin(ome_r*t)*sin(ome_t*t);
% q11(t)=r0*sin(ome_t*t);
% q22(t)=r0*sin(pi/2*sin(ome_t*t));

%w=4*pi;

% q11(t)=r0*sin(w*t);
% q22(t)=r0*sin(pi/2*sin(w*t));


dq_1(t)=diff(q11(t) ,t);
dq_2(t)=diff(q22(t) ,t);

ddq_1(t)=diff(q11(t),t,2);
ddq_2(t)=diff(q22(t),t,2);

Cx1=q1-a/2*cos(q3)+sqrt(3)/6*a*sin(q3);
Cx2=q2-a/2*sin(q3)-sqrt(3)/6*a*cos(q3);
Cx3=q1+a/2*cos(q3)+sqrt(3)/6*a*sin(q3);
Cx4=q2+a/2*sin(q3)-sqrt(3)/6*a*cos(q3);
Cx5=q1-sqrt(3)/3*a*sin(q3);
Cx6=q2+sqrt(3)/3*a*cos(q3);

CX=[Cx1;Cx2;Cx3;Cx4;Cx5;Cx6];

Ct1=-b/2+l1*cos(q4)+l2*cos(q5);
Ct2=-sqrt(3)*b/6+l1*sin(q4)+l2*sin(q5);
Ct3=b/2+l1*cos(q6)+l2*cos(q7);
Ct4=-sqrt(3)*b/6+l1*sin(q6)+l2*sin(q7);
Ct5=l1*cos(q8)+l2*cos(q9);
Ct6=sqrt(3)*b/3+l1*sin(q8)+l2*sin(q9);

CT=[Ct1;Ct2;Ct3;Ct4;Ct5;Ct6];

X=[q1 q2 q3].';
DX=[dq1 dq2 dq3].';
DDX=[ddq1 ddq2 ddq3].';
T=[q4 q5 q6 q7 q8 q9].';
JX=jacobian(CX,X);
JT=jacobian(CT,T);

DT=(JT)\JX*DX;
DDT=(JT)\(jacobian(JX*DX,X)*DX+JX*DDX-jacobian(JT*DT,T)*DT);

    Mt=[0,0,0,1,0,0,0,0,0;
        0,0,0,0,0,1,0,0,0;
        0,0,0,0,0,0,0,1,0].';


N=30;
time=linspace(0,0.25,N);
% time=load('time05.mat');
% time=time.time05;
%%

for i=1:N
    tt=time(i);
    q_1=q11(tt);
    q_2=q22(tt);
    q_3=0;

    CX_1=double(subs(CX,[q1 q2 q3],[q_1 q_2 q_3]));
    [q_4,q_5]=solve(CT(1:2)==CX_1(1:2),[q4 q5]);
    [q_6,q_7]=solve(CT(3:4)==CX_1(3:4),[q6 q7]);
    [q_8,q_9]=solve(CT(5:6)==CX_1(5:6),[q8 q9]);

    q_44(i)=q_4(2);
    q_55(i)=q_5(2);
    q_66(i)=q_6(2);
    q_77(i)=q_7(2);
    q_88(i)=q_8(2);
    q_99(i)=q_9(2);

    
    
    
    dq_1_1=double(dq_1(tt));
    dq_2_2=double (dq_2(tt));
    ddq_1_1=double (ddq_1(tt));
    ddq_2_2=double (ddq_2(tt));
%     
    D=double(subs(DT,[q1 q2 q3 q4 q5 q6 q7 q8 q9 dq1 dq2 dq3],[q_1 q_2 q_3 q_4(2) q_5(2) q_6(2) q_7(2) q_8(2) q_9(2) dq_1_1 dq_2_2 0]));
    DD=double(subs(DDT,[q1 q2 q3 q4 q5 q6 q7 q8 q9 dq1 dq2 dq3 ddq1 ddq2 ddq3],[q_1 q_2 q_3 q_4(2) q_5(2) q_6(2) q_7(2) q_8(2) q_9(2) dq_1_1 dq_2_2 0 ddq_1_1 ddq_2_2 0]));


    M = Mfunc(A_1,A_2,Iz1,Iz2,Ize,Rho1,Rho2,l1,l2,me,q_4(2),q_5(2),q_6(2),q_7(2),q_8(2),q_9(2));
    H = Hfunc(A_2,Rho2,D(1),D(2),D(3),D(4),D(5),D(6),l1,l2,q_4(2),q_5(2),q_6(2),q_7(2),q_8(2),q_9(2));
    A1 = A1func(a,l1,l2,q_3,q_4(2),q_5(2),q_6(2),q_7(2),q_8(2),q_9(2));

    MTT=double([Mt -A1.']);

    Ans1(i,:)=double((MTT)\(M*[ddq_1_1;ddq_2_2;0;DD]+H))

    tt
    i
end

toc
%%

figure
plot(time,Ans1(:,1),'color','#A2142F','linewidth',1.5)
hold on
plot(time,Ans1(:,2),'color','#4DBEEE','linewidth',1.5)
hold on
plot(time,Ans1(:,3),'color','#2E515E','linewidth',1.5)

grid minor
xlabel('{Time [s]} ','Interpreter','latex')
ylabel('{Torque [N.m] }','Interpreter','latex')
legend('${\tau_1}$','${\tau_2}$','${\tau_3}$','Interpreter','latex','Location','best')

%%
figure
plot(time,Ans1(:,4),'color','#A2142F','linewidth',1.5)
figure
plot(time,Ans1(:,5),'color','#4DBEEE','linewidth',1.5)
figure
plot(time,Ans1(:,6),'color','#2E515E','linewidth',1.5)
figure
plot(time,Ans1(:,7),'color','#2E515E','linewidth',1.5)
figure
plot(time,Ans1(:,8),'color','#2E515E','linewidth',1.5)
figure
plot(time,Ans1(:,9),'color','#2E515E','linewidth',1.5)

%%


%     x=optimvar('x',2);
% %     CT=subs(CT,[t11 t12 t13],[x(1) x(2) x(3)]);
%     eq1=(2*cos(x(1)))/5 + (2*cos(x(2)))/5 - 9/20==CX_1(1);
%     eq2=(2*sin(x(1)))/5 + (2*sin(x(2)))/5 - (3*3^(1/2))/20==CX_1(2);
% 
%     prob=eqnproblem;
%     prob.Equations.eq1=eq1;
%     prob.Equations.eq2=eq2;
% 
% %     show(prob)
%     x0.x=[0 0];
%     [sol,fval,exitflag] = solve(prob,x0);
% q_41=double(sol.x(1))
% q_51=double(sol.x(2))
% clear x
% 
%     x=optimvar('x',2);
% %     CT=subs(CT,[t11 t12 t13],[x(1) x(2) x(3)]);
%     eq1=(2*cos(x(1)))/5 + (2*cos(x(2)))/5 + 9/20==CX_1(3);
%     eq2=(2*sin(x(1)))/5 + (2*sin(x(2)))/5 - (3*3^(1/2))/20==CX_1(4);
% 
%     prob=eqnproblem;
%     prob.Equations.eq1=eq1;
%     prob.Equations.eq2=eq2;
% 
% %     show(prob)
%     x0.x=[0 0];
%     [sol,fval,exitflag] = solve(prob,x0);
% q_61=double(sol.x(1))
% q_71=double(sol.x(2))
% clear x
% 
%     x=optimvar('x',2);
% %     CT=subs(CT,[t11 t12 t13],[x(1) x(2) x(3)]);
%     eq1=(2*cos(x(1)))/5 + (2*cos(x(2)))/5==CX_1(5);
%     eq2=(2*sin(x(1)))/5 + (2*sin(x(2)))/5 + (3*3^(1/2))/10==CX_1(6);
% 
%     prob=eqnproblem;
%     prob.Equations.eq1=eq1;
%     prob.Equations.eq2=eq2;
% 
% %     show(prob)
%     x0.x=[0 0];
%     [sol,fval,exitflag] = solve(prob,x0);
% q_81=double(sol.x(1))
% q_91=double(sol.x(2))
% clear x
%%
syms x y z
r2_bar=[x y z].';
for j=1:size(q_55.')
ro1_prime=[-b/2+l1*cos(q_44(j)) -sqrt(3)/6*b+l1*sin(q_44(j)) 0].';
Qz2=[cos(q_55(j)) -sin(q_55(j)) 0;sin(q_55(j)) cos(q_55(j)) 0;0 0 1];
Qze=[cos(q_3) -sin(q_3) 0;sin(q_3) cos(q_3) 0;0 0 1];
r2(x,y,z)=ro1_prime+Qz2*r2_bar;
re1_band=[a/2 sqrt(3)/6*a 0].';
reodprime1=Qze*re1_band;
ro1dprime=r2(l2,0,0);
re1=ro1dprime+reodprime1;

rex(j)=re1(1);
rey(j)=re1(2);
end
%
figure
plot(rex,rey)
