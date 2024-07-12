clc
clear 
close all
%%
% q1=theta1   q2=theta2 

addpath('GFunction')

syms t x y z a b
syms q4 q5 q6 q7 q8 q9 q3 q1 q2 q3
syms dq4 dq5 dq6 dq7 dq8 dq9 dq3 dq1 dq2 dq3
syms l1 l2 m1 m2 m3 g Rho1 A_1 E1 I1 Rho2 A_2 E2 I2 Iz1 Iz2 me Ize
tor1=0;



%% Link 1
r1_bar=[x y z].';
ro1=[-b/2 -sqrt(3)/6*b 0].';
Qz1=[cos(q4) -sin(q4) 0;sin(q4) cos(q4) 0;0 0 1];
r1(x,y,z)=ro1+Qz1*r1_bar;

v1(x,y,z)=diff(r1,q4)*dq4;

%% Link 2
r2_bar=[x y z].';
ro1_prime=[-b/2+l1*cos(q4) -sqrt(3)/6*b+l1*sin(q4) 0].';
Qz2=[cos(q5) -sin(q5) 0;sin(q5) cos(q5) 0;0 0 1];
r2(x,y,z)=ro1_prime+Qz2*r2_bar;

v2(x,y,z)=diff(r2,q4)*dq4+diff(r2,q5)*dq5;

%% Link 3
r3_bar=[x y z].';
ro2=[-b/2 -sqrt(3)/6*b 0].';
Qz3=[cos(q6) -sin(q6) 0;sin(q6) cos(q6) 0;0 0 1];
r3(x,y,z)=ro2+Qz3*r3_bar;
v3(x,y,z)=diff(r3,q6)*dq6;
%% Link 4
r4_bar=[x y z].';
ro2_prime=[b/2+l1*cos(q6) -sqrt(3)/6*b+l1*sin(q6) 0].';
Qz4=[cos(q7) -sin(q7) 0;sin(q7) cos(q7) 0;0 0 1];
r4(x,y,z)=ro2_prime+Qz4*r4_bar;
v4(x,y,z)=diff(r4,q6)*dq6+diff(r4,q7)*dq7;
%% Link 5
r5_bar=[x y z].';
ro3=[0 +sqrt(3)/3*b 0].';
Qz5=[cos(q8) -sin(q8) 0;sin(q8) cos(q8) 0;0 0 1];
r5(x,y,z)=ro3+Qz5*r5_bar;
v5(x,y,z)=diff(r5,q8)*dq8;
%% Link 6
r6_bar=[x y z].';
ro3_prime=[l1*cos(q8) sqrt(3)/3*b+l1*sin(q8) 0].';
Qz6=[cos(q9) -sin(q9) 0;sin(q9) cos(q9) 0;0 0 1];
r6(x,y,z)=ro3_prime+Qz6*r6_bar;
v4(x,y,z)=diff(r6,q8)*dq8+diff(r6,q9)*dq9;
%% Mojri
r_e=[q1 q2 0].';
rend_bar=[x y z].';
Qze=[cos(q3) -sin(q3) 0;sin(q3) cos(q3) 0;0 0 1];
rend(x,y,z)=r_e+Qze*rend_bar;
vend(x,y,z)=diff(rend,q1)*dq1+diff(rend,q2)*dq2+diff(rend,q3)*dq3;
%% Kinematic Cons. p66
re1_band=[a/2 sqrt(3)/6*a 0].'; %37-3
reodprime1=Qze*re1_band; %38-3
ro1dprime=r2(l2,0,0);
re1=ro1dprime+reodprime1;

re2_band=[-a/2 sqrt(3)/6*a 0].';
reodprime2=Qze*re2_band;
ro2dprime=r4(l2,0,0);
re2=ro2dprime+reodprime2;

re3_band=[0 -sqrt(3)/3*a 0].';
reodprime3=Qze*re3_band;
ro3dprime=r6(l2,0,0);
re3=ro3dprime+reodprime3; %p.72
%%
T1=Rho1*l1*(A_1*l1^2/6+0.5*Iz1)*dq4^2;
T2=0.5*Rho2*l2*l1^2*A_2*dq4^2+0.5*Rho2*l1*l2^2*A_2*cos(q4-q5)*dq4*dq5+...
    +Rho2*l2*(A_2*l2^2/6+0.5*Iz2)*dq5^2;
T3=Rho1*l1*(A_1*l1^2/6+0.5*Iz1)*dq6^2;
T4=0.5*Rho2*l2*l1^2*A_2*dq6^2+0.5*Rho2*l1*l2^2*A_2*cos(q6-q7)*dq6*dq7+...
    +Rho2*l2*(A_2*l2^2/6+0.5*Iz2)*dq7^2;
T5=Rho1*l1*(A_1*l1^2/6+0.5*Iz1)*dq8^2;
T6=0.5*Rho2*l2*l1^2*A_2*dq8^2+0.5*Rho2*l1*l2^2*A_2*cos(q8-q9)*dq8*dq9+...
    +Rho2*l2*(A_2*l2^2/6+0.5*Iz2)*dq9^2;
Te=0.5*me*(dq1^2+dq2^2)+0.5*Ize*dq3^2;
T=T1+T2+T3+T4+T5+T6+Te;

%% Energy

E=T;
V=0;

q=[q1 q2 q3 q4 q5 q6 q7 q8 q9].';
dq=[dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8 dq9 ].';
% page83(108-3) (115-3)
Mass_matrix=jacobian(jacobian(T,dq),dq);
H_matrix=(jacobian(jacobian(T,dq),q))*dq-jacobian(T,q).'+jacobian(V,q).';
M1=[tor1;0;tor1;0;0;0;0;0];


C1=[q1;q2]-re1(1:2);
C2=[q1;q2]-re2(1:2);
C3=[q1;q2]-re3(1:2);
A2=[C1;C2;C3];

% A=jacobian(A2,q)*dq;
% 
% A1=jacobian(A,dq).';

A1=jacobian(A2,q);

AE=simplify(A1.'*A1);



A1dot=A1;
for n=1:size(A1,1)
       A1dot(n,:)= (jacobian(A1(n,:),q)*dq).';
end

M=simplify(Mass_matrix);
H_matrix=simplify(H_matrix);
A1=simplify(A1);
A1dot=simplify(A1dot);

%%

Mfunc=matlabFunction(M,'File','GFunction\Mfunc');
Hfunc=matlabFunction(H_matrix,'File','GFunction\Hfunc');
A1func=matlabFunction(A1,'File','GFunction\A1func');
A1dotfunc=matlabFunction(A1dot,'File','GFunction\A1dotfunc');
Efunc=matlabFunction(simplify(E),'File','GFunction\Efunc');
AEfunc=matlabFunction((AE),'File','GFunction\AEfunc');

