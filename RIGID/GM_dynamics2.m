function dZ=GM_dynamics2(t,Z)

dZ=Z;

global A_1 A_2 Rho1 Rho2 l1 l2 Iz1 Iz2 Ize me a
 

q=Z(1:9); dq=Z(10:18);

q1=q(1); q2=q(2); q3=q(3); q4=q(4); 
q5=q(5); q6=q(6); q7=q(7); q8=q(8); 
q9=q(9);

r0=0.1;
y0=0.25;
ome_t=4*pi;
ome_r=12*pi;

q(1)=r0*sin(ome_r*t)*cos(ome_t*t);
q(2)=y0+r0*sin(ome_r*t)*sin(ome_t*t);

dq(1)=ome_r*r0*cos(ome_r*t)*cos(ome_t*t) - ome_t*r0*sin(ome_r*t)*sin(ome_t*t);
dq(2)=ome_r*r0*cos(ome_r*t)*sin(ome_t*t) + ome_t*r0*cos(ome_t*t)*sin(ome_r*t);


dq3=dq(3); 
dq4=dq(4);  
dq5=dq(5); 
dq6=dq(6); 
dq7=dq(7); 
dq8=dq(8); 
dq9=dq(9); 



M = Mfunc(A_1,A_2,Iz1,Iz2,Ize,Rho1,Rho2,l1,l2,me,q4,q5,q6,q7,q8,q9);
H = Hfunc(A_2,Rho2,dq4,dq5,dq6,dq7,dq8,dq9,l1,l2,q4,q5,q6,q7,q8,q9);
A1 = A1func(a,l1,l2,q3,q4,q5,q6,q7,q8,q9).';
A1dot = A1dotfunc(a,dq3,dq4,dq5,dq6,dq7,dq8,dq9,l1,l2,q3,q4,q5,q6,q7,q8,q9).';



FTT=[0,0,0,1,0,1,0,1,0].';


F=FTT-H;


F_G=F-A1'*(pinv(A1*(pinv(M)*A1'))*(A1dot*dq+A1*(pinv(M)*F)));

ddq=pinv(M)*F_G;


dZ=[dq;ddq];

% size(FTT)
% size(A1)
% size(M)
% size(ddq)
% size(H)

dz1=pinv([FTT -A1.'])*(M*ddq+H);

dZ=[dq;ddq;dz1];
end