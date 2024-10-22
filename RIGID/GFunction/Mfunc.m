function M = Mfunc(A_1,A_2,Iz1,Iz2,Ize,Rho1,Rho2,l1,l2,me,q4,q5,q6,q7,q8,q9)
%Mfunc
%    M = Mfunc(A_1,A_2,Iz1,Iz2,Ize,Rho1,Rho2,L1,L2,ME,Q4,Q5,Q6,Q7,Q8,Q9)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    30-Jul-2022 10:20:08

t2 = Iz2.*3.0;
t3 = l1.^2;
t4 = l2.^2;
t5 = -q5;
t6 = -q7;
t7 = -q9;
t9 = Iz1./2.0;
t8 = A_2.*t4;
t10 = q4+t5;
t11 = q6+t6;
t12 = q8+t7;
t16 = A_2.*Rho2.*l2.*t3;
t17 = (A_1.*t3)./6.0;
t13 = cos(t10);
t14 = cos(t11);
t15 = cos(t12);
t18 = t2+t8;
t19 = t9+t17;
t20 = (Rho2.*l2.*t18)./3.0;
t21 = (Rho2.*l1.*t8.*t13)./2.0;
t22 = (Rho2.*l1.*t8.*t14)./2.0;
t23 = (Rho2.*l1.*t8.*t15)./2.0;
t24 = Rho1.*l1.*t19.*2.0;
t25 = t16+t24;
M = reshape([me,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,me,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,Ize,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t25,t21,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t21,t20,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t25,t22,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t22,t20,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t25,t23,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t23,t20],[9,9]);
