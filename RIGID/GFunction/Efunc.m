function out1 = Efunc(A_1,A_2,Iz1,Iz2,Ize,Rho1,Rho2,dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,dq9,l1,l2,me,q4,q5,q6,q7,q8,q9)
%Efunc
%    OUT1 = Efunc(A_1,A_2,Iz1,Iz2,Ize,Rho1,Rho2,DQ1,DQ2,DQ3,DQ4,DQ5,DQ6,DQ7,DQ8,DQ9,L1,L2,ME,Q4,Q5,Q6,Q7,Q8,Q9)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    30-Jul-2022 10:20:14

t2 = dq4.^2;
t3 = dq6.^2;
t4 = dq8.^2;
t5 = l1.^2;
t6 = l2.^2;
t7 = Iz1./2.0;
t8 = Iz2./2.0;
t9 = (A_1.*t5)./6.0;
t10 = (A_2.*t6)./6.0;
t11 = t7+t9;
t12 = t8+t10;
out1 = (Ize.*dq3.^2)./2.0+(me.*(dq1.^2+dq2.^2))./2.0+Rho1.*l1.*t2.*t11+Rho1.*l1.*t3.*t11+Rho1.*l1.*t4.*t11+Rho2.*dq5.^2.*l2.*t12+Rho2.*dq7.^2.*l2.*t12+Rho2.*dq9.^2.*l2.*t12+(A_2.*Rho2.*l2.*t2.*t5)./2.0+(A_2.*Rho2.*l2.*t3.*t5)./2.0+(A_2.*Rho2.*l2.*t4.*t5)./2.0+(A_2.*Rho2.*dq4.*dq5.*l1.*t6.*cos(q4-q5))./2.0+(A_2.*Rho2.*dq6.*dq7.*l1.*t6.*cos(q6-q7))./2.0+(A_2.*Rho2.*dq8.*dq9.*l1.*t6.*cos(q8-q9))./2.0;