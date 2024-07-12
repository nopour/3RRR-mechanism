function H_matrix = Hfunc(A_2,Rho2,dq4,dq5,dq6,dq7,dq8,dq9,l1,l2,q4,q5,q6,q7,q8,q9)
%Hfunc
%    H_matrix = Hfunc(A_2,Rho2,DQ4,DQ5,DQ6,DQ7,DQ8,DQ9,L1,L2,Q4,Q5,Q6,Q7,Q8,Q9)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    30-Jul-2022 10:20:10

t2 = l2.^2;
t3 = -q5;
t4 = -q7;
t5 = -q9;
t6 = q4+t3;
t7 = q6+t4;
t8 = q8+t5;
t9 = sin(t6);
t10 = sin(t7);
t11 = sin(t8);
H_matrix = [0.0;0.0;0.0;(A_2.*Rho2.*dq5.^2.*l1.*t2.*t9)./2.0;A_2.*Rho2.*dq4.^2.*l1.*t2.*t9.*(-1.0./2.0);(A_2.*Rho2.*dq7.^2.*l1.*t2.*t10)./2.0;A_2.*Rho2.*dq6.^2.*l1.*t2.*t10.*(-1.0./2.0);(A_2.*Rho2.*dq9.^2.*l1.*t2.*t11)./2.0;A_2.*Rho2.*dq8.^2.*l1.*t2.*t11.*(-1.0./2.0)];
