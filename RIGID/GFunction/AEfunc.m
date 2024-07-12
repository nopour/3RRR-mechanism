function AE = AEfunc(a,l1,l2,q3,q4,q5,q6,q7,q8,q9)
%AEfunc
%    AE = AEfunc(A,L1,L2,Q3,Q4,Q5,Q6,Q7,Q8,Q9)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    30-Jul-2022 10:20:15

t2 = cos(q4);
t3 = cos(q5);
t4 = cos(q6);
t5 = cos(q7);
t6 = cos(q8);
t7 = cos(q9);
t8 = sin(q4);
t9 = sin(q5);
t10 = sin(q6);
t11 = sin(q7);
t12 = sin(q8);
t13 = sin(q9);
t14 = l1.^2;
t15 = l2.^2;
t28 = -q3;
t29 = -q4;
t30 = -q5;
t31 = -q7;
t32 = -q8;
t33 = -q9;
t34 = sqrt(3.0);
t35 = pi./6.0;
t16 = l1.*t2;
t17 = l1.*t4;
t18 = l2.*t3;
t19 = l1.*t6;
t20 = l2.*t5;
t21 = l2.*t7;
t22 = l1.*t8;
t23 = l1.*t10;
t24 = l2.*t9;
t25 = l1.*t12;
t26 = l2.*t11;
t27 = l2.*t13;
t42 = q4+t30;
t43 = q3+t32;
t44 = q3+t33;
t45 = q6+t31;
t46 = q8+t33;
t55 = q3+t29+t35;
t56 = q3+t30+t35;
t57 = q6+t28+t35;
t58 = q7+t28+t35;
t36 = -t16;
t37 = -t17;
t38 = -t18;
t39 = -t19;
t40 = -t20;
t41 = -t21;
t47 = cos(t42);
t48 = cos(t45);
t49 = cos(t46);
t50 = sin(t43);
t51 = sin(t44);
t59 = cos(t55);
t60 = cos(t56);
t61 = cos(t57);
t62 = cos(t58);
t52 = l1.*l2.*t47;
t53 = l1.*l2.*t48;
t54 = l1.*l2.*t49;
t63 = (a.*l1.*t34.*t50)./3.0;
t64 = (a.*l2.*t34.*t51)./3.0;
t65 = (a.*l1.*t34.*t59)./3.0;
t66 = (a.*l1.*t34.*t61)./3.0;
t67 = (a.*l2.*t34.*t60)./3.0;
t68 = (a.*l2.*t34.*t62)./3.0;
t69 = -t66;
t70 = -t68;
AE = reshape([3.0,0.0,0.0,t22,t24,t23,t26,t25,t27,0.0,3.0,0.0,t36,t38,t37,t40,t39,t41,0.0,0.0,a.^2,t65,t67,t69,t70,t63,t64,t22,t36,t65,t14,t52,0.0,0.0,0.0,0.0,t24,t38,t67,t52,t15,0.0,0.0,0.0,0.0,t23,t37,t69,0.0,0.0,t14,t53,0.0,0.0,t26,t40,t70,0.0,0.0,t53,t15,0.0,0.0,t25,t39,t63,0.0,0.0,0.0,0.0,t14,t54,t27,t41,t64,0.0,0.0,0.0,0.0,t54,t15],[9,9]);
