function A1 = A1func(a,l1,l2,q3,q4,q5,q6,q7,q8,q9)
%A1func
%    A1 = A1func(A,L1,L2,Q3,Q4,Q5,Q6,Q7,Q8,Q9)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    30-Jul-2022 10:20:11

t2 = sqrt(3.0);
t3 = pi./3.0;
t4 = q3+t3;
A1 = reshape([-1.0,0.0,-1.0,0.0,-1.0,0.0,0.0,-1.0,0.0,-1.0,0.0,-1.0,a.*t2.*cos(q3-t3).*(-1.0./3.0),(a.*t2.*cos(q3+pi./6.0))./3.0,a.*t2.*cos(t4).*(-1.0./3.0),a.*t2.*sin(t4).*(-1.0./3.0),(a.*t2.*cos(q3))./3.0,(a.*t2.*sin(q3))./3.0,-l1.*sin(q4),l1.*cos(q4),0.0,0.0,0.0,0.0,-l2.*sin(q5),l2.*cos(q5),0.0,0.0,0.0,0.0,0.0,0.0,-l1.*sin(q6),l1.*cos(q6),0.0,0.0,0.0,0.0,-l2.*sin(q7),l2.*cos(q7),0.0,0.0,0.0,0.0,0.0,0.0,-l1.*sin(q8),l1.*cos(q8),0.0,0.0,0.0,0.0,-l2.*sin(q9),l2.*cos(q9)],[6,9]);
