syms x x_dot theta1 theta1_dot theta2 theta2_dot p11 p12 p13 p14 p15 p16 p21 p22 p23 p24 p25 p26 p31 p32 p33 p34 p35 p36 p41 p42 p43 p44 p45 p46 p51 p52 p53 p54 p55 p56 p61 p62 p63 p64 p65 p66

g = 10;
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;

A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
AT = transpose(A);
B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
BT = transpose(B);
R = 1;

P_1 = [p11 p12 p13 p14 p15 p16; p21 p22 p23 p24 p25 p26; p31 p32 p33 p34 p35 p36; p41 p42 p43 p44 p45 p46; p51 p52 p53 p54 p55 p56; p61 p62 p63 p64 p65 p66];
Q = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];

%ric = AT*P + P*A - P*B*(R^-1)*BT*P == -Q;
%solve(ric,P)

%solves ricatti for P, then find control law K
[X, L, G] = icare(A,B,Q,R);
P = X;
e = eig(P);
K = -(R^-1)*BT*P


%check for stability
A_new = A+B*K;
transpose(A_new);
P_lyap = lyap(A_new,Q);
eig(P_lyap);

