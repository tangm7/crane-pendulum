function dx = crane_diffeq_Nonlinear_LQG_ConstantReference_fxn(x,t,u,K_optimal,L_optimal,Udstddev, Vstddev)
% Note: K_optimal passed in isn't actually used, it is embeddd within the u anonymous function

xhat = [ x(7);  x(8);  x(9);  x(10);  x(11);  x(12) ]; % Estimated States
F = u(xhat);

% note: F = u(x) moved further down to use xhat

g = 10;
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;

X = [x(1);  x(2);  x(3);  x(4);  x(5);  x(6);  x(7);  x(8);  x(9);  x(10);  x(11);  x(12)]; % this is actually Xe now, the errors


A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];


% A_nonlinear = A;
A_nonlinear = ...
[0, 1, 0,                                                                                                      0, 0,                                                                                                      0;
0, 0, 0,                 -(l1*m1*sin(x(3))*x(4))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2), 0,                 -(l2*m2*sin(x(5))*x(6))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2);
0, 0, 0,                                                                                                      1, 0,                                                                                                      0;
0, 0, 0,         -(m1*cos(x(3))*sin(x(3))*x(4))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2), 0, -(l2*m2*cos(x(3))*sin(x(5))*x(6))/(l1*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));
0, 0, 0,                                                                                                      0, 0,                                                                                                      1;
0, 0, 0, -(l1*m1*cos(x(5))*sin(x(3))*x(4))/(l2*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2)), 0,         -(m2*cos(x(5))*sin(x(5))*x(6))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2)];
 

% B_nonlinear = B;
B_nonlinear = [...
                                                                  0;
              1/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2);
                                                                  0;
cos(x(3))/(l1*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));
                                                                  0;
cos(x(5))/(l2*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));];
 



 E_nonlinear = [...
 
                                                                                                                                                                                      0;
                                                                         -(g*m1*cos(x(3))*sin(x(3)) + g*m2*cos(x(5))*sin(x(5)))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2);
                                                                                                                                                                                      0;
-(- g*m2*sin(x(3))*cos(x(5))^2 + g*m2*cos(x(3))*sin(x(5))*cos(x(5)) + M*g*sin(x(3)) + g*m1*sin(x(3)) + g*m2*sin(x(3)))/(l1*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));
                                                                                                                                                                                      0;
-(- g*m1*sin(x(5))*cos(x(3))^2 + g*m1*cos(x(5))*sin(x(3))*cos(x(3)) + M*g*sin(x(5)) + g*m1*sin(x(5)) + g*m2*sin(x(5)))/(l2*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2))];

% E_nonlinear = zeros(size(E_nonlinear));


%%%%%%%%%%%%%%%%%%%%%%%%%% Pick one C matrix of the three below
%%% x(t) is output vector
C_xonly = [1 0 0 0 0 0];

% %%%%%%%%%%%%%%%%%%%%%%%%%%% Pick one respective observer of the three below to use
% p_xonly = [-10 -20 -30 -40 -50 -60];
% L_xonly = place(A',C_xonly',p_xonly).';
% L_xonly = L_optimal;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Pick one respective observer of the three below to use
% Ac = [A_nonlinear                             -B_nonlinear*K_optimal;
%      L_optimal*C_xonly           A-L_optimal*C_xonly - B*K_optimal];

Ac = [A_nonlinear              zeros(size((A-L_optimal*C_xonly)));
    L_optimal*C_xonly                      A-L_optimal*C_xonly];

%%%%%%%%%%%%%%%%%%%%%%% B for inputs
Bu = [B_nonlinear;   B];

%%%%%%%%%%%%%%%%%%%%%%%% B is the same for all
B_D = eye(6); % process noise directly affects dynamics
Bc = [B_D;    zeros(size(B_D))]; % combined B

%%%%%%%%%%%%%%%%%%%%%%%% Measurement Noise
Dc = [zeros(size(L_optimal));  L_optimal];

%%%%%%%%%%%%%%%%%%%%%%% add E for nonlinear system
Ec = [E_nonlinear; zeros(size(E_nonlinear))];

%%%%%%%%%%%%%%%%%%%%%%% Define Noise
Ud = Udstddev*randn(6,1); % process noise has standard deviation of .1 (equates to about 5 degrees for angles, .1 meters = 4 in, or .1 m/s)

V = Vstddev*randn(1,1); % measurement noise has standard deviation of .1 (equates to about 5 degrees for angles, .1 meters = 4 in, or .1 m/s) % CORRECTED: V is 1x1

%%%%%%%%%%%%%%%%%%%%%%%%% Pick one respective (twelve state) state equation
Xdot = Ac * X + Bu*F + Bc*Ud + Dc*V + Ec;



dx(1,1) = Xdot(1);
dx(2,1) = Xdot(2);

dx(3,1) = Xdot(3);

dx(4,1) = Xdot(4);
dx(5,1) = Xdot(5);
dx(6,1) = Xdot(6);


dx(7,1) = Xdot(7);
dx(8,1) = Xdot(8);

dx(9,1) = Xdot(9);

dx(10,1) = Xdot(10);
dx(11,1) = Xdot(11);
dx(12,1) = Xdot(12);

end