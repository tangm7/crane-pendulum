function dx = crane_diffeq_nonlinear_fxn(x,t,u)


% crane_diffeq_linear_fxn(x,t,u,g,m1,L1,l1,b1,m2,L2,l2,b2,Jhat0,Jhat2)
% tau1 = u(x); % motor input torque
% tau2 = 0; % no disturbance torques on pendulum arm

% For the state:             x  = (   q1           q2            q1_dot          q2_dot       )
% Derivative of state:    x_dot = (   q1_dot      q2_dot     q1_doubledot    q2_doubledot  )

F = u(x);

g = 10;
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;

X = [x(1);  x(2);  x(3);  x(4);  x(5);  x(6)];









A = ...
[0, 1, 0,                                                                                                      0, 0,                                                                                                      0;
0, 0, 0,                 -(l1*m1*sin(x(3))*x(4))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2), 0,                 -(l2*m2*sin(x(5))*x(6))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2);
0, 0, 0,                                                                                                      1, 0,                                                                                                      0;
0, 0, 0,         -(m1*cos(x(3))*sin(x(3))*x(4))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2), 0, -(l2*m2*cos(x(3))*sin(x(5))*x(6))/(l1*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));
0, 0, 0,                                                                                                      0, 0,                                                                                                      1;
0, 0, 0, -(l1*m1*cos(x(5))*sin(x(3))*x(4))/(l2*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2)), 0,         -(m2*cos(x(5))*sin(x(5))*x(6))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2)];
 
 
B = [...
                                                                  0;
              1/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2);
                                                                  0;
cos(x(3))/(l1*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));
                                                                  0;
cos(x(5))/(l2*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));];
 




 E = [...
 
                                                                                                                                                                                      0;
                                                                         -(g*m1*cos(x(3))*sin(x(3)) + g*m2*cos(x(5))*sin(x(5)))/(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2);
                                                                                                                                                                                      0;
-(- g*m2*sin(x(3))*cos(x(5))^2 + g*m2*cos(x(3))*sin(x(5))*cos(x(5)) + M*g*sin(x(3)) + g*m1*sin(x(3)) + g*m2*sin(x(3)))/(l1*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2));
                                                                                                                                                                                      0;
-(- g*m1*sin(x(5))*cos(x(3))^2 + g*m1*cos(x(5))*sin(x(3))*cos(x(3)) + M*g*sin(x(5)) + g*m1*sin(x(5)) + g*m2*sin(x(5)))/(l2*(- m1*cos(x(3))^2 - m2*cos(x(5))^2 + M + m1 + m2))];


% A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
% B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

% Xdot = A*X + B*F;
Xdot = A * X   +   B*F   +   E;



dx(1,1) = Xdot(1);
dx(2,1) = Xdot(2);

dx(3,1) = Xdot(3);

dx(4,1) = Xdot(4);
dx(5,1) = Xdot(5);
dx(6,1) = Xdot(6);

end