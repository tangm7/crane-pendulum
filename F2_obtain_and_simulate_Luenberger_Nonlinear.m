clear, clc;
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
R = .0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find L that makes the system VERY DETECTABLE/OBSERVABLE

%%% x(t) is output vector
C_xonly = [1 0 0 0 0 0];

%%% q1(t) q2(t) is output vector
C_q1_q2 = [0 0   1 0   1 0];

%%% x(t) theta2(t) is output vector
C_x_q2 = [1 0   0 0   1 0];

%%% x(t) theta1(t) theta2(t) is output vector
C_x_q1_q2 = [1 0   1 0   1 0];

syms L1 L2 L3 L4 L5 L6
Larray = [L1;  L2;  L3;  L4;  L5;  L6]; % L gain matrix

display(" Finding eigenvalues of A - Larray * C ")
% [Ve, De] = eig(A - Larray*C_xonly) % Note: pick a C matrix to do
% [L_xonly, S_xonly, E_xonly] = lqr(A-Larray*C_xonly, B, Q, R);

% p_xonly = [-10 -20 -30 -40 -50 -60]; % place very negative, 100 times further left than the usual poles
% p_all = [-10 -11 -12 -13 -14 -15];

p_all = [-10 -20 -30 -40 -50 -60];

p_xonly = [-10 -20 -30 -40 -50 -60]; % place very negative, 10 times further left than the usual poles
L_xonly = place(A',C_xonly',p_xonly).';

% p_q1_q2 = [-.1 -.2 -.3 -.4 -.5 -.6]; % not possible because output vectors are not observable!
% L_q1_q2 = place(A',C_q1_q2',p_all).'

% p_x_q2 = [-7.5 -7.6 -7.7 -7.8 -7.9 -8];
p_x_q2 = [-5.1 -5.2 -5.3 -5.4 -5.5 -5.6];
L_x_q2 = place(A',C_x_q2',p_x_q2).';

p_x_q1_q2 = [-5.1 -5.2 -5.3 -5.4 -5.5 -5.6];
L_x_q1_q2 = place(A',C_x_q1_q2',p_x_q1_q2).';






tspan = 0:.004:1;

%%%%  x    xdot    q1            q1d           q2           q2d
% x0 = [0;    0;  deg2rad(45);  deg2rad(0);  deg2rad(45);  deg2rad(0)]; % Define initial condition for Xe, the errors


%%%   x(1);  x(2);  x(3);         x(4);          x(5);        x(6)           x(7);  x(8);  x(9);  x(10);  x(11);  x(12)
% xc0 = [0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);        0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);  ]; % Define combined initial condition for X, the state and errors
% xc0 = [0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);        0;    0;  0;  deg2rad(0);  0;  deg2rad(0);  ]; % Define combined initial condition for X, the state and errors

% xc0 = [0;    0;  0;  0;  0;  0;        0;    0;  0;  0;  0;  0;  ]; % Define combined initial condition for X, the state and errors
xc0 = [0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);        0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);  ]; % Define combined initial condition for X, the state and errors

% u = @(x,t) -K*x
u = @(x,t) 1; % as per assignment, force is unit step input. Note control_input (see below) is using ones function right now.

[t,x] = ode45(  @(t,x)crane_diffeq_Nonlinear_Luenberger_fxn(x,t,u)    ,   tspan,   xc0);

size_x = size(x)
number_timesteps = size_x(1)
control_input = ones( number_timesteps ); % for plotting later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot results

figure

subplot(7,1,1) % Plot Cart Position
title(['Scenario ', 1])
plot(t, (x(:,1)), '.')
hold on
plot(t, (x(:,7)), '.')
hold off
grid on
legend('x','xhat')
xlabel('Time (s)') 
ylabel('(m)')



subplot(7,1,2) % Plot Cart Velocity
plot(t, (x(:,2)), 'LineWidth', 3)
hold on
plot(t, (x(:,8)), '.')
hold off
grid on
legend('xdot','xdothat')
xlabel('Time (s)') 
ylabel('(m/s)') 



subplot(7,1,3) % Plot pendulum 1 angle
% plot(t, (x(:,3)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,3))  ), 'LineWidth', 3)
hold on
plot(t, rad2deg(  (x(:,9))  ), '.')
hold off
grid on
legend('theta1','theta1hat')
xlabel('Time (s)') 
% ylabel('(rad)')
ylabel('(deg)')




subplot(7,1,4) % Plot pendulum 1 angular velocity
% plot(t, (x(:,4)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,4))  ), 'LineWidth', 3)
hold on
plot(t, rad2deg(  (x(:,10)) ), '.')
hold off
grid on
legend('theta1dot','theta1dothat')
xlabel('Time (s)') 
% ylabel('(rad/s)') 
ylabel('(deg/s)') 


subplot(7,1,5) % Plot pendulum 2 position
% plot(t, (x(:,5)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,5))  ), 'LineWidth', 3)
hold on
plot(t, rad2deg(  (x(:,11)) ), '.')
hold off
grid on
legend('theta2','theta2hat')
xlabel('Time (s)') 
% ylabel('(rad)') 
ylabel('(deg)') 



subplot(7,1,6) % Plot pendulum 2 angular velocity
% plot(t, (x(:,6)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,6))  ), 'LineWidth', 3)
hold on
plot(t, rad2deg(  (x(:,12)) ), '.')
hold off
grid on
legend('theta2dot','theta2dothat')
xlabel('Time (s)') 
% ylabel('(rad/s)') 
ylabel('(deg/s)') 


subplot(7,1,7) % Plot Cart Position
plot(t, (control_input), 'LineWidth', 3)
grid on
legend('u: control_input')
xlabel('Time (s)') 
ylabel('(N)')


