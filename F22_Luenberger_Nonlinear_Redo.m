clear, clc;
syms x x_dot theta1 theta1_dot theta2 theta2_dot p11 p12 p13 p14 p15 p16 p21 p22 p23 p24 p25 p26 p31 p32 p33 p34 p35 p36 p41 p42 p43 p44 p45 p46 p51 p52 p53 p54 p55 p56 p61 p62 p63 p64 p65 p66

g = 10;
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;

A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find K that makes the control law OPTIMAL
Q = diag([.1         .1      1000000        .1          1000000           .1]);
R = .0001;
[K, S, E] = lqr(A, B, Q, R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define measurement noise characteristics
Udstddev = .1; % no noise for part F, but let's put something here so Ricatti eqn can be solved
Vstddev = .1; % no noise for part F, but let's put something here so Ricatti eqn can be solved

Sigma_D = Udstddev*Udstddev*eye(6); % covariance matrix: no correlation so identity matrix, variance equals std dev squared
Sigma_V = Vstddev*Vstddev*eye(1); % covariance matrix: no correlation so identity matrix, variance equals std dev squared 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find L that makes the system OPTIMALLY DETECTABLE/OBSERVABLE
%%% x(t) is output vector
C_xonly = [1 0 0 0 0 0];

% p_xonly = [-10 -20 -30 -40 -50 -60]; % place very negative, 10 times further left than the usual poles
% L_xonly = place(A',C_xonly',p_xonly).';

% % % Sigma_D = 0.1*eyes(6);
% % % Sigma_V = 1;
L_xonly = (lqr(A',C_xonly',Sigma_D,Sigma_V)).'; % place OPTIMALLY


%%%%%%% Check out how negative the poles were (optional) %%%%%%%%%
eig(A'-(C_xonly')*(L_xonly')  )




tspan = 0:.004:20;

%%%%  x    xdot    q1            q1d           q2           q2d
% x0 = [0;    0;  deg2rad(45);  deg2rad(0);  deg2rad(45);  deg2rad(0)]; % Define initial condition for Xe, the errors


%%%   x(1);  x(2);  x(3);         x(4);          x(5);        x(6)           x(7);  x(8);  x(9);  x(10);  x(11);  x(12)
% xc0 = [0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);        0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);  ]; % Define combined initial condition for X, the state and errors
% xc0 = [0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);        0;    0;  0;  deg2rad(0);  0;  deg2rad(0);  ]; % Define combined initial condition for X, the state and errors

% xc0 = [0;    0;  0;  0;  0;  0;        0;    0;  0;  0;  0;  0;  ]; % Define combined initial condition for X, the state and errors
xc0 = [0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);        0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0);  ]; % Define combined initial condition for X, the state and errors

% u = @(xhat,t) -K*xhat
% u = @(x,t) 1; % as per assignment, force is unit step input. Note control_input (see below) is using ones function right now.

[t,x] = ode45(  @(t,x)crane_diffeq_Nonlinear_Luenberger_fxn_Redo(x,t,K,L_xonly,Udstddev, Vstddev)    ,   tspan,   xc0);

size_x = size(x)
number_timesteps = size_x(1)

% wr = 0;
% control_input = -K*(x(:,7:12)' - wr); % for plotting later, control input was based on estimated states
control_input = ones(size(t));
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

