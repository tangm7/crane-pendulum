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

P_1 = [p11 p12 p13 p14 p15 p16; p21 p22 p23 p24 p25 p26; p31 p32 p33 p34 p35 p36; p41 p42 p43 p44 p45 p46; p51 p52 p53 p54 p55 p56; p61 p62 p63 p64 p65 p66];

Qall = 10000;
% Q = diag([Qall Qall Qall Qall Qall Qall]);
% Q = diag([100 Qall Qall Qall Qall Qall]);

%%%%      x       xdot      q1         q1d           q2              q2d
Q = diag([.1         .1      1000000        .1          1000000           .1]);

% Q = diag([100 100 100 100 100 100]);
% Q = [1 0 0 0 0 0; 
%     0 1 0 0 0 0; 
%     0 0 1 0 0 0; 
%     0 0 0 1 0 0; 
%     0 0 0 0 1 0; 
%     0 0 0 0 0 1];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K, S, E] = lqr(A, B, Q, R);

%%%%%%%%%%
C = [1 1 1 1 1 1]; % perfect sensors, linear state feedback
D = 0; % no feedforward terms

tspan = 0:.004:20;

%%%%  x    xdot    q1            q1d           q2           q2d
x0 = [0;    0;  deg2rad(15);  deg2rad(0);  deg2rad(15);  deg2rad(0)];

u = @(x,t) -K*x

[t,x] = ode45(  @(t,x)crane_diffeq_nonlinear_fxn(x,t,u)    ,   tspan,   x0);

wr = 0;
control_input = -K*(x' - wr); % for plotting later

% u = @(x,t) K*x;
% sys = ss(A,B,C,D)
% [Kbuiltin,S,e] = lqr(sys,Q,R)
% lsim(sys,u,t,x0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot results

figure

subplot(7,1,1) % Plot Cart Position
title(['Scenario ', 1])
plot(t, (x(:,1)), 'LineWidth', 3)
grid on
legend('x (cart position)')
xlabel('Time (s)') 
ylabel('(m)')



subplot(7,1,2) % Plot Cart Velocity
plot(t, (x(:,2)), 'LineWidth', 3)
grid on
legend('xdot (cart velocity)')
xlabel('Time (s)') 
ylabel('(m/s)') 



subplot(7,1,3) % Plot pendulum 1 angle
% plot(t, (x(:,3)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,3))  ), 'LineWidth', 3)
grid on
legend('theta1 (pendulum 1 angle)')
xlabel('Time (s)') 
% ylabel('(rad)')
ylabel('(deg)')




subplot(7,1,4) % Plot pendulum 1 angular velocity
% plot(t, (x(:,4)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,4))  ), 'LineWidth', 3)
grid on
legend('theta1dot (pendulum 1 angular velocity)')
xlabel('Time (s)') 
% ylabel('(rad/s)') 
ylabel('(deg/s)') 


subplot(7,1,5) % Plot pendulum 2 position
% plot(t, (x(:,5)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,5))  ), 'LineWidth', 3)
grid on
legend('theta2 (pendulum 2 position)')
xlabel('Time (s)') 
% ylabel('(rad)') 
ylabel('(deg)') 



subplot(7,1,6) % Plot pendulum 2 angular velocity
% plot(t, (x(:,6)), 'LineWidth', 3)
plot(t, rad2deg(  (x(:,6))  ), 'LineWidth', 3)
grid on
legend('theta2dot (pendulum 2 angular velocity)')
xlabel('Time (s)') 
% ylabel('(rad/s)') 
ylabel('(deg/s)') 


subplot(7,1,7) % Plot Cart Position
plot(t, (control_input), 'LineWidth', 3)
grid on
legend('u: control_input')
xlabel('Time (s)') 
ylabel('(N)')


