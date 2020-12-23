function dx = crane_diffeq_linear_Luenberger_fxn(x,t,u)


F = u(x);

g = 10;
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;

X = [x(1);  x(2);  x(3);  x(4);  x(5);  x(6);  x(7);  x(8);  x(9);  x(10);  x(11);  x(12)]; % this is actually Xe now, the errors

A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -(M+m1)*g/(M*l1) 0 -m2*g/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -(M+m2)*g/(M*l2) 0];
B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];




%%%%%%%%%%%%%%%%%%%%%%%%%% Pick one C matrix of the three below
%%% x(t) is output vector
% C_xonly = [1 0 0 0 0 0];

%%% x(t) theta2(t) is output vector
% C_x_q2 = [1 0   0 0   1 0];

%%% x(t) theta1(t) theta2(t) is output vector
C_x_q1_q2 = [1 0   1 0   1 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Pick one respective observer of the three below to use
Udstddev = .1; % no noise for part F, but let's put something here so Ricatti eqn can be solved
Vstddev = .1; % no noise for part F, but let's put something here so Ricatti eqn can be solved
Sigma_D = Udstddev*Udstddev*eye(6); % covariance matrix: no correlation so identity matrix, variance equals std dev squared
Sigma_V = Vstddev*Vstddev*eye(1); % covariance matrix: no correlation so identity matrix, variance equals std dev squared 

% p_xonly = [-10 -20 -30 -40 -50 -60];
% L_xonly = place(A',C_xonly',p_xonly).';
% L_xonly = (lqr(A',C_xonly',Sigma_D,Sigma_V)).'; % place OPTIMALLY

% p_x_q2 = [-5.1 -5.2 -5.3 -5.4 -5.5 -5.6];
% L_x_q2 = place(A',C_x_q2',p_x_q2).';
% L_x_q2 = (lqr(A',C_x_q2',Sigma_D,Sigma_V)).'; % place OPTIMALLY

% p_x_q1_q2 = [-5.1 -5.2 -5.3 -5.4 -5.5 -5.6];
% L_x_q1_q2 = place(A',C_x_q1_q2',p_x_q1_q2).';
L_x_q1_q2 = (lqr(A',C_x_q1_q2',Sigma_D,Sigma_V)).'; % place OPTIMALLY



%%%%%%%%%%%%%%%%%%%%%%%%%%% Pick one respective observer of the three below to use

% Ac_xonly = [      A                zeros(6,6);
%             L_xonly*C_xonly    A - L_xonly*C_xonly;];
                

% Ac_x_q2 = [      A                zeros(6,6);
%             L_x_q2*C_x_q2    A - L_x_q2*C_x_q2;];


Ac_x_q1_q2 = [      A                zeros(6,6);
            L_x_q1_q2*C_x_q1_q2    A - L_x_q1_q2*C_x_q1_q2;];
        
%%%%%%%%%%%%%%%%%%%%%%%% B is the same for all
Bc = [B; B]; % combined B

%%%%%%%%%%%%%%%%%%%%%%%%% Pick one respective (twelve state) state equation
% Xdot = Ac_xonly * X + Bc*F;
% Xdot = Ac_x_q2 * X + Bc*F;
Xdot = Ac_x_q1_q2 * X + Bc*F;



% Xdot = A*X + B*F;



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