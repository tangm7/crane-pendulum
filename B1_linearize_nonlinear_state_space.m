clc, clear;

% Do not edit this
syms Atransformation(t,alpha,r,d); % Robogrok convention

Atransformation = [cos(t) -sin(t)*cos(alpha) sin(t)*sin(alpha) r*cos(t); 
    sin(t) cos(t)*cos(alpha) -cos(t)*sin(alpha) r*sin(t);0 sin(alpha) cos(alpha) d; 0 0 0 1 ]; % Robogrok, r

%MT
syms q1(t) q2(t) x(t)
syms tau1 tau2 F
syms M m1 m2 l1 l2 g 

% plug in values of 1st row of DH table
s=struct('t',0,       'alpha',0,        'r',x(t),        'd',0); % Robogrok
A1 = subs(Atransformation,s); % T 0->1 matrix

%2nd row of DH table
s=struct('t',-pi/2-q1(t),       'alpha',0,        'r',l1,        'd',0); % Robogrok
A2 = subs(Atransformation,s); % T 1->2 matrix

%3rd row of DH table
s=struct('t',-pi/2-q2(t),       'alpha',0,        'r',l2,        'd',0); % Robogrok
A3 = subs(Atransformation,s); % T 2->3 matrix


T0a=simplify(A1)
Ta1=simplify(A2)
T01=simplify(T0a*Ta1)

Ta2=simplify(A3)
T02=simplify(T0a*Ta2)



z0a = T0a(1:3,3)

o0a = T0a(1:3,4)
o01 = T01(1:3,4)
o02 = T02(1:3,4)

Jva = [   [1;0;0]    zeros(3,1)    zeros(3,1)]
Jwa = [ zeros(3,1)   zeros(3,1)    zeros(3,1)]

Jv1 = [   [1;0;0]     cross(z0a,o01-o0a)    zeros(3,1)   ]
Jw1 = [ zeros(3,1)        z0a               zeros(3,1)   ]

Jv2 = [   [1;0;0]     zeros(3,1)         cross(z0a,o02-o0a)       ] % fixed
Jw2 = [   zeros(3,1)        zeros(3,1)            z0a    ] % fixed


qdot = [diff(x(t),t); -diff(q1(t),t); -diff(q2(t),t) ]

v0a = Jva * qdot
v01 = Jv1 * qdot
v02 = Jv2 * qdot

w0a = Jwa * qdot
w01 = Jw1 * qdot
w02 = Jw2 * qdot

R01 = T01(1:3,1:3)
R02 = T02(1:3,1:3)

I1 = zeros(3,3)
I2 = zeros(3,3)



Ta = 1/2 * M * (v0a.') * v0a

T1 = 1/2 * m1 * (v01.') * v01 + 1/2 * (w01.') * R01*I1*(R01.') * w01
T2 = 1/2 * m2 * (v02.') * v02 + 1/2 * (w02.') * R02*I2*(R02.') * w02


T = simplify( Ta + T1 + T2 )



G = [0;g;0]

Va = 0
V1 = m1*(G.')*o01 + m1*(G.')*l1*[0;1;0]
V2 = m2*(G.')*o02 + m2*(G.')*l2*[0;1;0]

V = Va + V1 + V2



L = T - V

eqn1 = (    collect(    simplify(    diff( diff(L, diff(x(t), t)), t)  - diff(L, x(t))      ), [q1(t) diff(q1(t), t) diff(q1(t), t, t) q2(t) diff(q2(t), t) diff(q2(t), t, t)    x(t)  diff(x(t),t)   diff(x(t),t,t)   ]    )    == F       )
eqn2 = (    collect(    simplify(    diff( diff(L, diff(q1(t), t)), t) - diff(L, q1(t))     ), [q1(t) diff(q1(t), t) diff(q1(t), t, t) q2(t) diff(q2(t), t) diff(q2(t), t, t)    x(t)  diff(x(t),t)   diff(x(t),t,t)   ]    )    == 0    )
eqn3 = (    collect(    simplify(    diff( diff(L, diff(q2(t), t)), t) - diff(L, q2(t))     ), [q1(t) diff(q1(t), t) diff(q1(t), t, t) q2(t) diff(q2(t), t) diff(q2(t), t, t)    x(t)  diff(x(t),t)   diff(x(t),t,t)   ]    )    == 0    )

display("equations of motion")
latex(eqn1)
latex(eqn2)
latex(eqn3)






% Lsimple = (    collect(    simplify(    L     ), [1/2*m1 1/2*m2 q1(t) diff(q1(t), t) diff(q1(t), t, t) q2(t) diff(q2(t), t) diff(q2(t), t, t)    x(t)  diff(x(t),t)   diff(x(t),t,t)   ]    )    == tau2    )
% 
% latex(Lsimple)


disp("latex equations of energy expressions:")
latex(Ta)
latex(T1)
latex(T2)

% latex(Va)
latex(V1)
latex(V2)



T2
T2_ramzi = 1/2*m2*diff(x(t),t)^2 - m2*l2*diff(x(t),t)*diff(q2(t),t)*cos(q2(t)) + 1/2*m2*l2^2 * diff(q2(t),t)^2

L_ramzi_origin = 1/2 * M*diff(x(t),t)^2 + 1/2*m1*diff(x(t),t)^2 + 1/2*m1*l1^2*diff(q1(t),t)^2  - m1*l1*diff(x(t),t) * diff(q1(t),t) * cos(q1(t)) + 1/2 * m2 * diff(x(t),t)^2 +  1/2 * m2 * l2^2 * diff(q2(t),t)^2 - m2 * l2 * diff(x(t),t) * diff(q2(t),t) * cos(q2(t)) - m1*g*l1*( 1-cos(q1(t)) ) - m2*g*l2*( 1- cos(q2(t)))
L_ramzi_retype = 1/2 * M*diff(x(t),t)^2 + 1/2*m1*diff(x(t),t)^2 + 1/2*m1*l1^2*diff(q1(t),t)^2  - m1*l1*diff(x(t),t) * diff(q1(t),t) * cos(q1(t)) + 1/2 * m2 * diff(x(t),t)^2 +  1/2 * m2 * l2^2 * diff(q2(t),t)^2 - m2 * l2 * diff(x(t),t) * diff(q2(t),t) * cos(q2(t)) - m1*g*l1*( 1-cos(q1(t)) ) - m2*g*l2*( 1- cos(q2(t)))

display("checked ramzi equations were typed correctly")
L_ramzi_origin-L_ramzi_retype

simplify(T2-T2_ramzi)


double(   subs(T1, [q1(t) diff(q1(t), t) diff(q1(t), t, t) q2(t) diff(q2(t), t) diff(q2(t), t, t)    x(t)  diff(x(t),t)   diff(x(t),t,t)   m1  m2   l1  l2], [-1.5708  .1 0   0 0 0   0 0 0   1 1   1 1])   )


latex(simplify(L))

display("compare L expressions")
simplify( L - L_ramzi_origin)




% display("convert equations of motion to thetas")
% 
% string_eqn1 = latex(eqn1)
% string_eqn2 = latex(eqn2)
% string_eqn3 = latex(eqn3)
% string_L    = latex(simplify(L))
% 
% oldarray = ['q_{1}','q_{2}']
% newarray = ['\\theta_1','\\theta_2']
% 
% newChr = strrep(string_eqn1,oldarray(1),newarray(1))

% replaced_eqn1 = regexprep(string_eqn1,{'q_{1}','q_{2}'},{'\\theta_1','\\theta_2'})

% replaced_eqnasdf = regexprep('qwer{1}',{'qwer{1}','q_{2}'},{'asdfasdf','\\theta_2'})

display("Now we need to solve for our second derivative terms, 3 equations, 3 unknowns. We can use matrix for this system of equations")

z = [diff(x(t),t,t);   diff(q1(t),t,t);    diff(q2(t),t,t)]


G = [ (M+m1+m2)           -l1*m1*cos(q1(t))    -l2*m2*cos(q2(t));
    -l1*m1*cos(q1(t))     l1^2*m1                     0         ;
    -l2*m2*cos(q2(t))           0                     l2^2*m2   ;]



b = [ -l1*m1*sin(q1(t))*diff(q1(t),t)^2 - l2*m2*sin(q2(t))*diff(q2(t),t)^2 + F ;
%       -l1*m1*sin(q1(t))*diff(q1(t),t)^2 - l2*m2*sin(q2(t))*diff(q2(t),t)^2 + F
    -g*l1*m1*sin(q1(t));
    -g*l2*m2*sin(q2(t));]

matrixeqns = G*z - b == 0
matrixeqnsLHS = G*z - b

eqn1_matrix = matrixeqns(1)
eqn2_matrix = matrixeqns(2)
eqn3_matrix = matrixeqns(3)

display("check we have typed in matrix equations correctly")
simplify( (eqn1-F) - eqn1_matrix)
isAlways(      (  diff( diff(L, diff(x(t), t)), t)  - diff(L, x(t)) - F  )      ==     matrixeqnsLHS(1)   ) 
simplify(      (  diff( diff(L, diff(x(t), t)), t)  - diff(L, x(t)) - F  )       -     matrixeqnsLHS(1)   )
simplify(eqn2 - eqn2_matrix)
simplify(eqn3 - eqn3_matrix)



second_derivatives = simplify( G\b ) % solve for second derivatives


latex(second_derivatives)


display("next we need to factor out the expressions for the second derivatives into parts multiplied by each state variable")

Ad = -(   -m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2   ); % note negative in front of whole thing
A24 = l1*m1*sin(q1(t))*diff(q1(t),t);
A26 = l2*m2*sin(q2(t))*diff(q2(t),t);

% E2 = g*m1*cos(q1(t))*sin(q1(t)) + g*m2*cos(q2(t))*sin(q2(t))   - F;
E2 = g*m1*cos(q1(t))*sin(q1(t)) + g*m2*cos(q2(t))*sin(q2(t));

B2 = -1;


x_doubledot_matrix_style = 1/Ad * (  A24 * diff(q1(t),t)   + A26*diff(q2(t),t)    + E2  +    B2*F    );

simplify(   x_doubledot_matrix_style - second_derivatives(1)   ) % verify xdoubledot equation is correct


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A44 = l1*m1*cos(q1(t))*sin(q1(t))*diff(q1(t),t)
A46 = l2*m2*cos(q1(t))*sin(q2(t))*diff(q2(t),t)
% E4 = M*g*sin(q1(t)) - F*cos(q1(t)) + g*m1*sin(q1(t)) + g*m2*sin(q1(t)) - g*m2*cos(q2(t))^2*sin(q1(t)) +  g*m2*cos(q1(t)) * cos(q2(t)) * sin(q2(t))
E4 = M*g*sin(q1(t)) + g*m1*sin(q1(t)) + g*m2*sin(q1(t)) - g*m2*cos(q2(t))^2*sin(q1(t)) +  g*m2*cos(q1(t)) * cos(q2(t)) * sin(q2(t))

B4 = -cos(q1(t))

theta1_doubledot_matrix_style = 1/(Ad*l1) * (E4      +  A44*diff(q1(t))      +   A46*diff(q2(t),t)     +     B4*F     )

simplify(   theta1_doubledot_matrix_style - second_derivatives(2)   ) % verify matrix style theta1_doubledot equation is correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A64 = l1*m1*cos(q2(t))*sin(q1(t))*diff(q1(t),t);
A66 = l2*m2*cos(q2(t))*sin(q2(t))*diff(q2(t),t);
% E6 = M*g*sin(q2(t)) - F*cos(q2(t)) + g*m1*sin(q2(t)) + g*m2*sin(q2(t)) - g*m1*cos(q1(t))^2*sin(q2(t)) + g*m1*cos(q1(t))*cos(q2(t))*sin(q1(t));
E6 = M*g*sin(q2(t)) + g*m1*sin(q2(t)) + g*m2*sin(q2(t)) - g*m1*cos(q1(t))^2*sin(q2(t)) + g*m1*cos(q1(t))*cos(q2(t))*sin(q1(t));

B6 = -cos(q2(t))

theta2_doubledot_matrix_style = 1/(Ad*l2) * (E6       + A64 * diff(q1(t),t)      +  A66*diff(q2(t),t)    +    B6*F        );

simplify( theta2_doubledot_matrix_style - second_derivatives(3) )% verify matrix style theta2_doubledot equation is correct




A = [0    1    0        0          0       0;
    0    0    0    1/Ad*A24        0    1/Ad*A26;
    0    0    0         1          0       0;
    0    0    0    1/(Ad*l1)*A44   0    1/(Ad*l1)*A46;
    0    0    0         0          0       1;
    0    0    0    1/(Ad*l2)*A64   0    1/(Ad*l2)*A66;]

B = [0;   1/Ad*B2;   0;   1/(Ad*l1)*B4;   0;   1/(Ad*l2)*B6]
u = F

E = [0;   1/Ad*E2;   0;   1/(Ad*l1)*E4;    0;   1/(Ad*l2)*E6]

X = [x(t);   diff(x(t),t);   q1(t);   diff(q1(t),t);   q2(t);   diff(q2(t),t)]
X_dot = A * X   +   B*u   +   E


display("now let's verify the nonlinear state space equations")


simplify( X_dot(2) - second_derivatives(1)  )
simplify( X_dot(4) - second_derivatives(2)  )
simplify( X_dot(6) - second_derivatives(3)  )

latex( [ diff(x(t),t);   diff(x(t),t,t);   diff(q1(t),t);   diff(q1(t),t,t);   diff(q2(t),t);   diff(q2(t),t,t)] )

latex(simplify(X_dot))

latex(A)

latex(X)

latex(B)

latex(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display("Now all we need to do is linearize the nonlinear state space equations")


f1 = X_dot(1)
f2 = X_dot(2)
f3 = X_dot(3)
f4 = X_dot(4)
f5 = X_dot(5)
f6 = X_dot(6)


x_e = 0;
xd_e = 0;
xdd_e = 0; % shouldn't need this

q1_e = 0;
q1d_e = 0;
q1dd_e = 0; % shouldn't need this

q2_e = 0;
q2d_e = 0;
q2dd_e = 0; % shouldn't need this

A_linear = subs(   [   diff(f1,X(1))    diff(f1,X(2))    diff(f1,X(3))    diff(f1,X(4))    diff(f1,X(5))    diff(f1,X(6));
                       diff(f2,X(1))    diff(f2,X(2))    diff(f2,X(3))    diff(f2,X(4))    diff(f2,X(5))    diff(f2,X(6));
                       diff(f3,X(1))    diff(f3,X(2))    diff(f3,X(3))    diff(f3,X(4))    diff(f3,X(5))    diff(f3,X(6));
                       diff(f4,X(1))    diff(f4,X(2))    diff(f4,X(3))    diff(f4,X(4))    diff(f4,X(5))    diff(f4,X(6));
                       diff(f5,X(1))    diff(f5,X(2))    diff(f5,X(3))    diff(f5,X(4))    diff(f5,X(5))    diff(f5,X(6));
                       diff(f6,X(1))    diff(f6,X(2))    diff(f6,X(3))    diff(f6,X(4))    diff(f6,X(5))    diff(f6,X(6));]  , ... 
           [q1(t) diff(q1(t), t) diff(q1(t), t, t) q2(t) diff(q2(t), t) diff(q2(t), t, t)    x(t)  diff(x(t),t)   diff(x(t),t,t)   ], ... % old
           [q1_e q1d_e q1dd_e q2_e q2d_e q2dd_e    x_e  xd_e   xdd_e   ]   )  % new
   

B_linear = subs(    [  diff(f1,F);
                       diff(f2,F);
                       diff(f3,F);
                       diff(f4,F);
                       diff(f5,F);
                       diff(f6,F); ]   , ...
           [q1(t) diff(q1(t), t) diff(q1(t), t, t) q2(t) diff(q2(t), t) diff(q2(t), t, t)    x(t)  diff(x(t),t)   diff(x(t),t,t)   ], ... % old
           [q1_e q1d_e q1dd_e q2_e q2d_e q2dd_e    x_e  xd_e   xdd_e   ]   )  % new


    
latex(A_linear)
latex(B_linear)


