clear, clc;
syms q1(t) q2(t) x(t)
syms tau1 tau2 F
syms M m1 m2 l1 l2 g 
A = ...
[0, 1, 0,                                                                                                      0, 0,                                                                                                      0;
0, 0, 0,                 -(l1*m1*sin(q1(t))*diff(q1(t), t))/(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2), 0,                 -(l2*m2*sin(q2(t))*diff(q2(t), t))/(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2);
0, 0, 0,                                                                                                      1, 0,                                                                                                      0;
0, 0, 0,         -(m1*cos(q1(t))*sin(q1(t))*diff(q1(t), t))/(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2), 0, -(l2*m2*cos(q1(t))*sin(q2(t))*diff(q2(t), t))/(l1*(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2));
0, 0, 0,                                                                                                      0, 0,                                                                                                      1;
0, 0, 0, -(l1*m1*cos(q2(t))*sin(q1(t))*diff(q1(t), t))/(l2*(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2)), 0,         -(m2*cos(q2(t))*sin(q2(t))*diff(q2(t), t))/(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2)]




B = [...
                                                                  0;
              1/(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2);
                                                                  0;
cos(q1(t))/(l1*(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2));
                                                                  0;
cos(q2(t))/(l2*(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2));]




 E = [...
 
                                                                                                                                                                                      0;
                                                                         -(g*m1*cos(q1(t))*sin(q1(t)) + g*m2*cos(q2(t))*sin(q2(t)))/(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2);
                                                                                                                                                                                      0;
-(- g*m2*sin(q1(t))*cos(q2(t))^2 + g*m2*cos(q1(t))*sin(q2(t))*cos(q2(t)) + M*g*sin(q1(t)) + g*m1*sin(q1(t)) + g*m2*sin(q1(t)))/(l1*(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2));
                                                                                                                                                                                      0;
-(- g*m1*sin(q2(t))*cos(q1(t))^2 + g*m1*cos(q2(t))*sin(q1(t))*cos(q1(t)) + M*g*sin(q2(t)) + g*m1*sin(q2(t)) + g*m2*sin(q2(t)))/(l2*(- m1*cos(q1(t))^2 - m2*cos(q2(t))^2 + M + m1 + m2))]