% a = 0;
% b = 0;
% r1 = 10;
% r2 = 15;
% d = 6;
% dir_x = 3;
% dir_y = 4;
syms a b r1 r2 d dir_x dir_y
syms x1 y1 x2 y2 t
eq1 = (x1-a)^2 + (y1-b)^2 - r1^2;
eq2 = (x2-a)^2 + (y2-b)^2 - r2^2;
eq3 = (x1-x2)^2 + (y1-y2)^2 - d^2;
eq4 = (x2-x1) - dir_x*t;
eq5 = (y2-y1) - dir_y*t;

ret = solve(eq1,eq2,eq3,eq4,eq5,'x1','y1','x2','y2','t')
ret
ret.x1
ret.y1
ret.x2
ret.y2
ret.t
