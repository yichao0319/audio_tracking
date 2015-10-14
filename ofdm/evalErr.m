function f=evalErr(z,anchor,n,d)

% z: variables, including power and source position
% anchor: anchor information
% n: the number of sources beside the anchor
% d: training data

Cx=anchor(1);
Cy=anchor(2);
CHeight=anchor(3);

x=d(1);
y=d(2);
c=d(3);

P0=z(1);
P=z(2:3:3*n-1);
sx=z(3:3:3*n);
sy=z(4:3:3*n+1);

dista2=(x-Cx)^2+(y-Cy)^2+CHeight^2;
dist2=(x-sx).^2+(y-sy).^2+CHeight^2;
f=(c-P0/dista2-sum(P./dist2))^2;
