clear all;

anchor=[0,0,0.1];
z=[6,2,0.2,1.4,1,0.2,2.0];
n=2;
d=[1,1,5];
f0=evalErr(z,anchor,n,d);
disp(f0);
f2=(d(3)-evalEDP(z,anchor,n,d))^2;
disp(f2);

f=zeros(1,3*n+1);
delta=1e-8;
for i=1:3*n+1
    z1=z;
    z1(i)=z1(i)+delta;
    f(i)=evalErr(z1,anchor,n,d);    
end
gest=(f-f0)/delta;
disp(gest);

g=calcGradient(z,anchor,n,d);
disp(g);