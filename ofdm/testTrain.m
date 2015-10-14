clear all;

% generate ground truth
anchor=[0,0,0.1];
P0=5;

n=1;
s1=[0;1];
P1=1;
% s2=[0;2];
% P2=1;

% z_true=[P0,P1,s1.',P2,s2.'];
z_true=[P0,P1,s1.'];

% generate data
N=1000;
x=zeros(1,N);
y=zeros(1,N);
c=zeros(1,N);
d=zeros(N,3);
dev=0.01;
for i=1:N
    x(i)=rand*2+1;
    y(i)=rand*2;
    c(i)=evalEDP(z_true,anchor,n,[x(i),y(i)])+dev*randn;
end

z_est=trainGD(anchor,x,y,c,5000,n);
disp(z_true);
disp(z_est);

