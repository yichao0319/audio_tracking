clear;

% load models
r1=load('LSModel.mat');
r2=load('NNModelN5.mat');
r3=load('KNNModel.mat');

% set the parameters
Ntest=1000;
Fs=44100;
sampleInterval=0.04;
fmin=17000;
fmax=19500;
vs=346;
Np=2;
results=zeros(4,Ntest);
v=[0,0];
weights=zeros(Ntest,2);
length=zeros(Ntest,2);

% testing
for i=1:Ntest
    w=[1 0.4+0.6*rand];
    R=[1 1.05+0.3*rand];
    weights(i,:)=w;
    length(i,:)=R;
    temp=genFMCWPeakStat(Fs,sampleInterval,fmin,fmax,vs,R,v,w,Np,1);
    % baseline
    % results(1,i)=getSawShapeFMCWDist(Fs,sampleInterval,fmin,fmax,vs,R,v,w,Np,1);
    % linear regression
    results(2,i)=(temp-r1.mu)./r1.sigma*r1.coef+r1.bias;
    % neural network
    results(3,i)=r2.net(transpose((temp-r2.mu)./r2.sigma))+r2.bias; 
    % knn
    idx=knnsearch(r3.features,((temp-r3.mu)./r3.sigma),'K',1);
    results(4,i)=mean(r3.labels(idx));
end

meanErr=mean(abs(results-1),2);
maxErr=[max(abs(results(1,:)-1));max(abs(results(2,:)-1));...
    max(abs(results(3,:)-1));max(abs(results(4,:)-1))];
disp(meanErr);
disp(maxErr);
