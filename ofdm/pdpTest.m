clear all;

% parameters
Fs=44100;
Nfft=128;
Ncp=32;
fc=18000;
nSamp=20;
Ns=50;

% get channel gains and distances
Np=9;
% posiX=0.88+(0:0.1:0.7);
% posiY=0.48;
posiX=0.41:0.1:1.21;
posiY=zeros(1,Np);
p=[posiX;posiY];
CHeight=0.095;
C1=[0;0];

% fileIndex=1435875000+[620,666,715,763,836,878,960,1037];
% fileIndex=1436191000+[1324,1274,1218,800,746,692,611,561,511];
fileIndex=1436194000+[166,215,267,327,384,431,478,538,592];
h1=cell(1,Np);
h1Mag= cell(1,Np);
d=zeros(1,Np);
h1MagMean=zeros(1,Np);
for i=1:Np
    h1{i}=ofdm_rx(int2str(fileIndex(i)),Fs,Nfft,Ncp,fc,nSamp,Ns);
    h1Mag{i}=abs(h1{i}).^2;
    h1MagMean(i)=mean(h1Mag{i});
    d(i)=normCorr(CHeight,p(:,i)-C1);
end
res=[d;h1MagMean];

% % get channel gains and distances
% fileIndex=1435874000+[-57,27,84,141,204,291,379,533,606,664,708,759,856,935,1002];
% Np=numel(fileIndex);
% h1=cell(1,Np);
% h1Mag= cell(1,Np);
% p=zeros(2,Np);
% d=zeros(1,Np);
% h1MagMean=zeros(1,Np);
% C1=[0;0];
% CHeight=0.085;
% for i=1:Np
%     [h1{i},p(:,i)]=ofdm_rx(int2str(fileIndex(i)),Fs,Nfft,Ncp,fc,nSamp,Ns);
%     h1Mag{i}=abs(h1{i}).^2;
%     h1MagMean(i)=mean(h1Mag{i});
%     d(i)=normCorr(CHeight,p(:,i)-C1);
% end
% res=[d;h1MagMean];


