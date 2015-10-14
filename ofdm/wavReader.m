fileIndex='1436305545';
fileSuffix='_point_0_audio_audio.wav';
fileName=strcat(fileIndex,fileSuffix);
[y,Fs]=audioread(fileName);

Ly=numel(y);
plot((1:Ly),abs(y))

% % modify the wav file
% fileIndex='1435875763m';
% fileSuffix='_point_0_audio_audio.wav';
% fileName=strcat(fileIndex,fileSuffix);
% y=y(1e5:numel(y));
% audiowrite(fileName,y,Fs,'BitsPerSample', 16);

% downconvert
T=numel(y);
fc=18000;
Ts=1/Fs;
y=y.'.*exp(-1i*2*pi*fc*(1:T)*Ts);
y=lowPassFilterByFFT(y,Fs,2000,500);
plot((1:Ly),abs(y));