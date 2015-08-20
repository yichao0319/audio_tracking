Fs = 44100;
r = audiorecorder(Fs,24,1);
recordblocking(r,5);
rcv_packet = getaudiodata(r,'int16');
%plot(rcv_packet);
wavwrite(rcv_packet,Fs,'chirp_record'); 

