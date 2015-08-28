Fs = 44100;
r = audiorecorder(Fs, 16, 1);
recordblocking(r, 3);
rcv_packet = getaudiodata(r, 'int16');
plot(rcv_packet);
wavwrite(rcv_packet, Fs, '0928.pc1.exp1.wav');
