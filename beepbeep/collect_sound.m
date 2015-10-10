Fs = 44100;
r = audiorecorder(Fs, 16, 1);
recordblocking(r, 20);
rcv_packet = getaudiodata(r, 'int16');
plot(rcv_packet);
wavwrite(rcv_packet, Fs, '0901.exp5.pc1.wav');
