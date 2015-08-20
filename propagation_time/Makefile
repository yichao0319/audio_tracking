CC=gcc
CFLAGS=#-c -Wall
AUDIOFLAGS=-Wimplicit -framework AudioToolbox -framework CoreFoundation -lm

all: tx rx test_audio_tx test_audio_rx 

tx: sender.c
	$(CC) $(CFLAGS) -o $@ $<

rx: receiver.c
	$(CC) $(CFLAGS) $(AUDIOFLAGS) -o $@ $< 

test_audio_tx: test_audio_tx.c
	$(CC) $(CFLAGS) $(AUDIOFLAGS) -o $@ $< 

test_audio_rx: test_audio_rx.c
	$(CC) $(CFLAGS) $(AUDIOFLAGS) -o $@ $< 

clean:
	rm tx rx test_audio_tx test_audio_rx
