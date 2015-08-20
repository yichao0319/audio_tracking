CC=clang++
CFLAGS=#-c -Wall
LDFLAGS=`pkg-config --cflags --libs opencv`

ALL=camera_calibration camera_tracking

all: $(ALL)

camera_calibration: camera_calibration.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

camera_calibration.o: camera_calibration.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c $<

camera_tracking: camera_tracking.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

camera_tracking.o: camera_tracking.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c $<

clean:
	rm -rf *.o core* $(ALL)
