#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>          /* See NOTES */
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <unistd.h>
#include <string.h>

#include <inttypes.h>
#include <math.h>
#include <time.h>

// #include "inc/fmod.h"
// #pragma comment (lib, "lib/libfmod.dylib")

#define SRV_PORT 9930
#define SRV_IP   "10.147.124.83"
#define BUFLEN   1
#define NPACK    100

#ifdef __MACH__
#include <sys/time.h>
#include <mach/mach_time.h>
//clock_gettime is not implemented on OSX
#define CLOCK_REALTIME 0 
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

void diep(const char *str) {
    fprintf(stderr, "\n>> ERROR: %s\n", str);
    exit(1);
}

uint64_t get_current_time_with_ms (void)
{
    uint64_t        ns;
    time_t          s;  // Seconds
    struct timespec spec;

    clock_gettime(CLOCK_REALTIME, &spec);

    s  = spec.tv_sec;
    ns = spec.tv_nsec;
    // printf("Current time: %llu ns since the Epoch\n", ns);
    return ns;
}

// void play_wave() 
// {
//     FSOUND_SAMPLE* handle;
//     FSOUND_Init (44100, 32, 0);

//     FSOUND_Close();
// }

void senderFunc() {
    struct sockaddr_in si_other;
    int s, i;
    int slen = sizeof(si_other);
    char buf[BUFLEN];

    uint64_t        start;
    uint64_t        end;
    uint64_t        elapsed;
    uint64_t        elapsedNano;
    static mach_timebase_info_data_t    sTimebaseInfo;

    uint64_t        start2;
    uint64_t        end2;
    uint64_t        elapsed2;

    if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }
    start = mach_absolute_time();


    if((s=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP)) == -1) {
        diep("socket");
    }

    memset((char *) &si_other, 0, sizeof(si_other));
    si_other.sin_family = AF_INET;
    si_other.sin_port = htons(SRV_PORT);
    if (inet_aton(SRV_IP, &si_other.sin_addr)==0) {
        diep("inet_aton() failed\n");
    }

    for (i=0; i < NPACK; i++) {
        // printf("Sending packet %d\n", i);
        // sprintf(buf, "This is packet %d\n", i);
        // start = mach_absolute_time();
        // start2 = get_current_time_with_ms();
        if (sendto(s, buf, BUFLEN, 0, (struct sockaddr *)&si_other, slen) == -1){
            diep("sendto()");
        }
        // end = mach_absolute_time();
        // end2 = get_current_time_with_ms();
        // elapsed = end - start;
        // elapsedNano = elapsed * sTimebaseInfo.numer / sTimebaseInfo.denom;
        // elapsed2 = end2 - start2;
        // printf("elapsed: %lluns v.s. %lluns\n", elapsedNano, elapsed2);
        // printf("elapsed: %lluns\n", elapsedNano);

        end = mach_absolute_time();
        printf("  time: %llu ns, elapsed=%llu ns\n", end, end-start);
        start = end;

        while(1) {
            if((mach_absolute_time() - start) > 95000) {
                break;
            }
        }
    }

    close(s);
}

int main()
{
    senderFunc();
    return 0;
}
