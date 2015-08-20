#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>          /* See NOTES */
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <unistd.h>
#include <string.h>

#include <sys/time.h>
#ifdef __MACH__
#include <mach/mach_time.h>
#endif

#include "AudioToolbox/AudioToolbox.h"

#define BUFLEN 512
#define NPACK 100
#define PORT 9930

void diep(const char *str) {
    fprintf(stderr, "\n>> ERROR: %s\n", str);
    exit(1);
}

int main(void)
{
    struct sockaddr_in si_me, si_other;
    int s, i;
    int slen = sizeof(si_other);
    char buf[BUFLEN];

    uint64_t        start;
    uint64_t        end;
    uint64_t        elapsed;
    uint64_t        elapsedNano;
    static mach_timebase_info_data_t    sTimebaseInfo;
    if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }
    start = mach_absolute_time();
    


    if((s=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1) {
        diep("socket");
    }

    memset((char *) &si_me, 0, sizeof(si_me));
    si_me.sin_family = AF_INET;
    si_me.sin_port = htons(PORT);
    si_me.sin_addr.s_addr = htonl(INADDR_ANY);
    if (bind(s, (struct sockaddr *)&si_me, sizeof(si_me))==-1) {
        diep("bind");
    }

    for (i = 0; i < NPACK; i ++) {
        if (recvfrom(s, buf, BUFLEN, 0, (struct sockaddr *)&si_other, (unsigned int *)&slen)==-1) {
            diep("recvfrom()");
        }
        end = mach_absolute_time();
        // printf("Received packet from %s:%d\nData: %s\n", 
        //     inet_ntoa(si_other.sin_addr), ntohs(si_other.sin_port), buf);
        printf("  time: %llu ns, elapsed=%llu ns\n", end, end-start);
        start = end;
    }
    close(s);
    return 0;
}

