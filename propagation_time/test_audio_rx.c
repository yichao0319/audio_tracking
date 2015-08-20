#include "AudioToolbox/AudioToolbox.h"

static const int kNumberBuffers = 3;                            

typedef struct AQRecorderState {
    AudioStreamBasicDescription  mDataFormat;
    AudioQueueRef                mQueue;
    AudioQueueBufferRef          mBuffers[kNumberBuffers];
    AudioFileID                  mAudioFile;
    UInt32                       bufferByteSize;
    SInt64                       mCurrentPacket;
    bool                         mIsRunning;
} AQRecorderState;


/* call back function */
static void HandleInputBuffer (
    void                                *aqData,             // 1
    AudioQueueRef                       inAQ,                // 2
    AudioQueueBufferRef                 inBuffer,            // 3
    const AudioTimeStamp                *inStartTime,        // 4
    UInt32                              inNumPackets,        // 5
    const AudioStreamPacketDescription  *inPacketDesc        // 6
) {
    AQRecorderState *pAqData = (AQRecorderState *)aqData;

    if (inNumPackets == 0 && pAqData->mDataFormat.mBytesPerPacket != 0) {
        inNumPackets = inBuffer->mAudioDataByteSize / pAqData->mDataFormat.mBytesPerPacket;
    }

    // Writing an Audio Queue Buffer to Disk
    if( AudioFileWritePackets (
            pAqData->mAudioFile,
            false,
            inBuffer->mAudioDataByteSize,
            inPacketDesc,
            pAqData->mCurrentPacket,
            &inNumPackets,
            inBuffer->mAudioData) == noErr) {

        pAqData->mCurrentPacket += inNumPackets;
    }

    if (pAqData->mIsRunning == 0) {
        return;
    }

    // Enqueuing an Audio Queue Buffer
    AudioQueueEnqueueBuffer (
        pAqData->mQueue,
        inBuffer,
        0,          // not used for recording
        NULL        // not used for recording
    );

}


/* Write a Function to Derive Recording Audio Queue Buffer Size */
void DeriveBufferSize (
    AudioQueueRef                audioQueue,
    AudioStreamBasicDescription  ASBDescription,
    Float64                      seconds,
    UInt32                       *outBufferSize
) {
    static const int maxBufferSize = 0x50000;   // upper bound
    int maxPacketSize = ASBDescription.mBytesPerPacket; // for CBR

    // for VBR
    if (maxPacketSize == 0) {
        UInt32 maxVBRPacketSize = sizeof(maxPacketSize);
        AudioQueueGetProperty (
                audioQueue,
                kAudioQueueProperty_MaximumOutputPacketSize,
                // in Mac OS X v10.5, instead use
                //   kAudioConverterPropertyMaximumOutputPacketSize
                &maxPacketSize,
                &maxVBRPacketSize
        );
    }

 
    Float64 numBytesForTime = ASBDescription.mSampleRate * maxPacketSize * seconds;

    UInt32 tmp = numBytesForTime < maxBufferSize ? numBytesForTime : maxBufferSize;
    *outBufferSize = tmp;

}


OSStatus SetMagicCookieForFile (
    AudioQueueRef inQueue,
    AudioFileID   inFile
) {
    OSStatus result = noErr;
    UInt32 cookieSize;


    if(AudioQueueGetPropertySize (
        inQueue,
        kAudioQueueProperty_MagicCookie,
        &cookieSize ) == noErr
    ) {
        char* magicCookie = (char *) malloc (cookieSize);
        if(AudioQueueGetProperty (
            inQueue,
            kAudioQueueProperty_MagicCookie,
            magicCookie,
            &cookieSize) == noErr
        ) {
            result = AudioFileSetProperty (
                        inFile,
                        kAudioFilePropertyMagicCookieData,
                        cookieSize,
                        magicCookie
                    );
        }
        free (magicCookie);
    }
    return result;
}


int main(int argc, char const *argv[])
{
    char filePath[100] = "./test.wav";

    // Set Up an Audio Format for Recording
    AQRecorderState aqData;
    aqData.mDataFormat.mFormatID         = kAudioFormatLinearPCM;
    aqData.mDataFormat.mSampleRate       = 44100.0;
    aqData.mDataFormat.mChannelsPerFrame = 1;
    aqData.mDataFormat.mBitsPerChannel   = 16;
    aqData.mDataFormat.mBytesPerPacket   =
        aqData.mDataFormat.mBytesPerFrame = 2;
        // aqData.mDataFormat.mChannelsPerFrame * sizeof (SInt16);
    printf(">> %ld\n", sizeof(SInt16));
    aqData.mDataFormat.mFramesPerPacket  = 1;
    AudioFileTypeID fileType             = kAudioFileAIFFType;
    aqData.mDataFormat.mFormatFlags =
        kLinearPCMFormatFlagIsBigEndian
        | kLinearPCMFormatFlagIsSignedInteger
        | kLinearPCMFormatFlagIsPacked;


    // Create a Recording Audio Queue
    AudioQueueNewInput (
        &aqData.mDataFormat,
        HandleInputBuffer,
        &aqData,
        NULL,
        kCFRunLoopCommonModes,
        0,
        &aqData.mQueue
    );


    // Getting the Full Audio Format from an Audio Queue
    UInt32 dataFormatSize = sizeof (aqData.mDataFormat);
    AudioQueueGetProperty (
        aqData.mQueue,
        kAudioQueueProperty_StreamDescription,
        &aqData.mDataFormat,
        &dataFormatSize
    );


    // Create an Audio File
    CFURLRef audioFileURL =
        CFURLCreateFromFileSystemRepresentation (
            NULL,
            (const UInt8 *) filePath,
            strlen (filePath),
            false
        );

    AudioFileCreateWithURL (
        audioFileURL,
        fileType,
        &aqData.mDataFormat,
        kAudioFileFlags_EraseFile,
        &aqData.mAudioFile
    );


    // Set an Audio Queue Buffer Size
    DeriveBufferSize (
        aqData.mQueue,
        aqData.mDataFormat,
        0.5,
        &aqData.bufferByteSize
    );


    // Prepare a Set of Audio Queue Buffers
    for (int i = 0; i < kNumberBuffers; ++i) {
        AudioQueueAllocateBuffer (                       // 2
            aqData.mQueue,
            aqData.bufferByteSize,
            &aqData.mBuffers[i]
        );

        AudioQueueEnqueueBuffer (
            aqData.mQueue,
            aqData.mBuffers[i],
            0,
            NULL
        );
    }


    // Record Audio
    aqData.mCurrentPacket = 0;
    aqData.mIsRunning = true;
    AudioQueueStart (
        aqData.mQueue,
        NULL    //  NULL to indicate that the audio queue should start recording immediately.
    );

    
    // Wait, on user interface thread, until user stops the recording
    char tmp;
    scanf("%c", &tmp);

    AudioQueueStop (
        aqData.mQueue,
        true
    );
    aqData.mIsRunning = false;


    // Clean Up After Recording
    AudioQueueDispose (
        aqData.mQueue,
        true
    );
    AudioFileClose (aqData.mAudioFile);


    return 0;
}
