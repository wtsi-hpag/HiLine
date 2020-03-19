/*
Copyright (c) 2020 Ed Harry, Wellcome Sanger Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#define SAM_Buffer_Size 3925
#define Number_Of_SAM_Buffers_Per_Queue 16
#define Number_Of_SAM_Buffer_Queues 16

enum
read_type
{
    typeUnpaired,
    typeUnmapped,
    typeSecondary,
    typeQCFail,
    typeDuplicate,
    typeSupplementary,
    typeInvalidReference,
    typeBelowMinMapQ,
    typeTooFarFromRestrictionSite,
    typeDump,
    typeSelfCircle,
    typeDanglingEnd,
    typeSameFragmentAndStrand,
    typeReLigated,
    typeFF,
    typeFR,
    typeRF,
    typeRR
};

struct
sam_buffer
{
    u08 sam[SAM_Buffer_Size];
    u08 isGood;
    u08 isReverse;
    u08 isRead1;
    u32 homeIndex;
    u32 lineLength;
    s32 distanceFromRestrictionSite;
    read_type type;
    u64 samIndex;
    sequence *sequence;
    wavl_node *fragment;
    u64 readStart;
    process_queue *processQueue;
    process_queue *writeQueue;
    sam_buffer *prev;
};

struct
single_sam_buffer_queue
{
    u32 queueLength;
    u32 pad;
    mutex rwMutex;
    sam_buffer *front;
    sam_buffer *rear;
};

struct
sam_buffer_queue
{
    single_sam_buffer_queue **queues;
    threadSig index;
    u32 pad;
};

global_variable
struct sam_buffer_queue *
SAM_Buffer_Queue;

global_function
void
InitialiseSingleSAMBufferQueue(single_sam_buffer_queue *queue)
{
    InitialiseMutex(queue->rwMutex);
    queue->queueLength = 0;
}

global_function
void
AddSingleSAMBufferToQueue(single_sam_buffer_queue *queue, sam_buffer *buffer);

global_function
void
InitialiseSAMBufferQueue(memory_arena *arena, sam_buffer_queue *queue, process_queue *processQueue, process_queue *writeQueue)
{
    queue->queues = PushArrayP(arena, single_sam_buffer_queue *, Number_Of_SAM_Buffer_Queues);
    queue->index = 0;

    ForLoop(Number_Of_SAM_Buffer_Queues)
    {
        queue->queues[index] = PushStructP(arena, single_sam_buffer_queue);
        InitialiseSingleSAMBufferQueue(queue->queues[index]);

        ForLoop2(Number_Of_SAM_Buffers_Per_Queue)
        {
            sam_buffer *buffer = PushStructP(arena, sam_buffer);
            buffer->homeIndex = index;
            buffer->processQueue = processQueue;
            buffer->writeQueue = writeQueue;
            AddSingleSAMBufferToQueue(queue->queues[index], buffer);
        }
    }
}

global_function
void
AddSingleSAMBufferToQueue(single_sam_buffer_queue *queue, sam_buffer *buffer)
{
    LockMutex(queue->rwMutex);
    buffer->prev = 0;

    switch (queue->queueLength)
    {
        case 0:
            queue->front = buffer;
            queue->rear = buffer;
            break;

        default:
            queue->rear->prev = buffer;
            queue->rear = buffer;
    }

    ++queue->queueLength;
    UnlockMutex(queue->rwMutex);
}

global_function
void
AddSAMBufferToQueue(sam_buffer_queue *queue, sam_buffer *buffer)
{
    single_sam_buffer_queue *singleQueue = queue->queues[buffer->homeIndex];
    AddSingleSAMBufferToQueue(singleQueue, buffer);
}

global_function
sam_buffer *
TakeSingleSAMBufferFromQueue(single_sam_buffer_queue *queue)
{
    LockMutex(queue->rwMutex);
    sam_buffer *buffer = queue->front;

    switch (queue->queueLength)
    {
        case 0:
            break;

        case 1:
            queue->front = 0;
            queue->rear = 0;
            queue->queueLength = 0;
            break;

        default:
            queue->front = buffer->prev;
            --queue->queueLength;
    }

    UnlockMutex(queue->rwMutex);

    return(buffer);
}

global_function
single_sam_buffer_queue *
GetSingleSAMBufferQueue(sam_buffer_queue *queue)
{
    u32 index = __atomic_fetch_add(&queue->index, 1, 0) % Number_Of_SAM_Buffer_Queues;
    return(queue->queues[index]);
}

global_function
sam_buffer *
TakeSAMBufferFromQueue(sam_buffer_queue *queue)
{
    return(TakeSingleSAMBufferFromQueue(GetSingleSAMBufferQueue(queue)));
}

global_function
sam_buffer *
TakeSAMBufferFromQueue_Wait(sam_buffer_queue *queue)
{
    sam_buffer *buffer = 0;
    while (!buffer)
    {
        buffer = TakeSAMBufferFromQueue(queue);
    }
    return(buffer);
}
