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

#define Number_Of_Restriction_Site_Regexs_Per_Queue 32
#define Number_Of_Restriction_Site_Regex_Queues 32

struct
restriction_site_regex
{
    u32 pad;
    u32 homeIndex;
    hs_database_t *database;
    hs_scratch_t *scratch;
    memory_arena *arena;
    restriction_site_regex *prev;
};

struct
single_restriction_site_regex_queue
{
    u32 queueLength;
    u32 pad;
    mutex rwMutex;
    restriction_site_regex *front;
    restriction_site_regex *rear;
};

struct
restriction_site_regex_queue
{
    single_restriction_site_regex_queue **queues;
    threadSig index;
    u32 pad;
};

global_function
void
InitialiseSingleRestrictionSiteRegexQueue(single_restriction_site_regex_queue *queue)
{
    InitialiseMutex(queue->rwMutex);
    queue->queueLength = 0;
}

global_function
void
AddSingleRestrictionSiteRegexToQueue(single_restriction_site_regex_queue *queue, restriction_site_regex *regex);

global_function
u32
InitialiseRegex(restriction_site_regex *regex, char *pattern);

global_function
u32
InitialiseRestrictionSiteRegexQueue(memory_arena *arena, restriction_site_regex_queue *queue, char *pattern, u32 individualMemCapacity = MegaByte(1))
{
    u32 result = 1;
    
    queue->queues = PushArrayP(arena, single_restriction_site_regex_queue *, Number_Of_Restriction_Site_Regex_Queues);
    queue->index = 0;

    ForLoop(Number_Of_Restriction_Site_Regex_Queues)
    {
        queue->queues[index] = PushStructP(arena, single_restriction_site_regex_queue);
        InitialiseSingleRestrictionSiteRegexQueue(queue->queues[index]);

        ForLoop2(Number_Of_Restriction_Site_Regexs_Per_Queue)
        {
            restriction_site_regex *regex = PushStructP(arena, restriction_site_regex);
            regex->arena = 0;
            if ((result = InitialiseRegex(regex, pattern)))
            {
                regex->homeIndex = index;
                AddSingleRestrictionSiteRegexToQueue(queue->queues[index], regex);
                regex->arena = PushStructP(arena, memory_arena);
                CreateMemoryArenaP(regex->arena, individualMemCapacity);
            }
            else break;
        }

        if (!result) break;
    }

    return(result);
}

global_function
u32
InitialiseRegex(restriction_site_regex *regex, char *pattern)
{
    u32 result = 1;

    hs_compile_error_t *compile_err;
    if (hs_compile(pattern, HS_FLAG_SOM_LEFTMOST | HS_FLAG_CASELESS, HS_MODE_BLOCK, NULL, &regex->database, &compile_err) != HS_SUCCESS)
    {
        PrintError("Unable to compile pattern \'%s\': %s", pattern, compile_err->message);
        hs_free_compile_error(compile_err);
        result = 0;
    }

    if (result && hs_alloc_scratch(regex->database, &regex->scratch) != HS_SUCCESS)
    {
        PrintError("Unable to allocate regex scratch space for pattern \'%s\'", pattern);
        hs_free_database(regex->database);
        result = 0;
    }

    return(result);
}

global_function
void
AddSingleRestrictionSiteRegexToQueue(single_restriction_site_regex_queue *queue, restriction_site_regex *regex)
{
    LockMutex(queue->rwMutex);
    regex->prev = 0;

    switch (queue->queueLength)
    {
        case 0:
            queue->front = regex;
            queue->rear = regex;
            break;

        default:
            queue->rear->prev = regex;
            queue->rear = regex;
    }

    ++queue->queueLength;
    UnlockMutex(queue->rwMutex);
}

global_function
void
AddRestrictionSiteRegexToQueue(restriction_site_regex_queue *queue, restriction_site_regex *regex)
{
    single_restriction_site_regex_queue *singleQueue = queue->queues[regex->homeIndex];
    AddSingleRestrictionSiteRegexToQueue(singleQueue, regex);
}

global_function
restriction_site_regex *
TakeSingleRestrictionSiteRegexFromQueue(single_restriction_site_regex_queue *queue)
{
    LockMutex(queue->rwMutex);
    restriction_site_regex *regex = queue->front;

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
            queue->front = regex->prev;
            --queue->queueLength;
    }

    UnlockMutex(queue->rwMutex);

    return(regex);
}

global_function
single_restriction_site_regex_queue *
GetSingleRestrictionSiteRegexQueue(restriction_site_regex_queue *queue)
{
    u32 index = __atomic_fetch_add(&queue->index, 1, 0) % Number_Of_Restriction_Site_Regex_Queues;
    return(queue->queues[index]);
}

global_function
restriction_site_regex *
TakeRestrictionSiteRegexFromQueue(restriction_site_regex_queue *queue)
{
    return(TakeSingleRestrictionSiteRegexFromQueue(GetSingleRestrictionSiteRegexQueue(queue)));
}

global_function
restriction_site_regex *
TakeRestrictionSiteRegexFromQueue_Wait(restriction_site_regex_queue *queue)
{
    restriction_site_regex *regex = 0;
    while (!regex)
    {
        regex = TakeRestrictionSiteRegexFromQueue(queue);
    }
    return(regex);
}

global_function
void
FreeSingleRestrictionSiteRegexQueue(single_restriction_site_regex_queue *queue)
{
    LockMutex(queue->rwMutex);

    for (   restriction_site_regex *regex = queue->front;
            regex;
            regex = regex->prev )
    {
        if (regex->scratch) hs_free_scratch(regex->scratch);
        if (regex->database) hs_free_database(regex->database);
        if (regex->arena) FreeMemoryArenaP(regex->arena);
    }

    UnlockMutex(queue->rwMutex);
}

global_function
void
FreeRestrictionSiteRegexQueue(restriction_site_regex_queue *queue)
{
    ForLoop(Number_Of_Restriction_Site_Regex_Queues)
    {
        FreeSingleRestrictionSiteRegexQueue(queue->queues[index]);
    }
}
