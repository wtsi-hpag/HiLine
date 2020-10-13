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

#define PY_SSIZE_T_CLEAN
#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>
#include <Python.h>

#include "Header.h"

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wreserved-id-macro"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wextra-semi-stmt"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconditional-uninitialized"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"
#pragma clang diagnostic pop

#define String_(x) #x
#define String(x) String_(x)
#define ErrorLine String(Error at line __LINE__)

#define TestForPythonError \
do \
{ \
    if ((PyErr_Occurred())) \
    { \
        PyErr_SetString(PyExc_Exception, ErrorLine); \
        return(0); \
    } \
} while (0)

global_variable
u08
Status_Marco_Expression_Sponge = 0;

global_variable
char
Message_Buffer[1024];

global_variable
s32
Error_Handle;

#define PrintError(message, ...) \
{ \
    u32 halfSize = sizeof(Message_Buffer) >> 1; \
    stbsp_snprintf(Message_Buffer, halfSize, message, ##__VA_ARGS__); \
    u64 n = (u64)stbsp_snprintf(Message_Buffer + halfSize, halfSize, "%s\n", Message_Buffer); \
    write(Error_Handle, Message_Buffer + halfSize, n); \
    PyErr_SetString(PyExc_Exception, Message_Buffer + halfSize); \
} \
Status_Marco_Expression_Sponge = 0

global_variable
s32
Status_Handle;

#define PrintStatus(message, ...) \
{ \
    u32 halfSize = sizeof(Message_Buffer) >> 1; \
    stbsp_snprintf(Message_Buffer, halfSize, message, ##__VA_ARGS__); \
    u64 n = (u64)stbsp_snprintf(Message_Buffer + halfSize, halfSize, "%s\n", Message_Buffer); \
    write(Status_Handle, Message_Buffer + halfSize, n); \
} \
Status_Marco_Expression_Sponge = 0

global_variable
memory_arena
Working_Set;

struct
sequence
{
    char *name;
    u32 length;
    u08 frontOverhang;
    u08 endOverhang;
    u08 frontOverhangIsFwd;
    u08 endOverhangIsFwd;
    sequence *next;
};

struct
sequence_hash_table_node
{
    sequence *sequence;
    sequence_hash_table_node *next;
};

struct
sequence_hash_table
{
    sequence_hash_table_node **table;
    u32 size;
    u32 pad;
};

global_function
sequence_hash_table *
CreateSequenceHashTable(memory_arena *arena, u32 nEntries)
{
    u32 size = NextPrime((u32)((f32)nEntries * 1.3f));

    sequence_hash_table *table = PushStructP(arena, sequence_hash_table);
    table->size = size;
    table->table = PushArrayP(arena, sequence_hash_table_node *, size);
    memset(table->table, 0, size * sizeof(sequence_hash_table_node *));

    return(table);
}

global_function
void
AddSequenceToHashTable(memory_arena *arena, sequence_hash_table *table, sequence *sequence)
{
#define SequenceHashTableSeed 0xf20b503a1896576e
    u32 name32[16];
    PushStringIntoIntArray((u32 *)name32, ArrayCount(name32), (u08 *)sequence->name);

    u32 hash = FastHash32(name32, sizeof(name32), SequenceHashTableSeed) % table->size;
    sequence_hash_table_node *node = table->table[hash];
    sequence_hash_table_node *prevNode = 0;

    while (node)
    {
        prevNode = node;
        node = node->next;
    }

    sequence_hash_table_node *newNode = PushStructP(arena, sequence_hash_table_node);
    newNode->sequence = sequence;
    newNode->next = 0;

    (prevNode ? prevNode->next : table->table[hash]) = newNode;
}

global_function
sequence *
GetSequenceFromHashTable(sequence_hash_table *table, char *name, char nameTerm = '\t')
{
    u32 name32[16];
    PushStringIntoIntArray((u32 *)name32, ArrayCount(name32), (u08 *)name, (u08)nameTerm);

    u32 hash = FastHash32(name32, sizeof(name32), SequenceHashTableSeed) % table->size;
    sequence_hash_table_node *node = table->table[hash];
    sequence_hash_table_node *prevNode = 0;
    sequence *result = 0;

    while (node)
    {
        prevNode = node;
        sequence *test = node->sequence;
        if (AreStringsEqual(name, nameTerm, test->name, '\0'))
        {
            result = test;
            break;
        }
        node = node->next;
    }

    return(result);
}

struct
write_buffer
{
    u08 *buffer;
    u64 size;
};

global_function
write_buffer *
CreateWriteBuffer(memory_arena *arena, u64 bufferSize)
{
    write_buffer *buffer = PushStructP(arena, write_buffer);
    buffer->buffer = PushArrayP(arena, u08, bufferSize);
    buffer->size = 0;
    
    return(buffer);
}

global_variable
s32
Sam_Out;

global_variable
u08
Global_Write_Error_Flag = 0;

global_function
void
WriteFunction(void *in)
{
    write_buffer *buffer = (write_buffer *)in;
    if (write(Sam_Out, (const void *)buffer->buffer, buffer->size) != buffer->size)
    {
        PrintError("Error writing to output");
        Global_Write_Error_Flag = 1;
    }
}

struct
read_buffer
{
    u08 *buffer;
    u64 size;
};

struct
read_pool
{
    thread_pool *pool;
    s32 handle;
    u32 bufferPtr;
    read_buffer *buffers[2];
};

global_function
read_pool *
CreateReadPool(memory_arena *arena)
{
    read_pool *pool = PushStructP(arena, read_pool);
    pool->pool = ThreadPoolInit(arena, 1);

#define ReadBufferSize MegaByte(16)
    pool->bufferPtr = 0;
    pool->buffers[0] = PushStructP(arena, read_buffer);
    pool->buffers[0]->buffer = PushArrayP(arena, u08, ReadBufferSize);
    pool->buffers[0]->size = 0;
    pool->buffers[1] = PushStructP(arena, read_buffer);
    pool->buffers[1]->buffer = PushArrayP(arena, u08, ReadBufferSize);
    pool->buffers[1]->size = 0;

    return(pool);
}

global_function
void
FillReadBuffer(void *in)
{
    read_pool *pool = (read_pool *)in;
    read_buffer *buffer = pool->buffers[pool->bufferPtr];
    buffer->size = (u64)read(pool->handle, buffer->buffer, ReadBufferSize);
}

global_function
read_buffer *
GetNextReadBuffer(read_pool *readPool)
{
    FenceIn(ThreadPoolWait(readPool->pool));
    read_buffer *buffer = readPool->buffers[readPool->bufferPtr];
    readPool->bufferPtr = (readPool->bufferPtr + 1) & 1;
    ThreadPoolAddTask(readPool->pool, FillReadBuffer, readPool);
    return(buffer);
}

global_function
void
CloseReadPool(read_pool *pool)
{
    ThreadPoolWait(pool->pool);
    ThreadPoolDestroy(pool->pool);
}

global_function
PyObject *
Aligner_Main(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *objParams;
    static const char *kwlist[] = {"params", 0};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char **)kwlist, &objParams))
    {
        return(0);
    }

    PyObject *objSamInput, *objSamOutput, *objNThreads, *objInfo, *objError;

    {
        struct
            inputParam
            {
                PyObject **obj;
                char *name;
            };

        inputParam params[] = {
            {
                &objSamInput,
                (char *)"samInput"
            },
            {
                &objSamOutput,
                (char *)"samOutput"
            },
            {
                &objNThreads,
                (char *)"nThreads"
            },
            {
                &objInfo,
                (char *)"info"
            },
            {
                &objError,
                (char *)"error"
            }
        };

        ForLoop(ArrayCount(params))
        {
            inputParam *param = params + index;
            PyObject *objName = Py_BuildValue("s", param->name);
            *(param->obj) = PyObject_GetAttr(objParams, objName);
            Py_DECREF(objName);
            if (!(*(param->obj)))
            {
                PrintError("No param attribute \'%s\' found", param->name);
                return(0);
            }
        }
    }

    s32 exitCode = EXIT_SUCCESS;

    if ((Status_Handle = (s32)PyObject_AsFileDescriptor(objInfo)) < 0)
    {
        PyErr_SetString(PyExc_Exception, "Cannot open the info/status handle");
        return(0);
    }

    if ((Error_Handle = (s32)PyObject_AsFileDescriptor(objError)) < 0)
    {
        PyErr_SetString(PyExc_Exception, "Cannot open the error handle");
        return(0);
    }

    u32 nThreads = (u32)PyLong_AsUnsignedLong(objNThreads);
    TestForPythonError;
    Py_DECREF(objNThreads);

    CreateMemoryArena(Working_Set, MegaByte(512));
    thread_pool *writePool = ThreadPoolInit(&Working_Set, 1);
    read_pool *readPool = CreateReadPool(&Working_Set);

    s32 samIn;
    if ((samIn = PyObject_AsFileDescriptor(objSamInput)) > 0)
    {
        readPool->handle = samIn;
        if ((Sam_Out = PyObject_AsFileDescriptor(objSamOutput)) > 0)
        {
            Py_BEGIN_ALLOW_THREADS;
            u08 samLine[KiloByte(16)];
            u32 linePtr = 0;
            u08 goodRead = 1;

            write_buffer *writeBuffers[2];
            writeBuffers[0] = CreateWriteBuffer(&Working_Set, ReadBufferSize);
            writeBuffers[1] = CreateWriteBuffer(&Working_Set, ReadBufferSize);
            u08 writeBufferFlip = 0;
            
            write_buffer *currentWriteBuffer = writeBuffers[writeBufferFlip];

            u08 headerMode = 1;
            u64 numHeaderLines = 0;
            u32 numSequences = 0;
            sequence *sequences = 0;
            sequence *tailSequence = 0;
            sequence_hash_table *sequenceHashTable;

            u64 totalRead = 0;
            char printNBuffers[2][16] = {{0}};
            u08 printNBufferPtr = 0;

            read_buffer *readBuffer = GetNextReadBuffer(readPool);
            do
            {
                readBuffer = GetNextReadBuffer(readPool);
                
                for (   u64 bufferIndex = 0;
                        bufferIndex < readBuffer->size;
                        ++bufferIndex )
                {
                    u08 character = readBuffer->buffer[bufferIndex];
                    samLine[linePtr++] = character;

                    if (character == '\n')
                    {
                        if (headerMode && samLine[0] != '@')
                        {
                            sequenceHashTable = CreateSequenceHashTable(&Working_Set, numSequences);
                            PrintStatus("Created digestion fragment hash-table with capacity %u", sequenceHashTable->size);
                            
                            TraverseLinkedList(sequences, sequence) AddSequenceToHashTable(&Working_Set, sequenceHashTable, node);
                            headerMode = 0;

                            PrintStatus("Inserted %u sequences into hash-table", numSequences);
                        }

                        if (headerMode && samLine[1] == 'S' && samLine[2] == 'Q')
                        {
                            ++numHeaderLines;
                            u32 nameStart = 3;
                            while (!(samLine[nameStart] == '\t' && samLine[nameStart + 1] == 'S' && samLine[nameStart + 2] == 'N' && samLine[nameStart + 3] == ':')) ++nameStart;
                            nameStart += 4;

                            u08 *name = samLine + nameStart;
                            u32 nameLength = 1;
                            while (*name != '\t' && *name != '\n')
                            {
                                ++name;
                                ++nameLength;
                            }

                            name = PushArray(Working_Set, u08, nameLength);
                            ForLoop((nameLength - 1))
                            {
                                name[index] = samLine[nameStart + index];
                            }
                            name[nameLength - 1] = 0;

                            u08 *nameTmp = name + nameLength - 1;
                            u32 endLength = 0;
                            while (*(--nameTmp) != '_') ++endLength;
                            u32 end = StringToInt(name + nameLength - 2, endLength - 1);
                            u08 endOverhangIsFwd = name[nameLength - 2] == '+' ? 1 : 0;

                            u08 *nameTmpTmp = nameTmp;
                            u32 startLength = 0;
                            while (*(--nameTmpTmp) != '_') ++startLength;
                            u32 start = StringToInt(nameTmp - 1, startLength - 1);
                            u08 startOverhangIsFwd = nameTmp[-1] == '+' ? 1 : 0;

                            u32 lengthStart = 3;
                            while (!(samLine[lengthStart] == '\t' && samLine[lengthStart + 1] == 'L' && samLine[lengthStart + 2] == 'N' && samLine[lengthStart + 3] == ':')) ++lengthStart;
                            lengthStart += 4;

                            u08 *length = samLine + lengthStart;
                            u32 lengthLength = 0;
                            while (*length != '\t' && *length != '\n')
                            {
                                ++length;
                                ++lengthLength;
                            }
                            u32 seqLength = StringToInt(length, lengthLength);

                            sequence *newSequence = PushStruct(Working_Set, sequence);
                            newSequence->name = (char *)name;
                            newSequence->next = 0;
                            newSequence->length = seqLength;
                            newSequence->frontOverhang = (u08)start;
                            newSequence->endOverhang = (u08)end;
                            newSequence->frontOverhangIsFwd = startOverhangIsFwd;
                            newSequence->endOverhangIsFwd = endOverhangIsFwd;
                            ++numSequences;

                            if (!sequences) sequences = newSequence;
                            else tailSequence->next = newSequence;
                            tailSequence = newSequence;
                        }

                        if (!headerMode)
                        {
                            u32 copyPtr = 0;

                            while (samLine[copyPtr++] != '\t') {} // read name
                            u32 flagsLength = 0;
                            while (samLine[copyPtr++] != '\t') ++flagsLength; // flags
                            u08 readIsRev = (StringToInt(samLine + copyPtr - 1, flagsLength) & 16) ? 1 : 0;

                            sequence *seq;
                            if ((seq = GetSequenceFromHashTable(sequenceHashTable, (char *)(samLine + copyPtr))))
                            {
                                while (samLine[copyPtr++] != '\t') {} // ref name

                                u32 posLen = 0;
                                while (samLine[copyPtr++] != '\t') ++posLen; // pos
                                u32 pos = StringToInt(samLine + copyPtr - 1, posLen);

                                if (pos)
                                {
                                    pos -= 1;

                                    while (samLine[copyPtr++] != '\t') {} // mapq
                                    u32 startSamLength = copyPtr;

                                    if (samLine[startSamLength] != '*')
                                    {
                                        u32 cigarLen = 0;
                                        while (samLine[copyPtr++] != '\t') ++cigarLen; // cigar
                                        u32 otherStart = copyPtr;

                                        u32 otherLength = 2;
                                        while (samLine[copyPtr++] != '\t') ++otherLength; // next read name
                                        while (samLine[copyPtr++] != '\t') ++otherLength; // next pos
                                        while (samLine[copyPtr++] != '\t') ++otherLength; // temp length

                                        u32 seqStart = copyPtr;
                                        u32 seqLength = 0;
                                        while (samLine[copyPtr++] != '\t') ++seqLength; // seq

                                        u32 qualStart = copyPtr;
                                        while (samLine[copyPtr++] != '\t') {} // qual

                                        u32 firstRefCigar = 0;
                                        {
                                            u08 *cigar = samLine + startSamLength;
                                            u32 numLen = 0;
                                            u08 cigOp;
                                            while ((cigOp = *(cigar++)) != '\t')
                                            {
                                                u08 process = 1;
                                                u08 ref = 0;
                                                switch (cigOp)
                                                {
                                                    case 'M':
                                                    case 'D':
                                                    case 'N':
                                                    case '=':
                                                    case 'X':
                                                        ref = 1;
                                                    case 'I':
                                                    case 'S':
                                                    case 'H':
                                                    case 'P':
                                                        break;

                                                    default:
                                                        ++numLen;
                                                        process = 0;
                                                }

                                                if (ref) break;

                                                if (process)
                                                {
                                                    firstRefCigar += (numLen ? StringToInt(cigar - 1, numLen) : 1);
                                                    numLen = 0;
                                                }
                                            }
                                        }

                                        u32 frontOverhang = readIsRev == seq->frontOverhangIsFwd ? 0 : seq->frontOverhang;
                                        u32 endOverhang = readIsRev == seq->endOverhangIsFwd ? 0 : seq->endOverhang;

                                        u32 startOffset = firstRefCigar + frontOverhang;
                                        startOffset = startOffset > pos ? startOffset - pos : 0;
                                        pos = startOffset ? 0 : pos - frontOverhang;

                                        u32 referenceLength = seq->length;
                                        referenceLength = referenceLength > frontOverhang ? referenceLength - frontOverhang : 0;
                                        referenceLength = referenceLength > endOverhang ? referenceLength - endOverhang : 0;
                                        referenceLength = referenceLength > pos ? referenceLength - pos : 0;

                                        u08 newSamLine[sizeof(samLine)];
                                        u32 newLinePtr = (u32)stbsp_snprintf((char *)newSamLine, startSamLength + 1, "%s", samLine);

                                        if (referenceLength)
                                        {
                                            u32 seqStartOffset = 0;
                                            u32 seqEndOffset = 0;

                                            {
                                                u08 *cigar = samLine + startSamLength;
                                                u32 numLen = 0;
                                                u08 cigOp;
                                                while ((cigOp = *(cigar++)) != '\t')
                                                {
                                                    u08 process = 1;
                                                    u08 ref = 0;
                                                    u08 qu = 0;
                                                    switch (cigOp)
                                                    {
                                                        case 'I':
                                                        case 'S':
                                                        case 'M':
                                                        case '=':
                                                        case 'X':
                                                            qu = 1;
                                                            if (cigOp == 'I' || cigOp == 'S') break;
                                                        case 'D':
                                                        case 'N':
                                                            ref = 1;
                                                        case 'H':
                                                        case 'P':
                                                            break;

                                                        default:
                                                            ++numLen;
                                                            process = 0;
                                                    }

                                                    if (process)
                                                    {
                                                        u32 n = numLen ? StringToInt(cigar - 1, numLen) : 1;
                                                        numLen = 0;

                                                        if (startOffset)
                                                        {
                                                            if (startOffset > n)
                                                            {
                                                                startOffset -= n;
                                                                if (qu) seqStartOffset += n;
                                                                n = 0;
                                                            }
                                                            else
                                                            {
                                                                if (qu) seqStartOffset += startOffset;
                                                                n -= startOffset;
                                                                startOffset = 0;
                                                            }
                                                        }

                                                        if (!startOffset)
                                                        {
                                                            if (referenceLength > n) referenceLength -= n;
                                                            else
                                                            {
                                                                if (qu) seqEndOffset += (n - referenceLength);
                                                                n = referenceLength;
                                                                referenceLength = 0;
                                                            }

                                                            if (n) newLinePtr += (u32)stbsp_snprintf((char *)newSamLine + newLinePtr, 6, "%u%c", n, cigOp);
                                                        }
                                                    }
                                                }
                                            }

                                            newLinePtr += (u32)stbsp_snprintf((char *)newSamLine + newLinePtr, otherLength + 2, "\t%s", samLine + otherStart);

                                            u32 totalOffset = seqStartOffset + seqEndOffset;
                                            if (totalOffset < seqLength)
                                            {
                                                newLinePtr += (u32)stbsp_snprintf((char *)newSamLine + newLinePtr, seqLength - totalOffset + 2, "\t%s", samLine + seqStart + seqStartOffset);
                                                newLinePtr += (u32)stbsp_snprintf((char *)newSamLine + newLinePtr, seqLength - totalOffset + 2, "\t%s", samLine + qualStart + seqStartOffset);
                                            }
                                            else
                                            {
                                                //newLinePtr += (u32)stbsp_snprintf((char *)newSamLine + newLinePtr, 5, "\t*\t*");
                                                goodRead = 0;
                                            }
                                        }
                                        else
                                        {
                                            //newLinePtr += (u32)stbsp_snprintf((char *)newSamLine + newLinePtr, 5, "\t*\t*");
                                            goodRead = 0;
                                        }

                                        if (goodRead)
                                        {
                                            newLinePtr += (u32)stbsp_snprintf((char *)newSamLine + newLinePtr, linePtr - copyPtr + 2, "\t%s", samLine + copyPtr);

                                            stbsp_snprintf((char *)samLine, newLinePtr + 1, "%s", newSamLine);
                                            linePtr = newLinePtr;
                                        }
                                    }
                                }
                            }
                        }

                        if (goodRead)
                        {
                            if ((u64)linePtr > (ReadBufferSize - currentWriteBuffer->size))
                            {
                                if (Global_Write_Error_Flag)
                                {
                                    exitCode = EXIT_FAILURE;
                                    goto End;
                                }

                                FenceIn(ThreadPoolWait(writePool));

                                ThreadPoolAddTask(writePool, WriteFunction, currentWriteBuffer);

                                writeBufferFlip = (writeBufferFlip + 1) & 1;
                                currentWriteBuffer = writeBuffers[writeBufferFlip];
                                currentWriteBuffer->size = 0;
                            }

                            ForLoop(linePtr) currentWriteBuffer->buffer[currentWriteBuffer->size++] = samLine[index];
                        }
                        linePtr = 0;
                        goodRead = 1;

#define Log2_Print_Interval 14
                        if (!(++totalRead & ((1 << Log2_Print_Interval) - 1)))
                        {
                            u08 currPtr = printNBufferPtr;
                            u08 otherPtr = (currPtr + 1) & 1;
                            stbsp_snprintf(printNBuffers[currPtr], sizeof(printNBuffers[currPtr]), "%$" PRIu64, headerMode ? numHeaderLines : totalRead - numHeaderLines);

                            if (strcmp(printNBuffers[currPtr], printNBuffers[otherPtr]))
                            {
                                PrintStatus("%s %s processed", printNBuffers[currPtr], headerMode ? "header lines" : "reads");
                            }

                            printNBufferPtr = otherPtr;
                        }
                    }
                }
            } while (readBuffer->size);

            ThreadPoolAddTask(writePool, WriteFunction, currentWriteBuffer);
            ThreadPoolWait(writePool);
            
            Py_END_ALLOW_THREADS;
        }
        else
        {	
            PrintError("Could not open sam output");
            exitCode = EXIT_FAILURE;
        }
    }
    else
    {	
        PrintError("Could not open sam input");
        exitCode = EXIT_FAILURE;
    }

End:
    CloseReadPool(readPool);
    ThreadPoolWait(writePool);
    ThreadPoolDestroy(writePool);
    FreeMemoryArena(Working_Set);
    if (exitCode == EXIT_FAILURE || Global_Write_Error_Flag) return(0);
    Py_RETURN_NONE;
}

#define ModuleName "_Aligner_Main"
#define ModuleDescription "HiLine alignment read trimmer"

global_variable
PyMethodDef
AlignerMethods[] =
{
    {ModuleName, (PyCFunction) Aligner_Main, METH_VARARGS | METH_KEYWORDS, ModuleDescription},
    {NULL, NULL, 0, NULL}
};

global_variable
struct
PyModuleDef
AlignerModule =
{
    PyModuleDef_HEAD_INIT,
    ModuleName,
    ModuleDescription,
    -1,
    AlignerMethods
};

PyMODINIT_FUNC
PyInit__Aligner()
{
    return(PyModule_Create(&AlignerModule));
}
