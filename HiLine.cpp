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
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>
#include <Python.h>
#include "numpy/arrayobject.h"

#include "Header.h"
#include "WAVLTree.cpp"
#include "LineBufferQueue.cpp"
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#include "hs.h"
#pragma clang diagnostic pop

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

global_variable
char *
Name;

global_variable
char *
Version;

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

#include "RestrictionSiteRegexQueue.cpp"

global_variable
memory_arena
Working_Set;

global_variable
thread_pool *
Thread_Pool;

struct
restriction_site
{
	u32 pad;
	u32 restrictionLocation;
	char *pattern;
	restriction_site_regex_queue queue;
};

struct
restriction_sites
{
	u32 longestLen;
	u32 num;
	restriction_site *sites;
};

global_variable
restriction_sites *
Restriction_Sites = 0;

global_function
u32
CompileRestrictionSiteRegex(restriction_site *site)
{
	return(InitialiseRestrictionSiteRegexQueue(&Working_Set, &site->queue, site->pattern));
}

global_function
void
FreeRestrictionSites(restriction_sites *sites)
{
	if (sites)
	{
		ForLoop(sites->num)
		{
			restriction_site *site = sites->sites + index;
			FreeRestrictionSiteRegexQueue(&site->queue);
		}
	}
}

struct
scan_context_data
{
	u64 offset;
	memory_arena *arena;
	char *name;
	line_buffer_value **headValue;
	line_buffer_value **tailValue;
};

global_function
s32
EventHandler(u32 id, unsigned long long from, unsigned long long to, u32 flags, void *ctx)
{
	(void) id;
	(void) to;
	(void) flags;

	scan_context_data *data = (scan_context_data *)ctx;

	line_buffer_value *value = PushStructP(data->arena, line_buffer_value);
	value->name = data->name;
	value->value = data->offset + (u64)from;
	value->next = 0;

	if (!(*(data->headValue))) *(data->headValue) = value;
	else (*(data->tailValue))->next = value;

	*(data->tailValue) = value;

	return(0);
}

global_variable
threadSig
Global_Regex_Scan_Error_Flag = 0;

global_function
void
ScanForRestrictionSite(void *in)
{
	line_buffer *buffer = (line_buffer *)in;
	char *seq = (char *)buffer->line;
	u32 lineLen = buffer->lineLength;

	scan_context_data data;
	data.name = buffer->name;
	data.headValue = &buffer->headValue;
	data.tailValue = &buffer->tailValue;

	ForLoop(Restriction_Sites->num)
	{
		restriction_site *site = Restriction_Sites->sites + index;

		restriction_site_regex *regex = TakeRestrictionSiteRegexFromQueue_Wait(&site->queue);
		u64 offset = (u64)site->restrictionLocation + buffer->startCoord;

		data.arena = regex->arena;
		data.offset = offset;

		if (hs_scan(regex->database, seq, lineLen, 0, regex->scratch, EventHandler, &data) != HS_SUCCESS)
		{
			PrintError("Unable to scan input buffer");
			Global_Regex_Scan_Error_Flag = 1;
		}

		AddRestrictionSiteRegexToQueue(&site->queue, regex);
	}

	AddLineBufferToQueue(Line_Buffer_Queue, buffer);
}

struct
sequence
{
	char *name;
	u64 length;
	wavl_tree *tree;
	sequence *next;
	u32 index;
	u32 pad;
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

global_variable
sequence_hash_table *
Sequence_Hash_Table;

struct
create_tree_data
{
	line_buffer_value **values;
	sequence *sequence;
	memory_arena *arena;
	u32 numValueLists;
	u32 pad;
};

global_function
void
CreateTree_Thread(void *in)
{
	create_tree_data *data = (create_tree_data *)in;
	ForLoop(data->numValueLists)
	{
		TraverseLinkedList(data->values[index], line_buffer_value)
		{
			WavlTreeInsertValue(data->arena, data->sequence->tree, node->value);
		}
	}
	WavlTreeFreeze(data->sequence->tree);
}

struct sam_buffer;

struct
process_queue
{
	mutex lock;
	binary_semaphore *hasItems;
	sam_buffer *head;
	sam_buffer *tail;
	void (*function)(sam_buffer *buffer);
	volatile u64 nextSAMIndex;
	volatile u64 isClosed;
};

#include "SAMLineBufferQueue.cpp"

global_function
sam_buffer *
PullSAMRecordFromProcessQueue(process_queue *queue)
{
	BinarySemaphoreWait(queue->hasItems);
	LockMutex(queue->lock);

	sam_buffer *buffer = queue->head;
	if (buffer) queue->head = buffer->prev;
	if (queue->head || !queue->tail) BinarySemaphorePost(queue->hasItems);
	else queue->tail = 0;

	UnlockMutex(queue->lock);

	return(buffer);
}

global_function
void
ProcessQueueJob(void *in)
{
	process_queue *queue = (process_queue *)in;
	sam_buffer *buffer;
	while ((buffer = PullSAMRecordFromProcessQueue(queue))) queue->function(buffer);
	queue->function(0);
	queue->isClosed = 1;
}

global_function
process_queue *
CreateProcessQueue(memory_arena *arena, void (*function)(sam_buffer*), thread_pool *pool)
{
	process_queue *queue = PushStructP(arena, process_queue);

	queue->hasItems = PushStructP(arena, binary_semaphore);
	BinarySemaphoreInit(queue->hasItems, 0);
	queue->head = 0;
	queue->tail = 0;
	InitialiseMutex(queue->lock);
	queue->function = function;
	queue->nextSAMIndex = 0;
	queue->isClosed = 0;

	ThreadPoolAddTask(pool, ProcessQueueJob, queue);

	return(queue);
}

global_function
void
AddSAMRecordToProcessQueue(process_queue *queue, sam_buffer *buffer, u32 ordered = 1)
{
	while (ordered && queue->nextSAMIndex != buffer->samIndex) {}

	LockMutex(queue->lock);

	if (!queue->head) queue->head = buffer;
	if (queue->tail) queue->tail->prev = buffer;
	queue->tail = buffer;
	if (buffer) buffer->prev = 0;

	if (ordered && buffer) queue->nextSAMIndex = buffer->samIndex + 1;

	BinarySemaphorePost(queue->hasItems);

	UnlockMutex(queue->lock);
}

global_function
void
CloseProcessQueue(process_queue *queue, u32 wait = 1)
{
	AddSAMRecordToProcessQueue(queue, 0, 0);

	if (wait)
	{
		while (!queue->isClosed) {}
	}
}

global_variable
threadSig
Global_SAM_Write_Error_Flag = 0;

global_variable
u32
Min_Map_Quality;

struct
output_handles
{
	union
	{
		struct
		{
			s32 unPaired;
			s32 unMapped;
			s32 secondary;
			s32 qcFail;
			s32 duplicate;
			s32 supplementary;
			s32 invalidReferenceName;
			s32 belowMinMapq;
			s32 tooFarFromRestrictionSite;
			s32 dump;
			s32 selfCircle;
			s32 danglingEnd;
			s32 sameFragmentAndStrand;
			s32 reLigated;
			s32 FF;
			s32 FR;
			s32 RF;
			s32 RR;
		};
		s32 handles[18];
	};
};

global_variable
output_handles *
Output_Handles;

global_variable
output_handles *
Unique_Output_Handles;

global_function
output_handles *
CreateUniqueOutputHandles(memory_arena *arena, output_handles *handles)
{
	output_handles *result = PushStructP(arena, output_handles);
	ForLoop(ArrayCount(result->handles))
	{
		result->handles[index] = -1;
	}

	auto IsIn = [result](s32 test)->u32
	{
		u32 isIn = 0;
		ForLoop(ArrayCount(result->handles))
		{
			if (test == result->handles[index])
			{
				isIn = 1;
				break;
			}
		}
		return(isIn);
	};

	u32 resultIndex = 0;
	ForLoop(ArrayCount(handles->handles))
	{
		if (!IsIn(handles->handles[index])) result->handles[resultIndex++] = handles->handles[index]; 
	}

	return(result);
}

global_variable
char *
Handle_Names[ArrayCount(Output_Handles->handles)] =
{
	"unpaired",
	"unmapped",
	"secondary",
	"qcfail",
	"duplicate",
	"supplementary",
	"invalidreferencename",
	"belowminmapq",
	"toofarfromrestrictionsite",
	"dump",
	"selfcircle",
	"danglingend",
	"samefragmentandstrand",
	"religated",
	"ff",
	"fr",
	"rf",
	"rr"
};

global_variable
char *
Handle_Short_Names[ArrayCount(Output_Handles->handles)] =
{
	"UP",
	"UM",
	"SE",
	"QC",
	"DP",
	"SP",
	"IR",
	"MQ",
	"RS",
	"DM",
	"SC",
	"DE",
	"SS",
	"RL",
	"FF",
	"FR",
	"RF",
	"RR"
};

struct
linked_list_node
{
	union
	{
		volatile u64 value;
		struct
		{
			volatile s32 signedValue;
			u32 pad;
		};
	};
	linked_list_node *next;
};

struct
linked_list
{
	linked_list_node *head;
	linked_list_node *tail;
};

#define DeclareLinkedListFunction(type) \
	global_function \
	void \
	PushOntoLinkedList(type value, linked_list *list, memory_arena *arena) \
{ \
	linked_list_node *node = PushStructP(arena, linked_list_node); \
	node->value = value; \
	node->next = 0; \
	if (!list->head) list->head = node; \
	if (list->tail) list->tail->next = node; \
	list->tail = node; \
}

DeclareLinkedListFunction(volatile u64)
DeclareLinkedListFunction(volatile s32)

global_function
void
CreatePGLineAndTags(sam_buffer *buffer, s32 handle)
{
	char nameBuff[128];
	u08 *namePtr = (u08 *)nameBuff;
	u32 index = buffer->lineLength - 1;

	buffer->sam[index++] = '\t';

	buffer->sam[index++] = 'm';
	buffer->sam[index++] = 'q';
	buffer->sam[index++] = ':';
	
	char numBuff[16];
	stbsp_snprintf(numBuff, sizeof(numBuff), "%u", Min_Map_Quality);

	namePtr = (u08 *)numBuff;
	while (*namePtr) buffer->sam[index++] = *(namePtr++);
	buffer->sam[index++] = '\t';

	ForLoop2((Restriction_Sites->num - 2))
	{
		buffer->sam[index++] = 'r';
		buffer->sam[index++] = 's';
		buffer->sam[index++] = ':';

		namePtr = (u08 *)Restriction_Sites->sites[index2].pattern;
		while (*namePtr) buffer->sam[index++] = *(namePtr++);

		stbsp_snprintf(numBuff, sizeof(numBuff), "%u", Restriction_Sites->sites[index2].restrictionLocation);

		buffer->sam[index++] = ',';
		namePtr = (u08 *)numBuff;
		while (*namePtr) buffer->sam[index++] = *(namePtr++);

		buffer->sam[index++] = '\t';
	}

	ForLoop2(ArrayCount(Output_Handles->handles))
	{
		if (Output_Handles->handles[index2] == handle)
		{
			buffer->sam[index++] = 'h';
			buffer->sam[index++] = 'c';
			buffer->sam[index++] = ':';

			namePtr = (u08 *)Handle_Short_Names[index2];
			while (*namePtr) buffer->sam[index++] = *(namePtr++);

			buffer->sam[index++] = '\t';
		}
	}

	--index;
	buffer->sam[index++] = '\n';

	buffer->lineLength = index;
}

global_function
void
AddTags(sam_buffer *buffer)
{
	u32 index = buffer->lineLength - 1;
	buffer->sam[index++] = '\t';

	buffer->sam[index++] = 'h';
	buffer->sam[index++] = 'c';
	buffer->sam[index++] = ':';
	buffer->sam[index++] = 'Z';
	buffer->sam[index++] = ':';

	u08 *namePtr = (u08 *)Handle_Short_Names[(s32)buffer->type];
	while (*namePtr) buffer->sam[index++] = *(namePtr++);

	if (buffer->fragment)
	{
		buffer->sam[index++] = '\t';
		ForLoop2(buffer->fragment->tagLength)
		{
			buffer->sam[index++] = buffer->fragment->tag[index2];
		}
	}

	buffer->sam[index++] = '\n';
	buffer->lineLength = index;
}

global_function
void
WriteSamRecordToHandle(sam_buffer *buffer, s32 handle, u08 addPG = 0)
{
	if (handle > 0)
	{
		u32 oldLength = 0;
		if (addPG)
		{
			oldLength = buffer->lineLength;
			CreatePGLineAndTags(buffer, handle);
		}

		if (buffer->sam[0] != '@') AddTags(buffer);

		Global_SAM_Write_Error_Flag = write(handle, buffer->sam, buffer->lineLength) != (ssize_t)buffer->lineLength;

		if (oldLength) buffer->lineLength = oldLength;

		if (Global_SAM_Write_Error_Flag)
		{
			buffer->sam[64] = 0;
			PrintError("Failed to write SAM record:\n%s", buffer->sam);
		}
	}
}

global_function
void
WriteSamRecord(sam_buffer *buffer)
{
	if (buffer)
	{
		if (buffer->sam[0] == '@')
		{
			u08 addPG = 0;
			if (buffer->sam[1] == 'P' && buffer->sam[2] == 'G')
			{
				u32 index = 4;
				while ((!(buffer->sam[index - 4] == '\t' && buffer->sam[index - 3] == 'P' && buffer->sam[index - 2] == 'N' && buffer->sam[index - 1] == ':')) && index < buffer->lineLength) ++index;
				
				if (index < buffer->lineLength && AreStringsEqual(Name, 0, (char *)(buffer->sam + index), '_'))
				{
					index += StringLength((u08 *)Name);

					while (buffer->sam[index + 1] != '\n')
					{
						buffer->sam[index] = buffer->sam[index + 1];
						++index;
					}
					
					buffer->sam[index++] = '\n';
					buffer->lineLength = index;

					addPG = 1;
				}
			}

			ForLoop(ArrayCount(Unique_Output_Handles->handles))
			{
				WriteSamRecordToHandle(buffer, Unique_Output_Handles->handles[index], addPG);
			}
		}
		else
		{
			s32 handle;
			switch (buffer->type)
			{
				case typeUnpaired:
					handle = Output_Handles->unPaired;
					break;

				case typeUnmapped:
					handle = Output_Handles->unMapped;
					break;

				case typeSecondary:
					handle = Output_Handles->secondary;
					break;

				case typeQCFail:
					handle = Output_Handles->qcFail;
					break;

				case typeDuplicate:
					handle = Output_Handles->duplicate;
					break;

				case typeSupplementary:
					handle = Output_Handles->supplementary;
					break;

				case typeInvalidReference:
					handle = Output_Handles->invalidReferenceName;
					break;

				case typeBelowMinMapQ:
					handle = Output_Handles->belowMinMapq;
					break;

				case typeTooFarFromRestrictionSite:
					handle = Output_Handles->tooFarFromRestrictionSite;
					break;

				case typeDump:
					handle = Output_Handles->dump;
					break;

				case typeSelfCircle:
					handle = Output_Handles->selfCircle;
					break;

				case typeDanglingEnd:
					handle = Output_Handles->danglingEnd;
					break;

				case typeSameFragmentAndStrand:
					handle = Output_Handles->sameFragmentAndStrand;
					break;

				case typeReLigated:
					handle = Output_Handles->reLigated;
					break;

				case typeFF:
					handle = Output_Handles->FF;
					break;

				case typeFR:
					handle = Output_Handles->FR;
					break;

				case typeRF:
					handle = Output_Handles->RF;
					break;

				case typeRR:
					handle = Output_Handles->RR;
					break;
			}

			WriteSamRecordToHandle(buffer, handle);
		}

		AddSAMBufferToQueue(SAM_Buffer_Queue, buffer);
	}
}

struct
stats
{
	u64 inputHeaderLines;
	u64 inputSAMLines;
	volatile u64 invalidReferenceName;
	volatile u64 unPaired;
	volatile u64 supplementary;
	volatile u64 duplicate;
	volatile u64 qcFail;
	volatile u64 secondary;
	volatile u64 unMapped;
	volatile u64 belowMinMapq;
	volatile u64 tooFarFromRestrictionSite;
	volatile u64 goodRead;
	volatile u64 selfCircle;
	volatile u64 danglingEnd;
	volatile u64 sameFragmentAndStrand;
	volatile u64 reLigated;
	volatile u64 FF;
	volatile u64 FR;
	volatile u64 RF;
	volatile u64 RR;
	volatile u64 *selfCirclePairSequenceVector;
	volatile u64 *danglingEndPairSequenceVector;
	volatile u64 *sameFragmentAndStrandPairSequenceVector;
	volatile u64 *reLigatedPairSequenceVector;
	volatile u64 **FFPairSequenceMatrix;
	volatile u64 **FRPairSequenceMatrix;
	volatile u64 **RFPairSequenceMatrix;
	volatile u64 **RRPairSequenceMatrix;
	volatile u64 FFintraSequence;
	volatile u64 FRintraSequence;
	volatile u64 RFintraSequence;
	volatile u64 RRintraSequence;
	linked_list FFIntraSequenceReadSepValue;
	linked_list FFIntraSequenceFragmentSepValue;
	linked_list FFRestrictionSiteDistanceValue;
	linked_list FRIntraSequenceReadSepValue;
	linked_list FRIntraSequenceFragmentSepValue;
	linked_list FRRestrictionSiteDistanceValue;
	linked_list RFIntraSequenceReadSepValue;
	linked_list RFIntraSequenceFragmentSepValue;
	linked_list RFRestrictionSiteDistanceValue;
	linked_list RRIntraSequenceReadSepValue;
	linked_list RRIntraSequenceFragmentSepValue;
	linked_list RRRestrictionSiteDistanceValue;
	linked_list selfCircleRestrictionSiteDistanceValue;
	linked_list danglingEndRestrictionSiteDistanceValue;
	linked_list sameFragmentAndStrandRestrictionSiteDistanceValue;
	linked_list reLigatedRestrictionSiteDistanceValue;
};

global_variable
stats *
Stats;

global_function
void
ProcessPairs(sam_buffer *buffer)
{
	static sam_buffer *read1 = 0;
	static sam_buffer *read2 = 0;

	if (buffer)
	{
		if (buffer->isGood)
		{
			if (buffer->isRead1)
			{
				if (read1)
				{
					read1->type = typeDump;
					AddSAMRecordToProcessQueue(read1->writeQueue, read1, 0);
				}
				read1 = buffer;
			}
			else
			{
				if (read2)
				{
					read2->type = typeDump;
					AddSAMRecordToProcessQueue(read2->writeQueue, read2, 0);
				}
				read2 = buffer;
			}

			if (read1 && read2 && AreStringsEqual((char *)read1->sam, '\t', (char *)read2->sam, '\t'))
			{
				sam_buffer *originalRead1 = read1;
				sam_buffer *originalRead2 = read2;

				u64 pairSeparation = 0;
				u64 fragmentSeparation = 0;
				if (read1->sequence == read2->sequence)
				{
					if (read1->readStart > read2->readStart)
					{
						sam_buffer *tmp = read1;
						read1 = read2;
						read2 = tmp;
					}
					pairSeparation = read2->readStart - read1->readStart;
					if (read1->fragment->next) fragmentSeparation = read2->fragment->value - read1->fragment->next->value;

					// This preserves original pair-order
					/*{
					  read1 = originalRead1;
					  read2 = originalRead2;
					  }*/
				}

				if (read1->fragment == read2->fragment)
				{
					if (read1->isReverse && !read2->isReverse)
					{
						__atomic_add_fetch(&Stats->selfCircle, 1, 0);
						read1->type = typeSelfCircle;
						read2->type = typeSelfCircle;
						PushOntoLinkedList(read1->distanceFromRestrictionSite, &Stats->selfCircleRestrictionSiteDistanceValue, &Working_Set);
						PushOntoLinkedList(read2->distanceFromRestrictionSite, &Stats->selfCircleRestrictionSiteDistanceValue, &Working_Set);
						__atomic_add_fetch(Stats->selfCirclePairSequenceVector + read1->sequence->index, 1, 0);
					}
					else if (!read1->isReverse && read2->isReverse)
					{
						__atomic_add_fetch(&Stats->danglingEnd, 1, 0);
						read1->type = typeDanglingEnd;
						read2->type = typeDanglingEnd;
						PushOntoLinkedList(read1->distanceFromRestrictionSite, &Stats->danglingEndRestrictionSiteDistanceValue, &Working_Set);
						PushOntoLinkedList(read2->distanceFromRestrictionSite, &Stats->danglingEndRestrictionSiteDistanceValue, &Working_Set);
						__atomic_add_fetch(Stats->danglingEndPairSequenceVector + read1->sequence->index, 1, 0);
					}
					else
					{
						__atomic_add_fetch(&Stats->sameFragmentAndStrand, 1, 0);
						read1->type = typeSameFragmentAndStrand;
						read2->type = typeSameFragmentAndStrand;
						PushOntoLinkedList(read1->distanceFromRestrictionSite, &Stats->sameFragmentAndStrandRestrictionSiteDistanceValue, &Working_Set);
						PushOntoLinkedList(read2->distanceFromRestrictionSite, &Stats->sameFragmentAndStrandRestrictionSiteDistanceValue, &Working_Set);
						__atomic_add_fetch(Stats->sameFragmentAndStrandPairSequenceVector + read1->sequence->index, 1, 0);
					}
				}
				else if (read1->fragment->next == read2->fragment)
				{
					__atomic_add_fetch(&Stats->reLigated, 1, 0);
					read1->type = typeReLigated;
					read2->type = typeReLigated;
					PushOntoLinkedList(read1->distanceFromRestrictionSite, &Stats->reLigatedRestrictionSiteDistanceValue, &Working_Set);
					PushOntoLinkedList(read2->distanceFromRestrictionSite, &Stats->reLigatedRestrictionSiteDistanceValue, &Working_Set);
					__atomic_add_fetch(Stats->reLigatedPairSequenceVector + read1->sequence->index, 1, 0);
				}
				else
				{
					sam_buffer *read1_cap = read1;
					sam_buffer *read2_cap = read2;
					auto ProcessValidPairStats = [read1_cap, read2_cap, pairSeparation, fragmentSeparation](volatile u64 **matrix, volatile u64 *sepCount, linked_list *readSepList, linked_list *fragmentSepList, linked_list *siteList)
					{
						sam_buffer *read1 = read1_cap;
						sam_buffer *read2 = read2_cap;

						u32 id1 = read1->sequence->index;
						u32 id2 = read2->sequence->index;
						u32 min = Min(id1, id2);
						u32 max = Max(id1, id2);
						id1 = min;
						id2 = max - min;
						__atomic_add_fetch(matrix[id1] + id2, 1, 0);

						if (!id2)
						{
							PushOntoLinkedList(pairSeparation, readSepList, &Working_Set);
							PushOntoLinkedList(fragmentSeparation, fragmentSepList, &Working_Set);
							__atomic_add_fetch(sepCount, 1, 0);
						}

						PushOntoLinkedList(read1->distanceFromRestrictionSite, siteList, &Working_Set);
						PushOntoLinkedList(read2->distanceFromRestrictionSite, siteList, &Working_Set);
					};

					if (!read1->isReverse && !read2->isReverse)
					{
						__atomic_add_fetch(&Stats->FF, 1, 0);
						read1->type = typeFF;
						read2->type = typeFF;

						ProcessValidPairStats(Stats->FFPairSequenceMatrix, &Stats->FFintraSequence, &Stats->FFIntraSequenceReadSepValue, &Stats->FFIntraSequenceFragmentSepValue, &Stats->FFRestrictionSiteDistanceValue);
					}
					else if (!read1->isReverse && read2->isReverse)
					{
						__atomic_add_fetch(&Stats->FR, 1, 0);
						read1->type = typeFR;
						read2->type = typeFR;

						ProcessValidPairStats(Stats->FRPairSequenceMatrix, &Stats->FRintraSequence, &Stats->FRIntraSequenceReadSepValue, &Stats->FRIntraSequenceFragmentSepValue, &Stats->FRRestrictionSiteDistanceValue);
					}
					else if (read1->isReverse && !read2->isReverse)
					{
						__atomic_add_fetch(&Stats->RF, 1, 0);
						read1->type = typeRF;
						read2->type = typeRF;

						ProcessValidPairStats(Stats->RFPairSequenceMatrix, &Stats->RFintraSequence, &Stats->RFIntraSequenceReadSepValue, &Stats->RFIntraSequenceFragmentSepValue, &Stats->RFRestrictionSiteDistanceValue);
					}
					else
					{
						__atomic_add_fetch(&Stats->RR, 1, 0);
						read1->type = typeRR;
						read2->type = typeRR;

						ProcessValidPairStats(Stats->RRPairSequenceMatrix, &Stats->RRintraSequence, &Stats->RRIntraSequenceReadSepValue, &Stats->RRIntraSequenceFragmentSepValue, &Stats->RRRestrictionSiteDistanceValue);
					}
				}

				read1 = originalRead1;
				read2 = originalRead2;

				AddSAMRecordToProcessQueue(read1->writeQueue, read1, 0);
				AddSAMRecordToProcessQueue(read2->writeQueue, read2, 0);

				read1 = read2 = 0;
			}
		}
		else
		{
			AddSAMRecordToProcessQueue(buffer->writeQueue, buffer, 0);
		}
	}
	else
	{
		if (read1)
		{
			read1->type = typeDump;
			AddSAMRecordToProcessQueue(read1->writeQueue, read1, 0);
		}

		if (read2)
		{
			read2->type = typeDump;
			AddSAMRecordToProcessQueue(read2->writeQueue, read2, 0);
		}
	}
}

global_function
void
ProcessCigar(u08 *cigar, u32 *referenceDelta_out, u32 *queryDelta_out)
{
	u32 refDelta = 0;
	u32 queryDelta = 0;

	u32 numLen = 0;
	while (*cigar != '\t')
	{
		u08 cigChar = *cigar++;
		u32 process = 1;
		u32 delRef = 0;
		u32 delQu = 0;

		switch (cigChar)
		{
			case 'M':
				{
					delQu = 1;
					delRef = 1;
				} break;

			case 'I':
				{
					delQu = 1;
				} break;

			case 'D':
				{
					delRef = 1;
				} break;

			case 'N':
				{
					delRef = 1;
				} break;

			case 'S':
				{
					delQu = 1;
				} break;

			case 'H':
				{

				} break;

			case 'P':
				{

				} break;

			case '=':
				{
					delQu = 1;
					delRef = 1;
				} break;

			case 'X':
				{
					delQu = 1;
					delRef = 1;
				} break;

			default:
				++numLen;
				process = 0;
		}

		if (process)
		{
			u32 n = numLen ? StringToInt(cigar - 1, numLen) : 1;
			numLen = 0;

			if (delRef) refDelta += n;
			if (delQu) queryDelta += n;
		}
	}

	*referenceDelta_out = refDelta;
	*queryDelta_out = queryDelta;
}

global_function
void
ProcessSamRecord(void *in)
{
	sam_buffer *buffer = (sam_buffer *)in;
	buffer->sequence = 0;
	buffer->fragment = 0;
	u08 *line = buffer->sam;

	u32 goodRead = 1;

	// skip qname
	while (*line++ != '\t') {}

	// get flags
	u32 len = 1;
	u32 flags;
	while (*++line != '\t') ++len;
	flags = StringToInt(line, len);
	{
		// unpaired
		if (!(flags & 0x1))
		{
			goodRead = 0;
			__atomic_add_fetch(&Stats->unPaired, 1, 0);
			buffer->type = typeUnpaired;
		}

		// unmapped
		if (flags & 0x4)
		{
			goodRead = 0;
			__atomic_add_fetch(&Stats->unMapped, 1, 0);
			buffer->type = typeUnmapped;
		}

		// secondary
		if (flags & 0x100)
		{
			goodRead = 0;
			__atomic_add_fetch(&Stats->secondary, 1, 0);
			buffer->type = typeSecondary;
		}

		// qcfail
		if (flags & 0x200)
		{
			goodRead = 0;
			__atomic_add_fetch(&Stats->qcFail, 1, 0);
			buffer->type = typeQCFail;
		}

		// duplicate
		if (flags & 0x400)
		{
			goodRead = 0;
			__atomic_add_fetch(&Stats->duplicate, 1, 0);
			buffer->type = typeDuplicate;
		}

		// supplementary
		if (flags & 0x800)
		{
			goodRead = 0;
			__atomic_add_fetch(&Stats->supplementary, 1, 0);
			buffer->type = typeSupplementary;
		}
	}

	if (goodRead)
	{
		if ((buffer->sequence = GetSequenceFromHashTable(Sequence_Hash_Table, (char *)++line)))
		{
			// skip past rname
			while (*line++ != '\t') {}

			// get reference start
			len = 1;
			while (*++line != '\t') ++len;
			u64 refStart = StringToInt64(line, len) - 1;

			// get mapq
			len = 0;
			while (*++line != '\t') ++len;
			if (StringToInt(line, len) >= Min_Map_Quality)
			{
				u32 readLen;
				u32 referenceDelta;
				ProcessCigar(line + 1, &referenceDelta, &readLen);

				u32 isReverse = flags & 0x10;

				u64 readStart = isReverse ? refStart + (u64)referenceDelta - 1 : refStart;
				wavl_node *fragment = WavlTreeFindInterval(buffer->sequence->tree, readStart);

				u32 distanceFromRestrictionSite = (u32)(readStart - fragment->value);
				u32 closerToNextFragment = 0;
				if (fragment->next)
				{
					u32 nextDis = (u32)(fragment->next->value - readStart);
					if (nextDis < distanceFromRestrictionSite)
					{
						distanceFromRestrictionSite = nextDis;
						closerToNextFragment = 1;
					}
				}

				if (distanceFromRestrictionSite < (readLen << 1))
				{
					buffer->distanceFromRestrictionSite = ((s32)distanceFromRestrictionSite) * (closerToNextFragment ? -1 : 1);
					buffer->readStart = readStart;
					buffer->fragment = fragment;
					buffer->isReverse = isReverse ? 1 : 0;
					buffer->isRead1 = flags & 0x40 ? 1 : 0;

					__atomic_add_fetch(&Stats->goodRead, 1, 0);
				}
				else
				{
					goodRead = 0;
					__atomic_add_fetch(&Stats->tooFarFromRestrictionSite, 1, 0);
					buffer->type = typeTooFarFromRestrictionSite;
				}
			}
			else
			{
				goodRead = 0;
				__atomic_add_fetch(&Stats->belowMinMapq, 1, 0);
				buffer->type = typeBelowMinMapQ;
			}
		}
		else
		{
			goodRead = 0;
			__atomic_add_fetch(&Stats->invalidReferenceName, 1, 0);
			buffer->type = typeInvalidReference;
		}
	}

	buffer->isGood = goodRead ? 1 : 0;
	AddSAMRecordToProcessQueue(buffer->processQueue, buffer);
}

global_function
PyObject *
HiLine_Main(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyObject *objParams;
	static const char *kwlist[] = {"params", 0};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char **)kwlist, &objParams))
	{
		return(0);
	}

	PyObject *objSamInput, *objSamOutput, *objFastaInput, *objRestrictionSites, *objNThreads, *objMinMapQ,
			 *objInfo, *objError, *objPName, *objPVersion;

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
				&objRestrictionSites,
				(char *)"restrictionSites"
			},
			{
				&objFastaInput,
				(char *)"fastaInput"
			},
			{
				&objNThreads,
				(char *)"nThreads"
			},
			{
				&objMinMapQ,
				(char *)"minMapQ"
			},
			{
				&objInfo,
				(char *)"info"
			},
			{
				&objError,
				(char *)"error"
			},
			{
				&objPName,
				(char *)"name"
			},
			{
				&objPVersion,
				(char *)"version"
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

	Name = (char *)PyUnicode_1BYTE_DATA(objPName);
	Version = (char *)PyUnicode_1BYTE_DATA(objPVersion);

	PyObject *objName = Py_BuildValue("s", "name");
	char *samPath = (char *)PyUnicode_1BYTE_DATA(PyObject_GetAttr(objSamInput, objName));
	char *fastaPath = (char *)PyUnicode_1BYTE_DATA(PyObject_GetAttr(objFastaInput, objName));
	Py_DECREF(objName);

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
	if (nThreads < 3)
	{
		PrintError("At-least 3 threads are required at a minimum");
		return(0);
	}
	TestForPythonError;
	Py_DECREF(objNThreads);

	CreateMemoryArena(Working_Set, MegaByte(512));
	Thread_Pool = ThreadPoolInit(&Working_Set, nThreads);
	Line_Buffer_Queue = PushStruct(Working_Set, line_buffer_queue);
	InitialiseLineBufferQueue(&Working_Set, Line_Buffer_Queue);

	{
		Output_Handles = PushStruct(Working_Set, output_handles);
		ForLoop(ArrayCount(Handle_Names))
		{
			PyObject *objName = Py_BuildValue("s", Handle_Names[index]);
			PyObject *objAttrib = PyObject_GetAttr(objSamOutput, objName);
			if (!objAttrib)
			{
				PrintError("No output handle \'%s\' found", Handle_Names[index]);
				return(0);
			}
			Output_Handles->handles[index] = PyObject_AsFileDescriptor(objAttrib);
			if (PyErr_Occurred())
			{
				return(0);
			}
			TestForPythonError;
			Py_DECREF(objName);
			Py_DECREF(objAttrib);
		}
		Unique_Output_Handles = CreateUniqueOutputHandles(&Working_Set, Output_Handles);
	}

	Restriction_Sites = PushStruct(Working_Set, restriction_sites);
	PyObject *objRestrictionSites_fast = PySequence_Fast(objRestrictionSites, "restrictionSites needs to be a list");
	TestForPythonError;
	Py_DECREF(objRestrictionSites);

	Restriction_Sites->num = (u32)PySequence_Fast_GET_SIZE(objRestrictionSites_fast) + 2;
	Restriction_Sites->sites = PushArray(Working_Set, restriction_site, Restriction_Sites->num);

	{
		u32 longestLen = 0;
		PyObject *objLocName = Py_BuildValue("s", "loc");
		ForLoop((Restriction_Sites->num - 2))
		{
			PyObject *objRestrictionSite = PySequence_Fast_GET_ITEM(objRestrictionSites_fast, (Py_ssize_t)index);

			PyObject *objLoc = PyObject_GetAttr(objRestrictionSite, objLocName);
			if (!objLoc)
			{
				PyErr_SetString(PyExc_Exception, "No restrictionSite attribute \'loc\' found");
				return(0);
			}
			u32 loc = (u32)PyLong_AsUnsignedLong(objLoc);
			TestForPythonError;
			Py_DECREF(objLoc);

			PyObject *objPattern = PyObject_GetAttr(objRestrictionSite, Py_BuildValue("s", "pattern"));
			if (!objPattern)
			{
				PyErr_SetString(PyExc_Exception, "No restrictionSite attribute \'pattern\' found");
				return(0);
			}
			char *pattern = (char *)PyUnicode_1BYTE_DATA(objPattern);
			TestForPythonError;

			Restriction_Sites->sites[index].pattern = pattern;
			Restriction_Sites->sites[index].restrictionLocation = loc;

			u32 len = 0;
			{
				u08 mode = 1;
				char *str = Restriction_Sites->sites[index].pattern;
				while (*str != '\0')
				{
					switch (*str)
					{
						case '[':
							mode = 0;
							break;

						case ']':
							mode = 1;
							++len;
							break;

						default:
							if (mode) ++len;
					}

					++str;
				}
			}

			if (len < Restriction_Sites->sites[index].restrictionLocation)
			{
				PrintError("Restriction site \'%s\' of length %u cannot cut at site %u", Restriction_Sites->sites[index].pattern, len, Restriction_Sites->sites[index].restrictionLocation);
				exitCode = EXIT_FAILURE;
				return(0);
			}

			longestLen = Max(longestLen, len);

			if (!CompileRestrictionSiteRegex(Restriction_Sites->sites + index))
			{
				exitCode = EXIT_FAILURE;
				return(0);
			}
		}
		Py_DECREF(objLocName);

		Restriction_Sites->sites[Restriction_Sites->num - 2].pattern = "[ATGC]N";
		Restriction_Sites->sites[Restriction_Sites->num - 2].restrictionLocation = 1;
		Restriction_Sites->sites[Restriction_Sites->num - 1].pattern = "N[ATGC]";
		Restriction_Sites->sites[Restriction_Sites->num - 1].restrictionLocation = 1;

		if (!CompileRestrictionSiteRegex(Restriction_Sites->sites + Restriction_Sites->num - 2) || !CompileRestrictionSiteRegex(Restriction_Sites->sites + Restriction_Sites->num - 1))
		{
			exitCode = EXIT_FAILURE;
			return(0);
		}

		Py_DECREF(objRestrictionSites_fast);
		Restriction_Sites->longestLen = Max(longestLen, 1);
	}

	Min_Map_Quality = (u32)PyLong_AsUnsignedLong(objMinMapQ);
	if (PyErr_Occurred())
	{
		return(0);
	}
	TestForPythonError;
	Py_DECREF(objMinMapQ);

	PyObject *returnObj;
	{
		s32 samFile;
		if ((samFile = PyObject_AsFileDescriptor(objSamInput)) > 0)
		{
#define ReadBufferSize MegaByte(16)
			u08 *readBuffer = PushArray(Working_Set, u08, ReadBufferSize);
			s64 bytesRead;

			s32 fastaFile;
			if ((fastaFile = PyObject_AsFileDescriptor(objFastaInput)) > 0)
			{
				sequence *sequences = 0;
				u32 numSequences = 0;

				Py_BEGIN_ALLOW_THREADS;
				PrintStatus("Reading \'%s\'...", fastaPath);
				{
					sequence *tailSequence = 0;

					line_buffer *buffer = 0;
					u32 bufferPtr = 0;
					u08 character;
					u64 startCoord = 0;
					u32 overlap = Restriction_Sites->longestLen - 1;
					u08 *overlapBuffer = PushArray(Working_Set, u08, overlap);
					overlapBuffer[0] = 0;
					char *currName = 0;

					char nameBuffer[128];
					u32 nameBufferPtr = 0;

					enum enumMode {modeName, modeComment, modeSeq};
					enumMode mode = modeSeq;

					do
					{
						bytesRead = read(fastaFile, readBuffer, ReadBufferSize);
						if (bytesRead < 0)
						{
							PrintError("Error reading from \'%s\'", samPath);
							exitCode = EXIT_FAILURE;
							goto EndFastaRead;
						}

						for (	u64 bufferIndex = 0;
								bufferIndex < (u64)bytesRead;
								++bufferIndex )
						{
							character = readBuffer[bufferIndex];

							if (Global_Regex_Scan_Error_Flag)
							{
								exitCode = EXIT_FAILURE;
								goto EndFastaRead;
							}

							if (mode == modeName || mode == modeComment)
							{
								if (character == '\n')
								{
									mode = modeSeq;

									char *newName = PushArray(Working_Set, char, nameBufferPtr + 1);
									ForLoop(nameBufferPtr)
									{
										newName[index] = nameBuffer[index];
									}
									newName[nameBufferPtr + 1] = '\0';

									sequence *newSequence = PushStruct(Working_Set, sequence);
									newSequence->name = newName;
									newSequence->tree = InitialiseWavlTree(&Working_Set);
									newSequence->next = 0;
									newSequence->index = numSequences++;

									if (!sequences) sequences = newSequence;
									else tailSequence->next = newSequence;

									sequence *prevSequence = tailSequence;
									tailSequence = newSequence;

									currName = newName;
									overlapBuffer[0] = 0;
									bufferPtr = 0;
									buffer = 0;

									if (prevSequence)
									{
										ThreadPoolWait(Thread_Pool);

										create_tree_data *createTreeData = PushStructP(Restriction_Sites->sites->queue.queues[0]->front->arena, create_tree_data);
										createTreeData->arena = PushStruct(Working_Set, memory_arena);
										CreateMemoryArenaP(createTreeData->arena, MegaByte(1));
										prevSequence->length = startCoord + overlap;
										createTreeData->sequence = prevSequence;
										createTreeData->numValueLists = Number_Of_Line_Buffer_Queues * Number_Of_Line_Buffers_Per_Queue;
										createTreeData->values = PushArrayP(Restriction_Sites->sites->queue.queues[0]->front->arena, line_buffer_value *, createTreeData->numValueLists);

										u32 index2 = 0;
										ForLoop(Number_Of_Line_Buffer_Queues)
										{
											single_line_buffer_queue *queue = Line_Buffer_Queue->queues[index];
											line_buffer *lineBuffer = queue->front;

											while (lineBuffer)
											{
												createTreeData->values[index2++] = lineBuffer->headValue;
												lineBuffer->headValue = 0;
												lineBuffer->tailValue = 0;
												lineBuffer = lineBuffer->prev;
											}
										}

										ThreadPoolAddTask(Thread_Pool, CreateTree_Thread, createTreeData);
									}

									startCoord = 0;
									continue;
								}

								if (character < 33) mode = modeComment;

								if (mode == modeComment) continue;

								if (nameBufferPtr == sizeof(nameBuffer))
								{
									PrintError("Seq name \'%s\' too long", nameBuffer);
									exitCode = EXIT_FAILURE;
									goto EndFastaRead;
								}

								nameBuffer[nameBufferPtr++] = (char)character;

								continue;
							}

							if (character == '\n') continue;
							if (character == '>')
							{
								mode = modeName;
								nameBufferPtr = 0;

								if (buffer)
								{
									buffer->lineLength = bufferPtr;
									buffer->startCoord = startCoord;
									startCoord += (u64)(bufferPtr - overlap);
									buffer->name = currName;
									ThreadPoolAddTask(Thread_Pool, ScanForRestrictionSite, buffer);
								}

								continue;
							}

							if (!buffer)
							{
								buffer = TakeLineBufferFromQueue_Wait(Line_Buffer_Queue);
								if (overlapBuffer[0])
								{
									ForLoop(overlap)
									{
										buffer->line[bufferPtr++] = overlapBuffer[index];
									}
								}
							}
							buffer->line[bufferPtr++] = character;

							if (bufferPtr == Line_Buffer_Size)
							{
								buffer->lineLength = bufferPtr;
								buffer->startCoord = startCoord;
								buffer->name = currName;

								ForLoop(overlap)
								{
									overlapBuffer[index] = buffer->line[bufferPtr - overlap + index];
								}

								ThreadPoolAddTask(Thread_Pool, ScanForRestrictionSite, buffer);

								startCoord += (u64)(bufferPtr - overlap);
								bufferPtr = 0;
								buffer = 0;
							}
						}
					} while (bytesRead);

					if (bufferPtr)
					{
						buffer->lineLength = bufferPtr;
						buffer->startCoord = startCoord;
						startCoord += (u64)(bufferPtr - overlap);
						buffer->name = currName;

						ThreadPoolAddTask(Thread_Pool, ScanForRestrictionSite, buffer);
					}

					{
						ThreadPoolWait(Thread_Pool);

						create_tree_data *createTreeData = PushStructP(Restriction_Sites->sites->queue.queues[0]->front->arena, create_tree_data);
						createTreeData->arena = PushStruct(Working_Set, memory_arena);
						CreateMemoryArenaP(createTreeData->arena, MegaByte(1));
						tailSequence->length = startCoord + overlap;
						createTreeData->sequence = tailSequence;
						createTreeData->numValueLists = Number_Of_Line_Buffer_Queues * Number_Of_Line_Buffers_Per_Queue;
						createTreeData->values = PushArrayP(Restriction_Sites->sites->queue.queues[0]->front->arena, line_buffer_value *, createTreeData->numValueLists);

						u32 index2 = 0;
						ForLoop(Number_Of_Line_Buffer_Queues)
						{
							single_line_buffer_queue *queue = Line_Buffer_Queue->queues[index];
							line_buffer *lineBuffer = queue->front;

							while (lineBuffer)
							{
								createTreeData->values[index2++] = lineBuffer->headValue;
								lineBuffer->headValue = 0;
								lineBuffer->tailValue = 0;
								lineBuffer = lineBuffer->prev;
							}
						}

						ThreadPoolAddTask(Thread_Pool, CreateTree_Thread, createTreeData);
					}

EndFastaRead:
					if (exitCode == EXIT_FAILURE) goto End;
				}

				{
					u64 totalBP = 0;

					Sequence_Hash_Table = CreateSequenceHashTable(&Working_Set, numSequences);
					TraverseLinkedList(sequences, sequence)
					{
						AddSequenceToHashTable(&Working_Set, Sequence_Hash_Table, node);
						totalBP += node->length;
					}

					ThreadPoolWait(Thread_Pool);

					u64 totalRestrictionSites = 0;
					u64 treeId = 0;
					TraverseLinkedList(sequences, sequence)
					{
						totalRestrictionSites += (node->tree->numIntervals - 1);

						for (	wavl_node *fragment = WavlTreeFindInterval(node->tree, 0);
								fragment;
								fragment = fragment->next )
						{

							char tagBuffer[256];
							u32 tagLength = 0;

							tagBuffer[tagLength++] = 'r';
							tagBuffer[tagLength++] = 'f';
							tagBuffer[tagLength++] = ':';
							tagBuffer[tagLength++] = 'B';
							tagBuffer[tagLength++] = ':';
							tagBuffer[tagLength++] = 'I';
							tagBuffer[tagLength++] = ',';

							tagLength += (u32)stbsp_snprintf(((char *)tagBuffer) + tagLength, sizeof(tagBuffer) - tagLength - 1, "%u", (u32)(treeId++ & 0xffffffff));
							tagBuffer[tagLength++] = ',';
							tagLength += (u32)stbsp_snprintf(((char *)tagBuffer) + tagLength, sizeof(tagBuffer) - tagLength - 1, "%u", (u32)(fragment->value & 0xffffffff));
							tagBuffer[tagLength++] = ',';
							u64 upstreamFragmentDistance = fragment->next ? fragment->next->value : node->length;
							tagLength += (u32)stbsp_snprintf(((char *)tagBuffer) + tagLength, sizeof(tagBuffer) - tagLength - 1, "%u", (u32)(upstreamFragmentDistance & 0xffffffff));

							fragment->tag = PushArray(Working_Set, u08, tagLength);
							ForLoop(tagLength)
							{
								fragment->tag[index] = tagBuffer[index];
							}
							fragment->tagLength = tagLength;
						}
					}

					PrintStatus("%u sequences scanned for restriction sites, %$.2f sites/sequence, %$.2f BP/site",
							numSequences, (f64)totalRestrictionSites / (f64)numSequences, (f64)totalBP / (f64)totalRestrictionSites);
				}

				{

					ThreadPoolDestroy(Thread_Pool);
					u32 nWorkerThreads = nThreads - 2;
					if (nWorkerThreads > 2) nWorkerThreads = 2;
					Thread_Pool = ThreadPoolInit(&Working_Set, nWorkerThreads);

					thread_pool *processPool = ThreadPoolInit(&Working_Set, 2);

					process_queue *processQueue = CreateProcessQueue(&Working_Set, ProcessPairs, processPool);
					process_queue *writeQueue = CreateProcessQueue(&Working_Set, WriteSamRecord, processPool);

					SAM_Buffer_Queue = PushStruct(Working_Set, sam_buffer_queue);
					InitialiseSAMBufferQueue(&Working_Set, SAM_Buffer_Queue, processQueue, writeQueue);

					Stats = PushStruct(Working_Set, stats);
					memset(Stats, 0, sizeof(stats));

					{
						auto CreateMatrix = [numSequences](volatile u64 ***matrixPtr)
						{
							*matrixPtr = PushArray(Working_Set, volatile u64 *, numSequences);
							ForLoop(numSequences)
							{
								(*matrixPtr)[index] = PushArray(Working_Set, volatile u64, numSequences - index);
								memset((void *)(*matrixPtr)[index], 0, sizeof(u64) * (numSequences - index));
							}
						};

						CreateMatrix(&Stats->FFPairSequenceMatrix);
						CreateMatrix(&Stats->FRPairSequenceMatrix);
						CreateMatrix(&Stats->RFPairSequenceMatrix);
						CreateMatrix(&Stats->RRPairSequenceMatrix);
					}

					{
						auto CreateVector = [numSequences](volatile u64 **vectorPtr)
						{
							*vectorPtr = PushArray(Working_Set, volatile u64, numSequences);
							memset((void *)(*vectorPtr), 0, sizeof(u64) * numSequences);
						};

						CreateVector(&Stats->selfCirclePairSequenceVector);
						CreateVector(&Stats->danglingEndPairSequenceVector);
						CreateVector(&Stats->sameFragmentAndStrandPairSequenceVector);
						CreateVector(&Stats->reLigatedPairSequenceVector);
					}

					sam_buffer *buffer = 0;
					u32 bufferPtr = 0;

					u64 samIndex = 0;
					u64 nHeaderLines = 0;
					u32 doneHeader = 0;
					u64 totalRead = 0;

					char printNBuffers[2][16] = {{0}};
					u08 printNBufferPtr = 0;
					
					PrintStatus("Reading from \'%s\'...", samPath);
					do
					{
						bytesRead = read(samFile, readBuffer, ReadBufferSize);
						if (bytesRead < 0)
						{
							PrintError("Error reading from \'%s\'", samPath);
							exitCode = EXIT_FAILURE;
							goto EndSamRead;
						}

						for (	u64 bufferIndex = 0;
								bufferIndex < (u64)bytesRead;
								++bufferIndex )
						{
							if (Global_SAM_Write_Error_Flag)
							{
								exitCode = EXIT_FAILURE;
								goto EndSamRead;
							}

							u08 character = readBuffer[bufferIndex];

							if (!buffer) buffer = TakeSAMBufferFromQueue_Wait(SAM_Buffer_Queue);

							if (bufferPtr == SAM_Buffer_Size)
							{
								PrintError("Buffer error, a SAM record is larger than %lu bytes", SAM_Buffer_Size);
								exitCode = EXIT_FAILURE;
								goto EndSamRead;
							}

							buffer->sam[bufferPtr++] = character;

							if (character == '\n')
							{
								if (!doneHeader && buffer->sam[0] != '@')
								{
									doneHeader = 1;
									nHeaderLines = samIndex;
									samIndex = 0;
								}

								buffer->lineLength = bufferPtr;
								buffer->samIndex = samIndex++;
								if (doneHeader)
								{
									ThreadPoolAddTask(Thread_Pool, ProcessSamRecord, buffer);
								}
								else
								{
									AddSAMRecordToProcessQueue(buffer->writeQueue, buffer);
								}

#define Log2_Print_Interval 14
								if (!(++totalRead & ((1 << Log2_Print_Interval) - 1)))
								{
									u08 currPtr = printNBufferPtr;
									u08 otherPtr = (currPtr + 1) & 1;
									stbsp_snprintf(printNBuffers[currPtr], sizeof(printNBuffers[currPtr]), "%$" PRIu64, totalRead);

									if (strcmp(printNBuffers[currPtr], printNBuffers[otherPtr]))
									{
										PrintStatus("%s SAM lines processed", printNBuffers[currPtr]);
									}
									
									printNBufferPtr = otherPtr;
								}

								buffer = 0;
								bufferPtr = 0;
							}
						}

					} while (bytesRead);

					Stats->inputHeaderLines = nHeaderLines;
					Stats->inputSAMLines = samIndex;

EndSamRead:
					ThreadPoolWait(Thread_Pool);

					FenceIn(CloseProcessQueue(processQueue));
					FenceIn(CloseProcessQueue(writeQueue, 0));

					ThreadPoolWait(processPool);
					ThreadPoolDestroy(processPool);
				}
				Py_END_ALLOW_THREADS;

				if (exitCode == EXIT_SUCCESS)
				{
					PrintStatus("");
					PrintStatus("%$" PRIu64 "/%$" PRIu64 " SAM header/record lines read", Stats->inputHeaderLines, Stats->inputSAMLines);
					PrintStatus("%$" PRIu64 " good reads", Stats->goodRead);
					PrintStatus("%$" PRIu64 " unpaired reads", Stats->unPaired);
					PrintStatus("%$" PRIu64 " supplementary reads", Stats->supplementary);
					PrintStatus("%$" PRIu64 " duplicate reads", Stats->duplicate);
					PrintStatus("%$" PRIu64 " qc-failing reads", Stats->qcFail);
					PrintStatus("%$" PRIu64 " secondary reads", Stats->secondary);
					PrintStatus("%$" PRIu64 " unmapped reads", Stats->unMapped);
					PrintStatus("%$" PRIu64 " reads below minimum mapping quality (%u)", Stats->belowMinMapq, Min_Map_Quality);
					PrintStatus("%$" PRIu64 " reads too far from a restriction site", Stats->tooFarFromRestrictionSite);
					PrintStatus("%$" PRIu64 " reads with an invalid reference name", Stats->invalidReferenceName);

					PrintStatus("");
					u64 totalValidPairs = Stats->FF + Stats->FR + Stats->RF + Stats->RR;
					PrintStatus("%$" PRIu64 " valid HiC pairs", totalValidPairs);
					PrintStatus("%$.2f%% FF", 100.0 * (f64)Stats->FF / (f64)totalValidPairs);
					PrintStatus("%$.2f%% FR", 100.0 * (f64)Stats->FR / (f64)totalValidPairs);
					PrintStatus("%$.2f%% RF", 100.0 * (f64)Stats->RF / (f64)totalValidPairs);
					PrintStatus("%$.2f%% RR", 100.0 * (f64)Stats->RR / (f64)totalValidPairs);
					PrintStatus("%$" PRIu64 " intra-sequence pairs", Stats->FFintraSequence + Stats->FRintraSequence + Stats->RFintraSequence + Stats->RRintraSequence);
					PrintStatus("%$" PRIu64 " self-circle pairs", Stats->selfCircle);
					PrintStatus("%$" PRIu64 " dangling-end pairs", Stats->danglingEnd);
					PrintStatus("%$" PRIu64 " re-ligated pairs", Stats->reLigated);
					PrintStatus("%$" PRIu64 " error pairs (same fragment and strand)", Stats->sameFragmentAndStrand);
					PrintStatus("");

					{
						npy_intp dims[1] = {18};
						PyArrayObject *typeCounts = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);

						{
							u64 *out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (0 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->invalidReferenceName);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (1 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->unPaired);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (2 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->supplementary);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (3 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->duplicate);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (4 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->qcFail);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (5 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->secondary);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (6 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->unMapped);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (7 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->belowMinMapq);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (8 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->tooFarFromRestrictionSite);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (9 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->goodRead);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (10 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->selfCircle);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (11 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->danglingEnd);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (12 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->sameFragmentAndStrand);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (13 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->reLigated);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (14 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->FF);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (15 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->FR);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (16 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->RF);

							out = (u64 *)((u08 *)PyArray_DATA(typeCounts) + (17 * PyArray_STRIDE(typeCounts, 0)));
							*out = (u64)(Stats->RR);
						}

						PyObject *seqNames = PyTuple_New((Py_ssize_t)numSequences);
						dims[0] = numSequences;
						PyArrayObject *seqLengths = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *restrictionSitesPerSeq = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						{
							u32 index = 0;
							TraverseLinkedList(sequences, sequence)
							{
								PyTuple_SET_ITEM(seqNames, (Py_ssize_t)index, Py_BuildValue("s", node->name));
								u64 *out = (u64 *)((u08 *)PyArray_DATA(seqLengths) + (index * PyArray_STRIDE(seqLengths, 0)));
								*out = node->length;
								out = (u64 *)((u08 *)PyArray_DATA(restrictionSitesPerSeq) + (index++ * PyArray_STRIDE(restrictionSitesPerSeq, 0)));
								*out = node->tree->numIntervals - 1;
							}
						}

						PyArrayObject *selfCirclePairsPerSeq = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *danglingEndPairsPerSeq = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *sameFragmentAndStrandPairsPerSeq = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *reLigatedPairsPerSeq = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						{
							ForLoop(numSequences)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(selfCirclePairsPerSeq) + (index * PyArray_STRIDE(selfCirclePairsPerSeq, 0)));
								*out = (u64)Stats->selfCirclePairSequenceVector[index];

								out = (u64 *)((u08 *)PyArray_DATA(danglingEndPairsPerSeq) + (index * PyArray_STRIDE(danglingEndPairsPerSeq, 0)));
								*out = (u64)Stats->danglingEndPairSequenceVector[index];

								out = (u64 *)((u08 *)PyArray_DATA(sameFragmentAndStrandPairsPerSeq) + (index * PyArray_STRIDE(sameFragmentAndStrandPairsPerSeq, 0)));
								*out = (u64)Stats->sameFragmentAndStrandPairSequenceVector[index];

								out = (u64 *)((u08 *)PyArray_DATA(reLigatedPairsPerSeq) + (index * PyArray_STRIDE(reLigatedPairsPerSeq, 0)));
								*out = (u64)Stats->reLigatedPairSequenceVector[index];
							}
						}

						PyObject *FFPerSeq = PyTuple_New((Py_ssize_t)numSequences);
						PyObject *FRPerSeq = PyTuple_New((Py_ssize_t)numSequences);
						PyObject *RFPerSeq = PyTuple_New((Py_ssize_t)numSequences);
						PyObject *RRPerSeq = PyTuple_New((Py_ssize_t)numSequences);
						{
							ForLoop(numSequences)
							{
								dims[0] = numSequences - index;
								PyArrayObject *arrayFF = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
								PyArrayObject *arrayFR = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
								PyArrayObject *arrayRF = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
								PyArrayObject *arrayRR = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);

								ForLoop2((numSequences - index))
								{
									u64 *out = (u64 *)((u08 *)PyArray_DATA(arrayFF) + (index2 * PyArray_STRIDE(arrayFF, 0)));
									*out = (u64)Stats->FFPairSequenceMatrix[index][index2];

									out = (u64 *)((u08 *)PyArray_DATA(arrayFR) + (index2 * PyArray_STRIDE(arrayFR, 0)));
									*out = (u64)Stats->FRPairSequenceMatrix[index][index2];

									out = (u64 *)((u08 *)PyArray_DATA(arrayRF) + (index2 * PyArray_STRIDE(arrayRF, 0)));
									*out = (u64)Stats->RFPairSequenceMatrix[index][index2];

									out = (u64 *)((u08 *)PyArray_DATA(arrayRR) + (index2 * PyArray_STRIDE(arrayRR, 0)));
									*out = (u64)Stats->RRPairSequenceMatrix[index][index2];
								}

								PyTuple_SET_ITEM(FFPerSeq, (Py_ssize_t)index, (PyObject *)arrayFF);
								PyTuple_SET_ITEM(FRPerSeq, (Py_ssize_t)index, (PyObject *)arrayFR);
								PyTuple_SET_ITEM(RFPerSeq, (Py_ssize_t)index, (PyObject *)arrayRF);
								PyTuple_SET_ITEM(RRPerSeq, (Py_ssize_t)index, (PyObject *)arrayRR);
							}
						}

						dims[0] = 2 * Stats->selfCircle;
						PyArrayObject *selfCircleRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->selfCircleRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(selfCircleRestrictionSiteDistances) + (index++ * PyArray_STRIDE(selfCircleRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = 2 * Stats->danglingEnd;
						PyArrayObject *danglingEndRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->danglingEndRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(danglingEndRestrictionSiteDistances) + (index++ * PyArray_STRIDE(danglingEndRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = 2 * Stats->sameFragmentAndStrand;
						PyArrayObject *sameFragmentAndStrandRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->sameFragmentAndStrandRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(sameFragmentAndStrandRestrictionSiteDistances) + (index++ * PyArray_STRIDE(sameFragmentAndStrandRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = 2 * Stats->reLigated;
						PyArrayObject *reLigatedRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->reLigatedRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(reLigatedRestrictionSiteDistances) + (index++ * PyArray_STRIDE(reLigatedRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = 2 * Stats->FF;
						PyArrayObject *FFRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->FFRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(FFRestrictionSiteDistances) + (index++ * PyArray_STRIDE(FFRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = 2 * Stats->FR;
						PyArrayObject *FRRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->FRRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(FRRestrictionSiteDistances) + (index++ * PyArray_STRIDE(FRRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = 2 * Stats->RF;
						PyArrayObject *RFRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->RFRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(RFRestrictionSiteDistances) + (index++ * PyArray_STRIDE(RFRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = 2 * Stats->RR;
						PyArrayObject *RRRestrictionSiteDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_INT32);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->RRRestrictionSiteDistanceValue.head, linked_list_node)
							{
								s32 *out = (s32 *)((u08 *)PyArray_DATA(RRRestrictionSiteDistances) + (index++ * PyArray_STRIDE(RRRestrictionSiteDistances, 0)));
								*out = (s32)node->signedValue;
							}
						}

						dims[0] = Stats->FFintraSequence;
						PyArrayObject *FFIntraSequenceReadDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *FFIntraSequenceFragmentDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->FFIntraSequenceReadSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(FFIntraSequenceReadDistances) + (index++ * PyArray_STRIDE(FFIntraSequenceReadDistances, 0)));
								*out = (u64)node->value;
							}
							index = 0;
							TraverseLinkedList(Stats->FFIntraSequenceFragmentSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(FFIntraSequenceFragmentDistances) + (index++ * PyArray_STRIDE(FFIntraSequenceFragmentDistances, 0)));
								*out = (u64)node->value;
							}
						}

						dims[0] = Stats->FRintraSequence;
						PyArrayObject *FRIntraSequenceReadDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *FRIntraSequenceFragmentDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->FRIntraSequenceReadSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(FRIntraSequenceReadDistances) + (index++ * PyArray_STRIDE(FRIntraSequenceReadDistances, 0)));
								*out = (u64)node->value;
							}
							index = 0;
							TraverseLinkedList(Stats->FRIntraSequenceFragmentSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(FRIntraSequenceFragmentDistances) + (index++ * PyArray_STRIDE(FRIntraSequenceFragmentDistances, 0)));
								*out = (u64)node->value;
							}
						}

						dims[0] = Stats->RFintraSequence;
						PyArrayObject *RFIntraSequenceReadDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *RFIntraSequenceFragmentDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->RFIntraSequenceReadSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(RFIntraSequenceReadDistances) + (index++ * PyArray_STRIDE(RFIntraSequenceReadDistances, 0)));
								*out = (u64)node->value;
							}
							index = 0;
							TraverseLinkedList(Stats->RFIntraSequenceFragmentSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(RFIntraSequenceFragmentDistances) + (index++ * PyArray_STRIDE(RFIntraSequenceFragmentDistances, 0)));
								*out = (u64)node->value;
							}
						}

						dims[0] = Stats->RRintraSequence;
						PyArrayObject *RRIntraSequenceReadDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						PyArrayObject *RRIntraSequenceFragmentDistances = (PyArrayObject *)PyArray_SimpleNew(1, (const npy_intp*)dims, NPY_UINT64);
						{
							u64 index = 0;
							TraverseLinkedList(Stats->RRIntraSequenceReadSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(RRIntraSequenceReadDistances) + (index++ * PyArray_STRIDE(RRIntraSequenceReadDistances, 0)));
								*out = (u64)node->value;
							}
							index = 0;
							TraverseLinkedList(Stats->RRIntraSequenceFragmentSepValue.head, linked_list_node)
							{
								u64 *out = (u64 *)((u08 *)PyArray_DATA(RRIntraSequenceFragmentDistances) + (index++ * PyArray_STRIDE(RRIntraSequenceFragmentDistances, 0)));
								*out = (u64)node->value;
							}
						}

						returnObj = PyTuple_Pack(28, typeCounts, seqNames, seqLengths, restrictionSitesPerSeq, selfCirclePairsPerSeq, danglingEndPairsPerSeq, sameFragmentAndStrandPairsPerSeq, reLigatedPairsPerSeq, FFPerSeq, FRPerSeq, RFPerSeq, RRPerSeq, selfCircleRestrictionSiteDistances, danglingEndRestrictionSiteDistances, sameFragmentAndStrandRestrictionSiteDistances, reLigatedRestrictionSiteDistances, FFRestrictionSiteDistances, FRRestrictionSiteDistances, RFRestrictionSiteDistances, RRRestrictionSiteDistances, FFIntraSequenceReadDistances, FRIntraSequenceReadDistances, RFIntraSequenceReadDistances, RRIntraSequenceReadDistances, FFIntraSequenceFragmentDistances, FRIntraSequenceFragmentDistances, RFIntraSequenceFragmentDistances, RRIntraSequenceFragmentDistances);
					}
				}
			}
			else
			{
				PrintError("Could not open \'%s\'", fastaPath);
				exitCode = EXIT_FAILURE;
				goto End;
			}

			if (exitCode == EXIT_FAILURE) goto End;
		}
		else
		{
			PrintError("Could not open \'%s\'", samPath);
			exitCode = EXIT_FAILURE;
			goto End;
		}

		if (exitCode == EXIT_FAILURE) goto End;
	}

End:
	if (Thread_Pool)
	{
		ThreadPoolWait(Thread_Pool);
		ThreadPoolDestroy(Thread_Pool);
	}

	FreeRestrictionSites(Restriction_Sites);
	FreeMemoryArena(Working_Set);

	if (exitCode == EXIT_FAILURE) return(0);
	return(returnObj);
}

#define ModuleName "_HiLine_Main"
#define ModuleDescription "A HiC alignment and classification pipeline"

global_variable
PyMethodDef
HiLineMethods[] =
{
	{ModuleName, (PyCFunction) HiLine_Main, METH_VARARGS | METH_KEYWORDS, ModuleDescription},
	{NULL, NULL, 0, NULL}
};

global_variable
struct
PyModuleDef
HiLineModule =
{
	PyModuleDef_HEAD_INIT,
	ModuleName,
	ModuleDescription,
	-1,
	HiLineMethods
};

PyMODINIT_FUNC
PyInit__HiLine()
{
	PyObject *obj = PyModule_Create(&HiLineModule);
	import_array();
	return(obj);
}
