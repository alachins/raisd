/*  
 *  RAiSD, Raised Accuracy in Sweep Detection
 *
 *  Copyright January 2017 by Nikolaos Alachiotis and Pavlos Pavlidis
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  For any other enquiries send an email to
 *  Nikolaos Alachiotis (n.alachiotis@gmail.com)
 *  Pavlos Pavlidis (pavlidisp@gmail.com)  
 *  
 */

#include "RAiSD.h"

RSDChunk_t * RSDChunk_new(void)
{
	RSDChunk_t * ch = NULL;
	ch = (RSDChunk_t *) malloc(sizeof(RSDChunk_t));
	assert(ch!=NULL);

	ch->chunkID = -1;
	//posPosition (not init here)
	ch->seqPosition = NULL;
	
	ch->chunkMemSize = CHUNK_MEMSIZE_AND_INCREMENT;

	ch->chunkSize = 0;
	ch->sitePosition = NULL;
	ch->derivedAlleleCount = NULL;
	ch->patternID = NULL;	

	return ch;
}

void RSDChunk_free(RSDChunk_t * ch, int numberOfSamples)
{
	assert(ch!=NULL);

	if(numberOfSamples!=-1)
		MemoryFootprint += sizeof(fpos_t)*numberOfSamples;
	
	if(ch->seqPosition!=NULL)
		free(ch->seqPosition);

	MemoryFootprint += sizeof(float)*ch->chunkMemSize;

	if(ch->sitePosition!=NULL)
		free(ch->sitePosition);

	MemoryFootprint += sizeof(int)*ch->chunkMemSize;
	
	if(ch->derivedAlleleCount!=NULL)
		free(ch->derivedAlleleCount);

	MemoryFootprint += sizeof(int)*ch->chunkMemSize;

	if(ch->patternID!=NULL)
		free(ch->patternID);

	free(ch);
}

void RSDChunk_init(RSDChunk_t * RSDChunk, int numberOfSamples)
{
	RSDChunk->seqPosition = (fpos_t*)malloc(sizeof(fpos_t)*numberOfSamples);
	assert(RSDChunk->seqPosition!=NULL);

	RSDChunk->sitePosition = (float*)malloc(sizeof(float)*RSDChunk->chunkMemSize);
	assert(RSDChunk->sitePosition != NULL);

	RSDChunk->derivedAlleleCount = (int*)malloc(sizeof(int)*RSDChunk->chunkMemSize);
	assert(RSDChunk->derivedAlleleCount != NULL);

	RSDChunk->patternID = (int*)malloc(sizeof(int)*RSDChunk->chunkMemSize);
	assert(RSDChunk->patternID != NULL);

	RSDChunk->derAll1CntTotal = 0;
	RSDChunk->derAllNCntTotal = 0;
}

void RSDChunk_reset(RSDChunk_t * RSDChunk)
{
	int i;

	if(RSDChunk->chunkID==-1)
	{
		for(i=0;i<RSDChunk->chunkMemSize;i++)
		{
			RSDChunk->sitePosition[i] = 0.0;
			RSDChunk->derivedAlleleCount[i] = 0;
			RSDChunk->patternID[i] = 0;
		}

		RSDChunk->chunkSize = 0;
	}
	else
	{
		assert(RSDChunk->chunkID>=0);

		if(RSDChunk->chunkSize <= WINDOW_SIZE)
			return;

		int i_offset = RSDChunk->chunkSize - WINDOW_SIZE + 1;
		assert(i_offset>=0);

		RSDChunk->chunkSize = WINDOW_SIZE - 1;

		for(i=0;i<RSDChunk->chunkSize;i++)
		{
			RSDChunk->sitePosition[i] = RSDChunk->sitePosition[i_offset+i];
			RSDChunk->derivedAlleleCount[i] = RSDChunk->derivedAlleleCount[i_offset+i];
			RSDChunk->patternID[i] = RSDChunk->patternID[i_offset+i];
		}
	
		for(;i<RSDChunk->chunkMemSize;i++)
		{
			RSDChunk->sitePosition[i] = 0.0;
			RSDChunk->derivedAlleleCount[i] = 0;
			RSDChunk->patternID[i] = 0;
		}
	}	
}
