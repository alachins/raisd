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
	ch = (RSDChunk_t *) rsd_malloc(sizeof(RSDChunk_t));
	assert(ch!=NULL);

	ch->chunkID = -1;
	//posPosition (not init here)
	ch->seqPosition = NULL;
	
	ch->chunkMemSize = CHUNK_MEMSIZE_AND_INCREMENT;

	ch->chunkSize = 0;

	ch->sitePosition = NULL;
	ch->derivedAlleleCount = NULL;
	ch->patternID = NULL;	
	
	ch->chunkData = NULL;

	return ch;
}

void RSDChunk_free(RSDChunk_t * ch, int64_t numberOfSamples)
{
	assert(ch!=NULL);
	assert(numberOfSamples!=0);

	//if(numberOfSamples!=-1)
	//	MemoryFootprint += sizeof(fpos_t)*((unsigned long)numberOfSamples);
	
	if(ch->seqPosition!=NULL)
		free(ch->seqPosition);

	//MemoryFootprint += 3*sizeof(int)*((unsigned long)ch->chunkMemSize);

	if(ch->sitePosition!=NULL)
		free(ch->sitePosition);
	
	if(ch->derivedAlleleCount!=NULL)
		free(ch->derivedAlleleCount);

	if(ch->patternID!=NULL)
		free(ch->patternID);

	if(ch->chunkData!=NULL)
		free(ch->chunkData);

	free(ch);
}

void RSDChunk_init(RSDChunk_t * RSDChunk, int64_t numberOfSamples, int64_t createPatternPoolMask)
{
	RSDChunk->seqPosition = (fpos_t*)rsd_malloc(sizeof(fpos_t)*((unsigned long)numberOfSamples));
	assert(RSDChunk->seqPosition!=NULL);

#ifdef _MLT

	assert(0); // sitePosition was changed to double in version 2.5. MLT is in a proof-of-concept state.

	if(createPatternPoolMask==1)
	{
		RSDChunk->sitePosition = (float*)rsd_malloc(sizeof(float)*((unsigned long)RSDChunk->chunkMemSize));
		assert(RSDChunk->sitePosition != NULL);

		RSDChunk->derivedAlleleCount = (int*)rsd_malloc(sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
		assert(RSDChunk->derivedAlleleCount != NULL);

		RSDChunk->patternID = (int*)rsd_malloc(sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
		assert(RSDChunk->patternID != NULL);
	}
	else
	{
		RSDChunk->chunkData = (float*)rsd_malloc(sizeof(float)*((unsigned long)RSDChunk->chunkMemSize*3));
		assert(RSDChunk->chunkData!=NULL);
	}
#else
	assert(createPatternPoolMask==0 || createPatternPoolMask==1);

	RSDChunk->sitePosition = (double*)rsd_malloc(sizeof(double)*((unsigned long)RSDChunk->chunkMemSize));
	assert(RSDChunk->sitePosition != NULL);

	RSDChunk->derivedAlleleCount = (int*)rsd_malloc(sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
	assert(RSDChunk->derivedAlleleCount != NULL);

	RSDChunk->patternID = (int*)rsd_malloc(sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
	assert(RSDChunk->patternID != NULL);
#endif

	RSDChunk->derAll1CntTotal = 0;
	RSDChunk->derAllNCntTotal = 0;
}

void RSDChunk_reset(RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);

	int i;

	if(RSDChunk->chunkID==-1)
	{

#ifdef _MLT
		if(RSDCommandLine->createPatternPoolMask==1)
		{
			for(i=0;i<RSDChunk->chunkMemSize;i++)
			{
				RSDChunk->sitePosition[i] = 0.0;
				RSDChunk->derivedAlleleCount[i] = 0;
				RSDChunk->patternID[i] = 0;
			}
		}
		else
		{
			for(i=0;i<RSDChunk->chunkMemSize*3;i++)
				RSDChunk->chunkData[i] = 0.0f;
		}
#else
		for(i=0;i<RSDChunk->chunkMemSize;i++)
		{
			RSDChunk->sitePosition[i] = 0.0;
			RSDChunk->derivedAlleleCount[i] = 0;
			RSDChunk->patternID[i] = 0;
		}
#endif

		RSDChunk->chunkSize = 0;
	}
	else
	{
		assert(RSDChunk->chunkID>=0);

		if(RSDChunk->chunkSize <= RSDCommandLine->windowSize)
			return;

		int64_t i_offset = RSDChunk->chunkSize - RSDCommandLine->windowSize + 1;
		assert(i_offset>=0);

		RSDChunk->chunkSize = RSDCommandLine->windowSize - 1;

#ifdef _MLT
		if(RSDCommandLine->createPatternPoolMask==1)
		{
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
		else
		{
			for(i=0;i<RSDChunk->chunkSize;i++)
			{
				RSDChunk->chunkData[i*3+0] = RSDChunk->chunkData[(i_offset+i)*3+0];
				RSDChunk->chunkData[i*3+1] = RSDChunk->chunkData[(i_offset+i)*3+1];
				RSDChunk->chunkData[i*3+2] = RSDChunk->chunkData[(i_offset+i)*3+2];
			}

			for(;i<RSDChunk->chunkMemSize;i++)
			{
				RSDChunk->chunkData[i*3+0] = 0.0f;
				RSDChunk->chunkData[i*3+1] = 0.0f;
				RSDChunk->chunkData[i*3+2] = 0.0f;
			}
		}
#else
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
#endif
	}

}
