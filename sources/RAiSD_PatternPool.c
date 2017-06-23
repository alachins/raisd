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

RSDPatternPool_t * RSDPatternPool_new(void)
{
	RSDPatternPool_t * pp = NULL;
	pp = (RSDPatternPool_t *) malloc(sizeof(RSDPatternPool_t));
	assert(pp!=NULL);

	pp->memorySize = PATTERNPOOL_SIZE;
	pp->maxSize = -1;
	pp->patternSize = -1;
	pp->dataSize = -1;
	pp->incomingSite = NULL;
	pp->incomingSiteCompact = NULL;
	pp->incomingSitePosition = -1.0;
	pp->poolData = NULL;
	pp->exchangeBuffer = NULL;

	return pp;
}

void RSDPatternPool_free(RSDPatternPool_t * pp, int numberOfSamples)
{
	assert(pp!=NULL);

	MemoryFootprint += (numberOfSamples+1);
	free(pp->incomingSite);

	MemoryFootprint += sizeof(uint64_t)*pp->patternSize;
	free(pp->incomingSiteCompact);

	MemoryFootprint += sizeof(uint64_t)*(pp->patternSize*pp->maxSize);
	free(pp->poolData);

	MemoryFootprint += sizeof(int)*pp->maxSize;
	free(pp->poolDataAlleleCount);

	MemoryFootprint += sizeof(int)*pp->maxSize;
	free(pp->poolDataPatternCount);

	MemoryFootprint += sizeof(uint64_t)*pp->patternSize;
	free(pp->exchangeBuffer);

	free(pp);	
}

void RSDPatternPool_init (RSDPatternPool_t * RSDPatternPool, RSDCommandLine_t * RSDCommandLine, int numberOfSamples)
{
	assert(RSDCommandLine!=NULL);

	int wordLength = sizeof(uint64_t)*8;

	int wordsPerSNP = numberOfSamples%wordLength==0?numberOfSamples/wordLength:numberOfSamples/wordLength+1;
	RSDPatternPool->patternSize = wordsPerSNP;

	float maxSNPsInPool = ((float)(RSDPatternPool->memorySize)*1024.0*1024.0)/((wordsPerSNP*8.0)+(8.0)); // +8 for the allelecount and the patterncount
	RSDPatternPool->maxSize = (int) maxSNPsInPool;

	RSDPatternPool->incomingSite = (char*) malloc(sizeof(char)*(numberOfSamples+1));
	assert(RSDPatternPool->incomingSite!=NULL);

	RSDPatternPool->incomingSiteCompact = (uint64_t*) malloc(sizeof(uint64_t)*RSDPatternPool->patternSize);
	assert(RSDPatternPool->incomingSiteCompact!=NULL);

	RSDPatternPool->poolData = (uint64_t*) malloc(sizeof(uint64_t)*(RSDPatternPool->patternSize*RSDPatternPool->maxSize));
	assert(RSDPatternPool->poolData!=NULL);

	RSDPatternPool->poolDataAlleleCount = (int*) malloc(sizeof(int)*RSDPatternPool->maxSize);
	assert(RSDPatternPool->poolDataAlleleCount!=NULL);

	RSDPatternPool->poolDataPatternCount = (int*) malloc(sizeof(int)*RSDPatternPool->maxSize);
	assert(RSDPatternPool->poolDataPatternCount!=NULL);

	RSDPatternPool->exchangeBuffer = (uint64_t*) malloc(sizeof(uint64_t)*RSDPatternPool->patternSize);
	assert(RSDPatternPool->exchangeBuffer!=NULL);
}

void RSDPatternPool_resize (RSDPatternPool_t * RSDPatternPool, int setSamples, FILE * fpOut)
{
	if(setSamples<=(int)(RSDPatternPool->patternSize*sizeof(uint64_t)*8))
		return;
	
	int wordLength = sizeof(uint64_t)*8;

	int wordsPerSNP = setSamples%wordLength==0?setSamples/wordLength:setSamples/wordLength+1;
	RSDPatternPool->patternSize = wordsPerSNP;

	float maxSNPsInPool = ((float)(RSDPatternPool->memorySize)*1024.0*1024.0)/((wordsPerSNP*8.0)+(8.0)); // +8 for the allelecount and the patterncount
	RSDPatternPool->maxSize = (int) maxSNPsInPool;

	RSDPatternPool->incomingSiteCompact = (uint64_t*) malloc(sizeof(uint64_t)*RSDPatternPool->patternSize);
	assert(RSDPatternPool->incomingSiteCompact!=NULL);

	RSDPatternPool->poolData = (uint64_t*) malloc(sizeof(uint64_t)*(RSDPatternPool->patternSize*RSDPatternPool->maxSize));
	assert(RSDPatternPool->poolData!=NULL);

	RSDPatternPool->exchangeBuffer = (uint64_t*) malloc(sizeof(uint64_t)*RSDPatternPool->patternSize);
	assert(RSDPatternPool->exchangeBuffer!=NULL);

	fprintf(fpOut, "The pattern structure has been resized to %d patterns (max. capacity) and approx. %d MB memory footprint.\n", RSDPatternPool->maxSize, RSDPatternPool->memorySize);	
	fflush(fpOut);	

	fprintf(stdout, "The pattern structure has been resized to %d patterns (max. capacity) and approx. %d MB memory footprint.\n", RSDPatternPool->maxSize, RSDPatternPool->memorySize);	
	fflush(stdout);
}

void RSDPatternPool_print(RSDPatternPool_t * RSDPatternPool, FILE * fpOut)
{
	fprintf(fpOut, "\nA pattern structure of %d patterns (max. capacity) and approx. %d MB memory footprint has been created.\n", RSDPatternPool->maxSize, RSDPatternPool->memorySize);	
	fflush(fpOut);
}

void RSDPatternPool_exhangePatterns (RSDPatternPool_t * RSDPatternPool, int pID_a, int pID_b)
{
	if(pID_a==pID_b)
		return;

	memcpy(RSDPatternPool->exchangeBuffer, &(RSDPatternPool->poolData[pID_a*RSDPatternPool->patternSize]), RSDPatternPool->patternSize*8);
	int alleleCount =  RSDPatternPool->poolDataAlleleCount[pID_a];
	int patternCount =  0;//RSDPatternPool->poolDataPatternCount[pID_a]; // Not copying the patterncount because it changes per chunk
 
	memcpy(&(RSDPatternPool->poolData[pID_a*RSDPatternPool->patternSize]), &(RSDPatternPool->poolData[pID_b*RSDPatternPool->patternSize]), RSDPatternPool->patternSize*8);
	RSDPatternPool->poolDataAlleleCount[pID_a] = RSDPatternPool->poolDataAlleleCount[pID_b];
	RSDPatternPool->poolDataPatternCount[pID_a] = 0;//RSDPatternPool->poolDataPatternCount[pID_b];

	memcpy(&(RSDPatternPool->poolData[pID_b*RSDPatternPool->patternSize]), RSDPatternPool->exchangeBuffer, RSDPatternPool->patternSize*8);
	RSDPatternPool->poolDataAlleleCount[pID_b] = alleleCount;
	RSDPatternPool->poolDataPatternCount[pID_b] = patternCount;

}

void RSDPatternPool_reset (RSDPatternPool_t * RSDPatternPool, int numberOfSamples, int setSamples, RSDChunk_t * RSDChunk)
{
	int i;

	if(RSDChunk->chunkID==-1)
	{
		if(numberOfSamples>=setSamples)
		{
			for(i=0;i<numberOfSamples+1;i++)
				RSDPatternPool->incomingSite[i]='\0';
		}
		else
		{
			for(i=0;i<setSamples+1;i++) 
				RSDPatternPool->incomingSite[i]='\0';
		}

		RSDPatternPool->incomingSiteAlleleCount = 0;

		for(i=0;i<RSDPatternPool->patternSize;i++)
			RSDPatternPool->incomingSiteCompact[i] = 0ull;

		RSDPatternPool->incomingSitePosition = -1.0;

		for(i=0;i<RSDPatternPool->patternSize*RSDPatternPool->dataSize;i++)
			RSDPatternPool->poolData[i] = 0ull;

		for(i=0;i<RSDPatternPool->dataSize;i++)
		{
			RSDPatternPool->poolDataAlleleCount[i] = 0;
			RSDPatternPool->poolDataPatternCount[i] = 0;
		}

		RSDPatternPool->dataSize = 0;
	}
	else
	{
		if(numberOfSamples>=setSamples)
		{
			for(i=0;i<numberOfSamples+1;i++)
				RSDPatternPool->incomingSite[i]='\0';
		}
		else
		{
			for(i=0;i<setSamples+1;i++) 
				RSDPatternPool->incomingSite[i]='\0';
		}

		RSDPatternPool->incomingSiteAlleleCount = 0;

		for(i=0;i<RSDPatternPool->patternSize;i++)
			RSDPatternPool->incomingSiteCompact[i] = 0ull;

		RSDPatternPool->incomingSitePosition = -1.0;


		// End part of chunk rellocated to the beginning
		int j;
		int pLoc = 0;

		for(i=0;i<RSDChunk->chunkSize;i++)
		{
			int pID_a = RSDChunk->patternID[i];
			int pID_b = pLoc;
		
			if(pID_a >= pID_b)
			{
				RSDPatternPool_exhangePatterns (RSDPatternPool, pID_a, pID_b);

				for(j=i;j<RSDChunk->chunkSize;j++)
				{
					if(RSDChunk->patternID[j]==pID_a)
					{
						RSDChunk->patternID[j] = pID_b;
						RSDPatternPool->poolDataPatternCount[pLoc]++;
					}
					else
					{
						if(RSDChunk->patternID[j]==pID_b)
							RSDChunk->patternID[j] = pID_a;
					}
				}

				pLoc++;
			}
		}

		for(i=pLoc*RSDPatternPool->patternSize;i<RSDPatternPool->patternSize*RSDPatternPool->dataSize;i++)
			RSDPatternPool->poolData[i] = 0ull;

		for(i=pLoc;i<RSDPatternPool->dataSize;i++)
		{
			RSDPatternPool->poolDataAlleleCount[i] = 0;
			RSDPatternPool->poolDataPatternCount[i] = 0;
		}

		RSDPatternPool->dataSize = pLoc; // The new patternpool corresponds to the number of patterns found in the WINDOW_SIZE snps of the chunk.
	}

}

int RSDPatternPool_pushSNP (RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, int numberOfSamples)
{
	if(RSDPatternPool->incomingSitePosition==-1.0)
		return 0;

	int i, j, lcnt=0, acnt=0;

	//init compact
	for(i=0;i<RSDPatternPool->patternSize;i++)
		RSDPatternPool->incomingSiteCompact[i] = 0ull;

	//compress snp
	i = 0;
	for(j=0;j<numberOfSamples;j++)
	{
		if(RSDPatternPool->incomingSite[j]!='0') // ambiguous
			RSDPatternPool->incomingSite[j] = '1';

		uint8_t b = RSDPatternPool->incomingSite[j]-48;
		RSDPatternPool->incomingSiteCompact[i] = (RSDPatternPool->incomingSiteCompact[i]<<1) | b;
	
		lcnt++;
		if(lcnt==64)
		{	
			acnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompact[i]);
			lcnt=0;
			i++;			
		}	
	}
	if(lcnt!=64 && lcnt!=0)
	{	
		RSDPatternPool->incomingSiteCompact[i] = RSDPatternPool->incomingSiteCompact[i] << (64 - lcnt); // get rid of stray bits
		RSDPatternPool->incomingSiteCompact[i] = RSDPatternPool->incomingSiteCompact[i] >> (64 - lcnt);

		acnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompact[i]);
	}
	
	i = 0;
	if(RSDPatternPool->dataSize==0)
	{
		memcpy(&(RSDPatternPool->poolData[i]), RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize*8);
		RSDPatternPool->poolDataAlleleCount[i] =  RSDPatternPool->incomingSiteAlleleCount;
		RSDPatternPool->poolDataPatternCount[i] = 1;
		RSDPatternPool->dataSize++;
	}
	else
	{
		int match = 0;

		for(i=0;i<RSDPatternPool->dataSize;i++) 
		{
			if(!snpv_cmp(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize))
			{	
				match=1;
				break;
			}

			if(!snpv_cmp(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize))
			{
				match=1;		
				break;
			}
		}

		if(match)
		{
			RSDPatternPool->poolDataPatternCount[i] += 1;
		}
		else
		{	
			memcpy(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize*8);
			RSDPatternPool->poolDataAlleleCount[i] =  RSDPatternPool->incomingSiteAlleleCount;
			RSDPatternPool->poolDataPatternCount[i] = 1;
			RSDPatternPool->dataSize++;
		}
	}

	RSDChunk->chunkSize++;
	
	if(RSDChunk->chunkSize>RSDChunk->chunkMemSize)
	{
		RSDChunk->chunkMemSize += CHUNK_MEMSIZE_AND_INCREMENT;
		RSDChunk->sitePosition = realloc(RSDChunk->sitePosition, sizeof(float)*RSDChunk->chunkMemSize);
		RSDChunk->derivedAlleleCount = realloc(RSDChunk->derivedAlleleCount, sizeof(int)*RSDChunk->chunkMemSize);
		RSDChunk->patternID = realloc(RSDChunk->patternID, sizeof(int)*RSDChunk->chunkMemSize);
	}

	RSDChunk->sitePosition[RSDChunk->chunkSize-1] = (float) RSDPatternPool->incomingSitePosition;
	RSDChunk->derivedAlleleCount[RSDChunk->chunkSize-1] = RSDPatternPool->incomingSiteAlleleCount;
	RSDChunk->patternID[RSDChunk->chunkSize-1] = i;

	RSDChunk->derAll1CntTotal += RSDPatternPool->incomingSiteAlleleCount==1?1:0;
	RSDChunk->derAllNCntTotal += RSDPatternPool->incomingSiteAlleleCount==numberOfSamples-1?1:0;

	// pattern pool full check
	int poolFull = 0;		

	if(RSDPatternPool->dataSize == RSDPatternPool->maxSize)
		poolFull = 1;

	return poolFull;
}
