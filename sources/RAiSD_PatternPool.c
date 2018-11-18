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

void RSDPatternPool_partialReset (RSDPatternPool_t * RSDPatternPool);

#ifdef _PTIMES
struct timespec requestStartOoC;
struct timespec requestEndOoC;
double TotalOoCTime;
#endif

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
	pp->incomingSiteCompactMask = NULL;
	pp->incomingSiteCompactWithMissing = 0;
	pp->incomingSitePosition = -1.0;
	pp->createPatternPoolMask = 0;
	pp->patternPoolMaskMode = 0;
	pp->poolData = NULL;
	pp->poolDataMask = NULL;
	pp->poolDataAlleleCount = NULL;
	pp->poolDataPatternCount = NULL;
	pp->poolDataWithMissing = NULL;
	pp->poolDataMaskCount = NULL;
	pp->poolDataAppliedMaskCount = NULL;
	pp->exchangeBuffer = NULL;	
	pp->hashMap=NULL;

	return pp;
}

void RSDPatternPool_free(RSDPatternPool_t * pp, int64_t numberOfSamples)
{
	assert(pp!=NULL);

	MemoryFootprint += (numberOfSamples+1);
	free(pp->incomingSite);

	MemoryFootprint += sizeof(uint64_t)*((unsigned long)pp->patternSize);
	free(pp->incomingSiteCompact);

	MemoryFootprint += sizeof(uint64_t)*((unsigned long)(pp->patternSize*pp->maxSize));
	free(pp->poolData);

	if(pp->createPatternPoolMask==1)
	{
		MemoryFootprint += sizeof(uint64_t)*((unsigned long)pp->patternSize);
		free(pp->incomingSiteCompactMask);

		MemoryFootprint += sizeof(uint64_t)*((unsigned long)(pp->patternSize*pp->maxSize));
		free(pp->poolDataMask);

		MemoryFootprint += sizeof(int)*((unsigned long)pp->maxSize);
		free(pp->poolDataWithMissing);

		MemoryFootprint += sizeof(int)*((unsigned long)pp->maxSize);
		free(pp->poolDataMaskCount);

		MemoryFootprint += sizeof(int)*((unsigned long)pp->maxSize);
		free(pp->poolDataAppliedMaskCount);
	}

	MemoryFootprint += sizeof(int)*((unsigned long)pp->maxSize);
	free(pp->poolDataAlleleCount);

	MemoryFootprint += sizeof(int)*((unsigned long)pp->maxSize);
	free(pp->poolDataPatternCount);

	MemoryFootprint += sizeof(uint64_t)*((unsigned long)pp->patternSize);
	free(pp->exchangeBuffer);

	if(pp->hashMap!=NULL)
		RSDHashMap_free(pp->hashMap);

	free(pp);	
}

void RSDPatternPool_init (RSDPatternPool_t * RSDPatternPool, RSDCommandLine_t * RSDCommandLine, int64_t numberOfSamples)
{
	assert(RSDCommandLine!=NULL);

	RSDPatternPool->createPatternPoolMask = (uint64_t)RSDCommandLine->createPatternPoolMask;
	RSDPatternPool->patternPoolMaskMode = (uint64_t)RSDCommandLine->patternPoolMaskMode;

	int wordLength = sizeof(uint64_t)*8;

	int wordsPerSNP = numberOfSamples%wordLength==0?(int)(numberOfSamples/wordLength):(int)(numberOfSamples/wordLength+1);
	RSDPatternPool->patternSize = wordsPerSNP;

	float maxSNPsInPool_num = ((float)RSDPatternPool->memorySize)*1024.0f*1024.0f;
	float maxSNPsInPool = maxSNPsInPool_num/(wordsPerSNP*8.0f+8.0f); // +8 for the allelecount and the patterncount
	RSDPatternPool->maxSize = (int) maxSNPsInPool;

	RSDPatternPool->incomingSite = (char*) malloc(sizeof(char)*((unsigned long)(numberOfSamples+1)));
	assert(RSDPatternPool->incomingSite!=NULL);

	RSDPatternPool->incomingSiteCompact = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)RSDPatternPool->patternSize));
	assert(RSDPatternPool->incomingSiteCompact!=NULL);

	RSDPatternPool->poolData = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)(RSDPatternPool->patternSize*RSDPatternPool->maxSize)));
	assert(RSDPatternPool->poolData!=NULL);

	if(RSDPatternPool->createPatternPoolMask==1)
	{
		RSDPatternPool->incomingSiteCompactMask = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)RSDPatternPool->patternSize));
		assert(RSDPatternPool->incomingSiteCompactMask!=NULL);

		RSDPatternPool->poolDataMask = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)(RSDPatternPool->patternSize*RSDPatternPool->maxSize)));
		assert(RSDPatternPool->poolDataMask!=NULL);

		RSDPatternPool->poolDataWithMissing = (int*) malloc(sizeof(int)*((unsigned long)RSDPatternPool->maxSize));
		assert(RSDPatternPool->poolDataWithMissing!=NULL);

		RSDPatternPool->poolDataMaskCount = (int*) malloc(sizeof(int)*((unsigned long)RSDPatternPool->maxSize));
		assert(RSDPatternPool->poolDataMaskCount!=NULL);

		RSDPatternPool->poolDataAppliedMaskCount = (int*) malloc(sizeof(int)*((unsigned long)RSDPatternPool->maxSize));
		assert(RSDPatternPool->poolDataAppliedMaskCount!=NULL);
	}

	RSDPatternPool->poolDataAlleleCount = (int*) malloc(sizeof(int)*((unsigned long)RSDPatternPool->maxSize));
	assert(RSDPatternPool->poolDataAlleleCount!=NULL);

	RSDPatternPool->poolDataPatternCount = (int*) malloc(sizeof(int)*((unsigned long)RSDPatternPool->maxSize));
	assert(RSDPatternPool->poolDataPatternCount!=NULL);

	RSDPatternPool->exchangeBuffer = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)RSDPatternPool->patternSize));
	assert(RSDPatternPool->exchangeBuffer!=NULL);

#ifdef _HM
	if(RSDPatternPool->createPatternPoolMask==0)
		RSDPatternPool->hashMap = RSDHashMap_new();
	else
		assert(0);
#endif

}

void RSDPatternPool_resize (RSDPatternPool_t * RSDPatternPool, int64_t setSamples, FILE * fpOut)
{
	if(setSamples<=(int)(((unsigned long)RSDPatternPool->patternSize)*sizeof(uint64_t)*8))
		return;
	
	int wordLength = sizeof(uint64_t)*8;

	int wordsPerSNP = setSamples%wordLength==0?(int)(setSamples/wordLength):(int)(setSamples/wordLength+1);
	RSDPatternPool->patternSize = wordsPerSNP;

	int prevMaxSize = RSDPatternPool->maxSize;
	float maxSNPsInPool_num = ((float)RSDPatternPool->memorySize)*1024.0f*1024.0f;
	float maxSNPsInPool = maxSNPsInPool_num/(wordsPerSNP*8.0f+8.0f); // +8 for the allelecount and the patterncount
	RSDPatternPool->maxSize = (int) maxSNPsInPool;
	assert(RSDPatternPool->maxSize<=prevMaxSize);

	free(RSDPatternPool->incomingSiteCompact);
	RSDPatternPool->incomingSiteCompact = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)RSDPatternPool->patternSize));
	assert(RSDPatternPool->incomingSiteCompact!=NULL);

	free(RSDPatternPool->poolData);
	RSDPatternPool->poolData = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)(RSDPatternPool->patternSize*RSDPatternPool->maxSize)));
	assert(RSDPatternPool->poolData!=NULL);

	if(RSDPatternPool->createPatternPoolMask==1)
	{
		free(RSDPatternPool->incomingSiteCompactMask);
		RSDPatternPool->incomingSiteCompactMask = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)RSDPatternPool->patternSize));
		assert(RSDPatternPool->incomingSiteCompactMask!=NULL);

		free(RSDPatternPool->poolDataMask);
		RSDPatternPool->poolDataMask = NULL;
		RSDPatternPool->poolDataMask = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)(RSDPatternPool->patternSize*RSDPatternPool->maxSize)));
		assert(RSDPatternPool->poolDataMask!=NULL);
	}

	RSDPatternPool_partialReset (RSDPatternPool);

	free(RSDPatternPool->exchangeBuffer);
	RSDPatternPool->exchangeBuffer = (uint64_t*) malloc(sizeof(uint64_t)*((unsigned long)RSDPatternPool->patternSize));
	assert(RSDPatternPool->exchangeBuffer!=NULL);

	fprintf(fpOut, " The pattern structure has been resized to %d patterns (max. capacity) and approx. %d MB memory footprint.\n", RSDPatternPool->maxSize, PATTERNPOOL_SIZE_MASK_FACTOR*RSDPatternPool->memorySize);	
	fflush(fpOut);	

	fprintf(stdout, " The pattern structure has been resized to %d patterns (max. capacity) and approx. %d MB memory footprint.\n", RSDPatternPool->maxSize, PATTERNPOOL_SIZE_MASK_FACTOR*RSDPatternPool->memorySize);	
	fflush(stdout);
}

void RSDPatternPool_print(RSDPatternPool_t * RSDPatternPool, FILE * fpOut)
{
	fprintf(fpOut, "\n A pattern structure of %d patterns (max. capacity) and approx. %d MB memory footprint has been created.\n", RSDPatternPool->maxSize, PATTERNPOOL_SIZE_MASK_FACTOR*RSDPatternPool->memorySize);	
	fflush(fpOut);
}

void RSDPatternPool_exchangePatterns (RSDPatternPool_t * RSDPatternPool, int pID_a, int pID_b)
{
	if(pID_a==pID_b)
		return;

	memcpy(RSDPatternPool->exchangeBuffer, &(RSDPatternPool->poolData[pID_a*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
	int alleleCount =  RSDPatternPool->poolDataAlleleCount[pID_a];
	int patternCount =  0; //RSDPatternPool->poolDataPatternCount[pID_a]; // Not copying the patterncount because it changes per chunk		

	memcpy(&(RSDPatternPool->poolData[pID_a*RSDPatternPool->patternSize]), &(RSDPatternPool->poolData[pID_b*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
	RSDPatternPool->poolDataAlleleCount[pID_a] = RSDPatternPool->poolDataAlleleCount[pID_b];
	RSDPatternPool->poolDataPatternCount[pID_a] = 0; //RSDPatternPool->poolDataPatternCount[pID_b];	

	memcpy(&(RSDPatternPool->poolData[pID_b*RSDPatternPool->patternSize]), RSDPatternPool->exchangeBuffer, (unsigned long)(RSDPatternPool->patternSize*8));
	RSDPatternPool->poolDataAlleleCount[pID_b] = alleleCount;
	RSDPatternPool->poolDataPatternCount[pID_b] = patternCount;	
	
	if(RSDPatternPool->createPatternPoolMask==1)
	{
		int withMissing =   RSDPatternPool->poolDataWithMissing[pID_a];
		RSDPatternPool->poolDataWithMissing[pID_a] = RSDPatternPool->poolDataWithMissing[pID_b];
		RSDPatternPool->poolDataWithMissing[pID_b] = withMissing;

		int maskCount = RSDPatternPool->poolDataMaskCount[pID_a];
		RSDPatternPool->poolDataMaskCount[pID_a] = RSDPatternPool->poolDataMaskCount[pID_b];
		RSDPatternPool->poolDataMaskCount[pID_b] = maskCount;

		int appliedMaskCount = RSDPatternPool->poolDataAppliedMaskCount[pID_a];
		RSDPatternPool->poolDataAppliedMaskCount[pID_a] = RSDPatternPool->poolDataAppliedMaskCount[pID_b];
		RSDPatternPool->poolDataAppliedMaskCount[pID_b] = appliedMaskCount;

		memcpy(RSDPatternPool->exchangeBuffer, &(RSDPatternPool->poolDataMask[pID_a*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
		memcpy(&(RSDPatternPool->poolDataMask[pID_a*RSDPatternPool->patternSize]), &(RSDPatternPool->poolDataMask[pID_b*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
		memcpy(&(RSDPatternPool->poolDataMask[pID_b*RSDPatternPool->patternSize]), RSDPatternPool->exchangeBuffer, (unsigned long)(RSDPatternPool->patternSize*8));
	}
}

void RSDPatternPool_exchangePatternsFractions (RSDPatternPool_t * RSDPatternPool, int pID_a, int pID_b)
{
	if(pID_a==pID_b)
		return;

	RSDHashMap_t * RSDHashMap = RSDPatternPool->hashMap;	

	memcpy(RSDPatternPool->exchangeBuffer, &(RSDHashMap->poolDataFractions[pID_a*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
	int alleleCount =  RSDPatternPool->poolDataAlleleCount[pID_a];
	int patternCount =  0; //RSDPatternPool->poolDataPatternCount[pID_a]; // Not copying the patterncount because it changes per chunk		

	memcpy(&(RSDHashMap->poolDataFractions[pID_a*RSDPatternPool->patternSize]), &(RSDHashMap->poolDataFractions[pID_b*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
	RSDPatternPool->poolDataAlleleCount[pID_a] = RSDPatternPool->poolDataAlleleCount[pID_b];
	RSDPatternPool->poolDataPatternCount[pID_a] = 0; //RSDPatternPool->poolDataPatternCount[pID_b];	

	memcpy(&(RSDHashMap->poolDataFractions[pID_b*RSDPatternPool->patternSize]), RSDPatternPool->exchangeBuffer, (unsigned long)(RSDPatternPool->patternSize*8));
	RSDPatternPool->poolDataAlleleCount[pID_b] = alleleCount;
	RSDPatternPool->poolDataPatternCount[pID_b] = patternCount;	
	
	if(RSDPatternPool->createPatternPoolMask==1)
	{
		int withMissing =   RSDPatternPool->poolDataWithMissing[pID_a];
		RSDPatternPool->poolDataWithMissing[pID_a] = RSDPatternPool->poolDataWithMissing[pID_b];
		RSDPatternPool->poolDataWithMissing[pID_b] = withMissing;

		int maskCount = RSDPatternPool->poolDataMaskCount[pID_a];
		RSDPatternPool->poolDataMaskCount[pID_a] = RSDPatternPool->poolDataMaskCount[pID_b];
		RSDPatternPool->poolDataMaskCount[pID_b] = maskCount;

		int appliedMaskCount = RSDPatternPool->poolDataAppliedMaskCount[pID_a];
		RSDPatternPool->poolDataAppliedMaskCount[pID_a] = RSDPatternPool->poolDataAppliedMaskCount[pID_b];
		RSDPatternPool->poolDataAppliedMaskCount[pID_b] = appliedMaskCount;

		memcpy(RSDPatternPool->exchangeBuffer, &(RSDPatternPool->poolDataMask[pID_a*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
		memcpy(&(RSDPatternPool->poolDataMask[pID_a*RSDPatternPool->patternSize]), &(RSDPatternPool->poolDataMask[pID_b*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
		memcpy(&(RSDPatternPool->poolDataMask[pID_b*RSDPatternPool->patternSize]), RSDPatternPool->exchangeBuffer, (unsigned long)(RSDPatternPool->patternSize*8));
	}
}

void RSDPatternPool_partialReset (RSDPatternPool_t * RSDPatternPool)
{
	int i;
	for(i=0;i<RSDPatternPool->patternSize;i++)
	{
		RSDPatternPool->incomingSiteCompact[i] = 0ull;

		if(RSDPatternPool->createPatternPoolMask==1)
			RSDPatternPool->incomingSiteCompactMask[i] = 0ull;
	}

	for(i=0;i<RSDPatternPool->patternSize*RSDPatternPool->dataSize;i++)
	{
		RSDPatternPool->poolData[i] = 0ull;

		if(RSDPatternPool->createPatternPoolMask==1)
			RSDPatternPool->poolDataMask[i] = 0ull;
	}	
}

void RSDPatternPool_reset (RSDPatternPool_t * RSDPatternPool, int64_t numberOfSamples, int64_t setSamples, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);

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

		RSDPatternPool->incomingSiteDerivedAlleleCount = 0;

		for(i=0;i<RSDPatternPool->patternSize;i++)
		{
			RSDPatternPool->incomingSiteCompact[i] = 0ull;

			if(RSDPatternPool->createPatternPoolMask==1)
				RSDPatternPool->incomingSiteCompactMask[i] = 0ull;
		}
		RSDPatternPool->incomingSitePosition = -1.0;

#ifdef _HM
		for(i=0;i<RSDPatternPool->patternSize*RSDPatternPool->maxSize;i++)
		{
			RSDPatternPool->poolData[i] = 0ull;

			if(RSDPatternPool->createPatternPoolMask==1)
				RSDPatternPool->poolDataMask[i] = 0ull;
		}

		for(i=0;i<RSDPatternPool->maxSize;i++)
		{
			RSDPatternPool->poolDataAlleleCount[i] = 0;
			RSDPatternPool->poolDataPatternCount[i] = 0;

			if(RSDPatternPool->createPatternPoolMask==1)
			{
				RSDPatternPool->poolDataWithMissing[i] = 0;
				RSDPatternPool->poolDataMaskCount[i] = 0;
				RSDPatternPool->poolDataAppliedMaskCount[i] = 0;
			}
		}

		for(i=0;i<RSDPatternPool->hashMap->addressListSize;i++)
			RSDPatternPool->hashMap->addressListEntrySize[i] = 0;

		RSDPatternPool->hashMap->addressListEntryFull = 0;
#else

		for(i=0;i<RSDPatternPool->patternSize*RSDPatternPool->dataSize;i++)
		{
			RSDPatternPool->poolData[i] = 0ull;

			if(RSDPatternPool->createPatternPoolMask==1)
				RSDPatternPool->poolDataMask[i] = 0ull;
		}

		for(i=0;i<RSDPatternPool->dataSize;i++)
		{
			RSDPatternPool->poolDataAlleleCount[i] = 0;
			RSDPatternPool->poolDataPatternCount[i] = 0;

			if(RSDPatternPool->createPatternPoolMask==1)
			{
				RSDPatternPool->poolDataWithMissing[i] = 0;
				RSDPatternPool->poolDataMaskCount[i] = 0;
				RSDPatternPool->poolDataAppliedMaskCount[i] = 0;
			}
		}
		
#endif
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

		RSDPatternPool->incomingSiteDerivedAlleleCount = 0;

		for(i=0;i<RSDPatternPool->patternSize;i++)
		{
			RSDPatternPool->incomingSiteCompact[i] = 0ull;

			if(RSDPatternPool->createPatternPoolMask==1)
				RSDPatternPool->incomingSiteCompactMask[i] = 0ull;

		}
		RSDPatternPool->incomingSitePosition = -1.0;


		// End part of chunk relocated to the beginning
		int j = 0;
		int pLoc = 0;

#ifdef _MLT
#ifdef _HM  // HM and MLT
		RSDHashMap_t * RSDHashMap = RSDPatternPool->hashMap;

		for(i=1;i<RSDHashMap->addressListSize;i++)
			RSDHashMap->addressListEntrySize[i] = 0;

		for(i=1;i<RSDHashMap->addressListSize;i++)
		{
			int base = (i-1)*(int)RSDHashMap->addressListEntryMaxSize;

			for(j=0;j<RSDChunk->chunkSize;j++)
			{
				int pID = (int)RSDChunk->chunkData[j*3+2];
				
				if(pID>=base && pID<base+(int)RSDHashMap->addressListEntryMaxSize) 
				{	
					int pID_a = pID;
					int pID_b = base+(int)RSDHashMap->addressListEntrySize[i];
	
					if(pID_a >= pID_b)
					{
						RSDPatternPool_exchangePatternsFractions (RSDPatternPool, pID_a, pID_b);
				
						int k;
						for(k=j;k<RSDChunk->chunkSize;k++)
						{
							if((int)RSDChunk->chunkData[k*3+2]==pID_a)
							{
								RSDChunk->chunkData[k*3+2] = (float)pID_b;
								RSDPatternPool->poolDataPatternCount[pID_b]++;
							}
							else
							{
								if((int)RSDChunk->chunkData[k*3+2]==pID_b)
									RSDChunk->chunkData[k*3+2] = (float)pID_a;
							}
						}

						RSDHashMap->addressListEntrySize[i]++;
					}
				}
			}
		}

		RSDPatternPool->hashMap->addressListEntryFull = 0;
		
#else  // MLT, no HM, only for runs without mask
		if(RSDCommandLine->createPatternPoolMask==1)
		{
			for(i=0;i<RSDChunk->chunkSize;i++)
			{
				int pID_a = RSDChunk->patternID[i];
				int pID_b = pLoc;
		
				if(pID_a >= pID_b)
				{
					RSDPatternPool_exchangePatterns (RSDPatternPool, pID_a, pID_b);

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
		}
		else
		{
			for(i=0;i<RSDChunk->chunkSize;i++)
			{
				int pID_a = (int)RSDChunk->chunkData[i*3+2];
				int pID_b = pLoc;
		
				if(pID_a >= pID_b)
				{
					RSDPatternPool_exchangePatterns (RSDPatternPool, pID_a, pID_b);

					for(j=i;j<RSDChunk->chunkSize;j++)
					{
						if(((int)RSDChunk->chunkData[j*3+2])==pID_a)
						{
							RSDChunk->chunkData[j*3+2] = (float)pID_b;
							RSDPatternPool->poolDataPatternCount[pLoc]++;
						}
						else
						{
							if(((int)RSDChunk->chunkData[j*3+2])==pID_b)
								RSDChunk->chunkData[j*3+2] = (float)pID_a;
						}
					}

					pLoc++;
				}
			}			
		}
#endif
#else
#ifdef _HM // HM, no MLT
		RSDHashMap_t * RSDHashMap = RSDPatternPool->hashMap;

		for(i=1;i<RSDHashMap->addressListSize;i++)
			RSDHashMap->addressListEntrySize[i] = 0;

		for(i=1;i<RSDHashMap->addressListSize;i++)
		{
			int base = (i-1)*(int)RSDPatternPool->hashMap->addressListEntryMaxSize;

			for(j=0;j<RSDChunk->chunkSize;j++)
			{
				int pID = RSDChunk->patternID[j];
				
				if(pID>=base && pID<base+(int)RSDHashMap->addressListEntryMaxSize) 
				{	
					int pID_a = pID;
					int pID_b = base+(int)RSDHashMap->addressListEntrySize[i];
	
					if(pID_a >= pID_b)
					{
						RSDPatternPool_exchangePatternsFractions (RSDPatternPool, pID_a, pID_b);
				
						int k;
						for(k=j;k<RSDChunk->chunkSize;k++)
						{
							if(RSDChunk->patternID[k]==pID_a)
							{
								RSDChunk->patternID[k] = pID_b;
								RSDPatternPool->poolDataPatternCount[pID_b]++;
							}
							else
							{
								if(RSDChunk->patternID[k]==pID_b)
									RSDChunk->patternID[k] = pID_a;
							}
						}

						RSDHashMap->addressListEntrySize[i]++;
					}
				}
			}
		}

		RSDPatternPool->hashMap->addressListEntryFull = 0;

#else // Default: no HM, no MLT

		for(i=0;i<RSDChunk->chunkSize;i++)
		{
			int pID_a = RSDChunk->patternID[i];
			int pID_b = pLoc;
	
			if(pID_a >= pID_b)
			{
				RSDPatternPool_exchangePatterns (RSDPatternPool, pID_a, pID_b);

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
#endif
#endif

		for(i=pLoc*RSDPatternPool->patternSize;i<RSDPatternPool->patternSize*RSDPatternPool->dataSize;i++)
		{
			RSDPatternPool->poolData[i] = 0ull;
			
			if(RSDPatternPool->createPatternPoolMask==1)
				RSDPatternPool->poolDataMask[i] = 0ull;
		}

		for(i=pLoc;i<RSDPatternPool->dataSize;i++)
		{
			RSDPatternPool->poolDataAlleleCount[i] = 0;
			RSDPatternPool->poolDataPatternCount[i] = 0;

			if(RSDPatternPool->createPatternPoolMask==1)
			{
				RSDPatternPool->poolDataWithMissing[i] = 0;
				RSDPatternPool->poolDataMaskCount[i] = 0;
				RSDPatternPool->poolDataAppliedMaskCount[i] = 0;
			}
		}

		RSDPatternPool->dataSize = pLoc; // The new pattern pool corresponds to the number of patterns found in the WINDOW_SIZE snps of the chunk.
	}
}

int RSDPatternPool_pushSNP (RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, int64_t numberOfSamples, RSDCommandLine_t * RSDCommandLine)
{
#ifdef _PTIMES
	clock_gettime(CLOCK_REALTIME, &requestStartOoC);
#endif

	if(RSDPatternPool->incomingSitePosition<=-1.0)
		return 0;

	assert(RSDCommandLine!=NULL);

	int i, j, lcnt=0, acnt=0, vcnt=0, avcnt=0, match=0;

	//init compact
	for(i=0;i<RSDPatternPool->patternSize;i++)
	{
		RSDPatternPool->incomingSiteCompact[i] = 0ull;

		if(RSDPatternPool->createPatternPoolMask==1)
			RSDPatternPool->incomingSiteCompactMask[i] = 0ull;
	}

	//compress snp
	i = 0;

	// with valid mask (-M)
	if(RSDPatternPool->createPatternPoolMask==1)
	{
		RSDPatternPool->incomingSiteCompactWithMissing = 0;
		for(j=0;j<numberOfSamples;j++)
		{
			switch(RSDPatternPool->incomingSite[j])
			{
				case 'N':
					assert(RSDPatternPool->incomingSite[j]=='N');

					RSDPatternPool->incomingSite[j] = '0';

					uint8_t bN = (uint8_t)(RSDPatternPool->incomingSite[j]-48);
					RSDPatternPool->incomingSiteCompact[i] = (RSDPatternPool->incomingSiteCompact[i]<<1) | bN;

					uint8_t bNValid = 0;
					RSDPatternPool->incomingSiteCompactMask[i] = (RSDPatternPool->incomingSiteCompactMask[i]<<1) | bNValid;

					RSDPatternPool->incomingSiteCompactWithMissing = 1;

					break;

				default:
					assert(RSDPatternPool->incomingSite[j]=='0' || RSDPatternPool->incomingSite[j]=='1');

					uint8_t b = (uint8_t)(RSDPatternPool->incomingSite[j]-48);
					RSDPatternPool->incomingSiteCompact[i] = (RSDPatternPool->incomingSiteCompact[i]<<1) | b;

					uint8_t bValid = 1;
					RSDPatternPool->incomingSiteCompactMask[i] = (RSDPatternPool->incomingSiteCompactMask[i]<<1) | bValid;

					break;
			}

			lcnt++;
			if(lcnt==64)
			{	
				acnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompact[i]);
				vcnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompactMask[i]);
				avcnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompact[i]&RSDPatternPool->incomingSiteCompactMask[i]);
				lcnt=0;
				i++;			
			}	
		}
	}
	else
	{	// default or with N imputation (-M 1)
		for(j=0;j<numberOfSamples;j++)
		{
			assert(RSDPatternPool->incomingSite[j]=='0' || RSDPatternPool->incomingSite[j]=='1');

			uint8_t b = (uint8_t)(RSDPatternPool->incomingSite[j]-48);
			RSDPatternPool->incomingSiteCompact[i] = (RSDPatternPool->incomingSiteCompact[i]<<1) | b;
	
			lcnt++;
			if(lcnt==64)
			{	
				acnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompact[i]);
				lcnt=0;
				i++;			
			}	
		}
	}

	if(lcnt!=64 && lcnt!=0)
	{	
		RSDPatternPool->incomingSiteCompact[i] = RSDPatternPool->incomingSiteCompact[i] << (64 - lcnt); // get rid of stray bits
		RSDPatternPool->incomingSiteCompact[i] = RSDPatternPool->incomingSiteCompact[i] >> (64 - lcnt);

		acnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompact[i]);

		if(RSDPatternPool->createPatternPoolMask==1)
		{
			RSDPatternPool->incomingSiteCompactMask[i] = RSDPatternPool->incomingSiteCompactMask[i] << (64 - lcnt); // get rid of stray bits
			RSDPatternPool->incomingSiteCompactMask[i] = RSDPatternPool->incomingSiteCompactMask[i] >> (64 - lcnt);

			vcnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompactMask[i]);
			avcnt += rsd_popcnt_u64(RSDPatternPool->incomingSiteCompact[i]&RSDPatternPool->incomingSiteCompactMask[i]);		
		}
	}

	assert(i==RSDPatternPool->patternSize-1);


#ifdef _HM
	RSDHashMap_t * RSDHashMap = RSDPatternPool->hashMap;
	assert(RSDHashMap!=NULL);

	match = 0;
	i = 0;

	RSDHashMap_setMainKey (RSDHashMap, (int64_t)acnt);
	RSDHashMap_setSecondaryKey (RSDHashMap, numberOfSamples-((int64_t)acnt));
	i = RSDHashMap_scanPatternPoolFractions (RSDHashMap, RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize, (int)numberOfSamples, &match);

	if(match)
	{
		RSDPatternPool->poolDataPatternCount[i] += 1;
	}
	else
	{	
		memcpy(&(RSDHashMap->poolDataFractions[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, (unsigned long)(RSDPatternPool->patternSize*8));
		RSDPatternPool->poolDataAlleleCount[i] =  RSDPatternPool->incomingSiteDerivedAlleleCount;
		RSDPatternPool->poolDataPatternCount[i] = 1;
		RSDPatternPool->dataSize++;	
	}
#else
	match = 0;
	i = 0;

	if(RSDPatternPool->dataSize==0)
	{
		memcpy(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, (unsigned long)(RSDPatternPool->patternSize*8));

		if(RSDPatternPool->createPatternPoolMask==1)
		{
			memcpy(&(RSDPatternPool->poolDataMask[i]), RSDPatternPool->incomingSiteCompactMask, (unsigned long)(RSDPatternPool->patternSize*8));
			RSDPatternPool->poolDataWithMissing[i] = (int)RSDPatternPool->incomingSiteCompactWithMissing;

			if(RSDPatternPool->poolDataWithMissing[i]==1)
			{
				assert(vcnt<numberOfSamples);
				assert(avcnt<vcnt);
			}
			else
			{
				assert(vcnt==numberOfSamples);
				assert(avcnt==RSDPatternPool->incomingSiteDerivedAlleleCount);
			}
			assert(acnt<vcnt);

			RSDPatternPool->poolDataMaskCount[i] = vcnt;
			RSDPatternPool->poolDataAppliedMaskCount[i] = avcnt;
		}
		RSDPatternPool->poolDataAlleleCount[i] =  RSDPatternPool->incomingSiteDerivedAlleleCount;
		RSDPatternPool->poolDataPatternCount[i] = 1;
		RSDPatternPool->dataSize++;
	}
	else
	{
		match = 0;
		if(RSDPatternPool->createPatternPoolMask==1)
		{
			if(RSDPatternPool->patternPoolMaskMode==0)
			{
				for(i=0;i<RSDPatternPool->dataSize;i++) 
				{
					if(!snpv_cmp(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize))
					{	
						if(!snpv_cmp(&(RSDPatternPool->poolDataMask[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompactMask, RSDPatternPool->patternSize))
						{	
							match=1;
							break;
						}
					}
				}

				if(match==0)
				for(i=0;i<RSDPatternPool->dataSize;i++) 
				{
					if(!snpv_cmp(&(RSDPatternPool->poolDataMask[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompactMask, RSDPatternPool->patternSize))
					{	
						
						if(!isnpv_cmp_with_mask(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, &(RSDPatternPool->poolDataMask[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompactMask, RSDPatternPool->patternSize, (int)numberOfSamples))
						{	
							match=1;
							break;
						}
					}
				}
			}
			else
			{
				for(i=0;i<RSDPatternPool->dataSize;i++) 
				{

					if(!snpv_cmp_cross_masks(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), 
								   RSDPatternPool->incomingSiteCompact, 
								 &(RSDPatternPool->poolDataMask[i*RSDPatternPool->patternSize]), 
								   RSDPatternPool->incomingSiteCompactMask, 
								   RSDPatternPool->patternSize))
					{	
						match=1;
						break;
					}
				}

				if(match==0)
				for(i=0;i<RSDPatternPool->dataSize;i++) 
				{

					if(!isnpv_cmp_cross_masks(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), 
								   RSDPatternPool->incomingSiteCompact, 
								 &(RSDPatternPool->poolDataMask[i*RSDPatternPool->patternSize]), 
								   RSDPatternPool->incomingSiteCompactMask, 
								   RSDPatternPool->patternSize))
					{	
						match=1;
						break;
					}
				}
			}
		}
		else
		{
			for(i=0;i<RSDPatternPool->dataSize;i++) 
			{
				if(!snpv_cmp(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize))
				{	
					match=1;
					break;
				}
			}
			if(match==0)
			for(i=0;i<RSDPatternPool->dataSize;i++) 
			{
				if(!isnpv_cmp(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, RSDPatternPool->patternSize, (int)numberOfSamples))
				{	
					match=1;
					break;
				}	
			}	

		}

		if(match)
		{
			RSDPatternPool->poolDataPatternCount[i] += 1;
		}
		else
		{	
			memcpy(&(RSDPatternPool->poolData[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompact, (unsigned long)(RSDPatternPool->patternSize*8));

			if(RSDPatternPool->createPatternPoolMask==1)
			{
				memcpy(&(RSDPatternPool->poolDataMask[i*RSDPatternPool->patternSize]), RSDPatternPool->incomingSiteCompactMask, (unsigned long)(RSDPatternPool->patternSize*8));
				RSDPatternPool->poolDataWithMissing[i] = (int)RSDPatternPool->incomingSiteCompactWithMissing;

				if(RSDPatternPool->poolDataWithMissing[i]==1)
				{
					assert(vcnt<numberOfSamples);
					assert(avcnt<vcnt);
				}
				else
				{
					assert(vcnt==numberOfSamples);
					assert(avcnt==RSDPatternPool->incomingSiteDerivedAlleleCount);
				}
				assert(acnt<vcnt);				

				RSDPatternPool->poolDataMaskCount[i] = vcnt;
				RSDPatternPool->poolDataAppliedMaskCount[i] = avcnt;	
			}

			RSDPatternPool->poolDataAlleleCount[i] =  RSDPatternPool->incomingSiteDerivedAlleleCount;
			RSDPatternPool->poolDataPatternCount[i] = 1;
			RSDPatternPool->dataSize++;
		}
	}
#endif
	RSDChunk->chunkSize++;
	
	if(RSDChunk->chunkSize>RSDChunk->chunkMemSize)
	{
		RSDChunk->chunkMemSize += CHUNK_MEMSIZE_AND_INCREMENT;

#ifdef _MLT
		if(RSDCommandLine->createPatternPoolMask==1)
		{
			RSDChunk->sitePosition = realloc(RSDChunk->sitePosition, sizeof(float)*((unsigned long)RSDChunk->chunkMemSize));
			RSDChunk->derivedAlleleCount = realloc(RSDChunk->derivedAlleleCount, sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
			RSDChunk->patternID = realloc(RSDChunk->patternID, sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
		}
		else
		{
			RSDChunk->chunkData = realloc(RSDChunk->chunkData, sizeof(float)*((unsigned long)RSDChunk->chunkMemSize*3));
			assert(RSDChunk->chunkData!=NULL);
		}
#else
		RSDChunk->sitePosition = realloc(RSDChunk->sitePosition, sizeof(float)*((unsigned long)RSDChunk->chunkMemSize));
		RSDChunk->derivedAlleleCount = realloc(RSDChunk->derivedAlleleCount, sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
		RSDChunk->patternID = realloc(RSDChunk->patternID, sizeof(int)*((unsigned long)RSDChunk->chunkMemSize));
#endif
	}

#ifdef _MLT
	if(RSDCommandLine->createPatternPoolMask==1)
	{
		RSDChunk->sitePosition[RSDChunk->chunkSize-1] = (float) RSDPatternPool->incomingSitePosition;
		RSDChunk->derivedAlleleCount[RSDChunk->chunkSize-1] = RSDPatternPool->incomingSiteDerivedAlleleCount;
		RSDChunk->patternID[RSDChunk->chunkSize-1] = i;
	}
	else
	{
		RSDChunk->chunkData[(RSDChunk->chunkSize-1)*3+0] =  (float) RSDPatternPool->incomingSitePosition;
		RSDChunk->chunkData[(RSDChunk->chunkSize-1)*3+1] =  (float) RSDPatternPool->incomingSiteDerivedAlleleCount;
		RSDChunk->chunkData[(RSDChunk->chunkSize-1)*3+2] =  (float) i;
	}
#else
	RSDChunk->sitePosition[RSDChunk->chunkSize-1] = (float) RSDPatternPool->incomingSitePosition;
	RSDChunk->derivedAlleleCount[RSDChunk->chunkSize-1] = RSDPatternPool->incomingSiteDerivedAlleleCount;
	RSDChunk->patternID[RSDChunk->chunkSize-1] = i;
#endif
	RSDChunk->derAll1CntTotal += RSDPatternPool->incomingSiteDerivedAlleleCount==1?1:0;
	RSDChunk->derAllNCntTotal += RSDPatternPool->incomingSiteDerivedAlleleCount==numberOfSamples-1?1:0;

	// pattern pool full check
	int poolFull = 0;

#ifdef _HM
	if(RSDHashMap->addressListEntryFull==1)
		poolFull = 1;
#else
	if(RSDPatternPool->dataSize == RSDPatternPool->maxSize)
		poolFull = 1;
#endif	

#ifdef _PTIMES
	clock_gettime(CLOCK_REALTIME, &requestEndOoC);
	double OoCTime = (requestEndOoC.tv_sec-requestStartOoC.tv_sec)+ (requestEndOoC.tv_nsec-requestStartOoC.tv_nsec) / BILLION;
	TotalOoCTime += OoCTime;
#endif	

	return poolFull;
}

int RSDPatternPool_imputeIncomingSite (RSDPatternPool_t * RSDPatternPool, int64_t setSamples)
{
	double prob1 = ((double)RSDPatternPool->incomingSiteDerivedAlleleCount)/((double)RSDPatternPool->incomingSiteTotalAlleleCount);
	int64_t i;
	int flag = 0;
	for(i=0;i<setSamples;i++)
	{
		if(RSDPatternPool->incomingSite[i]=='N')
		{
			flag = 1;

			int val = rand();
			double rval = ((double)val) / ((double)RAND_MAX);
			if(rval<=prob1)
			{
				RSDPatternPool->incomingSite[i] = '1';
				RSDPatternPool->incomingSiteDerivedAlleleCount++;
				RSDPatternPool->incomingSiteTotalAlleleCount++;	
			}
			else
			{
				RSDPatternPool->incomingSite[i] = '0';
				RSDPatternPool->incomingSiteTotalAlleleCount++;
			}
		}		
	}
	return flag;
}

void RSDPatternPool_assessMissing  (RSDPatternPool_t * RSDPatternPool, int64_t numberOfSamples)
{
	if(RSDPatternPool->createPatternPoolMask==0)
		return;

	int i;
	for(i=0;i<RSDPatternPool->dataSize;i++)
	{
		int j;
		int acnt = 0;
		int vcnt = 0;
		for(j=0;j<RSDPatternPool->patternSize;j++)
		{
			acnt += rsd_popcnt_u64(RSDPatternPool->poolData[i*RSDPatternPool->patternSize+j]);
			vcnt += rsd_popcnt_u64(RSDPatternPool->poolDataMask[i*RSDPatternPool->patternSize+j]);
		}

		if(RSDPatternPool->poolDataWithMissing[i]==1)
			assert(vcnt<numberOfSamples);
		else
			assert(vcnt==numberOfSamples);

		assert(acnt<vcnt);

		assert(RSDPatternPool->poolDataMaskCount[i]==vcnt);

		if(RSDPatternPool->poolDataWithMissing[i]==1)
			assert(RSDPatternPool->poolDataMaskCount[i]-RSDPatternPool->poolDataAppliedMaskCount[i]<numberOfSamples-RSDPatternPool->poolDataAlleleCount[i]);
		else
			assert(RSDPatternPool->poolDataMaskCount[i]-RSDPatternPool->poolDataAppliedMaskCount[i]==numberOfSamples-RSDPatternPool->poolDataAlleleCount[i]);	
		
		assert(RSDPatternPool->poolDataAppliedMaskCount[i]==RSDPatternPool->poolDataAlleleCount[i]);
	}	
}
