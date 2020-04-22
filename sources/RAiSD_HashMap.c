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

RSDHashMap_t * 	RSDHashMap_new 	(void)
{
	RSDHashMap_t * hm = NULL;
	hm = (RSDHashMap_t *)rsd_malloc(sizeof(RSDHashMap_t));
	assert(hm!=NULL);

	hm->addressListSize = 0;
	hm->addressList = NULL;
	hm->addressListEntryMaxSize = 0;
	hm->addressListEntrySize = NULL;
	hm->addressListEntryFull = 0;

	hm->mainKey = 0;
	hm->mainAddress = NULL;
	hm->secondaryKey = 0;
	hm->secondaryAddress = NULL;

	hm->poolDataFractions = NULL;

	return hm;	
}

void RSDHashMap_free (RSDHashMap_t * hm)
{
	assert(hm!=NULL);
	
	if(hm->addressList!=NULL)
		free(hm->addressList);

	hm->addressListSize = 0;
	hm->addressList = NULL;

	if(hm->addressListEntrySize!=NULL)
		free(hm->addressListEntrySize);

	hm->addressListEntryMaxSize = 0;
	hm->addressListEntrySize = NULL;
	hm->addressListEntryFull = 0;

	if(hm->poolDataFractions!=NULL)
		free(hm->poolDataFractions);

	hm->poolDataFractions = NULL;

	hm->mainKey = 0;
	hm->secondaryKey = 0;
	hm->mainAddress = NULL;
	hm->secondaryAddress = NULL;
		
	free(hm);
}

void RSDHashMap_init (RSDHashMap_t * RSDHashMap, int64_t setSamples, int maxSize, int patternSize)
{
	assert(RSDHashMap!=NULL);

	RSDHashMap->addressListSize = setSamples;
	RSDHashMap->addressList = (uint64_t**)rsd_malloc(sizeof(uint64_t*)*(unsigned long)RSDHashMap->addressListSize);
	assert(RSDHashMap->addressList!=NULL);

	RSDHashMap->addressListEntrySize = (uint64_t*)rsd_malloc(sizeof(uint64_t)*(unsigned long)RSDHashMap->addressListSize);
	assert(RSDHashMap->addressListEntrySize!=NULL);

	RSDHashMap->addressList[0]=NULL;
	RSDHashMap->addressListEntrySize[0] = 0;

	int curMaxSize = maxSize - maxSize%(setSamples-1);
	RSDHashMap->addressListEntryMaxSize = (unsigned long)curMaxSize / (unsigned long)(setSamples-1);
	assert(RSDHashMap->addressListEntryMaxSize>=1);

	RSDHashMap->poolDataFractions = (uint64_t*) rsd_malloc(sizeof(uint64_t)*((unsigned long)(patternSize*maxSize)));
	assert(RSDHashMap->poolDataFractions!=NULL);
	
	int i=0;
	for(i=1;i<RSDHashMap->addressListSize;i++)
	{
		RSDHashMap->addressList[i] = RSDHashMap->poolDataFractions+(unsigned long)(i-1)*RSDHashMap->addressListEntryMaxSize*(unsigned long)patternSize;
		RSDHashMap->addressListEntrySize[i] = 0;
	}	
}

void RSDHashMap_setMainKey (RSDHashMap_t * RSDHashMap, int64_t mainKey)
{
	assert(RSDHashMap!=NULL);
	assert(mainKey>0 && mainKey<RSDHashMap->addressListSize);

	RSDHashMap->mainKey =  mainKey; 
	RSDHashMap->mainAddress = RSDHashMap->addressList[mainKey];
}

void RSDHashMap_setSecondaryKey (RSDHashMap_t * RSDHashMap, int64_t secondaryKey)
{
	assert(RSDHashMap!=NULL);
	assert(secondaryKey>0 && secondaryKey<RSDHashMap->addressListSize);

	assert(secondaryKey==RSDHashMap->addressListSize-RSDHashMap->mainKey);

	RSDHashMap->secondaryKey =  secondaryKey; 
	RSDHashMap->secondaryAddress = RSDHashMap->addressList[secondaryKey];
}

int RSDHashMap_scanPatternPoolFractions (RSDHashMap_t * RSDHashMap, uint64_t * incomingSiteCompact, int patternSize, int numberOfSamples, int * match)
{
	assert(RSDHashMap!=NULL);
	assert(incomingSiteCompact!=NULL);
	assert(patternSize>=1);

	uint64_t * poolData = NULL;
	uint64_t key=0, size=0;
	uint64_t i = 0;
	int isMatch = 0;

	// mainKey fraction scan
	poolData = RSDHashMap->mainAddress;
	assert(poolData!=NULL);

	key = (uint64_t)RSDHashMap->mainKey;
	size = RSDHashMap->addressListEntrySize[key];

	if(size>0)
	{
		for(i=0;i<size;i++)
		{
			if(!snpv_cmp(&(poolData[i*(unsigned long)patternSize]), incomingSiteCompact, patternSize))
			{	
				isMatch=1;
				break;
			}
		}
		if(isMatch==1)
		{
			*match = isMatch;
			return (int)((key-1)*RSDHashMap->addressListEntryMaxSize+i);
		}
	}

	// secondaryKey fraction scan
	poolData = RSDHashMap->secondaryAddress;
	assert(poolData!=NULL);

	key = (uint64_t)RSDHashMap->secondaryKey;
	size = RSDHashMap->addressListEntrySize[key];

	if(size>0)
	{
		for(i=0;i<size;i++)
		{
			if(!isnpv_cmp(&(poolData[i*(unsigned long)patternSize]), incomingSiteCompact, patternSize, numberOfSamples))
			{	
				isMatch=1;
				break;
			}
		}
		if(isMatch==1)
		{
			*match = isMatch;
			return (int)((key-1)*RSDHashMap->addressListEntryMaxSize+i);
		}
	}

	// return no match for mainKey
	assert(isMatch==0);
	*match = isMatch;

	key = (uint64_t)RSDHashMap->mainKey;
	RSDHashMap->addressListEntrySize[key]++;
	
	RSDHashMap->addressListEntryFull = RSDHashMap->addressListEntrySize[key]==RSDHashMap->addressListEntryMaxSize?1:0;

	return (int)((key-1)*RSDHashMap->addressListEntryMaxSize+(RSDHashMap->addressListEntrySize[key]-1));
}
	

