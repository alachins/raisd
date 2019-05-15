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

RSDLutMap_t * RSDLutMap_new (void)
{
	RSDLutMap_t * lm = (RSDLutMap_t *)rsd_malloc(sizeof(RSDLutMap_t));
	assert(lm!=NULL);
	
	lm->lutMapSize = 0;

	lm->lutMapMem = NULL;
	lm->lutMap = NULL;

	lm->lutMapMemC = NULL;
	lm->lutMapC = NULL;

	return lm;
}

void RSDLutMap_free (RSDLutMap_t * lm)
{
	if(lm->lutMapMem!=NULL)
	{
		free(lm->lutMapMem);
		lm->lutMapMem = NULL;
	}

	if(lm->lutMap!=NULL)
	{
		free(lm->lutMap);
		lm->lutMap = NULL;
	}

	if(lm->lutMapMemC!=NULL)
	{
		free(lm->lutMapMemC);
		lm->lutMapMemC = NULL;
	}

	if(lm->lutMapC!=NULL)
	{
		free(lm->lutMapC);
		lm->lutMapC = NULL;
	}

	free(lm);
}

void RSDLutMap_reset (RSDLutMap_t * lm)
{
	int i=0, j=0;
	for(i=0;i<lm->lutMapSize;i++)
		for(j=0;j<LUTMAP_WIDTH;j++)
		{	
			lm->lutMap[i][j]=0;
			lm->lutMapC[i][j]=0;
		}
}

void RSDLutMap_init (RSDLutMap_t * lm, int64_t numberOfSamples)
{
	assert(lm!=NULL);

	lm->lutMapSize =  numberOfSamples%LUTMAP_GROUPSIZE==0?numberOfSamples/LUTMAP_GROUPSIZE:numberOfSamples/LUTMAP_GROUPSIZE+1;

	lm->lutMapMem = (uint8_t*)calloc((unsigned long)lm->lutMapSize*LUTMAP_WIDTH, sizeof(uint8_t));
	assert(lm->lutMapMem!=NULL);

	lm->lutMap = (uint8_t**)rsd_malloc(sizeof(uint8_t*)*(unsigned long)lm->lutMapSize);
	assert(lm->lutMap);

	int i=0;
	for(i=0;i<lm->lutMapSize;i++)
		lm->lutMap[i] = &(lm->lutMapMem[i*LUTMAP_WIDTH]);

	lm->lutMapMemC = (uint8_t*)calloc((unsigned long)lm->lutMapSize*LUTMAP_WIDTH, sizeof(uint8_t));
	assert(lm->lutMapMemC!=NULL);

	lm->lutMapC = (uint8_t**)rsd_malloc(sizeof(uint8_t*)*(unsigned long)lm->lutMapSize);
	assert(lm->lutMapC);

	for(i=0;i<lm->lutMapSize;i++)
		lm->lutMapC[i] = &(lm->lutMapMemC[i*LUTMAP_WIDTH]);
	
}

int RSDLutMap_scan (RSDLutMap_t * lm, uint64_t * query) // This returns 1 when pattern definitely new.
{
	assert(lm!=NULL);
	assert(query!=NULL);

	int i=0, y=0;

	uint64_t tquery = query[0];
	for(i=0;i<lm->lutMapSize;i++)
	{
		tquery = i%LUTMAP_INTERVAL==0?query[y++]:tquery;

		uint8_t lutMapAddr = tquery & 0xff;

		if(lm->lutMap[i][lutMapAddr]==0)
			return 1;

		tquery >>= LUTMAP_GROUPSIZE;		
	}

	return 0;
}

int RSDLutMap_scanC (RSDLutMap_t * lm, uint64_t * query, int64_t patternSize, int64_t numberOfSamples) // This returns 1 when pattern definitely new.
{
	assert(lm!=NULL);
	assert(query!=NULL);

	int i=0, y=0;

	int fGroups = ((int)patternSize-1)*LUTMAP_INTERVAL;

	uint64_t tquery = ~query[0];
	for(i=0;i<fGroups;i++)
	{
		tquery = i%LUTMAP_INTERVAL==0?~query[y++]:tquery;

		uint8_t lutMapAddr = tquery & 0xff;

		if(lm->lutMap[i][lutMapAddr]==0)
			return 1;

		tquery >>= LUTMAP_GROUPSIZE;		
	}

	int lSamples = (int)numberOfSamples - ((int)patternSize-1)*64;
	int lGroups = lSamples%LUTMAP_GROUPSIZE==0?lSamples/LUTMAP_GROUPSIZE:lSamples/LUTMAP_GROUPSIZE+1;

	// last iters
	tquery = ~query[patternSize-1];

	int shiftLast = 64-lSamples;
	tquery = tquery << shiftLast;
	tquery = tquery >> shiftLast;	


	for(i=0;i<lGroups;i++)
	{
		uint8_t lutMapAddr = tquery & 0xff;

		if(lm->lutMap[i][lutMapAddr]==0)
			return 1;

		tquery >>= LUTMAP_GROUPSIZE;		
	}

	return 0;
}

void RSDLutMap_update (RSDLutMap_t * lm, uint64_t * query) 
{
	assert(lm!=NULL);
	assert(query!=NULL);

	int i=0, y=0;

	uint64_t tquery = query[0];
	for(i=0;i<lm->lutMapSize;i++)
	{
		tquery = i%LUTMAP_INTERVAL==0?query[y++]:tquery;
		uint8_t lutMapAddr = tquery & 0xff;

		lm->lutMap[i][lutMapAddr]=1;

		tquery >>= LUTMAP_GROUPSIZE;		
	}
}



