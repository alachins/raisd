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

static RSDTreeNode_t * (*RSDTreeNode_new) (RSDTreeMap_t * RSDTreeMap);
void RSDTreeNode_free (RSDTreeNode_t * RSDTreeNode);
RSDTreeNode_t * RSDTreeNodePool_getNode (RSDTreeNodePool_t * RSDTreeNodePool);
RSDTreeNode_t * RSDTreeNode_new_0 (RSDTreeMap_t * RSDTreeMap);
RSDTreeNode_t * RSDTreeNode_new_1 (RSDTreeMap_t * RSDTreeMap);
RSDTreeNodePool_t * RSDTreeNodePool_new (void);
void RSDTreeNodePool_free (RSDTreeNodePool_t * np);

RSDTreeNode_t * RSDTreeNodePool_getNode (RSDTreeNodePool_t * RSDTreeNodePool)
{
	assert(RSDTreeNodePool!=NULL);

	RSDTreeNode_t * RSDTreeNode = NULL;

	RSDTreeNodePool->nextNodeDimY++;

	if(RSDTreeNodePool->nextNodeDimY==RSDTreeNodePool->nodeMatrixSizeY)
	{
		RSDTreeNodePool->nodeMatrixSizeX++;

		RSDTreeNodePool->nodeMatrix = (RSDTreeNode_t **)rsd_realloc(RSDTreeNodePool->nodeMatrix, sizeof(RSDTreeNode_t *)*((unsigned long)RSDTreeNodePool->nodeMatrixSizeX));
		assert(RSDTreeNodePool->nodeMatrix!=NULL);
	
		RSDTreeNodePool->nodeMatrix[RSDTreeNodePool->nodeMatrixSizeX-1] = (RSDTreeNode_t *)rsd_malloc(sizeof(RSDTreeNode_t)*((unsigned long)RSDTreeNodePool->nodeMatrixSizeY));
		assert(RSDTreeNodePool->nodeMatrix[RSDTreeNodePool->nodeMatrixSizeX-1]!=NULL);

		RSDTreeNodePool->nextNodeDimY = 0;	
	}

	RSDTreeNode = &(RSDTreeNodePool->nodeMatrix[RSDTreeNodePool->nodeMatrixSizeX-1][RSDTreeNodePool->nextNodeDimY]);
	assert(RSDTreeNode!=NULL);

	return RSDTreeNode;
}

RSDTreeNode_t * RSDTreeNode_new_0 (RSDTreeMap_t * RSDTreeMap)
{
#ifdef _TM_NODE_POOL
	RSDTreeNode_t * tn = RSDTreeNodePool_getNode (RSDTreeMap->nodePool[0]);
	assert(tn!=NULL);

	RSDTreeMap->totalMemory += sizeof(RSDTreeNode_t);
#else
	RSDTreeNode_t * tn = (RSDTreeNode_t *)rsd_malloc(sizeof(RSDTreeNode_t));
	assert(tn!=NULL);

	RSDTreeMap->totalMemory += sizeof(RSDTreeNode_t);
#endif

	RSDTreeMap->totalNodes++;

	tn->childNode[0] = NULL;
	tn->childNode[1] = NULL;

#ifndef _TM_PATTERN_ARRAY
	tn->patternID = -1;
#endif

	return tn;
}

RSDTreeNode_t * RSDTreeNode_new_1 (RSDTreeMap_t * RSDTreeMap)
{
#ifdef _TM_NODE_POOL
	RSDTreeNode_t * tn = RSDTreeNodePool_getNode (RSDTreeMap->nodePool[1]);
	assert(tn!=NULL);

	RSDTreeMap->totalMemory += sizeof(RSDTreeNode_t);
#else
	RSDTreeNode_t * tn = (RSDTreeNode_t *)rsd_malloc(sizeof(RSDTreeNode_t));
	assert(tn!=NULL);

	RSDTreeMap->totalMemory += sizeof(RSDTreeNode_t);
#endif

	RSDTreeMap->totalNodes++;

	tn->childNode[0] = NULL;
	tn->childNode[1] = NULL;

#ifndef _TM_PATTERN_ARRAY
	tn->patternID = -1;
#endif

	return tn;
}

void RSDTreeNode_free (RSDTreeNode_t * RSDTreeNode)
{
	if(RSDTreeNode==NULL)
		return;

	RSDTreeNode_free (RSDTreeNode->childNode[0]);
	RSDTreeNode_free (RSDTreeNode->childNode[1]);

	free(RSDTreeNode);
	RSDTreeNode = NULL;
}

RSDTreeNodePool_t * RSDTreeNodePool_new (void)
{
	RSDTreeNodePool_t * np = NULL;

	np = (RSDTreeNodePool_t *)rsd_malloc(sizeof(RSDTreeNodePool_t));
	assert(np!=NULL);

	np->nodeMatrixSizeX = 1;
	np->nodeMatrixSizeY = TREEMAP_NODEPOOL_CHUNKSIZE;

	np->nextNodeDimX = -1;
	np->nextNodeDimY = -1;

	np->nodeMatrix = (RSDTreeNode_t **)rsd_malloc(sizeof(RSDTreeNode_t *)*((unsigned long)np->nodeMatrixSizeX));
	assert(np->nodeMatrix!=NULL);

	np->nodeMatrix[0] = (RSDTreeNode_t *)rsd_malloc(sizeof(RSDTreeNode_t)*((unsigned long)np->nodeMatrixSizeY));
	assert(np->nodeMatrix[0]!=NULL);
	
	return np;	
}

RSDTreeMap_t * RSDTreeMap_new (void)
{
	RSDTreeMap_t * tm = (RSDTreeMap_t *)rsd_malloc(sizeof(RSDTreeMap_t));
	assert(tm!=NULL);
	
	tm->totalPatterns = 0;
	tm->totalNodes = 0;
	tm->totalMemory = 0;

	RSDTreeNode_new = &RSDTreeNode_new_0; 

#ifdef _TM_NODE_POOL
	tm->nodePool[0] = RSDTreeNodePool_new ();
#ifdef _TM_NODE_POOL_DOUBLE
	tm->nodePool[1] = RSDTreeNodePool_new ();
#endif
	tm->rootNode[0] = RSDTreeNode_new (tm); // left subtree, starting with first bit 0
#ifdef _TM_NODE_POOL_DOUBLE
	RSDTreeNode_new = &RSDTreeNode_new_1;
#endif
	tm->rootNode[1] = RSDTreeNode_new (tm);  // right subtree, starting with first bit 1
#else
	tm->rootNode[0] = RSDTreeNode_new (tm); // left subtree, starting with first bit 0
	tm->rootNode[1] = RSDTreeNode_new (tm);  // right subtree, starting with first bit 1
#endif
	// _TM_PATTERN_ARRAY
	tm->maxPatterns = 0;
	tm->patternID = NULL;	

	return tm;
}

void RSDTreeNodePool_free (RSDTreeNodePool_t * np)
{
	assert(np!=NULL);

	int i;

	if(np->nodeMatrix != NULL)
	{
		for(i=0;i<np->nodeMatrixSizeX;i++)
		{
			if(np->nodeMatrix[i]!=NULL)
			{
				free(np->nodeMatrix[i]);
				np->nodeMatrix[i] = NULL;
			}
		}
	
		free (np->nodeMatrix);
		np->nodeMatrix = NULL;
	}

	free(np);
	np=NULL;
}

void RSDTreeMap_free (RSDTreeMap_t * RSDTreeMap)
{
	assert(RSDTreeMap!=NULL);

#ifdef _TM_NODE_POOL
	RSDTreeNodePool_free (RSDTreeMap->nodePool[0]);
#ifdef _TM_NODE_POOL_DOUBLE
	RSDTreeNodePool_free (RSDTreeMap->nodePool[1]);
#endif
#else
	RSDTreeNode_free (RSDTreeMap->rootNode[0]);
	RSDTreeNode_free (RSDTreeMap->rootNode[1]);
#endif

#ifdef _TM_PATTERN_ARRAY
	if(RSDTreeMap->patternID!=NULL)
	{
		free(RSDTreeMap->patternID);
		RSDTreeMap->patternID = NULL;

		RSDTreeMap->maxPatterns = 0;
	}
#endif

	free(RSDTreeMap);
	RSDTreeMap = NULL;
}

#ifdef _TM_PATTERN_ARRAY
static inline int64_t RSDTreeMap_setPatternID (RSDTreeMap_t * RSDTreeMap, RSDTreeNode_t * RSDTreeNode)
{
	assert(RSDTreeNode!=NULL);

	RSDTreeMap->totalPatterns++;
	
	if(RSDTreeMap->totalPatterns>RSDTreeMap->maxPatterns)
	{
		RSDTreeMap->maxPatterns += TREEMAP_REALLOC_INCR;
		RSDTreeMap->patternID = (RSDTreeNode_t **) rsd_realloc(RSDTreeMap->patternID, sizeof(RSDTreeNode_t *)*(unsigned long)RSDTreeMap->maxPatterns);
		assert(RSDTreeMap->patternID!=NULL);
	}

	RSDTreeMap->patternID[RSDTreeMap->totalPatterns-1] = RSDTreeNode;

	return RSDTreeMap->totalPatterns-1;	
}

static inline int64_t RSDTreeMap_getPatternID (RSDTreeMap_t * RSDTreeMap, RSDTreeNode_t * RSDTreeNode)
{
	int i;
	for(i=0;i!=RSDTreeMap->totalPatterns;++i)
		if(RSDTreeMap->patternID[i]==RSDTreeNode)
			return i;

	return -1;
}
#endif

int64_t RSDTreeMap_matchSNP (RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples)
{
	assert(RSDTreeMap!=NULL);
	assert(RSDPatternPool_g!=NULL);
	assert(numberOfSamples>=1);

	RSDPatternPool_t * RSDPatternPool = (RSDPatternPool_t *) RSDPatternPool_g;

	int64_t i=0;
	uint8_t curBit = (uint8_t)(RSDPatternPool->incomingSite[i]-48);
	assert(curBit==0 || curBit==1);

	assert(RSDTreeMap->rootNode[curBit]!=NULL);

	RSDTreeNode_t * RSDTreeNode = RSDTreeMap->rootNode[curBit];

	for(i=0;i<RSDPatternPool->patternSize-1;i++)
	{
		uint64_t patternPart = RSDPatternPool->incomingSiteCompact[i];

		int j;
		for(j=0;j!=64;++j)
		{
			uint64_t tpatternPart = patternPart >> (63-j); 
			curBit = tpatternPart & 0x00000001;

			if(RSDTreeNode->childNode[curBit]==NULL)	
				return -1;	

			RSDTreeNode = RSDTreeNode->childNode[curBit];				
		}		
	}

	assert(i==RSDPatternPool->patternSize-1);

	{
		uint64_t patternPart = RSDPatternPool->incomingSiteCompact[i];
		int j;
		for(j=0;j!=numberOfSamples-(RSDPatternPool->patternSize-1)*64;++j)
		{
			uint64_t tpatternPart = patternPart >> ((numberOfSamples-(RSDPatternPool->patternSize-1)*64)-1-j);
			curBit = tpatternPart & 0x00000001;

			if(RSDTreeNode->childNode[curBit]==NULL)
				return -1;

			RSDTreeNode = RSDTreeNode->childNode[curBit];			
		}
	} 

#ifdef _TM_PATTERN_ARRAY
	return RSDTreeMap_getPatternID (RSDTreeMap, RSDTreeNode);
#else
	return RSDTreeNode->patternID;
#endif
}

int64_t RSDTreeMap_matchSNPC (RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples)
{
	assert(RSDTreeMap!=NULL);
	assert(RSDPatternPool_g!=NULL);
	assert(numberOfSamples>=1);

	RSDPatternPool_t * RSDPatternPool = (RSDPatternPool_t *) RSDPatternPool_g;

	int64_t i=0;
	uint8_t curBit = (uint8_t)(RSDPatternPool->incomingSite[i]-48);
	assert(curBit==0 || curBit==1);

	curBit = 1 - curBit;

	assert(RSDTreeMap->rootNode[curBit]!=NULL);

	RSDTreeNode_t * RSDTreeNode = RSDTreeMap->rootNode[curBit];

	for(i=0;i<RSDPatternPool->patternSize-1;i++)
	{
		uint64_t patternPart = RSDPatternPool->incomingSiteCompact[i];

		int j;
		for(j=0;j<64;j++)
		{
			uint64_t tpatternPart = patternPart >> (63-j); 
			curBit = tpatternPart & 0x00000001;

			curBit = 1 - curBit;

			if(RSDTreeNode->childNode[curBit]==NULL)
				return -1;

			RSDTreeNode = RSDTreeNode->childNode[curBit];
		}
	}

	assert(i==RSDPatternPool->patternSize-1);

	{
		uint64_t patternPart = RSDPatternPool->incomingSiteCompact[i];
		int j;
		for(j=0;j<numberOfSamples-(RSDPatternPool->patternSize-1)*64;j++)
		{
			uint64_t tpatternPart = patternPart >> ((numberOfSamples-(RSDPatternPool->patternSize-1)*64)-1-j);
			curBit = tpatternPart & 0x00000001;

			curBit = 1 - curBit;

			if(RSDTreeNode->childNode[curBit]==NULL)
				return -1;

			RSDTreeNode = RSDTreeNode->childNode[curBit];
		}
	} 

#ifdef _TM_PATTERN_ARRAY
	return RSDTreeMap_getPatternID (RSDTreeMap, RSDTreeNode);
#else
	return RSDTreeNode->patternID;
#endif
}

int64_t RSDTreeMap_updateTree (RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples)
{
	assert(RSDTreeMap!=NULL);
	assert(RSDPatternPool_g!=NULL);
	assert(numberOfSamples>=1);

	RSDPatternPool_t * RSDPatternPool = (RSDPatternPool_t *) RSDPatternPool_g;

	int64_t i=0;

	uint8_t curBit = (uint8_t)(RSDPatternPool->incomingSite[i]-48);
	assert(curBit==0 || curBit==1);

	assert(RSDTreeMap->rootNode[curBit]!=NULL);

	RSDTreeNode_t * RSDTreeNode = RSDTreeMap->rootNode[curBit];

	RSDTreeNode_new = &RSDTreeNode_new_0;

#ifdef _TM_NODE_POOL_DOUBLE
	if(curBit==1)
		RSDTreeNode_new = &RSDTreeNode_new_1;
#endif

	for(i=0;i<RSDPatternPool->patternSize-1;i++)
	{
		uint64_t patternPart = RSDPatternPool->incomingSiteCompact[i];

		int j;
		for(j=0;j!=64;++j)
		{
			uint64_t tpatternPart = patternPart >> (63-j);
			curBit = tpatternPart & 0x00000001;

			if(RSDTreeNode->childNode[curBit]==NULL)
				RSDTreeNode->childNode[curBit] = RSDTreeNode_new(RSDTreeMap);

			RSDTreeNode = RSDTreeNode->childNode[curBit];
		}
	}

	assert(i==RSDPatternPool->patternSize-1);

	{
		uint64_t patternPart = RSDPatternPool->incomingSiteCompact[i];

		int j;
		for(j=0;j<numberOfSamples-(RSDPatternPool->patternSize-1)*64;j++)
		{
			uint64_t tpatternPart = patternPart >> ((numberOfSamples-(RSDPatternPool->patternSize-1)*64)-1-j);
			curBit = tpatternPart & 0x00000001;

			if(RSDTreeNode->childNode[curBit]==NULL)
				RSDTreeNode->childNode[curBit] = RSDTreeNode_new(RSDTreeMap);

			RSDTreeNode = RSDTreeNode->childNode[curBit];
		}
	}

#ifdef _TM_PATTERN_ARRAY
	return RSDTreeMap_setPatternID (RSDTreeMap, RSDTreeNode);
#else
	assert(RSDTreeNode->patternID==-1);
	RSDTreeNode->patternID = RSDTreeMap->totalPatterns++;
	return RSDTreeNode->patternID;
#endif
}


int64_t RSDTreeMap_updateTreeInit (RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples, uint64_t * pattern)
{
	assert(RSDTreeMap!=NULL);
	assert(RSDPatternPool_g!=NULL);
	assert(numberOfSamples>=1);

	RSDPatternPool_t * RSDPatternPool = (RSDPatternPool_t *) RSDPatternPool_g;

	int64_t i=0;

	uint64_t tpatternPart = pattern[0] >> 63; 

	if(numberOfSamples<64)
		tpatternPart = pattern[0] >> (numberOfSamples-1);
	 
	uint8_t curBit = tpatternPart & 0x00000001;

	assert(RSDTreeMap->rootNode[curBit]!=NULL);

	RSDTreeNode_t * RSDTreeNode = RSDTreeMap->rootNode[curBit];

	RSDTreeNode_new = &RSDTreeNode_new_0;

#ifdef _TM_NODE_POOL_DOUBLE
	if(curBit==1)
		RSDTreeNode_new = &RSDTreeNode_new_1;
#endif

	for(i=0;i<RSDPatternPool->patternSize-1;i++)
	{
		uint64_t patternPart = pattern[i];

		int j;
		for(j=0;j<64;j++)
		{
			tpatternPart = patternPart >> (63-j);
			curBit = tpatternPart & 0x00000001;

			if(RSDTreeNode->childNode[curBit]==NULL)
				RSDTreeNode->childNode[curBit] = RSDTreeNode_new(RSDTreeMap);

			RSDTreeNode = RSDTreeNode->childNode[curBit];	
		}
	}

	assert(i==RSDPatternPool->patternSize-1);

	{
		uint64_t patternPart = pattern[i];
		int j;
		for(j=0;j<numberOfSamples-(RSDPatternPool->patternSize-1)*64;j++)
		{
			tpatternPart = patternPart >> ((numberOfSamples-(RSDPatternPool->patternSize-1)*64)-1-j);
			curBit = tpatternPart & 0x00000001;		

			if(RSDTreeNode->childNode[curBit]==NULL)
				RSDTreeNode->childNode[curBit] = RSDTreeNode_new(RSDTreeMap);
		
			RSDTreeNode = RSDTreeNode->childNode[curBit];
		}
	}

#ifdef _TM_PATTERN_ARRAY
	return RSDTreeMap_setPatternID (RSDTreeMap, RSDTreeNode);
#else
	assert(RSDTreeNode->patternID==-1);
	RSDTreeNode->patternID = RSDTreeMap->totalPatterns++;
	return RSDTreeNode->patternID;
#endif
}
