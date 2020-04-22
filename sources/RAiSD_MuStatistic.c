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

//float getPatternCounts (int * pCntVec, int offset, int * patternID, int p0, int p1, int p2, int p3, int * pcntl, int * pcntr, int * pcntexll, int * pcntexlr);
float 	getPatternCounts (int winMode, RSDMuStat_t * RSDMuStat, int sizeL, int sizeR, int * patternID, int p0, int p1, int p2, int p3, int * pcntl, int * pcntr, int * pcntexll, int * pcntexlr);
float 	getCrossLD 			(RSDPatternPool_t * pp, int p0, int p1, int p2, int p3, int samples);
float 	getRegionLD 			(RSDPatternPool_t * pp, int p0, int p1, int samples);
float 	pwLD 				(RSDPatternPool_t * pp, int p1, int p2, int samples);
void	(*RSDMuStat_scanChunk) 		(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
void 	RSDMuStat_scanChunkBinary	(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
void 	RSDMuStat_scanChunkWithMask	(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
void 	(*RSDMuStat_storeOutput) 	(RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu);
void 	RSDMuStat_output2FileSimple 	(RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu);
void 	RSDMuStat_output2FileFull 	(RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu);
void 	RSDMuStat_output2BufferSimple 	(RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu);
void 	RSDMuStat_output2BufferFull 	(RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu);




RSDMuStat_t * RSDMuStat_new (void)
{
	RSDMuStat_t * mu = NULL;

	mu = (RSDMuStat_t *)rsd_malloc(sizeof(RSDMuStat_t));
	assert(mu!=NULL);

	mu->reportFP = NULL;
	
	mu->windowSize = DEFAULT_WINDOW_SIZE;
	
	mu->pCntVec = NULL; 

	mu->muVarMax = 0.0f; 
	mu->muVarMaxLoc = 0.0;

	mu->muSfsMax = 0.0f; 
	mu->muSfsMaxLoc = 0.0;

	mu->muLdMax = 0.0f;
	mu->muLdMaxLoc = 0.0;

	mu->muMax = 0.0f; 
	mu->muMaxLoc = 0.0;

#ifdef _MUMEM
	mu->muReportBufferSize = 1;
	mu->muReportBuffer = (float*)rsd_malloc(sizeof(float)*7);
	assert(mu->muReportBuffer!=NULL);
#endif

	mu->exclTableSize = 0;
	mu->exclTableChromName = NULL;
	mu->exclTableRegionStart = NULL;
	mu->exclTableRegionStop = NULL;

	mu->excludeRegionsTotal = 0;
	mu->excludeRegionStart = NULL;
	mu->excludeRegionStop = NULL;

	mu->bufferMemMaxSize = 0;
	mu->buffer0Data = NULL; // windowCenter
	mu->buffer1Data = NULL; // windowStart	
	mu->buffer2Data = NULL; // windowEnd
	mu->buffer3Data = NULL; // muVAR
	mu->buffer4Data = NULL; // muSFS
	mu->buffer5Data = NULL; // muLD
	mu->buffer6Data = NULL; // mu
	mu->currentScoreIndex = -1;

	return mu;
}

void RSDMuStat_free (RSDMuStat_t * mu)
{
	assert(mu!=NULL);

	if(mu->reportFP!=NULL)
	{
		fclose(mu->reportFP);
		mu->reportFP = NULL;
	}
	
	//MemoryFootprint += sizeof(int)*((unsigned long)(mu->windowSize*4));
	if(mu->pCntVec!=NULL)
		free(mu->pCntVec);

#ifdef _MUMEM	
	if(mu->muReportBuffer!=NULL)
	{
		//MemoryFootprint += sizeof(float)*((unsigned long)(mu->muReportBufferSize*7));

		free(mu->muReportBuffer);
		mu->muReportBuffer =  NULL;
	}
#endif

	if(mu->exclTableSize!=0)
	{
		if(mu->exclTableChromName!=NULL)
		{
			int i;
			for(i=0;i<mu->exclTableSize;i++)
			{
				if(mu->exclTableChromName[i]!=NULL)
				{
					free(mu->exclTableChromName[i]);
					mu->exclTableChromName[i] = NULL;
				}
			}
		
			free(mu->exclTableChromName);
			mu->exclTableChromName = NULL;
		}

		if(mu->exclTableRegionStart!=NULL)
		{
			free(mu->exclTableRegionStart);
			mu->exclTableRegionStart = NULL;
		}

		if(mu->exclTableRegionStop!=NULL)
		{
			free(mu->exclTableRegionStop);
			mu->exclTableRegionStop = NULL;
		}		
	}

	if(mu->excludeRegionStart!=NULL)
	{
		free(mu->excludeRegionStart);
		mu->excludeRegionStart = NULL;
	}

	if(mu->excludeRegionStop!=NULL)
	{
		free(mu->excludeRegionStop);
		mu->excludeRegionStop = NULL;
	}

	if(mu->buffer0Data!=NULL)
	{
		free(mu->buffer0Data);
		mu->buffer0Data = NULL;
	}

	if(mu->buffer1Data!=NULL)
	{
		free(mu->buffer1Data);
		mu->buffer1Data = NULL;
	}

	if(mu->buffer2Data!=NULL)
	{
		free(mu->buffer2Data);
		mu->buffer2Data = NULL;
	}

	if(mu->buffer3Data!=NULL)
	{
		free(mu->buffer3Data);
		mu->buffer3Data = NULL;
	}

	if(mu->buffer4Data!=NULL)
	{
		free(mu->buffer4Data);
		mu->buffer4Data = NULL;
	}

	if(mu->buffer5Data!=NULL)
	{
		free(mu->buffer5Data);
		mu->buffer5Data = NULL;
	}

	if(mu->buffer6Data!=NULL)
	{
		free(mu->buffer6Data);
		mu->buffer6Data = NULL;
	}

	free(mu);
}

void RSDMuStat_init (RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);

	RSDMuStat->muVarMax = 0.0f; 
	RSDMuStat->muVarMaxLoc = 0.0;

	RSDMuStat->muSfsMax = 0.0f; 
	RSDMuStat->muSfsMaxLoc = 0.0;

	RSDMuStat->muLdMax = 0.0f;
	RSDMuStat->muLdMaxLoc = 0.0;

	RSDMuStat->muMax = 0.0f; 
	RSDMuStat->muMaxLoc = 0.0;

	RSDMuStat->windowSize = RSDCommandLine->windowSize;

	if(RSDMuStat->pCntVec==NULL)
	{
		RSDMuStat->pCntVec = (int *)rsd_malloc(sizeof(int)*((unsigned long)(RSDMuStat->windowSize*4)));
		assert(RSDMuStat->pCntVec!=NULL);
	}

	if(RSDCommandLine->createPatternPoolMask==1)
		RSDMuStat_scanChunk = &RSDMuStat_scanChunkWithMask;
	else
		RSDMuStat_scanChunk = &RSDMuStat_scanChunkBinary;
		
}

void RSDMuStat_setReportName (RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine, FILE * fpOut)
{
	strcpy(RSDMuStat->reportName, "RAiSD_Report.");
	strcat(RSDMuStat->reportName, RSDCommandLine->runName);

	if(RSDCommandLine->splitOutput==1)
		return;
	
	RSDMuStat->reportFP = fopen(RSDMuStat->reportName, "r");
	if(RSDMuStat->reportFP!=NULL && RSDCommandLine->overwriteOutput==0)
	{
		fprintf(fpOut, "\nERROR: Output report file %s exists. Use -f to overwrite it.\n\n", RSDMuStat->reportName);
		fprintf(stderr, "\nERROR: Output report file %s exists. Use -f to overwrite it.\n\n", RSDMuStat->reportName);
		exit(0);	
	}

	if(RSDMuStat->reportFP!=NULL)
		fclose(RSDMuStat->reportFP);

	RSDMuStat->reportFP = fopen(RSDMuStat->reportName, "w");
	assert(RSDMuStat->reportFP!=NULL);
}

void RSDMuStat_loadExcludeTable (RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDMuStat!=NULL);
	assert(RSDCommandLine!=NULL);

	if(!strcmp(RSDCommandLine->excludeRegionsFile, "\0"))
		return;

	FILE * fp = fopen(RSDCommandLine->excludeRegionsFile, "r");
	assert(fp!=NULL);

	char tstring[STRING_SIZE];

	int rcnt = fscanf(fp, "%s", tstring);
	assert(rcnt==1);

	while(rcnt==1)
	{	
		RSDMuStat->exclTableSize++;
		
		RSDMuStat->exclTableChromName = rsd_realloc(RSDMuStat->exclTableChromName, sizeof(char*)*((unsigned long)RSDMuStat->exclTableSize));
		assert(RSDMuStat->exclTableChromName!=NULL);
		
		RSDMuStat->exclTableChromName[RSDMuStat->exclTableSize-1] = (char*)rsd_malloc(sizeof(char)*STRING_SIZE);
		assert(RSDMuStat->exclTableChromName[RSDMuStat->exclTableSize-1]!=NULL);

		RSDMuStat->exclTableRegionStart = rsd_realloc(RSDMuStat->exclTableRegionStart, sizeof(int64_t)*((unsigned long)RSDMuStat->exclTableSize));
		assert(RSDMuStat->exclTableRegionStart!=NULL);

		RSDMuStat->exclTableRegionStop = rsd_realloc(RSDMuStat->exclTableRegionStop, sizeof(int64_t)*((unsigned long)RSDMuStat->exclTableSize));
		assert(RSDMuStat->exclTableRegionStop!=NULL);

		strncpy(RSDMuStat->exclTableChromName[RSDMuStat->exclTableSize-1], tstring, STRING_SIZE);

		rcnt = fscanf(fp, "%s", tstring);
		assert(rcnt==1);

		RSDMuStat->exclTableRegionStart[RSDMuStat->exclTableSize-1] = (int64_t)strtoull(tstring, NULL, 0);

		rcnt = fscanf(fp, "%s", tstring);
		assert(rcnt==1);

		RSDMuStat->exclTableRegionStop[RSDMuStat->exclTableSize-1] = (int64_t)strtoull(tstring, NULL, 0);

		rcnt = fscanf(fp, "%s", tstring);
	}

	fclose(fp);
	fp=NULL;
}

void RSDMuStat_excludeRegion (RSDMuStat_t * RSDMuStat, RSDDataset_t * RSDDataset)
{
	assert(RSDMuStat!=NULL);
	assert(RSDDataset!=NULL);

	if(RSDMuStat->excludeRegionsTotal!=0 || RSDMuStat->excludeRegionStart!=NULL || RSDMuStat->excludeRegionStop!=NULL)
	{
		if(RSDMuStat->excludeRegionStart!=NULL)
			free(RSDMuStat->excludeRegionStart);

		if(RSDMuStat->excludeRegionStop!=NULL)
			free(RSDMuStat->excludeRegionStop);
	}

	RSDMuStat->excludeRegionsTotal = 0;
	RSDMuStat->excludeRegionStart = NULL;
	RSDMuStat->excludeRegionStop = NULL;

	if(RSDMuStat->exclTableSize==0)
		return;

	int i;

	for(i=0;i<RSDMuStat->exclTableSize;i++)
	{
		if(!strcmp(RSDMuStat->exclTableChromName[i], RSDDataset->setID))
		{
			RSDMuStat->excludeRegionsTotal++;
			RSDMuStat->excludeRegionStart = rsd_realloc(RSDMuStat->excludeRegionStart, sizeof(int64_t)*((unsigned long)RSDMuStat->excludeRegionsTotal));
			assert(RSDMuStat->excludeRegionStart!=NULL);
			RSDMuStat->excludeRegionStart[RSDMuStat->excludeRegionsTotal-1] = RSDMuStat->exclTableRegionStart[i];
			RSDMuStat->excludeRegionStop = rsd_realloc(RSDMuStat->excludeRegionStop, sizeof(int64_t)*((unsigned long)RSDMuStat->excludeRegionsTotal));
			assert(RSDMuStat->excludeRegionStop!=NULL);
			RSDMuStat->excludeRegionStop[RSDMuStat->excludeRegionsTotal-1] = RSDMuStat->exclTableRegionStop[i];
		}
	}
	
}

void RSDMuStat_setReportNamePerSet (RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine, FILE * fpOut, RSDDataset_t * RSDDataset, RSDCommonOutliers_t * RSDCommonOutliers)
{
	assert(fpOut!=NULL);

	if(RSDCommandLine->splitOutput==0)
		return;

	if(RSDMuStat->reportFP!=NULL)
	{
		fclose(RSDMuStat->reportFP);
		RSDMuStat->reportFP = NULL;
	}

	char tstring[STRING_SIZE];
	strcpy(tstring, RSDMuStat->reportName);
	strcat(tstring, ".");
	strcat(tstring, RSDDataset->setID);
	
	strcpy(RSDMuStat->reportFPFileName, tstring);

	if(RSDCommandLine->createMPlot==1)
	{
		fprintf(RAiSD_ReportList_FP, "%s\n", tstring);
		fflush(RAiSD_ReportList_FP);
	}	

	RSDMuStat->reportFP = fopen(tstring, "r");
	if(RSDMuStat->reportFP!=NULL && RSDCommandLine->overwriteOutput==0)
	{
		fprintf(fpOut, "\nERROR: Output report file %s exists. Use -f to overwrite it.\n\n", tstring);
		fprintf(stderr, "\nERROR: Output report file %s exists. Use -f to overwrite it.\n\n", tstring);
		exit(0);	
	}

	if(RSDMuStat->reportFP!=NULL)
	{
		fclose(RSDMuStat->reportFP);
		RSDMuStat->reportFP = NULL;
	}

	RSDMuStat->reportFP = fopen(tstring, "w");
	assert(RSDMuStat->reportFP!=NULL);

	if(!strcmp(RSDCommonOutliers->reportFilenameRAiSD, "\0"))
	{
		strncpy(RSDCommonOutliers->reportFilenameRAiSD, tstring, STRING_SIZE);
		RSDCommonOutliers->positionIndexRAiSD = 1;
		RSDCommonOutliers->scoreIndexRAiSD = 7;
	}
} 

float pwLD (RSDPatternPool_t * pp, int p1, int p2, int samples)
{
	float ld = 0.0f;

	int i, cntL=0, cntR=0, cntLR=0;
	for(i=0;i<pp->patternSize;i++)
	{
		cntL += rsd_popcnt_u64(pp->poolData[p1*pp->patternSize+i]);
		cntR += rsd_popcnt_u64(pp->poolData[p2*pp->patternSize+i]);
		cntLR += rsd_popcnt_u64(pp->poolData[p1*pp->patternSize+i] & pp->poolData[p2*pp->patternSize+i]);
	}

	float f1 = ((float)cntL)/((float)samples);
	float f2 = ((float)cntR)/((float)samples);
	float f12 = ((float)cntLR)/((float)samples);

	float num = (f12 - f1*f2)*(f12 - f1*f2);
	float den =  f1 * (1.0f-f1) * f2 * (1.0f-f2);

	ld = num/den;

	return ld;
}

float getRegionLD (RSDPatternPool_t * pp, int p0, int p1, int samples)
{
	double totalLD = 0.0;
	int i, j;
	for(i=p0;i<=p1-1;i++)
	{
		for(j=p0+1;j<=p1;j++)
			totalLD += (double) pwLD(pp, i, j, samples);
	}
	return (float)totalLD;
}

float getCrossLD (RSDPatternPool_t * pp, int p0, int p1, int p2, int p3, int samples)
{
	double totalLD = 0.0;
	int i, j;
	for(i=p0;i<=p1;i++)
	{
		for(j=p2;j<=p3;j++)
			totalLD += (double) pwLD(pp, i, j, samples);
	}
	return (float)totalLD;
}

/*#ifdef _REF
float getPatternCounts (int * pCntVec, int offset, int * patternID, int p0, int p1, int p2, int p3, int * pcntl, int * pcntr, int * pcntexll, int * pcntexlr)
{
	int i, j;

	int list_left_size = 0;
	int * list_left = &(pCntVec[0]); 
	int * list_left_cnt = &(pCntVec[offset]); 
	list_left[0] = patternID[p0];
	list_left_cnt[0] = 1;
	list_left_size++;

	// Left subwindow
	for(i=p0+1;i<=p1;i++)
	{
		int match = 0;
		for(j=0;j<list_left_size;j++)
		{
			if(list_left[j]==patternID[i])
			{
				match = 1;
				list_left_cnt[j]++;
				
			}
		}

		if(match==0)
		{
			list_left_size++;
			list_left[list_left_size-1] = patternID[i];
			list_left_cnt[list_left_size-1] = 1;
		}
	}
		
	*pcntl = list_left_size; 

	int checksum = 0;
	for(i=0;i<list_left_size;i++)
		checksum += list_left_cnt[i];

	assert(checksum==p1-p0+1);


	int list_right_size = 0;
	int * list_right = &(pCntVec[2*offset]); 
	int * list_right_cnt = &(pCntVec[3*offset]); 
	list_right[0] = patternID[p2];
	list_right_cnt[0] = 1;
	list_right_size++;	

	for(i=p2+1;i<=p3;i++)
	{
		int match = 0;
		for(j=0;j<list_right_size;j++)
		{
			if(list_right[j]==patternID[i])
			{
				match = 1;
				list_right_cnt[j]++;
			}
		}

		if(match==0)
		{
			list_right_size++;
			list_right[list_right_size-1] = patternID[i];
			list_right_cnt[list_right_size-1] = 1;
		}
	}
		
	*pcntr = list_right_size;


	int excl_left = list_left_size;
	for(i=0;i<list_left_size;i++)
	{
		for(j=0;j<list_right_size;j++)
		{
			if(list_left[i] == list_right[j])
			{
				excl_left--;
			}
		}
	}
	*pcntexll = excl_left;

	int excl_right = list_right_size;
	for(i=0;i<list_right_size;i++)
	{
		for(j=0;j<list_left_size;j++)
		{
			if(list_right[i] == list_left[j])
			{
				excl_right--;
			}
		}
	}
	*pcntexlr = excl_right;
	assert(list_left_size - excl_left == list_right_size - excl_right);


	// Exclusive SNPs left 
	int excntsnpsl = 0;
	for(i=p0;i<=p1;i++)
	{
		int match = 0;
		for(j=p2;j<=p3;j++)
		{
			if(patternID[i]==patternID[j])
			{
				match = 1;
				j=p3+1;
				
			}
		}

		if(match==0)
		{
			excntsnpsl++;
		}
	}
	*pcntexll = (*pcntexll) * excntsnpsl;

	int excntsnpsr = 0;
	for(i=p2;i<=p3;i++)
	{
		int match = 0;
		for(j=p0;j<=p1;j++)
		{
			if(patternID[i]==patternID[j])
			{
				match = 1;
				j=p1+1;
				
			}
		}

		if(match==0)
		{
			excntsnpsr++;
		}
	}
	*pcntexlr = (*pcntexlr) * excntsnpsr;

	return 0;
}
#endif
*/

float getPatternCounts (int winMode, RSDMuStat_t * RSDMuStat, int sizeL, int sizeR, int * patternID, int p0, int p1, int p2, int p3, int * pcntl, int * pcntr, int * pcntexll, int * pcntexlr)
{
	int i, j;

	if(winMode!=WIN_MODE_FXD)
	{ // TODO: Realloc must be called when needed. Not per call. Fix this before release.
		if(RSDMuStat->pCntVec!=NULL)
		{
			free(RSDMuStat->pCntVec);
			RSDMuStat->pCntVec = NULL;
		}

		if(RSDMuStat->pCntVec==NULL)
		{
			RSDMuStat->pCntVec = (int *)rsd_malloc(sizeof(int)*((unsigned long)(RSDMuStat->windowSize*4)));
			assert(RSDMuStat->pCntVec!=NULL);
		}
	}

	int * pCntVec = RSDMuStat->pCntVec;
	assert(pCntVec!=NULL);

	int list_left_size = 0;
	int * list_left = &(pCntVec[0]); 
	int * list_left_cnt = &(pCntVec[sizeL]); 
	list_left[0] = patternID[p0];
	list_left_cnt[0] = 1;
	list_left_size++;

	// Left subwindow
	for(i=p0+1;i<=p1;i++)
	{
		int match = 0;
		for(j=0;j<list_left_size;j++)
		{
			if(list_left[j]==patternID[i])
			{
				match = 1;
				list_left_cnt[j]++;
				
			}
		}

		if(match==0)
		{
			list_left_size++;
			list_left[list_left_size-1] = patternID[i];
			list_left_cnt[list_left_size-1] = 1;
		}
	}
		
	*pcntl = list_left_size; 

	int checksum = 0;
	for(i=0;i<list_left_size;i++)
		checksum += list_left_cnt[i];

	assert(checksum==p1-p0+1);


	int list_right_size = 0;
	int * list_right = &(pCntVec[2*sizeL]); 
	int * list_right_cnt = &(pCntVec[2*sizeL+sizeR]); 
	list_right[0] = patternID[p2];
	list_right_cnt[0] = 1;
	list_right_size++;	

	for(i=p2+1;i<=p3;i++)
	{
		int match = 0;
		for(j=0;j<list_right_size;j++)
		{
			if(list_right[j]==patternID[i])
			{
				match = 1;
				list_right_cnt[j]++;
			}
		}

		if(match==0)
		{
			list_right_size++;
			list_right[list_right_size-1] = patternID[i];
			list_right_cnt[list_right_size-1] = 1;
		}
	}
		
	*pcntr = list_right_size;


	int excl_left = list_left_size;
	for(i=0;i<list_left_size;i++)
	{
		for(j=0;j<list_right_size;j++)
		{
			if(list_left[i] == list_right[j])
			{
				excl_left--;
			}
		}
	}
	*pcntexll = excl_left;

	int excl_right = list_right_size;
	for(i=0;i<list_right_size;i++)
	{
		for(j=0;j<list_left_size;j++)
		{
			if(list_right[i] == list_left[j])
			{
				excl_right--;
			}
		}
	}
	*pcntexlr = excl_right;
	assert(list_left_size - excl_left == list_right_size - excl_right);


	/* Exclusive SNPs left */
	/* */
	int excntsnpsl = 0;
	for(i=p0;i<=p1;i++)
	{
		int match = 0;
		for(j=p2;j<=p3;j++)
		{
			if(patternID[i]==patternID[j])
			{
				match = 1;
				j=p3+1;
				
			}
		}

		if(match==0)
		{
			excntsnpsl++;
		}
	}
	*pcntexll = (*pcntexll) * excntsnpsl;

	int excntsnpsr = 0;
	for(i=p2;i<=p3;i++)
	{
		int match = 0;
		for(j=p0;j<=p1;j++)
		{
			if(patternID[i]==patternID[j])
			{
				match = 1;
				j=p1+1;
				
			}
		}

		if(match==0)
		{
			excntsnpsr++;
		}
	}
	*pcntexlr = (*pcntexlr) * excntsnpsr;

	if(winMode!=WIN_MODE_FXD)
	{
		free(RSDMuStat->pCntVec);
		RSDMuStat->pCntVec = NULL;
	}

	return 0;
}

#ifdef _HW
// FPL2018 MuStatistic implementation (DAER)
void pru1_var_scan (float * vec_a_in, float * vec_a_tmp, float * window_loc, int iterations, int window_size, uint64_t region_length)
{
	int i;
	float muVar = 0.0f;

	// Processing
	for(i=0;i<iterations;i++)
	{
		muVar = vec_a_in[i+window_size-1] - vec_a_in[i];
		muVar /= region_length;
		muVar /= window_size;

		vec_a_tmp[i] = muVar;
	
		window_loc[i] = (vec_a_in[i+window_size-1] + vec_a_in[i]) / 2.0f;
	}
}

void pru2_sfs_scan (int * vec_b_in, float * vec_b_tmp, int iterations, int window_size, int sample_size)
{
	int i;
	float muSFS = 0.0f;

	// Preprocessing
	int dCnt1 = 0, dCntN = 0;
	for(i=0;i<window_size - 0;i++)
	{
		dCnt1 += (vec_b_in[i]==1);
		dCntN += (vec_b_in[i]==sample_size-1);
	}

	// Processing
	for(i=0;i<iterations;i++)
	{
		dCnt1 -= (vec_b_in[i-1]==1);
		dCnt1 += (vec_b_in[i+window_size-1]==1);

		dCntN -= (vec_b_in[i-1]==sample_size-1);
		dCntN += (vec_b_in[i+window_size-1]==sample_size-1);


		muSFS = (float)dCnt1 + (float)dCntN; 

		if(dCnt1+dCntN==0) 
			muSFS = 0.000001f;

		muSFS /= (float)window_size;

		vec_b_tmp[i] = muSFS;
	}	
}

void pru3_ld_scan (int * vec_c_in, float * vec_c_tmp, int iterations, int window_size, int sample_size)
{
	int i, j;
	float muLD = 0.0f;

	// Preprocessing
	int * tmp_vec_a = (int *)rsd_malloc(sizeof(int)*(window_size/2)); // patterns left
	int * tmp_vec_b = (int *)rsd_malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_c = (int *)rsd_malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_d = (int *)rsd_malloc(sizeof(int)*(window_size/2));

	int * tmp_vec_a_cnt = (int *)rsd_malloc(sizeof(int)*(window_size/2)); // patterns left
	int * tmp_vec_b_cnt = (int *)rsd_malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_c_cnt = (int *)rsd_malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_d_cnt = (int *)rsd_malloc(sizeof(int)*(window_size/2));
	
	for(i=0;i<(window_size/2);i++)
	{
		tmp_vec_a[i] = -1;
		tmp_vec_b[i] = -1;
		tmp_vec_c[i] = -1;
		tmp_vec_d[i] = -1;

		tmp_vec_a_cnt[i] = 0;
		tmp_vec_b_cnt[i] = 0;
		tmp_vec_c_cnt[i] = 0;
		tmp_vec_d_cnt[i] = 0;
	}

	// Temp Registers
	int reg1 = vec_c_in[0];
	int reg2 = vec_c_in[0 + (window_size/2) - 1];
	int reg3 = vec_c_in[0 + (window_size/2)];
	int reg4 = vec_c_in[0 + (window_size/2) + (window_size/2) - 1];

	// Left
	tmp_vec_a[0] = vec_c_in[0];
	tmp_vec_a_cnt[0] = 1;

	for( i=1 ; i<= 0 + (window_size/2) - 1; i++)
	{
		int match = 0;
		for(j=0;j<(window_size/2);j++)
		{
			if(tmp_vec_a[j]==vec_c_in[i])
			{
				match = 1;
				tmp_vec_a_cnt[j]++;
				
			}
		}

		if(match==0)
		{
			for(j=0;j<(window_size/2);j++)
				if(tmp_vec_a[j]==-1)
					break;

			
			tmp_vec_a[j] = vec_c_in[i];
			tmp_vec_a_cnt[j] = 1;
		}
	}

	// Right
	tmp_vec_b[0] = vec_c_in[(window_size/2)];
	tmp_vec_b_cnt[0] = 1;


	for( i=(window_size/2)+1 ; i<= (window_size/2) + (window_size/2) - 1; i++)
	{
		int match = 0;
		for(j=0;j<(window_size/2);j++)
		{
			if(tmp_vec_b[j]==vec_c_in[i])
			{
				match = 1;
				tmp_vec_b_cnt[j]++;
				
			}
		}

		if(match==0)
		{
			for(j=0;j<(window_size/2);j++)
				if(tmp_vec_b[j]==-1)
					break;

			
			tmp_vec_b[j] = vec_c_in[i];
			tmp_vec_b_cnt[j] = 1;
		}
	}

	// Excl Left, Excl Right
	int left_sz = 0;
	for(j=0;j<window_size/2;j++)
		if(tmp_vec_a[j]!=-1)
			left_sz++;

	int right_sz = 0;
	for(j=0;j<window_size/2;j++)
		if(tmp_vec_b[j]!=-1)
			right_sz++;

	int excl_left = left_sz;
	for(i=0;i<window_size/2;i++)
	{
		int match = 0;
		for(j=0;j<window_size/2;j++)
		{
			if((tmp_vec_a[i] == tmp_vec_b[j]) && tmp_vec_a[i]!=-1)
			{
				excl_left--;
			}
		}
	}

	int excl_right = right_sz;
	for(i=0;i<window_size/2;i++)
	{
		int match = 0;
		for(j=0;j<window_size/2;j++)
		{
			if((tmp_vec_a[i] == tmp_vec_b[j]) && tmp_vec_a[i]!=-1)
			{
				excl_right--;
			}
		}
	}

	
	muLD = (((float)excl_left)+((float)excl_right)) / ((float)(left_sz * right_sz));
	vec_c_tmp[0] = muLD;

	int ii;
	// Processing
	int match = 0;
	for(i=1;i<iterations;i++)
	{
		int snp_l_first = i;
		int snp_l_last = i + (window_size/2) - 1;

		int snp_r_first = i + (window_size/2);
		int snp_r_last = i + (window_size/2) + (window_size/2) - 1;


		// Remove previous Left
		for(j=0;j<(window_size/2);j++)
		{
			if(tmp_vec_a[j]==reg1)
			{
				tmp_vec_a_cnt[j]--;

				if(tmp_vec_a_cnt[j]==0)
				{
					tmp_vec_a[j] = -1;
					left_sz--;
				}				
			}
		}
		reg1 = vec_c_in[snp_l_first];
		

		// Add new Left
		match = 0;
		for(j=0;j<(window_size/2);j++)
		{
			if(tmp_vec_a[j]==vec_c_in[snp_l_last])
			{
				match = 1;
				tmp_vec_a_cnt[j]++;
				
			}
		}

		if(match==0)
		{
			for(j=0;j<(window_size/2);j++)
				if(tmp_vec_a[j]==-1)
					break;

			
			tmp_vec_a[j] = vec_c_in[snp_l_last];
			tmp_vec_a_cnt[j] = 1;
			left_sz++;
		}
		reg2 =  vec_c_in[snp_l_last];

		// Remove previous right
		for(j=0;j<(window_size/2);j++)
		{
			if(tmp_vec_b[j]==reg2)
			{
				tmp_vec_b_cnt[j]--;

				if(tmp_vec_b_cnt[j]==0)
				{
					tmp_vec_b[j] = -1;
					right_sz--;
				}				
			}
		}
		reg3 = vec_c_in[snp_r_first];

		// Add new Right
		match = 0;
		for(j=0;j<(window_size/2);j++)
		{
			if(tmp_vec_b[j]==vec_c_in[snp_r_last])
			{
				match = 1;
				tmp_vec_b_cnt[j]++;
				
			}
		}

		if(match==0)
		{
			for(j=0;j<(window_size/2);j++)
				if(tmp_vec_b[j]==-1)
					break;

			
			tmp_vec_b[j] = vec_c_in[snp_r_last];
			tmp_vec_b_cnt[j] = 1;
			right_sz++;
		}
		reg4 =  vec_c_in[snp_r_last];
		
		// Calculate excl counters		
		excl_left = left_sz;
		for(ii=0;ii<window_size/2;ii++)
		{
			int match = 0;
			for(j=0;j<window_size/2;j++)
			{
				if((tmp_vec_a[ii] == tmp_vec_b[j]) && tmp_vec_a[ii]!=-1)
				{
					excl_left--;
				}
			}
		}

		excl_right = right_sz;
		for(ii=0;ii<window_size/2;ii++)
		{
			int match = 0;
			for(j=0;j<window_size/2;j++)
			{
				if((tmp_vec_a[ii] == tmp_vec_b[j]) && tmp_vec_a[ii]!=-1)
				{
					excl_right--;
				}
			}
		}

		
		muLD = (((float)excl_left)+((float)excl_right)) / ((float)(left_sz * right_sz));
		vec_c_tmp[i] = muLD;
	}

	// Free
	free(tmp_vec_a);
	free(tmp_vec_b);
	free(tmp_vec_c);
	free(tmp_vec_d);
}

void pru4_mu_calc (float * vec_a_tmp, float * vec_b_tmp, float * vec_c_tmp, float * vec_out, int iterations)
{
	int i;
	for(i=0;i<iterations;i++)
		vec_out[i]  = vec_a_tmp[i] * vec_b_tmp[i] * vec_c_tmp[i]; 
}

void RSDMuStat_scanChunk (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	// Exec parameters
	int iterations = RSDChunk->chunkSize - RSDMuStat->windowSize + 1;
	int vec_size = RSDChunk->chunkSize;
	int window_size = RSDMuStat->windowSize;
	uint64_t  region_length = RSDDataset->setRegionLength;
	int sample_size = RSDDataset->setSamples;

	// Input vectors
	float * vec_a_in = (float *)rsd_malloc(sizeof(float)*vec_size);
	int * vec_b_in = (int *)rsd_malloc(sizeof(int)*vec_size);
	int * vec_c_in = (int *)rsd_malloc(sizeof(int)*vec_size);

	memcpy(vec_a_in, RSDChunk->sitePosition, sizeof(float)*vec_size);
	memcpy(vec_b_in, RSDChunk->derivedAlleleCount, sizeof(int)*vec_size);
	memcpy(vec_c_in, RSDChunk->patternID, sizeof(int)*vec_size);

	// Intermediate vectors
	float * vec_a_tmp = (float *)rsd_malloc(sizeof(float)*vec_size);
	float * vec_b_tmp = (float *)rsd_malloc(sizeof(float)*vec_size);
	float * vec_c_tmp = (float *)rsd_malloc(sizeof(float)*vec_size);

	// Output vectors
	float * window_loc = (float*)rsd_malloc(sizeof(float)*vec_size);
	float * vec_out = (float *)rsd_malloc(sizeof(float)*vec_size);

	// Processing Units
	pru1_var_scan (vec_a_in, vec_a_tmp, window_loc, iterations, window_size, region_length);
	pru2_sfs_scan (vec_b_in, vec_b_tmp, iterations, window_size, sample_size);
	pru3_ld_scan (vec_c_in, vec_c_tmp, iterations, window_size, sample_size);
	pru4_mu_calc (vec_a_tmp, vec_b_tmp, vec_c_tmp, vec_out, iterations);

	int i, size = RSDChunk->chunkSize;

	float windowCenter = 0.0f;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	for(i=0;i<size-RSDMuStat->windowSize+1;i++)
	{

		windowCenter = window_loc[i];
		mu =  vec_out[i];//muVar * muSfs * muLd;

		// MuVar Max
		if (muVar > RSDMuStat->muVarMax)
		{
			RSDMuStat->muVarMax = muVar;
			RSDMuStat->muVarMaxLoc = windowCenter;
		}

		// MuSfs Max
		if (muSfs > RSDMuStat->muSfsMax)
		{
			RSDMuStat->muSfsMax = muSfs;
			RSDMuStat->muSfsMaxLoc = windowCenter;
		}

		// MuLd Max
		if (muLd > RSDMuStat->muLdMax)
		{
			RSDMuStat->muLdMax = muLd;
			RSDMuStat->muLdMaxLoc = windowCenter;
		}

		// Mu Max
		if (mu > RSDMuStat->muMax)
		{
			RSDMuStat->muMax = mu;
			RSDMuStat->muMaxLoc = windowCenter;
		}
	
		fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", windowCenter, muVar, muSfs, muLd, mu);
	}

	// Free
	free(vec_a_in);
	free(vec_b_in); 
	free(vec_c_in); 
	free(vec_a_tmp);
	free(vec_b_tmp); 
	free(vec_c_tmp);
}
// End of FPL2018 MuStatistic implementation (DAER)
#else
#ifdef _MLT
// TRETS2019 MuStatistic implementation (DAER2+HAM+MLT)
void RSDMuStat_scanChunkBinary (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i = 0, j = 0, 
	    size = (int)RSDChunk->chunkSize,
	    dCnt1 = 0,
	    dCntN = 0,
	    prevDerivedAllele = -1,
            patternMemLeftSize = 0,
	    patternMemLeft[DEFAULT_WINDOW_SIZE],
	    patternMemRightSize = 0,
	    patternMemRight[DEFAULT_WINDOW_SIZE],
            patternMemExclusiveSize = 0,
	    patternMemLeftExclusive[DEFAULT_WINDOW_SIZE],
	    patternMemRightExclusive[DEFAULT_WINDOW_SIZE],
	    rotbufIndex = 0,
	    firstLeftID_prev=-1, 
	    firstRightID_prev=-1,
	    patternsLeftExclusive = 0, 
	    patternsRightExclusive = 0, 
	    snpsLeftExclusive = 0, 
	    snpsRightExclusive = 0,
	    pcntl = 0, 
	    pcntr = 0, 
	    pcntexll=0, 
	    pcntexlr=0,
	    sitePositionL=-1,
	    derivedAlleleL=-1,
	    patternIDL=-1,
	    sitePositionF=-1,
	    derivedAlleleF=-1,
	    firstLeftID=-1,
	    firstRightID=-1;

	float windowCenter = 0.0f, 
	      windowStart = 0.0f, 
	      windowEnd = 0.0f,
	      muVar = 0.0f,
	      muSfs = 0.0f,
	      muLd = 0.0f,
	      mu = 0.0f,
	      rotbuf[DEFAULT_WINDOW_SIZE*3]; // circular queue

#ifdef _MUMEM
	int muReportBufferIndex = -1;

	if(size>RSDMuStat->muReportBufferSize)
	{
		RSDMuStat->muReportBufferSize = size;
		RSDMuStat->muReportBuffer = rsd_realloc(RSDMuStat->muReportBuffer, sizeof(float)*(unsigned long)RSDMuStat->muReportBufferSize*7);
		assert(RSDMuStat->muReportBuffer);	
	}
#endif

	// Circular queue init
	memcpy(rotbuf, RSDChunk->chunkData, DEFAULT_WINDOW_SIZE*sizeof(float)*3);

	// SFS preprocessing
	for(i=0;i<RSDMuStat->windowSize;i++)
	{
		dCnt1 += (((int)RSDChunk->chunkData[i*3+1])==1);
		dCntN += (((int)RSDChunk->chunkData[i*3+1])==RSDDataset->setSamples-1);
	}
		
	// LD preprocessing: init
	for(i=0;i<DEFAULT_WINDOW_SIZE/2;i++)
	{
		patternMemLeft[i*2+0] = -1; // ID 
		patternMemLeft[i*2+1] = 0; // count
	
		patternMemLeftExclusive[i*2+0] = -1;
		patternMemLeftExclusive[i*2+1] = 0;

		patternMemRight[i*2+0] = -1;
		patternMemRight[i*2+1] = 0;

		patternMemRightExclusive[i*2+0] = -1;		
		patternMemRightExclusive[i*2+1] = 0;
	}

	// LD preprocessing: patternMemLeft
	for(i=0;i<RSDMuStat->windowSize/2;i++)
	{
		int k, match=0;
		for(k=0;k<patternMemLeftSize;k++)
		{
			if(patternMemLeft[k*2+0]==(int)RSDChunk->chunkData[i*3+2])
			{
				match=1;
				patternMemLeft[k*2+1]++;
			}	
		}
		if(match==0)
		{	k = patternMemLeftSize;
			patternMemLeft[k*2+0] = (int)RSDChunk->chunkData[i*3+2];
			patternMemLeft[k*2+1]++;
			patternMemLeftSize++;
		}		
	}

	// LD preprocessing: patternMemRight
	for(;i<RSDMuStat->windowSize;i++)
	{
		int k, match=0;
		for(k=0;k<patternMemRightSize;k++)
		{
			if(patternMemRight[k*2+0]==(int)RSDChunk->chunkData[i*3+2])
			{
				match=1;
				patternMemRight[k*2+1]++;
			}	
		}
		if(match==0)
		{	k = patternMemRightSize;
			patternMemRight[k*2+0] = (int)RSDChunk->chunkData[i*3+2];
			patternMemRight[k*2+1]++;
			patternMemRightSize++;
		}
	}

	// LD preprocessing: patternMemLeftExclusive
	patternMemExclusiveSize = 0;
	for(i=0;i<patternMemLeftSize;i++)
	{
		int k, match = 0;
		for(k=0;k<patternMemRightSize;k++)
		{
			if(patternMemLeft[i*2+0]==patternMemRight[k*2+0])
			{
				match=1;
				k = patternMemRightSize;
			}
		}
		if(match==0)
		{
			k = patternMemExclusiveSize;
			patternMemLeftExclusive[k*2+0] = patternMemLeft[i*2+0];
			patternMemLeftExclusive[k*2+1] += patternMemLeft[i*2+1];
			patternMemExclusiveSize++;
		}		
	}

	// LD preprocessing: patternMemRightExclusive
	patternMemExclusiveSize = 0;
	for(i=0;i<patternMemRightSize;i++)
	{
		int k, match = 0;
		for(k=0;k<patternMemLeftSize;k++)
		{
			if(patternMemRight[i*2+0]==patternMemLeft[k*2+0])
			{
				match=1;
				k = patternMemLeftSize;
			}
		}
		if(match==0)
		{
			k = patternMemExclusiveSize;
			patternMemRightExclusive[k*2+0] = patternMemRight[i*2+0];
			patternMemRightExclusive[k*2+1] += patternMemRight[i*2+1];
			patternMemExclusiveSize++;
		}
		
	}

	
	// Pattern count preprocessing
	for(i=0;i<DEFAULT_WINDOW_SIZE/2;i++)
	{
		pcntl += (patternMemLeft[i*2+1]!=0?1:0);
		pcntr += (patternMemRight[i*2+1]!=0?1:0);
		patternsLeftExclusive += (patternMemLeftExclusive[i*2+1]!=0?1:0);
		patternsRightExclusive += (patternMemRightExclusive[i*2+1]!=0?1:0);
		snpsLeftExclusive += patternMemLeftExclusive[i*2+1];
		snpsRightExclusive += patternMemRightExclusive[i*2+1]; 		
	}

	// Iteration 0
	i=(int)RSDMuStat->windowSize-1;
	{
		// mem access
		sitePositionL = (int)RSDChunk->chunkData[i*3]; 
		derivedAlleleL = (int)RSDChunk->chunkData[i*3+1];
		patternIDL = (int)RSDChunk->chunkData[i*3+2];

		// rotbuf access
		sitePositionF = (int)rotbuf[rotbufIndex*3]; 
		derivedAlleleF = (int)rotbuf[rotbufIndex*3+1];
		firstLeftID = (int)rotbuf[rotbufIndex*3+2]; 
		
		// we need to load the previous sliding step's right window's first pattern ID
		int rotbufIndex_temp = rotbufIndex + (DEFAULT_WINDOW_SIZE/2);

		if(rotbufIndex>=(DEFAULT_WINDOW_SIZE/2)-1)
			rotbufIndex_temp -= (DEFAULT_WINDOW_SIZE-1);

		firstRightID = (int)rotbuf[rotbufIndex_temp*3+2];

		// rotbuf update
		rotbuf[rotbufIndex*3] = (float)sitePositionL;		
		rotbuf[rotbufIndex*3+1] = (float)derivedAlleleL;
		rotbuf[rotbufIndex*3+2] = (float)patternIDL; 	
	
		rotbufIndex++;
		if(rotbufIndex==RSDMuStat->windowSize-1)
			rotbufIndex = 0;
		
		// Window center (bp)
		windowCenter = (sitePositionF + sitePositionL) / 2.0f;
		windowStart = sitePositionF;
		windowEnd = sitePositionL;

		// Var
		muVar = sitePositionL - sitePositionF;
		muVar /= RSDDataset->setRegionLength;
		muVar /= RSDMuStat->windowSize;		

		pcntexll = patternsLeftExclusive*snpsLeftExclusive;
		pcntexlr = patternsRightExclusive*snpsRightExclusive;

		// Mu_Sfs
		muSfs = dCnt1+dCntN==0?0.000001f:(((float)dCnt1)+((float)dCntN));		
		muSfs /= (float)RSDMuStat->windowSize;		
	
		// Mu_Ld
		muLd = pcntexll + pcntexlr==0?0.000001f:((((float)pcntexll)+((float)pcntexlr))/((float)(pcntl*pcntr)));
	
		// Mu
		mu =  muVar * muSfs * muLd;

		// Set for the next iteration
		prevDerivedAllele = derivedAlleleF; 
		firstLeftID_prev = firstLeftID; 
		firstRightID_prev = firstRightID;	

		
		// MuVar Max
		RSDMuStat->muVarMax = muVar;
		RSDMuStat->muVarMaxLoc = windowCenter;
		
		// MuSfs Max
		RSDMuStat->muSfsMax = muSfs;
		RSDMuStat->muSfsMaxLoc = windowCenter;
		
		// MuLd Max
		RSDMuStat->muLdMax = muLd;
		RSDMuStat->muLdMaxLoc = windowCenter;
		
		// Mu Max
		RSDMuStat->muMax = mu;
		RSDMuStat->muMaxLoc = windowCenter;
#ifdef _MUMEM
		muReportBufferIndex++;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+0] = windowCenter;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+1] = windowStart;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+2] = windowEnd;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+3] = muVar;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+4] = muSfs;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+5] = muLd;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+6] = mu;
#else		
		if(RSDCommandLine->fullReport==1)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	
		else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
#endif
	}
	
	// Scan chunk: all sliding window steps except for the first
	for(i=(int)RSDMuStat->windowSize;i<size;i++)
	{
		// mem access
		sitePositionL = (int)RSDChunk->chunkData[i*3]; 
		derivedAlleleL = (int)RSDChunk->chunkData[i*3+1];
		patternIDL = (int)RSDChunk->chunkData[i*3+2];

		// rotbuf access
		sitePositionF = (int)rotbuf[rotbufIndex*3]; 
		derivedAlleleF = (int)rotbuf[rotbufIndex*3+1];
		firstLeftID = (int)rotbuf[rotbufIndex*3+2]; 
		
		// we need to load the previous sliding step's right window's first pattern ID
		int rotbufIndex_temp = rotbufIndex + (DEFAULT_WINDOW_SIZE/2);

		if(rotbufIndex>=(DEFAULT_WINDOW_SIZE/2)-1)
			rotbufIndex_temp -= (DEFAULT_WINDOW_SIZE-1);

		firstRightID = (int)rotbuf[rotbufIndex_temp*3+2];

		// rotbuf update
		rotbuf[rotbufIndex*3] = (float)sitePositionL;		
		rotbuf[rotbufIndex*3+1] = (float)derivedAlleleL;
		rotbuf[rotbufIndex*3+2] = (float)patternIDL; 	
	
		rotbufIndex++;
		if(rotbufIndex==RSDMuStat->windowSize-1)
			rotbufIndex = 0;
		
		// Window center (bp)
		windowCenter = (sitePositionF + sitePositionL) / 2.0f;
		windowStart = sitePositionF;
		windowEnd = sitePositionL;

		// Var
		muVar = sitePositionL - sitePositionF;
		muVar /= RSDDataset->setRegionLength;
		muVar /= RSDMuStat->windowSize;		
		
		// SFS update
		dCnt1 -= (prevDerivedAllele==1);
		dCnt1 += (derivedAlleleL==1);

		dCntN -= (prevDerivedAllele==RSDDataset->setSamples-1);
		dCntN += (derivedAlleleL==RSDDataset->setSamples-1);

		// LD update
		int A = firstLeftID_prev;
		int B = firstRightID_prev;
		int C = patternIDL;

		int * LMem = patternMemLeft;
		int * RMem = patternMemRight;
		int * xLMem = patternMemLeftExclusive;
		int * xRMem = patternMemRightExclusive;
		
		int A_Le_ind=-1, A_Le_cnt=0;
		int B_Le_ind=-1, B_Le_cnt=0;
		int C_Le_ind=-1;//, C_Le_cnt=0;

		int A_Ri_ind=-1, A_Ri_cnt=0;
		int B_Ri_ind=-1, B_Ri_cnt=0;
		int C_Ri_ind=-1, C_Ri_cnt=0;

		int A_xLe_ind=-1, A_xLe_cnt=0;
		int B_xLe_ind=-1, B_xLe_cnt=0;
		int C_xLe_ind=-1, C_xLe_cnt=0;

		int A_xRi_ind=-1, A_xRi_cnt=0;
		int B_xRi_ind=-1, B_xRi_cnt=0;
		int C_xRi_ind=-1, C_xRi_cnt=0;

		int F_Le_vec[DEFAULT_WINDOW_SIZE/2], F_Le_ind=-1; // stacks
		int F_Ri_vec[DEFAULT_WINDOW_SIZE/2], F_Ri_ind=-1;
		int F_xLe_vec[DEFAULT_WINDOW_SIZE/2], F_xLe_ind=-1;
		int F_xRi_vec[DEFAULT_WINDOW_SIZE/2], F_xRi_ind=-1;

		for(j=0;j<DEFAULT_WINDOW_SIZE/2;j++)
		{
			// Mem read
			int curID_Le = LMem[j*2+0]; // ID
			int curSZ_Le = LMem[j*2+1]; // Count

			int curID_Ri = RMem[j*2+0];
			int curSZ_Ri = RMem[j*2+1];

			int curID_xLe = xLMem[j*2+0];
			int curSZ_xLe = xLMem[j*2+1];

			int curID_xRi = xRMem[j*2+0];
			int curSZ_xRi = xRMem[j*2+1];


			// Comp Le
			if(curID_Le==A)
			{	
				A_Le_ind = j;
				A_Le_cnt = curSZ_Le;
			}

			if(curID_Le==B)
			{	
				B_Le_ind = j;
				B_Le_cnt = curSZ_Le;
			}

			if(curID_Le==C)
			{	
				C_Le_ind = i;
				//C_Le_cnt = curSZ_Le;
			}

			if(curID_Le==-1)
			{
				F_Le_ind++;
				F_Le_vec[F_Le_ind] = j;
			}

			// Comp Ri
			if(curID_Ri==A)
			{	
				A_Ri_ind = j;
				A_Ri_cnt = curSZ_Ri;
			}

			if(curID_Ri==B)
			{	
				B_Ri_ind = j;
				B_Ri_cnt = curSZ_Ri;
			}

			if(curID_Ri==C)
			{	
				C_Ri_ind = j;
				C_Ri_cnt = curSZ_Ri;
			}

			if(curID_Ri==-1)
			{
				F_Ri_ind++;
				F_Ri_vec[F_Ri_ind] = j;
			}


			// Comp xLe
			if(curID_xLe==A)
			{	
				A_xLe_ind = j;
				A_xLe_cnt = curSZ_xLe;
			}

			if(curID_xLe==B)
			{	
				B_xLe_ind = j;
				B_xLe_cnt = curSZ_xLe;
			}

			if(curID_xLe==C)
			{	
				C_xLe_ind = j;
				C_xLe_cnt = curSZ_xLe;
			}

			if(curID_xLe==-1)
			{
				F_xLe_ind++;
				F_xLe_vec[F_xLe_ind] = j;
			}

			// Comp Ri
			if(curID_xRi==A)
			{	
				A_xRi_ind = j;
				A_xRi_cnt = curSZ_xRi;
			}

			if(curID_xRi==B)
			{	
				B_xRi_ind = j;
				B_xRi_cnt = curSZ_xRi;
			}

			if(curID_xRi==C)
			{	
				C_xRi_ind = j;
				C_xRi_cnt = curSZ_xRi;
			}

			if(curID_xRi==-1)
			{
				F_xRi_ind++;
				F_xRi_vec[F_xRi_ind] = j;
			}

		}

		assert(A_Le_ind!=-1);
		assert(B_Ri_ind!=-1);

		/************** Left **************/
		if(A!=B)
		{
			// Removing A
			if(A_Le_cnt==1)
			{
				assert(A_Le_ind>=0);

				LMem[A_Le_ind*2+0] = -1;
				LMem[A_Le_ind*2+1] = 0;

				F_Le_ind++;
				F_Le_vec[F_Le_ind] = A_Le_ind;
				
				A_Le_ind = -1;
				A_Le_cnt = 0;

				pcntl--;
			}
			else // A_Le_cnt!=1
			{
				A_Le_cnt--;
				LMem[A_Le_ind*2+1] = A_Le_cnt;				
			}

			// Adding B
			if(B_Le_ind==-1)
			{
				assert(F_Le_ind>=0);

				B_Le_ind = F_Le_vec[F_Le_ind];
				F_Le_ind--;

				B_Le_cnt = 1;

				assert(B_Le_ind>=0);
				
				LMem[B_Le_ind*2+0] = B;
				LMem[B_Le_ind*2+1] = B_Le_cnt;
				
				pcntl++;
			}
			else
			{
				B_Le_cnt++;
				LMem[B_Le_ind*2+1] = B_Le_cnt;
			}
		}
		// Checking C
		if(C_Le_ind!=-1 && C==A)
		{
			C_Le_ind = A_Le_ind;
			//C_Le_cnt = A_Le_cnt;
		}

		if(C==B)
		{
			C_Le_ind = B_Le_ind;
			//C_Le_cnt = B_Le_cnt;
		}
		/************** End Left **************/

		/************** Right **************/
		if(B!=C)
		{
			// Removing B
			if(B_Ri_cnt==1)
			{
				assert(B_Ri_ind>=0);

				RMem[B_Ri_ind*2+0] = -1;
				RMem[B_Ri_ind*2+1] = 0;

				F_Ri_ind++;
				F_Ri_vec[F_Ri_ind] = B_Ri_ind;
				
				B_Ri_ind = -1;
				B_Ri_cnt = 0;

				pcntr--;
			}
			else // B_Ri_cnt!=1
			{
				B_Ri_cnt--;
				RMem[B_Ri_ind*2+1] = B_Ri_cnt;				
			}

			// Adding C
			if(C_Ri_ind==-1)
			{
				assert(F_Ri_ind>=0);

				C_Ri_ind = F_Ri_vec[F_Ri_ind];
				F_Ri_ind--;

				C_Ri_cnt = 1;

				assert(C_Ri_ind>=0);
				
				RMem[C_Ri_ind*2+0] = C;
				RMem[C_Ri_ind*2+1] = C_Ri_cnt;

				pcntr++;
			}
			else
			{
				C_Ri_cnt++;
				RMem[C_Ri_ind*2+1] = C_Ri_cnt;
			}
		}
		// Checking A
		if(A_Ri_ind!=-1 && A==B)
		{
			A_Ri_ind = B_Ri_ind;
			A_Ri_cnt = B_Ri_cnt;
		}
		if(A==C)
		{
			A_Ri_ind = C_Ri_ind;
			A_Ri_cnt = C_Ri_cnt;
		}
		/************** End Right **************/

		/************** xLeft **************/			
		if(A!=B)
		{
			// Removing A 
			if(A_xLe_ind!=-1)
			{
				if(A_xLe_cnt==1) // remove A
				{
					assert(A_xLe_ind>=0);

					xLMem[A_xLe_ind*2+0] = -1;
					xLMem[A_xLe_ind*2+1] = 0;

					F_xLe_ind++;
					F_xLe_vec[F_xLe_ind] = A_xLe_ind;

					snpsLeftExclusive -= A_xLe_cnt; 
		
					A_xLe_ind = -1;
					A_xLe_cnt = 0;

					patternsLeftExclusive--;
				}
				else // A_xLe_cnt!=1
				{
					if(A_Ri_ind==-1)
					{
						A_xLe_cnt--;
						xLMem[A_xLe_ind*2+1] = A_xLe_cnt;
						snpsLeftExclusive--;
					}
					else
					{
						assert(A_xLe_ind>=0);

						xLMem[A_xLe_ind*2+0] = -1;
						xLMem[A_xLe_ind*2+1] = 0;

						F_xLe_ind++;
						F_xLe_vec[F_xLe_ind] = A_xLe_ind;

						snpsLeftExclusive -= A_xLe_cnt; 
		
						A_xLe_ind = -1;
						A_xLe_cnt = 0;

						patternsLeftExclusive--;
					}				
				}
			}
		}

		// Adding B
		if(B_Ri_ind==-1)
		{
			assert(F_xLe_ind>=0);

			B_xLe_ind = F_xLe_vec[F_xLe_ind];
			F_xLe_ind--;

			B_xLe_cnt = B_Le_cnt;

			assert(B_xLe_ind>=0);
		
			xLMem[B_xLe_ind*2+0] = B;
			xLMem[B_xLe_ind*2+1] = B_xLe_cnt;

			snpsLeftExclusive += B_xLe_cnt; 

			patternsLeftExclusive++;
		}

		if(C_xLe_ind!=-1 && A!=C) // Remove C
		{
			xLMem[C_xLe_ind*2+0] = -1;
			xLMem[C_xLe_ind*2+1] = 0;

			F_xLe_ind++;
			F_xLe_vec[F_xLe_ind] = C_xLe_ind;

			snpsLeftExclusive -= C_xLe_cnt; 

			C_xLe_ind = -1;
			C_xLe_cnt = 0;

			patternsLeftExclusive--;
		}
		/************** End xLeft **************/

		/************** xRight **************/
		if(B!=C)
		{
			// Removing B
			if(B_xRi_ind!=-1)
			{
				assert(B_xRi_ind>=0);

				xRMem[B_xRi_ind*2+0] = -1;
				xRMem[B_xRi_ind*2+1] = 0;

				F_xRi_ind++;
				F_xRi_vec[F_xRi_ind] = B_xRi_ind;

				snpsRightExclusive -= B_xRi_cnt;
			
				B_xRi_ind = -1;
				B_xRi_cnt = 0;

				patternsRightExclusive--;
			}
		}

		// Adding C
		if(C_Le_ind==-1)
		{
			if(C_xRi_ind==-1)
			{
				assert(F_xRi_ind>=0);

				C_xRi_ind = F_xRi_vec[F_xRi_ind];
				F_xRi_ind--;

				C_xRi_cnt = C_Ri_cnt;

				assert(C_xRi_ind>=0);
			
				xRMem[C_xRi_ind*2+0] = C;
				xRMem[C_xRi_ind*2+1] = C_xRi_cnt;

				patternsRightExclusive++;
				snpsRightExclusive += C_xRi_cnt;
			}
			else
			{
				xRMem[C_xRi_ind*2+1] = C_Ri_cnt; // C_Ri_cnt equals the previous value of xRMem[C_xRi_ind*2+1] + 1
				snpsRightExclusive++;
			}
			
		}
		else
		{
			if(C_xRi_ind!=-1) // removing C
			{
				assert(C_xRi_ind>=0);

				xRMem[C_xRi_ind*2+0] = -1;
				xRMem[C_xRi_ind*2+1] = 0;

				F_xRi_ind++;
				F_xRi_vec[F_xRi_ind] = C_xRi_ind;

				snpsRightExclusive -= C_xRi_cnt;
			
				C_xRi_ind = -1;
				C_xRi_cnt = 0;

				patternsRightExclusive--;
			}				
		}

		// Checking to add A
		if(A_Ri_ind!=-1 && A_Le_ind==-1 && A!=C)
		{
			assert(F_xRi_ind>=0);

			A_xRi_ind = F_xRi_vec[F_xRi_ind];
			F_xRi_ind--;

			A_xRi_cnt = A_Ri_cnt;

			assert(A_xRi_ind>=0);
		
			xRMem[A_xRi_ind*2+0] = A;
			xRMem[A_xRi_ind*2+1] = A_xRi_cnt;

			patternsRightExclusive++;
			snpsRightExclusive += A_xRi_cnt;
		}
		/************** End Right **************/		

		pcntexll = patternsLeftExclusive*snpsLeftExclusive;
		pcntexlr = patternsRightExclusive*snpsRightExclusive;

		// Mu_Sfs
		muSfs = dCnt1+dCntN==0?0.000001f:(((float)dCnt1)+((float)dCntN));		
		muSfs /= (float)RSDMuStat->windowSize;		
	
		// Mu_Ld
		muLd = pcntexll + pcntexlr==0?0.000001f:((((float)pcntexll)+((float)pcntexlr))/((float)(pcntl*pcntr)));
	
		// Mu
		mu =  muVar * muSfs * muLd;

		// Set for the next iteration
		prevDerivedAllele = derivedAlleleF; 
		firstLeftID_prev = firstLeftID; 
		firstRightID_prev = firstRightID;
		
		// MuVar Max
		RSDMuStat->muVarMaxLoc = muVar > RSDMuStat->muVarMax?windowCenter:RSDMuStat->muVarMaxLoc;
		RSDMuStat->muVarMax = muVar > RSDMuStat->muVarMax?muVar:RSDMuStat->muVarMax;

		// MuSfs Max
		RSDMuStat->muSfsMaxLoc = muSfs > RSDMuStat->muSfsMax?windowCenter:RSDMuStat->muSfsMaxLoc;
		RSDMuStat->muSfsMax = muSfs > RSDMuStat->muSfsMax?muSfs:RSDMuStat->muSfsMax;

		// MuLd Max
		RSDMuStat->muLdMaxLoc = muLd > RSDMuStat->muLdMax?windowCenter:RSDMuStat->muLdMaxLoc;
		RSDMuStat->muLdMax = muLd > RSDMuStat->muLdMax?muLd:RSDMuStat->muLdMax;

		// Mu Max
		RSDMuStat->muMaxLoc = mu > RSDMuStat->muMax?windowCenter:RSDMuStat->muMaxLoc;
		RSDMuStat->muMax = mu > RSDMuStat->muMax?mu:RSDMuStat->muMax;
		
#ifdef _MUMEM
		muReportBufferIndex++;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+0] = windowCenter;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+1] = windowStart;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+2] = windowEnd;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+3] = muVar;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+4] = muSfs;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+5] = muLd;
		RSDMuStat->muReportBuffer[muReportBufferIndex*7+6] = mu;
#else		
		if(RSDCommandLine->fullReport==1)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	
		else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
#endif
	}

#ifdef _MUMEM
	if(RSDCommandLine->fullReport==1)
	{
		for(i=0;i<=muReportBufferIndex;i++)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)RSDMuStat->muReportBuffer[i*7+0], 
												   (double)RSDMuStat->muReportBuffer[i*7+1], 
												   (double)RSDMuStat->muReportBuffer[i*7+2], 
												   (double)RSDMuStat->muReportBuffer[i*7+3], 
												   (double)RSDMuStat->muReportBuffer[i*7+4],
												   (double)RSDMuStat->muReportBuffer[i*7+5],
												   (double)RSDMuStat->muReportBuffer[i*7+6]);
	}
	else
	{
		for(i=0;i<=muReportBufferIndex;i++)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)RSDMuStat->muReportBuffer[i*7+0], 
								     (double)RSDMuStat->muReportBuffer[i*7+6]);
	}
#endif
}
// End of TRETS2019 MuStatistic implementation (DAER2+HAM+MLT)
#else
#ifdef _EXP1
// Start of experimental sliding-window implementation 1: one-sided expansion, single evaluation
void RSDMuStat_scanChunkBinary (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, j, size = (int)RSDChunk->chunkSize;

	double windowCenter = 0.0, windowStart = 0.0, windowEnd = 0.0;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	int64_t sfsSlack = RSDCommandLine->sfsSlack;

	int doExtra = 0; // min 0

	int64_t initWinSize = RSDMuStat->windowSize;

	for(i=0;i<size-RSDMuStat->windowSize+1-doExtra*2;i++)
	{
		for(int k=0;k<=doExtra*2;k=k+2)
		{
			RSDMuStat->windowSize = initWinSize + k;

			assert(RSDMuStat->windowSize%2==0);
			
			// SNP window range
			int snpf = i;
			int snpl = (int)(snpf + RSDMuStat->windowSize - 1);
			
			int winlsnpf = snpf;
			int winlsnpl = (int)(winlsnpf + RSDMuStat->windowSize/2 - 1);

			int winrsnpf = winlsnpl + 1;
			int winrsnpl = snpl;

			int totalSNPsL = winlsnpl - winlsnpf + 1;
			int totalSNPsR = winrsnpl - winrsnpf + 1;

			// Window center (bp)
			windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
			windowStart = RSDChunk->sitePosition[snpf];
			windowEnd = RSDChunk->sitePosition[snpl];
			/**/

			float isValid = 1.0;
	
			for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
				if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
					isValid = 0.0;

			// Mu_Var
			muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
			muVar /= RSDDataset->setRegionLength;
			muVar /= RSDMuStat->windowSize;
			//muVar *= RSDDataset->setSNPs;
			muVar *= RSDDataset->preLoadedsetSNPs; // v2.4

			// Mu_SFS
			int dCnt1 = 0, dCntN = 0;
			for(j=0;j<RSDMuStat->windowSize - 0;j++)
			{
				dCnt1 += (RSDChunk->derivedAlleleCount[i+j]<=sfsSlack);
				dCntN += (RSDChunk->derivedAlleleCount[i+j]>=RSDDataset->setSamples-sfsSlack);
			}

			float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
			muSfs = (float)dCnt1 + (float)dCntN*facN; 

			if(dCnt1+dCntN==0) 
				muSfs = 0.0000000001f;

			muSfs *= RSDDataset->muVarDenom; 
			muSfs /= (float)RSDMuStat->windowSize;

			// Mu_Ld
			int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;

			float tempTest = getPatternCounts (WIN_MODE_OSE, RSDMuStat, totalSNPsL, totalSNPsR, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
			muLd = tempTest;

			if(pcntexll + pcntexlr==0) 
			{
				muLd = 0.0000000001f;
			}
			else
			{
				muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
			}

			// Mu
			mu =  muVar * muSfs * muLd * isValid;

			// MuVar Max
			if (muVar > RSDMuStat->muVarMax)
			{
				RSDMuStat->muVarMax = muVar;
				RSDMuStat->muVarMaxLoc = windowCenter;
			}

			// MuSfs Max
			if (muSfs > RSDMuStat->muSfsMax)
			{
				RSDMuStat->muSfsMax = muSfs;
				RSDMuStat->muSfsMaxLoc = windowCenter;
			}

			// MuLd Max
			if (muLd > RSDMuStat->muLdMax)
			{
				RSDMuStat->muLdMax = muLd;
				RSDMuStat->muLdMaxLoc = windowCenter;
			}

			// Mu Max
			if (mu > RSDMuStat->muMax)
			{
				RSDMuStat->muMax = mu;
				RSDMuStat->muMaxLoc = windowCenter;
			}	
		}

		RSDMuStat->windowSize = initWinSize ;

		if(RSDCommandLine->fullReport==1)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
	}
}
// End of experimental sliding-window implementation 1: one-sided expansion
#endif

#ifdef _EXP2
// Start of experimental sliding-window implementation 2: floating window center, single evaluation
void RSDMuStat_scanChunkBinary (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, j, size = (int)RSDChunk->chunkSize;

	double windowCenter = 0.0, windowStart = 0.0, windowEnd = 0.0;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	int64_t sfsSlack = RSDCommandLine->sfsSlack;

	int slackVal = 5;

	for(i=0;i<size-RSDMuStat->windowSize+1;i++)
	{
		// SNP window range
		int snpf = i;
		int snpl = (int)(snpf + RSDMuStat->windowSize - 1);

		int curMiddleFxd = (int)(snpf + RSDMuStat->windowSize/2);		

		for(int curMiddle=curMiddleFxd-slackVal;curMiddle<=curMiddleFxd+slackVal;curMiddle++) 
		{
			int winlsnpf = snpf;
			int winlsnpl = curMiddle-1;

			int winrsnpf = curMiddle;
			int winrsnpl = snpl;
		
			int totalSNPsL = winlsnpl - winlsnpf + 1;
			int totalSNPsR = winrsnpl - winrsnpf + 1;

			// Window center (bp)
			windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
			windowStart = RSDChunk->sitePosition[snpf];
			windowEnd = RSDChunk->sitePosition[snpl];

			float isValid = 1.0;

			for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
				if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
					isValid = 0.0;

			// Mu_Var
			muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
			muVar /= RSDDataset->setRegionLength;
			muVar /= RSDMuStat->windowSize;
			//muVar *= RSDDataset->setSNPs;
			muVar *= RSDDataset->preLoadedsetSNPs; // v2.4

			// Mu_SFS
			int dCnt1 = 0, dCntN = 0;
			for(j=0;j<RSDMuStat->windowSize - 0;j++)
			{
				dCnt1 += (RSDChunk->derivedAlleleCount[i+j]<=sfsSlack);
				dCntN += (RSDChunk->derivedAlleleCount[i+j]>=RSDDataset->setSamples-sfsSlack);
			}

			float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
			muSfs = (float)dCnt1 + (float)dCntN*facN; 

			if(dCnt1+dCntN==0) 
				muSfs = 0.0000000001f;

			muSfs *= RSDDataset->muVarDenom; 
			muSfs /= (float)RSDMuStat->windowSize;

			// Mu_Ld
			int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;

			float tempTest = getPatternCounts (WIN_MODE_FC, RSDMuStat, totalSNPsL, totalSNPsR, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
			muLd = tempTest;

			if(pcntexll + pcntexlr==0) 
			{
				muLd = 0.0000000001f;
			}
			else
			{
				muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
			}

			// Mu
			mu =  muVar * muSfs * muLd * isValid;

			// MuVar Max
			if (muVar > RSDMuStat->muVarMax)
			{
				RSDMuStat->muVarMax = muVar;
				RSDMuStat->muVarMaxLoc = windowCenter;
			}

			// MuSfs Max
			if (muSfs > RSDMuStat->muSfsMax)
			{
				RSDMuStat->muSfsMax = muSfs;
				RSDMuStat->muSfsMaxLoc = windowCenter;
			}

			// MuLd Max
			if (muLd > RSDMuStat->muLdMax)
			{
				RSDMuStat->muLdMax = muLd;
				RSDMuStat->muLdMaxLoc = windowCenter;
			}

			// Mu Max
			if (mu > RSDMuStat->muMax)
			{
				RSDMuStat->muMax = mu;
				RSDMuStat->muMaxLoc = windowCenter;
			}				
		}		

		if(RSDCommandLine->fullReport==1)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
	}
}
// End of experimental sliding-window implementation 2: floating window center
#endif

#ifdef _EXP3
// Start of experimental sliding-window implementation 3: double-sided reduction, all-to-all evaluation
void RSDMuStat_scanChunkBinary (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, j, size = (int)RSDChunk->chunkSize;

	double windowCenter = 0.0, windowStart = 0.0, windowEnd = 0.0;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	int64_t sfsSlack = RSDCommandLine->sfsSlack;

	int doSearch = 10; // min 1

	for(i=0;i<size-RSDMuStat->windowSize+1;i++)
	{
		int hSize = (int)(RSDMuStat->windowSize/2);
		
		// SNP window range
		int snpf = i;
		int snpl = (int)(snpf + RSDMuStat->windowSize - 1);	

		int curMiddle = (int)(snpf + RSDMuStat->windowSize/2);		

		for(int m=snpf;m<=curMiddle-(hSize-doSearch+1);m++)
		{
			for(int n=curMiddle+(hSize-doSearch);n<=snpl;n++)
			{
				int winlsnpf = m;
				int winlsnpl = curMiddle-1;

				int winrsnpf = curMiddle;
				int winrsnpl = n;

				int totalSNPsL = winlsnpl - winlsnpf + 1;
				int totalSNPsR = winrsnpl - winrsnpf + 1;

				// Window center (bp)
				windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
				windowStart = RSDChunk->sitePosition[snpf];
				windowEnd = RSDChunk->sitePosition[snpl];

				float isValid = 1.0;

				for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
					if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
						isValid = 0.0;

				// Mu_Var
				muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
				muVar /= RSDDataset->setRegionLength;
				muVar /= RSDMuStat->windowSize;
				//muVar *= RSDDataset->setSNPs;
				muVar *= RSDDataset->preLoadedsetSNPs; // v2.4	
		
				// Mu_SFS
				int dCnt1 = 0, dCntN = 0;
				for(j=0;j<RSDMuStat->windowSize - 0;j++)
				{
					dCnt1 += (RSDChunk->derivedAlleleCount[i+j]<=sfsSlack);
					dCntN += (RSDChunk->derivedAlleleCount[i+j]>=RSDDataset->setSamples-sfsSlack);
				}

				float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
				muSfs = (float)dCnt1 + (float)dCntN*facN; 

				if(dCnt1+dCntN==0) 
					muSfs = 0.0000000001f;

				muSfs *= RSDDataset->muVarDenom; 
				muSfs /= (float)RSDMuStat->windowSize;
		
				// Mu_Ld
				int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;

				float tempTest = getPatternCounts (WIN_MODE_DSR, RSDMuStat, totalSNPsL, totalSNPsR, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
				muLd = tempTest;

				if(pcntexll + pcntexlr==0) 
				{
					muLd = 0.0000000001f;
				}
				else
				{
					muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
				}	

				// Mu
				mu =  muVar * muSfs * muLd * isValid;

				// MuVar Max
				if (muVar > RSDMuStat->muVarMax)
				{
					RSDMuStat->muVarMax = muVar;
					RSDMuStat->muVarMaxLoc = windowCenter;
				}

				// MuSfs Max
				if (muSfs > RSDMuStat->muSfsMax)
				{
					RSDMuStat->muSfsMax = muSfs;
					RSDMuStat->muSfsMaxLoc = windowCenter;
				}

				// MuLd Max
				if (muLd > RSDMuStat->muLdMax)
				{
					RSDMuStat->muLdMax = muLd;
					RSDMuStat->muLdMaxLoc = windowCenter;
				}

				// Mu Max
				if (mu > RSDMuStat->muMax)
				{
					RSDMuStat->muMax = mu;
					RSDMuStat->muMaxLoc = windowCenter;
				}				
			}
		}			
	
		if(RSDCommandLine->fullReport==1)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
	}
}
// End of experimental sliding-window implementation 3: double-sided expansion, all-to-all evaluation
#endif

#ifdef _EXP4
// Start of experimental sliding-window implementation 4: all three
void RSDMuStat_scanChunkBinary (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, j, size = (int)RSDChunk->chunkSize;

	double windowCenter = 0.0, windowStart = 0.0, windowEnd = 0.0;

	float muVar = 0.0f, muVarM1 = 0.0f, muVarM2 = 0.0f, muVarM3 = 0.0f;
	float muSfs = 0.0f, muSfsM1 = 0.0f, muSfsM2 = 0.0f, muSfsM3 = 0.0f;
	float muLd = 0.0f, muLdM1 = 0.0f, muLdM2 = 0.0f, muLdM3 = 0.0f;
	float mu = 0.0f, muM1 = 0.0f, muM2 = 0.0f, muM3 = 0.0f;

	double windowCenter_muM1 = 0.0, windowStart_muM1 = 0.0, windowEnd_muM1 = 0.0;
	double windowCenter_muVarM1 = 0.0;//, windowStart_muVarM1 = 0.0, windowEnd_muVarM1 = 0.0;
	double windowCenter_muSfsM1 = 0.0;//, windowStart_muSfsM1 = 0.0, windowEnd_muSfsM1 = 0.0;
	double windowCenter_muLdM1 = 0.0;//, windowStart_muLdM1 = 0.0, windowEnd_muLdM1 = 0.0;

	double windowCenter_muM2 = 0.0, windowStart_muM2 = 0.0, windowEnd_muM2 = 0.0;
	double windowCenter_muVarM2 = 0.0;//, windowStart_muVarM2 = 0.0, windowEnd_muVarM2 = 0.0;
	double windowCenter_muSfsM2 = 0.0;//, windowStart_muSfsM2 = 0.0, windowEnd_muSfsM2 = 0.0;
	double windowCenter_muLdM2 = 0.0;//, windowStart_muLdM2 = 0.0, windowEnd_muLdM2 = 0.0;

	double windowCenter_muM3 = 0.0, windowStart_muM3 = 0.0, windowEnd_muM3 = 0.0;
	double windowCenter_muVarM3 = 0.0;//, windowStart_muVarM3 = 0.0, windowEnd_muVarM3 = 0.0;
	double windowCenter_muSfsM3 = 0.0;//, windowStart_muSfsM3 = 0.0, windowEnd_muSfsM3 = 0.0;
	double windowCenter_muLdM3 = 0.0;//, windowStart_muLdM3 = 0.0, windowEnd_muLdM3 = 0.0;


	int64_t sfsSlack = RSDCommandLine->sfsSlack;
		
	int64_t initWinSize = RSDMuStat->windowSize;

	int doExtra = (int)(initWinSize/2); // min 0
	int slackVal = (int)(initWinSize/2-1);
	int doSearch = (int)(initWinSize/2); // min 1	

	for(i=0;i<size-RSDMuStat->windowSize+1-doExtra*2;i++)
	{
		// mode 1
		{
			muVarM1 = 0.0f;
			muSfsM1 = 0.0f;
			muLdM1 = 0.0f;
  			muM1 = 0.0f;

			for(int k=0;k<=doExtra*2;k=k+2)
			{
				RSDMuStat->windowSize = initWinSize + k;

				assert(RSDMuStat->windowSize%2==0);
			
				// SNP window range
				int snpf = i;
				int snpl = (int)(snpf + RSDMuStat->windowSize - 1);
			
				int winlsnpf = snpf;
				int winlsnpl = (int)(winlsnpf + RSDMuStat->windowSize/2 - 1);

				int winrsnpf = winlsnpl + 1;
				int winrsnpl = snpl;

				int totalSNPsL = winlsnpl - winlsnpf + 1;
				int totalSNPsR = winrsnpl - winrsnpf + 1;

				// Window center (bp)
				windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
				windowStart = RSDChunk->sitePosition[snpf];
				windowEnd = RSDChunk->sitePosition[snpl];

				float isValid = 1.0;
	
				for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
					if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
						isValid = 0.0;

				// Mu_Var
				muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
				muVar /= RSDDataset->setRegionLength;
				muVar /= RSDMuStat->windowSize;
				//muVar *= RSDDataset->setSNPs;
				muVar *= RSDDataset->preLoadedsetSNPs; // v2.4

				// Mu_SFS
				int dCnt1 = 0, dCntN = 0;
				for(j=0;j<RSDMuStat->windowSize - 0;j++)
				{
					dCnt1 += (RSDChunk->derivedAlleleCount[i+j]<=sfsSlack);
					dCntN += (RSDChunk->derivedAlleleCount[i+j]>=RSDDataset->setSamples-sfsSlack);
				}

				float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
				muSfs = (float)dCnt1 + (float)dCntN*facN; 

				if(dCnt1+dCntN==0) 
					muSfs = 0.0000000001f;

				muSfs *= RSDDataset->muVarDenom; 
				muSfs /= (float)RSDMuStat->windowSize;

				// Mu_Ld
				int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;

				float tempTest = getPatternCounts (WIN_MODE_OSE, RSDMuStat, totalSNPsL, totalSNPsR, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
				muLd = tempTest;

				if(pcntexll + pcntexlr==0) 
				{
					muLd = 0.0000000001f;
				}
				else
				{
					muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
				}

				// Mu
				mu =  muVar * muSfs * muLd * isValid;

				if(mu > muM1)
				{
					muM1 = mu;
					windowCenter_muM1 = windowCenter;
					windowStart_muM1 = windowStart;
					windowEnd_muM1 = windowEnd;					
				}

				if(muVar > muVarM1)
				{
					muVarM1 = muVar;
					windowCenter_muVarM1 = windowCenter;
					//windowStart_muVarM1 = windowStart;
					//windowEnd_muVarM1 = windowEnd;					
				}

				if(muSfs > muSfsM1)
				{
					muSfsM1 = muSfs;
					windowCenter_muSfsM1 = windowCenter;
					//windowStart_muSfsM1 = windowStart;
					//windowEnd_muSfsM1 = windowEnd;					
				}

				if(muLd > muLdM1)
				{
					muLdM1 = muLd;
					windowCenter_muLdM1 = windowCenter;
					//windowStart_muLdM1 = windowStart;
					//windowEnd_muLdM1 = windowEnd;					
				}				
			}

			RSDMuStat->windowSize = initWinSize ;
		}

		// mode 2
		{
			muVarM2 = 0.0f;
			muSfsM2 = 0.0f;
			muLdM2 = 0.0f;
  			muM2 = 0.0f;

			// SNP window range
			int snpf = i;
			int snpl = (int)(snpf + RSDMuStat->windowSize - 1);

			int curMiddleFxd = (int)(snpf + RSDMuStat->windowSize/2);		

			for(int curMiddle=curMiddleFxd-slackVal;curMiddle<=curMiddleFxd+slackVal;curMiddle++) 
			{
				int winlsnpf = snpf;
				int winlsnpl = curMiddle-1;

				int winrsnpf = curMiddle;
				int winrsnpl = snpl;
		
				int totalSNPsL = winlsnpl - winlsnpf + 1;
				int totalSNPsR = winrsnpl - winrsnpf + 1;

				// Window center (bp)
				windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
				windowStart = RSDChunk->sitePosition[snpf];
				windowEnd = RSDChunk->sitePosition[snpl];

				float isValid = 1.0;

				for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
					if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
						isValid = 0.0;

				// Mu_Var
				muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
				muVar /= RSDDataset->setRegionLength;
				muVar /= RSDMuStat->windowSize;
				//muVar *= RSDDataset->setSNPs;
				muVar *= RSDDataset->preLoadedsetSNPs; // v2.4

				// Mu_SFS
				int dCnt1 = 0, dCntN = 0;
				for(j=0;j<RSDMuStat->windowSize - 0;j++)
				{
					dCnt1 += (RSDChunk->derivedAlleleCount[i+j]<=sfsSlack);
					dCntN += (RSDChunk->derivedAlleleCount[i+j]>=RSDDataset->setSamples-sfsSlack);
				}

				float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
				muSfs = (float)dCnt1 + (float)dCntN*facN; 

				if(dCnt1+dCntN==0) 
					muSfs = 0.0000000001f;

				muSfs *= RSDDataset->muVarDenom; 
				muSfs /= (float)RSDMuStat->windowSize;

				// Mu_Ld
				int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;

				float tempTest = getPatternCounts (WIN_MODE_FC, RSDMuStat, totalSNPsL, totalSNPsR, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
				muLd = tempTest;

				if(pcntexll + pcntexlr==0) 
				{
					muLd = 0.0000000001f;
				}
				else
				{
					muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
				}

				// Mu
				mu =  muVar * muSfs * muLd * isValid;

				if(mu > muM2)
				{
					muM2 = mu;
					windowCenter_muM2 = windowCenter;
					windowStart_muM2 = windowStart;
					windowEnd_muM2 = windowEnd;					
				}

				if(muVar > muVarM2)
				{
					muVarM2 = muVar;
					windowCenter_muVarM2 = windowCenter;
					//windowStart_muVarM2 = windowStart;
					//windowEnd_muVarM2 = windowEnd;					
				}

				if(muSfs > muSfsM2)
				{
					muSfsM2 = muSfs;
					windowCenter_muSfsM2 = windowCenter;
					//windowStart_muSfsM2 = windowStart;
					//windowEnd_muSfsM2 = windowEnd;					
				}

				if(muLd > muLdM2)
				{
					muLdM2 = muLd;
					windowCenter_muLdM2 = windowCenter;
					//windowStart_muLdM2 = windowStart;
					//windowEnd_muLdM2 = windowEnd;					
				}					
			}
		}

		// Mode 3
		{
			muVarM3 = 0.0f;
			muSfsM3 = 0.0f;
			muLdM3 = 0.0f;
  			muM3 = 0.0f;

			int hSize = (int)(RSDMuStat->windowSize/2);
		
			// SNP window range
			int snpf = i;
			int snpl = (int)(snpf + RSDMuStat->windowSize - 1);	

			int curMiddle = (int)(snpf + RSDMuStat->windowSize/2);		

			for(int m=snpf;m<=curMiddle-(hSize-doSearch+1);m++)
			{
				for(int n=curMiddle+(hSize-doSearch);n<=snpl;n++)
				{
					int winlsnpf = m;
					int winlsnpl = curMiddle-1;

					int winrsnpf = curMiddle;
					int winrsnpl = n;

					int totalSNPsL = winlsnpl - winlsnpf + 1;
					int totalSNPsR = winrsnpl - winrsnpf + 1;

					// Window center (bp)
					windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
					windowStart = RSDChunk->sitePosition[snpf];
					windowEnd = RSDChunk->sitePosition[snpl];

					float isValid = 1.0;

					for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
						if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
							isValid = 0.0;

					// Mu_Var
					muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
					muVar /= RSDDataset->setRegionLength;
					muVar /= RSDMuStat->windowSize;
					//muVar *= RSDDataset->setSNPs;
					muVar *= RSDDataset->preLoadedsetSNPs; // v2.4	
		
					// Mu_SFS
					int dCnt1 = 0, dCntN = 0;
					for(j=0;j<RSDMuStat->windowSize - 0;j++)
					{
						dCnt1 += (RSDChunk->derivedAlleleCount[i+j]<=sfsSlack);
						dCntN += (RSDChunk->derivedAlleleCount[i+j]>=RSDDataset->setSamples-sfsSlack);
					}

					float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
					muSfs = (float)dCnt1 + (float)dCntN*facN; 

					if(dCnt1+dCntN==0) 
						muSfs = 0.0000000001f;

					muSfs *= RSDDataset->muVarDenom; 
					muSfs /= (float)RSDMuStat->windowSize;
		
					// Mu_Ld
					int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;

					float tempTest = getPatternCounts (WIN_MODE_DSR, RSDMuStat, totalSNPsL, totalSNPsR, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
					muLd = tempTest;

					if(pcntexll + pcntexlr==0) 
					{
						muLd = 0.0000000001f;
					}
					else
					{
						muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
					}	

					// Mu
					mu =  muVar * muSfs * muLd * isValid;

					if(mu > muM3)
					{
						muM3 = mu;
						windowCenter_muM3 = windowCenter;
						windowStart_muM3 = windowStart;
						windowEnd_muM3 = windowEnd;					
					}

					if(muVar > muVarM3)
					{
						muVarM3 = muVar;
						windowCenter_muVarM3 = windowCenter;
						//windowStart_muVarM3 = windowStart;
						//windowEnd_muVarM3 = windowEnd;
					
					}

					if(muSfs > muSfsM3)
					{
						muSfsM3 = muSfs;
						windowCenter_muSfsM3 = windowCenter;
						//windowStart_muSfsM3 = windowStart;
						//windowEnd_muSfsM3 = windowEnd;					
					}

					if(muLd > muLdM3)
					{
						muLdM3 = muLd;
						windowCenter_muLdM3 = windowCenter;
						//windowStart_muLdM3 = windowStart;
						//windowEnd_muLdM3 = windowEnd;					
					}				
				}
			}
		}

		// MuVar Max over M1, M2, M3
		if(muVarM1 > muVarM2)
		{
			muVar = muVarM1;
			windowCenter = windowCenter_muVarM1;
			//windowStart = windowStart_muVarM1;
			//windowEnd = windowEnd_muVarM1;
		}
		else
		{
			muVar = muVarM2;
			windowCenter = windowCenter_muVarM2;
			//windowStart = windowStart_muVarM2;
			//windowEnd = windowEnd_muVarM2;
		}

		if(muVarM3 > muVar)
		{
			muVar = muVarM3;
			windowCenter = windowCenter_muVarM3;
			//windowStart = windowStart_muVarM3;
			//windowEnd = windowEnd_muVarM3;
		}

		// MuVar Max
		if (muVar > RSDMuStat->muVarMax)
		{
			RSDMuStat->muVarMax = muVar;
			RSDMuStat->muVarMaxLoc = windowCenter;
		}


		// MuSfs Max over M1, M2, M3
		if(muSfsM1 > muSfsM2)
		{
			muSfs = muSfsM1;
			windowCenter = windowCenter_muSfsM1;
			//windowStart = windowStart_muVarM1;
			//windowEnd = windowEnd_muVarM1;
		}
		else
		{
			muSfs = muSfsM2;
			windowCenter = windowCenter_muSfsM2;
			//windowStart = windowStart_muVarM2;
			//windowEnd = windowEnd_muVarM2;
		}

		if(muSfsM3 > muSfs)
		{
			muSfs = muSfsM3;
			windowCenter = windowCenter_muSfsM3;
			//windowStart = windowStart_muVarM3;
			//windowEnd = windowEnd_muVarM3;
		}		

		// MuSfs Max
		if (muSfs > RSDMuStat->muSfsMax)
		{
			RSDMuStat->muSfsMax = muSfs;
			RSDMuStat->muSfsMaxLoc = windowCenter;
		}

		// MuLd Max over M1, M2, M3
		if(muLdM1 > muLdM2)
		{
			muLd = muLdM1;
			windowCenter = windowCenter_muLdM1;
			//windowStart = windowStart_muVarM1;
			//windowEnd = windowEnd_muVarM1;
		}
		else
		{
			muLd = muLdM2;
			windowCenter = windowCenter_muLdM2;
			//windowStart = windowStart_muVarM2;
			//windowEnd = windowEnd_muVarM2;
		}

		if(muLdM3 > muLd)
		{
			muLd = muLdM3;
			windowCenter = windowCenter_muLdM3;
			//windowStart = windowStart_muVarM3;
			//windowEnd = windowEnd_muVarM3;
		}	

		// MuLd Max
		if (muLd > RSDMuStat->muLdMax)
		{
			RSDMuStat->muLdMax = muLd;
			RSDMuStat->muLdMaxLoc = windowCenter;
		}

		// Mu Max over M1, M2, M3
		if(muM1 > muM2)
		{
			mu = muM1;
			windowCenter = windowCenter_muM1;
			windowStart = windowStart_muM1;
			windowEnd = windowEnd_muM1;
		}
		else
		{
			mu = muM2;
			windowCenter = windowCenter_muM2;
			windowStart = windowStart_muM2;
			windowEnd = windowEnd_muM2;
		}

		if(muM3 > mu)
		{
			mu = muM3;
			windowCenter = windowCenter_muM3;
			windowStart = windowStart_muM3;
			windowEnd = windowEnd_muM3;
		}

		// Mu Max
		if (mu > RSDMuStat->muMax)
		{
			RSDMuStat->muMax = mu;
			RSDMuStat->muMaxLoc = windowCenter;
		}				

		if(RSDCommandLine->fullReport==1)
		{
			//assert(0);
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	}
		else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
	}
}
// End of experimental sliding-window implementation 4: all three
#endif

void RSDMuStat_output2FileSimple (RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu)
{
	windowStart = windowStart+0.0; // make compilers happy
	windowEnd = windowEnd+0.0;
	muVar = muVar+0.0;
	muSfs = muSfs+0.0;
	muLd = muLd+0.0;

	fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", windowCenter, mu);	
}

void RSDMuStat_output2FileFull (RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu)
{
	fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", windowCenter, windowStart, windowEnd, muVar, muSfs, muLd, mu);
}

void RSDMuStat_output2BufferSimple (RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu)
{
	windowStart = windowStart+0.0; // make compilers happy
	windowEnd = windowEnd+0.0;
	muVar = muVar+0.0;
	muSfs = muSfs+0.0;
	muLd = muLd+0.0;

	if(RSDMuStat->currentScoreIndex+1>=RSDMuStat->bufferMemMaxSize)
	{
		RSDMuStat->bufferMemMaxSize+=SCOREBUFFER_REALLOC_INCR;

		RSDMuStat->buffer0Data = rsd_realloc(RSDMuStat->buffer0Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer0Data!=NULL);

		RSDMuStat->buffer6Data = rsd_realloc(RSDMuStat->buffer6Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer6Data!=NULL);
	}

	RSDMuStat->buffer0Data[RSDMuStat->currentScoreIndex]=windowCenter;
	RSDMuStat->buffer6Data[RSDMuStat->currentScoreIndex]=mu;	
}

void RSDMuStat_output2BufferFull (RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu)
{
	if(RSDMuStat->currentScoreIndex+1>=RSDMuStat->bufferMemMaxSize)
	{
		RSDMuStat->bufferMemMaxSize+=SCOREBUFFER_REALLOC_INCR;

		RSDMuStat->buffer0Data = rsd_realloc(RSDMuStat->buffer0Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer0Data!=NULL);

		RSDMuStat->buffer1Data = rsd_realloc(RSDMuStat->buffer1Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer1Data!=NULL);

		RSDMuStat->buffer2Data = rsd_realloc(RSDMuStat->buffer2Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer2Data!=NULL);

		RSDMuStat->buffer3Data = rsd_realloc(RSDMuStat->buffer3Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer3Data!=NULL);

		RSDMuStat->buffer4Data = rsd_realloc(RSDMuStat->buffer4Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer4Data!=NULL);

		RSDMuStat->buffer5Data = rsd_realloc(RSDMuStat->buffer5Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer5Data!=NULL);

		RSDMuStat->buffer6Data = rsd_realloc(RSDMuStat->buffer6Data, sizeof(double)*(unsigned long)RSDMuStat->bufferMemMaxSize);
		assert(RSDMuStat->buffer6Data!=NULL);
	}

	RSDMuStat->buffer0Data[RSDMuStat->currentScoreIndex]=windowCenter;
	RSDMuStat->buffer1Data[RSDMuStat->currentScoreIndex]=windowStart;
	RSDMuStat->buffer2Data[RSDMuStat->currentScoreIndex]=windowEnd;
	RSDMuStat->buffer3Data[RSDMuStat->currentScoreIndex]=muVar;
	RSDMuStat->buffer4Data[RSDMuStat->currentScoreIndex]=muSfs;
	RSDMuStat->buffer5Data[RSDMuStat->currentScoreIndex]=muLd;
	RSDMuStat->buffer6Data[RSDMuStat->currentScoreIndex]=mu;
}

inline void RSDMuStat_storeOutputConfigure (RSDCommandLine_t * RSDCommandLine)
{
	if(RSDCommandLine->gridSize==-1)
	{
		if(RSDCommandLine->fullReport==1)
		{
			RSDMuStat_storeOutput = &RSDMuStat_output2FileFull;
		}
		else
		{
			RSDMuStat_storeOutput = &RSDMuStat_output2FileSimple;
		}
	}
	else
	{
		if(RSDCommandLine->fullReport==1)
		{
			RSDMuStat_storeOutput = &RSDMuStat_output2BufferFull;
		}
		else
		{
			RSDMuStat_storeOutput = &RSDMuStat_output2BufferSimple;
		}
	}
}

void RSDMuStat_writeBuffer2File (RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine)
{
	if(RSDCommandLine->gridSize==-1)
		return;

	int i = 0;

	double * windowCenter = (double*)rsd_malloc(sizeof(double)*(unsigned long)RSDCommandLine->gridSize);
	assert(windowCenter!=NULL);

	double firstPos = RSDMuStat->buffer0Data[0];
	double lastPos = RSDMuStat->buffer0Data[RSDMuStat->currentScoreIndex];
	double step = (double)(lastPos - firstPos)/(RSDCommandLine->gridSize-1.0);

	for(i=0;i<RSDCommandLine->gridSize-1;i++)
	{
		windowCenter[i] = firstPos + i*step;
	}
	windowCenter[RSDCommandLine->gridSize-1] = lastPos;

	//gsl_interp_accel * acc = gsl_interp_accel_alloc ();
    	//gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, RSDMuStat->currentScoreIndex+1);
	//gsl_spline_init (spline, RSDMuStat->buffer0Data, RSDMuStat->buffer6Data, RSDMuStat->currentScoreIndex+1);
	//val = gsl_spline_eval (spline, windowCenter, acc);	

	double * mu = (double*)rsd_malloc(sizeof(double)*(unsigned long)RSDCommandLine->gridSize);
	assert(mu!=NULL);

   	gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
	gsl_interp * interpolation = gsl_interp_alloc (gsl_interp_linear,(size_t)RSDMuStat->currentScoreIndex+1);
      	gsl_interp_init(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer6Data, (size_t)RSDMuStat->currentScoreIndex+1);

	for(i=0;i<RSDCommandLine->gridSize;i++)
	{
		mu[i] =  gsl_interp_eval(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer6Data, windowCenter[i], accelerator);
	}
     	gsl_interp_free (interpolation);
	gsl_interp_accel_free(accelerator);

	if(RSDCommandLine->fullReport==1)
	{
		double * muVAR = (double*)rsd_malloc(sizeof(double)*(unsigned long)RSDCommandLine->gridSize);
		assert(muVAR!=NULL);

		accelerator =  gsl_interp_accel_alloc();
		interpolation = gsl_interp_alloc (gsl_interp_linear,(size_t)RSDMuStat->currentScoreIndex+1);
		gsl_interp_init(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer3Data, (size_t)RSDMuStat->currentScoreIndex+1);

		for(i=0;i<RSDCommandLine->gridSize;i++)
		{
			muVAR[i] =  gsl_interp_eval(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer3Data, windowCenter[i], accelerator);
		}
		gsl_interp_free (interpolation);
		gsl_interp_accel_free(accelerator);

		double * muSFS = (double*)rsd_malloc(sizeof(double)*(unsigned long)RSDCommandLine->gridSize);
		assert(muSFS!=NULL);

		accelerator =  gsl_interp_accel_alloc();
		interpolation = gsl_interp_alloc (gsl_interp_linear,(size_t)RSDMuStat->currentScoreIndex+1);
		gsl_interp_init(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer4Data, (size_t)RSDMuStat->currentScoreIndex+1);

		for(i=0;i<RSDCommandLine->gridSize;i++)
		{
			muSFS[i] =  gsl_interp_eval(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer4Data, windowCenter[i], accelerator);
		}
		gsl_interp_free (interpolation);
		gsl_interp_accel_free(accelerator);

		double * muLD = (double*)rsd_malloc(sizeof(double)*(unsigned long)RSDCommandLine->gridSize);
		assert(muLD!=NULL);

		accelerator =  gsl_interp_accel_alloc();
		interpolation = gsl_interp_alloc (gsl_interp_linear,(size_t)RSDMuStat->currentScoreIndex+1);
		gsl_interp_init(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer5Data, (size_t)RSDMuStat->currentScoreIndex+1);

		for(i=0;i<RSDCommandLine->gridSize;i++)
		{
			muLD[i] =  gsl_interp_eval(interpolation, RSDMuStat->buffer0Data, RSDMuStat->buffer5Data, windowCenter[i], accelerator);
		}
		gsl_interp_free (interpolation);
		gsl_interp_accel_free(accelerator);

		for(i=1;i<RSDCommandLine->gridSize;i++)
		{
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", windowCenter[i], 
												   0.0, 
												   0.0, 	
												   muVAR[i],
												   muSFS[i],
												   muLD[i],
												   mu[i]);

		}

		free(muVAR);
		free(muSFS);
		free(muLD);
	}
	else
	{
		for(i=1;i<RSDCommandLine->gridSize;i++)
		{
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", windowCenter[i], mu[i]);

		}
	}     
	
	free(windowCenter);	
	free(mu);

	fflush(RSDMuStat->reportFP);
}

#ifdef _REF
// Default SW implementation
void RSDMuStat_scanChunkBinary (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, j, size = (int)RSDChunk->chunkSize;

	double windowCenter = 0.0, windowStart = 0.0, windowEnd = 0.0;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	int64_t sfsSlack = RSDCommandLine->sfsSlack;
	
	RSDMuStat_storeOutputConfigure(RSDCommandLine);

	// Mu_SFS initialization
	int dCnt1 = 0, dCntN = 0;
	for(i=0;i<RSDMuStat->windowSize - 0;i++)
	{
		dCnt1 += (RSDChunk->derivedAlleleCount[i]<=sfsSlack);
		dCntN += (RSDChunk->derivedAlleleCount[i]>=RSDDataset->setSamples-sfsSlack);
	}

	for(i=0;i<size-RSDMuStat->windowSize+1;i++)
	{
		// SNP window range
		int snpf = i;
		int snpl = (int)(snpf + RSDMuStat->windowSize - 1);

		int winlsnpf = snpf;
		int winlsnpl = (int)(winlsnpf + RSDMuStat->windowSize/2 - 1);

		int winrsnpf = winlsnpl + 1;
		int winrsnpl = snpl;		
		
		// Window center (bp)
		windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
		windowStart = RSDChunk->sitePosition[snpf];
		windowEnd = RSDChunk->sitePosition[snpl];

		float isValid = 1.0;		

		for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
			if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
				isValid = 0.0;

		// Mu_Var
		muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
		muVar /= RSDDataset->setRegionLength;
		muVar /= RSDMuStat->windowSize;
		//muVar *= RSDDataset->setSNPs;
		muVar *= RSDDataset->preLoadedsetSNPs; // v2.4

		// Mu_SFS
		if(snpf>=1)
		{
			dCnt1 -= (RSDChunk->derivedAlleleCount[snpf-1]<=sfsSlack);
			dCnt1 += (RSDChunk->derivedAlleleCount[snpl]<=sfsSlack);

			dCntN -= (RSDChunk->derivedAlleleCount[snpf-1]>=RSDDataset->setSamples-sfsSlack);
			dCntN += (RSDChunk->derivedAlleleCount[snpl]>=RSDDataset->setSamples-sfsSlack);
		}
		float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
		muSfs = (float)dCnt1 + (float)dCntN*facN; 

		if(dCnt1+dCntN==0) 
			muSfs = 0.0000000001f;

		muSfs *= RSDDataset->muVarDenom; 
		muSfs /= (float)RSDMuStat->windowSize;

		// Mu_Ld
		int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;
		
		float tempTest = getPatternCounts (WIN_MODE_FXD, RSDMuStat, (int)RSDMuStat->windowSize, (int)RSDMuStat->windowSize, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
		muLd = tempTest;

		if(pcntexll + pcntexlr==0) 
		{
			muLd = 0.0000000001f;
		}
		else
		{
			muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
		}

		// Mu
		mu =  muVar * muSfs * muLd * isValid;

		// MuVar Max
		if (muVar > RSDMuStat->muVarMax)
		{
			RSDMuStat->muVarMax = muVar;
			RSDMuStat->muVarMaxLoc = windowCenter;
		}

		// MuSfs Max
		if (muSfs > RSDMuStat->muSfsMax)
		{
			RSDMuStat->muSfsMax = muSfs;
			RSDMuStat->muSfsMaxLoc = windowCenter;
		}

		// MuLd Max
		if (muLd > RSDMuStat->muLdMax)
		{
			RSDMuStat->muLdMax = muLd;
			RSDMuStat->muLdMaxLoc = windowCenter;
		}

		// Mu Max
		if (mu > RSDMuStat->muMax)
		{
			RSDMuStat->muMax = mu;
			RSDMuStat->muMaxLoc = windowCenter;
		}

		RSDMuStat->currentScoreIndex++;

		RSDMuStat_storeOutput (RSDMuStat, (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);

		//if(RSDCommandLine->fullReport==1)
		//	fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	//else
		//	fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
	}	
}
// End of Default SW implementation
#endif
#endif

void RSDMuStat_scanChunkWithMask (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, j, size = (int)RSDChunk->chunkSize;

	double windowCenter = 0.0, windowStart = 0.0, windowEnd = 0.0;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	int64_t sfsSlack = RSDCommandLine->sfsSlack;

	RSDMuStat_storeOutputConfigure(RSDCommandLine);

	// Mu_SFS initialization
	int dCnt1 = 0, dCntN = 0;
	for(i=0;i<RSDMuStat->windowSize - 0;i++)
	{
		int pID = RSDChunk->patternID[i];
		
		if(RSDPatternPool->poolDataWithMissing[pID]==1)
		{
			dCnt1 += (RSDPatternPool->poolDataAppliedMaskCount[pID]<=sfsSlack);
			dCntN += (RSDPatternPool->poolDataAppliedMaskCount[pID]>=RSDPatternPool->poolDataMaskCount[pID]-sfsSlack);
		}
		else // This can be omitted - test more first
		{
			dCnt1 += (RSDChunk->derivedAlleleCount[i]<=sfsSlack);
			dCntN += (RSDChunk->derivedAlleleCount[i]>=RSDDataset->setSamples-sfsSlack);			
		}
	}

	// TODO if needed 
	//float slack = 0.25;
	//int64_t exclRegSize = RSDMuStat->excludeRegionStop - RSDMuStat->excludeRegionStart;
	//int64_t exclRegSlack = (int64_t) (slack * exclRegSize);
	//RSDMuStat->excludeRegionStart -= exclRegSlack;
	//RSDMuStat->excludeRegionStop += exclRegSlack;

	for(i=0;i<size-RSDMuStat->windowSize+1;i++)
	{
		// SNP window range
		int snpf = i;
		int snpl = (int)(snpf + RSDMuStat->windowSize - 1);

		int winlsnpf = snpf;
		int winlsnpl = (int)(winlsnpf + RSDMuStat->windowSize/2 - 1);

		int winrsnpf = winlsnpl + 1;
		int winrsnpl = snpl;
		
		// Window center (bp)
		windowCenter = round((RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0);
		windowStart = RSDChunk->sitePosition[snpf];
		windowEnd = RSDChunk->sitePosition[snpl];

		float isValid = 1.0;

		for(j=0;j<RSDMuStat->excludeRegionsTotal;j++)
			if((windowEnd>=RSDMuStat->excludeRegionStart[j]) && (windowStart<=RSDMuStat->excludeRegionStop[j]))
				isValid = 0.0;

		// Mu_Var
		muVar = (float)(RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf]);
		muVar /= RSDDataset->setRegionLength;
		muVar /= RSDMuStat->windowSize;
		//muVar *= RSDDataset->setSNPs;
		muVar *= RSDDataset->preLoadedsetSNPs; // v2.4

		// Mu_SFS
		if(snpf>=1)
		{
			int pID = RSDChunk->patternID[snpf-1];
			
			if(RSDPatternPool->poolDataWithMissing[pID]==1)
			{
				dCnt1 -= (RSDPatternPool->poolDataAppliedMaskCount[pID]<=sfsSlack);
				dCntN -= (RSDPatternPool->poolDataAppliedMaskCount[pID]>=RSDPatternPool->poolDataMaskCount[pID]-sfsSlack);
			}
			else // This can be omitted - test more first
			{
				dCnt1 -= (RSDChunk->derivedAlleleCount[snpf-1]<=sfsSlack);
				dCntN -= (RSDChunk->derivedAlleleCount[snpf-1]>=RSDDataset->setSamples-sfsSlack);
			}

			pID = RSDChunk->patternID[snpl];

			if(RSDPatternPool->poolDataWithMissing[pID]==1)
			{
				dCnt1 += (RSDPatternPool->poolDataAppliedMaskCount[pID]<=sfsSlack);
				dCntN += (RSDPatternPool->poolDataAppliedMaskCount[pID]>=RSDPatternPool->poolDataMaskCount[pID]-sfsSlack);
			}
			else // This can be omitted - test more first
			{
				dCnt1 += (RSDChunk->derivedAlleleCount[snpl]<=sfsSlack);
				dCntN += (RSDChunk->derivedAlleleCount[snpl]>=RSDDataset->setSamples-sfsSlack);
			}
		}
		float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
		muSfs = (float)dCnt1 + (float)dCntN*facN; 

		if(dCnt1+dCntN==0) 
			muSfs = 0.0000000001f;

		muSfs *= RSDDataset->muVarDenom; 
		muSfs /= (float)RSDMuStat->windowSize;

		// Mu_Ld
		int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;
		
		float tempTest = getPatternCounts (WIN_MODE_FXD, RSDMuStat, (int)RSDMuStat->windowSize, (int)RSDMuStat->windowSize, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
		muLd = tempTest;

		if(pcntexll + pcntexlr==0) 
		{
			muLd = 0.0000000001f;
		}
		else
		{
			muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
		}

		// Mu
		mu =  muVar * muSfs * muLd * isValid;

		// MuVar Max
		if (muVar > RSDMuStat->muVarMax)
		{
			RSDMuStat->muVarMax = muVar;
			RSDMuStat->muVarMaxLoc = windowCenter;
		}

		// MuSfs Max
		if (muSfs > RSDMuStat->muSfsMax)
		{
			RSDMuStat->muSfsMax = muSfs;
			RSDMuStat->muSfsMaxLoc = windowCenter;
		}

		// MuLd Max
		if (muLd > RSDMuStat->muLdMax)
		{
			RSDMuStat->muLdMax = muLd;
			RSDMuStat->muLdMaxLoc = windowCenter;
		}

		// Mu Max
		if (mu > RSDMuStat->muMax)
		{
			RSDMuStat->muMax = mu;
			RSDMuStat->muMaxLoc = windowCenter;
		}

		RSDMuStat->currentScoreIndex++;

		RSDMuStat_storeOutput (RSDMuStat, (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);

		//if(RSDCommandLine->fullReport==1)
		//	fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	//else
		//	fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	

	}
}
#endif
