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

float 	getPatternCounts 		(int * pCntVec, int offset, int * patternID, int p0, int p1, int p2, int p3, int * pcntl, int * pcntr, int * pcntexll, int * pcntexlr);
float 	getCrossLD 			(RSDPatternPool_t * pp, int p0, int p1, int p2, int p3, int samples);
float 	getRegionLD 			(RSDPatternPool_t * pp, int p0, int p1, int samples);
float 	pwLD 				(RSDPatternPool_t * pp, int p1, int p2, int samples);
void	(*RSDMuStat_scanChunk) 		(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
void 	RSDMuStat_scanChunkBinary	(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
void 	RSDMuStat_scanChunkWithMask	(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);

RSDMuStat_t * RSDMuStat_new (void)
{
	RSDMuStat_t * mu = NULL;

	mu = (RSDMuStat_t *)malloc(sizeof(RSDMuStat_t));
	assert(mu!=NULL);

	mu->reportFP = NULL;
	
	mu->windowSize = WINDOW_SIZE;
	
	mu->pCntVec = NULL; 

	mu->muVarMax = 0.0f; 
	mu->muVarMaxLoc = 0.0f;

	mu->muSfsMax = 0.0f; 
	mu->muSfsMaxLoc = 0.0f;

	mu->muLdMax = 0.0f;
	mu->muLdMaxLoc = 0.0f;

	mu->muMax = 0.0f; 
	mu->muMaxLoc = 0.0f;

#ifdef _MUMEM
	mu->muReportBufferSize = 1;
	mu->muReportBuffer = (float*)malloc(sizeof(float)*7);
	assert(mu->muReportBuffer!=NULL);
#endif
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
	
	MemoryFootprint += sizeof(int)*((unsigned long)(mu->windowSize*4));
	if(mu->pCntVec!=NULL)
		free(mu->pCntVec);

#ifdef _MUMEM	
	if(mu->muReportBuffer!=NULL)
	{
		MemoryFootprint += sizeof(float)*((unsigned long)(mu->muReportBufferSize*7));

		free(mu->muReportBuffer);
		mu->muReportBuffer =  NULL;
	}
#endif

	free(mu);
}

void RSDMuStat_init (RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);

	RSDMuStat->muVarMax = 0.0f; 
	RSDMuStat->muVarMaxLoc = 0.0f;

	RSDMuStat->muSfsMax = 0.0f; 
	RSDMuStat->muSfsMaxLoc = 0.0f;

	RSDMuStat->muLdMax = 0.0f;
	RSDMuStat->muLdMaxLoc = 0.0f;

	RSDMuStat->muMax = 0.0f; 
	RSDMuStat->muMaxLoc = 0.0f;

	if(RSDMuStat->pCntVec==NULL)
	{
		RSDMuStat->pCntVec = (int *)malloc(sizeof(int)*((unsigned long)(RSDMuStat->windowSize*4)));
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

void RSDMuStat_setReportNamePerSet (RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine, FILE * fpOut, RSDDataset_t * RSDDataset)
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

	return 0;
}

#ifdef _HW
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
	int * tmp_vec_a = (int *)malloc(sizeof(int)*(window_size/2)); // patterns left
	int * tmp_vec_b = (int *)malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_c = (int *)malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_d = (int *)malloc(sizeof(int)*(window_size/2));

	int * tmp_vec_a_cnt = (int *)malloc(sizeof(int)*(window_size/2)); // patterns left
	int * tmp_vec_b_cnt = (int *)malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_c_cnt = (int *)malloc(sizeof(int)*(window_size/2));
	int * tmp_vec_d_cnt = (int *)malloc(sizeof(int)*(window_size/2));
	
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
	float * vec_a_in = (float *)malloc(sizeof(float)*vec_size);
	int * vec_b_in = (int *)malloc(sizeof(int)*vec_size);
	int * vec_c_in = (int *)malloc(sizeof(int)*vec_size);

	memcpy(vec_a_in, RSDChunk->sitePosition, sizeof(float)*vec_size);
	memcpy(vec_b_in, RSDChunk->derivedAlleleCount, sizeof(int)*vec_size);
	memcpy(vec_c_in, RSDChunk->patternID, sizeof(int)*vec_size);

	// Intermediate vectors
	float * vec_a_tmp = (float *)malloc(sizeof(float)*vec_size);
	float * vec_b_tmp = (float *)malloc(sizeof(float)*vec_size);
	float * vec_c_tmp = (float *)malloc(sizeof(float)*vec_size);

	// Output vectors
	float * window_loc = (float*)malloc(sizeof(float)*vec_size);
	float * vec_out = (float *)malloc(sizeof(float)*vec_size);

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
#else
#ifdef _MLT
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
	    patternMemLeft[WINDOW_SIZE],
	    patternMemRightSize = 0,
	    patternMemRight[WINDOW_SIZE],
            patternMemExclusiveSize = 0,
	    patternMemLeftExclusive[WINDOW_SIZE],
	    patternMemRightExclusive[WINDOW_SIZE],
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
	      rotbuf[WINDOW_SIZE*3]; // circular queue

#ifdef _MUMEM
	int muReportBufferIndex = -1;

	if(size>RSDMuStat->muReportBufferSize)
	{
		RSDMuStat->muReportBufferSize = size;
		RSDMuStat->muReportBuffer = realloc(RSDMuStat->muReportBuffer, sizeof(float)*(unsigned long)RSDMuStat->muReportBufferSize*7);
		assert(RSDMuStat->muReportBuffer);	
	}
#endif

	// Circular queue init
	memcpy(rotbuf, RSDChunk->chunkData, WINDOW_SIZE*sizeof(float)*3);

	// SFS preprocessing
	for(i=0;i<RSDMuStat->windowSize;i++)
	{
		dCnt1 += (((int)RSDChunk->chunkData[i*3+1])==1);
		dCntN += (((int)RSDChunk->chunkData[i*3+1])==RSDDataset->setSamples-1);
	}
		
	// LD preprocessing: init
	for(i=0;i<WINDOW_SIZE/2;i++)
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
	for(i=0;i<WINDOW_SIZE/2;i++)
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
		int rotbufIndex_temp = rotbufIndex + (WINDOW_SIZE/2);

		if(rotbufIndex>=(WINDOW_SIZE/2)-1)
			rotbufIndex_temp -= (WINDOW_SIZE-1);

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
		int rotbufIndex_temp = rotbufIndex + (WINDOW_SIZE/2);

		if(rotbufIndex>=(WINDOW_SIZE/2)-1)
			rotbufIndex_temp -= (WINDOW_SIZE-1);

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

		int F_Le_vec[WINDOW_SIZE/2], F_Le_ind=-1; // stacks
		int F_Ri_vec[WINDOW_SIZE/2], F_Ri_ind=-1;
		int F_xLe_vec[WINDOW_SIZE/2], F_xLe_ind=-1;
		int F_xRi_vec[WINDOW_SIZE/2], F_xRi_ind=-1;

		for(j=0;j<WINDOW_SIZE/2;j++)
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
#else
void RSDMuStat_scanChunkBinary (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, size = (int)RSDChunk->chunkSize;

	float windowCenter = 0.0f, windowStart = 0.0f, windowEnd = 0.0f;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	// Mu_SFS initialization
	int dCnt1 = 0, dCntN = 0;
	for(i=0;i<RSDMuStat->windowSize - 0;i++)
	{
		dCnt1 += (RSDChunk->derivedAlleleCount[i]==1);
		dCntN += (RSDChunk->derivedAlleleCount[i]==RSDDataset->setSamples-1);
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
		windowCenter = (RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0f;
		windowStart = RSDChunk->sitePosition[snpf];
		windowEnd = RSDChunk->sitePosition[snpl];

		// Mu_Var
		muVar = RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf];
		muVar /= RSDDataset->setRegionLength;
		muVar /= RSDMuStat->windowSize; 

		// Mu_SFS
		if(snpf>=1)
		{
			dCnt1 -= (RSDChunk->derivedAlleleCount[snpf-1]==1);
			dCnt1 += (RSDChunk->derivedAlleleCount[snpl]==1);

			dCntN -= (RSDChunk->derivedAlleleCount[snpf-1]==RSDDataset->setSamples-1);
			dCntN += (RSDChunk->derivedAlleleCount[snpl]==RSDDataset->setSamples-1);
		}
		float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
		muSfs = (float)dCnt1 + (float)dCntN*facN; 

		if(dCnt1+dCntN==0) 
			muSfs = 0.000001f;

		muSfs /= (float)RSDMuStat->windowSize;

		// Mu_Ld
		int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;
		
		float tempTest = getPatternCounts (RSDMuStat->pCntVec, (int)RSDMuStat->windowSize, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
		muLd = tempTest;

		if(pcntexll + pcntexlr==0) 
		{
			muLd = 0.000001f;
		}
		else
		{
			muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
		}


		// Mu
		mu =  muVar * muSfs * muLd;

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

		if(RSDCommandLine->fullReport==1)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	
	}
}
#endif

void RSDMuStat_scanChunkWithMask (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommandLine!=NULL);
	assert(RSDPatternPool!=NULL);

	int i, size = (int)RSDChunk->chunkSize;

	float windowCenter = 0.0f, windowStart = 0.0f, windowEnd = 0.0f;

	float muVar = 0.0f;
	float muSfs = 0.0f;
	float muLd = 0.0f;
	float mu = 0.0f;

	// Mu_SFS initialization
	int dCnt1 = 0, dCntN = 0;
	for(i=0;i<RSDMuStat->windowSize - 0;i++)
	{
		int pID = RSDChunk->patternID[i];
		
		if(RSDPatternPool->poolDataWithMissing[pID]==1)
		{
			
			dCnt1 += (RSDPatternPool->poolDataAppliedMaskCount[pID]==1);
			dCntN += (RSDPatternPool->poolDataAppliedMaskCount[pID]==RSDPatternPool->poolDataMaskCount[pID]-1);
		}
		else // This can be ommitted - test more first
		{
			dCnt1 += (RSDChunk->derivedAlleleCount[i]==1);
			dCntN += (RSDChunk->derivedAlleleCount[i]==RSDDataset->setSamples-1);			
		}
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
		windowCenter = (RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0f;
		windowStart = RSDChunk->sitePosition[snpf];
		windowEnd = RSDChunk->sitePosition[snpl];

		// Mu_Var
		muVar = RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf];
		muVar /= RSDDataset->setRegionLength;
		muVar /= RSDMuStat->windowSize; 

		// Mu_SFS
		if(snpf>=1)
		{
			int pID = RSDChunk->patternID[snpf-1];
			
			if(RSDPatternPool->poolDataWithMissing[pID]==1)
			{
				dCnt1 -= (RSDPatternPool->poolDataAppliedMaskCount[pID]==1);
				dCntN -= (RSDPatternPool->poolDataAppliedMaskCount[pID]==RSDPatternPool->poolDataMaskCount[pID]-1);
			}
			else // This can be ommitted - test more first
			{
				dCnt1 -= (RSDChunk->derivedAlleleCount[snpf-1]==1);
				dCntN -= (RSDChunk->derivedAlleleCount[snpf-1]==RSDDataset->setSamples-1);
			}

			pID = RSDChunk->patternID[snpl];

			if(RSDPatternPool->poolDataWithMissing[pID]==1)
			{
				dCnt1 += (RSDPatternPool->poolDataAppliedMaskCount[pID]==1);
				dCntN += (RSDPatternPool->poolDataAppliedMaskCount[pID]==RSDPatternPool->poolDataMaskCount[pID]-1);
			}
			else // This can be ommitted - test more first
			{
				dCnt1 += (RSDChunk->derivedAlleleCount[snpl]==1);
				dCntN += (RSDChunk->derivedAlleleCount[snpl]==RSDDataset->setSamples-1);
			}
		}
		float facN = 1.0;//(float)RSDChunk->derAll1CntTotal/(float)RSDChunk->derAllNCntTotal;
		muSfs = (float)dCnt1 + (float)dCntN*facN; 

		if(dCnt1+dCntN==0) // Adapt this
			muSfs = 0.000001f;

		muSfs /= (float)RSDMuStat->windowSize;

		// Mu_Ld
		int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0;
		
		float tempTest = getPatternCounts (RSDMuStat->pCntVec, (int)RSDMuStat->windowSize, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr);
		muLd = tempTest;

		if(pcntexll + pcntexlr==0) // Adapt this
		{
			muLd = 0.000001f;
		}
		else
		{
			muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
		}

		// Mu
		mu =  muVar * muSfs * muLd;

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

		if(RSDCommandLine->fullReport==1)
			fprintf(RSDMuStat->reportFP, "%.0f\t%.0f\t%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", (double)windowCenter, (double)windowStart, (double)windowEnd, (double)muVar, (double)muSfs, (double)muLd, (double)mu);	else
			fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\n", (double)windowCenter, (double)mu);	

	}
}
#endif
