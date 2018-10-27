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
