#include "RAiSD.h"

RSDMuStat_t * RSDMuStat_new (void)
{
	RSDMuStat_t * mu = NULL;

	mu = (RSDMuStat_t *)malloc(sizeof(RSDMuStat_t));
	assert(mu!=NULL);
	
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
	
	MemoryFootprint += sizeof(int)*(mu->windowSize*4);
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
		RSDMuStat->pCntVec = (int *)malloc(sizeof(int)*(RSDMuStat->windowSize*4));
		assert(RSDMuStat->pCntVec!=NULL);
	}
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
		fclose(RSDMuStat->reportFP);

	char tstring[STRING_SIZE];
	strcpy(tstring, RSDMuStat->reportName);
	strcat(tstring, ".");
	strcat(tstring, RSDDataset->setID);

	RSDMuStat->reportFP = fopen(tstring, "r");
	if(RSDMuStat->reportFP!=NULL && RSDCommandLine->overwriteOutput==0)
	{
		fprintf(fpOut, "\nERROR: Output report file %s exists. Use -f to overwrite it.\n\n", tstring);
		fprintf(stderr, "\nERROR: Output report file %s exists. Use -f to overwrite it.\n\n", tstring);
		exit(0);	
	}

	if(RSDMuStat->reportFP!=NULL)
		fclose(RSDMuStat->reportFP);

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

float getRegionLD (RSDPatternPool_t * pp, int * patternID, int p0, int p1, int samples)
{
	double totalLD = 0.0f;
	int i, j;
	for(i=p0;i<=p1-1;i++)
	{
		for(j=p0+1;j<=p1;j++)
			totalLD += (double) pwLD(pp, i, j, samples);
	}
	return totalLD;
}

float getCrossLD (RSDPatternPool_t * pp, int * patternID, int p0, int p1, int p2, int p3, int samples)
{
	double totalLD = 0.0f;
	int i, j;
	for(i=p0;i<=p1;i++)
	{
		for(j=p2;j<=p3;j++)
			totalLD += (double) pwLD(pp, i, j, samples);
	}
	return totalLD;
}

float getPatternCounts (int * pCntVec, int offset, int * patternID, int p0, int p1, int p2, int p3, int * pcntl, int * pcntr, int * pcntexll, int * pcntexlr, float * pcntexllf, float * pcntexlrf, int * exsnpcntleft, int * exsnpcntright, RSDPatternPool_t * RSDPatternPool, int samples)
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
		
	*pcntl = list_left_size; // number of patterns in left sudwindow

	//float ldl = getRegionLD (RSDPatternPool, patternID, p0, p1, samples);


	

	int checksum = 0;
	for(i=0;i<list_left_size;i++)
		checksum += list_left_cnt[i];

	assert(checksum==p1-p0+1);


	int list_right_size = 0;
	int * list_right = &(pCntVec[2*offset]); //(int *) malloc(sizeof(int));
	int * list_right_cnt = &(pCntVec[3*offset]); //(int *) malloc(sizeof(int));
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
			//list_right = realloc(list_right, sizeof(int)*list_right_size);
			//list_right_cnt = realloc(list_right_cnt, sizeof(int)*list_right_size);
			list_right[list_right_size-1] = patternID[i];
			list_right_cnt[list_right_size-1] = 1;
		}
	}
		
	*pcntr = list_right_size;

	//printf("\nstart %d end %d: ", p0, p3);
	//for(i=p0;i<=p3;i++)
	//{
	//	printf("%d ", patternID[i]);
	//}
	//fflush(stdout);
	//assert(0);

	int excl_pat_left_total = 0;
	int * excl_pat_left = (int *) malloc(sizeof(int)*WINDOW_SIZE);
	int * excl_pat_left_cnt = (int *) malloc(sizeof(int)*WINDOW_SIZE);
	int excl_pat_rght_total = 0;
	int * excl_pat_rght = (int *) malloc(sizeof(int)*WINDOW_SIZE);
	int * excl_pat_rght_cnt = (int *) malloc(sizeof(int)*WINDOW_SIZE);
	int shared_pat_total = 0;
	int * shared_pat = (int *) malloc(sizeof(int)*WINDOW_SIZE);
	int * shared_pat_cnt = (int *) malloc(sizeof(int)*WINDOW_SIZE);

	for(i=0;i<WINDOW_SIZE;i++)
	{
		excl_pat_left[i] = 0;
		excl_pat_left_cnt [i] = 0;
		excl_pat_rght[i] = 0;
		excl_pat_rght_cnt [i] = 0;
		shared_pat[i] = 0;
		shared_pat_cnt [i] = 0;
	}

	int max_cnt_of_excl_snps_for_one_excl_pat = 0;
	for(i=0;i<list_left_size;i++)
	{
		int match = 0;
		for(j=0;j<list_right_size;j++)
		{
			if(list_left[i] == list_right[j])
			{
				match = 1;
				j = list_right_size;
			}
		}
		if(match==0)
		{
			excl_pat_left [excl_pat_left_total] = list_left[i];
			excl_pat_left_cnt [excl_pat_left_total++] = list_left_cnt[i];

			if(list_left_cnt[i]>=max_cnt_of_excl_snps_for_one_excl_pat)
				max_cnt_of_excl_snps_for_one_excl_pat = list_left_cnt[i];
		}
	}


	int max_cnt_of_excl_snps_for_one_excl_pat_rght = 0;
	for(i=0;i<list_right_size;i++)
	{
		int match = 0;
		for(j=0;j<list_left_size;j++)
		{
			if(list_left[j] == list_right[i])
			{
				match = 1;
				j = list_left_size;
			}
		}
		if(match==0)
		{
			excl_pat_rght [excl_pat_rght_total] = list_right[i];
			excl_pat_rght_cnt [excl_pat_rght_total++] = list_right_cnt[i];

			if(list_right_cnt[i]>=max_cnt_of_excl_snps_for_one_excl_pat_rght)
				max_cnt_of_excl_snps_for_one_excl_pat_rght = list_right_cnt[i];
		}
	}


	int max_cnt_of_shared_snps_for_one_shared_pat = 0;
	for(i=0;i<list_left_size;i++)
	{
		int match = 0;
		for(j=0;j<list_right_size;j++)
		{
			if(list_left[i] == list_right[j])
			{
				shared_pat [shared_pat_total] = list_left[i];
				shared_pat_cnt [shared_pat_total++] = list_left_cnt[i] + list_right_cnt[i] ;

				if(list_left_cnt[i] + list_right_cnt[j] >=max_cnt_of_shared_snps_for_one_shared_pat)
					max_cnt_of_shared_snps_for_one_shared_pat = list_left_cnt[i] + list_right_cnt[j];
				
				j = list_right_size;
			}
		}
	}

	float tempTest = 0.0f;
	
	tempTest = ((float)max_cnt_of_excl_snps_for_one_excl_pat * (float)max_cnt_of_excl_snps_for_one_excl_pat_rght);// / ((float)max_cnt_of_shared_snps_for_one_shared_pat+1.0);
	//if(tempTest<(float)max_cnt_of_excl_snps_for_one_excl_pat_rght)
	//	tempTest = (float)max_cnt_of_excl_snps_for_one_excl_pat_rght;
//	tempTest /= ((*pcntl)*(*pcntr));
//	tempTest = (float)max_cnt_of_shared_snps_for_one_shared_pat+0.00001;

	free(excl_pat_left);
	free(excl_pat_left_cnt);
	free(excl_pat_rght);
	free(excl_pat_rght_cnt);
	free(shared_pat);
	free(shared_pat_cnt);
	

/*
	if(list_left_size==PCNT)
	{
	//printf("\n N\t");
	//for(i=p0;i<=p1;i++)
	//	printf("%d ", patternID[i]);

	//printf(" P-l %d LD %f\n", list_left_size, ldl);	

	if(ldl<minld)
		minld = ldl;

	if(ldl>maxld)
		maxld = ldl;


	cntld++;
	accumld += ldl;

	}*/
	

	checksum = 0;
	for(i=0;i<list_right_size;i++)
		checksum += list_right_cnt[i];

	assert(checksum==p3-p2+1);


	//int * list_left_backup = (int *) malloc(sizeof(int)*list_left_size);
	//for(i=0;i<list_left_size;i++)
	//	list_left_backup [i] = list_left[i];

	//int * list_right_backup = (int *) malloc(sizeof(int)*list_right_size);
	//for(i=0;i<list_right_size;i++)
	//	list_right_backup [i] = list_right[i];
	

	int excl_left = list_left_size;
	for(i=0;i<list_left_size;i++)
	{
		int match = 0;
		for(j=0;j<list_right_size;j++)
		{
			if(list_left[i] == list_right[j])
			{
				excl_left--;
				//list_left_backup [i] = -1;
			}
		}
	}
	*pcntexll = excl_left;
	assert(excl_pat_left_total==excl_left);

	int excl_right = list_right_size;
	for(i=0;i<list_right_size;i++)
	{
		int match = 0;
		for(j=0;j<list_left_size;j++)
		{
			if(list_right[i] == list_left[j])
			{
				excl_right--;
				//list_right_backup [i] = -1;
			}
		}
	}
	*pcntexlr = excl_right;
	assert(excl_pat_rght_total==excl_right);
	assert(list_left_size - excl_left == list_right_size - excl_right);


	/* Exclusive SNPs left */
	/* */
	int excntsnpsl = 0;
	for(i=p0+1;i<=p1;i++)
	{
		int match = 0;
		for(j=p2+1;j<=p3;j++)
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
	//*pcntexll = (*pcntexll) * excntsnpsl;
	//*exsnpcntleft = excntsnpsl;
	int excntsnpsr = 0;
	for(i=p2+1;i<=p3;i++)
	{
		int match = 0;
		for(j=p0+1;j<=p1;j++)
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
	}/* */
	//*exsnpcntright = excntsnpsr;
	//*pcntexlr = (*pcntexlr) * excntsnpsr;

	// New Test - (Max Number of Exclusive SNPs left + Max NUmber of Exclusive SNPs right ) / shared SNPs or Patterns + 1
	

	
		
	//*pcntexll =  (*pcntexll) + (*pcntexlr);
	//*pcntexlr =  excntsnpsl + excntsnpsr;


	return tempTest;

	//*pcntexll = 1.0;
	//fprintf(pattest, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", *pcntl, *pcntr, *pcntexll, *pcntexlr, *pcntl-*pcntexll, *pcntr - *pcntexlr, excntsnpsl,excntsnpsr);

	//*pcntexll = ((*pcntexll)+(*pcntexlr));//*(excntsnpsl+excntsnpsr);
	//*pcntexlr = 0;
	//int pairs = ((WINDOW_SIZE/2) * ((WINDOW_SIZE/2)-1))/2;
 
	//float crossld = getCrossLD (RSDPatternPool, patternID, p0, p1, p2, p3, samples);
	//float ldl = getRegionLD (RSDPatternPool, patternID, p0, p1, samples);
	//float ldr = getRegionLD (RSDPatternPool, patternID, p2, p3, samples);
	//printf("%d %d %d %d %d %d -> %f %f %f: %f %f\n", *pcntl, *pcntr, *pcntexll, *pcntexlr, excntsnpsl,excntsnpsr,  ldl/pairs, ldr/pairs, crossld/(WINDOW_SIZE*WINDOW_SIZE/4), (float)(*pcntexll)/(float)(*pcntl), (float)(*pcntexlr)/(float)(*pcntr), (crossld/(WINDOW_SIZE*WINDOW_SIZE/4))/(((float)(*pcntexlr)/(float)(*pcntr))*((float)(*pcntexll)/(float)(*pcntl))));
	

	//float ldl = getRegionLD (RSDPatternPool, patternID, p0, p1, samples);
	
	//float ldr = getRegionLD (RSDPatternPool, patternID, p2, p3, samples);

/* *
	//if((*pcntexll)+(*pcntexlr)==PCNT)
	{
	//printf("\n N\t");
	//for(i=p0;i<=p1;i++)
	//	printf("%d ", patternID[i]);

	printf(" Pats %d LD %f\n", (*pcntexll)+(*pcntexlr), crossld);	

	if(crossld<minld)
		minld = crossld;

	if(crossld>maxld)
		maxld = crossld;


	cntld++;
	accumld += crossld;

	}/* */



	//*pcntexlr = excntsnpsr;
	//if((*pcntexll)==0)
	//	(*pcntexll)++;

	//if((*pcntexlr)==0)
	//	(*pcntexlr)++;

	//*pcntexllf = ((float)excntsnpsl/WINDOW_SIZE) *  ((float)(*pcntexll)/(*pcntl));// * excntsnpsl;
	//*pcntexlrf = ((float)excntsnpsr/WINDOW_SIZE) * ((float)(*pcntexlr)/(*pcntr));//* excntsnpsr;
	//printf("%d %d %d  %f - %d %f\n", excntsnpsl, (*pcntexll), (*pcntl), *pcntexllf, excntsnpsr, *pcntexlrf);
	
	//int commonSNPs = excl_left + excl_right;//list_left_size - excl_left; // number of common patterns
	//int dif = list_left_size>list_right_size?list_left_size-list_right_size:list_right_size-list_left_size;
	//float indicator = (float)(excl_left + excl_right)/(float)(dif);
	//printf("Ind %f LD %f\n", indicator, crossld);
	//commonSNPs = (int)indicator;
/*	printf("\n\nL: ");
	for(i=p0;i<=p1;i++)
		printf("%d ", patternID[i]);

	printf("\n\nPatL: %d / %d\n", excl_left, list_left_size);
	for(i=0;i<list_left_size;i++)
		printf("%d ", list_left[i]);
	

	printf("\n\nR: ");
	for(i=p2;i<=p3;i++)
		printf("%d ", patternID[i]);

	printf("\n\nPatR: %d / %d\n", excl_right, list_right_size);
	for(i=0;i<list_right_size;i++)
		printf("%d ", list_right[i]);
*/

/* */
/*	if(commonSNPs==PCNT)
	{
	//printf("\n N\t");
	//for(i=p0;i<=p1;i++)
	//	printf("%d ", patternID[i]);

	//printf(" P-l %d LD %f\n", list_left_size, ldl);	

	if(crossld<minld)
		minld = crossld;

	if(crossld>maxld)
		maxld = crossld;


	cntld++;
	accumld += crossld;

	} *//**/


	//fflush(stdout);
	//assert(0);

/*	// New Test: exclusive pattern * corresponding snp count
	int exclusive_snps_left = p1-p0+1;

	for(i=p0;i<=p1;i++)
	{
		for(j=0;j<list_right_size;j++)	
		{
			if(patternID[i]==list_right[i])
			{
				exclusive_snps_left--;
				j=list_right_size;
			}
		}
	}


	int exclusive_snps_right = p3-p2+1;

	for(i=p2;i<=p3;i++)
	{
		for(j=0;j<list_left_size;j++)	
		{
			if(patternID[i]==list_left[i])
			{
				exclusive_snps_right--;
				j=list_left_size;
			}
		}
	}
	

	//ldmap[exclusive_snps_left][exclusive_snps_right] += crossld;
	*pcntexll = exclusive_snps_left;
	*pcntexlr = exclusive_snps_right;
*/
	//*pcntexll = (int)((float)exclusive_snps_left/(float)excl_left);
	//*pcntexlr = (int)((float)exclusive_snps_right/(float)excl_right);


	//*pcntexll = (int)((float)excl_left*(float)exclusive_snps_left); // kinda okay
	//*pcntexlr = (int)((float)excl_right*(float)exclusive_snps_right); // kinda okay

	//*pcntexll = excl_left + excl_right; //(int)((float)((float)excl_left*(float)exclusive_snps_left) + (float)((float)excl_right*(float)exclusive_snps_right));
	//*pcntexlr = (exclusive_snps_left + exclusive_snps_right);


	//for(i=0;i<list_left_size;i++)
	//	if(list_left_backup [i]!=-1)
	//		(*exsnpcntleft) += list_left_cnt[i];


	//for(i=0;i<list_right_size;i++)
	//	if(list_right_backup [i]!=-1)
	//		(*exsnpcntright) += list_right_cnt[i];


	//printf("%d %d\n", *exsnpcntleft, *exsnpcntright);
	//assert((*exsnpcntleft)!=0);
	//assert((*exsnpcntright)!=0);

	//free(list_left);
	//free(list_left_backup);
	//free(list_left_cnt);

	//free(list_right);
	//free(list_right_backup);
	//free(list_right_cnt);


}

void RSDMuStat_scanChunk (RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine)
{
	int i, size = RSDChunk->chunkSize;

	float windowCenter = 0.0f;

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

	minld = 1000000000000000.0f;
	maxld = 0.0f;
	accumld = 0.0f;
	cntld = 0;

	int j;
	//ldmap = (double **)malloc(sizeof(double*)*WINDOW_SIZE);
	//for(i=0;i<WINDOW_SIZE;i++)
	//{
	//	ldmap[i] = (double*)malloc(sizeof(double)*WINDOW_SIZE);
	//	for(j=0;j<WINDOW_SIZE;j++)
	//		ldmap[i][j] = 0.0f;
	//}

	for(i=0;i<size-RSDMuStat->windowSize+1;i++)
	{
		//printf("window at %d\n", i);
		//fflush(stdout);

		// SNP window range
		int snpf = i;
		int snpl = snpf + RSDMuStat->windowSize - 1;

		int winlsnpf = snpf;
		int winlsnpl = winlsnpf + RSDMuStat->windowSize/2 - 1;

		int winrsnpf = winlsnpl + 1;
		int winrsnpl = snpl;
		
		// Window center (bp)
		windowCenter = (RSDChunk->sitePosition[snpf] + RSDChunk->sitePosition[snpl]) / 2.0f;

		// Mu_Var
		//printf("%d: %f %f\n", i, RSDChunk->sitePosition[snpf], RSDChunk->sitePosition[snpl]);
		//fflush(stdout);
		muVar = RSDChunk->sitePosition[snpl] - RSDChunk->sitePosition[snpf];
		//muVar /= RSDCommandLine->regionLength; // TODO: maybe remove this
		//muVar *= RSDDataset->setSNPs;                                       //     <<< ---- Check why we need this because it generates issues with the multistep parsing approach
		//muVar /= RSDCommandLine->regionLength;
		muVar /= RSDDataset->setRegionLength;
		muVar /= RSDMuStat->windowSize; // This doesnot affect the results at all.
		// TODO: maybe divide by window size in snps / to make it comparable with different window sizes
		// This is the number of expected SNPs if the win snps were uniformly distributed.

		//printf("VAR parameters: %d %d %d", RSDDataset->setSNPs, RSDCommandLine->regionLength, RSDMuStat->windowSize);

		//fflush(stdout);
		//assert(0);

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
		// TODO: weight vasi total arithmo of singletons X kai arithmo K apo N-1 sites kai weight sta N-1 to X/K
		// This has been tested, gives worse performance

		if(dCnt1+dCntN==0) // TODO: Handle this properly
			muSfs = 0.000001f;

		muSfs /= (float)RSDMuStat->windowSize;
		//muSfs *= (float)RSDMuStat->windowSize;
		//muSfs /= RSDCommandLine->regionLength;


		/*if((int)windowCenter==20070)
		{
			printf("SFS parameters: %f %f %d %d", muSfs, (float)RSDMuStat->windowSize, dCnt1, dCntN);

			fflush(stdout);
			assert(0);
		}*/



		// Mu_Ld
		int pcntl = 0, pcntr = 0, pcntexll=0, pcntexlr=0, exsnpcntleft = 0, exsnpcntright = 0;
		float pcntexllf=0.0, pcntexlrf=0.0;
		float tempTest = getPatternCounts (RSDMuStat->pCntVec, RSDMuStat->windowSize, RSDChunk->patternID, winlsnpf, winlsnpl, winrsnpf, winrsnpl, &pcntl, &pcntr, &pcntexll, &pcntexlr, &pcntexllf, &pcntexlrf , &exsnpcntleft, &exsnpcntright, RSDPatternPool, RSDDataset->setSamples);
		//printf("v: %f\n", tempTest);
		muLd = tempTest;

/* */
		if(pcntexll + pcntexlr==0) // TODO: Handle this properly
		{
			muLd = 0.000001f;
		}
		else
		{
			muLd = (((float)pcntexll)+((float)pcntexlr)) / ((float)(pcntl * pcntr));
			//muLd = (((float)pcntexll)*((float)pcntexlr)) / ((float)(pcntl * pcntr));
		}



		//if((int)windowCenter==44850)
		//{
		//	printf("LD parameters: %f %d %d %d %d", muLd, pcntexll, pcntexlr, pcntl, pcntr);
//
		//	fflush(stdout);
			//assert(0);
		//}




//		muLd *= tempTest;
/* */		//muLd = (((float)pcntexll)*((float)pcntexlr)) / ((float)(pcntl * pcntr));

		//muLd = ((((float)(((float)pcntexll) * ((float)pcntexlr))) / ((float)(pcntl * pcntr))));
		//muLd = ((float)(pcntl * pcntr))/((float)(((float)pcntexll) * ((float)pcntexlr)));
		//muLd = 1.0 / ((float)(pcntl * pcntr));

		//muLd = 

		///muLd *= (float)RSDMuStat->windowSize;
		//muLd /= RSDCommandLine->regionLength;

		//muLd = ((float)(((float)pcntexllf) * ((float)pcntexlrf))) / ((float)(pcntl * pcntr));
		//muLd = ((float)(((float)pcntexll)*((float)pcntexlr)))*(((float)(pcntl + pcntr))) / ((float)(pcntl * pcntr));
		//muLd = ((float)(((float)pcntexll/(float)pcntl) + ((float)pcntexlr/(float)pcntr))) / ((float)(pcntl * pcntr));
		//muLd = ((float)(((float)pcntexll) * ((float)pcntexlr))) / ((float)(pcntl * pcntr));


		//muLd = ((float)(pcntl * pcntr)) / ((float)(((float)pcntexll) + ((float)pcntexlr))) ; 
		 
/*		int dif = pcntl>pcntr?pcntl-pcntr:pcntr-pcntl;
		assert(dif>=0);
	
		//dif++;
		float diff = (float)dif + 0.0000001;
	//printf("dif %d\n", dif);

		muLd = (((float)((float)pcntexll + (float)pcntexlr))/diff) / ((float)(pcntl * pcntr)); 
*/
		//muLd = (float)(pcntl + pcntr) / ((float)pcntexll * (float)pcntexlr);
		//muLd = (4.0f / (((float)pcntexll + (float)pcntexlr) * (WINDOW_SIZE/2) *(WINDOW_SIZE/2)    )) / ((float)(pcntl + pcntr)/((float)(pcntl * pcntr)*(WINDOW_SIZE/2)*((WINDOW_SIZE/2)-1)));

		//muLd *= RSDDataset->setSNPs;
		//muLd /= RSDCommandLine->regionLength;

		//muLd = 1.0 / ((float)(pcntl * pcntr)) * ((float)((float)pcntexll * (float)pcntexlr)); // REVERSE: [* -> +] kai [+ -> *]


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

		//FILE * fp = fopen("results.txt", "a");
		//fprintf(fp, "%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", windowCenter, muVar, muSfs, muLd, mu);
		//fclose(fp);
		//printf("%f -> %f\n", windowCenter, mu);

	
		fprintf(RSDMuStat->reportFP, "%.0f\t%.3e\t%.3e\t%.3e\t%.3e\n", windowCenter, muVar, muSfs, muLd, mu);
	}

	//printf("\nMin Cross LD  %f\nMax Cross LD %f\nAvg Cross %f\n\n", minld, maxld, accumld/(double)cntld);

	//printf("\n\n");
/*	for(i=0;i<WINDOW_SIZE;i++)
	{
		for(j=0;j<WINDOW_SIZE;j++)
			printf("%.0f ", ldmap[i][j]);

		printf("\n");

		free(ldmap[i]);
	}
	free(ldmap);
*/

	//fflush(stdout);
	//assert(0);

	//printf(" [VAR %.0f %f ] ", muVarMaxLoc, muVarMax);
	//printf(" [SFS %.0f %f ] ", muSfsMaxLoc, muSfsMax);
	//printf(" [LD  %.0f %f ] ", muLdMaxLoc, muLdMax);
	//printf(" [MU  %.0f %f ] ", muMaxLoc, muMax);

	//MuVar_Accum += DIST (RSDMuStat->muVarMaxLoc, 50000.0);
 	//MuSfs_Accum += DIST (RSDMuStat->muSfsMaxLoc, 50000.0);
 	//MuLd_Accum += DIST (RSDMuStat->muLdMaxLoc, 50000.0);
 	//Mu_Accum += DIST (RSDMuStat->muMaxLoc, 50000.0);
	
	//printf("last i index %d\n", i);

	//printf("\n\nEnd of scanning chunk %d\n", RSDChunk->chunkID);
}
