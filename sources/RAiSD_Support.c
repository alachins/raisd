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

void	ignoreLineSpaces	(FILE *fp, char *ent);
int 	flagMatch		(FILE *fp, char flag[], int flaglength, char tmp);
void 	RSD_printTime 		(FILE * fp1, FILE * fp2);
void 	RSD_printMemory 	(FILE * fp1, FILE * fp2);
int 	RSD_Rscript_check	(void);
void 	RSD_Rscript_generate 	(void);
void 	RSD_Rscript_remove 	(void);

char POPCNT_U16_LUT [0x1u << 16];

unsigned long long rdtsc(void)
{
	unsigned a, d;

	__asm__ volatile("rdtsc" : "=a" (a), "=d" (d));

	return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

#ifndef _INTRINSIC_POPCOUNT
int popcount_u32_iterative (unsigned int n)
{
	int count=0;    
    
	while(n)
	{
		count += n & 0x1u ;    
		n >>= 1 ;
	}
	
	assert(count<=16);
	assert(count>=0);

	return count;
}

void popcount_u64_init (void)
{
	unsigned int i;    
	for (i = 0; i < (0x1u<<16); i++)
	{
		int t = popcount_u32_iterative(i);
		POPCNT_U16_LUT[i] = (char)t;
	}	
	
}
#endif

inline int rsd_popcnt_u64 (uint64_t input)
{
#ifndef _INTRINSIC_POPCOUNT
	return POPCNT_U16_LUT[input & 0xffffu] + POPCNT_U16_LUT[ (input>>16) & 0xffffu] + POPCNT_U16_LUT[ (input>>32) & 0xffffu] + POPCNT_U16_LUT[ (input>>48) & 0xffffu]; 
#else
	return _mm_popcnt_u64(input);
#endif
}

int snpv_cmp (uint64_t * A, uint64_t * B, int size)
{
	assert(A!=NULL);
	assert(B!=NULL);
	assert(size>=1);
	
	int i;
	for(i=0;i!=size;++i)
		if(A[i]!=B[i])
		{	
			return 1;
		}
	
	return 0;
}

int snpv_cmp_cross_masks (uint64_t * A, uint64_t * B, uint64_t * mA, uint64_t * mB, int size)
{
	assert(A!=NULL);
	assert(B!=NULL);
	assert(size>=1);
	
	int i;
	for(i=0;i!=size;++i)
		if((A[i]&mB[i])!=(B[i]&mA[i]))
		{	
			return 1;
		}

	return 0;
}

int isnpv_cmp_cross_masks (uint64_t * A, uint64_t * B, uint64_t * mA, uint64_t * mB, int size)
{	
	assert(A!=NULL);
	assert(B!=NULL);
	assert(size>=1);
	
	int i;
	for(i=0;i!=size;++i)
		if(((~A[i])&mA[i]&mB[i])!=(B[i]&mA[i]))
		{	
			return 1;
		}

	return 0;
}

int isnpv_cmp (uint64_t * A, uint64_t * B, int size, int numberOfSamples)
{
	assert(A!=NULL);
	assert(B!=NULL);
	assert(size>=1);
	
	int i;
	for(i=0;i!=size-1;++i)
	{
		if((~A[i])!=B[i])
		{
			return 1;
		}
	}
	uint64_t temp = ~B[size-1];
	int shiftLast = 64-(numberOfSamples-(size-1)*64);
	temp = temp << shiftLast;
	temp = temp >> shiftLast;
	if(temp!=A[size-1])
		return 1;

	return 0;
}

int isnpv_cmp_with_mask (uint64_t * A, uint64_t * B, uint64_t * mA, uint64_t * mB, int size, int numberOfSamples)
{
	assert(A!=NULL);
	assert(B!=NULL);
	assert(size>=1);
	
	int i;
	for(i=0;i!=size-1;++i)
	{
		if(((~A[i])&mA[i])!=B[i])
		{
			return 1;
		}
	}

	uint64_t temp = (~B[size-1])&mB[size-1];
	int shiftLast = 64-(numberOfSamples-(size-1)*64);
	temp = temp << shiftLast;
	temp = temp >> shiftLast;
	if(temp!=A[size-1])
		return 1;

	return 0;
}

int getGXLocation_vcf (char * string, char * GX)
{
	assert(string!=NULL);
	
	char tstring [STRING_SIZE];
	int Location = 0;
	int i, length = (int)strlen(string);
	int startIndex = 0;
	int endIndex = 0;
	for(i=0;i<length;i++)
	{
		if(string[i]==':')
		{
			endIndex =  i;
			memcpy(tstring, &string[startIndex], (size_t)(endIndex-startIndex));
			tstring[endIndex]='\0';
			if(!strcmp(tstring, GX))
				return Location;
			
			Location++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(tstring, &string[startIndex], (size_t)endIndex);
	tstring[endIndex]='\0';
	if(!strcmp(tstring, GX))
		return Location;

	return -1;	
}

int diploidyCheck(char * data)
{
	unsigned int i, count=0;	
	for(i=0;i<strlen(data);i++)
		count += (data[i]==','? 1:0);
	
	int check = count==2?1:0;
	return check;
}

void getGPProbs (char * data, double *p00, double *p01, double * p11, int isLik)
{
	char tstring[STRING_SIZE];
	unsigned int i;
	int tstr_i=0;
	
	for(i=0;i<strlen(data);i++)
	{
		if(data[i]==',')
		{
			tstring[tstr_i] = '\0';
			tstr_i=0;

			if(*p00<0.0)
				*p00 = atof(tstring);
			else
				*p01 = atof(tstring);			
			
		}
		else
			tstring[tstr_i++] = data[i];
	}

	tstring[tstr_i] = '\0';
	*p11 = atof(tstring);	

	if(isLik==0)
		assert(*p00>=0.0 && *p01>=0.0 && *p11>=0.0);
	else
		assert(*p00<=0.0 && *p01<=0.0 && *p11<=0.0);
}

void reconGT (char * data)
{
	assert(diploidyCheck(data)==1);

	double p00 = -1.0, p01 = -1.0, p11 = -1.0;
	getGPProbs(data, &p00, &p01, &p11, 0);

	double checksum = p00+p01+p11;
	assert(checksum>=0.999 && checksum<=1.001);

	int val = rand();
	double rval = ((double)val) / ((double)RAND_MAX);
	
	if(rval<p00)
		strcpy (data, "0/0");
	else
	{
		if(rval>=p00 && rval<p00+p01)
			strcpy (data, "0/1");
		else
			strcpy (data, "1/1");	
	}
}

int getGXData_vcf (char * string, int location, char * data)
{
	assert(string!=NULL);
	assert(location>=0);
	assert(data!=NULL);	

	int tLocation = 0;
	int i, length = (int)strlen(string);
	int startIndex = 0;
	int endIndex = 0;
	for(i=0;i<length;i++)
	{
		if(string[i]==':')
		{
			endIndex =  i;
			memcpy(data, &string[startIndex], (size_t)(endIndex-startIndex));
			data[endIndex]='\0';

			if(tLocation==location)
				return 1;

			tLocation++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(data, &string[startIndex], (size_t)endIndex);
	data[endIndex]='\0';
	if(tLocation==location)
		return 1;

	return 0;		
}

void getGTData_vcf (char * string, int locationGT, int locationGP, int locationGL, char * data) // from sample
{
	assert(locationGT!=-1 || (locationGP!=-1 && locationGL!=-1));

	if(locationGT!=-1)
	{
		int ret = getGXData_vcf(string, locationGT, data);
		assert(ret==1);
	}
	else
	{
		int ret = getGXData_vcf(string, locationGL, data);
		assert(ret==1);

		double p00 = 0.0, p01 = 0.0, p11 = 0.0; // likelihoods
		getGPProbs(data, &p00, &p01, &p11, 1); 

		if(p00+p01+p11!=0.0)
		{
			ret = getGXData_vcf(string, locationGP, data);
			assert(ret==1);
		
			reconGT (data);	
		}
		else
			strcpy (data, "./.");
			
	}
}

void dataShuffleKnuth(char * data, int startIndex, int endIndex)
{
	if(startIndex == endIndex)
		return;

	int i, index;
	char tmp;

	for (i = endIndex; i > startIndex; i--)
	{
		index = startIndex + (rand() % (i - startIndex + 1));

		tmp = data[index];
		data[index] = data[i];
		data[i] = tmp;
	}
}

int getGTAlleles_vcf (char * string, char * stateVector, int statesTotal, char * sampleData, int * derivedAlleleCount, int * totalAlleleCount, int ploidy)
{	
	assert(statesTotal>=2);
	assert(stateVector!=NULL);

	int i, j=0, index=0, start=0, end=0, len = (int)strlen(string), skipSNP=0;
	
	for(i=0;i<len;i++)
	{	
		if(string[i]>=48 && string[i]<=57)
		{
			if(string[i]!='0' && string[i]!='1')
			{
				fprintf(stderr, "\n ERROR: Invalid character (%c) found!\n\n", string[i]);
				exit(0);
			}

			index = string[i]-48;
	
			assert(index==0 || index==1);

			(*totalAlleleCount)++;

			(*derivedAlleleCount)+=index;

			assert(index<statesTotal);

			sampleData[j++] = string[i];			

			sampleData[j] = '\0';
		}
		else
		{
			if(string[i]=='.')
			{
				sampleData[j++] = 'N';
				sampleData[j] = '\0';
				skipSNP=1;
	
				if(ploidy>1 && len==1)
				{
					assert(!strcmp(string, "."));

					int cnt = ploidy-1;
					while(cnt--!=0)
					{
						sampleData[j++] = 'N';
						sampleData[j] = '\0';
					}
				}
			}

			if(string[i]=='/')
			{
				end++;
			}

			if(string[i]=='|')
			{
				dataShuffleKnuth(sampleData, start, end);
				start = j;
				end = j;
			}			
		}
	}

	dataShuffleKnuth(sampleData, start, end);

	return skipSNP;
}


float * putInSortVector(int * size, float * vector, float value)
{
	(*size)++;
	vector = realloc(vector, sizeof(float)*((size_t)(*size)));
	vector[(*size)-1] = 0.0f;

	if(*size==1)
	{
		vector[0] = value;
		return vector;
	}

	int i;
	for(i=0;i<*size;i++)
	{
		if(value >= vector[i])
		{
			int pos = i;
			for(i=*size-1;i>pos;i--)
				vector[i] = vector[i-1];

			vector[pos] = value;
			i = *size + 1;
		}		
	}

	return vector;
}

float DIST (float a, float b)
{
	if(a>=b)
		return a-b;
	else
		return b-a;
}

char alleleMask_binary (char c, int * isDerived, int * isValid, FILE * fpOut)
{
	*isDerived = 0;
	*isValid = 0;

	switch(c)
	{
		case '0':
			*isDerived = 0;
			*isValid = 1;
			return c;

		case '1':
			*isDerived = 1;
			*isValid = 1;
			return c;

		/*case '-':
			*isDerived = 0;
			*isValid = 1;
			return 'N';

		case 'N':
			*isDerived = 0;
			*isValid = 1;
			return 'N';*/

		default:
			fprintf(fpOut, "ERROR: Unrecognized character %c\n\n",c);
			fprintf(stderr, "ERROR: Unrecognized character %c\n\n",c);
			exit(0);

	}	
}

int monomorphic_check (int incomingSiteDerivedAlleleCount, int setSamples, int64_t * cnt, int skipSNP)
{
	if(skipSNP==1) // to avoid double counting
		return 0;

	int check = 1;

	if(incomingSiteDerivedAlleleCount==0 || incomingSiteDerivedAlleleCount==setSamples)
	{
		check = 0;
 		++(*cnt);
	}

	return check;	
}

int strictPolymorphic_check (int incomingSiteDerivedAlleleCount, int incomingSiteTotalAlleleCount, int64_t * cnt, int skipSNP)
{
	if(skipSNP==1) // to avoid double counting
		return 0;

	int check = 1;

	if(incomingSiteDerivedAlleleCount==incomingSiteTotalAlleleCount)
	{
		check = 0;
		++(*cnt);
	}

	return check;	
}

int maf_check (int ac, int at, double maf, int64_t * cnt, int skipSNP)
{
	if(skipSNP==1) // to avoid double counting
		return 0;

	if(ac<=0||at<=0)
		return 0;

	int check = 1; // default: do not discard

	double aaf = ((double)ac)/((double)at);
	double adf = ((double)at-ac)/((double)at);

	if(aaf<maf || adf<maf)
	{
		check = 0;
		++(*cnt);
	}

	return check;
}

void ignoreLineSpaces(FILE *fp, char *ent)
{
	while(*ent==' '|| *ent == 9) // horizontal tab
		*ent = (char)fgetc(fp);  
}

int flagMatch(FILE *fp, char flag[], int flaglength, char tmp)
{
	int counter = 0;
	while(counter < flaglength)
	{
		if(tmp != flag[counter])
		  {
		    break;
		  }
		tmp = (char)fgetc(fp);

		++counter;
	}
	
	return counter;
}

void RSD_printTime (FILE * fp1, FILE * fp2)
{
	clock_gettime(CLOCK_REALTIME, &requestEnd);
	double TotalTime = (requestEnd.tv_sec-requestStart.tv_sec)+ (requestEnd.tv_nsec-requestStart.tv_nsec) / BILLION;


	fprintf(fp1, "\n\n");
	fprintf(fp1, " Total execution time %.5f seconds\n", TotalTime);

	fprintf(fp2, "\n\n");
	fprintf(fp2, " Total execution time %.5f seconds\n", TotalTime);

#ifdef _PTIMES
	fprintf(fp1, " Total OoC time %.5f seconds\n", TotalOoCTime);
	fprintf(fp2, " Total OoC time %.5f seconds\n", TotalOoCTime);
#endif
}

void RSD_printMemory (FILE * fp1, FILE * fp2)
{
	fprintf(fp1, " Total memory footprint %.0f kbytes\n", MemoryFootprint/1024.0);
	fprintf(fp1, "\n");

	fprintf(fp2, " Total memory footprint %.0f kbytes\n", MemoryFootprint/1024.0);
	fprintf(fp2, "\n");
}

void RSD_printSiteReportLegend (FILE * fp, int64_t imputePerSNP, int64_t createPatternPoolMask)
{
	if(fp==NULL)
		return;

	if(imputePerSNP==0 && createPatternPoolMask==0) // M=0
		fprintf(fp, "\n Index: Name | Sites = SNPs + Discarded | Discarded = HeaderCheckFailed + MAFCheckFailed + WithMissing + Monomorphic\n");
	else
	{
		if(imputePerSNP==1) // M=1
			fprintf(fp, "\n Index: Name | Sites = SNPs + Discarded | Discarded = HeaderCheckFailed + MAFCheckFailed + PotentiallyMonomorphicSites | Imputed\n");
		else // M=2,3
			fprintf(fp, "\n Index: Name | Sites = SNPs + Discarded | Discarded = HeaderCheckFailed + MAFCheckFailed + PotentiallyMonomorphicSites\n");
	}
	fflush(fp);
}

#ifdef _ZLIB
int gzscanf (gzFile fp, char * string)
{
	assert(string!=NULL);
	char t = (char) gzgetc(fp);

	while(t==' ' || t=='\n' || t==EOF || t==13 || t=='\t')
	{
		if(t==EOF)
			return 0;

		t = (char) gzgetc(fp);
	}

	int i=0;
	string[0]='\0';
	while(t!=' ' && t!='\n' && t!=EOF && t!=13 && t!='\t')
	{
		string[i++]=t;
		//assert(t!=' ' && t!='\n' && t!=EOF && t!=13 && t!='\t');
		t = (char) gzgetc(fp);
	}
	gzungetc(t,fp);
	//assert(t==' ' || t=='\n' || t==EOF || t==13 || t=='\t');

	string[i]='\0';

	return 1;
}
#endif

