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

int getGTLocation_vcf (char * string)
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
			if(!strcmp(tstring, "GT"))
				return Location;
			
			Location++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(tstring, &string[startIndex], (size_t)endIndex);
	tstring[endIndex]='\0';
	if(!strcmp(tstring, "GT"))
		return Location;

	return -1;	
}

void getGTData_vcf (char * string, int location, char * data) // from sample
{
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
				return;
			
			tLocation++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(data, &string[startIndex], (size_t)endIndex);
	data[endIndex]='\0';
	if(tLocation==location)
		return;	

	assert(0);
	return;
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

int getGTAlleles_vcf (char * string, char * stateVector, int statesTotal, char * sampleData, int * derivedAlleleCount, int * totalAlleleCount)
{	
	int i, j=0, index=0, start=0, end=0, len = (int)strlen(string), skipSNP=0;

	assert(stateVector!=NULL);

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

int maf_check (int ac, int at, double maf)
{
	if(ac<=0||at<=0)
		return 0;

	int check = 1; // default: do not discard

	double aaf = ((double)ac)/((double)at);
	double adf = ((double)at-ac)/((double)at);

	if(aaf<maf || adf<maf)
		check = 0;
	
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
}

void RSD_printMemory (FILE * fp1, FILE * fp2)
{
	fprintf(fp1, " Total memory footprint %.0f kbytes\n", MemoryFootprint/1024.0);
	fprintf(fp1, "\n");

	fprintf(fp2, " Total memory footprint %.0f kbytes\n", MemoryFootprint/1024.0);
	fprintf(fp2, "\n");
}

