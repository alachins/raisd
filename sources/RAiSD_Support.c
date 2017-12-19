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

unsigned long long rdtsc(void)
{
	unsigned a, d;

	__asm__ volatile("rdtsc" : "=a" (a), "=d" (d));

	return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

double gettime(void)
{
	struct timeval ttime;
	gettimeofday(&ttime , NULL);
	return ttime.tv_sec + ttime.tv_usec * 0.000001;
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

	return count;
}

void popcount_u64_init (void)
{
	unsigned int i;    
	for (i = 0; i < (0x1u<<16); i++)
		POPCNT_U16_LUT[i] = popcount_u32_iterative(i);
	
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
		{	return 1;
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
	uint8_t temp = ~B[size-1];
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
	int i, length = strlen(string);
	int startIndex = 0;
	int endIndex = 0;
	for(i=0;i<length;i++)
	{
		if((string[i]==':'))
		{
			endIndex =  i;
			memcpy(tstring, &string[startIndex], endIndex-startIndex);
			tstring[endIndex]='\0';
			if(!strcmp(tstring, "GT"))
				return Location;
			
			Location++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(tstring, &string[startIndex], endIndex);
	tstring[endIndex]='\0';
	if(!strcmp(tstring, "GT"))
		return Location;

	return -1;	
}

void getGTData_vcf (char * string, int location, char * data) // from sample
{
	int tLocation = 0;
	int i, length = strlen(string);
	int startIndex = 0;
	int endIndex = 0;
	for(i=0;i<length;i++)
	{
		if((string[i]==':'))
		{
			endIndex =  i;
			memcpy(data, &string[startIndex], endIndex-startIndex);
			data[endIndex]='\0';

			if(tLocation==location)
				return;
			
			tLocation++;
			startIndex = i+1;
		}
	}
	endIndex =  i;
	memcpy(data, &string[startIndex], endIndex);
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
		//printf("%d\n", index);

		tmp = data[index];
		data[index] = data[i];
		data[i] = tmp;
	}
}

void getGTAlleles_vcf (char * string, char * stateVector, int statesTotal, char * sampleData, int * alleleCount)
{	
	int i, j=0, index=0, start=0, end=0, len = strlen(string);

	for(i=0;i<len;i++)
	{	
		if(string[i]>=48 && string[i]<=57)
		{
			index = string[i]-48;

			if(index>0)
				(*alleleCount)++;
			
			assert(index<statesTotal);
			
			//FSM
			//sampleData[j++] = stateVector[index]; 
			assert(stateVector!=NULL);

			//ISM
			if(string[i]=='0')
				sampleData[j++] = '0';
			else
				sampleData[j++] = '1';

			sampleData[j] = '\0';
		}
		else
		{
  			if(string[i]=='.')
			{
				sampleData[j++] = 'N';
				sampleData[j] = '\0';
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
}


float * putInSortVector(int * size, float * vector, float value)
{
	(*size)++;
	vector = realloc(vector, sizeof(float)*(*size));
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



