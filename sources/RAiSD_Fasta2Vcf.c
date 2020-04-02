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

RSDFastaStates_t * 	RSDFastaStates_new 				(void);
void 			RSDFastaStates_reset 				(RSDFastaStates_t * RSDFastaStates);
int 			RSDFasteStates_matchState 			(RSDFastaStates_t * RSDFastaStates, char state);
void 			RSDFastaStates_addState 			(RSDFastaStates_t * RSDFastaStates, char state, int index);
void 			RSDFastaStates_free 				(RSDFastaStates_t * RSDFastaStates);
void 			RSDFastaStates_getStates 			(RSDFastaStates_t * RSDFastaStates, char * site, int64_t size);
void 			RSDFastaStates_print 				(RSDFastaStates_t * RSDFastaStates, int pos, char outgroupState, char altState, FILE * fpOut);
void 			RSDFastaStates_printImpute 			(RSDFastaStates_t * RSDFastaStates, int pos, char outgroupState, FILE * fpOut);
void 			RSDFastaStates_filterAmbiguousIngroup 		(RSDFastaStates_t * RSDFastaStates, char * tsite, int64_t numberOfSamples);
void			RSDFastaStates_filterAmbiguousOutgroup 		(char * outgroupSequence, int64_t sequenceLength);
void 			RSDFastaStates_filterAmbiguousOutgroup2 	(char * outgroupSequence, int64_t sequenceLength);
void 			RSDFastaStates_gaps2UN 				(char * tsite, int numberOfSamples);
void 			RSDFastaStates_resetImpute 			(RSDFastaStates_t * RSDFastaStates);
void 			RSDFastaStates_imputeProbsEqualize 		(RSDFastaStates_t * RSDFastaStates);
void 			RSDFastaStates_imputeStatesCompute 		(RSDFastaStates_t * RSDFastaStates);
char 			RSDFastaStates_getNewState 			(RSDFastaStates_t * RSDFastaStates);
void 			RSDFastaStates_imputeUN 			(RSDFastaStates_t * RSDFastaStates, char * tsite, int numberOfSamples);
void 			RSDFastaStates_imputeGAP 			(RSDFastaStates_t * RSDFastaStates, char * tsite, int numberOfSamples);
void 			RSDFastaStates_imputeOutgroupUN 		(RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int posIndex, int mode);
void 			RSDFastaStates_imputeOutgroupGAP 		(RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int posIndex, int mode);
int 			RSDFastaStates_matchOutgroupStateIngroup	(RSDFastaStates_t * RSDFastaStates, char outgroupState);
void 			RSDFastaStates_noninformativeOutgroup 		(RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int posIndex, int mode);
void 			RSDFastaStates_fixOutgroup 			(RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int pos);
char 			RSDFastaStates_getREF				(RSDFastaStates_t * RSDFastaStates, char outState);
char 			RSDFastaStates_getALT				(RSDFastaStates_t * RSDFastaStates, char refState);
int 			isACGT 						(char input);
int 			isAmbiguousDNACharacter				(char input);
int 			isValidDNACharacter				(char input);
char 			getRandACGT 					(void);
int 			seqOver 					(char state);
int 			skipSequence 					(FILE * fp);
void 			checkSequence 					(FILE * fp, double * percUN, double * percGAP, double * percACGT);
char 			getStateVCF					(char input, char refState, char altState);

char			rsd2upper					(char in);

int isACGT (char input)
{
	if(input==AD || input==ad ||
	   input==CY || input==cy ||
	   input==GU || input==gu ||
	   input==TH || input==th) return 1;

	return 0;	
}

int isAmbiguousDNACharacter(char input)
{
	if(input=='X' || input=='K' || 
	   input=='M' || input=='R' || 
	   input=='Y' || input=='S' || 
	   input=='W' || input=='B' || 
	   input=='V' || input=='H' || 
	   input=='D' || input=='.') return 1;

	if(input=='x' || input=='k' || 
	   input=='m' || input=='r' || 
	   input=='y' || input=='s' || 
	   input=='w' || input=='b' || 
	   input=='v' || input=='h' || 
	   input=='d' || input=='.') return 1;

	return 0;
}

int isValidDNACharacter(char input)
{
	if(input==AD || input==ad ||
	   input==CY || input==cy ||
	   input==GU || input==gu ||
	   input==TH || input==th ||
	   input==UN || input==un ||
	   input==GAP) return 1;

	return isAmbiguousDNACharacter(input);
} 

char getRandACGT (void)
{
	char acgt[4]={'A', 'C', 'G', 'T'};
	int rval = rand()%4;
	return acgt[rval];
}


int seqOver (char state)
{
	if(state=='>' || state==EOF)
		return 1;

	return 0;
}

int skipSequence (FILE * fp)
{
	assert(fp!=NULL);
	
	char tchar = (char)fgetc(fp);
	
	while(!seqOver(tchar))
	{
		tchar = (char)fgetc(fp);
		
		if(tchar==EOF)
			return EOF;
	}

	if(tchar=='>')
		ungetc(tchar, fp);
	
	return 0;
}

inline char rsd2upper (char in)
{
	char out = in;
	
	if(in>=97 && in <=122)
		out = in-32;

	return out;
}

void checkSequence (FILE * fp, double * percUN, double * percGAP, double * percACGT)
{
	assert(fp!=NULL);
	
	char tchar = (char)fgetc(fp);

	*percUN = 0.0;
	*percGAP = 0.0;
	*percACGT = 0.0;

	double totalLength = 0;

	while(!seqOver(tchar))
	{
		tchar = (char)fgetc(fp);
		
		if(tchar==EOF)
		{	
			//(*percACGT) = (*percACGT)/totalLength; 
			//(*percUN) = (*percUN)/totalLength; 
			//(*percGAP) = (*percGAP)/totalLength;
 			
			return;			
		}	

		if(isACGT(tchar))
			(*percACGT)+=1.0;

		if(isAmbiguousDNACharacter(tchar) || tchar==UN || tchar==un)
			(*percUN)+=1.0;

		if(tchar==GAP)
			(*percGAP)+=1.0;

		totalLength += 1.0;
	}

	if(tchar=='>')
	{
		ungetc(tchar, fp);
		totalLength -= 1.0;

		//(*percACGT) = (*percACGT)/totalLength; 
		//(*percUN) = (*percUN)/totalLength; 
		//(*percGAP) = (*percGAP)/totalLength;
	}
	
	return;
}

char getStateVCF(char input, char refState, char altState)
{

	// change1
	
	/**if(input==refState )
		return '0';

	if(input==altState || isACGT(input))
		return '1';
**/
	if(input==altState)
		return '1';

	if(input==refState || isACGT(input))
		return '0';

	if(input==UN)
		return '.';

	if(input==GAP)
		return '.';

	printf("Error: %c (%c %c)\n", input, refState, altState);
	fflush(stdout);

	assert(0);
	return '.';
}

RSDFastaStates_t * RSDFastaStates_new (void)
{
	RSDFastaStates_t * RSDFastaStates = (RSDFastaStates_t *)rsd_malloc(sizeof(RSDFastaStates_t));
	assert(RSDFastaStates!=NULL);

	RSDFastaStates->statesMax = 0;
	RSDFastaStates->statesTotal = 0;
	RSDFastaStates->statesChar = NULL;
	RSDFastaStates->statesCount = NULL;

	RSDFastaStates->imputeStatesMax = 0;
	RSDFastaStates->imputeStatesTotal = 0;
	RSDFastaStates->imputeStatesChar = NULL;
	RSDFastaStates->imputeStatesProb = NULL;

	return RSDFastaStates;
}

void RSDFastaStates_reset (RSDFastaStates_t * RSDFastaStates)
{
	assert(RSDFastaStates!=NULL);
	
	int i;
	for(i=0;i<RSDFastaStates->statesMax;i++)
	{
		RSDFastaStates->statesChar[i] = '/';
		RSDFastaStates->statesCount[i] = 0;
	}
	RSDFastaStates->statesTotal = 0;
}

int RSDFasteStates_matchState (RSDFastaStates_t * RSDFastaStates, char state)
{
	assert(RSDFastaStates!=NULL);
	
	int i;
	for(i=0;i<RSDFastaStates->statesTotal;i++)
	{
		if(state==RSDFastaStates->statesChar[i])
		{
			return i;
		}
	}
	return -1;
}

void RSDFastaStates_addState (RSDFastaStates_t * RSDFastaStates, char state, int index)
{
	assert(RSDFastaStates!=NULL);
	
	if(index==-1)
	{
		RSDFastaStates->statesTotal++;

		if(RSDFastaStates->statesTotal>RSDFastaStates->statesMax)
		{
			RSDFastaStates->statesChar = (char *)rsd_realloc(RSDFastaStates->statesChar, sizeof(char)*(unsigned long)RSDFastaStates->statesTotal);
			assert(RSDFastaStates->statesChar!=NULL);

			RSDFastaStates->statesCount = (int *)rsd_realloc(RSDFastaStates->statesCount, sizeof(int)*(unsigned long)RSDFastaStates->statesTotal);
			assert(RSDFastaStates->statesCount!=NULL);

			RSDFastaStates->statesMax = RSDFastaStates->statesTotal;
		}
		RSDFastaStates->statesChar[RSDFastaStates->statesTotal-1] = state;
		RSDFastaStates->statesCount[RSDFastaStates->statesTotal-1] = 1;		
	}
	else
	{
		assert(index<RSDFastaStates->statesTotal);
		assert(index<RSDFastaStates->statesMax);
		assert(state==RSDFastaStates->statesChar[index]);

		RSDFastaStates->statesCount[index]++;
	}

}

void RSDFastaStates_free (RSDFastaStates_t * RSDFastaStates)
{
	assert(RSDFastaStates!=NULL);

	if(RSDFastaStates->statesChar!=NULL)
		free(RSDFastaStates->statesChar);
	
	if(RSDFastaStates->statesCount!=NULL)
		free(RSDFastaStates->statesCount);

	if(RSDFastaStates->imputeStatesChar!=NULL)
		free(RSDFastaStates->imputeStatesChar);
	
	if(RSDFastaStates->imputeStatesProb!=NULL)
		free(RSDFastaStates->imputeStatesProb);

	free(RSDFastaStates);	
}

void RSDFastaStates_getStates (RSDFastaStates_t * RSDFastaStates, char * site, int64_t size)
{
	assert(RSDFastaStates!=NULL);
	assert(site!=NULL);
	assert(RSDFastaStates->statesTotal==0);

	int i=0, ind=-1;	
	
	for(i=0;i<size;i++)
	{
		ind=RSDFasteStates_matchState (RSDFastaStates, site[i]);
		RSDFastaStates_addState(RSDFastaStates, site[i], ind);	
	}

	int totalCount = 0;
	for(i=0;i<RSDFastaStates->statesTotal;i++)
		totalCount += RSDFastaStates->statesCount[i];

	assert(totalCount==size);
}

void RSDFastaStates_print (RSDFastaStates_t * RSDFastaStates, int pos, char outgroupState, char altState, FILE * fpOut)
{
	assert(RSDFastaStates!=NULL);
	assert(fpOut!=NULL);

	int i;
	fprintf(fpOut, "[%d][%c][%c]:", pos, outgroupState, altState);

	for(i=0;i<RSDFastaStates->statesTotal;i++)
	{
		fprintf(fpOut, "[\'%c\',%d]", RSDFastaStates->statesChar[i], RSDFastaStates->statesCount[i]);
	}
	fprintf(fpOut, "\n");
	fflush(fpOut);
}

void RSDFastaStates_printImpute (RSDFastaStates_t * RSDFastaStates, int pos, char outgroupState, FILE * fpOut)
{
	assert(RSDFastaStates!=NULL);
	assert(fpOut!=NULL);

	int i;
	fprintf(fpOut, "[%d][%c]:", pos, outgroupState);

	for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
	{
		fprintf(fpOut, "[\'%c\',%f]", RSDFastaStates->imputeStatesChar[i], RSDFastaStates->imputeStatesProb[i]);
	}
	fprintf(fpOut, "\n");
	fflush(fpOut);
}

void RSDFastaStates_filterAmbiguousIngroup (RSDFastaStates_t * RSDFastaStates, char * tsite, int64_t numberOfSamples)
{
	assert(RSDFastaStates!=NULL);
	assert(tsite!=NULL);

	int i;
	for(i=0;i<numberOfSamples;i++)
	{
		if(isAmbiguousDNACharacter(tsite[i])==1)
			tsite[i] = UN;
	}

	RSDFastaStates_reset (RSDFastaStates);
	RSDFastaStates_getStates (RSDFastaStates, tsite, numberOfSamples);
}

void RSDFastaStates_filterAmbiguousOutgroup (char * outgroupSequence, int64_t sequenceLength)
{
	assert(outgroupSequence!=NULL);

	int i;
	for(i=0;i<sequenceLength;i++)
	{
		if(isAmbiguousDNACharacter(outgroupSequence[i])==1)
			outgroupSequence[i] = UN;
	}		
}

void RSDFastaStates_filterAmbiguousOutgroup2 (char * outgroupSequence, int64_t sequenceLength)
{
	if(outgroupSequence==NULL)
		return;

	int i;
	for(i=0;i<sequenceLength;i++)
	{
		if(isAmbiguousDNACharacter(outgroupSequence[i])==1)
			outgroupSequence[i] = UN;
	}		
}

void RSDFastaStates_gaps2UN (char * tsite, int numberOfSamples)
{
	return;

	assert(tsite!=NULL);

	int i;
	for(i=0;i<numberOfSamples;i++)
	{
		if(tsite[i]==GAP)
			tsite[i] = UN;
	}
}

void RSDFastaStates_resetImpute (RSDFastaStates_t * RSDFastaStates)
{
	assert(RSDFastaStates!=NULL);

	int i;
	for(i=0;i<RSDFastaStates->imputeStatesMax;i++)
	{
		RSDFastaStates->imputeStatesChar[i] = '/';
		RSDFastaStates->imputeStatesProb[i] = 0.0;
	}
	RSDFastaStates->imputeStatesTotal = 0;	
}

void RSDFastaStates_imputeProbsEqualize (RSDFastaStates_t * RSDFastaStates)
{
	assert(RSDFastaStates!=NULL);
	
	if(RSDFastaStates->imputeStatesTotal==0) 
	{
		if(RSDFastaStates->imputeStatesMax<4)
		{
			RSDFastaStates->imputeStatesMax = 4;

			RSDFastaStates->imputeStatesChar = (char *)rsd_realloc(RSDFastaStates->imputeStatesChar, sizeof(char)*(unsigned long)RSDFastaStates->imputeStatesMax);
			assert(RSDFastaStates->imputeStatesChar!=NULL);

			RSDFastaStates->imputeStatesProb = (double *)rsd_realloc(RSDFastaStates->imputeStatesProb, sizeof(double)*(unsigned long)RSDFastaStates->imputeStatesMax);
			assert(RSDFastaStates->imputeStatesProb!=NULL);
		}

		RSDFastaStates_resetImpute (RSDFastaStates);

		RSDFastaStates->imputeStatesTotal = 4;

		RSDFastaStates->imputeStatesChar[0] = AD;
		RSDFastaStates->imputeStatesChar[1] = TH;
		RSDFastaStates->imputeStatesChar[2] = GU;
		RSDFastaStates->imputeStatesChar[3] = CY;

		RSDFastaStates->imputeStatesProb[0] = 0.25;
		RSDFastaStates->imputeStatesProb[1] = 0.25;
		RSDFastaStates->imputeStatesProb[2] = 0.25;
		RSDFastaStates->imputeStatesProb[3] = 0.25;
	}	
}

void RSDFastaStates_imputeStatesCompute (RSDFastaStates_t * RSDFastaStates)
{
	assert(RSDFastaStates!=NULL);
	assert(RSDFastaStates->imputeStatesTotal==0);

	int i;
	double totalACGT=0;

	for(i=0;i<RSDFastaStates->statesTotal;i++)
	{
		if(isACGT(RSDFastaStates->statesChar[i])==1)
		{
			RSDFastaStates->imputeStatesTotal++;

			if(RSDFastaStates->imputeStatesTotal>RSDFastaStates->imputeStatesMax)
			{
				RSDFastaStates->imputeStatesChar = (char *)rsd_realloc(RSDFastaStates->imputeStatesChar, sizeof(char)*(unsigned long)RSDFastaStates->imputeStatesTotal);
				assert(RSDFastaStates->imputeStatesChar!=NULL);

				RSDFastaStates->imputeStatesProb = (double *)rsd_realloc(RSDFastaStates->imputeStatesProb, sizeof(double)*(unsigned long)RSDFastaStates->imputeStatesTotal);
				assert(RSDFastaStates->imputeStatesProb!=NULL);

				RSDFastaStates->imputeStatesMax = RSDFastaStates->imputeStatesTotal;
			}
			RSDFastaStates->imputeStatesChar[RSDFastaStates->imputeStatesTotal-1] = RSDFastaStates->statesChar[i];
			RSDFastaStates->imputeStatesProb[RSDFastaStates->imputeStatesTotal-1] = (double)RSDFastaStates->statesCount[i];
			totalACGT += (double)RSDFastaStates->statesCount[i];
		}			
	}

	for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
		RSDFastaStates->imputeStatesProb[i] /= totalACGT;

	RSDFastaStates_imputeProbsEqualize (RSDFastaStates);	
}

char RSDFastaStates_getNewState (RSDFastaStates_t * RSDFastaStates)
{
	assert(RSDFastaStates!=NULL);
	assert(RSDFastaStates->imputeStatesTotal!=0);

	int randT = rand();
	double randVal = ((double)randT)/((double)RAND_MAX);

	int i;
	double low = 0.0, high = 0.0;
	for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
	{
		high = low + RSDFastaStates->imputeStatesProb[i];
		
		if(randVal>=low && randVal<high)
			return RSDFastaStates->imputeStatesChar[i];

		low = high;
	}	

	return RSDFastaStates->imputeStatesChar[RSDFastaStates->imputeStatesTotal-1];
}

void RSDFastaStates_imputeUN (RSDFastaStates_t * RSDFastaStates, char * tsite, int numberOfSamples)
{
	assert(RSDFastaStates!=NULL);
	assert(tsite!=NULL);

	RSDFastaStates_resetImpute (RSDFastaStates);
	RSDFastaStates_imputeStatesCompute (RSDFastaStates);

	int i;
	for(i=0;i<numberOfSamples;i++)
	{
		if(tsite[i]==UN)
		{
			tsite[i] = RSDFastaStates_getNewState (RSDFastaStates);
			tsite[i] = getRandACGT();
		}
	}	

	RSDFastaStates_reset (RSDFastaStates);
	RSDFastaStates_getStates (RSDFastaStates, tsite, numberOfSamples);
}

void RSDFastaStates_imputeGAP (RSDFastaStates_t * RSDFastaStates, char * tsite, int numberOfSamples)
{
	assert(RSDFastaStates!=NULL);
	assert(tsite!=NULL);

	RSDFastaStates_resetImpute (RSDFastaStates);
	RSDFastaStates_imputeStatesCompute (RSDFastaStates);

	int i;
	for(i=0;i<numberOfSamples;i++)
	{
		if(tsite[i]==GAP)
		{
			tsite[i] = RSDFastaStates_getNewState (RSDFastaStates);
		}
	}
	RSDFastaStates_reset (RSDFastaStates);
	RSDFastaStates_getStates (RSDFastaStates, tsite, numberOfSamples);	
}

void RSDFastaStates_imputeOutgroupUN (RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int posIndex, int mode)
{
	assert(RSDFastaStates!=NULL);
	assert(RSDDataset!=NULL);
	assert(posIndex>=0);
	assert(mode==MAJORITY || mode==PROBABILITY);

	if(RSDDataset->outgroupSequence[posIndex]!=UN)
		return;

	RSDFastaStates_resetImpute (RSDFastaStates);
	RSDFastaStates_imputeStatesCompute (RSDFastaStates);

	int i;
	if(mode==MAJORITY)
	{
		int maxIndex = -1;
		double maxVal = 0.0;
		for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
		{
			if(RSDFastaStates->imputeStatesProb[i]>maxVal)
			{
				maxVal = RSDFastaStates->imputeStatesProb[i];
				maxIndex = i;
			}
		}
	
		RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[maxIndex];
		return;		
	}

	if(mode==PROBABILITY)
	{
		int randT = rand();
		double randVal = ((double)randT)/((double)RAND_MAX);

		double low = 0.0, high = 0.0;
		for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
		{
			high = low + RSDFastaStates->imputeStatesProb[i];
		
			if(randVal>=low && randVal<high)
			{
				RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[i];
				return;
			}

			low = high;
		}	
		RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[RSDFastaStates->imputeStatesTotal-1];
		return;	
	}

	return;
}

void RSDFastaStates_imputeOutgroupGAP (RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int posIndex, int mode)
{
	assert(RSDFastaStates!=NULL);
	assert(RSDDataset!=NULL);
	assert(posIndex>=0);
	assert(mode==MAJORITY || mode==PROBABILITY);

	if(RSDDataset->outgroupSequence[posIndex]!=GAP)
		return;	

	RSDFastaStates_resetImpute (RSDFastaStates);
	RSDFastaStates_imputeStatesCompute (RSDFastaStates);

	//RSDFastaStates_printImpute (RSDFastaStates, 0, 'x', stdout);

	int i;
	if(mode==MAJORITY)
	{
		int maxIndex = -1;
		double maxVal = 0.0;
		for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
		{
			if(RSDFastaStates->imputeStatesProb[i]>maxVal)
			{
				maxVal = RSDFastaStates->imputeStatesProb[i];
				maxIndex = i;
			}
		}
		RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[maxIndex];
		return;		
	}

	if(mode==PROBABILITY)
	{
		int randT = rand();
		double randVal = ((double)randT)/((double)RAND_MAX);

		double low = 0.0, high = 0.0;
		for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
		{
			high = low + RSDFastaStates->imputeStatesProb[i];
		
			if(randVal>=low && randVal<high)
			{
				RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[i];
				return;
			}

			low = high;
		}	
		RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[RSDFastaStates->imputeStatesTotal-1];
		return;	
	}
	return;
}

int RSDFastaStates_matchOutgroupStateIngroup (RSDFastaStates_t * RSDFastaStates, char outgroupState)
{
	assert(RSDFastaStates!=NULL);

	int i;
	
	for(i=0;i<RSDFastaStates->statesTotal;i++)
	{
		if(outgroupState==RSDFastaStates->statesChar[i])
			return 1;
	}
	
	return 0;
}

void RSDFastaStates_noninformativeOutgroup (RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int posIndex, int mode)
{
	assert(RSDFastaStates!=NULL);
	assert(RSDDataset!=NULL);
	assert(posIndex>=0);
	assert(mode==MAJORITY || mode==PROBABILITY);

	char outgroupState = RSDDataset->outgroupSequence[posIndex];

	if(outgroupState!=GAP && outgroupState!=UN && RSDFastaStates_matchOutgroupStateIngroup (RSDFastaStates, outgroupState)==1) 
		return;

	if(RSDDataset->outgroupSequence2!=NULL)
	{
		outgroupState = RSDDataset->outgroupSequence2[posIndex];

		if(outgroupState!=GAP && outgroupState!=UN && RSDFastaStates_matchOutgroupStateIngroup (RSDFastaStates, outgroupState)==1)
		{ 
			RSDDataset->outgroupSequence[posIndex] = outgroupState;
			return;
		}
	}

	RSDFastaStates_resetImpute (RSDFastaStates);
	RSDFastaStates_imputeStatesCompute (RSDFastaStates);

	int i;
	if(mode==MAJORITY)
	{
		int maxIndex = -1;
		double maxVal = 0.0;
		for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
		{
			if(RSDFastaStates->imputeStatesProb[i]>maxVal)
			{
				maxVal = RSDFastaStates->imputeStatesProb[i];
				maxIndex = i;
			}
		}
		RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[maxIndex];
		return;		
	}

	if(mode==PROBABILITY)
	{
		int randT = rand();
		double randVal = ((double)randT)/((double)RAND_MAX);

		double low = 0.0, high = 0.0;
		for(i=0;i<RSDFastaStates->imputeStatesTotal;i++)
		{
			high = low + RSDFastaStates->imputeStatesProb[i];
		
			if(randVal>=low && randVal<high)
			{
				RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[i];
				return;
			}

			low = high;
		}	
		RSDDataset->outgroupSequence[posIndex] = RSDFastaStates->imputeStatesChar[RSDFastaStates->imputeStatesTotal-1];
		return;	
	}
	return;
}

void RSDFastaStates_fixOutgroup (RSDFastaStates_t * RSDFastaStates, RSDDataset_t * RSDDataset, int pos)
{
	assert(RSDFastaStates!=NULL);
	assert(RSDDataset!=NULL);
	assert(pos>=0);

	int i;

	if(isACGT(RSDDataset->outgroupSequence[pos])) // or gap
	{		
		for(i=0;i<RSDFastaStates->statesTotal;i++)
		{
			if(RSDFastaStates->statesChar[i]==RSDDataset->outgroupSequence[pos])
				return;
		}

		RSDFastaStates_resetImpute (RSDFastaStates);
		RSDFastaStates_imputeStatesCompute (RSDFastaStates);

		RSDDataset->outgroupSequence[pos] = RSDFastaStates->imputeStatesChar[rand()%RSDFastaStates->imputeStatesTotal]; // this is wrong rand draw from the valid states
	}

	assert(isACGT(RSDDataset->outgroupSequence[pos]));
}

char RSDFastaStates_getREF(RSDFastaStates_t * RSDFastaStates, char outState)
{
	assert(RSDFastaStates!=NULL);

	char retOutState = outState;
	int i=0, acgtStates = 0;
	if(!isACGT(retOutState))
	{
		for(i=0;i<RSDFastaStates->statesTotal;i++)
		{
			if(isACGT(RSDFastaStates->statesChar[i]))
				acgtStates++;
		}
		
		if(acgtStates!=0)
		{
			int sval=-1, rval = rand() % acgtStates;
			for(i=0;i<RSDFastaStates->statesTotal;i++)
			{
				if(isACGT(RSDFastaStates->statesChar[i]))
				{
					sval++;
					if(sval==rval)
						retOutState = RSDFastaStates->statesChar[i];
				}
			}
		}
		else
		{
			retOutState = getRandACGT();
		}			
	}	

	return retOutState;
}

char RSDFastaStates_getALT(RSDFastaStates_t * RSDFastaStates, char refState)
{
	assert(RSDFastaStates!=NULL);
		
	int i=0;

	if(isACGT(refState))
	{
		for(i=0;i<RSDFastaStates->statesTotal;i++)
		{
			if(RSDFastaStates->statesChar[i]!=refState && isACGT(RSDFastaStates->statesChar[i]))
				return RSDFastaStates->statesChar[i];
		}

		char randState = getRandACGT ();
		while(randState==refState)
			randState = getRandACGT ();

		return randState;
	}

	assert(0);	

	return '.';
}

void RSDDataset_convertFasta2VCF (RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut)
{
	assert(RSDDataset!=NULL);
	assert(RSDCommandLine!=NULL);
	assert(fpOut!=NULL);
	
	FILE * fp = fopen(RSDCommandLine->inputFileName, "r");
	assert(fp!=NULL);

	char fileNameNew[STRING_SIZE];
	strncpy(fileNameNew, RSDCommandLine->inputFileName, STRING_SIZE);
	strcat(fileNameNew, ".vcf");

	fprintf(fpOut, "\nMESSAGE: Converting input file %s to vcf (%s)", RSDCommandLine->inputFileName, fileNameNew);
	fflush(fpOut);

	fprintf(stdout, "\nMESSAGE: Converting input file %s to vcf (%s)", RSDCommandLine->inputFileName, fileNameNew);
	fflush(stdout);

	FILE * fpNew = fopen(fileNameNew, "r");

	if(fpNew!=NULL) // fxd vcf file already exists
	{
		if(RSDCommandLine->overwriteOutput==0)
		{
			fprintf(fpOut, "\n\nERROR: Converted file %s exists.\n        Use -f to overwrite it or give it as input with -I.\n\n", fileNameNew);
			fflush(fpOut);

			fprintf(stderr, "\n\nERROR: Converted file %s exists.\n       Use -f to overwrite it or give it as input with -I.\n\n", fileNameNew);
			exit(0);
		}
		else
		{
			fclose(fpNew);
		}
	}

	fpNew = fopen(fileNameNew, "w");
	assert(fpNew!=NULL);

	fprintf(fpNew, "##fileformat=VCFv4.2\n");
	fprintf(fpNew, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

	char tstring[STRING_SIZE], seqName[STRING_SIZE];
	tstring[0]='\0';
	seqName[0]='\0';
	int64_t numberOfSamples = 0, sequenceLength = 0, sequenceLength2 = 0;

	int i=-1, j=0, outgroup=0, firstSeqOutgroup=0, outgroup2=0; 

	fpos_t outseqPosition;
	fpos_t outseqPosition2;
	fpos_t * seqPosition = NULL;

	strcpy(RSDDataset->outgroupName, RSDCommandLine->outgroupName);
	strcpy(RSDDataset->outgroupName2, RSDCommandLine->outgroupName2);
	
	int secondaryOutgroup=0;

	if(!strcmp(RSDDataset->outgroupName,"\0")) // if primary outgroup not provided
	{
		if(!strcmp(RSDDataset->outgroupName2,"\0")) // if secondary outgroup not provided
		{
			firstSeqOutgroup = 1;
		}
		else
		{
			strcpy(RSDDataset->outgroupName, RSDDataset->outgroupName2);
			strcpy(RSDDataset->outgroupName2, "\0");
		}
	}
	else
	{
		if(strcmp(RSDDataset->outgroupName2,"\0")) // if secondary outgroup is given
			secondaryOutgroup=1;
	}

	int rcnt = fscanf(fp, "%s", tstring);
	assert(rcnt==1);

	fgetpos(fp, &outseqPosition2); // just to make the compiler happy, will either be overwritten or never used			

	// check if outgroups exist and set ptrs
	while(rcnt!=EOF)
	{
		strcpy(seqName, &tstring[1]);

		if(firstSeqOutgroup==1 && !strcmp(RSDDataset->outgroupName, "\0"))
		{
			strcpy(RSDDataset->outgroupName, seqName);
			assert(secondaryOutgroup==0);
		}
	
		if(!strcmp(RSDDataset->outgroupName, seqName))
		{
			outgroup=1;

			tstring[0] = (char)fgetc(fp); // \n newline
			assert(tstring[0]=='\n');

			fgetpos(fp, &outseqPosition);
		}

		if(!strcmp(RSDDataset->outgroupName2, seqName))
		{
			outgroup2=1;

			tstring[0] = (char)fgetc(fp); // \n newline
			assert(tstring[0]=='\n');

			fgetpos(fp, &outseqPosition2);			
		}
		
		rcnt = skipSequence(fp);

		if(rcnt!=EOF)
			rcnt = fscanf(fp, "%s", tstring);
	}

	rewind(fp);

	if(outgroup==1) // primary outgroup found
	{
		if(outgroup2==0) // secondary outgroup not found 
		{
			fprintf(fpOut, "\nMESSAGE: Using sequence %s as outgroup", RSDDataset->outgroupName);
			fflush(fpOut);

			fprintf(stdout, "\nMESSAGE: Using sequence %s as outgroup", RSDDataset->outgroupName);
			fflush(stdout);

			if(secondaryOutgroup==1)
			{
				fprintf(fpOut, "\nMESSAGE: Secondary outgroup sequence %s not found", RSDDataset->outgroupName2);
				fflush(fpOut);

				fprintf(stdout, "\nMESSAGE: Secondary outgroup sequence %s not found", RSDDataset->outgroupName2);
				fflush(stdout);

				strcpy(RSDDataset->outgroupName2, "\0");
			}
		}
		else // secondary outgroup found
		{
			fprintf(fpOut, "\nMESSAGE: Using sequence %s as primary outgroup", RSDDataset->outgroupName);
			fprintf(fpOut, "\nMESSAGE: Using sequence %s as secondary outgroup", RSDDataset->outgroupName2);
			fflush(fpOut);

			fprintf(stdout, "\nMESSAGE: Using sequence %s as primary outgroup", RSDDataset->outgroupName);
			fprintf(stdout, "\nMESSAGE: Using sequence %s as secondary outgroup", RSDDataset->outgroupName2);
			fflush(stdout);
		}
	}
	else // primary outgroup not found
	{
		fprintf(fpOut, "\nMESSAGE: Outgroup sequence %s not found", RSDDataset->outgroupName);
		fflush(fpOut);

		fprintf(stdout, "\nMESSAGE: Outgroup sequence %s not found", RSDDataset->outgroupName);
		fflush(stdout);
		
		if(outgroup2==0) // secondary outgroup not found
		{
			if(secondaryOutgroup==1)
			{
				fprintf(fpOut, "\nMESSAGE: Outgroup sequence %s not found", RSDDataset->outgroupName2);
				fflush(fpOut);

				fprintf(stdout, "\nMESSAGE: Outgroup sequence %s not found", RSDDataset->outgroupName2);
				fflush(stdout);	

				strcpy(RSDDataset->outgroupName2, "\0");	
			}

			rcnt = fscanf(fp, "%s", tstring);
			assert(rcnt==1);

			strcpy(seqName, &tstring[1]);

			strcpy(RSDDataset->outgroupName, seqName);
			firstSeqOutgroup = 1;
		
			tstring[0] = (char)fgetc(fp); // \n newline
			assert(tstring[0]=='\n');

			fgetpos(fp, &outseqPosition);

			fprintf(fpOut, "\nMESSAGE: Using sequence %s as outgroup", RSDDataset->outgroupName);
			fflush(fpOut);

			fprintf(stdout, "\nMESSAGE: Using sequence %s as outgroup", RSDDataset->outgroupName);
			fflush(stdout);
		}
		else // secondary outgroup found
		{
			fprintf(fpOut, "\nMESSAGE: Using only secondary outgroup sequence %s", RSDDataset->outgroupName2);
			fflush(fpOut);

			fprintf(stdout, "\nMESSAGE: Using only secondary outgroup sequence %s", RSDDataset->outgroupName2);
			fflush(stdout);	

			outseqPosition = outseqPosition2;

			strcpy(RSDDataset->outgroupName, RSDDataset->outgroupName2);
			strcpy(RSDDataset->outgroupName2, "\0");
		}
	}
	
	assert(strcmp(RSDDataset->outgroupName, "\0"));

	// fetch primary outgroup sequence
	fsetpos(fp, &outseqPosition);

	tstring[0] = (char)(fgetc(fp));

	while(!seqOver(tstring[0]))
	{
		if(tstring[0]!='\n')
		{
			tstring[0] = rsd2upper(tstring[0]);

			sequenceLength++;

			RSDDataset->outgroupSequence = rsd_realloc(RSDDataset->outgroupSequence, sizeof(char)*(unsigned long)(sequenceLength+1));
			assert(RSDDataset->outgroupSequence!=NULL);

			RSDDataset->outgroupSequence[sequenceLength-1] = tstring[0];
			RSDDataset->outgroupSequence[sequenceLength] = '\0';				

			tstring[0] = (char)(fgetc(fp));

			if(tstring[0]!='\n' && !isValidDNACharacter(tstring[0]))
			{
				fprintf(fpOut, "\n\nERROR: Invalid character (%c) found in sequence %s (outgroup)\n\n", tstring[0], RSDDataset->outgroupName);
				fflush(fpOut);

				fprintf(stderr, "\n\nERROR: Invalid character (%c) found in sequence %s (outgroup)\n\n", tstring[0], RSDDataset->outgroupName);
				exit(0);
			}
		}
		else
			tstring[0] = (char)(fgetc(fp));	
	}

	assert(sequenceLength!=0); // proceed based on expected number of sites equal to the outgroup size
	assert(RSDDataset->outgroupSequence!=NULL);

	if(strcmp(RSDDataset->outgroupName2, "\0"))
	{
		rewind(fp);

		// fetch secondary outgroup sequence
		fsetpos(fp, &outseqPosition2);

		tstring[0] = (char)(fgetc(fp));

		while(!seqOver(tstring[0]))
		{
			if(tstring[0]!='\n')
			{
				tstring[0] = rsd2upper(tstring[0]);

				sequenceLength2++;

				RSDDataset->outgroupSequence2 = rsd_realloc(RSDDataset->outgroupSequence2, sizeof(char)*(unsigned long)(sequenceLength2+1));
				assert(RSDDataset->outgroupSequence2!=NULL);

				RSDDataset->outgroupSequence2[sequenceLength2-1] = tstring[0];
				RSDDataset->outgroupSequence2[sequenceLength2] = '\0';				

				tstring[0] = (char)(fgetc(fp));

				if(tstring[0]!='\n' && !isValidDNACharacter(tstring[0]))
				{
					fprintf(fpOut, "\n\nERROR: Invalid character (%c) found in sequence %s (outgroup)\n\n", tstring[0], RSDDataset->outgroupName2);
					fflush(fpOut);

					fprintf(stderr, "\n\nERROR: Invalid character (%c) found in sequence %s (outgroup)\n\n", tstring[0], RSDDataset->outgroupName2);
					exit(0);
				}
			}
			else
				tstring[0] = (char)(fgetc(fp));	
		}
	
		assert(sequenceLength2!=0);
		assert(sequenceLength==sequenceLength2);
		assert(RSDDataset->outgroupSequence2!=NULL);
	}

	// get sequence pointers
	rewind(fp);

	rcnt = fscanf(fp, "%s", tstring);
	assert(rcnt==1);

	//char fileNameNew2[STRING_SIZE];
	//strncpy(fileNameNew2, RSDCommandLine->inputFileName, STRING_SIZE);
	//strcat(fileNameNew2, ".seqQual");

	while(rcnt!=EOF)
	{
		strcpy(seqName, &tstring[1]);

		tstring[0] = (char)fgetc(fp); // \n newline
		assert(tstring[0]=='\n');
		
		outgroup=0;
		if(!strcmp(RSDDataset->outgroupName, seqName))
			outgroup=1;

		outgroup2=0;
		if(!strcmp(RSDDataset->outgroupName2, seqName))
			outgroup2=1;

		if((outgroup==0 && outgroup2==0) || (outgroup==1 && firstSeqOutgroup==1) ) // if its ingroup
		{
			numberOfSamples++;
				
			seqPosition = (fpos_t*)rsd_realloc(seqPosition, sizeof(fpos_t)*((unsigned long)numberOfSamples));
			assert(seqPosition!=NULL);

			fgetpos(fp, &seqPosition[numberOfSamples-1]);

			fprintf(fpNew, "\t%s", seqName);

			// add sequence filtering here
			//double percUN = 0.0, percGAP=0.0, percACGT=0.0;
			//checkSequence (fp, &percUN, &percGAP, &percACGT);

			rcnt = skipSequence(fp);			
		}
		else // skip outgroup sequence without any processing only when outgroup is given by the user
			rcnt = skipSequence(fp);

		if(rcnt!=EOF)
			rcnt = fscanf(fp, "%s", tstring);
	}

	fprintf(fpNew, "\n");

	assert(numberOfSamples!=0);
	assert(seqPosition!=NULL);

	// fetch data in column-major order
	//char fileNameSiteReport[STRING_SIZE];
	//strncpy(fileNameSiteReport, RSDCommandLine->inputFileName, STRING_SIZE);
	//strcat(fileNameSiteReport, ".siteReport");

	//FILE * fpSiteReport = fopen(fileNameSiteReport, "w");
	//assert(fpSiteReport!=NULL);

	char * tsite = (char*)rsd_malloc(sizeof(char)*(unsigned long)(numberOfSamples+1));
	assert(tsite!=NULL);

	RSDFastaStates_t * RSDFastaStates = RSDFastaStates_new ();

	for(j=0;j<sequenceLength;j++)	
	{
		for(i=0;i<numberOfSamples;i++)		
		{
			fsetpos(fp, &seqPosition[i]);

			tstring[0] = rsd2upper((char)fgetc(fp));

			while(tstring[0]=='\n' || tstring[0]=='\t' || tstring[0]==' ')
				tstring[0] = rsd2upper((char)fgetc(fp));

			if(!isValidDNACharacter(tstring[0]))
			{
				fprintf(fpOut, "\n\nERROR: Invalid character (%c) found in sequence %s\n\n", tstring[0], seqName);
				fflush(fpOut);

				fprintf(stderr, "\n\nERROR: Invalid character (%c) found in sequence %s\n\n", tstring[0], seqName);
				exit(0);
			}

			tsite[i] = tstring[0]; 
			tsite[i+1] = '\0';
			fgetpos(fp, &seqPosition[i]);
		}

		RSDFastaStates_reset (RSDFastaStates);
		RSDFastaStates_getStates (RSDFastaStates, tsite, numberOfSamples); // get all states at site, including gaps, Ns, and ambiguous characters.
		
		//RSDFastaStates_print (RSDFastaStates, j+1, RSDDataset->outgroupSequence[j],'X', fpSiteReport);

		RSDFastaStates_filterAmbiguousIngroup (RSDFastaStates, tsite, numberOfSamples); // ingroup ambiguous -> N
		RSDFastaStates_filterAmbiguousOutgroup (RSDDataset->outgroupSequence, sequenceLength); // outgroup ambiguous -> N
		RSDFastaStates_filterAmbiguousOutgroup2 (RSDDataset->outgroupSequence2, sequenceLength); // outgroup2 ambiguous -> N

		RSDFastaStates_noninformativeOutgroup (RSDFastaStates, RSDDataset, j, MAJORITY);

		//RSDFastaStates_imputeUN (RSDFastaStates, tsite, numberOfSamples); // probability-based impute N
		//RSDFastaStates_imputeGAP (RSDFastaStates, tsite, numberOfSamples); // probability-based impute GAP
		//RSDFastaStates_print (RSDFastaStates, j+1,  RSDDataset->outgroupSequence[j],'X', fpSiteReport);

		char refState = RSDDataset->outgroupSequence[j]; 
		char altState = RSDFastaStates_getALT(RSDFastaStates, refState);

		//RSDFastaStates_print (RSDFastaStates, j+1, refState, altState, fpSiteReport);

		fprintf(fpNew, "%s\t", RSDCommandLine->chromNameVCF);
		fprintf(fpNew, "%d\t", j+1);
		fprintf(fpNew, ".\t");
		fprintf(fpNew, "%c\t", refState);
		fprintf(fpNew, "%c\t", altState);
		fprintf(fpNew, ".\tPASS\t.\tGT");

		for(i=0;i<numberOfSamples;i++)
		{
			tsite[i] = getStateVCF(tsite[i], refState, altState);
			fprintf(fpNew, "\t%c", tsite[i]);
		}
		fprintf(fpNew, "\n");	

		RSDFastaStates_reset (RSDFastaStates);
		RSDFastaStates_getStates (RSDFastaStates, tsite, numberOfSamples);

		//RSDFastaStates_print (RSDFastaStates, j+1, refState, altState, fpSiteReport);
		//fprintf(fpSiteReport, "\n");
	}

	if(RSDCommandLine->fasta2vcfMode==FASTA2VCF_CONVERT_n_PROCESS)
	{
		fprintf(fpOut, "\nMESSAGE: Processing continues using file %s", fileNameNew);
		fflush(fpOut);

		fprintf(stdout, "\nMESSAGE: Processing continues using file %s", fileNameNew);
		fflush(stdout);

		strncpy(RSDCommandLine->inputFileName, fileNameNew, STRING_SIZE);
		strcpy(RSDDataset->inputFileFormat, "vcf");
	}

	fprintf(fpOut, "\n\n");
	fprintf(stdout, "\n\n");

	if(tsite!=NULL)
		free(tsite);

	if(seqPosition!=NULL)
		free(seqPosition);

	RSDFastaStates_free (RSDFastaStates);
	
	fclose(fp);
	fclose(fpNew);
	//fclose(fpSiteReport);

	if(RSDCommandLine->fasta2vcfMode==FASTA2VCF_CONVERT_n_EXIT)
	{
		fprintf(fpOut, " Conversion completed!\n VCF output file %s generated!\n\n", fileNameNew);
		fflush(fpOut);

		fprintf(stdout, " Conversion completed!\n VCF output file %s generated!\n\n", fileNameNew);
		fflush(stdout);

		RSDCommandLine_free(RSDCommandLine);
		RSDDataset_free(RSDDataset);

		if(RAiSD_Info_FP!=NULL)
			fclose(RAiSD_Info_FP);

		if(RAiSD_SiteReport_FP!=NULL)
			fclose(RAiSD_SiteReport_FP);

		if(RAiSD_ReportList_FP!=NULL)
			fclose(RAiSD_ReportList_FP);

		exit(0);
	}
}
