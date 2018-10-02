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

char	(*RSDDataset_goToNextSet) 		(RSDDataset_t * RSDDataset);
int	(*RSDDataset_getNumberOfSamples) 	(RSDDataset_t * RSDDataset);
int 	(*RSDDataset_getValidSampleList) 	(RSDDataset_t * RSDDataset);
int 	(*RSDDataset_getFirstSNP)		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
int 	(*RSDDataset_getNextSNP)		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
int 	RSDDataset_getFirstSNP_vcf 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
int 	RSDDataset_getNumberOfSamples_vcf 	(RSDDataset_t * RSDDataset);
int 	RSDDataset_getValidSampleSize 		(RSDDataset_t * RSDDataset);
void 	RSDDataset_validateSample		(RSDDataset_t * RSDDataset, char * sampleName, int sampleIndex);
int 	RSDDataset_getValidSampleList_vcf 	(RSDDataset_t * RSDDataset);
char 	RSDDataset_goToNextSet_vcf 		(RSDDataset_t * RSDDataset);
int 	RSDDataset_getFirstSNP_ms 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
int 	RSDDataset_getNumberOfSamples_ms 	(RSDDataset_t * RSDDataset);
int 	RSDDataset_getValidSampleList_ms 	(RSDDataset_t * RSDDataset);
char 	RSDDataset_goToNextSet_ms 		(RSDDataset_t * RSDDataset);
void 	RSDDataset_reportMissing 		(RSDDataset_t * RSDDataset, FILE * fpOut);
void 	RSDDataset_detectFormat			(RSDDataset_t * RSDDataset);


/***********/
/* Dataset */
/***********/

RSDDataset_t * RSDDataset_new(void)
{
	RSDDataset_t * d = NULL;

	d = (RSDDataset_t *) malloc(sizeof(RSDDataset_t));
	assert(d!=NULL);

	d->inputFilePtr = NULL;
	strcpy(d->inputFileFormat, "invalid");
	d->inputFileIsMBS = 0;
	d->numberOfSamples = -1;
	strcpy(d->setID, "-1");
	d->setSamples = 0;
	d->setSize = 0;
	d->setSNPs = 0;
	d->setProgress = -1;
	d->setParsingMode = MULTI_STEP_PARSING;
	d->setRegionLength = 0ull;
	d->sampleFilePtr = NULL;
	d->sampleValidListSize = ALL_SAMPLES_VALID;
	d->sampleValidList = NULL;
	d->sampleValid = NULL;
	d->numberOfSamplesVCF = 0;
	d->setSitesDiscarded = 0;
	d->setSitesDiscardedHeaderCheckFailed=0;
	d->setSitesDiscardedWithMissing = 0;
	d->setSitesDiscardedMafCheckFailed = 0;
	d->setSitesDiscardedStrictPolymorphicCheckFailed = 0;
	d->setSitesDiscardedMonomorphic = 0;
	d->setSitesImputedTotal = 0;

	return d;
}

void RSDDataset_free(RSDDataset_t * d)
{
	assert(d!=NULL);
	
	if(d->inputFilePtr!=NULL)
		fclose(d->inputFilePtr);

	if(d->sampleValidList!=NULL)
	{
		free(d->sampleValidList);
	}

	if(d->sampleValid!=NULL)
		free(d->sampleValid);

	free(d);
}

void RSDDataset_detectFormat(RSDDataset_t * RSDDataset)
{
	FILE * fp = RSDDataset->inputFilePtr;
	
	char tmp;

	int cnt = 0;

	tmp=(char)fgetc(fp);

	char macsFLAG[] = "COMMAND";

	int flaglength = 7,j;

	char vcfFLAG[] = "##fileformat=VCF";

	int vcf_flaglength = 16;

	char sfFLAG[]="position";

	int sf_flaglength = 8;	  

	int format = -999;

	while(tmp!=EOF)
	{

	  
		if(tmp=='/')
		{
			tmp = (char)fgetc(fp);			

			if(tmp=='/' && cnt < 3)
			{
			    	format = MS_FORMAT;
				strcpy(RSDDataset->inputFileFormat, "ms");
			    	break; 
		  	}
			else
	  			tmp = (char)fgetc(fp);				
		}
		else
		{
			if(tmp=='>')
		  	{
			    	format = FASTA_FORMAT;
				strcpy(RSDDataset->inputFileFormat, "fasta");
				assert(0);
			    	break;
		  	}
			else
			{
			  	j = flagMatch(fp, macsFLAG, flaglength, tmp);
			  
			  	if(j == flaglength)
			    	{
					strcpy(RSDDataset->inputFileFormat, "macs");
					format = MACS_FORMAT;
					assert(0);
					break;
				}
			  else
			    {
			      // unset
			      
			      fseek(fp, -j-1, SEEK_CUR);

			      tmp = (char)fgetc(fp);
			      
			      j = flagMatch(fp, vcfFLAG, vcf_flaglength, tmp);
			      
			      if(j == vcf_flaglength)
				{
					strcpy(RSDDataset->inputFileFormat, "vcf");
					format = VCF_FORMAT;
					break;
				}
			      else
				{
				  // unset
				  
				  fseek(fp, -j-1, SEEK_CUR);
				  
				  tmp = (char)fgetc( fp );
				  
				  j = flagMatch(fp, sfFLAG, sf_flaglength, tmp);
				  
				  if( j == sf_flaglength)
				    {
				      tmp = (char)fgetc(fp);
				      
				      ignoreLineSpaces(fp, &tmp);
				      
				      int j1 = flagMatch(fp, "x", 1, tmp);
				      
				      if( j1 == 1)
					{
					  tmp = (char)fgetc( fp );
					  
					  ignoreLineSpaces(fp, &tmp);

					  int j2 = flagMatch(fp, "n", 1, tmp);

					  if( j2 == 1)
					    {	
						strcpy(RSDDataset->inputFileFormat, "sf");				      
					      	format = SF_FORMAT;
					      	break;
					    }
					}
				    
				    }
				  else
				    {
				      tmp = (char)fgetc(fp);
				      
				      if(tmp == 10 || tmp == 13)
					cnt = 0;
				      
				      if(tmp != ' ')
					cnt++;
				      
				    }
				}
			    }
			}
		}		
	}
	
	if(format == -999)
	  	format = OTHER_FORMAT;
}

void RSDDataset_init (RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut)
{
	assert(RSDDataset!=NULL);
	assert(RSDCommandLine!=NULL);

	RSDDataset->inputFilePtr = fopen(RSDCommandLine->inputFileName, "r");
	assert(RSDDataset->inputFilePtr!=NULL);

	fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));

	RSDDataset_detectFormat(RSDDataset);

	fclose(RSDDataset->inputFilePtr);

	// Last command line check for length value, required with ms files.
	if(RSDCommandLine->regionLength==0ull && !strcmp(RSDDataset->inputFileFormat, "ms"))
	{
		fprintf(stderr, "\nERROR: Missing required input parameter -L\n\n");
		exit(0);
	}

	RSDDataset->inputFilePtr = fopen(RSDCommandLine->inputFileName, "r");
	assert(RSDDataset->inputFilePtr!=NULL);

	if(RSDCommandLine->printSampleList==1)
	{
		char sampleFileName[STRING_SIZE];
		strcpy(sampleFileName, "RAiSD_SampleList.");
		strcat(sampleFileName, RSDCommandLine->runName);
		strcat(sampleFileName, "\0");
	
		RSDDataset->sampleFilePtr = fopen(sampleFileName, "w");
		assert(RSDDataset->sampleFilePtr!=NULL);
	}
	else
	{
		if(strcmp(RSDCommandLine->sampleFileName, "\0"))
		{
			RSDDataset->sampleFilePtr = fopen(RSDCommandLine->sampleFileName, "r");
			assert(RSDDataset->sampleFilePtr!=NULL);
		}

	} 
	
	RSDDataset_initParser (RSDDataset, fpOut, RSDCommandLine);

	if(strcmp(RSDCommandLine->sampleFileName, "\0")) // gets in here when sample file is provided as input
		RSDDataset_getValidSampleList (RSDDataset);

	RSDDataset_getNumberOfSamples (RSDDataset);	

	if (RSDDataset->sampleFilePtr!=NULL)
	{
		fclose(RSDDataset->sampleFilePtr);
		RSDDataset->sampleFilePtr = NULL;

		if(RSDCommandLine->printSampleList==1)
			exit(0);
	}

	if(RSDDataset->numberOfSamples<MIN_NUMBER_OF_SAMPLES)
	{
		fprintf(stderr, "\nERROR: Number of samples (%d) less than minimum (%d)\n\n", RSDDataset->numberOfSamples, MIN_NUMBER_OF_SAMPLES);
		exit(0);
	}	
}

void RSDDataset_reportMissing (RSDDataset_t * RSDDataset, FILE * fpOut)
{
	int i;
	for(i=0;i<RSDDataset->sampleValidListSize;i++)
		if(strlen(RSDDataset->sampleValidList[i])!=0)
			fprintf(fpOut, "\nMISSING SAMPLE: %s", RSDDataset->sampleValidList[i]);	

}

void RSDDataset_print (RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut)
{
	assert(RSDDataset!=NULL);

	if(fpOut!=NULL)
	{
		fprintf(fpOut, " Samples: %d", RSDDataset->numberOfSamples);
			
		if(RSDDataset->sampleValidListSize!=ALL_SAMPLES_VALID)
		{
			fprintf(fpOut, " [Total: %d, Not found: %d, Requested %d]", RSDDataset->numberOfSamplesVCF, RSDDataset->sampleValidListSize-RSDDataset->numberOfSamples, RSDDataset->sampleValidListSize);	
		RSDDataset_reportMissing (RSDDataset, fpOut);	

		}

		fprintf(fpOut, "\n");
		if(!strcmp(RSDDataset->inputFileFormat, "ms"))
			fprintf(fpOut, " Region:  %lu bp\n", RSDCommandLine->regionLength);

		fprintf(fpOut, " Format:  %s\n", RSDDataset->inputFileFormat);

		fflush(fpOut);
	}

	if(!strcmp(RSDDataset->inputFileFormat, "invalid"))
	{
		fprintf(fpOut, "\nERROR: Unrecognized input file format\n\n");
		exit(0);
	}
	if(!strcmp(RSDDataset->inputFileFormat, "sf"))
	{
		fprintf(fpOut, "\nERROR: The SweepFinder file format is not supported\n\n");
		exit(0);
	}		
}

void RSDDataset_setPosition (RSDDataset_t * RSDDataset, int * setIndex)
{
	assert(RSDDataset!=NULL);

	(*setIndex)++;

	fsetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));

	RSDDataset->setSize = 0;
	RSDDataset->setSNPs = 0;
	RSDDataset->setProgress = 0;
	RSDDataset->setSamples = -1;
	RSDDataset->setParsingMode = MULTI_STEP_PARSING;
}


/*************/
/* ms format */
/*************/

char RSDDataset_goToNextSet_ms (RSDDataset_t * RSDDataset)
{
	assert(RSDDataset->inputFilePtr!=NULL);

	fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
	char tchar=(char)fgetc(RSDDataset->inputFilePtr);

	while(tchar!=EOF)
	{
		if(tchar=='/')
		{
			tchar=(char)fgetc(RSDDataset->inputFilePtr);
			if(tchar=='/')
			{
				int setID = atoi(RSDDataset->setID);
				setID++;
				sprintf(RSDDataset->setID, "%d", setID);
				return tchar;
			}				
		}

		fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
		tchar=(char)fgetc(RSDDataset->inputFilePtr);
	}

	return EOF;
}

int RSDDataset_getValidSampleList_ms (RSDDataset_t * RSDDataset)
{
	assert(RSDDataset!=NULL);

	return 0;
}

int RSDDataset_getNumberOfSamples_ms (RSDDataset_t * RSDDataset)
{
	char tstring[STRING_SIZE];
	int rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // ms/mssel/msHOT/command(if mbs)
	assert(rcnt==1);
	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // samples if !mbs
	assert(rcnt==1);
	
	if(RSDDataset->inputFileIsMBS==1)
	{
		rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // samples if mbs
		assert(rcnt==1);
	}	

	int val = atoi(tstring);
	RSDDataset->numberOfSamples = val;
	assert(RSDDataset->numberOfSamples>=1);
	
	return RSDDataset->numberOfSamples;
}

void RSDDataset_getSetRegionLength_ms (RSDDataset_t * RSDDataset, uint64_t length)
{
	RSDDataset->setRegionLength = length;
}

int RSDDataset_getFirstSNP_ms (RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut)
{
	RSDDataset_getSetRegionLength_ms (RSDDataset, length);

	// This function does not return until a SNP is found or the set is over.

	char tstring[STRING_SIZE];
	int setDone = 0;

	RSDPatternPool->incomingSitePosition = -1.0;

	fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));

	int rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	assert(rcnt==1); // //

	tstring[2] = '\0'; // for all ms-like keep just // (required for mbs)

	if(strcmp(tstring, "//"))
	{
		fsetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
		RSDDataset->setSize = 0;
		RSDDataset->setSNPs = 0;
		RSDDataset->setProgress = 0;
		RSDDataset->setSamples = 0;
		setDone = 1;
		return setDone; 
	}

	if(RSDDataset->inputFileIsMBS==1)
	{	
		tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
		while(tstring[0]!='\n')
			tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // skip line
	}

	fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	assert(rcnt==1); // segsites:

	if(strcmp(tstring, "segsites:"))
	{	
		fsetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
		RSDDataset->setSize = 0;
		RSDDataset->setSNPs = 0;
		RSDDataset->setProgress = 0;
		RSDDataset->setSamples = 0;
		setDone = 1;
		return setDone;
	}

	RSDDataset->setSamples = RSDDataset->numberOfSamples; // Always equal for ms. It can differ in vcf.

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	assert(rcnt==1); // number of sites

	int vali = atoi(tstring);
	RSDDataset->setSize = RSDDataset->inputFileIsMBS==1?vali-1:vali;
	RSDDataset->setProgress = 1; // first site loaded
	RSDDataset->setSNPs = 0; // Init to 0 since we dont know yet whether this is a polymorphic one or not

	if(vali<MIN_SET_SNPS)
	{	
		fsetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
		//RSDDataset->setSize = 0;
		RSDDataset->setSNPs = 0;
		RSDDataset->setProgress = 0;
		RSDDataset->setSamples = 0;
		setDone = 1;
		return setDone;
	}

	RSDDataset->setParsingMode = MULTI_STEP_PARSING;

	if(RSDDataset->setSize<=RSDPatternPool->maxSize)
		RSDDataset->setParsingMode = SINGLE_STEP_PARSING;

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	assert(rcnt==1); // positions:

	if(strcmp(tstring, "positions:"))
	{	
		fsetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
		RSDDataset->setSize = 0;
		RSDDataset->setSNPs = 0;
		RSDDataset->setProgress = 0;
		RSDDataset->setSamples = 0;
		setDone = 1;
		return setDone;
	}

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	assert(rcnt==1); // first SNP position

	double valf = atof(tstring);
	valf = RSDDataset->inputFileIsMBS==1?valf:valf*(double)length;
	//valf = valf * (double)length;
	RSDPatternPool->incomingSitePosition = valf;

	fgetpos(RSDDataset->inputFilePtr, &(RSDChunk->posPosition));

	int i;
	for(i=1;i<RSDDataset->setSize;i++) // parse/ignore remaining SNP positions
	{	
		rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
		assert(rcnt==1);
	}

	tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
	while(tstring[0]!='\n')
		tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // get to the end of the line

	RSDPatternPool->incomingSiteDerivedAlleleCount = 0;
	RSDPatternPool->incomingSiteTotalAlleleCount = 0;
	int drvAll_incr = 0, totAll_incr=0;
	for(i=0;i<RSDDataset->numberOfSamples;i++)
	{
		tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // load first bit on each line
		
		RSDPatternPool->incomingSite[i] = alleleMask_binary(tstring[0], &drvAll_incr, &totAll_incr, fpOut);
		
		RSDPatternPool->incomingSiteDerivedAlleleCount += drvAll_incr;

		RSDPatternPool->incomingSiteTotalAlleleCount += totAll_incr;

		fgetpos(RSDDataset->inputFilePtr, &(RSDChunk->seqPosition[i]));

		if(RSDDataset->setParsingMode == MULTI_STEP_PARSING)
		{
			while(tstring[0]!='\n')
				tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // skip line
		}
		else // SINGLE_STEP_PARSING
		{	
			int siteIndex = 0;
			while(tstring[0]!='\n')
			{	
				siteIndex++;
				tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // parse line

				if(tstring[0]!='\n')
				{
					uint64_t tval = RSDPatternPool->poolData[siteIndex*RSDPatternPool->patternSize+i/64];
					char c = alleleMask_binary(tstring[0], &drvAll_incr, &totAll_incr, fpOut);
					uint8_t b = (uint8_t)(c-48); 
					tval = (tval<<1) | b;
					RSDPatternPool->poolData[siteIndex*RSDPatternPool->patternSize+i/64] = tval;
				}
			}	
		}		
	}
	RSDPatternPool->incomingSite[RSDDataset->numberOfSamples] = '\0'; 

	if(RSDPatternPool->incomingSiteDerivedAlleleCount!=0 && RSDPatternPool->incomingSiteDerivedAlleleCount!=RSDDataset->numberOfSamples && maf_check(RSDPatternPool->incomingSiteDerivedAlleleCount, RSDPatternPool->incomingSiteTotalAlleleCount, maf, &RSDDataset->setSitesDiscardedMafCheckFailed, 0)==1)
		RSDDataset->setSNPs++;

	while((RSDPatternPool->incomingSiteDerivedAlleleCount==0||RSDPatternPool->incomingSiteDerivedAlleleCount==RSDDataset->numberOfSamples || maf_check(RSDPatternPool->incomingSiteDerivedAlleleCount, RSDPatternPool->incomingSiteTotalAlleleCount, maf, &RSDDataset->setSitesDiscardedMafCheckFailed, 0)!=1) && setDone==0) // keep loading until first SNP is found
		setDone = RSDDataset_getNextSNP_ms (RSDDataset, RSDPatternPool, RSDChunk, RSDCommandLine, length, maf, fpOut);

	return setDone;
}

int RSDDataset_getNextSNP_ms (RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut)
{
	// This function does not return until a SNP is found or the set is over.

	char tstring[STRING_SIZE];
	int setDone = 0;

	RSDPatternPool->incomingSitePosition = -1.0;

	fsetpos(RSDDataset->inputFilePtr, &(RSDChunk->posPosition));

	int rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	assert(rcnt==1);

	double valf = atof(tstring);
	valf = RSDDataset->inputFileIsMBS==1?valf:valf*(double)length;

	RSDPatternPool->incomingSitePosition = valf;

	fgetpos(RSDDataset->inputFilePtr, &(RSDChunk->posPosition));

	RSDDataset->setProgress++;

	if(RSDDataset->setProgress==RSDDataset->setSize)
		setDone = 1;

	int drvAll_incr = 0, totAll_incr=0;

	if(RSDDataset->setParsingMode == MULTI_STEP_PARSING)
	{
		int i;

		RSDPatternPool->incomingSiteDerivedAlleleCount = 0;
		RSDPatternPool->incomingSiteTotalAlleleCount = 0;
		for(i=0;i<RSDDataset->numberOfSamples;i++)
		{
			fsetpos(RSDDataset->inputFilePtr, &(RSDChunk->seqPosition[i]));

			tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);

			RSDPatternPool->incomingSite[i] = alleleMask_binary(tstring[0], &drvAll_incr, &totAll_incr, fpOut);
		
			RSDPatternPool->incomingSiteDerivedAlleleCount += drvAll_incr;

			RSDPatternPool->incomingSiteTotalAlleleCount += totAll_incr;

			fgetpos(RSDDataset->inputFilePtr, &(RSDChunk->seqPosition[i]));
		}
	}
	else // SINGLE_STEP_PARSING
	{
		int i;

		memcpy(RSDPatternPool->incomingSiteCompact, &(RSDPatternPool->poolData[(RSDDataset->setProgress-1)*RSDPatternPool->patternSize]), (unsigned long)(RSDPatternPool->patternSize*8));
		
		for(i=RSDDataset->numberOfSamples-1;i!=-1;--i)
		{
			uint8_t v = RSDPatternPool->incomingSiteCompact[i/64] & 1;
			RSDPatternPool->incomingSiteCompact[i/64] =  RSDPatternPool->incomingSiteCompact[i/64] >> 1;
			RSDPatternPool->incomingSite[i] = (char)(v + 48);
		}
		RSDPatternPool->incomingSite[RSDDataset->numberOfSamples] = '\0';

		RSDPatternPool->incomingSiteDerivedAlleleCount = 0;
		RSDPatternPool->incomingSiteTotalAlleleCount = 0;

		int tlen = (int)strlen(RSDPatternPool->incomingSite);
		for(i=0;i<tlen;i++)
		{
			RSDPatternPool->incomingSite[i] = alleleMask_binary(RSDPatternPool->incomingSite[i], &drvAll_incr, &totAll_incr, fpOut);
			RSDPatternPool->incomingSiteDerivedAlleleCount += drvAll_incr;
			RSDPatternPool->incomingSiteTotalAlleleCount += totAll_incr;
			
		}		
	}

	if(RSDPatternPool->incomingSiteDerivedAlleleCount!=0 && RSDPatternPool->incomingSiteDerivedAlleleCount!=RSDDataset->numberOfSamples && maf_check(RSDPatternPool->incomingSiteDerivedAlleleCount, RSDPatternPool->incomingSiteTotalAlleleCount, maf, &RSDDataset->setSitesDiscardedMafCheckFailed, 0)==1)
		RSDDataset->setSNPs++;
	else
		RSDPatternPool->incomingSitePosition = -1.0; // This is only useful if the last site of the set is not a SNP.

	while((RSDPatternPool->incomingSiteDerivedAlleleCount==0||RSDPatternPool->incomingSiteDerivedAlleleCount==RSDDataset->numberOfSamples || maf_check(RSDPatternPool->incomingSiteDerivedAlleleCount, RSDPatternPool->incomingSiteTotalAlleleCount, maf, &RSDDataset->setSitesDiscardedMafCheckFailed, 0)!=1) && setDone==0) // keep loading until first SNP is found
		setDone = RSDDataset_getNextSNP_ms (RSDDataset, RSDPatternPool, RSDChunk, RSDCommandLine, length, maf, fpOut);

	return setDone;
}

/**************/
/* vcf format */
/**************/

char RSDDataset_goToNextSet_vcf (RSDDataset_t * RSDDataset)
{
	assert(RSDDataset->inputFilePtr!=NULL);

	fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
	char tchar=(char)fgetc(RSDDataset->inputFilePtr);

	char tstring[STRING_SIZE];
	while(tchar!=EOF)
	{
		if(tchar=='\n')
		{
			int rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
			assert(rcnt==1);

			if(strcmp(tstring, RSDDataset->setID))
			{
				strcpy(RSDDataset->setID, tstring);
				return tchar;
			}
		}

		fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
		tchar=(char)fgetc(RSDDataset->inputFilePtr);
	}

	return EOF;
}

int RSDDataset_getValidSampleList_vcf (RSDDataset_t * RSDDataset)
{
	char tstring[STRING_SIZE];

	int rcnt = fscanf(RSDDataset->sampleFilePtr, "%s", tstring);

	RSDDataset->sampleValidListSize=0;
	while(rcnt==1)
	{
		RSDDataset->sampleValidListSize++;

		RSDDataset->sampleValidList = realloc(RSDDataset->sampleValidList, sizeof(char**)*((unsigned long)RSDDataset->sampleValidListSize));
		assert(RSDDataset->sampleValidList);

		RSDDataset->sampleValidList[RSDDataset->sampleValidListSize-1] = (char*)malloc(sizeof(char)*STRING_SIZE);
		assert(RSDDataset->sampleValidList[RSDDataset->sampleValidListSize-1]!=NULL);

		strcpy(RSDDataset->sampleValidList[RSDDataset->sampleValidListSize-1], tstring);		

		rcnt = fscanf(RSDDataset->sampleFilePtr, "%s", tstring);
	}

	return 0;
}

void RSDDataset_validateSample(RSDDataset_t * RSDDataset, char * sampleName, int sampleIndex)
{
	RSDDataset->sampleValid[sampleIndex] = SAMPLE_IS_NOT_VALID; // assume that is not valid

	int i;
	for(i=0;i<RSDDataset->sampleValidListSize;i++)
	{
		if(!strcmp(RSDDataset->sampleValidList[i], sampleName))
		{
			RSDDataset->sampleValid[sampleIndex] = SAMPLE_IS_VALID;
			strcpy(RSDDataset->sampleValidList[i], "\0");
			return;
		}
	}
}

int RSDDataset_getValidSampleSize (RSDDataset_t * RSDDataset)
{
	int sampleSize = 0;

	int i;
	for(i=0;i<RSDDataset->numberOfSamplesVCF;i++) // CHECK: we need all samples in the vcf here
		if(RSDDataset->sampleValid[i]==SAMPLE_IS_VALID)
			sampleSize++;

	return sampleSize;
}

int RSDDataset_getNumberOfSamples_vcf (RSDDataset_t * RSDDataset)
{
	char tstring[STRING_SIZE];

	int rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	while(rcnt==1 && strcmp(tstring, "#CHROM"))
		rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);

	if(rcnt!=1)
	{
		RSDDataset->numberOfSamples = 0;
		assert(RSDDataset->numberOfSamples>=1);
		return RSDDataset->numberOfSamples;
	}

	if(rcnt==1 && !strcmp(tstring, "#CHROM"))
	{
		int i;
		for(i=0;i<8 && rcnt==1;i++)
			rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
		
		assert(rcnt==1); // check for incomplete vcf header line 

		int sampleCntr = 0;
		int sampleIndex = -1;

		tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
		char tchar=tstring[0]; 
		char sampleName [STRING_SIZE];
	
		while(tstring[0]!='\n' && tstring[0]!=EOF && tstring[0]!='\r')
		{
			rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
			assert(rcnt==1);

			if(RSDDataset->sampleFilePtr!=NULL)
			{
				if(sampleCntr==0)
				{
					strcpy(sampleName, tstring);
				}
				else
				{
					sampleName[0] = tchar;
					sampleName[1] = '\0';
					strcat(sampleName, tstring);
				}
				
				fprintf(RSDDataset->sampleFilePtr, "%s\n", sampleName);
			}
		
			sampleCntr++; // this counts the number of samples in the file, dont care about validity yet

			sampleIndex = sampleCntr-1;
			RSDDataset->sampleValid = realloc(RSDDataset->sampleValid, sizeof(int)*((unsigned long)sampleCntr));
			assert(RSDDataset->sampleValid!=NULL);
			RSDDataset->sampleValid[sampleIndex]=SAMPLE_IS_VALID; // default operation: all samples are valid

			if(RSDDataset->sampleValidListSize!=ALL_SAMPLES_VALID)
				RSDDataset_validateSample(RSDDataset, sampleName, sampleIndex);

			tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
			tchar = tstring[0];

			while(tstring[0]==' ' || tstring[0]=='\t')
			{
				tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
				tchar = tstring[0];
			}
		}
		assert(tstring[0]=='\n' || tstring[0]=='\r');
		ungetc(tstring[0], RSDDataset->inputFilePtr);		

		RSDDataset->numberOfSamplesVCF = sampleCntr;
		RSDDataset->numberOfSamples = RSDDataset_getValidSampleSize (RSDDataset);

		return RSDDataset->numberOfSamples;
	}

	assert(0);
	return -1;	
}

void RSDDataset_getSetRegionLength_vcf (RSDDataset_t * RSDDataset)
{
	double setRegionSize=0.0, setSizeTmp=0.0;
	int rcnt=0, setDone = 0;
	char tstring[STRING_SIZE], chromosome[STRING_SIZE];

	fpos_t 	setPosition;
	fgetpos(RSDDataset->inputFilePtr, &setPosition); // start position

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", chromosome);
	assert(rcnt==1);

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
	assert(rcnt==1);

	setRegionSize = (double)atof(tstring);

	tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
	while(tstring[0]!='\n' && tstring[0]!=EOF)
		tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // get to the end of the line

	while(setDone==0)
	{
		rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
		if(rcnt!=1)
			setDone = 1;		

		if(setDone==0 && !strcmp(tstring, chromosome))
		{
			assert(rcnt==1); 

			rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring);
			assert(rcnt==1);

			setSizeTmp = (double)atof(tstring);

			assert(setSizeTmp>=setRegionSize);
			setRegionSize = setSizeTmp;

			tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
			while(tstring[0]!='\n' && tstring[0]!=EOF)
				tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // get to the end of the line

			if(tstring[0]==EOF)
				setDone = 1;
		}
		else
		{
			setDone = 1;
		}	
	}

	fsetpos(RSDDataset->inputFilePtr, &setPosition); // restore start position
	RSDDataset->setRegionLength = (uint64_t)setRegionSize;
}

int RSDDataset_getFirstSNP_vcf (RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut)
{
	if(length!=0ull)
		RSDDataset->setRegionLength = length;
	else
		RSDDataset_getSetRegionLength_vcf (RSDDataset);

	int setDone = 0;
	while(!setDone && RSDPatternPool->incomingSitePosition<=-1.0) 
		setDone = RSDDataset_getNextSNP(RSDDataset, RSDPatternPool, RSDChunk, RSDCommandLine, RSDDataset->setRegionLength, maf, fpOut);

	return setDone;
}

int RSDDataset_getNextSNP_vcf (RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut)
{
	assert(RSDChunk!=NULL);

	char tstring[STRING_SIZE];
	int setDone = 0;

	fgetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));

	RSDPatternPool->incomingSitePosition = -1.0;

	int rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // chrom
	if(rcnt!=1) // hit eof
	{
		tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
		if(tstring[0]==EOF)
		{
			setDone = 1;
			return setDone;
		}
		else
			assert(rcnt==1);
	}
	assert(rcnt==1);

	if(strcmp(RSDDataset->setID, tstring))
	{	
		fsetpos(RSDDataset->inputFilePtr, &(RSDDataset->setPosition));
		setDone = 1;

		fseek(RSDDataset->inputFilePtr,-1,SEEK_CUR);

		return setDone; 
	}
	assert(!strcmp(RSDDataset->setID, tstring)); // make sure we are still in the same chrom	

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // position
	assert(rcnt==1);

	double vali = (double)atof(tstring);
	RSDPatternPool->incomingSitePosition = vali;

	if(vali>(double)length)
	{
		printf("%f vs %f\n", vali, (double)length);
		fflush(stdout);

		assert(vali<=(double)length);
	}

	RSDDataset->setSize++;
	RSDDataset->setProgress++; // site loaded
	RSDDataset->setSNPs += 0; // Init to 0 since we dont know yet whether this is a polymorphic one or not

	//fgetpos(RSDDataset->inputFilePtr, &(RSDChunk->posPosition)); // not used

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // id
	assert(rcnt==1);

	int successSum = 0;
	int alleleMapSize = 0;
	char alleleMap[STRING_SIZE];

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // ref
	assert(rcnt==1);
	
	if(strlen(tstring)==1)
	{
		successSum++;
		alleleMap[alleleMapSize++] = tstring[0];
		alleleMap[alleleMapSize] = '\0';		
	}

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // alt
	assert(rcnt==1);

	if(strlen(tstring)==1)
	{
		successSum++;
		alleleMap[alleleMapSize++] = tstring[0];
		alleleMap[alleleMapSize] = '\0';
	}

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // qual
	assert(rcnt==1);
	
	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // filter
	assert(rcnt==1);

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // info
	assert(rcnt==1);

	rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // format
	assert(rcnt==1);

	int gtLoc = getGXLocation_vcf (tstring, "GT");
	
	if(gtLoc!=-1)
		successSum++;

	int gpLoc = -1, glLoc = -1;

	if(gtLoc==-1)
	{
		gpLoc = getGXLocation_vcf (tstring, "GP");
		glLoc = getGXLocation_vcf (tstring, "GL");

		if(gpLoc!=-1 && glLoc!=-1)
			successSum++;
	}
	
	int firstSNP = 0;
	int skipSNP = 0;

	if(successSum>=3)
	{
		successSum = 1;
	
		int i, sampleSize = 0;
		char data[STRING_SIZE], data2[STRING_SIZE];
		
		if(RSDDataset->setSamples==-1)
			firstSNP = 1;

		if(firstSNP)
			RSDDataset->setSamples = RSDDataset->numberOfSamples;

		RSDPatternPool->incomingSiteDerivedAlleleCount = 0;
		RSDPatternPool->incomingSiteTotalAlleleCount = 0;

		for(i=0;i<RSDDataset->numberOfSamplesVCF;i++) 
		{
			rcnt = fscanf(RSDDataset->inputFilePtr, "%s", tstring); // format
			assert(rcnt==1);

			if(RSDDataset->sampleValid[i]==SAMPLE_IS_VALID)
			{
				getGTData_vcf(tstring, gtLoc, gpLoc, glLoc, data);
				skipSNP = getGTAlleles_vcf (data, alleleMap, alleleMapSize, data2, &RSDPatternPool->incomingSiteDerivedAlleleCount, &RSDPatternPool->incomingSiteTotalAlleleCount, (int)RSDCommandLine->ploidy);

				if(RSDCommandLine->imputePerSNP==1 || RSDCommandLine->createPatternPoolMask==1 || (!skipSNP))
				{
					sampleSize += strlen(data2);

					//resize
					if((sampleSize>RSDDataset->setSamples))
					{
						if(firstSNP)
						{
							RSDDataset->setSamples = sampleSize;
							RSDPatternPool->incomingSite = realloc(RSDPatternPool->incomingSite, sizeof(char)*((unsigned long)(RSDDataset->setSamples+1)));
							assert(RSDPatternPool->incomingSite!=NULL);
						}
						else
						{
							// ERROR with snp larger than the first
							fprintf(fpOut, "\n\nERROR: Wrong SNP size (L) found!\n\n\n");
							fprintf(stderr, "\n\nERROR: Wrong SNP size (L) found!\n\n\n");
							exit(0);
						}
					}
			
					memcpy(&(RSDPatternPool->incomingSite[((unsigned long)sampleSize)-strlen(data2)]), data2, strlen(data2));
				}
				else
				{
					//skipline
					tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
					while(tstring[0]!='\n')
						tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // get to the end of the line

					RSDDataset->setSitesDiscardedWithMissing++;

					assert(skipSNP==1);
								
					break;
				}
			}			
		}
		RSDPatternPool->incomingSite[sampleSize] = '\0';

		if(RSDCommandLine->imputePerSNP==1 || RSDCommandLine->createPatternPoolMask==1)
		{
			skipSNP = strictPolymorphic_check(RSDPatternPool->incomingSiteDerivedAlleleCount, RSDPatternPool->incomingSiteTotalAlleleCount, &RSDDataset->setSitesDiscardedStrictPolymorphicCheckFailed, 0)==1?0:1;

			if(RSDCommandLine->imputePerSNP==1 && (!skipSNP))
			{
				RSDDataset->setSitesImputedTotal+=RSDPatternPool_imputeIncomingSite (RSDPatternPool, RSDDataset->setSamples);
				assert(RSDPatternPool->incomingSiteTotalAlleleCount==RSDDataset->setSamples);
				skipSNP=0; // assume snp after impute
			}	
		}
	}
	else
	{
		assert(skipSNP==0);
		skipSNP=1; // no snp if header-check failed

		successSum = 0;
		RSDPatternPool->incomingSiteDerivedAlleleCount = 0;
		RSDPatternPool->incomingSiteTotalAlleleCount = 0;
		RSDPatternPool->incomingSite[0] = '\0';

		//skipline
		tstring[0] = (char)fgetc(RSDDataset->inputFilePtr);
		while(tstring[0]!='\n')
			tstring[0] = (char)fgetc(RSDDataset->inputFilePtr); // get to the end of the line

		RSDDataset->setSitesDiscardedHeaderCheckFailed++;		
	}

	// snp checks
	skipSNP = monomorphic_check(RSDPatternPool->incomingSiteDerivedAlleleCount, (int)RSDDataset->setSamples, &RSDDataset->setSitesDiscardedMonomorphic, skipSNP)==1?0:1;
	skipSNP = maf_check(RSDPatternPool->incomingSiteDerivedAlleleCount, RSDPatternPool->incomingSiteTotalAlleleCount, maf, &RSDDataset->setSitesDiscardedMafCheckFailed, skipSNP)==1?0:1;
	
	if(skipSNP)
	{
		assert(skipSNP==1);

		RSDDataset->setSitesDiscarded++;

		RSDPatternPool->incomingSitePosition = -1.0; // This is only useful if the last site of the set is not a SNP

		if(firstSNP)
		{
			RSDDataset->setSamples=-1;
			free(RSDPatternPool->incomingSite);
			RSDPatternPool->incomingSite =  NULL;
			RSDPatternPool->incomingSite = (char*)malloc(sizeof(char)*((unsigned long)(RSDDataset->numberOfSamples+1)));
			assert(RSDPatternPool->incomingSite!=NULL);
		}
	}
	else
	{
		RSDDataset->setSNPs++;
	}
	
	return setDone;
}

/**********/
/* Parser */
/**********/

void RSDDataset_initParser (RSDDataset_t * RSDDataset, FILE * fpOut, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDDataset!=NULL);

	if(!strcmp(RSDDataset->inputFileFormat, "invalid"))
	{
		fprintf(fpOut, "\nERROR: unexpected error during parser initialization\n\n");
		fprintf(stderr, "\nERROR: unexpected error during parser initialization\n\n");
		exit(0);
	}

	if(!strcmp(RSDDataset->inputFileFormat, "ms"))
	{
		RSDDataset_goToNextSet = &RSDDataset_goToNextSet_ms;
		RSDDataset_getValidSampleList = &RSDDataset_getValidSampleList_ms;
		RSDDataset_getNumberOfSamples = &RSDDataset_getNumberOfSamples_ms;		
		RSDDataset_getFirstSNP = &RSDDataset_getFirstSNP_ms;
		RSDDataset_getNextSNP = &RSDDataset_getNextSNP_ms;

		RSDDataset->inputFileIsMBS = RSDCommandLine->mbs==1?1:RSDDataset->inputFileIsMBS;

		return;
	}

	if(!strcmp(RSDDataset->inputFileFormat, "vcf"))
	{
		RSDDataset_goToNextSet = &RSDDataset_goToNextSet_vcf;
		RSDDataset_getValidSampleList = &RSDDataset_getValidSampleList_vcf;
		RSDDataset_getNumberOfSamples = &RSDDataset_getNumberOfSamples_vcf;
		RSDDataset_getFirstSNP = &RSDDataset_getFirstSNP_vcf;
		RSDDataset_getNextSNP = &RSDDataset_getNextSNP_vcf;

		return;
	}

	assert(0);
	return;	
}

void RSDDataset_resetSiteCounters (RSDDataset_t * RSDDataset)
{
	RSDDataset->setSitesDiscarded = 0;
	RSDDataset->setSitesDiscardedHeaderCheckFailed = 0;
	RSDDataset->setSitesDiscardedWithMissing = 0;
	RSDDataset->setSitesDiscardedMafCheckFailed = 0;
	RSDDataset->setSitesDiscardedStrictPolymorphicCheckFailed = 0;
	RSDDataset->setSitesDiscardedMonomorphic = 0;
	RSDDataset->setSitesImputedTotal = 0;
}

void RSDDataset_printSiteReport (RSDDataset_t * RSDDataset, FILE * fp, int setIndex, int64_t imputePerSNP, int64_t createPatternPoolMask)
{
	assert(fp!=NULL);

	assert(RSDDataset->setSize==RSDDataset->setSNPs+RSDDataset->setSitesDiscarded);

	if(imputePerSNP==0 && createPatternPoolMask==0) // M=0
	{
		assert(RSDDataset->setSitesDiscarded==RSDDataset->setSitesDiscardedHeaderCheckFailed+RSDDataset->setSitesDiscardedMafCheckFailed+RSDDataset->setSitesDiscardedWithMissing+RSDDataset->setSitesDiscardedMonomorphic);

		fprintf(fp,"\n %d: %s | %d = %d + %d | %d = %d + %d + %d + %d", setIndex, RSDDataset->setID, (int)RSDDataset->setSize, (int)RSDDataset->setSNPs, (int)RSDDataset->setSitesDiscarded, (int)RSDDataset->setSitesDiscarded, (int)RSDDataset->setSitesDiscardedHeaderCheckFailed, (int)RSDDataset->setSitesDiscardedMafCheckFailed, (int)RSDDataset->setSitesDiscardedWithMissing, (int)RSDDataset->setSitesDiscardedMonomorphic);
	}
	else
	{
		assert(RSDDataset->setSitesDiscarded==RSDDataset->setSitesDiscardedHeaderCheckFailed+RSDDataset->setSitesDiscardedMafCheckFailed+RSDDataset->setSitesDiscardedStrictPolymorphicCheckFailed);

		if(imputePerSNP==1)
			fprintf(fp,"\n %d: %s | %d = %d + %d | %d = %d + %d + %d | %d ", setIndex, RSDDataset->setID, (int)RSDDataset->setSize, (int)RSDDataset->setSNPs, (int)RSDDataset->setSitesDiscarded, (int)RSDDataset->setSitesDiscarded, (int)RSDDataset->setSitesDiscardedHeaderCheckFailed, (int)RSDDataset->setSitesDiscardedMafCheckFailed, (int)RSDDataset->setSitesDiscardedStrictPolymorphicCheckFailed, (int)RSDDataset->setSitesImputedTotal);
		else
			fprintf(fp,"\n %d: %s | %d = %d + %d | %d = %d + %d + %d ", setIndex, RSDDataset->setID, (int)RSDDataset->setSize, (int)RSDDataset->setSNPs, (int)RSDDataset->setSitesDiscarded, (int)RSDDataset->setSitesDiscarded, (int)RSDDataset->setSitesDiscardedHeaderCheckFailed, (int)RSDDataset->setSitesDiscardedMafCheckFailed, (int)RSDDataset->setSitesDiscardedStrictPolymorphicCheckFailed);
	}
	fflush(fp);
}

