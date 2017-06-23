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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <linux/limits.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <immintrin.h>
#include <time.h>
#include <sys/time.h>

double minld;
double maxld;
double accumld;
int cntld;
double ** ldmap;

/*Testing*/
uint64_t selectionTarget;
double MuVar_Accum;
double MuSfs_Accum;
double MuLd_Accum;
double Mu_Accum;
uint64_t selectionTargetDThreshold;
double MuVar_Success;
double MuSfs_Success;
double MuLd_Success;
double Mu_Success;
double fpr_loc;
int scr_svec_sz;
float * scr_svec;
float tpr_thres;
double tpr_scr;
int setIndexValid;
/**/

double global_pos;

#define STRING_SIZE 10000 
#define PATTERNPOOL_SIZE 1 // MBs
#define CHUNK_MEMSIZE_AND_INCREMENT 1024 // sites
#define MULTI_STEP_PARSING 0
#define SINGLE_STEP_PARSING 1 // set this to 0 to deactivate completely
#define WINDOW_SIZE 10
#define MIN_SET_SNPS 10
#define MIN_NUMBER_OF_SAMPLES 3
#define ALL_SAMPLES_VALID -1 // when -1, all samples are assumed valid and will be processed
#define SAMPLE_IS_VALID 1
#define SAMPLE_IS_NOT_VALID 0

#define OTHER_FORMAT -1
#define MS_FORMAT 0
#define FASTA_FORMAT 1
#define MACS_FORMAT 2
#define VCF_FORMAT 3
#define SF_FORMAT 4

#define MIN_MU_VAL 0.00000001

double StartTime;
double FinishTime;
double MemoryFootprint;

FILE * RAiSD_Info_FP;

// RAiSD_Support.c
unsigned long long 	rdtsc			(void);
double 			gettime			(void);
int			snpv_cmp 		(uint64_t * A, uint64_t * B, int size);
int			isnpv_cmp 		(uint64_t * A, uint64_t * B, int size, int numberOfSamples);
int			getGTLocation_vcf 	(char * string);
void			getGTData_vcf 		(char * string, int location, char * data);
void			getGTAlleles_vcf	(char * string, char * stateVector, int statesTotal, char * sampleData, int * alleleCount);
int			rsd_popcnt_u64		(uint64_t input);
float 			DIST 			(float a, float b);
float * 		putInSortVector		(int * size, float * vector, float value);

#ifndef _INTRINSIC_POPCOUNT
char POPCNT_U16_LUT [0x1u << 16];
int 			popcount_u32_iterative	(unsigned int n);
void 			popcount_u64_init 	(void);
#endif

// RAiSD_CommandLine.c
typedef struct
{
	char 		runName[STRING_SIZE]; // Flag: N
	char 		inputFileName[STRING_SIZE]; // Flag: I
	uint64_t	regionLength; // Flag: L
	int		overwriteOutput; // Flag: f
	int		splitOutput; // Flag: s 
	int		setSeparator; // Flag: t
	int		printSampleList; // Flag: p
	char 		sampleFileName[STRING_SIZE]; // Flag: S

} RSDCommandLine_t;

RSDCommandLine_t * 	RSDCommandLine_new			(void);
void 			RSDCommandLine_free			(RSDCommandLine_t * RSDCommandLine);
void 			RSDCommandLine_init			(RSDCommandLine_t * RSDCommandLine);
void 			RSDCommandLine_load			(RSDCommandLine_t * RSDCommandLine, int argc, char ** argv);
void 			RSDCommandLine_print			(int argc, char ** argv, FILE * fpOut);

// RAiSD_Chunk.c
typedef struct
{
	int 		chunkID; // index
	
	fpos_t		posPosition;
	fpos_t * 	seqPosition;

	int 		chunkMemSize; // preallocated

	float * 	sitePosition;
	int * 		derivedAlleleCount;
	int * 		patternID;	

	int		chunkSize; // number of snps

	int		derAll1CntTotal;
	int		derAllNCntTotal; 

} RSDChunk_t;

RSDChunk_t *	RSDChunk_new 			(void);
void 		RSDChunk_free			(RSDChunk_t * ch, int numberOfSamples);
void 		RSDChunk_init			(RSDChunk_t * RSDChunk, int numberOfSamples);
void		RSDChunk_reset			(RSDChunk_t * RSDChunk);

// RAiSD_PatternPool.c
typedef struct
{
	int 		memorySize; // MBs
	
	int 		maxSize; // maximum number of patterns that can be stored

	int 		patternSize; // 64-bit words

	int 		dataSize; // number of patterns stored in the pool

	char * 		incomingSite;
	int		incomingSiteAlleleCount;
	uint64_t * 	incomingSiteCompact;
	double		incomingSitePosition;

	uint64_t * 	poolData; // pattern data
	int *		poolDataAlleleCount; // number of derived alleles per pattern
	int *		poolDataPatternCount; // number of specific pattern occurences

	uint64_t * 	exchangeBuffer; // used to exchange patterns between location in the pool

} RSDPatternPool_t;

RSDPatternPool_t * 	RSDPatternPool_new		(void);
void 			RSDPatternPool_free		(RSDPatternPool_t * pp, int numberOfSamples);
void 			RSDPatternPool_init 		(RSDPatternPool_t * RSDPatternPool, RSDCommandLine_t * RSDCommandLine, int numberOfSamples);
void 			RSDPatternPool_print		(RSDPatternPool_t * RSDPatternPool, FILE * fpOut);
void 			RSDPatternPool_reset 		(RSDPatternPool_t * RSDPatternPool, int numberOfSamples, int setSamples, RSDChunk_t * RSDChunk);
int			RSDPatternPool_pushSNP		(RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, int numberOfSamples);
void			RSDPatternPool_resize 		(RSDPatternPool_t * RSDPatternPool, int setSamples, FILE * fpOut);
void 			RSDPatternPool_exhangePatterns 	(RSDPatternPool_t * RSDPatternPool, int pID_a, int pID_b);

// RAiSD_Dataset.c
typedef struct
{
	FILE * 		inputFilePtr;
	char 		inputFileFormat[STRING_SIZE];

	fpos_t 		setPosition; // set position, first "/" in ms or first line with dif chrom number in VCF
	int 		numberOfSamples; // -1 when unitialized

	int 		setParsingMode;
	
	char		setID[STRING_SIZE]; // for vcf this is the chrom number
	int		setSamples; // number of samples in the set - differs from numberOfSamples in vcf because of ploidy
	int 		setSize; // number of segsites, provided by the set header in ms, -1 for vcf - this can also include non-polymorphic sites
	int		setSNPs; // number of non-discarded sites
	int		setProgress; // number of loaded segsites

	uint64_t 	setRegionLength;

	FILE *		sampleFilePtr;
	int		sampleValidListSize; // number of names in the sampleValidList
	char **		sampleValidList; // list of names of valid samples
	int *		sampleValid; // this must be of size numberOfSamples with sampleValidListSize number of 1s or less, indicates which of the dataset samples are valid
	
	int		numberOfSamplesVCF; // this is the full number of samples in the vcf file
} RSDDataset_t;

RSDDataset_t * 	RSDDataset_new				(void);
void 		RSDDataset_free				(RSDDataset_t * RSDDataset);
void 		RSDDataset_init				(RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
void 		RSDDataset_print 			(RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
void 		RSDDataset_setPosition 			(RSDDataset_t * RSDDataset, int * setIndex);

void 		RSDDataset_initParser			(RSDDataset_t * RSDDataset, FILE * fpOut);
char 		(*RSDDataset_goToNextSet) 		(RSDDataset_t * RSDDataset);
int 		(*RSDDataset_getNumberOfSamples) 	(RSDDataset_t * RSDDataset);
int 		(*RSDDataset_getValidSampleList) 	(RSDDataset_t * RSDDataset);
int 		(*RSDDataset_getFirstSNP) 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, uint64_t length);
int 		(*RSDDataset_getNextSNP) 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, uint64_t length);
int 		RSDDataset_getNextSNP_ms 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, uint64_t length);
int 		RSDDataset_getNextSNP_vcf 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, uint64_t length);
void		RSDDataset_getSetRegionLength_ms	(RSDDataset_t * RSDDataset, uint64_t length);
void		RSDDataset_getSetRegionLength_vcf	(RSDDataset_t * RSDDataset);

// RAiSD_MuStatistic.c
typedef struct
{
	char 	reportName [STRING_SIZE];
	FILE *	reportFP;

	int	windowSize; // number of SNPs in each window

	int *	pCntVec;    // temp array used to count pattern occurences, contains (sequentially and at a stride): pattern count left, pattern count right, pattern count exclusive left, pattern count exclusive right
	
	float 	muVarMax; 
	float 	muVarMaxLoc;

	float 	muSfsMax; 
	float 	muSfsMaxLoc;

	float 	muLdMax; 
	float 	muLdMaxLoc;

	float 	muMax; 
	float 	muMaxLoc;

} RSDMuStat_t;

RSDMuStat_t * 	RSDMuStat_new 			(void);
void 		RSDMuStat_free 			(RSDMuStat_t * mu);
void 		RSDMuStat_init 			(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine);
void 		RSDMuStat_setReportName 	(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
void 		RSDMuStat_setReportNamePerSet 	(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine, FILE * fpOut, RSDDataset_t * RSDDataset);
void 		RSDMuStat_scanChunk 		(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);



