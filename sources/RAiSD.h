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
#include <inttypes.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#ifdef _ZLIB
#include "zlib.h"
#endif

#define MAJOR_VERSION 2
#define MINOR_VERSION 8
#define RELEASE_MONTH "April"
#define RELEASE_YEAR 2020

/*Testing*/
extern uint64_t selectionTarget;
extern double MuVar_Accum;
extern double MuSfs_Accum;
extern double MuLd_Accum;
extern double Mu_Accum;
extern uint64_t selectionTargetDThreshold;
extern double MuVar_Success;
extern double MuSfs_Success;
extern double MuLd_Success;
extern double Mu_Success;
extern double fpr_loc;
extern int scr_svec_sz;
extern float * scr_svec;
extern double tpr_thres;
extern double tpr_scr;
extern int setIndexValid;
/**/

#ifdef _PTIMES
extern struct timespec requestStartOoC;
extern struct timespec requestEndOoC;
extern double TotalOoCTime;

extern struct timespec requestStartMu;
extern struct timespec requestEndMu;
extern double TotalMuTime;
#endif

#define WIN_MODE_FXD 0 // this is the default one
#define WIN_MODE_OSE 1 // one-sided expansion
#define WIN_MODE_FC 2 // floating-center
#define WIN_MODE_DSR 3 // double-sided reduction

#define STRING_SIZE 8192
#define PATTERNPOOL_SIZE 1 // MBs without the mask (actual memfootprint approx. double)
#define	PATTERNPOOL_SIZE_MASK_FACTOR 2
#define CHUNK_MEMSIZE_AND_INCREMENT 1024 // sites
#define MULTI_STEP_PARSING 0
#define SINGLE_STEP_PARSING 1 // set this to 0 to deactivate completely
#define DEFAULT_WINDOW_SIZE 50
#define MIN_WINDOW_SIZE 6
#define MIN_NUMBER_OF_SAMPLES 3
#define ALL_SAMPLES_VALID -1 // when -1, all samples are assumed valid and will be processed
#define SAMPLE_IS_VALID 1
#define SAMPLE_IS_NOT_VALID 0
#define MAX_COMMANDLINE_FLAGS 100

#define OTHER_FORMAT -1
#define MS_FORMAT 0
#define FASTA_FORMAT 1
#define MACS_FORMAT 2
#define VCF_FORMAT 3
#define SF_FORMAT 4

#define MIN_MU_VAL 0.00000001

#define BILLION 1E9

#define RSDPLOT_BASIC_MU 0
#define	RSDPLOT_MANHATTAN 1
#define RSDPLOT_COMMONOUTLIERS 2

#define LUTMAP_WIDTH		256 // 65536
#define	LUTMAP_GROUPSIZE	8 // 16
#define LUTMAP_INTERVAL 	8 // 4

#define	TREEMAP_REALLOC_INCR	1000
#define	TREEMAP_NODEPOOL_CHUNKSIZE 10000

#define VCF_FILE_CHECK_PASS 1
#define VCF_FILE_CHECK_FAIL 0

#define	L_FLAG_INDEX 2
#define B_FLAG_INDEX 17

#define M_FLAG_INDEX 6
#define Y_FLAG_INDEX 5

#define CO_FLAG_INDEX 22

#define FASTA2VCF_CONVERT_n_PROCESS 0
#define FASTA2VCF_CONVERT_n_EXIT 1

#define GAP '-'
#define AD 'A'
#define CY 'C'
#define GU 'G'
#define TH 'T'
#define UN 'N'
#define ad 'a'
#define cy 'c'
#define gu 'g'
#define th 't'
#define un 'n'

#define MAJORITY 0
#define PROBABILITY 1

#define SCOREBUFFER_REALLOC_INCR 1024

#define SWEED_CO 0
#define RAiSD_CO 1

// RAiSD.c
extern struct timespec requestStart;
extern struct timespec requestEnd;

extern double StartTime;
extern double FinishTime;
extern double MemoryFootprint;
extern FILE * RAiSD_Info_FP;
extern FILE * RAiSD_SiteReport_FP;
extern FILE * RAiSD_ReportList_FP;

void 			RSD_header 			(FILE * fpOut);

// RAiSD_Support.c
unsigned long long 	rdtsc				(void);
double 			gettime				(void);
int			snpv_cmp 			(uint64_t * A, uint64_t * B, int size);
int			snpv_cmp_cross_masks		(uint64_t * A, uint64_t * B, uint64_t * mA, uint64_t * mB, int size);
int			isnpv_cmp_cross_masks		(uint64_t * A, uint64_t * B, uint64_t * mA, uint64_t * mB, int size);
int			isnpv_cmp 			(uint64_t * A, uint64_t * B, int size, int numberOfSamples);
int			isnpv_cmp_with_mask		(uint64_t * A, uint64_t * B, uint64_t * mA, uint64_t * mB, int size, int numberOfSamples);
int			getGXLocation_vcf 		(char * string, char * GX);
int 			getGXData_vcf 			(char * string, int location, char * data);
void			getGTData_vcf 			(char * string, int locationGT, int locationGP, int locationGL, char * data);
int			getGTAlleles_vcf		(char * string, char * stateVector, int statesTotal, char * sampleData, int * derivedAlleleCount, int * totalAlleleCount, int ploidy);
int			rsd_popcnt_u64			(uint64_t input);
double 			DIST 				(double a, double b);
float * 		putInSortVector			(int * size, float * vector, float value);
char 			alleleMask_binary 		(char c, int * isDerived, int *isValid, FILE * fpOut);
char 			alleleMask_fasta 		(char c, int * isDerived, int *isValid, FILE * fpOut, char outgroupState);
int 			monomorphic_check 		(int incomingSiteDerivedAlleleCount, int setSamples, int64_t * cnt, int skipSNP);
int 			maf_check 			(int ac, int at, double maf, int64_t * cnt, int skipSNP);
int 			strictPolymorphic_check 	(int incomingSiteDerivedAlleleCount, int incomingSiteTotalAlleleCount, int64_t * cnt, int skipSNP);
void 			dataShuffleKnuth		(char * data, int startIndex, int endIndex);
extern void		ignoreLineSpaces		(FILE * fp, char * ent);
extern int 		flagMatch			(FILE * fp, char flag[], int flaglength, char tmp);
extern void 		RSD_printTime 			(FILE * fp1, FILE * fp2);
extern void 		RSD_printMemory 		(FILE * fp1, FILE * fp2);
extern void 		RSD_countMemory 		(size_t newMemSz);
int 			diploidyCheck			(char * data);
void 			getGPProbs 			(char * data, double *p00, double *p01, double * p11, int isLik);
void 			reconGT 			(char * data);
void 			RSD_printSiteReportLegend 	(FILE * fp, int64_t imputePerSNP, int64_t createPatternPoolMask);
extern void *		rsd_malloc			(size_t size);
extern void *		rsd_realloc			(void * p, size_t size);
void 			VCFFileCheck 			(void * vRSDDataset, FILE * fpX, char * fileName, FILE * fpOut);
int 			VCFFileCheckAndReorder		(void * vRSDDataset, FILE * fp, char * fileName, int overwriteOutput, FILE * fpOut);
void 			printRAiSD 			(FILE * fpOut);

#ifndef _INTRINSIC_POPCOUNT
extern char	 	POPCNT_U16_LUT [0x1u << 16];
int 			popcount_u32_iterative	(unsigned int n);
void 			popcount_u64_init 	(void);
#endif

#ifdef _ZLIB
int 			gzscanf 		(gzFile fp, char * string);
#endif

// RAiSD_LinkedList.c
typedef struct RSDLinkedListNode
{
	struct 	RSDLinkedListNode * prv;
	struct 	RSDLinkedListNode * nxt;
	double 	snpPosition;
	char * 	snpData;
} RSDLinkedListNode_t;

RSDLinkedListNode_t * 	RSDLinkedList_new 		(double pos, char * snp);
RSDLinkedListNode_t * 	RSDLinkedList_free		(RSDLinkedListNode_t * listHead);
RSDLinkedListNode_t * 	RSDLinkedList_addNode 		(RSDLinkedListNode_t * listHead, double pos, char * snp);
int 			RSDLinkedList_getSize 		(RSDLinkedListNode_t * listHead);
void 			RSDLinkedList_appendToFile 	(RSDLinkedListNode_t * listHead, FILE * fp);

// RAiSD_CommandLine.c
typedef struct
{	// -a is also reserved for the rand gen's seed
	char 		runName[STRING_SIZE]; // Flag: N
	char 		inputFileName[STRING_SIZE]; // Flag: I
	uint64_t	regionLength; // Flag: L for ms, Flag: B for vcf
	uint64_t	regionSNPs; // Flag: B
	int		overwriteOutput; // Flag: f
	int		splitOutput; // Flag: s 
	int		setSeparator; // Flag: t
	int		printSampleList; // Flag: p
	char 		sampleFileName[STRING_SIZE]; // Flag: S
	double		maf; // Flag: m
	int64_t		mbs; // Flag: b
	int64_t		imputePerSNP; //Flag: i
	int64_t		createPatternPoolMask; //Flag: M
	int64_t		patternPoolMaskMode; //Flag: M
	int64_t		displayProgress; // Flag: O
	int64_t		fullReport; // Flag: R
	int64_t		createPlot; // Flag: P
	char 		manhattanThreshold[STRING_SIZE]; // Flag: P
	int64_t		createMPlot; // Flag: A
	double		muThreshold; // Flag: H
	int64_t		ploidy; // Flag: y
	int64_t		displayDiscardedReport; // Flag: D
	int64_t		windowSize; // Flag: w
	int64_t		sfsSlack; // Flag: c
	char		excludeRegionsFile[STRING_SIZE]; // Flag: X
	int64_t		orderVCF; // Flag: o
	char		outgroupName[STRING_SIZE]; // Flag: C
	char		outgroupName2[STRING_SIZE]; // Flag: C2
	char		chromNameVCF[STRING_SIZE]; // Flag: H
	int64_t		fasta2vcfMode; // E
	int64_t		gridSize; // G
	int64_t		createCOPlot; // CO
	char 		reportFilenameSweeD[STRING_SIZE]; // CO
	int		positionIndexSweeD; // CO
	int		scoreIndexSweeD; // CO
	char 		reportFilenameRAiSD[STRING_SIZE]; // CO
	int		positionIndexRAiSD; // CO
	int		scoreIndexRAiSD; // CO
	char		commonOutliersThreshold[STRING_SIZE]; // COT
	double		commonOutliersMaxDistance; // COD

} RSDCommandLine_t;

void 			RSDHelp 				(FILE * fp);
void 			RSDVersions				(FILE * fp);
RSDCommandLine_t * 	RSDCommandLine_new			(void);
void 			RSDCommandLine_free			(RSDCommandLine_t * RSDCommandLine);
void 			RSDCommandLine_init			(RSDCommandLine_t * RSDCommandLine);
void 			RSDCommandLine_load			(RSDCommandLine_t * RSDCommandLine, int argc, char ** argv);
void 			RSDCommandLine_print			(int argc, char ** argv, FILE * fpOut);
void 			RSDCommandLine_printWarnings 		(RSDCommandLine_t * RSDCommandLine, int argc, char ** argv, void * RSDDataset, FILE * fpOut);

// RAiSD_Chunk.c
typedef struct
{
	int64_t		chunkID; // index
	int64_t		chunkMemSize; // preallocated
	int64_t		chunkSize; // number of snps

	fpos_t		posPosition;
	fpos_t * 	seqPosition;

	double * 	sitePosition;
	int * 		derivedAlleleCount;
	int * 		patternID;

	float *		chunkData; // type casting for alleleCount and patternID

	int		derAll1CntTotal;
	int		derAllNCntTotal; 

} RSDChunk_t;

RSDChunk_t *	RSDChunk_new 			(void);
void 		RSDChunk_free			(RSDChunk_t * ch, int64_t numberOfSamples);
void 		RSDChunk_init			(RSDChunk_t * RSDChunk, int64_t numberOfSamples, int64_t createPatternPoolMask);
void		RSDChunk_reset			(RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine);

// RAiSD_HashMap.c
typedef struct
{
	int64_t		addressListSize; // = sampleSize for convenience
	uint64_t **	addressList;

	uint64_t	addressListEntryMaxSize;
	uint64_t *	addressListEntrySize;
	int64_t		addressListEntryFull;

	int64_t		mainKey;
	int64_t		secondaryKey;

	uint64_t *	mainAddress;
	uint64_t *	secondaryAddress;

	uint64_t *	poolDataFractions;	

} RSDHashMap_t;

RSDHashMap_t * 	RSDHashMap_new				(void);
void		RSDHashMap_init 			(RSDHashMap_t * RSDHashMap, int64_t numberOfSamples, uint64_t * poolData, int maxSize, int patternSize);
void 		RSDHashMap_free 			(RSDHashMap_t * hm);
void		RSDHashMap_setMainKey 			(RSDHashMap_t * RSDHashMap, int64_t mainKey);
void		RSDHashMap_setSecondaryKey 		(RSDHashMap_t * RSDHashMap, int64_t secondaryKey);
int 		RSDHashMap_scanPatternPoolFractions 	(RSDHashMap_t * RSDHashMap, uint64_t * incomingSiteCompact, int patternSize, int numberOfSamples, int * match);

// RAiSD_LutMap.c
typedef struct
{
	int64_t		lutMapSize;
	uint8_t *	lutMapMem;
	uint8_t ** 	lutMap;	
	uint8_t *	lutMapMemC;
	uint8_t **	lutMapC;

} RSDLutMap_t;

RSDLutMap_t *	RSDLutMap_new		(void);
void		RSDLutMap_init		(RSDLutMap_t * lm, int64_t numberOfSamples);
void 		RSDLutMap_free 		(RSDLutMap_t * lm);
void 		RSDLutMap_reset		(RSDLutMap_t * lm);
int		RSDLutMap_scan		(RSDLutMap_t * lm, uint64_t * query);
int		RSDLutMap_scanC		(RSDLutMap_t * lm, uint64_t * query, int64_t patternSize, int64_t numberOfSamples);
void 		RSDLutMap_update	(RSDLutMap_t * lm, uint64_t * query);

// RAiSD_TreeMap.c
typedef struct RSDTreeNode_t
{
	struct RSDTreeNode_t * 	childNode[2];
#ifndef _TM_PATTERN_ARRAY
	int64_t			patternID; 		
#endif
} RSDTreeNode_t;

typedef struct 
{
	int64_t nodeMatrixSizeX;
	int64_t nodeMatrixSizeY;

	int64_t	nextNodeDimX;
	int64_t	nextNodeDimY;
	
	RSDTreeNode_t ** nodeMatrix;

} RSDTreeNodePool_t;

typedef struct
{
	int64_t			totalPatterns;
	int64_t			totalNodes;
	size_t			totalMemory;

	RSDTreeNode_t *		rootNode[2];

	// _TM_NODE_POOL
	RSDTreeNodePool_t *	nodePool[2];

	// _TM_PATTERN_ARRAY
	int64_t			maxPatterns; // Thats how many patterns the allocated mem can accommodate.
	RSDTreeNode_t ** 	patternID; // This ptr index in this vec is the patternID.
	

} RSDTreeMap_t;

RSDTreeMap_t * 	RSDTreeMap_new 			(void);
void 		RSDTreeMap_free			(RSDTreeMap_t * RSDTreeMap);
int64_t		RSDTreeMap_matchSNP 		(RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples);
int64_t		RSDTreeMap_matchSNPC 		(RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples);
int64_t		RSDTreeMap_updateTree 		(RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples);
int64_t		RSDTreeMap_updateTreeInit 	(RSDTreeMap_t * RSDTreeMap, void * RSDPatternPool_g, int64_t numberOfSamples, uint64_t * pattern);

// RAiSD_PatternPool.c
typedef struct
{
	int 		memorySize; // MBs
	
	int 		maxSize; // maximum number of patterns that can be stored

	int 		patternSize; // 64-bit words

	int 		dataSize; // number of patterns stored in the pool

	char * 		incomingSite;
	int		incomingSiteDerivedAlleleCount;
	int		incomingSiteTotalAlleleCount;
	uint64_t * 	incomingSiteCompact;
	uint64_t *	incomingSiteCompactMask;
	int64_t		incomingSiteCompactWithMissing;
	double		incomingSitePosition;

	uint64_t	createPatternPoolMask; // 0 if no mask used, 1 if mask used
	uint64_t	patternPoolMaskMode; // 0 for ignoring allele pairs with N, 1 for treating N as extra state
	uint64_t * 	poolData; // pattern data
	uint64_t *	poolDataMask; // pattern data mask (for handling N)
	int *		poolDataAlleleCount; // number of derived alleles per pattern
	int *		poolDataPatternCount; // number of specific pattern occurences
	int *		poolDataWithMissing; // yes/no if missing/not per pattern
	int *		poolDataMaskCount; // popcnt the mask only
	int *		poolDataAppliedMaskCount; // popcnt data&mask

	uint64_t * 	exchangeBuffer; // used to exchange patterns between location in the pool

	RSDHashMap_t *	hashMap;
	RSDLutMap_t * 	lutMap;
	RSDTreeMap_t * 	treeMap;

} RSDPatternPool_t;

RSDPatternPool_t * 	RSDPatternPool_new			(void);
void 			RSDPatternPool_free			(RSDPatternPool_t * pp, int64_t numberOfSamples);
void 			RSDPatternPool_init 			(RSDPatternPool_t * RSDPatternPool, RSDCommandLine_t * RSDCommandLine, int64_t numberOfSamples);
void 			RSDPatternPool_print			(RSDPatternPool_t * RSDPatternPool, FILE * fpOut);
void 			RSDPatternPool_reset 			(RSDPatternPool_t * RSDPatternPool, int64_t numberOfSamples, int64_t setSamples, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine);
int			RSDPatternPool_pushSNP			(RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, int64_t numberOfSamples, RSDCommandLine_t * RSDCommandLine);
void			RSDPatternPool_resize 			(RSDPatternPool_t * RSDPatternPool, int64_t setSamples, FILE * fpOut);
void 			RSDPatternPool_exchangePatterns 	(RSDPatternPool_t * RSDPatternPool, int pID_a, int pID_b);
void 			RSDPatternPool_exchangePatternsFractions(RSDPatternPool_t * RSDPatternPool, int pID_a, int pID_b);
int			RSDPatternPool_imputeIncomingSite 	(RSDPatternPool_t * RSDPatternPool, int64_t setSamples);
void			RSDPatternPool_assessMissing 		(RSDPatternPool_t * RSDPatternPool, int64_t numberOfSamples);

// RAiSD_Dataset.c
typedef struct
{
	FILE * 		inputFilePtr;
#ifdef _ZLIB
	gzFile 		inputFilePtrGZ;
#endif
	char 		inputFileFormat[STRING_SIZE];
	char		setID[STRING_SIZE]; // for vcf this is the chrom number
	int 		inputFileIsMBS;
	int 		numberOfSamples; // -1 when unitialized

	fpos_t 		setPosition; // set position, first "/" in ms or first line with dif chrom number in VCF
#ifdef _ZLIB
	z_off_t		setPositionGZ;
#endif
	int64_t		setParsingMode;
	
	int64_t		setSamples; // number of samples in the set - differs from numberOfSamples in vcf because of ploidy
	int64_t 	setSize; // number of segsites, provided by the set header in ms, -1 for vcf - this can also include non-polymorphic sites
	int64_t		setSNPs; // number of non-discarded sites
	int64_t		setProgress; // number of loaded segsites
	uint64_t	preLoadedsetSNPs;

	uint64_t 	setRegionLength;

	FILE *		sampleFilePtr;
	char **		sampleValidList; // list of names of valid samples
	int		sampleValidListSize; // number of names in the sampleValidList
	int		numberOfSamplesVCF; // this is the full number of samples in the vcf file
	int *		sampleValid; // this must be of size numberOfSamples with sampleValidListSize number of 1s or less, indicates which of the dataset samples are valid

	int64_t		setSitesDiscarded;
	int64_t		setSitesDiscardedHeaderCheckFailed;
	int64_t		setSitesDiscardedWithMissing;
	int64_t		setSitesDiscardedMafCheckFailed;
	int64_t		setSitesDiscardedStrictPolymorphicCheckFailed;
	int64_t		setSitesDiscardedMonomorphic;
	int64_t		setSitesImputedTotal;
	
	double		muVarDenom;

	char *		outgroupSequence;
	char		outgroupName[STRING_SIZE];

	char *		outgroupSequence2;
	char		outgroupName2[STRING_SIZE];

} RSDDataset_t;

RSDDataset_t * 	RSDDataset_new				(void);
void 		RSDDataset_free				(RSDDataset_t * RSDDataset);
void 		RSDDataset_init				(RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
void 		RSDDataset_print 			(RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
void 		RSDDataset_setPosition 			(RSDDataset_t * RSDDataset, int * setIndex);

void 		RSDDataset_initParser			(RSDDataset_t * RSDDataset, FILE * fpOut, RSDCommandLine_t * RSDCommandLine, int isGZ);
extern char	(*RSDDataset_goToNextSet) 		(RSDDataset_t * RSDDataset);
extern int	(*RSDDataset_getNumberOfSamples) 	(RSDDataset_t * RSDDataset);
extern int 	(*RSDDataset_getValidSampleList) 	(RSDDataset_t * RSDDataset);
extern int 	(*RSDDataset_getFirstSNP) 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
extern int	(*RSDDataset_getNextSNP) 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
int 		RSDDataset_getNextSNP_ms 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
int 		RSDDataset_getNextSNP_vcf 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
void		RSDDataset_getSetRegionLength_ms	(RSDDataset_t * RSDDataset, uint64_t length);
void		RSDDataset_getSetRegionLength_vcf	(RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
void 		RSDDataset_printSiteReport 		(RSDDataset_t * RSDDataset, FILE * fp, int setIndex, int64_t imputePerSNP, int64_t createPatternPoolMask);
void 		RSDDataset_resetSiteCounters 		(RSDDataset_t * RSDDataset);
void		RSDDataset_calcMuVarDenom		(RSDDataset_t * RSDDataset);
#ifdef _ZLIB
void		RSDDataset_getSetRegionLength_vcf_gz	(RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
int 		RSDDataset_getNextSNP_vcf_gz 		(RSDDataset_t * RSDDataset, RSDPatternPool_t * RSDPatternPool, RSDChunk_t * RSDChunk, RSDCommandLine_t * RSDCommandLine, uint64_t length, double maf, FILE * fpOut);
#endif

// RAiSD_CommonOutliers.c
typedef struct
{
	char 		reportFilenameSweeD[STRING_SIZE];
	int		positionIndexSweeD;
	int		scoreIndexSweeD;

	int64_t		reportSizeSweeD;
	double	*	positionSweeD;
	double	* 	scoreSweeD;

	int64_t		coPointSizeSweeD;
	double	*	coPointPositionSweeD;
	double	* 	coPointScoreSweeD;

	char 		reportFilenameRAiSD[STRING_SIZE];
	int		positionIndexRAiSD;
	int		scoreIndexRAiSD;

	int64_t		reportSizeRAiSD;
	double	*	positionRAiSD;
	double	* 	scoreRAiSD;

	int64_t		coPointSizeRAiSD;
	double	*	coPointPositionRAiSD;
	double	* 	coPointScoreRAiSD;

	char 		report1Filename[STRING_SIZE];
	char 		report2Filename[STRING_SIZE];
	char 		common1Filename[STRING_SIZE];
	char 		common2Filename[STRING_SIZE];

} RSDCommonOutliers_t;

RSDCommonOutliers_t * 	RSDCommonOutliers_new 		(void);
void 			RSDCommonOutliers_free 		(RSDCommonOutliers_t * RSDCommonOutliers);
void 			RSDCommonOutliers_init 		(RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine);
void 			RSDCommonOutliers_process 	(RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine);

// RAiSD_MuStatistic.c
typedef struct
{
	char 		reportName [STRING_SIZE];
	FILE *		reportFP;
	char		reportFPFileName[STRING_SIZE];

	int64_t		windowSize; // number of SNPs in each window

	int *		pCntVec;    // temp array used to count pattern occurences, contains (sequentially and at a stride): pattern count left, pattern count right, pattern count exclusive left, pattern count exclusive right
	
	float 		muVarMax; 
	float 		muSfsMax; 
	float 		muLdMax; 
	float 		muMax;

	double 		muVarMaxLoc;
	double 		muSfsMaxLoc;
	double		muLdMaxLoc; 
	double 		muMaxLoc;

#ifdef _MUMEM
	int64_t		muReportBufferSize;
	float * 	muReportBuffer;
#endif

	int64_t		exclTableSize;
	char **		exclTableChromName;
	int64_t *	exclTableRegionStart;
	int64_t *	exclTableRegionStop;

	int64_t		excludeRegionsTotal;
	int64_t	*	excludeRegionStart;
	int64_t	*	excludeRegionStop;

	int64_t		bufferMemMaxSize;
	double *	buffer0Data; // windowCenter
	double *	buffer1Data; // windowStart	
	double *	buffer2Data; // windowEnd
	double *	buffer3Data; // muVAR
	double *	buffer4Data; // muSFS
	double *	buffer5Data; // muLD
	double *	buffer6Data; // mu
	int64_t		currentScoreIndex;

} RSDMuStat_t;

RSDMuStat_t * 	RSDMuStat_new 			(void);
void 		RSDMuStat_free 			(RSDMuStat_t * mu);
void 		RSDMuStat_init 			(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine);
void 		RSDMuStat_setReportName 	(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
void 		RSDMuStat_setReportNamePerSet 	(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine, FILE * fpOut, RSDDataset_t * RSDDataset, RSDCommonOutliers_t * RSDCommonOutliers);
extern void	(*RSDMuStat_scanChunk) 		(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
void 		RSDMuStat_scanChunkBinary	(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
void 		RSDMuStat_scanChunkWithMask	(RSDMuStat_t * RSDMuStat, RSDChunk_t * RSDChunk, RSDPatternPool_t * RSDPatternPool, RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine);
extern float   	getPatternCount			(RSDPatternPool_t * RSDPatternPool, int * pCntVec, int offset, int * patternID, int p0, int p1, int p2, int p3, int * pcntl, int * pcntr, int * pcntexll, int * pcntexlr);
void		RSDMuStat_loadExcludeTable 	(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine);
void		RSDMuStat_excludeRegion 	(RSDMuStat_t * RSDMuStat, RSDDataset_t * RSDDataset);
extern void 	(*RSDMuStat_storeOutput) 	(RSDMuStat_t * RSDMuStat, double windowCenter, double windowStart, double windowEnd, double muVar, double muSfs, double muLd, double mu);
void 		RSDMuStat_writeBuffer2File 	(RSDMuStat_t * RSDMuStat, RSDCommandLine_t * RSDCommandLine);

// RAiSD_Plot.c
void 		RSDPlot_printRscriptVersion 	(RSDCommandLine_t * RSDCommandLine, FILE * fpOut);
int 		RSDPlot_checkRscript 		(void);
void 		RSDPlot_createRscriptName 	(RSDCommandLine_t * RSDCommandLine, char * scriptName);
void 		RSDPlot_generateRscript 	(RSDCommandLine_t * RSDCommandLine, int mode);
void 		RSDPlot_removeRscript 		(RSDCommandLine_t * RSDCommandLine,int mode);
void 		RSDPlot_createPlot 		(RSDCommandLine_t * RSDCommandLine, RSDDataset_t * RSDDataset, RSDMuStat_t * RSDMuStat, RSDCommonOutliers_t * RSDCommonOutliers, int mode);
void 		RSDPlot_createReportListName 	(RSDCommandLine_t * RSDCommandLine, char * reportListName);

// RAiSD_Fasta2Vcf.c
typedef struct
{
	int		statesMax;
	int		statesTotal;
	char *		statesChar;
	int *		statesCount;

	int 		imputeStatesMax;
	int 		imputeStatesTotal;
	char * 		imputeStatesChar;
	double * 	imputeStatesProb;

} RSDFastaStates_t;

void RSDDataset_convertFasta2VCF (RSDDataset_t * RSDDataset, RSDCommandLine_t * RSDCommandLine, FILE * fpOut);

