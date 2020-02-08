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

void RSD_init (void);

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
double tpr_thres;
double tpr_scr;
int setIndexValid;
double StartTime;
double FinishTime;
double MemoryFootprint;
FILE * RAiSD_Info_FP;
FILE * RAiSD_SiteReport_FP;
FILE * RAiSD_ReportList_FP;
struct timespec requestStart;
struct timespec requestEnd;

#ifdef _PTIMES
struct timespec requestStartMu;
struct timespec requestEndMu;
double TotalMuTime;
#endif

void RSD_header (FILE * fpOut)
{
	if(fpOut==NULL)
		return;

	printRAiSD (fpOut);

	fprintf(fpOut, "\n");
	fprintf(fpOut, " RAiSD, Raised Accuracy in Sweep Detection\n");
	fprintf(fpOut, " This is version %d.%d (released in %s %d)\n", MAJOR_VERSION, MINOR_VERSION, RELEASE_MONTH, RELEASE_YEAR);
	fprintf(fpOut, " Copyright (C) 2017, and GNU GPL'd, by Nikolaos Alachiotis and Pavlos Pavlidis\n");
	fprintf(fpOut, " Contact n.alachiotis/pavlidisp at gmail.com\n");
	fprintf(fpOut, "\n");
}

void RSD_init (void)
{
	clock_gettime(CLOCK_REALTIME, &requestStart);
	//StartTime = gettime();

	FinishTime = 0.0;
	MemoryFootprint = 0.0;

	RAiSD_Info_FP = NULL;
	RAiSD_SiteReport_FP = NULL;
	RAiSD_ReportList_FP = NULL;

	/*Testing*/
	selectionTarget = 0ull;
	selectionTargetDThreshold = 0ull;
	MuVar_Accum = 0.0;
 	MuSfs_Accum = 0.0;
 	MuLd_Accum = 0.0;
 	Mu_Accum = 0.0;
 	MuVar_Success = 0.0;
	MuSfs_Success = 0.0;
	MuLd_Success = 0.0;
	Mu_Success = 0.0;
	fpr_loc = 0.0;
	scr_svec_sz = 0;
	scr_svec = NULL;
	tpr_thres = 0.0;
	tpr_scr = 0.0;
	setIndexValid = -1;
	/**/

#ifndef _INTRINSIC_POPCOUNT
	popcount_u64_init();
#endif	

	srand((unsigned int)time(NULL)); // if no seed given

#ifdef _PTIMES
	TotalOoCTime = 0.0;
	TotalMuTime = 0.0;
#endif

}

int main (int argc, char ** argv)
{
	RSD_init();

	RSDCommandLine_t * RSDCommandLine = RSDCommandLine_new();
	RSDCommandLine_init(RSDCommandLine);
	RSDCommandLine_load(RSDCommandLine, argc, argv);

	RSDPlot_generateRscript(RSDCommandLine, RSDPLOT_BASIC_MU);
	
	RSD_header(stdout);
	RSD_header(RAiSD_Info_FP);
	RSD_header(RAiSD_SiteReport_FP);

	RSDCommandLine_print(argc, argv, stdout);
	RSDCommandLine_print(argc, argv, RAiSD_Info_FP);
	RSDCommandLine_print(argc, argv, RAiSD_SiteReport_FP);

	RSD_printSiteReportLegend(RAiSD_SiteReport_FP, RSDCommandLine->imputePerSNP, RSDCommandLine->createPatternPoolMask);

	RSDDataset_t * RSDDataset = RSDDataset_new();
	RSDDataset_init(RSDDataset, RSDCommandLine, RAiSD_Info_FP); 
	RSDDataset_print(RSDDataset, RSDCommandLine, stdout);
	RSDDataset_print(RSDDataset, RSDCommandLine, RAiSD_Info_FP);

	RSDPlot_printRscriptVersion (RSDCommandLine, stdout);
	RSDPlot_printRscriptVersion (RSDCommandLine, RAiSD_Info_FP);

	RSDCommandLine_printWarnings (RSDCommandLine, argc, argv, (void*)RSDDataset, stdout);
	RSDCommandLine_printWarnings (RSDCommandLine, argc, argv, (void*)RSDDataset, RAiSD_Info_FP);

	RSDPatternPool_t * RSDPatternPool = RSDPatternPool_new();
	RSDPatternPool_init(RSDPatternPool, RSDCommandLine, RSDDataset->numberOfSamples);
	RSDPatternPool_print(RSDPatternPool, stdout);
	RSDPatternPool_print(RSDPatternPool, RAiSD_Info_FP);

	RSDChunk_t * RSDChunk = RSDChunk_new();
	RSDChunk_init(RSDChunk, RSDDataset->numberOfSamples, RSDCommandLine->createPatternPoolMask);

	RSDMuStat_t * RSDMuStat = RSDMuStat_new();
	RSDMuStat_setReportName (RSDMuStat, RSDCommandLine, RAiSD_Info_FP);
	RSDMuStat_loadExcludeTable (RSDMuStat, RSDCommandLine);		

	int setIndex = -1, setDone = 0, setsProcessedTotal=0;

	// Set processing
	while(RSDDataset_goToNextSet(RSDDataset)!=EOF) 
	{
		RSDDataset_setPosition (RSDDataset, &setIndex);

		if(setIndexValid!=-1 && setIndex!=setIndexValid)
		{
			char tchar = (char)fgetc(RSDDataset->inputFilePtr); // Note: this might generate a segfault with MakefileZLIB
			assert(tchar==tchar);
		}

		if(setIndexValid==-1 || setIndex == setIndexValid)
		{
			RSDMuStat_setReportNamePerSet (RSDMuStat, RSDCommandLine, RAiSD_Info_FP, RSDDataset);

			if(RSDCommandLine->setSeparator)
				fprintf(RSDMuStat->reportFP, "// %s\n", RSDDataset->setID);	
	
			RSDMuStat_init (RSDMuStat, RSDCommandLine);
			RSDMuStat_excludeRegion (RSDMuStat, RSDDataset);

			setDone = 0;

			RSDChunk->chunkID = -1;	
			RSDChunk_reset(RSDChunk, RSDCommandLine);

			RSDPatternPool_reset(RSDPatternPool, RSDDataset->numberOfSamples, RSDDataset->setSamples, RSDChunk, RSDCommandLine);	

			setDone = RSDDataset_getFirstSNP(RSDDataset, RSDPatternPool, RSDChunk, RSDCommandLine, RSDCommandLine->regionLength, RSDCommandLine->maf, RAiSD_Info_FP);
			if(setDone)
			{
				if(RSDCommandLine->displayProgress==1)
				fprintf(stdout, "\n %d: Set %s | sites %d | snps %d | region %lu - skipped", setIndex, RSDDataset->setID, (int)RSDDataset->setSize, (int)RSDDataset->setSNPs, RSDDataset->setRegionLength);
				fprintf(RAiSD_Info_FP, "\n %d: Set %s | sites %d | snps %d | region %lu - skipped", setIndex, RSDDataset->setID, (int)RSDDataset->setSize, (int)RSDDataset->setSNPs, RSDDataset->setRegionLength);

				if(RSDCommandLine->displayDiscardedReport==1)
					RSDDataset_printSiteReport (RSDDataset, RAiSD_SiteReport_FP, setIndex, RSDCommandLine->imputePerSNP, RSDCommandLine->createPatternPoolMask);

				RSDDataset_resetSiteCounters (RSDDataset);	

				continue;
			}
			RSDPatternPool_resize (RSDPatternPool, RSDDataset->setSamples, RAiSD_Info_FP);
#ifdef _HM
			RSDHashMap_free(RSDPatternPool->hashMap);
			RSDPatternPool->hashMap = RSDHashMap_new();
			RSDHashMap_init (RSDPatternPool->hashMap, RSDDataset->setSamples, RSDPatternPool->poolData, RSDPatternPool->maxSize, RSDPatternPool->patternSize);
#endif
#ifdef _LM
			RSDLutMap_free(RSDPatternPool->lutMap);
			RSDPatternPool->lutMap =  RSDLutMap_new();
			RSDLutMap_init (RSDPatternPool->lutMap, RSDDataset->setSamples);
#endif
			RSDPatternPool_pushSNP (RSDPatternPool, RSDChunk, RSDDataset->setSamples, RSDCommandLine);

			int sitesloaded = 0;
			int patternsloaded = 0;
			
			// Chunk processing
			while(!setDone) 
			{
				RSDChunk->chunkID++;

				int poolFull = 0;

				// SNP processing
				while(!poolFull && !setDone) 
				{
					setDone = RSDDataset_getNextSNP(RSDDataset, RSDPatternPool, RSDChunk, RSDCommandLine, RSDDataset->setRegionLength, RSDCommandLine->maf, RAiSD_Info_FP);
					poolFull = RSDPatternPool_pushSNP (RSDPatternPool, RSDChunk, RSDDataset->setSamples, RSDCommandLine); 
				}
	
				RSDPatternPool_assessMissing (RSDPatternPool, RSDDataset->setSamples);

				// version 2.4, for ms
				if(!strcmp(RSDDataset->inputFileFormat, "ms"))
					assert((uint64_t)RSDDataset->setSNPs==RSDDataset->preLoadedsetSNPs);

#ifdef _PTIMES
				clock_gettime(CLOCK_REALTIME, &requestStartMu);
#endif
				// Compute Mu statistic
				RSDMuStat_scanChunk (RSDMuStat, RSDChunk, RSDPatternPool, RSDDataset, RSDCommandLine);
#ifdef _PTIMES
				clock_gettime(CLOCK_REALTIME, &requestEndMu);
				double MuTime = (requestEndMu.tv_sec-requestStartMu.tv_sec)+ (requestEndMu.tv_nsec-requestStartMu.tv_nsec) / BILLION;
				TotalMuTime += MuTime;
#endif
				sitesloaded+=RSDChunk->chunkSize;
				patternsloaded += RSDPatternPool->dataSize;

				RSDChunk_reset(RSDChunk, RSDCommandLine);
				RSDPatternPool_reset(RSDPatternPool, RSDDataset->numberOfSamples, RSDDataset->setSamples, RSDChunk, RSDCommandLine);
			}

			// Dist and Succ
			if(selectionTarget!=0ull)
			{
				MuVar_Accum += DIST (RSDMuStat->muVarMaxLoc, (double)selectionTarget);
				if(selectionTargetDThreshold!=0ull)
				{
					if(DIST (RSDMuStat->muVarMaxLoc, (double)selectionTarget)<=selectionTargetDThreshold)
						MuVar_Success += 1.0;
				}

				MuSfs_Accum += DIST (RSDMuStat->muSfsMaxLoc, (double)selectionTarget);
				if(selectionTargetDThreshold!=0ull)
				{
					if(DIST (RSDMuStat->muSfsMaxLoc, (double)selectionTarget)<=selectionTargetDThreshold)
						MuSfs_Success += 1.0;
				}

				MuLd_Accum += DIST (RSDMuStat->muLdMaxLoc, (double)selectionTarget);
				if(selectionTargetDThreshold!=0ull)
				{
					if(DIST (RSDMuStat->muLdMaxLoc, (double)selectionTarget)<=selectionTargetDThreshold)
						MuLd_Success += 1.0;
				}

				Mu_Accum += DIST (RSDMuStat->muMaxLoc, (double)selectionTarget);
				if(selectionTargetDThreshold!=0ull)
				{
					if(DIST (RSDMuStat->muMaxLoc, (double)selectionTarget)<=selectionTargetDThreshold)
						Mu_Success += 1.0;
				}
			}

			if(RSDCommandLine->createPlot==1)
				RSDPlot_createPlot (RSDCommandLine, RSDDataset, RSDMuStat, RSDPLOT_BASIC_MU);

			setsProcessedTotal++;

			if(RSDCommandLine->displayProgress==1)
				fprintf(stdout, "\n %d: Set %s | sites %d | snps %d | region %lu - Var %.0f %.3e | SFS %.0f %.3e | LD %.0f %.3e | MuStat %.0f %.3e", setIndex, RSDDataset->setID, (int)RSDDataset->setSize, (int)RSDDataset->setSNPs, RSDDataset->setRegionLength, (double)RSDMuStat->muVarMaxLoc, (double)RSDMuStat->muVarMax, (double)RSDMuStat->muSfsMaxLoc, (double)RSDMuStat->muSfsMax, (double)RSDMuStat->muLdMaxLoc, (double)RSDMuStat->muLdMax, (double)RSDMuStat->muMaxLoc, (double)RSDMuStat->muMax);

			fflush(stdout);

			if(RSDCommandLine->displayDiscardedReport==1)
				RSDDataset_printSiteReport (RSDDataset, RAiSD_SiteReport_FP, setIndex, RSDCommandLine->imputePerSNP, RSDCommandLine->createPatternPoolMask);

			RSDDataset_resetSiteCounters (RSDDataset);			

			fprintf(RAiSD_Info_FP, "\n %d: Set %s | sites %d | snps %d | region %lu - Var %.0f %.3e | SFS %.0f %.3e | LD %.0f %.3e | MuStat %.0f %.3e", setIndex, RSDDataset->setID, (int)RSDDataset->setSize, (int)RSDDataset->setSNPs, RSDDataset->setRegionLength, (double)RSDMuStat->muVarMaxLoc, (double)RSDMuStat->muVarMax, (double)RSDMuStat->muSfsMaxLoc, (double)RSDMuStat->muSfsMax, (double)RSDMuStat->muLdMaxLoc, (double)RSDMuStat->muLdMax, (double)RSDMuStat->muMaxLoc, (double)RSDMuStat->muMax);

			fflush(RAiSD_Info_FP);

			// FPR/TPR
			if(fpr_loc>0.0)
				scr_svec = putInSortVector(&scr_svec_sz, scr_svec, RSDMuStat->muMax);

			if(tpr_thres>0.0)
				if(RSDMuStat->muMax>=(float)tpr_thres)
					tpr_scr += 1.0;

			if(setIndex == setIndexValid)
				break;
		}			
	}

	fprintf(stdout, "\n\n");
	fprintf(stdout, " Sets (total):     %d\n", setIndex+1);
	fprintf(stdout, " Sets (processed): %d\n", setsProcessedTotal);
	fprintf(stdout, " Sets (skipped):   %d\n", setIndex+1-setsProcessedTotal);

	fprintf(RAiSD_Info_FP, "\n\n");
	fprintf(RAiSD_Info_FP, " Sets (total):     %d\n", setIndex+1);
	fprintf(RAiSD_Info_FP, " Sets (processed): %d\n", setsProcessedTotal);
	fprintf(RAiSD_Info_FP, " Sets (skipped):   %d\n", setIndex+1-setsProcessedTotal);

	fflush(stdout);
	fflush(RAiSD_Info_FP);

	if(selectionTarget!=0ull)
	{
		fprintf(stdout, "\n");
		fprintf(stdout, " AVERAGE DISTANCE (Target %lu)\n mu-VAR\t%.3f\n mu-SFS\t%.3f\n mu-LD\t%.3f\n MuStat\t%.3f\n", selectionTarget, MuVar_Accum/setsProcessedTotal, MuSfs_Accum/setsProcessedTotal, MuLd_Accum/setsProcessedTotal, Mu_Accum/setsProcessedTotal);

		if(selectionTargetDThreshold!=0ull)
		{
			fprintf(stdout, "\n");
			fprintf(stdout, " SUCCESS RATE (Distance %lu)\n mu-VAR\t%.3f\n mu-SFS\t%.3f\n mu-LD\t%.3f\n MuStat\t%.3f\n", selectionTargetDThreshold, MuVar_Success/setsProcessedTotal, MuSfs_Success/setsProcessedTotal, MuLd_Success/setsProcessedTotal, Mu_Success/setsProcessedTotal);
		}

		fprintf(RAiSD_Info_FP, "\n");
		fprintf(RAiSD_Info_FP, " AVERAGE DISTANCE (Target %lu)\n mu-VAR\t%.3f\n mu-SFS\t%.3f\n mu-LD\t%.3f\n MuStat\t%.3f\n", selectionTarget, MuVar_Accum/setsProcessedTotal, MuSfs_Accum/setsProcessedTotal, MuLd_Accum/setsProcessedTotal, Mu_Accum/setsProcessedTotal);

		if(selectionTargetDThreshold!=0ull)
		{
			fprintf(RAiSD_Info_FP, "\n");
			fprintf(RAiSD_Info_FP, " SUCCESS RATE (Distance %lu)\n mu-VAR\t%.3f\n mu-SFS\t%.3f\n mu-LD\t%.3f\n MuStat\t%.3f\n", selectionTargetDThreshold, MuVar_Success/setsProcessedTotal, MuSfs_Success/setsProcessedTotal, MuLd_Success/setsProcessedTotal, Mu_Success/setsProcessedTotal);
		}
	}

	if(fpr_loc>0.0)
	{	
		fprintf(stdout, "\n");
		fprintf(stdout, " SORTED DATA (FPR %f)\n Size\t\t\t%d\n Highest Score\t\t%.9f\n Lowest Score\t\t%.9f\n FPR Threshold\t\t%.9f\n Threshold Location\t%d\n", fpr_loc, scr_svec_sz, (double)scr_svec[0],(double)scr_svec[scr_svec_sz-1], (double)scr_svec[(int)(scr_svec_sz*fpr_loc)], (int)(scr_svec_sz*fpr_loc));

		fprintf(RAiSD_Info_FP, "\n");
		fprintf(RAiSD_Info_FP, " SORTED DATA (FPR %f)\n Size\t\t\t%d\n Highest Score\t\t%.9f\n Lowest Score\t\t%.9f\n FPR Threshold\t\t%.9f\n Threshold Location\t%d\n", fpr_loc, scr_svec_sz, (double)scr_svec[0], (double)scr_svec[scr_svec_sz-1], (double)scr_svec[(int)(scr_svec_sz*fpr_loc)], (int)(scr_svec_sz*fpr_loc));
	}

	if(tpr_thres>0.0)
	{
		fprintf(stdout, "\n");
		fprintf(stdout, " SCORE COUNT (Threshold %f)\n TPR\t%f\n", tpr_thres, tpr_scr/setsProcessedTotal);

		fprintf(RAiSD_Info_FP, "\n");
		fprintf(RAiSD_Info_FP, " SCORE COUNT (Threshold %f)\n TPR\t%f\n", tpr_thres, tpr_scr/setsProcessedTotal);
	}

	fflush(stdout);
	fflush(RAiSD_Info_FP);

	RSDPlot_removeRscript(RSDCommandLine, RSDPLOT_BASIC_MU);

	if(RSDCommandLine->createMPlot==1)
	{
		RSDPlot_generateRscript(RSDCommandLine, RSDPLOT_MANHATTAN);
		RSDPlot_createPlot (RSDCommandLine, RSDDataset, RSDMuStat, RSDPLOT_MANHATTAN);		
		RSDPlot_removeRscript(RSDCommandLine, RSDPLOT_MANHATTAN);
	}

	RSDCommandLine_free(RSDCommandLine);
	RSDPatternPool_free(RSDPatternPool, RSDDataset->setSamples);
	RSDChunk_free(RSDChunk, RSDDataset->setSamples);
	RSDDataset_free(RSDDataset);
	RSDMuStat_free(RSDMuStat);

	if(scr_svec!=NULL)
		free(scr_svec);	
	
	RSD_printTime(stdout, RAiSD_Info_FP);
	RSD_printMemory(stdout, RAiSD_Info_FP);

	fclose(RAiSD_Info_FP);

	if(RAiSD_SiteReport_FP!=NULL)
		fclose(RAiSD_SiteReport_FP);

	return 0;
}
