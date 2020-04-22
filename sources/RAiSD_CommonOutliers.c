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

void 	RSDCommonOutliers_appendSorted 		(RSDCommonOutliers_t * RSDCommonOutliers, double pos, double sco, int mode);
void 	RSDCommonOutliers_writeToolReport 	(RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine, int mode);
int 	RSDCommonOutliers_findOutliers 		(RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine);

RSDCommonOutliers_t * RSDCommonOutliers_new (void)
{
	RSDCommonOutliers_t * co = NULL;

	co = rsd_malloc(sizeof(RSDCommonOutliers_t));
	assert(co!=NULL);

	strncpy(co->reportFilenameSweeD, "\0", STRING_SIZE);
	co->positionIndexSweeD = -1;
	co->scoreIndexSweeD = -1;

	co->reportSizeSweeD = -1;
	co->positionSweeD = NULL;
	co->scoreSweeD = NULL;

	co->coPointSizeSweeD = -1;
	co->coPointPositionSweeD = NULL;
	co->coPointScoreSweeD = NULL;

	strncpy(co->reportFilenameRAiSD, "\0", STRING_SIZE);
	co->positionIndexRAiSD = -1;
	co->scoreIndexRAiSD = -1;

	co->reportSizeRAiSD = -1;
	co->positionRAiSD = NULL;
	co->scoreRAiSD = NULL;

	co->coPointSizeRAiSD = -1;
	co->coPointPositionRAiSD = NULL;
	co->coPointScoreRAiSD = NULL;

	strncpy(co->report1Filename, "\0", STRING_SIZE);
	strncpy(co->report2Filename, "\0", STRING_SIZE);
	strncpy(co->common1Filename, "\0", STRING_SIZE);
	strncpy(co->common2Filename, "\0", STRING_SIZE);

	return co;
}

void RSDCommonOutliers_free (RSDCommonOutliers_t * RSDCommonOutliers)
{
	assert(RSDCommonOutliers!=NULL);

	if(RSDCommonOutliers->positionSweeD!=NULL)
		free(RSDCommonOutliers->positionSweeD);

	if(RSDCommonOutliers->scoreSweeD!=NULL)
		free(RSDCommonOutliers->scoreSweeD);

	if(RSDCommonOutliers->coPointPositionSweeD!=NULL)
		free(RSDCommonOutliers->coPointPositionSweeD);

	if(RSDCommonOutliers->coPointScoreSweeD!=NULL)
		free(RSDCommonOutliers->coPointScoreSweeD);

	if(RSDCommonOutliers->positionRAiSD!=NULL)
		free(RSDCommonOutliers->positionRAiSD);

	if(RSDCommonOutliers->scoreRAiSD!=NULL)
		free(RSDCommonOutliers->scoreRAiSD);

	if(RSDCommonOutliers->coPointPositionRAiSD!=NULL)
		free(RSDCommonOutliers->coPointPositionRAiSD);

	if(RSDCommonOutliers->coPointScoreRAiSD!=NULL)
		free(RSDCommonOutliers->coPointScoreRAiSD);	

	free(RSDCommonOutliers);
}

void RSDCommonOutliers_init (RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommonOutliers!=NULL);
	assert(RSDCommandLine!=NULL);

	strcpy(RSDCommonOutliers->reportFilenameSweeD, RSDCommandLine->reportFilenameSweeD);
	RSDCommonOutliers->positionIndexSweeD = RSDCommandLine->positionIndexSweeD;
	RSDCommonOutliers->scoreIndexSweeD = RSDCommandLine->scoreIndexSweeD;

	strcpy(RSDCommonOutliers->reportFilenameRAiSD, RSDCommandLine->reportFilenameRAiSD);
	RSDCommonOutliers->positionIndexRAiSD = RSDCommandLine->positionIndexRAiSD;
	RSDCommonOutliers->scoreIndexRAiSD = RSDCommandLine->scoreIndexRAiSD;
}

void RSDCommonOutliers_appendSorted (RSDCommonOutliers_t * RSDCommonOutliers, double pos, double sco, int mode)
{
	assert(RSDCommonOutliers!=NULL);
	assert(mode==SWEED_CO || mode==RAiSD_CO);

	int64_t i, j;

	if(mode==SWEED_CO)
	{
		RSDCommonOutliers->reportSizeSweeD++;

		RSDCommonOutliers->positionSweeD = rsd_realloc(RSDCommonOutliers->positionSweeD, sizeof(double)*(unsigned long)RSDCommonOutliers->reportSizeSweeD);
		assert(RSDCommonOutliers->positionSweeD!=NULL);

		RSDCommonOutliers->scoreSweeD = rsd_realloc(RSDCommonOutliers->scoreSweeD, sizeof(double)*(unsigned long)RSDCommonOutliers->reportSizeSweeD);
		assert(RSDCommonOutliers->scoreSweeD!=NULL);

		if(RSDCommonOutliers->reportSizeSweeD==1)
		{
			RSDCommonOutliers->positionSweeD[RSDCommonOutliers->reportSizeSweeD-1] = pos;
			RSDCommonOutliers->scoreSweeD[RSDCommonOutliers->reportSizeSweeD-1] = sco;
		}
		else
		{
			for(i=0;i<RSDCommonOutliers->reportSizeSweeD-2;i++)
			{
				if(sco>RSDCommonOutliers->scoreSweeD[i])
					break;
			}

			for(j=RSDCommonOutliers->reportSizeSweeD-1;j>i;j--)
			{
				RSDCommonOutliers->positionSweeD[j] = RSDCommonOutliers->positionSweeD[j-1];
				RSDCommonOutliers->scoreSweeD[j] = RSDCommonOutliers->scoreSweeD[j-1];
			}

			assert(j==i);

			RSDCommonOutliers->positionSweeD[i] = pos;
			RSDCommonOutliers->scoreSweeD[i] = sco;
		}
	}
	else // RAiSD_CO
	{
		RSDCommonOutliers->reportSizeRAiSD++;

		RSDCommonOutliers->positionRAiSD = rsd_realloc(RSDCommonOutliers->positionRAiSD, sizeof(double)*(unsigned long)RSDCommonOutliers->reportSizeRAiSD);
		assert(RSDCommonOutliers->positionRAiSD!=NULL);

		RSDCommonOutliers->scoreRAiSD = rsd_realloc(RSDCommonOutliers->scoreRAiSD, sizeof(double)*(unsigned long)RSDCommonOutliers->reportSizeRAiSD);
		assert(RSDCommonOutliers->scoreRAiSD!=NULL);

		if(RSDCommonOutliers->reportSizeRAiSD==1)
		{
			RSDCommonOutliers->positionRAiSD[RSDCommonOutliers->reportSizeRAiSD-1] = pos;
			RSDCommonOutliers->scoreRAiSD[RSDCommonOutliers->reportSizeRAiSD-1] = sco;
		}
		else
		{
			for(i=0;i<RSDCommonOutliers->reportSizeRAiSD-2;i++)
			{
				if(sco>RSDCommonOutliers->scoreRAiSD[i])
					break;
			}

			for(j=RSDCommonOutliers->reportSizeRAiSD-1;j>i;j--)
			{
				RSDCommonOutliers->positionRAiSD[j] = RSDCommonOutliers->positionRAiSD[j-1];
				RSDCommonOutliers->scoreRAiSD[j] = RSDCommonOutliers->scoreRAiSD[j-1];
			}

			assert(j==i);

			RSDCommonOutliers->positionRAiSD[i] = pos;
			RSDCommonOutliers->scoreRAiSD[i] = sco;
		}
	}
}

void RSDCommonOutliers_writeToolReport (RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine, int mode)
{
	assert(RSDCommonOutliers!=NULL);
	assert(mode==SWEED_CO || mode==RAiSD_CO);

	char filename [STRING_SIZE];
	int64_t reportSize = -1;
	if(mode==SWEED_CO)
	{	
		strcpy(filename, "RAiSD_SweeDReportSorted.");
		strcat(filename, RSDCommandLine->runName);
		strcpy(RSDCommonOutliers->report1Filename, filename);
		reportSize = RSDCommonOutliers->reportSizeSweeD;
	}
	else
	{
		strcpy(filename, "RAiSD_RAiSDReportSorted.");
		strcat(filename, RSDCommandLine->runName);
		strcpy(RSDCommonOutliers->report2Filename, filename);
		reportSize = RSDCommonOutliers->reportSizeRAiSD;
	}
	FILE * fp = fopen(filename, "w");
	assert(fp!=NULL);
	
	int i;

	for(i=0;i<reportSize;i++)
	{
		if(mode==SWEED_CO)
			fprintf(fp, "%f\t%e\n", RSDCommonOutliers->positionSweeD[i], RSDCommonOutliers->scoreSweeD[i]);
		else
			fprintf(fp, "%f\t%e\n", RSDCommonOutliers->positionRAiSD[i], RSDCommonOutliers->scoreRAiSD[i]);

	}

	fclose(fp);		
}

int RSDCommonOutliers_findOutliers (RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommonOutliers!=NULL);
	assert(RSDCommandLine!=NULL);

	double topP = (double)atof(RSDCommandLine->commonOutliersThreshold);
	sprintf(RSDCommandLine->commonOutliersThreshold, "%.2f", topP);

	fprintf(stdout, " Top\t\t\t%.2f\n", topP);
	fprintf(RAiSD_Info_FP, " Top\t\t\t%.2f\n", topP);

	int thres1Index = (int)(RSDCommonOutliers->reportSizeSweeD*topP);
	int thres2Index = (int)(RSDCommonOutliers->reportSizeRAiSD*topP);

	double thres1 = RSDCommonOutliers->scoreSweeD[thres1Index];
	double thres2 = RSDCommonOutliers->scoreRAiSD[thres2Index];

	fprintf(stdout, " SweeD threshold\t%f\n", thres1);
	fprintf(RAiSD_Info_FP, " SweeD threshold\t%f\n", thres1);

	fprintf(stdout, " RAiSD threshold\t%f\n", thres2);
	fprintf(RAiSD_Info_FP, " RAiSD threshold\t%f\n", thres2);

	RSDCommonOutliers->coPointSizeSweeD = 0;
	int i;
	for(i=0;i<RSDCommonOutliers->reportSizeSweeD;i++)
	{
		if(RSDCommonOutliers->scoreSweeD[i]>thres1)
		{
			RSDCommonOutliers->coPointSizeSweeD++;

			RSDCommonOutliers->coPointPositionSweeD = rsd_realloc(RSDCommonOutliers->coPointPositionSweeD, sizeof(double)*(unsigned long)RSDCommonOutliers->coPointSizeSweeD);
			assert(RSDCommonOutliers->coPointPositionSweeD!=NULL);

			RSDCommonOutliers->coPointScoreSweeD = rsd_realloc(RSDCommonOutliers->coPointScoreSweeD, sizeof(double)*(unsigned long)RSDCommonOutliers->coPointSizeSweeD);
			assert(RSDCommonOutliers->coPointScoreSweeD!=NULL);

			RSDCommonOutliers->coPointPositionSweeD[RSDCommonOutliers->coPointSizeSweeD-1] = RSDCommonOutliers->positionSweeD[i];
			RSDCommonOutliers->coPointScoreSweeD[RSDCommonOutliers->coPointSizeSweeD-1] = RSDCommonOutliers->scoreSweeD[i];
		}
	}

	RSDCommonOutliers->coPointSizeRAiSD = 0;
	for(i=0;i<RSDCommonOutliers->reportSizeRAiSD;i++)
	{
		if(RSDCommonOutliers->scoreRAiSD[i]>thres2)
		{
			RSDCommonOutliers->coPointSizeRAiSD++;

			RSDCommonOutliers->coPointPositionRAiSD = rsd_realloc(RSDCommonOutliers->coPointPositionRAiSD, sizeof(double)*(unsigned long)RSDCommonOutliers->coPointSizeRAiSD);
			assert(RSDCommonOutliers->coPointPositionRAiSD!=NULL);

			RSDCommonOutliers->coPointScoreRAiSD = rsd_realloc(RSDCommonOutliers->coPointScoreRAiSD, sizeof(double)*(unsigned long)RSDCommonOutliers->coPointSizeRAiSD);
			assert(RSDCommonOutliers->coPointScoreRAiSD!=NULL);

			RSDCommonOutliers->coPointPositionRAiSD[RSDCommonOutliers->coPointSizeRAiSD-1] = RSDCommonOutliers->positionRAiSD[i];
			RSDCommonOutliers->coPointScoreRAiSD[RSDCommonOutliers->coPointSizeRAiSD-1] = RSDCommonOutliers->scoreRAiSD[i];

		}
	}

	double prox = RSDCommandLine->commonOutliersMaxDistance;

	fprintf(stdout, " Max distance\t\t%.1f\n", prox);
	fprintf(RAiSD_Info_FP, " Max distance\t\t%.1f\n", prox);

	int j;

	char filename1 [STRING_SIZE];

	strcpy(filename1, "RAiSD_CommonOutlierPointsSweeD.");
	strcat(filename1, RSDCommandLine->runName);
	strcpy(RSDCommonOutliers->common1Filename, filename1);

	FILE * fp1 =  fopen(filename1, "w");
	assert(fp1!=NULL);

	char filename2 [STRING_SIZE];

	strcpy(filename2, "RAiSD_CommonOutlierPointsRAiSD.");
	strcat(filename2, RSDCommandLine->runName);
	strcpy(RSDCommonOutliers->common2Filename, filename2);

	FILE * fp2 =  fopen(filename2, "w");
	assert(fp2!=NULL);

	char filename3 [STRING_SIZE];

	strcpy(filename3, "RAiSD_CommonOutlierReport.");
	strcat(filename3, RSDCommandLine->runName);

	FILE * fp3 = fopen(filename3, "w");
	assert(fp3!=NULL);

	int flag1v = 0;
	int * flag1 = (int*)rsd_malloc(sizeof(int)*(unsigned long)RSDCommonOutliers->coPointSizeSweeD);
	assert(flag1!=NULL);

	int flag2v = 0;
	int * flag2 = (int*)rsd_malloc(sizeof(int)*(unsigned long)RSDCommonOutliers->coPointSizeRAiSD);
	assert(flag2!=NULL);

	for(i=0;i<RSDCommonOutliers->coPointSizeSweeD;i++)
		flag1[i] = 0;

	for(i=0;i<RSDCommonOutliers->coPointSizeRAiSD;i++)
		flag2[i] = 0;

	for(i=0;i<RSDCommonOutliers->coPointSizeSweeD;i++)
	{
		for(j=0;j<RSDCommonOutliers->coPointSizeRAiSD;j++)
		{
			if(fabs(RSDCommonOutliers->coPointPositionSweeD[i]-RSDCommonOutliers->coPointPositionRAiSD[j])<=prox)
			{
				if(!flag1[i])
				{	
					flag1[i] = 1;
					fprintf(fp1, "%f\t%e\n", RSDCommonOutliers->coPointPositionSweeD[i], RSDCommonOutliers->coPointScoreSweeD[i]);
					flag1v = 1;
				}
		
				if(!flag2[j])
				{
					flag2[j] = 1;
					fprintf(fp2, "%f\t%e\n", RSDCommonOutliers->coPointPositionRAiSD[j], RSDCommonOutliers->coPointScoreRAiSD[j]);
					flag2v = 1;
				}

				fprintf(fp3, "%f\t%f\n", RSDCommonOutliers->coPointPositionSweeD[i], RSDCommonOutliers->coPointPositionRAiSD[j]);
			}
		}
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	free(flag1);
	free(flag2);
	
	if(flag1v == 0 || flag2v == 0)
		return 0;
	
	return 1;	
}

void RSDCommonOutliers_process (RSDCommonOutliers_t * RSDCommonOutliers, RSDCommandLine_t * RSDCommandLine)
{
	assert(RSDCommonOutliers!=NULL);
	assert(RSDCommandLine!=NULL);

	if(!strcmp(RSDCommonOutliers->reportFilenameRAiSD, "\0"))
		return;

	fprintf(stdout, "\n");
	fprintf(stdout, " COMMON-OUTLIER ANALYSIS\n");

	fprintf(RAiSD_Info_FP, "\n");
	fprintf(RAiSD_Info_FP, " COMMON-OUTLIER ANALYSIS\n");

	FILE * fpReport = fopen(RSDCommonOutliers->reportFilenameSweeD, "r");
	assert(fpReport!=NULL);

	char tstring[STRING_SIZE];
	char tchar;

	int rcnt = fscanf(fpReport, "%s", tstring); 
	int rowIndex = 0, colIndex = 0, pair=0;
	double posSweeD = 0.0, scoSweeD = 0.0;

	RSDCommonOutliers->reportSizeSweeD = 0;

	while(rcnt!=EOF)
	{
		colIndex++;

		if(colIndex==RSDCommonOutliers->positionIndexSweeD && pair!=2)
		{
			posSweeD = (double)atof(tstring);
			pair++;
		}

		if(colIndex==RSDCommonOutliers->scoreIndexSweeD && pair!=2)
		{
			scoSweeD = (double)atof(tstring);
			pair++;
		}

		if(pair==2)
		{
			RSDCommonOutliers_appendSorted (RSDCommonOutliers, posSweeD, scoSweeD, SWEED_CO);
			pair = 0;
		}

		tchar = (char)fgetc(fpReport);
		if(tchar=='\n')
		{
			rowIndex++;
			colIndex=0;
			pair = 0;
		}
		else
			ungetc(tchar, fpReport);		

		rcnt = fscanf(fpReport, "%s", tstring); 
	}

	fclose(fpReport);

	fpReport = fopen(RSDCommonOutliers->reportFilenameRAiSD, "r");
	assert(fpReport!=NULL);

	rcnt = fscanf(fpReport, "%s", tstring);

	rowIndex = 0;
	colIndex = 0;
	pair=0;

	double posRAiSD = 0.0, scoRAiSD = 0.0;

	RSDCommonOutliers->reportSizeRAiSD = 0;

	while(rcnt!=EOF)
	{
		colIndex++;

		if(colIndex==RSDCommonOutliers->positionIndexRAiSD && pair!=2)
		{
			posRAiSD = (double)atof(tstring);
			pair++;
		}

		if(colIndex==RSDCommonOutliers->scoreIndexRAiSD && pair!=2)
		{
			scoRAiSD = (double)atof(tstring);
			pair++;
		}

		if(pair==2)
		{
			RSDCommonOutliers_appendSorted (RSDCommonOutliers, posRAiSD, scoRAiSD, RAiSD_CO);
			pair = 0;
		}

		tchar = (char)fgetc(fpReport);
		if(tchar=='\n')
		{
			rowIndex++;
			colIndex=0;
			pair = 0;
		}
		else
			ungetc(tchar, fpReport);		

		rcnt = fscanf(fpReport, "%s", tstring); 
	}

	fclose(fpReport);

	RSDCommonOutliers_writeToolReport (RSDCommonOutliers, RSDCommandLine, SWEED_CO);
	RSDCommonOutliers_writeToolReport (RSDCommonOutliers, RSDCommandLine, RAiSD_CO);

	int suc = RSDCommonOutliers_findOutliers (RSDCommonOutliers, RSDCommandLine);

	if(suc)
	{
		RSDPlot_generateRscript(RSDCommandLine, RSDPLOT_COMMONOUTLIERS);
		RSDPlot_createPlot (RSDCommandLine, NULL, NULL, RSDCommonOutliers, RSDPLOT_COMMONOUTLIERS);		
		RSDPlot_removeRscript(RSDCommandLine, RSDPLOT_COMMONOUTLIERS);
		
		strcpy(tstring, "RAiSD_SweeDReportSorted.");
		strcat(tstring, RSDCommandLine->runName);
		strcpy(RSDCommonOutliers->report1Filename, tstring);
		
		int ret = remove(tstring);
		ret = ret;
		assert(ret==0);	

		strcpy(tstring, "RAiSD_RAiSDReportSorted.");
		strcat(tstring, RSDCommandLine->runName);
		strcpy(RSDCommonOutliers->report1Filename, tstring);
		
		ret = remove(tstring);
		assert(ret==0);	
	}
	else
	{
		fprintf(stdout, "\nWARNING: No common outliers were found");
		fflush(stdout);

		fprintf(RAiSD_Info_FP, "\nWARNING: No common outliers were found");
		fflush(RAiSD_Info_FP);

		char outputFilename[STRING_SIZE];
		strncpy(outputFilename, "RAiSD_CommonOutlierPlot.", STRING_SIZE);
		strcat(outputFilename, RSDCommandLine->runName);

		FILE * fp = fopen(outputFilename, "r");
		if(fp!=NULL)
		{
			fprintf(stdout, "\nWARNING: %s was not updated\n\n", outputFilename);
			fflush(stdout);

			fprintf(RAiSD_Info_FP, "\nWARNING: %s was not updated\n\n", outputFilename);
			fflush(RAiSD_Info_FP);

			fclose(fp);
		}
		else
		{
			fprintf(stdout, "\nWARNING: %s was not created\n\n", outputFilename);
			fflush(stdout);

			fprintf(RAiSD_Info_FP, "\nWARNING: %s was not created\n\n", outputFilename);
			fflush(RAiSD_Info_FP);
		}
	}	
}
