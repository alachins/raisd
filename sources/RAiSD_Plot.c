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

int RSDPlot_checkRscript (void)
{
	int success = 1;

#ifdef _C1
	if (system("which Rscript > /dev/null 2>&1")) 
		success = 1;
	else
		success = 0;
#else
	FILE *fp;

	fp = fopen("Rscript.txt", "r");
	assert(fp==NULL);
	
	fp = popen("Rscript > /dev/null 2>Rscript.txt", "r");
	assert(fp!=NULL);

	int ret = pclose(fp);
	assert(ret!=-1);

	fp = fopen("Rscript.txt", "r");
	assert(fp!=NULL);

	char tstring[STRING_SIZE];
	ret = fscanf(fp, "%s", tstring);
	assert(ret==1);
	
	if(!strcmp(tstring, "Usage:"))
		success = 0;

	fclose(fp);

	ret = remove("Rscript.txt");
	assert(ret==0);	
#endif

	return success;	
}

void RSDPlot_printRscriptVersion (RSDCommandLine_t * RSDCommandLine, FILE * fpOut)
{
	if(RSDCommandLine->createPlot==0 && RSDCommandLine->createMPlot==0)
		return;

	fprintf(fpOut, " Rscript: ");

	FILE * fp;
	fp = fopen("RscriptVersion.txt", "r");
	assert(fp==NULL);
	
#ifdef _C1
	int ret = system("Rscript --version > /dev/null 2>RscriptVersion.txt");
	assert(ret!=-1);
	ret = ret;
#else
	fp = popen("Rscript --version > /dev/null 2>RscriptVersion.txt", "r");
	assert(fp!=NULL);

	int ret = pclose(fp);
	assert(ret!=-1);
	ret = ret;
#endif

	fp = fopen("RscriptVersion.txt", "r");
	assert(fp!=NULL);
	
	char tchar;

	tchar = (char)fgetc(fp);	

	while(tchar!=EOF)
	{
		fprintf(fpOut, "%c", tchar);
		tchar = (char)fgetc(fp);
	}

	fclose(fp);

	ret = remove("RscriptVersion.txt");
	assert(ret==0);	
	
	fflush(fpOut);
}

void RSDPlot_createRscriptName (RSDCommandLine_t * RSDCommandLine, char * scriptName)
{
	assert(scriptName!=NULL);

	strcpy(scriptName, "RSDPlot_");
	strcat(scriptName, RSDCommandLine->runName);
	strcat(scriptName, ".R");
}

void RSDPlot_createReportListName (RSDCommandLine_t * RSDCommandLine, char * reportListName)
{
	assert(reportListName!=NULL);

	strcpy(reportListName, "RAiSD_ReportList.");
	strcat(reportListName, RSDCommandLine->runName);
	strcat(reportListName, ".txt");
}

void RSDPlot_generateRscript (RSDCommandLine_t * RSDCommandLine, int mode)
{
	if(RSDCommandLine->createPlot==0 && RSDCommandLine->createMPlot==0 && RSDCommandLine->createCOPlot==0)
		return;

	char scriptName[STRING_SIZE]="RSDPlot.R";
	RSDPlot_createRscriptName (RSDCommandLine, scriptName);

	FILE * fp;
	fp = fopen(scriptName, "r");

	if(fp!=NULL)
	{
		if(RSDCommandLine->overwriteOutput==0)
		{
			fprintf(stderr, "\nERROR: Rscript file %s exists. Use -f to overwrite it.\n\n", scriptName);
			exit(0);
		}
		else
		{
			fclose(fp);
		}
	}
	
	fp = fopen(scriptName, "w");
	assert(fp!=NULL);

	if(mode==RSDPLOT_BASIC_MU)
	{
		const char * tscript = " \n \
args = commandArgs(trailingOnly=TRUE) \n \
RUNNAME <- args[1] \n \
ID <- args[2] \n \
output_file <- paste(\"RAiSD_Plot.\", RUNNAME,\".\",ID,\".pdf\", sep=\"\") \n \
RDPATH <- paste(\"RAiSD_Report.\",RUNNAME,\".\",ID, sep = \"\") \n \
rsd_data <- read.table(RDPATH, header=F) \n \
mup <- rsd_data[,1]/1000.0 \n \
mu1 <- rsd_data[,4] \n \
mu2 <- rsd_data[,5] \n \
mu3 <- rsd_data[,6] \n \
mu  <- rsd_data[,7] \n \
pdf(output_file, width=10, height=10) \n \
par(mfrow=c(2,2)) \n \
plot(mup, mu1, col=\"darkgray\", pch=16, cex = .6, ylab=bquote(mu ~ \"_var\"), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"_var curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
plot(mup, mu2, col=\"darkgray\", pch=16, cex = .6, ylab=bquote(mu ~ \"_sfs\"), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"_sfs curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
plot(mup, mu3, col=\"darkgray\", pch=16, cex = .6, ylab=bquote(mu ~ \"_ld\"), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"_ld curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
plot(mup, mu, col=\"black\", pch=16, cex = .6, ylab=bquote(mu), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
dev.off() \n";

		fprintf(fp,"%s", tscript);
	}

	if(mode==RSDPLOT_MANHATTAN)
	{
		const char * tscript = " \n \
args = commandArgs(trailingOnly=TRUE) \n \
library(qqman) \n \
RUNNAME <- args[1] \n \
THRESHOLD <- args[2] \n \
reportListN <- paste(\"RAiSD_ReportList.\", RUNNAME,\".txt\", sep=\"\") \n \
reportList <- read.table(reportListN)[,1] \n \
data <- data.frame(pos=c(), value=c(), chr=c()) \n \
for (i in 1:length(reportList[])) { \n \
d <- read.table(paste(reportList[i]), header=F, skip=0)[,1:7] \n \
tmp.dat <- data.frame(pos=d[,1], value=d[,7]*1, chr=rep(i, length(d[,1]))) \n \
data <- rbind(data, tmp.dat)} \n \
output_file <- paste(\"RAiSD_ManhattanPlot.\", RUNNAME,\".pdf\", sep=\"\") \n \
pdf(output_file) \n \
snp <- 1:dim(data)[1] \n \
mydf <- data.frame(snp, data) \n \
thres<-as.numeric(THRESHOLD) \n \
topQ<-thres*100 \n \
threshold <- quantile(x=data$value, probs = thres) \n \
title_msg <- paste(\"Manhattan plot for \", RUNNAME,\"\nthreshold=\",threshold,\" (\", topQ,\"%)\") \n \
manhattan(mydf, chr=\"chr\", bp=\"pos\", cex = 0.5, p=\"value\", snp=\"snp\", logp=F,  ylab=bquote(mu ~ \"statistic\"), genomewideline = FALSE, suggestiveline = FALSE,pos=0, col=c(\"blue2\", \"darkorange1\"), ylim=c(0.0, max(data$value, na.rm = TRUE)*1.4)) \n \
title(title_msg) \n \
abline(h=threshold, col=\"red\", lw=2) \n \
dev.off() \n";

		fprintf(fp,"%s", tscript);
	}

	if(mode==RSDPLOT_COMMONOUTLIERS)
	{
		const char * tscript = "\n \
 args = commandArgs(trailingOnly=TRUE) \n \
 d1 <- read.table( paste(args[3]), header=F, skip=0)[,1:2] \n \
 d2 <- read.table( paste(args[4]), header=F, skip=0)[,1:2] \n \
 d3 <- read.table( paste(args[5]), header=F, skip=0)[,1:2] \n \
 d4 <- read.table( paste(args[6]), header=F, skip=0)[,1:2] \n \
 pdf(paste(args[7]) , width=7, height=7) \n \
 par(mfrow=c(2,1)) \n \
 plot(d1[,1], d1[,2], col=\"darkgray\", pch=16, ylab=\"\", xlab=\"\") \n \
 points(d3[,1], d3[,2], col=\"red\", pch=19, cex=1.0) \n \
 mtext(side=1, text=\"Position\", 2) \n \
 mtext(side=2, text=\"SweeD\", 2) \n \
 title(paste(\"SweeD-RAiSD common outliers for \", args[1], \" ( top \", args[2],\" )\", sep=\"\")) \n \
 plot(d2[,1], d2[,2], col=\"darkgray\", pch=16, ylab=\"\", xlab=\"\") \n \
 points(d4[,1], d4[,2], col=\"red\", pch=19, cex=1.0) \n \
 mtext(side=1, text=\"Position\", 2) \n \
 mtext(side=2, text=\"RAiSD\", 2) \n \
 dev.off() \n";

		fprintf(fp,"%s", tscript);
	}

	fclose(fp);	
}

void RSDPlot_removeRscript (RSDCommandLine_t * RSDCommandLine, int mode)
{
	char scriptName[STRING_SIZE]="RSDPlot.R";
	char reportListName[STRING_SIZE]="RAiSD_ReportList.txt";
	RSDPlot_createRscriptName (RSDCommandLine, scriptName);

	FILE * fp;
	fp = fopen(scriptName, "r"); 
	if(fp!=NULL)
	{
		fclose(fp);
		int ret = remove(scriptName);
		assert(ret==0);	
		ret = ret;
	}
	fp=NULL;

	if(mode==RSDPLOT_MANHATTAN)
	{
		if(RAiSD_ReportList_FP!=NULL)
		{
			fclose(RAiSD_ReportList_FP);
			RSDPlot_createReportListName (RSDCommandLine, reportListName);
			int ret = remove(reportListName);
			assert(ret==0);
			ret = ret;	
		}
		RAiSD_ReportList_FP=NULL;
	}
}

void RSDPlot_createPlot (RSDCommandLine_t * RSDCommandLine, RSDDataset_t * RSDDataset, RSDMuStat_t * RSDMuStat, RSDCommonOutliers_t * RSDCommonOutliers, int mode)
{
	char tstring[STRING_SIZE];

	char scriptName[STRING_SIZE]="RSDPlot.R";
	RSDPlot_createRscriptName (RSDCommandLine, scriptName);

	if(mode==RSDPLOT_BASIC_MU)
	{
		strcpy(tstring, "Rscript ");
		strcat(tstring, scriptName);
		strcat(tstring, " ");	
		strcat(tstring, RSDCommandLine->runName);
		strcat(tstring, " ");
		strcat(tstring, RSDDataset->setID);
		strcat(tstring, " > /dev/null 2>&1");
	}

	if(mode==RSDPLOT_MANHATTAN)
	{
		strcpy(tstring, "Rscript ");
		strcat(tstring, scriptName);
		strcat(tstring, " ");	
		strcat(tstring, RSDCommandLine->runName);
		strcat(tstring, " ");
		strcat(tstring, RSDCommandLine->manhattanThreshold);
		strcat(tstring, " > /dev/null 2>&1");
	}

	if(mode==RSDPLOT_COMMONOUTLIERS)
	{
		char outputFilename[STRING_SIZE];
		strncpy(outputFilename, "RAiSD_CommonOutlierPlot.", STRING_SIZE);
		strcat(outputFilename, RSDCommandLine->runName);
		strcat(outputFilename, ".pdf");

		strcpy(tstring, "Rscript ");
		strcat(tstring, scriptName);
		strcat(tstring, " ");	
		strcat(tstring, RSDCommandLine->runName);
		strcat(tstring, " ");
		strcat(tstring, RSDCommandLine->commonOutliersThreshold);
		strcat(tstring, " ");
		strcat(tstring, RSDCommonOutliers->report1Filename);
		strcat(tstring, " ");
		strcat(tstring, RSDCommonOutliers->report2Filename);
		strcat(tstring, " ");
		strcat(tstring, RSDCommonOutliers->common1Filename);
		strcat(tstring, " ");
		strcat(tstring, RSDCommonOutliers->common2Filename);
		strcat(tstring, " ");
		strcat(tstring, outputFilename);
		strcat(tstring, " > /dev/null 2>&1");
	}

	if(RSDMuStat!=NULL)
	{
		if(RSDMuStat->reportFP!=NULL)
		{
			fclose(RSDMuStat->reportFP);
			RSDMuStat->reportFP = NULL;
		}
	}

#ifdef _C1
	int ret = system(tstring);
	ret = ret;
	assert(ret!=-1);
#else
	FILE *fp;

	fp = popen(tstring, "r");
	assert(fp!=NULL);

	int ret = pclose(fp);
	assert(ret!=-1);
#endif	
}
