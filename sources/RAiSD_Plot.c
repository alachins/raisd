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
	if(RSDCommandLine->createPlot==0)
		return;

	fprintf(fpOut, " Rscript: ");

	FILE * fp;
	fp = fopen("RscriptVersion.txt", "r");
	assert(fp==NULL);
	
#ifdef _C1
	int ret = system("Rscript --version > /dev/null 2>RscriptVersion.txt");
	assert(ret!=-1);
#else
	fp = popen("Rscript --version > /dev/null 2>RscriptVersion.txt", "r");
	assert(fp!=NULL);

	int ret = pclose(fp);
	assert(ret!=-1);
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

void RSDPlot_generateRscript (RSDCommandLine_t * RSDCommandLine)
{
	if(RSDCommandLine->createPlot==0)
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
plot(mup, mu1, col=\"darkgray\", pch=16, ylab=bquote(mu ~ \"_var\"), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"_var curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
plot(mup, mu2, col=\"darkgray\", pch=16, ylab=bquote(mu ~ \"_sfs\"), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"_sfs curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
plot(mup, mu3, col=\"darkgray\", pch=16, ylab=bquote(mu ~ \"_ld\"), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"_ld curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
plot(mup, mu, col=\"black\", pch=16, ylab=bquote(mu), xlab=\"Chromosome position (kb)\", main=bquote(mu ~ \"curve for\" ~ .(RUNNAME) ~ \".\" ~ .(ID))) \n \
dev.off() \n";

	fprintf(fp,"%s", tscript);

	fclose(fp);	
}

void RSDPlot_removeRscript (RSDCommandLine_t * RSDCommandLine)
{
	char scriptName[STRING_SIZE]="RSDPlot.R";
	RSDPlot_createRscriptName (RSDCommandLine, scriptName);

	FILE * fp;
	fp = fopen(scriptName, "r"); // TODO
	if(fp!=NULL)
	{
		fclose(fp);
		int ret = remove(scriptName);
		assert(ret==0);	
	}
	fp=NULL;
}

void RSDPlot_createPlot (RSDCommandLine_t * RSDCommandLine, RSDDataset_t * RSDDataset)
{
	char scriptName[STRING_SIZE]="RSDPlot.R";
	RSDPlot_createRscriptName (RSDCommandLine, scriptName);

	char tstring[STRING_SIZE];
	strcpy(tstring, "Rscript ");
	strcat(tstring, scriptName);
	strcat(tstring, " ");	
	strcat(tstring, RSDCommandLine->runName);
	strcat(tstring, " ");
	strcat(tstring, RSDDataset->setID);
	strcat(tstring, " > /dev/null 2>&1");

#ifdef _C1
	int ret = system(tstring);
	assert(ret!=-1);
#else
	FILE *fp;

	fp = popen(tstring, "r");
	assert(fp!=NULL);

	int ret = pclose(fp);
	assert(ret!=-1);
#endif	
}
