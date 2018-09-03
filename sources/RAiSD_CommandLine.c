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

void RSDHelp (FILE * fp)
{
	fprintf(fp, " This is RAiSD version 1.6, released in September 2018.\n\n");

	fprintf(fp, " RAiSD");

	fprintf(fp, "\n");
	fprintf(fp, "\t -n STRING\n");
	fprintf(fp, "\t -I STRING\n");
	fprintf(fp, "\t[-L INTEGER]\n");
	fprintf(fp, "\t[-h]\n");
	fprintf(fp, "\t[-v]\n");
	fprintf(fp, "\t[-f]\n");
	fprintf(fp, "\t[-s]\n");
	fprintf(fp, "\t[-t]\n");
	fprintf(fp, "\t[-p]\n");
	fprintf(fp, "\t[-S STRING]\n");
	fprintf(fp, "\t[-T INTEGER]\n");
	fprintf(fp, "\t[-d INTEGER]\n");
	fprintf(fp, "\t[-k FLOATING-POINT]\n");
	fprintf(fp, "\t[-l FLOATING-POINT]\n");
	fprintf(fp, "\t[-m FLOATING-POINT]\n");
	fprintf(fp, "\t[-b]\n");
	fprintf(fp, "\t[-a INTEGER]\n");
	fprintf(fp, "\t[-M 0|1|2|3]\n");
	fprintf(fp, "\t[-O]\n");
	fprintf(fp, "\t[-R]\n");
	fprintf(fp, "\t[-P]\n");

	fprintf(fp, "\n");	
	fprintf(fp, " -n\tProvides a unique run ID that is used to name the output files, i.e., the info file and the report(s).\n");
	fprintf(fp, " -I\tProvides the path to the input file, which can be either in ms or in vcf format.\n");
	fprintf(fp, " -L\tProvides the size of the region in basepairs for ms files. If known, it can be used for vcf as well, leading to faster processing.\n");
	fprintf(fp, " -h\tPrints this help message.\n");
	fprintf(fp, " -v\tPrints version information.\n");
	fprintf(fp, " -f\tOverwrites existing run files under the same run ID.\n");
	fprintf(fp, " -s\tGenerates a separate report file per set.\n");
	fprintf(fp, " -t\tRemoves the set separator symbol from the report(s).\n");
	fprintf(fp, " -p\tGenerates the output file RAiSD_Samples.STRING, where STRING is the run ID, comprising a list of samples in the input file (supported only with VCF).\n");
	fprintf(fp, " -S\tProvides the path to the list of samples to be processed (supported only with VCF).\n");
	fprintf(fp, " -T\tProvides the selection target (in basepairs) and calculates the average distance (over all datasets in the input file) between the selection target and the reported locations.\n");
	fprintf(fp, " -d\tProvides a maximum distance (in base pairs, from the selection target) to calculate success rate in terms of reported locations in the proximity of the target of selection (provided via -T).\n");
	fprintf(fp, " -k\tProvides the false positive rate (e.g., 0.05) to report the corresponding reported score after sorting the reported locations for all the datasets in the input file.\n");
	fprintf(fp, " -l\tProvides the threshold score, reported by a previous run using a false positive rate (e.g., 0.05, via -k) to report the true positive rate.\n");
	fprintf(fp, " -m\tProvides the threshold value for excluding SNPs with minor allele frequency < threshold (0.0-1.0).\n");
	fprintf(fp, " -b\tIndicates that the input file is in mbs format.\n");
	fprintf(fp, " -a\tProvides a seed for the random number generator.\n");	
	fprintf(fp, " -M\tIndicates the missing-data handling strategy (0: discards SNP (default), 1: imputes N per SNP, 2: represents N through a mask, 3: ignores allele pairs with N).\n");
	fprintf(fp, " -O\tShows progress on the display device (at snp set granularity).\n");
	fprintf(fp, " -R\tIncludes additional information (window start and end, and the mu-statistic factors for variation, SFS, and LD) in the report file.\n");
	fprintf(fp, " -P\tGenerates four plots (for the three mu-statistic factors and the final score) in one PDF file per set of SNPs in the input file using Rscript (activates -s, -t, and -R).\n");

	fprintf(fp, "\n");
}

void RSDVersions(FILE * fp)
{
	int releaseIndex = 0;
	int majorIndex = 1;
	int minorIndex = 0;

	fprintf(fp, " %d. RAiSD v%d.%d (Jun  9, 2017): first release\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Mar  7, 2018): -m to provide a MAF threshold\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Mar 28, 2018): -b to suppoert the mbs format\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Jul 18, 2018): -i to impute N per SNP, -a for rand seed\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Aug  3, 2018): -M to handle missing data with 4 strategies (removed -i)\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Aug  4, 2018): -R to include additional information in the report file\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Sep  3, 2018): -P to create plots per set of SNPs with Rscript\n", releaseIndex++, majorIndex, minorIndex++);

	majorIndex++;
}

RSDCommandLine_t * RSDCommandLine_new(void)
{
	RSDCommandLine_t * cl = NULL;
	cl = (RSDCommandLine_t *) malloc(sizeof(RSDCommandLine_t));
	assert(cl!=NULL);
	return cl;
}

void RSDCommandLine_free(RSDCommandLine_t * cl)
{
	assert(cl!=NULL);
	free(cl);
}

void RSDCommandLine_init(RSDCommandLine_t * RSDCommandLine)
{
	strncpy(RSDCommandLine->runName, "\0", STRING_SIZE);
	strncpy(RSDCommandLine->inputFileName, "\0", STRING_SIZE);
	RSDCommandLine->regionLength = 0ull;
	RSDCommandLine->overwriteOutput = 0;
	RSDCommandLine->splitOutput = 0;
	RSDCommandLine->setSeparator = 1;
	RSDCommandLine->printSampleList = 0;
	RSDCommandLine->maf = 0.0;
	strncpy(RSDCommandLine->sampleFileName, "\0", STRING_SIZE);
	RSDCommandLine->mbs = 0;
	RSDCommandLine->imputePerSNP = 0;
	RSDCommandLine->createPatternPoolMask = 0;
	RSDCommandLine->patternPoolMaskMode = 0;
	RSDCommandLine->displayProgress = 0;
	RSDCommandLine->fullReport = 0;
	RSDCommandLine->createPlot = 0;
	RSDCommandLine->muThreshold = 0.0;
}

void RSDCommandLine_load(RSDCommandLine_t * RSDCommandLine, int argc, char ** argv)
{
	int i;
	int info_exists = 0;
	char tstring [STRING_SIZE];
	for(i=1; i<argc; ++i)
	{
		if(!strcmp(argv[i], "-n")) 
		{ 
			if (i!=argc-1 && argv[i+1][0]!='-')
				strcpy(RSDCommandLine->runName, argv[++i]);
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}
			
			strcpy(tstring, "RAiSD_Info.");
			strcat(tstring, RSDCommandLine->runName);
			RAiSD_Info_FP = fopen(tstring, "r");
			if(RAiSD_Info_FP!=NULL)
			{
				fclose(RAiSD_Info_FP);
				info_exists = 1;
			}
			else
			{
				RAiSD_Info_FP = fopen(tstring, "w");
				assert(RAiSD_Info_FP!=NULL);
			}
			continue;
		}

		if(!strcmp(argv[i], "-I")) 
		{ 
			if (i!=argc-1 && argv[i+1][0]!='-')
				strcpy(RSDCommandLine->inputFileName, argv[++i]);
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-L")) 
		{ 
			if (i!=argc-1)
			{
				double len = atof(argv[++i]);
				RSDCommandLine->regionLength = (uint64_t) len;
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-help") || !strcmp(argv[i], "help") || !strcmp(argv[i], "-h")) 
		{ 
			RSD_header(stdout);
			RSDHelp(stdout);
			exit(0);
		}

		if(!strcmp(argv[i], "-v")) // version information
		{ 
			RSDVersions(stdout);
			exit(0);
		}

		if(!strcmp(argv[i], "-f")) // Force overwrite output report
		{ 
			RSDCommandLine->overwriteOutput = 1;
			continue;
		}

		if(!strcmp(argv[i], "-s")) // Split output reports per set
		{ 
			RSDCommandLine->splitOutput = 1;
			continue;
		}

		if(!strcmp(argv[i], "-t")) // Remove separator symbol
		{ 
			RSDCommandLine->setSeparator = 0;
			continue;
		}

		if(!strcmp(argv[i], "-p")) // Print sample list (currently only for VCF)
		{ 
			RSDCommandLine->printSampleList = 1;
			continue;
		}

		if(!strcmp(argv[i], "-S")) // Print sample list (currently only for VCF)
		{ 
			if (i!=argc-1 && argv[i+1][0]!='-')
				strcpy(RSDCommandLine->sampleFileName, argv[++i]);
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}
			continue;
		}

		if(!strcmp(argv[i], "-m")) // To provide a threshold for MAF
		{ 
			if (i!=argc-1)
			{
				RSDCommandLine->maf = (double)atof(argv[++i]);
				if(RSDCommandLine->maf<0.0 || RSDCommandLine->maf>1.0)
				{
					fprintf(stderr, "\nERROR: Invalid MAF value (valid: 0.0-1.0)\n\n");
					exit(0);
				}
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-b")) // To specify mbs format
		{ 
			RSDCommandLine->mbs = 1;
			continue;
		}	

		if(!strcmp(argv[i], "-M")) 
		{ 
			if (i!=argc-1)
			{	
				if((strlen(argv[i+1])!=1) || ((strlen(argv[i+1])==1)&&(argv[i+1][0]!='0' && argv[i+1][0]!='1' && argv[i+1][0]!='2' && argv[i+1][0]!='3')))// && argv[i+1][0]!='3')
				{
					fprintf(stderr, "\nERROR: Invalid argument after %s\n\n",argv[i]);
					exit(0);	
				}

				int mode = atoi(argv[++i]);
				switch(mode)
				{
					case 1:
						RSDCommandLine->imputePerSNP = 1;
						RSDCommandLine->createPatternPoolMask = 0;
						RSDCommandLine->patternPoolMaskMode = 0;
					break;
					case 2:
						RSDCommandLine->imputePerSNP = 0;
						RSDCommandLine->createPatternPoolMask = 1;
						RSDCommandLine->patternPoolMaskMode = 0;
					break;
					case 3:
						RSDCommandLine->imputePerSNP = 0;
						RSDCommandLine->createPatternPoolMask = 1;
						RSDCommandLine->patternPoolMaskMode = 1;
					break;
					default:
						RSDCommandLine->imputePerSNP = 0;
						RSDCommandLine->createPatternPoolMask = 0;
						RSDCommandLine->patternPoolMaskMode = 0;
					break;
				}				
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-O")) 
		{ 
			RSDCommandLine->displayProgress = 1;
			continue;
		}

		if(!strcmp(argv[i], "-R")) 
		{ 
			RSDCommandLine->fullReport = 1;
			continue;
		}


		if(!strcmp(argv[i], "-a")) 
		{ 
			if (i!=argc-1)
			{
				int seed = atoi(argv[++i]);
				srand((unsigned int)seed);
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-P")) 
		{ 
			RSDCommandLine->createPlot = 1;
			RSDCommandLine->splitOutput = 1; // activating output splitting
			RSDCommandLine->fullReport = 1; // activating full report
			RSDCommandLine->setSeparator = 0; // remove separator symbol

			if(RSDPlot_checkRscript()!=0)
			{
				fprintf(stderr, "\nERROR: Rscript is not installed, required by %s for plotting with R\n\n",argv[i]);
				exit(0);
			}			
			
			continue;
		}

		/*if(!strcmp(argv[i], "-P")) // To provide a threshold for plotting/reporting
		{ 
			if (i!=argc-1)
			{
				RSDCommandLine->plotThreshold = (double)atof(argv[++i]);
				if(RSDCommandLine->plotThreshold<0.0)
				{
					fprintf(stderr, "\nERROR: Invalid threshold value (valid: >=0.0)\n\n");
					exit(0);
				}
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}*/

		/* Testing */
		if(!strcmp(argv[i], "-T")) 
		{ 
			if (i!=argc-1)
			{
				double tar = atof(argv[++i]);
				selectionTarget = (uint64_t)tar;
			}			
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-d")) 
		{ 
			if (i!=argc-1)
			{
				double dist = atof(argv[++i]);
				selectionTargetDThreshold = (uint64_t)dist;
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-k")) 
		{ 
			if (i!=argc-1)
				fpr_loc = (double)atof(argv[++i]);
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-l")) 
		{ 
			if (i!=argc-1)
				tpr_thres = (double)atof(argv[++i]);
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		

		/*if(!strcmp(argv[i], "-set")) 
		{ 
			if (i!=argc-1)
				setIndexValid = atoi(argv[++i]);
			else
			{
				fprintf(stderr, "ERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}*/

		fprintf(stderr, "\nERROR: Unrecognized input parameter %s\n\n",argv[i]);
		exit(0);
	}

	// Checks
	if(info_exists==1) 
	{
		if(RSDCommandLine->overwriteOutput==0)
		{
			fprintf(stderr, "\nERROR: Output info file %s exists. Use -f to overwrite it.\n\n", tstring);
			exit(0);
		}
		else
		{
			RAiSD_Info_FP = fopen(tstring, "w");
			assert(RAiSD_Info_FP!=NULL);
		}
	}

	if(!strcmp(RSDCommandLine->runName, "\0"))
	{
		fprintf(stderr, "\nERROR: Missing required input parameter -n\n\n");
		exit(0);	
	}

	if(!strcmp(RSDCommandLine->inputFileName, "\0"))
	{
		fprintf(stderr, "\nERROR: Missing required input parameter -I\n\n");
		exit(0);	
	}
	else
	{
		FILE * inputFileExists = fopen(RSDCommandLine->inputFileName, "r");
		assert(inputFileExists!=NULL);
		fclose(inputFileExists);
	}

	if(RSDCommandLine->createPatternPoolMask==1)
		RSDCommandLine->imputePerSNP = 0;

}

void RSDCommandLine_print(int argc, char ** argv, FILE * fpOut)
{
	fprintf(fpOut, " Command: ");
	int i;
	for(i=0; i<argc; ++i)
	{
		fprintf(fpOut, "%s ", argv[i]);
	}
	fprintf(fpOut, "\n");

	fflush(fpOut);
}

