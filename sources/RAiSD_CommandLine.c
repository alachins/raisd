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
	fprintf(fp, " This is RAiSD version 1.3, released in July 2018.\n\n");

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
	fprintf(fp, "\t[-i]\n");
	fprintf(fp, "\t[-a INTEGER]\n");

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
	fprintf(fp, " -i\tImputes missing data per SNP (The default operation is to discard the SNP).\n");
	fprintf(fp, " -a\tProvides a seed for the random number generator.\n");


	fprintf(fp, "\n");
}

void RSDVersions(FILE * fp)
{
	int releaseIndex = 0;
	int majorIndex = 1;
	int minorIndex = 0;

	fprintf(fp, " %d. RAiSD v%d.%d (Jun  9, 2017): first release\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Mar  7, 2018): MAF threshold option\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Mar 28, 2018): mbs format with -b\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Jul 18, 2018): -i to impute N per SNP, -a for rand seed\n", releaseIndex++, majorIndex, minorIndex++);

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
		
		if(!strcmp(argv[i], "-i")) // To impute N per SNP
		{ 
			RSDCommandLine->imputePerSNP = 1;
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

