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

void flagCheck (char ** argv, int i, int * flagVector, int flagIndex);

void RSDHelp (FILE * fp)
{
	fprintf(fp, " This is RAiSD version %d.%d, released in %s %d.\n\n", MAJOR_VERSION, MINOR_VERSION, RELEASE_MONTH, RELEASE_YEAR);

	fprintf(fp, " RAiSD");

	fprintf(fp, "\n");
	fprintf(fp, "\t -n STRING\n");
	fprintf(fp, "\t -I STRING\n");

	fprintf(fp, "\n\t--- SNP and SAMPLE HANDLING\n\n");
	fprintf(fp, "\t[-L INTEGER]\n"); 	
	fprintf(fp, "\t[-S STRING]\n");
	fprintf(fp, "\t[-m FLOAT]\n");
	fprintf(fp, "\t[-M 0|1|2|3]\n");
	fprintf(fp, "\t[-y INTEGER]\n");
	fprintf(fp, "\t[-X STRING]\n");
	fprintf(fp, "\t[-B INTEGER INTEGER]\n");
	fprintf(fp, "\t[-o]\n");

	fprintf(fp, "\n\t--- SLIDING WINDOW and MU STATISTIC\n\n");
	fprintf(fp, "\t[-w INTEGER]\n");
	fprintf(fp, "\t[-c INTEGER]\n");

	fprintf(fp, "\n\t--- STANDARD OUTPUT and REPORTS\n\n");
	fprintf(fp, "\t[-f]\n");
	fprintf(fp, "\t[-s]\n");
	fprintf(fp, "\t[-t]\n");
	fprintf(fp, "\t[-p]\n");
	fprintf(fp, "\t[-O]\n");
	fprintf(fp, "\t[-R]\n");
	fprintf(fp, "\t[-P]\n");
	fprintf(fp, "\t[-D]\n");
	fprintf(fp, "\t[-A FLOAT]\n");

	fprintf(fp, "\n\t--- ACCURACY and SENSITIVITY EVALUATION\n\n");	
	fprintf(fp, "\t[-T INTEGER]\n");
	fprintf(fp, "\t[-d INTEGER]\n");
	fprintf(fp, "\t[-k FLOAT]\n");
	fprintf(fp, "\t[-l FLOAT]\n");

	fprintf(fp, "\n\t--- ADDITIONAL EXECUTION PARAMETERS\n\n");	
	fprintf(fp, "\t[-b]\n");
	fprintf(fp, "\t[-a INTEGER]\n");
	
	fprintf(fp, "\n\t--- HELP and VERSION NOTES\n\n");
	fprintf(fp, "\t[-h]\n");
	fprintf(fp, "\t[-v]\n");

	fprintf(fp, "\n\n");
	fprintf(fp, " DETAILED DESCRIPTION\n\n");	
	fprintf(fp, "\t-n\tProvides a unique run ID that is used to name the output files, i.e., the info file and the report(s).\n");
	fprintf(fp, "\t-I\tProvides the path to the input file, which can be either in ms or in vcf format.\n");

	fprintf(fp, "\n\t--- SNP and SAMPLE HANDLING\n\n");
	fprintf(fp, "\t-L\tProvides the size of the region in basepairs for ms files. See -B option for vcf files.\n");
	fprintf(fp, "\t-S\tProvides the path to the list of samples to be processed (supported only with VCF).\n");
	fprintf(fp, "\t-m\tProvides the threshold value for excluding SNPs with minor allele frequency < threshold (0.0-1.0).\n");
	fprintf(fp, "\t-M\tIndicates the missing-data handling strategy (0: discards SNP (default), 1: imputes N per SNP,\n\t\t2: represents N through a mask, 3: ignores allele pairs with N).\n");
	fprintf(fp, "\t-y\tProvides the ploidy (integer value), used to correctly represent missing data.\n");
	fprintf(fp, "\t-X\tProvides the path to a tab-delimited file that contains regions per chromosome (one per line) to be\n\t\texcluded from the analysis (Format: chromosome [tab] regionStart [tab] regionStop).\n");
	fprintf(fp, "\t-B\tProvides the chromosome size in basepairs (first INTEGER) and SNPs (second INTEGER) for vcf files that\n\t\tcontain a single chromosome. If -B is not provided, or the input vcf file contains multiple chromosomes,\n\t\tRAiSD will determine the respective values by parsing each chromosome in its entirety before processing,\n\t\twhich will lead to slightly longer overall execution time.\n");
	fprintf(fp, "\t-o\tEnables dataset check and ordering of the input vcf file (only unzipped vcf files are supported).\n");

	fprintf(fp, "\n\t--- SLIDING WINDOW and MU STATISTIC\n\n");
	fprintf(fp, "\t-w\tProvides the window size (integer value). The default value is 50 (empirically determined).\n");
	fprintf(fp, "\t-c\tProvides the slack for the SFS edges to be used for the calculation of mu_SFS. The default value is 1\n\t\t(singletons and S-1 snp class, where S is the sample size).\n");

	fprintf(fp, "\n\t--- STANDARD OUTPUT and REPORTS\n\n");
	fprintf(fp, "\t-f\tOverwrites existing run files under the same run ID.\n");
	fprintf(fp, "\t-s\tGenerates a separate report file per set.\n");
	fprintf(fp, "\t-t\tRemoves the set separator symbol from the report(s).\n");
	fprintf(fp, "\t-p\tGenerates the output file RAiSD_Samples.STRING, where STRING is the run ID, comprising a list of samples\n\t\tin the input file (supported only with VCF).\n");
	fprintf(fp, "\t-O\tShows progress on the display device (at snp set granularity).\n");
	fprintf(fp, "\t-R\tIncludes additional information in the report file(s), i.e., window start and end, and the mu-statistic\n\t\tfactors for variation, SFS, and LD.\n");
	fprintf(fp, "\t-P\tGenerates four plots (for the three mu-statistic factors and the final score) in one PDF file per set of\n\t\tSNPs in the input file using Rscript (activates -s, -t, and -R).\n");
	fprintf(fp, "\t-D\tGenerates a site report, e.g., total, discarded, imputed etc.\n");
	fprintf(fp, "\t-A\tProvides a probability value to be used for the quantile function in R and generates a Manhattan plot for\n\t\tthe final mu-statistic score using Rscript (activates -s, -t, and -R).\n");

	fprintf(fp, "\n\t--- ACCURACY and SENSITIVITY EVALUATION\n\n");
	fprintf(fp, "\t-T\tProvides the selection target (in basepairs) and calculates the average distance (over all datasets in the\n\t\tinput file) between the selection target and the reported locations.\n");
	fprintf(fp, "\t-d\tProvides a maximum distance from the selection target (in base pairs) to calculate success rates,\n\t\ti.e., reported locations in the proximity of the target of selection (provided via -T).\n");
	fprintf(fp, "\t-k\tProvides the false positive rate (e.g., 0.05) to report the corresponding reported score after sorting\n\t\tthe reported locations for all the datasets in the input file.\n");
	fprintf(fp, "\t-l\tProvides the threshold score, reported by a previous run using a false positive rate (e.g., 0.05, via -k)\n\t\tto report the true positive rate.\n");
	

	fprintf(fp, "\n\t--- ADDITIONAL EXECUTION PARAMETERS\n\n");
	fprintf(fp, "\t-b\tIndicates that the input file is in mbs format.\n");
	fprintf(fp, "\t-a\tProvides a seed for the random number generator.\n");	

	fprintf(fp, "\n\t--- HELP and VERSION NOTES\n\n");
	fprintf(fp, "\t-h\tPrints this help message.\n");
	fprintf(fp, "\t-v\tPrints version information.\n");

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
	fprintf(fp, " %d. RAiSD v%d.%d (Oct  2, 2018): -y for ploidy, -D for site report, fixed a bug in the plotting routine\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Dec 31, 2018): MakefileZLIB to parse VCF files in gzip file format (requires the zlib library)\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Apr 27, 2019): -w to set the window size (default 50), -c to set the SFS slack for the mu_SFS\n", releaseIndex++, majorIndex, minorIndex++);

	majorIndex++; minorIndex=0;

	fprintf(fp, " %d. RAiSD v%d.%d (May 15, 2019): -A to create Manhattan plots, scale factors for muVar and muSFS to yield comparable scores among different chromosomes\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Jan 21, 2020): Parser for unordered VCF files (e.g., extracted from DArTseq genotyping reports). The ordered VCF file is also created.\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Jan 22, 2020): Added missing field (discarded monomorphic sites) in the site report (Dataset.c file) for missing-data strategies M=1,2, or 3.\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Jan 23, 2020): -X to exclude regions per chromosome from the analysis.\n", releaseIndex++, majorIndex, minorIndex++);
	fprintf(fp, " %d. RAiSD v%d.%d (Jan 30, 2020): -B for chromosome length and SNP size. Fixed bug with the memory-reduction optimization for large chromosomes. -o to request vcf ordering and generation.\n", releaseIndex++, majorIndex, minorIndex++);

	// TODO: add message here for outoforder vcf
	// TODO: add message here for adding the monomorphic count in the site report
	// TODO: need to fix the inflated values bug
}

RSDCommandLine_t * RSDCommandLine_new(void)
{
	RSDCommandLine_t * cl = NULL;
	cl = (RSDCommandLine_t *) rsd_malloc(sizeof(RSDCommandLine_t));
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
	RSDCommandLine->regionSNPs = 0ull;
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
	RSDCommandLine->createMPlot = 0;
	strncpy(RSDCommandLine->manhattanThreshold, "\0", STRING_SIZE);
	RSDCommandLine->muThreshold = 0.0;
	RSDCommandLine->ploidy = 2; // default
	RSDCommandLine->displayDiscardedReport = 0;
	RSDCommandLine->windowSize = DEFAULT_WINDOW_SIZE;
	RSDCommandLine->sfsSlack = 1; // singletons, and S-1 snp class (S is the sample size)
	strncpy(RSDCommandLine->excludeRegionsFile, "\0", STRING_SIZE);
	RSDCommandLine->orderVCF = 0; 
}

void flagCheck (char ** argv, int i, int * flagVector, int flagIndex)
{
	if(flagVector[flagIndex]!=0)
	{
		fprintf(stderr, "\nERROR: Flag %s is given more than once!\n\n",argv[i]);
		exit(0);
	}

	flagVector[flagIndex]=1;
}

void RSDCommandLine_load(RSDCommandLine_t * RSDCommandLine, int argc, char ** argv)
{
	int i;
	int info_exists = 0;
	char tstring [STRING_SIZE], tstring2[STRING_SIZE];
	int * flagVector = (int*)calloc(MAX_COMMANDLINE_FLAGS, sizeof(int));
	assert(flagVector!=NULL);

	for(i=1; i<argc; ++i)
	{
		if(!strcmp(argv[i], "-n")) 
		{ 
			flagCheck (argv, i, flagVector, 0);

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
			flagCheck (argv, i, flagVector, 1);

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
			flagCheck (argv, i, flagVector, 2);

			if(flagVector[B_FLAG_INDEX])
			{
				fprintf(stderr, "\nERROR: Argument %s cannot be used if -B is already provided!\n\n",argv[i]);
				exit(0);
			}

			if (i!=argc-1 && argv[i+1][0]!='-')
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
			flagCheck (argv, i, flagVector, 3);

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
			flagCheck (argv, i, flagVector, 4);

			if (i!=argc-1 && argv[i+1][0]!='-')
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

		if(!strcmp(argv[i], "-y")) 
		{ 
			flagCheck (argv, i, flagVector, 5);

			if (i!=argc-1 && argv[i+1][0]!='-')
			{
				int ploidy = atoi(argv[++i]);
				if(ploidy<=0)
				{
					fprintf(stderr, "\nERROR: Invalid argument after %s\n\n",argv[i]);
					exit(0);
				}
				RSDCommandLine->ploidy = ploidy;
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}
			continue;
		}	

		if(!strcmp(argv[i], "-M")) 
		{ 
			flagCheck (argv, i, flagVector, 6);

			if (i!=argc-1 && argv[i+1][0]!='-')
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

		if(!strcmp(argv[i], "-D")) 
		{ 
			flagCheck (argv, i, flagVector, 7);

			RSDCommandLine->displayDiscardedReport = 1;
			
			tstring2[0]='\0';
			strcpy(tstring2, "RAiSD_SiteReport.");
			strcat(tstring2, RSDCommandLine->runName);

			RAiSD_SiteReport_FP = fopen(tstring2, "w");
			assert(RAiSD_SiteReport_FP!=NULL);

			continue;
		}


		if(!strcmp(argv[i], "-a")) 
		{ 
			flagCheck (argv, i, flagVector, 8);

			if (i!=argc-1 && argv[i+1][0]!='-')
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

		if(!strcmp(argv[i], "-w")) 
		{ 
			flagCheck (argv, i, flagVector, 9);

			if (i!=argc-1 && argv[i+1][0]!='-')
			{
				double windowSize = atof(argv[++i]);
				int64_t windowSizeInt = (int64_t)windowSize;
				if((windowSize<MIN_WINDOW_SIZE) || ((windowSizeInt&1)==1))
				{
					fprintf(stderr, "\nERROR: Invalid window size (valid: even number >= %d)\n\n", MIN_WINDOW_SIZE);
					exit(0);
				}
				RSDCommandLine->windowSize = (int64_t)windowSize;
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}


		if(!strcmp(argv[i], "-c")) 
		{ 
			flagCheck (argv, i, flagVector, 10);

			if (i!=argc-1 && argv[i+1][0]!='-')
			{
				double sfsSlack = atof(argv[++i]);
				if(sfsSlack<1)
				{
					fprintf(stderr, "\nERROR: Invalid sfs slack value (valid: >= 1)\n\n");
					exit(0);
				}
				RSDCommandLine->sfsSlack = (int64_t)sfsSlack;
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-A")) 
		{ 
			flagCheck (argv, i, flagVector, 11);

			double prob = 0.9999; // default

			if (i!=argc-1 && argv[i+1][0]!='-')
			{
				prob = (double)atof(argv[++i]);
				if(prob<0.0 || prob>1.0)
				{
					fprintf(stderr, "\nERROR: Invalid probability value (valid: 0.0-1.0)\n\n");
					exit(0);
				}
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			strcpy(RSDCommandLine->manhattanThreshold, argv[i]);

			RSDCommandLine->createMPlot = 1;
			RSDCommandLine->splitOutput = 1; // activating output splitting
			RSDCommandLine->fullReport = 1; // activating full report
			RSDCommandLine->setSeparator = 0; // remove separator symbol

			if(RSDPlot_checkRscript()!=0)
			{
				fprintf(stderr, "\nERROR: Rscript is not installed, required by %s for plotting with R\n\n",argv[i]);
				exit(0);
			}

			tstring2[0]='\0';
			strcpy(tstring2, "RAiSD_ReportList.txt"); 
			RSDPlot_createReportListName (RSDCommandLine, tstring2);

			RAiSD_ReportList_FP = fopen(tstring2, "w");
			assert(RAiSD_ReportList_FP!=NULL);	
			
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
			flagCheck (argv, i, flagVector, 12);

			if (i!=argc-1 && argv[i+1][0]!='-')
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
			flagCheck (argv, i, flagVector, 13);

			if (i!=argc-1 && argv[i+1][0]!='-')
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
			flagCheck (argv, i, flagVector, 14);

			if (i!=argc-1 && argv[i+1][0]!='-')
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
			flagCheck (argv, i, flagVector, 15);

			if (i!=argc-1 && argv[i+1][0]!='-')
				tpr_thres = (double)atof(argv[++i]);
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-X")) 
		{ 
			flagCheck (argv, i, flagVector, 16);

			if (i!=argc-1 && argv[i+1][0]!='-')
				strcpy(RSDCommandLine->excludeRegionsFile, argv[++i]);
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-B")) 
		{ 

			if(flagVector[L_FLAG_INDEX])
			{
				fprintf(stderr, "\nERROR: Argument %s cannot be used if -L is already provided!\n\n",argv[i]);
				exit(0);
			}
			
			flagCheck (argv, i, flagVector, 17);

			if ((i!=argc-1 && argv[i+1][0]!='-')&&(i!=argc-2 && argv[i+2][0]!='-'))
			{
				double len = atof(argv[++i]);
				RSDCommandLine->regionLength = (uint64_t) len;

				double snps = atof(argv[++i]);
				RSDCommandLine->regionSNPs = (uint64_t) snps;
			}
			else
			{
				fprintf(stderr, "\nERROR: Missing argument after %s\n\n",argv[i]);
				exit(0);	
			}

			continue;
		}

		if(!strcmp(argv[i], "-o")) 
		{ 
			RSDCommandLine->orderVCF = 1;
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
		assert(inputFileExists!=NULL); // file exists check
		fclose(inputFileExists);
	}

	if(strcmp(RSDCommandLine->excludeRegionsFile, "\0"))
	{
		FILE * excludeRegionsFileExists = fopen(RSDCommandLine->excludeRegionsFile, "r");
		assert(excludeRegionsFileExists!=NULL); // file exists check
		fclose(excludeRegionsFileExists);
	}

	if(RSDCommandLine->createPatternPoolMask==1)
		RSDCommandLine->imputePerSNP = 0;

	free(flagVector);
}

void RSDCommandLine_print(int argc, char ** argv, FILE * fpOut)
{
	if(fpOut==NULL)
		return;

	fprintf(fpOut, " Command: ");
	int i;
	for(i=0; i<argc; ++i)
	{
		fprintf(fpOut, "%s ", argv[i]);
	}
	fprintf(fpOut, "\n");

	fflush(fpOut);
}

void RSDCommandLine_printWarnings (RSDCommandLine_t * RSDCommandLine, int argc, char ** argv, void * RSDDataset, FILE * fpOut)
{
	if(fpOut==NULL)
		return;

	assert(RSDCommandLine!=NULL);
	assert(RSDDataset!=NULL);

	RSDDataset_t * myRSDDataset = (RSDDataset_t *)RSDDataset;
	
	int i=0;

	for(i=1; i<argc; ++i)
	{
		if(!strcmp(argv[i], "-L") && !strcmp(myRSDDataset->inputFileFormat, "vcf.gz")) 
		{ 
			fprintf(fpOut, "\n WARNING: Argument -L is not to be used with vcf.gz files and will be ignored!\n");
			fflush(fpOut);
		}

		if(!strcmp(argv[i], "-L") && !strcmp(myRSDDataset->inputFileFormat, "vcf")) 
		{ 
			fprintf(fpOut, "\n WARNING: Argument -L is not to be used with vcf files and will be ignored!\n");
			fflush(fpOut);
		}

		if(!strcmp(argv[i], "-B") && !strcmp(myRSDDataset->inputFileFormat, "ms")) 
		{ 
			fprintf(fpOut, "\n WARNING: Argument -B is not to be used with ms files and will be ignored!\n");
			fflush(fpOut);
		}
	}
}

