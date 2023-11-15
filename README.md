
[RAiSD: software to detect positive selection based on multiple signatures of a selective sweep and SNP vectors](https://www.nature.com/articles/s42003-018-0085-8)
===============================================

Authors: Nikolaos Alachiotis (n.alachiotis@gmail.com) and Pavlos Pavlidis (pavlidisp@gmail.com)

First release: 9/6/2017 

Last update: 6/8/2020       

Version: 2.9

About
-----

RAiSD (Raised Accuracy in Sweep Detection) is a stand-alone software implementation of the μ statistic for selective sweep detection. Unlike existing implementations, including our previously released tools (SweeD and OmegaPlus), RAiSD scans whole-genome SNP data based on a composite evaluation scheme that captures multiple sweep signatures at once. 

The main article describing RAiSD and the μ statistic is published in Communications Biology:

1. RAiSD detects positive selection based on multiple signatures of a selective sweep and SNP vectors ([PDF](https://www.nature.com/articles/s42003-018-0085-8.pdf))  
   URL: https://www.nature.com/articles/s42003-018-0085-8   

Other related publications:

2. Accelerated Inference of Positive Selection on Whole Genomes ([PDF](http://kalman.mee.tcd.ie/fpl2018/content/pdfs/FPL2018-43iDzVTplcpussvbfIaaHz/2IXjcULtzLwXxGVD0BRuki/4etDSBv306dBHUVjNtmSyT.pdf))  
URL: https://ieeexplore.ieee.org/abstract/document/8533493  

3. RAiSD-X: A Fast and Accurate FPGA System for the Detection of Positive Selection in Thousands of Genomes ([PDF](https://github.com/alachins/raisd/blob/master/publications/RAiSD-X_TRETS2019.pdf))  
URL: https://dl.acm.org/doi/10.1145/3364225

4. Population genomics insights into the recent evolution of SARS-CoV-2 ([PDF](https://www.biorxiv.org/content/10.1101/2020.04.21.054122v1.full.pdf))  
URL: https://doi.org/10.1101/2020.04.21.054122  


Changes to the μ statistic (May 2019)
-------------------------------------

As of RAiSD version 2.0, the μ-statistic factors `μ_VAR` and `μ_SFS` are modified as follows.

`μ_VAR`: The previous version of `μ_VAR` is multiplied by the total number of SNPs in the chromosome. This makes `μ_VAR` independent of the chromosome length, while its values are now distributed around 1.0 under neutrality. 

`μ_SFS`: The previous version of `μ_SFS` is multiplied by the (sample size) correcting factor of Watterson's θ: `Sum_i=1^(n-1) 1/i`, where `n` is the sample size. The values of `μ_SFS` are now distributed around 1.0 under neutrality. 


Download and Compile
--------------------

The following commands can be used to download and compile the source code. 

    $ mkdir RAiSD
    $ cd RAiSD
    $ wget https://github.com/alachins/raisd/archive/master.zip
    $ unzip master.zip
    $ cd raisd-master
    $ ./install-RAiSD.sh
    
The executable is placed in the path RAiSD/raisd-master/bin/release. A link to the executable is placed in the installation folder, i.e., raisd-master. 

As of version 2.7 (or later), the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) needs to be installed to compile and run RAiSD. The following error message (or similar) will appear if GSL is not installed and found in the expected /usr/include/gsl directory. In this case, please contact your system administrator to install it (RECOMMENDED). 
	
	Error message:
	Creating directories
	Beginning release build
	Compiling: sources/RAiSD.c -> build/release/RAiSD.o
	In file included from sources/RAiSD.c:22:0:
	sources/RAiSD.h:30:27: fatal error: gsl/gsl_errno.h: No such file or directory
 	#include <gsl/gsl_errno.h>
         	                  ^
	compilation terminated.
	make[1]: *** [build/release/RAiSD.o] Error 1
	make: *** [release] Error 2

As an alternative, you can install GSL locally (NOT RECOMMENDED) with the following command (in the raisd-master installation folder): 

  	 $ ./install-RAiSD-GSL.sh
	 
This will download the latest GSL version from https://ftp.gnu.org/gnu/gsl/, unzip it, install it in the raisd-master/gsl directory, and compile RAiSD. If the executable is not eventually created (this process might take a few minutes), you might need to change the version numbers in the install-RAiSD-GSL.sh script. A folder named gsl-X.Y is already created in the raisd-master directory, where X and Y are the GSL version numbers. Replace the version numbers in the 5th line of the install-RAiSD-GSL.sh script with X and Y, and execute the installation script again. 

Finally, and while still in the raisd-master directory, use the following commands to allow RAiSD to access the locally installed GSL library. Then, you could execute the ./install-RAiSD.sh to compile RAiSD. 
	
	$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/gsl/lib
 	$ export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)/gsl/lib

These commands should be used from within the installation directory (raisd-master) in every new terminal before running RAiSD. Otherwise, the following error message will appear when launching RAiSD:

	./RAiSD: error while loading shared libraries: libgsl.so.25: cannot open shared object file: No such file or directory
	


Test Run
--------

To verify that RAiSD is installed correctly, a test run can be done with the following commands. These commands are going to download a dataset (one of the many that we used for evaluation purposes) and execute RAiSD to process 1,000 simulated sets of SNPs (genomic region size = 100000 bp, weak bottleneck, selective sweep at the center of the region).
    
    $ wget 139.91.162.50/raisd_data/d1.tar.gz
    $ tar -xvzf d1.tar.gz
    $ ./RAiSD -n test_run -I d1/msselection1.out -L 100000
    
Upon completion, the output files RAiSD_Info.test_run and RAiSD_Report.test_run are generated. 
    
In-tool Help
------------

RAiSD outputs a quick-reference help message that provides a short description for each of the supported command-line flags with the following command. 

    $ RAiSD -h
    
The quick-reference help message generated by the current RAiSD release is the following.

   
      This is RAiSD version 2.8, released in April 2020.

      RAiSD
	      	-n STRING
	      	-I STRING

		--- SNP and SAMPLE HANDLING
		
		[-L INTEGER]
		[-S STRING]
		[-m FLOAT]
		[-M 0|1|2|3]
		[-y INTEGER]
		[-X STRING]
		[-B INTEGER INTEGER]
		[-o]

		--- SLIDING WINDOW and MU STATISTIC
		
		[-w INTEGER]
		[-c INTEGER]
		[-G INTEGER]

		--- STANDARD OUTPUT and REPORTS
		
		[-f]
		[-s]
		[-t]
		[-p]
		[-O]
		[-R]
		[-P]
		[-D]
		[-A FLOAT]

		--- ACCURACY and SENSITIVITY EVALUATION
		
		[-T INTEGER]
		[-d INTEGER]
		[-k FLOAT]
		[-l FLOAT]

		--- ADDITIONAL EXECUTION PARAMETERS
		
		[-b]
		[-a INTEGER]
		
		--- FASTA-to-VCF CONVERSION PARAMETERS

		[-C STRING]
		[-C2 STRING]
		[-H STRING]
		[-E STRING]
		
		--- COMMON-OUTLIER ANALYSIS

		[-CO STRING INTEGER INTEGER] or [-CO STRING INTEGER INTEGER STRING INTEGER INTEGER]
		[-COT FLOAT]
		[-COD INTEGER]

		--- HELP and VERSION NOTES
		
		[-h]
		[-v]


      DETAILED DESCRIPTION

		-n	Provides a unique run ID that is used to name the output files, i.e., the info file and the report(s).
		-I	Provides the path to the input file, which can be either in ms or in vcf format.

		--- SNP and SAMPLE HANDLING

		-L	Provides the size of the region in basepairs for ms files. See -B option for vcf files.
		-S	Provides the path to the list of samples to be processed (supported only with VCF).
		-m	Provides the threshold value for excluding SNPs with minor allele frequency < threshold (0.0-1.0).
		-M	Indicates the missing-data handling strategy (0: discards SNP (default), 1: imputes N per SNP, 2: represents N through a mask, 3: ignores allele pairs with N).
		-y	Provides the ploidy (integer value), used to correctly represent missing data.
		-X	Provides the path to a tab-delimited file that contains regions per chromosome (one per line) to be excluded from the analysis (Format: chromosome [tab] regionStart [tab] regionStop).
		-B	Provides the chromosome size in basepairs (first INTEGER) and SNPs (second INTEGER) for vcf files that contain a single chromosome. If -B is not provided, or the input vcf file contains multiple chromosomes, RAiSD will determine the respective values by parsing each chromosome in its entirety before processing, which will lead to slightly longer overall execution time.
		-o	Enables dataset check and ordering of the input vcf file (only unzipped vcf files are supported).

		--- SLIDING WINDOW and MU STATISTIC

		-w	Provides the window size (integer value). The default value is 50 (empirically determined).
		-c	Provides the slack for the SFS edges to be used for the calculation of mu_SFS. The default value is 1 (singletons and S-1 snp class, where S is the sample size).
		-G	Provides the grid size to specify the total number of evaluation points across the data. When used, RAiSD reports mu statistic scores at equidistant locations between the first and last SNPs.

		--- STANDARD OUTPUT and REPORTS

		-f	Overwrites existing run files under the same run ID.
		-s	Generates a separate report file per set.
		-t	Removes the set separator symbol from the report(s).
		-p	Generates the output file RAiSD_Samples.STRING, where STRING is the run ID, comprising a list of samples in the input file (supported only with VCF).
		-O	Shows progress on the display device (at snp set granularity).
		-R	Includes additional information in the report file(s), i.e., window start and end, and the mu-statistic factors for variation, SFS, and LD.
		-P	Generates four plots (for the three mu-statistic factors and the final score) in one PDF file per set of SNPs in the input file using Rscript (activates -s, -t, and -R).
		-D	Generates a site report, e.g., total, discarded, imputed etc.
		-A	Provides a probability value to be used for the quantile function in R and generates a Manhattan plot for the final mu-statistic score using Rscript (activates -s, -t, and -R).

		--- ACCURACY and SENSITIVITY EVALUATION

		-T	Provides the selection target (in basepairs) and calculates the average distance (over all datasets in the input file) between the selection target and the reported locations.
		-d	Provides a maximum distance from the selection target (in base pairs) to calculate success rates, i.e., reported locations in the proximity of the target of selection (provided via -T).
		-k	Provides the false positive rate (e.g., 0.05) to report the corresponding reported score after sorting the reported locations for all the datasets in the input file.
		-l	Provides the threshold score, reported by a previous run using a false positive rate (e.g., 0.05, via -k) to report the true positive rate.

		--- ADDITIONAL EXECUTION PARAMETERS

		-b	Indicates that the input file is in mbs format.
		-a	Provides a seed for the random number generator.
		
		--- FASTA-to-VCF CONVERSION PARAMETERS

		-C	Provides the outgroup to be used for the ancestral states (REF field in VCF). The first ingroup sequence is used if the outgroup is not given or found.
		-C2	Provides a second outgroup to be used for the ancestral states (REF field in VCF).
		-H	Provides the chromosome name (CHROM field in VCF) to overwrite default "chrom" string.
		-E	Converts input FASTA to VCF and terminates execution without further processing.
		
		--- COMMON-OUTLIER ANALYSIS

		-CO	Provides the report name (and column indices for positions and scores) to be used for common-outlier analysis. To perform a common-outlier analysis using RAiSD and SweeD, use -CO like this: "-CO SweeD_Report.SweeD-Run-Name 1 2". The SweeD report must not contain a header. If you have already run RAiSD on your data and only want to perform a common-outlier analysis, use -CO like this: "-CO SweeD_Report.SweeD-Run-Name 1 2 RAiSD_Report.RAiSD-Run-Name 1 X", where X is the index of the column you want to use depending on the RAiSD report. To use the mu-statistic, set Y=2 if RAiSD was invoked with the default parameters, or set Y=7 if -R was explicitly provided or implicitly activated through some other command-line parameter. Again, the RAiSD report must not contain a header.
		-COT	Provides the cut-off threshold for identifying top outliers per report (default: 0.05, i.e., top 5%). 
		-COD	Provides the maximum distance (in number of sites) between outlier points in the provided reports to identify matching outlier positions reported by RAiSD and SweeD. Based on the accuracy of the implemented methods in SweeD and RAiSD, we typically set -COD to a value between 100 and 400 sites (default: 1, i.e., exact match).

		--- HELP and VERSION NOTES

		-h	Prints this help message.
		-v	Prints version information.
		
		
    
    
Input File Formats
----------------------

The current RAiSD release can process SNP data in Hudson's ms or VCF (Variant Call Format) file formats. The d1 test dataset, previously used for the test run, is in Hudson's ms format. For the VCF format, refer to the respective Wikipedia entry (https://en.wikipedia.org/wiki/Variant_Call_Format).

RAiSD version 1.7 (or later) can also parse and process [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) VCF files.  
RAiSD version 1.8 (or later) can also parse compressed VCF files (including ANGSD-generated VCF files) in [GZ](https://www.gnu.org/software/gzip/) file format.                                                                     
RAiSD version 2.1 (or later) can also parse and process unordered VCF files (compressed unordered files are not yet supported).
RAiSD version 2.6 (or later) can also parse and process FASTA files.

Output Files
------------

RAiSD generates two output files, the RAiSD_Info and the RAiSD_Report, with the run name (provided via "-n") as file extension. 

#### RAiSD_Info.test_run

The first 20 lines of the RAiSD_Info.test_run file are shown below.

    RAiSD, Raised Accuracy in Sweep Detection
    Copyright (C) 2017, and GNU GPL'd, by Nikolaos Alachiotis and Pavlos Pavlidis
    Contact n.alachiotis/pavlidisp at gmail.com

    Command: ./RAiSD -n test_run -I d1/msselection1.out -L 100000 -f 
    Samples: 20
    Region:  100000 bp
    Format:  ms

    A pattern structure of 65536 patterns (max. capacity) and approx. 1 MB memory footprint has been created.

    0: Set 0 | sites 6300 | snps 6300 | region 100000 - Var 53090 8.400e-04 | SFS 12005 1.000e+00 | LD 38385 1.333e+00 | MuStat 49905 3.908e-04
    1: Set 1 | sites 6296 | snps 6296 | region 100000 - Var 49665 1.970e-03 | SFS 3745 1.000e+00 | LD 46945 2.000e+00 | MuStat 49685 8.320e-04
    2: Set 2 | sites 6109 | snps 6109 | region 100000 - Var 50835 1.890e-03 | SFS 42100 1.000e+00 | LD 89390 1.500e+00 | MuStat 50770 6.825e-04
    3: Set 3 | sites 6052 | snps 6052 | region 100000 - Var 49590 1.880e-03 | SFS 5545 1.000e+00 | LD 22460 1.500e+00 | MuStat 50150 1.102e-03
    4: Set 4 | sites 6184 | snps 6184 | region 100000 - Var 51920 3.020e-03 | SFS 26445 1.000e+00 | LD 55080 1.333e+00 | MuStat 52400 1.180e-03
    5: Set 5 | sites 5519 | snps 5519 | region 100000 - Var 49310 1.660e-03 | SFS 5475 1.000e+00 | LD 42490 1.333e+00 | MuStat 47165 8.925e-04
    6: Set 6 | sites 6052 | snps 6052 | region 100000 - Var 49640 1.700e-03 | SFS 9090 1.000e+00 | LD 61435 1.250e+00 | MuStat 48080 6.200e-04
    7: Set 7 | sites 6274 | snps 6274 | region 100000 - Var 50480 2.100e-03 | SFS 18910 1.000e+00 | LD 36695 1.333e+00 | MuStat 50490 8.505e-04

The RAiSD_Info file provides execution- and dataset-related information (command line, number of samples, region size, dataset format), as well as a result line per SNP set in the input file. Each per-set result line provides the following information:
a) set index, b) number of sites, c) number of SNPs, d) region size, e) best-score location and the respective score for each of the factors that form the μ statistic, denoted as VAR, SFS, and LD, and f) best-score location and the respective score for the μ statistic (MuStat).

The RAiSD_Info file additionally reports total execution time, memory footprint, and statistics about the total number of SNP sets in the input file (total, processed, skipped), as shown in the last 6 lines of the file.

    Sets (total):     1000
    Sets (processed): 1000
    Sets (skipped):   0

    Total execution time 32.03070 seconds
    Total memory footprint 1109 kbytes
    
#### RAiSD_Report.test_run

The first 20 lines of the RAiSD_Report.test_run file are shown below (using -R to include additional information in the report file).

    // 0
    430	20	840	1.640e-04	2.800e-01	2.580e+00	1.185e-04
    445	20	870	1.700e-04	2.600e-01	2.500e+00	1.105e-04
    450	20	880	1.720e-04	2.600e-01	2.077e+00	9.288e-05
    455	20	890	1.740e-04	2.800e-01	2.210e+00	1.077e-04
    460	20	900	1.760e-04	2.800e-01	2.071e+00	1.020e-04
    480	30	930	1.800e-04	2.800e-01	2.097e+00	1.057e-04
    490	30	950	1.840e-04	2.600e-01	1.667e+00	7.973e-05
    500	40	960	1.840e-04	2.600e-01	1.113e+00	5.324e-05
    515	50	980	1.860e-04	2.800e-01	8.798e-01	4.582e-05
    550	110	990	1.760e-04	2.800e-01	9.087e-01	4.478e-05
    565	120	1010	1.780e-04	2.800e-01	9.688e-01	4.828e-05
    585	150	1020	1.740e-04	2.800e-01	1.308e+00	6.371e-05
    595	170	1020	1.700e-04	2.800e-01	1.687e+00	8.031e-05
    600	170	1030	1.720e-04	2.800e-01	1.595e+00	7.683e-05
    605	180	1030	1.700e-04	3.000e-01	1.352e+00	6.897e-05
    625	190	1060	1.740e-04	3.000e-01	1.356e+00	7.077e-05
    650	240	1060	1.640e-04	3.000e-01	1.080e+00	5.315e-05
    665	260	1070	1.620e-04	3.000e-01	7.723e-01	3.753e-05
    670	270	1070	1.600e-04	2.800e-01	7.634e-01	3.420e-05
    
For each set of SNPs (separated by a line that contains the separator symbol "//" followed by the set index or chromosome name), the RAiSD_Report file contains a set of result lines, one per evaluated genomic location, with each result line providing the genomic location followed by the calculated μ statistic value (default) or the genomic location followed by the start and end positions of the window, the factors VAR, SFS, and LD, and finally the μ statistic (-R command-line flag). 

Required Input Parameters
-------------------------

The in-tool help message generates a list of parameters, as shown above. Those in brackets are optional.
    
    The basic execution mode does NOT require any FREE input parameters. 
   
In addition to the run name ("-n", used to name the output files accordingly) and the path to the input file ("-I") containing SNP data, the region length ("-L") is required only when simulated data are processed, i.e., ms files. The region length, which is required for calculating the μ statistic is extracted from the input file prior to processing when VCF files are analyzed. 

The following command lines provide two examples of command lines for processing ms and VCF files.

Hudson's ms format

    $ ./RAiSD -n ms_run -I input_file.ms -L 1000000

VCF format

    $ ./RAiSD -n vcf_run -I input_file.vcf

Optional Input Parameters
------------------------

In addition to the required set of parameters, several optional ones can be used. The in-tool help message provides a short description for each one of them. 

The -f parameter allows to overwrite existing output files. The default operation mode is to prevent execution in order to avoid accidental overwritting of existing output files.

Parameters -s and -t affect the way output files are generated. The -s parameter splits reports in separate files, one per SNP set, whereas the -t parameter removes the SNP set separator symbols "//".

Processing only a subset of the samples in a VCF file is possible with the use of parameters -p and -S. The former generates a list of all samples in the input file, stored in an output file named RAiSD_Samples.run_name, while the latter allows the user to provide a similar to RAiSD_Samples input file containing only the list of samples to process.

The optional parameters -T, -d, -k, and -l are used for evaluation purposes and are discussed in detail below.

Missing-data Strategies
-----------------------
RAiSD provides four different strategies to handle missing data, using the -M parameter followed by the strategy number.

    0: Discards all SNPs with missing data (default)
    1: Imputes N per SNP   
    2: Creates a mask for valid alleles and treats N as a third state 
    3: Creates a mask for valid alleles and ignores allele pairs with N

Evaluating Accuracy
-------------------
RAiSD provides a series of parameters that facilitate measuring performance when simulated datasets are analyzed. Given a simulated dataset, e.g., in ms format, that comprises several sets of SNPs, the -T parameter can be used in order to direct RAiSD to report accuracy, defined as the average distance between a known sweep location (provided via -T and the reported best-score locations). 

The following command line example for the previously conducted test run demonstrates the use of -T and the respective output. All provided datasets, like the d1 used in the test run, exhibit a selective sweep at the center of the genomic region they simulate. Therefore, RAiSD can be launched as follows:

    $ ./RAiSD -n test_run_accuracy -I d1/msselection1.out -L 100000 -T 50000
    
Note the additional output lines in the RAiSD_Info file, reporting accuracy. The values are in base pairs, averaged over the number of processed sets of SNPs.

    AVERAGE DISTANCE (Target 50000)
    mu-VAR	1022.500
    mu-SFS	26137.655
    mu-LD	18577.810
    MuStat	1253.685

Measuring Success Rate
-----------------------
The -d parameter can be used along with -T to direct RAiSD to report the success rate, defined as the percentage of sets with reported best-score location withing a distance (provided via -d, in base pairs) from the known target of selection (provided via -T, in base pairs). The respective command line for a distance threshold of 1% of the total region length, i.e., 1,000 bp, is shown below.

    $ ./RAiSD -n test_run_success_rate -I d1/msselection1.out -L 100000 -T 50000 -d 1000
    
Note the additional output lines in RAiSD_Info file, reporting success rate. The results show that 62.6% of the runs report locations within 1,000 bp from the known selection target (location 50,000).

    SUCCESS RATE (Distance 1000)
    mu-VAR	0.692
    mu-SFS	0.022
    mu-LD	0.019
    MuStat	0.626

Evaluating Sensitivity (True Positive Rate)
-------------------------------------------

The -k and -l parameters can be used to direct the tool to report sensitivity. This requires to conduct a run on neutral data first and sort all per-set best scores in order to define a threshold for a given False Positive Rate (FPR). RAiSD does this automatically when neutral data are processed, with the use of the -k parameter, providing an FPR. The following command line illustrates this for an FPR of 5%.

    $ ./RAiSD -n test_run_fpr -I d1/msneutral1.out -L 100000 -k 0.05
    
Note the additional lines that report the respective threshold for the chosen FPR value.

    SORTED DATA (FPR 0.050000)
    Size			1000
    Highest Score		0.000450000
    Lowest Score		0.000089600
    FPR Threshold		0.000221867
    Threshold Location	50
    
Thereafter, the FPR Threshold value can be given to RAiSD when selection data are processed (via the -l parameter) as shown below.

    $ ./RAiSD -n test_run_tpr -I d1/msselection1.out -L 100000 -l 0.000221867
    
This will report the respective TPR value as shown in the additional lines in RAiSD_Info file as shown below.

    SCORE COUNT (Threshold 0.000222)
    TPR	    0.995000
    
These results demonstrate a TPR of 99.5% when d1 is analyzed.

Generating μ-statistic plots
----------------------------

The -P parameter can be used to generate a set of four plots in a single PDF file per set of SNPs in the input dataset. This option requires R (https://www.r-project.org/) to be installed and Rscript to be in the default $PATH. By default, -P activates the generation of a separate report per set of SNPs (-s option), the removal of the set separator symbol from each report (-t option), and the inclusion of additional information in each report (-R option). The following command shows the use of -P for the test run.

    $ ./RAiSD -n test_run -I d1/msselection1.out -L 100000 -P
    
This command will generate RAiSD_Plot PDF files with the run name (provided via "-n") followed by the SNP set index (or name in VCF format) preceding the file extension.

RAiSD_Plot.test_run.0.pdf

An example of such a file can be found here: http://139.91.162.50/raisd_plots/examples/RAiSD_Plot.test_run.1.pdf

Processing VCF (ANGSD version)
------------------------------

As of version 1.7, RAiSD can process ANGSD (http://www.popgen.dk/angsd/index.php/ANGSD) VCF files. Such files (see a sample here: http://www.popgen.dk/angsd/index.php/Vcf) only include the GP and GL fields, rather than the GT field that is required by RAiSD. In the absence of GT, diploidy is assumed.  

Given a [GL:GP] entry [L1,L2,L3:P1,P2,P3], RAiSD constructs unphased genotype data on the fly as follows:  

    GT="./." for all-zero GL entries    
    GT="0/0" with probability P1   
    GT="0/1" with probability P2     
    GT="1/1" with probability P3        

Processing VCF in GZ file format
--------------------------------
Since version 1.8, RAiSD can process compressed VCF files in GZ file format. To activate this feature, compile the source code using the dedicated makefile MakefileZLIB as follows:

    $ make clean
    $ make -f MakefileZLIB 

    
The zlib library (https://zlib.net/) needs to be installed prior to compilation. In Ubuntu, the zlib library can be installed with the following command:

    $ sudo apt-get install zlib1g-dev

Processing unordered VCF
------------------------
As of version 2.1, RAiSD can parse and process unordered VCF files (unordered files in GZ file format are ***NOT*** yet supported). Unordered VCF files are derived from [DArTseq](https://www.diversityarrays.com/) genotyping reports. An ordered VCF file to be used in subsequent runs with RAiSD and/or other tools is automatically generated (extension .vcf.fxd). Same-chromosome SNPs are sorted according to their position (POS field in VCF format), and placed in the new VCF file chromosome by chromosome. The chromosome order can be retrieved from the RAiSD_Info file, and should be used to match plots with chromosomes when RAiSD is used to generate plots based on unordered VCF files. 

Generating the RAiSD_SiteReport file
------------------------------------

The -D option can be used to generate the RAiSD_SiteReport, i.e., a single file that includes a line per set of SNPs (ms format) or per CHROM (VCF format), providing a breakdown of the dataset in terms of total number of sites, number of SNPs used in the analysis, and total number of discarded sites. The total number of discarded sites is further broken down into site groups based on the reason they were discarded (failed check).  

When the default missing-data strategy is used (-M 0), the total number of discarded sites includes:  

    a) sites dropped due to a failed "header" check (only for VCF), 
    b) sites dropped because of the MAF check, 
    c) sites with missing data, and 
    d) monomorphic sites. 
    
When any other missing-data strategy is used (-M 1,2, or 3), the total number of discarded sites includes:

    a) sites dropped due to a failed "header" check (only for VCF), 
    b) sites dropped because of the MAF check, 
    c) monomorphic sites, and
    d) potentially monomorphic sites with missing data, i.e., sites with missing data and no variation. 
    
When missing data are imputed per SNP (-M 1), the RAiSD_SiteReport additionally includes the total number of sites that have been imputed before applying the rest of the checks (minor allele frequency, monomorphic sites, etc).

Changing the window size
------------------------

The window size for the sliding-window algorithm can be changed using the -w flag. The default value is set to 50 SNPs, which is empirically determined. Only even numbers can be used, due to a performance optimization employed in the computation of μ_LD. Using very large window sizes is expected to reduce the accuracy of the tool, and also leads to slower execution. Note that, RAiSD does not evaluate varying-size windows per genomic location. Therefore, μ statistic values obtained using different window sizes are not comparable. The following command shows the use of -w to set the window size to 20 SNPs. 


    $ ./RAiSD -n test_run -I d1/msselection1.out -L 100000 -w 20

Adapting the SFS edges for μ_SFS
--------------------------------

The SFS is expected to assume a U shape in the presence of a selective sweep, due to the fact that the number of singletons and SNPs with S-1 derived mutations (where S is the sample size) increases, whereas intermediate SNP classes shrink in size.
By default, the μ_SFS increases with the number of singletons and SNPs with S-1 derived mutations in a window. The -c flag allows to include SNPs with a larger number of derived mutations, e.g., doubletons or tripletons, in the calculation of μ_SFS values. For example, RAiSD can additionally include doubletons with the following example command. 

    $ ./RAiSD -n test_run -I d1/msselection1.out -L 100000 -c 2
    
When a number X that is larger than 1 is provided through -c for the SFS edges, all SNPs with a number of derived alleles less or equal than X and larger or equal than S-X are included in the calculation of μ_SFS values. 

Generating Manhattan plots
----------------------------

As of version 2.0, the -A parameter followed by a probability value (to be used for the quantile function in R) can be used to generate a Manhattan plot in a single PDF file for all chromosomes (sets of SNPs) in the input dataset. This option requires R (https://www.r-project.org/) to be installed and Rscript to be in the default $PATH. By default, -A activates the generation of a separate report per set of SNPs (-s option), the removal of the set separator symbol from each report (-t option), and the inclusion of additional information in each report (-R option). The following command shows the use of -A for the test run.

    $ ./RAiSD -n test_run -I d1/msselection1.out -L 100000 -A 0.995
    
This command will generate a single RAiSD_ManhattanPlot PDF file with the run name (provided via "-n") preceding the file extension.

RAiSD_ManhattanPlot.test_run.pdf

Note that, depending on the number of chromosomes (SNP sets) and the total number of SNPs in the input file, such a RAiSD run might require several minutes only to generate the RAiSD_ManhattanPlot PDF file. 

An example of such a file, obtained by a different run, can be found here: http://139.91.162.50/raisd_plots/examples/RAiSD_ManhattanPlot.test_run_2.pdf

Excluding regions from the analysis
-----------------------------------
As of version 2.3, the -X parameter can be used to provide a path to a tab-delimited text file that contains regions to be excluded from the analysis. The provided file can contain any number of lines, one per region per chromosome to be excluded.
An example file to exclude three regions, two in chr1 and one in chr2, is provided below.
	
	chr1	0	10
	chr1	20	30
	chr2	50	100

It is generally recommended to exclude the centromeres from the analysis, as they can lead to inflated scores and misleading results. 

FASTA-to-VCF conversion
-----------------------
As of version 2.6, RAiSD can also parse and process FASTA files. This is achieved by converting the input FASTA file to VCF format before processing. RAiSD generates the VCF file, which can be used in future runs to skip the conversion step. When processing FASTA files, up to two outgroups can be specified using parameters -C (primary outgroup) and -C2 (secondary outgroup). When no outgroup sequence is provided, the very first sequence in the alignment is used as outgroup. The outgroup sequence(s) must be included in the FASTA file. The -E parameter can be used to terminate execution right after the conversion is completed. The -H parameter can be used to provide a string for the CHROM field (default: "chrom"). 

The steps below describe the conversion process.

a) When no outgroups are provided (or found in the FASTA file), the first ingroup sequence is used. 
   In this case, this sequence is also included in the generated VCF. The outgroup sequences that are specified
   by the user are only used for the REF field and are not included in the sample list of the VCF file.

b) All ambiguous characters are replaced by N.

c) At each position x (VCF line, POS: x), if the primary outgroup state is not informative (gap, N, or 
   different than all ingroup states at position x), the secondary outgroup state is used for the REF field.
   If the secondary outgroup state is not informative as well, the ingroup dominant allele at position x is 
   used for the REF field. 

d) At each position x, the first ingroup state that differs from the REF state is used for the ALT field. If
   no valid ingroup state exists (all ingroup states are gaps, Ns, or the REF state), the ALT field is arbitrarily
   set to a valid DNA state that is not the REF state.

e) Per sample entry, all ALT states in the ingroup are set to 1, all remaining states are set to 0, and all gaps
   and Ns are represented by a dot (‘.’).

Grid size
---------
As of version 2.7, RAiSD can use a grid size parameter, provided via -G, to report mu, muVar, muSFS, and muLD scores. Note that the tool's inherent mu-statistic evaluation algorithm is driven by the number and location of SNPs, with denser SNP regions being evaluated more thoroughly whereas sparser ones being evaluated less exhaustively. It is generally recommended to provide a grid size value that is larger than the total number of SNPs in the dataset to be analyzed. When -G is provided with a grid size of N, RAiSD reports scores at N equidistant locations between the very first and very last sites in the dataset. 


Common outliers between RAiSD and SweeD/OmegaPlus
-------------------------------------------------

As of version 2.8, RAiSD can be used to detect common outliers with SweeD or OmegaPlus. To perform a common-outlier analysis between RAiSD and SweeD, for example, you can follow these steps:

	1) Execute SweeD to create report SweeD_Report.RUNNAME
	2) Remove the header line (if exists) from the generated SweeD_Report.RUNNAME file
	3) Execute RAiSD with the same grid size that was used for SweeD, and provide the SweeD_Report as input via -CO.
	
To use SweeD, refer to its manual, which is available at: https://github.com/alachins/sweed/blob/master/sweed3.0_manual.pdf.
RAiSD (step 3) can be used in two different modes: A) process the input dataset and then perform the common-outlier analysis, or B) only perform the common-outlier analysis (using a RAiSD_Report from a previous run). 

To use RAiSD in mode A, only the SweeD_Report needs to be provided via -CO as follows:

	RAiSD -n myrunname -I dataset.vcf -G SweeDgridSize -CO SweeD_Report.RUNNAME 1 2
	
In this mode, the -CO command-line flag receives 3 arguments, i.e., a path to the SweeD_Report file followed by the column indices for the positions and the scores, which, for SweeD reports, are 1 and 2, respectively.

When the RAiSD_Report already exists, generated by a previous RAiSD run with the same grid size used with SweeD (SweeDgridSize), RAiSD can be used in mode B as follows:

	RAiSD -n myrunname -CO SweeD_Report.RUNNAME 1 2 RAiSD_Report.RUNNAME2 1 2
	
In this case, the -CO command-line flag receives 3 additional arguments, i.e., a path to the RAiSD_Report file followed by the column indices for the positions and the scores. 

The optional -COT parameter allows to set the threshold for identifying top outliers per run. The default value is 0.05, which detects the top 5% scores as outliers in each of the two reports.

The optional -COD parameter allows some flexibility in matching common-outlier positions between the two runs by setting a maximum distance among RAiSD/SweeD common outliers. When the distance between two outlier positions, one in RAiSD_Report and one in SweeD_Report, is less or equal than the value provided via -COD, these points are identified as a common-outlier pair. The default COD value is 1 site, which requires an exact match between outlier positions. This is unlikely to occur given that RAiSD and SweeD employ different routines for computing the evaluated positions. Therefore, a larger number can be used, e.g., 100 sites.

The output of a common-outlier analysis with RAiSD is a set of four files:

	a) RAiSD_CommonOutlierReport.myrunname
	b) RAiSD_CommonOutlierPlot.myrunname.pdf
	c) RAiSD_CommonOutlierPointsRAiSD.myrunname
	d) RAiSD_CommonOutlierPointsSweeD.myrunname

The CommonOutlierReport contains all position pairs, with the first column referring to a position in SweeD_Report while the second column refers to a position in RAiSD_Report. The CommonOutlierPlot illustrates the results of both tools and the common-outlier positions (red points). The CommonOutlierPoints files contain outlier positions and their respective scores per tool/method. 

An example of a common-outlier analysis based on the reports obtained by the analyses of 1,601 SARS-CoV-2 genomes (available at the [GISAID database](https://platform.gisaid.org/)) using RAiSD and SweeD is provided in the common-outliers-SARS-CoV-2 directory. The command line used to perform the common-outlier analysis is the following:

	../RAiSD -n common-outliers -CO SweeD_Report.SARS-CoV-2 1 2 RAiSD_Report.SARS-CoV-2 1 2 -COT 0.05 -COD 400




Changelog
----------

	v1.0 (9/6/2017): first release
	v1.1 (7/3/2018): MAF threshold option	
	v1.2 (28/3/2018): mbs format with -b
	v1.3 (18/7/2018): -i to impute N per SNP, -a for rand seed
	v1.4 (3/8/2018): -M to handle missing data, -O to show progress on the display device
	v1.5 (4/8/2018): -R to include additional information in the report (reduced default to location and score)
	v1.6 (3/9/2018): -P to create plots per set of SNPs using Rscript
	v1.7 (2/10/2018): -y for ploidy, -D for site report, fixed a bug in the plotting routine
	v1.8 (31/12/2018): MakefileZLIB to parse VCF files in gzip file format (requires the zlib library)
	v1.9 (27/4/2019): -w to set the window size (default 50), -c to set the SFS slack for the μ_SFS
	v2.0 (15/5/2019): -A to create Manhattan plots, scale factors for muVar and muSFS to yield comparable scores among different chromosomes
	v2.1 (21/1/2020): Parser for unordered VCF files (e.g., DArTseq genotyping reports). Creates the ordered VCF file as well.
	v2.2 (22/1/2020): Added missing field (discarded monomorphic sites) in the site report for M={1,2,3}.
	v2.3 (23/1/2020): -X to exclude regions per chromosome from the analysis.
	v2.4 (30/1/2020): -B for chromosome length and SNP size. Fixed bug with the memory-reduction optimization for large chromosomes. -o to request vcf ordering and generation.
	v2.5 (8/2/2020): Fixed position bug due to typecasting. Some site positions were off by 1 bp.
	v2.6 (2/4/2020): Parses, converts to vcf, and analyzes fasta input files (-C/-C2 for outgroups, -H for chromosome name, -E for conversion-only mode).
	v2.7 (8/4/2020): -G parameter to specify the grid size
	v2.8 (22/4/2020): -CO, -COT, -COD parameters for common-outlier analysis between RAiSD and SweeD, install script for gsl
	v2.9 (6/8/2020): fixed bug in parsing one-character VCF sample names

Support
-------

We provide support for RAiSD through our generic [Google group for Population Genetics](https://groups.google.com/forum/#!forum/popgen-support).
    
    




