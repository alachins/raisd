
RAiSD: Software for selective sweep detection
===============================================

Authors: Nikolaos Alachiotis (n.alachiotis@gmail.com), Pavlos Pavlidis (ppavlidis@gmail.com)

Date: 9/6/2017

Version: 1.0

About
-----

RAiSD (Raised Accuracy in Sweep Detection) is a stand-alone software implementation of the μ statistic for selective sweep detection. Unlike existing implementations, including our previously released tools (SweeD and OmegaPlus), RAiSD scans whole-genome SNP data based on a composite evaluation scheme that captures multiple sweep signatures at once. 

Download and Compile
--------------------

The following commands can be used to download and compile the source code. 

    $ mkdir RAiSD
    $ cd RAiSD
    $ wget https://github.com/alachins/raisd/archive/master.zip
    $ unzip master.zip
    $ cd raisd-master
    $ make
    
The executable is placed in the path RAiSD/raisd-master/bin/release. A link to the executable is placed in the installation folder, i.e., raisd-master.

Test Run
--------

To check that RAiSD is installed correctly, a test run can be done with the following commands. These commands are going to download a dataset (one of the many that we used for evaluation purposes) and execute RAiSD.
    
    $ wget 139.91.162.50/raisd_data/d70.tar.gz
    $ tar -xvzf d70.tar.gz
    $ ./RAiSD -n test_run -I d70/msselection70.out
    
    
    


    
In-tool Help
------------
    $ RAiSD -h
    
Supported File Formats
----------------------
RAiSD can process SNP data in Hudson's ms or VCF (Variant Call Format) files.

Command-line Input Arguments
------------------------
The basic execution mode DOES NOT require any FREE input parameters. 









