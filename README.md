
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
    $ mkdir RAiSD
    $ cd RAiSD
    $ wget https://github.com/alachins/raisd/archive/master.zip
    $ unzip master.zip
    $ cd raisd-master
    $ make
    
In-tool Help
------------
    $ RAiSD -h
    
Supported File Formats
----------------------
RAiSD can process SNP data in Hudson's ms or VCF (Variant Call Format) files.

Command-line Input Arguments
------------------------
The basic execution mode DOES NOT require any FREE input parameters. 









