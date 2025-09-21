# gwas-BES

Scripts for identifying candidate SNPs suitable for base editing CRISPR system in GWAS data.

This tool was developed to detect GWAS SNPs across the genome that can be edited by the Base Editing CRISPR system.  

This was run on PYTHON 3.7.6 and it needs the following libraries to be installed:  
-Pandas 
-regular expression (re) 
-argparse

The algorithm includes two steps; 
1) determine the reference and major allele for each SNP; 
2) predict editable SNPs

Detailes:

STEP 1:run command: 

      python gwas.ref.maj.py --gwas <gwas sumstats> --ref <reference genpome assembly> --output <outfile name>

 
INPUT files:

--gwas; GWAS summary stats (GWAS sumstats) that  must have the column headers CHR SNP BP A1 A2. Its does not matter what 
delimiter is as the script recognises the delimiter.

NOTE: The gwas sumstats should be first filtered for only significat hits (P-value < 5*10^-8).  

--ref; Reference genome that  must match the GWAS sumstats (e.g. if the GWAS is hg19 then the assembly hg19 
should be provided)

NOTE: This script uses "g1000_eur.bim" so this file must be in the running folder.
NOTE: This script can also be used by its self if one wishes to get information about the SNPs (major/minor 
alleles and Ref/Alt alleles)

OUTPUT file:   

- The output will be a GWAS sumstats with three new columns including 'ref', 'major', and 'minor' appended to it.


STEP 2: run command: 

      python editBE_snp.py --gwas <gwas sumstats> --ref <reference genpome assembly> --output <outfile name>


INPUT files:

--gwas; A preproccessed gwas table which is the output of the setep 1.

--ref; Reference genome that must match the GWAS sumstats (see above).

OUTPUT file:

The output file will be the final result which is a table with the following headers:

SNP ; The gwas SNP 

chr : Chromosome number of the SNP	

BP ; Base pair location of SNP

A1 ; A1 column in the gwas	

A2 ; A1 column in the gwas	

major; Majore alle of the SNP	

ref ; Reference allele of the SNP

edit_strand ; The DNA strand that editing will occure 

edit_site ; a 5-nt site that editing take palce. This is the editing acticity window

#ofC/A	; Number of "C" or "A" nucleotides present in the activity window  

C/A_position ; The position of the C/A in the activity window	

gRNA ; Sequence of the protospacer. This is a 20nt long sequence spaning from PAM site to activity window

gRNA_gc_content	; peoportion of 'GC' in the gRAN sequence

editSite_priority ; efficiency-specificity of the editing of the target SNP, which is based on the number of C/A 
and their position in the activity window. This factor is classified into 5 groups, with A1 represents the higheset 
efficiency-specificity and D the lowest as follow:

A-1 ; there is 1 C/A in the site and the position is nucleotide 5 (5th nt in the protospacer)

A-2 ; there is 1 C/A in the site and the position is nucleotide 6 (6th nt in the protospacer)

A-3 ; there is 1 C/A in the site and the position is nucleotide 7 (7th nt in the protospacer)

A-4 ; there is 1 C/A in the site and the position is nucleotide 4 (4th nt in the protospacer)

A-5 ; there is 1 C/A in the site and the position is nucleotide 8 (8th nt in the protospacer)

B   ; there is 2 C/A in the site

C   ; there is 3 C/A in the site

D   ; there is 4 C/A in the site



END

