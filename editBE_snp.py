#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 09:30:45 2021

@author: em924
"""

""" search and find the editing window which is 20 nt upstream of PAM sequence (3nt) """
""" the script expects a pre-proccessed gwas table  which is the output of  the 'gwas.ref.maj.py' script"""
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description='proccessing gwas sumstats to find SNPs with editing poteintial.')
parser.add_argument('-g', '--gwas', help="input gwas summary statistics with atleaset these columns SNP,CHR,BP,A1,A2 (defult is a tab delimiter)",required=True)
parser.add_argument('-r', '--ref' ,help="reference genome",required=True)
parser.add_argument('-o', '--output' ,help="output file name",required=True)

args = parser.parse_args()
print(args)

#query=pd.read_csv(args.gwas.txt,sep= '\t')
query=pd.read_csv(args.gwas,sep=None,engine="python")
query2=query[((query['A1']=='C') & (query['A2']=='T')) | (query['A1']=='T') & (query['A2']=='C')]
query3=query[((query['A1']=='A') & (query['A2']=='G')) | (query['A1']=='G') & (query['A2']=='A')]
query4=query2.append(query3)
#query4.to_csv("acceptable_alleles_gwas.txt", sep='\t')

with open(args.ref,'r') as f:
        seq=f.read()

seq2=seq.split('>')[1:]

seq3={}

for i in seq2:
    i=i.strip().split('\n')
    seq3[i[0]]=''.join(i[1:])


def gc_content(gRNA):
    count=0
    gc=0
    for nt in gRNA:
        if nt=="G" or nt=="g" or nt=="C" or nt=="c" :
            gc+=1
        count+=1
    d=(gc*100)/count
    return d

def position(in_seq):
    pos=1
    dd=[]
    for i in in_seq:
        if i=='c' or i=='C':
            dd.append("c"+str(pos+3))
        pos+=1
    dd= ",".join(repr(e) for e in dd )
    return dd

def position_A(in_seq):
    pos=1
    dd=[]
    for i in in_seq:
        if i=='a' or i=='A':
            dd.append("a"+str(pos+3))
        pos+=1
    dd= ",".join(repr(e) for e in dd )
    return dd

def Rev_compl(st):
    complement= {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', 'N': 'N'}
    return "".join(complement[n] for n in reversed(st))

dic=seq3
edit_strand=[]
edit_site=[]
count_c=[]
pos_c=[]
gRNA=[]
gc=[]
feature=[]

""" A1,A2=Alt|Ref , MAJOR=major_allele """
for key,value in dic.items():
    for SNP,CHR,BP,A1,A2,MAJOR,REF in zip(query4['SNP'],query4['CHR'],query4['BP'],query4['A2'],query4['A1'],query4['major'],query4['ref']):

        #changing "C" to "T" on the POSITIVE strand of DNA to make a C->T conversion on + strand.
        if MAJOR=='C' or MAJOR=='not_found' and REF=='C':
               if key[3:] ==str(CHR):
                   pam_seq=value[BP+14:BP+19]
                   for match in re.finditer(r"(?=(gg))", pam_seq,re.IGNORECASE):
                       if match.start()==4:
                           site=value[BP:BP+5] if MAJOR==REF else 'C'+value[BP+1:BP+5]
                           gRNA_seq=value[BP-3:BP+17] if MAJOR==REF else value[BP-3:BP-1]+'C'+value[BP:BP+17]
                       elif match.start()==3:
                           site=value[BP-1:BP+4] if MAJOR==REF else 'C'+value[BP:BP+4]
                           gRNA_seq=value[BP-4:BP+16] if MAJOR==REF else value[BP-4:BP-1]+'C'+value[BP:BP+16]
                       elif match.start()==2:
                           site=value[BP-2:BP+3] if MAJOR==REF else value[BP-2:BP-1]+'C'+value[BP:BP+3]
                           gRNA_seq=value[BP-5:BP+15] if MAJOR==REF else value[BP-5:BP-1]+'C'+value[BP:BP+15]
                       elif match.start()==1:
                           site=value[BP-3:BP+2] if MAJOR==REF else value[BP-3:BP-1]+'C'+value[BP:BP+2]
                           gRNA_seq=value[BP-6:BP+14] if MAJOR==REF else value[BP-6:BP-1]+'C'+value[BP:BP+14]
                       elif match.start()==0:
                           site=value[BP-4:BP+1] if MAJOR==REF else value[BP-4:BP-1]+'C'+value[BP:BP+1]
                           gRNA_seq=value[BP-7:BP+13] if MAJOR==REF else value[BP-7:BP-1]+'C'+value[BP:BP+13]
                       number_c=site.count('C') + site.count('c')
                       if number_c>0:
                           pos_c.append(position(site))
                           gc.append(gc_content(gRNA_seq))
                           edit_strand.append('+')
                           edit_site.append(site)
                           count_c.append(number_c)
                           gRNA.append(gRNA_seq)
                           snp='\t'.join([SNP,str(CHR),str(BP),A1,A2,MAJOR,REF])
                           feature.append(snp)

        #changing "A" to "G" on the POSITIVE strand of DNA to make a A->G conversion on + strand.
        if MAJOR=='A' or MAJOR=='not_found' and REF=='A':
               if key[3:] ==str(CHR):
                   pam_seq=value[BP+14:BP+19]
                   for match in re.finditer(r"(?=(gg))", pam_seq,re.IGNORECASE):
                       if match.start()==4:
                           site=value[BP:BP+5] if MAJOR==REF else 'A'+value[BP+1:BP+5]
                           gRNA_seq=value[BP-3:BP+17] if MAJOR==REF else value[BP-3:BP-1]+'A'+value[BP:BP+17]
                       elif match.start()==3:
                           site=value[BP-1:BP+4] if MAJOR==REF else 'A'+value[BP:BP+4]
                           gRNA_seq=value[BP-4:BP+16] if MAJOR==REF else value[BP-4:BP-1]+'A'+value[BP:BP+16]
                       elif match.start()==2:
                           site=value[BP-2:BP+3] if MAJOR==REF else value[BP-2:BP-1]+'A'+value[BP:BP+3]
                           gRNA_seq=value[BP-5:BP+15] if MAJOR==REF else value[BP-5:BP-1]+'A'+value[BP:BP+15]
                       elif match.start()==1:
                           site=value[BP-3:BP+2] if MAJOR==REF else value[BP-3:BP-1]+'A'+value[BP:BP+2]
                           gRNA_seq=value[BP-6:BP+14] if MAJOR==REF else value[BP-6:BP-1]+'A'+value[BP:BP+14]
                       elif match.start()==0:
                           site=value[BP-4:BP+1] if MAJOR==REF else value[BP-4:BP-1]+'A'+value[BP:BP+1]
                           gRNA_seq=value[BP-7:BP+13] if MAJOR==REF else value[BP-7:BP-1]+'A'+value[BP:BP+13]
                       number_c=site.count('A') + site.count('a')
                       if number_c>0:
                           pos_c.append(position_A(site))
                           gc.append(gc_content(gRNA_seq))
                           edit_strand.append('+')
                           edit_site.append(site)
                           count_c.append(number_c)
                           gRNA.append(gRNA_seq)
                           snp='\t'.join([SNP,str(CHR),str(BP),A1,A2,MAJOR,REF])
                           feature.append(snp)

        #changing "A" to "G" on the NEGATIVE strand of DNA to make a T->C conversion on + strand.
        if MAJOR=='T' or MAJOR=='not_found' and REF=='T':
               if key[3:] ==str(CHR):
                   pam_seq=value[BP-20:BP-14]
                   for match in re.finditer(r"(?=(cc))", pam_seq,re.IGNORECASE):

                       if match.start()==0:
                           site=Rev_compl(value[BP-5:BP]) if MAJOR==REF else 'A'+Rev_compl(value[BP-5:BP-1])
                           gRNA_seq=Rev_compl(value[BP-17:BP+3]) if MAJOR==REF else Rev_compl(value[BP:BP+3])+'A'+Rev_compl(value[BP-17:BP-1])
                       elif match.start()==1:
                           site=Rev_compl(value[BP-4:BP+1]) if MAJOR==REF else Rev_compl(value[BP:BP+1])+'A'+Rev_compl(value[BP-4:BP-1])
                           gRNA_seq=Rev_compl(value[BP-16:BP+4]) if MAJOR==REF else Rev_compl(value[BP:BP+4])+'A'+Rev_compl(value[BP-16:BP-1])
                       elif match.start()==2:
                           site=Rev_compl(value[BP-3:BP+2]) if MAJOR==REF else Rev_compl(value[BP:BP+2])+'A'+Rev_compl(value[BP-3:BP-1])
                           gRNA_seq=Rev_compl(value[BP-15:BP+5]) if MAJOR==REF else Rev_compl(value[BP:BP+5])+'A'+Rev_compl(value[BP-15:BP-1])
                       elif match.start()==3:
                           site=Rev_compl(value[BP-2:BP+3]) if MAJOR==REF else Rev_compl(value[BP:BP+3])+'A'+Rev_compl(value[BP-2:BP-1])
                           gRNA_seq=Rev_compl(value[BP-14:BP+6]) if MAJOR==REF else Rev_compl(value[BP:BP+6])+'A'+Rev_compl(value[BP-14:BP-1])
                       elif match.start()==4:
                           site=Rev_compl(value[BP-1:BP+4]) if MAJOR==REF else Rev_compl(value[BP:BP+4])+'A'
                           gRNA_seq=Rev_compl(value[BP-13:BP+7]) if MAJOR==REF else Rev_compl(value[BP:BP+7])+'A'+Rev_compl(value[BP-13:BP-1])
                       number_c=site.count('A') + site.count('a')
                       if number_c>0:
                           pos_c.append(position_A(site))
                           gc.append(gc_content(gRNA_seq))
                           edit_strand.append('-')
                           edit_site.append(site)
                           count_c.append(number_c)
                           gRNA.append(gRNA_seq)
                           snp='\t'.join([SNP,str(CHR),str(BP),A1,A2,MAJOR,REF])
                           feature.append(snp)

        #changing "C" to "T" on the NEGATIVE strand of DNA to make a G->A conversion on + strand.
        if MAJOR=='G' or MAJOR=='not_found' and REF=='G':
               if key[3:] ==str(CHR):
                   pam_seq=value[BP-20:BP-14]
                   for match in re.finditer(r"(?=(cc))", pam_seq,re.IGNORECASE):
                       if match.start()==0:
                           site=Rev_compl(value[BP-5:BP]) if MAJOR==REF else 'C'+Rev_compl(value[BP-5:BP-1])
                           gRNA_seq=Rev_compl(value[BP-17:BP+3]) if MAJOR==REF else Rev_compl(value[BP:BP+3])+'C'+Rev_compl(value[BP-17:BP-1])
                       elif match.start()==1:
                           site=Rev_compl(value[BP-4:BP+1]) if MAJOR==REF else Rev_compl(value[BP:BP+1])+'C'+Rev_compl(value[BP-4:BP-1])
                           gRNA_seq=Rev_compl(value[BP-16:BP+4]) if MAJOR==REF else Rev_compl(value[BP:BP+4])+'C'+Rev_compl(value[BP-16:BP-1])
                       elif match.start()==2:
                           site=Rev_compl(value[BP-3:BP+2]) if MAJOR==REF else Rev_compl(value[BP:BP+2])+'C'+Rev_compl(value[BP-3:BP-1])
                           gRNA_seq=Rev_compl(value[BP-15:BP+5]) if MAJOR==REF else Rev_compl(value[BP:BP+5])+'C'+Rev_compl(value[BP-15:BP-1])
                       elif match.start()==3:
                           site=Rev_compl(value[BP-2:BP+3]) if MAJOR==REF else Rev_compl(value[BP:BP+3])+'C'+Rev_compl(value[BP-2:BP-1])
                           gRNA_seq=Rev_compl(value[BP-14:BP+6]) if MAJOR==REF else Rev_compl(value[BP:BP+6])+'C'+Rev_compl(value[BP-14:BP-1])
                       elif match.start()==4:
                           site=Rev_compl(value[BP-1:BP+4]) if MAJOR==REF else Rev_compl(value[BP:BP+4])+'C'
                           gRNA_seq=Rev_compl(value[BP-13:BP+7]) if MAJOR==REF else Rev_compl(value[BP:BP+7])+'C'+Rev_compl(value[BP-13:BP-1])
                       number_c=site.count('C') + site.count('c')
                       if number_c>0:
                           pos_c.append(position(site))
                           gc.append(gc_content(gRNA_seq))
                           edit_strand.append('-')
                           edit_site.append(site)
                           count_c.append(number_c)
                           gRNA.append(gRNA_seq)
                           snp='\t'.join([SNP,str(CHR),str(BP),A1,A2,MAJOR,REF])
                           feature.append(snp)

al=[None]*len(edit_site)
ind=0
for a,b,c,d,e,f,g in zip (feature,edit_strand,edit_site,count_c,pos_c,gRNA,gc):
        al[ind]='\t'.join([a,b,c,str(d),str(e),f,str(g)])
        ind+=1

""" prioritization of the editing sites"""
al2=[]
for item in al:
    item=item.split("\t")
    if (item[9]=='1') and ('c8' in item[10] or 'a8' in item[10]):
        item.append("A-5")
    elif (item[9]=='1') and ('c7' in item[10] or 'a7' in item[10]):
        item.append("A-3")
    elif (item[9]=='1') and ('c6' in item[10] or 'a6' in item[10]):
        item.append("A-2")
    elif (item[9]=='1') and ('c5' in item[10] or 'a5' in item[10]):
        item.append("A-1")
    elif (item[9]=='1') and ('c4' in item[10] or 'a4' in item[10]):
        item.append("A-4")

    else:
        if item[9]=='2':
            item.append("B")
        elif item[9]=='3':
            item.append("C")
        elif item[9]=='4':
            item.append("D")
        elif item[9]=='5':
            item.append("E")

    al2.append(item) 

for item in al2:
    item[10]=item[10].replace("'","")

""" write the results"""   
res=args.output 
with open(res,'w') as fe:
    fe.write("SNP\tchr\tBP\tA1\tA2\tmajor\tref\tedit_strand\tedit_site\t#ofC/A\tC/A_position\tgRNA\tgRNA_gc_content\teditSite_priority\n")
    for line in al2:
        for i in line:
            fe.write(i+"\t")
        fe.write('\n')
