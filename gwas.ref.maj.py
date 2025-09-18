#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Oct 14 09:30:45 2021

@author: em924
"""
""" search and find the Ref allele in GWAS summary statistics """
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description='Processing gwas sumstats to find SNPs with editing poteintial.')
parser.add_argument('-g', '--gwas', help="input gwas summary statistics with atleaset these columns SNP,CHR,BP,A1,A2,",required=True)
parser.add_argument('-r', '--ref' ,help="reference genome, the reference needs to match the gwas ref",required=True)
parser.add_argument('-o', '--output' ,help="output file name",required=True)

args = parser.parse_args()
print(args)

query=pd.read_csv(args.gwas,sep=None,engine="python")

with open(args.ref,'r') as f:
        seq=f.read()

seq2=seq.split('>')[1:]

seq3={}

for i in seq2:
    i=i.strip().split('\n')
    seq3[i[0]]=''.join(i[1:])

dic=seq3

query2=query
query2['ref']=""
n=len(query.columns)
print(n)
i=0
for CHR,SNP,BP,A1,A2 in zip(query['CHR'],query['SNP'],query['BP'],query['A1'],query['A2']):
    for key,seq in dic.items():
        if key.replace("chr","")==str(CHR):
            if seq[BP-1].upper()==str(A1):
                query2.iloc[i,n-1]=A1
            elif seq[BP-1].upper()==str(A2):
                query2.iloc[i,n-1]=A2
            else:
                query2.iloc[i,n-1]=A1
            i+=1
#query.to_csv('gwas_ref.txt',sep='\t',index=False)

#query=pd.read_csv('gwas_ref.txt',sep=None,engine="python")
g1k=pd.read_csv('g1000_eur.bim',sep=None,engine="python")

g1k2=pd.merge(g1k,query2['SNP'],on='SNP', how='inner')

query2['major']=""
query2['minor']=""
m=len(query2.columns)

i=0
for SNP in query2['SNP']:
    for SNP2,MAJOR,MINOR in zip(g1k2['SNP'],g1k2['MAJOR'],g1k2['MINOR']):
        if SNP==SNP2:
            query2.iloc[i,m-2]=MAJOR
            query2.iloc[i,m-1]=MINOR
            break
        else:
            query2.iloc[i,m-2]="not_found"
            query2.iloc[i,m-1]="not_found"
    i+=1

query2.to_csv(args.output,sep='\t',index=False)
