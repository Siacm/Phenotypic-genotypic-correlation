#!/usr/bin/python

###imported packages, uses python 3.0
import sys
import re
import getopt
import pandas as pd
import numpy as np
import csv
from io import StringIO
from collections import Counter, defaultdict

##################different functions#####################
def usage():
    print ("<merge_phenogeno.py> -i infile | -o outfile\n [Options: -h help | -a ST count | -b ST_gpattern count | -c gpattern_ST count | -d ST_allgenes count | -e allgenes_ST count | -f ST_gene count | -g gene_ST count")


def make_simpledict(dict, dictkey, dictvalue, df, index):
    for i in index:
        dict.setdefault(df[dictkey].iloc[i],[]).append(df[dictvalue].iloc[i])
    return dict
    
def dict_singlegene(index, dict, dictkey, str, df):
    for i in index:
        for col in df:
            if df[col].iloc[i] == str: #if column equals yes or maybe
                dict.setdefault(df[dictkey].iloc[i],[]).append(col)
    return dict

def dict_gpattern(index, dict, dictkey, str, df):
    for i in index:
        gpattern_list=[]
        for col in df:
            if df[col].iloc[i] == str: #if column equals yes or maybe
                gpattern_list.append(col)
                dict.setdefault(df[dictkey].iloc[i],[]).append(gpattern_list)
    return dict

def reverse_dict_simple(old_dict, new_dict):
    for olddict_key, olddict_value in old_dict.items():
        for x in olddict_value:
            new_dict.setdefault(x,[]).append(olddict_key)
    return new_dict

def reverse_dict(old_dict, new_dict):
    for olddict_key, olddict_value in old_dict.items():
        for x in olddict_value:
            new_dict.setdefault(repr(x),[]).append(olddict_key) #repr to make the list a key.
    return new_dict

def dict_loop(d):
    for key, value in d.items():
        outfile.write(key+"\t"+ str(len([item for item in value if item]))+"\n")
        
def counting_values(d):
    total_count=defaultdict(list)
    for key, value in d.items():
        for each_v in value:
            x=(repr(each_v))
            if x in total_count:
                total_count[x]+=1
            else:
                total_count[x]=1
        for k, v in total_count.items():
            outfile.write (key+"\t"+k+"\t"+str(v)+"\n")

def read_infile():
    #goes through each row to determine whether a column contains S,R and/or yes and assigns phenogeno value under column name "Concordance"
    infile_csvfile.loc[infile_csvfile.eq("S").any(1) & ~infile_csvfile.eq("yes").any(1)  & ~infile_csvfile.eq("maybe").any(1),"Concordance"] = "phenoSgenoS"
    infile_csvfile.loc[infile_csvfile.eq("S").any(1) & ~infile_csvfile.eq("yes").any(1) & infile_csvfile.eq("maybe").any(1) ,"Concordance"] = "phenoSgenoM"
    infile_csvfile.loc[infile_csvfile.eq("S").any(1) & infile_csvfile.eq("yes").any(1),"Concordance"] = "phenoSgenoR"
    infile_csvfile.loc[infile_csvfile.eq("R").any(1) & infile_csvfile.eq("yes").any(1),"Concordance"] = "phenoRgenoR"
    infile_csvfile.loc[infile_csvfile.eq("R").any(1) & ~infile_csvfile.eq("yes").any(1) & infile_csvfile.eq("maybe").any(1) ,"Concordance"] = "phenoRgenoM"
    infile_csvfile.loc[infile_csvfile.eq("R").any(1) & ~infile_csvfile.eq("yes").any(1) & ~infile_csvfile.eq("maybe").any(1) ,"Concordance"] = "phenoRgenoS"
    infile_csvfile.loc[infile_csvfile.eq("I").any(1) & infile_csvfile.eq("yes").any(1),"Concordance"] = "phenoIgenoR"
    infile_csvfile.loc[infile_csvfile.eq("I").any(1) & ~infile_csvfile.eq("yes").any(1) & infile_csvfile.eq("maybe").any(1) ,"Concordance"] = "phenoIgenoM"
    infile_csvfile.loc[infile_csvfile.eq("I").any(1) & ~infile_csvfile.eq("yes").any(1) & ~infile_csvfile.eq("maybe").any(1),"Concordance"] = "phenoIgenoS"

    #gets the row index of isolates with corresponding pheno geno pattern
    phenoSgenoS_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoSgenoS"].index
    phenoSgenoM_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoSgenoM"].index
    phenoSgenoR_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoSgenoR"].index
    phenoRgenoR_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoRgenoR"].index
    phenoRgenoM_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoRgenoM"].index
    phenoRgenoS_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoRgenoS"].index
    phenoIgenoR_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoIgenoR"].index
    phenoIgenoM_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoIgenoM"].index
    phenoIgenoS_index = infile_csvfile[infile_csvfile['Concordance'] == "phenoIgenoS"].index
    return(infile_csvfile,phenoSgenoS_index,phenoSgenoM_index,phenoSgenoR_index,phenoRgenoR_index,phenoRgenoM_index,phenoRgenoS_index,phenoIgenoR_index,phenoIgenoM_index,phenoIgenoS_index)

##########################function for overall counts
def main_counts():

    ID_ST_pSgS_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoSgenoS_index}
    ID_ST_pSgM_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoSgenoM_index}
    ID_ST_pSgR_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoSgenoR_index}
    ID_ST_pRgR_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoRgenoR_index}
    ID_ST_pRgM_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoRgenoM_index}
    ID_ST_pRgS_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoRgenoS_index}
    ID_ST_pIgR_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoIgenoR_index}
    ID_ST_pIgM_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoIgenoM_index}
    ID_ST_pIgS_dict={infile_csvfile['ID'].iloc[i]:infile_csvfile['Serovar_ST'].iloc[i] for i in phenoIgenoS_index}
    
    outfile.write("\n\nTotal Counts\n")
    outfile.write("There are "+str(len(ID_ST_pSgS_dict))+" isolates that are phenotypically and genotypically sensitive\n")
    outfile.write("There are "+str(len(ID_ST_pSgM_dict))+" isolates that are phenotypically sensitive and maybe genotypically resistant\n")
    outfile.write("There are "+str(len(ID_ST_pSgR_dict))+" isolates that are phenotypically sensitive but genotypically resistant\n")
    outfile.write("There are "+str(len(ID_ST_pRgR_dict))+" isolates that are phenotypically and genotypically resistant\n")
    outfile.write("There are "+str(len(ID_ST_pRgM_dict))+" isolates that are phenotypically resistant and maybe genotypically resistant\n")
    outfile.write("There are "+str(len(ID_ST_pRgS_dict))+" isolates that are phenotypically resistant but genotypically sensitive\n")
    outfile.write("There are "+str(len(ID_ST_pIgR_dict))+" isolates that are phenotypically intermediate and genotypically resistant\n")
    outfile.write("There are "+str(len(ID_ST_pIgM_dict))+" isolates that are phenotypically intermediate and maybe genotypically resistant\n")
    outfile.write("There are "+str(len(ID_ST_pIgS_dict))+" isolates that are phenotypically intermediate and genotypically sensitive\n\n\n")

    return (ID_ST_pSgS_dict, ID_ST_pSgM_dict, ID_ST_pSgR_dict, ID_ST_pRgR_dict, ID_ST_pRgM_dict, ID_ST_pRgS_dict, ID_ST_pIgR_dict, ID_ST_pIgM_dict, ID_ST_pIgS_dict)

##########################dictionary for Serovar and lab ID's
def ST_count():
    ST_ID_pSgS_dict=defaultdict(list)
    ST_ID_pSgM_dict=defaultdict(list)
    ST_ID_pSgR_dict=defaultdict(list)
    ST_ID_pRgR_dict=defaultdict(list)
    ST_ID_pRgM_dict=defaultdict(list)
    ST_ID_pRgS_dict=defaultdict(list)
    ST_ID_pIgR_dict=defaultdict(list)
    ST_ID_pIgM_dict=defaultdict(list)
    ST_ID_pIgS_dict=defaultdict(list)

    make_simpledict(ST_ID_pSgS_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoSgenoS_index)
    make_simpledict(ST_ID_pSgM_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoSgenoM_index)
    make_simpledict(ST_ID_pSgR_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoSgenoR_index)
    make_simpledict(ST_ID_pRgR_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoRgenoR_index)
    make_simpledict(ST_ID_pRgM_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoRgenoM_index)
    make_simpledict(ST_ID_pRgS_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoRgenoS_index)
    make_simpledict(ST_ID_pIgR_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoIgenoR_index)
    make_simpledict(ST_ID_pIgM_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoIgenoM_index)
    make_simpledict(ST_ID_pIgS_dict, 'Serovar_ST', 'ID', infile_csvfile, phenoIgenoS_index)

    return(ST_ID_pSgS_dict,ST_ID_pSgM_dict,ST_ID_pSgR_dict,ST_ID_pRgR_dict,ST_ID_pRgM_dict,ST_ID_pRgS_dict,ST_ID_pIgR_dict,ST_ID_pIgM_dict,ST_ID_pIgS_dict)

##########################dictionary for Serovar and gene_pattern

def  ST_gpattern_count():
    ST_gpattern_pSgM_dict=defaultdict(list)
    ST_gpattern_pSgR_dict=defaultdict(list)
    ST_gpattern_pRgR_dict=defaultdict(list)
    ST_gpattern_pRgM_dict=defaultdict(list)
    ST_gpattern_pIgR_dict=defaultdict(list)
    ST_gpattern_pIgM_dict=defaultdict(list)

    dict_gpattern(phenoSgenoM_index,ST_gpattern_pSgM_dict,"Serovar_ST", "maybe", infile_csvfile)
    dict_gpattern(phenoSgenoR_index,ST_gpattern_pSgR_dict,"Serovar_ST", "yes", infile_csvfile)
    dict_gpattern(phenoRgenoR_index,ST_gpattern_pRgR_dict,"Serovar_ST", "yes", infile_csvfile)
    dict_gpattern(phenoRgenoM_index,ST_gpattern_pRgM_dict,"Serovar_ST", "maybe", infile_csvfile)
    dict_gpattern(phenoIgenoR_index,ST_gpattern_pIgR_dict,"Serovar_ST", "yes", infile_csvfile)
    dict_gpattern(phenoIgenoM_index,ST_gpattern_pIgM_dict,"Serovar_ST", "maybe", infile_csvfile)
    
    return ST_gpattern_pSgM_dict,ST_gpattern_pSgR_dict,ST_gpattern_pRgR_dict,ST_gpattern_pRgM_dict,ST_gpattern_pIgR_dict,ST_gpattern_pIgM_dict

##########################dictionary for gene_pattern and serovars

def gpattern_ST_count():
    gpattern_ST_pSgM_dict=defaultdict(list)
    gpattern_ST_pSgR_dict=defaultdict(list)
    gpattern_ST_pRgR_dict=defaultdict(list)
    gpattern_ST_pRgM_dict=defaultdict(list)
    gpattern_ST_pIgR_dict=defaultdict(list)
    gpattern_ST_pIgM_dict=defaultdict(list)

    reverse_dict(ST_gpattern_pSgM_dict,gpattern_ST_pSgM_dict)
    reverse_dict(ST_gpattern_pSgR_dict,gpattern_ST_pSgR_dict)
    reverse_dict(ST_gpattern_pRgR_dict,gpattern_ST_pRgR_dict)
    reverse_dict(ST_gpattern_pRgM_dict,gpattern_ST_pRgM_dict)
    reverse_dict(ST_gpattern_pIgR_dict,gpattern_ST_pIgR_dict)
    reverse_dict(ST_gpattern_pIgM_dict,gpattern_ST_pIgM_dict)

    return(gpattern_ST_pSgM_dict,gpattern_ST_pSgR_dict,gpattern_ST_pRgR_dict,gpattern_ST_pRgM_dict,gpattern_ST_pIgR_dict,gpattern_ST_pIgM_dict)

##########################dictionary for Serovar and all genes

def ST_allgenes_count():
    ST_allgenes_pSgS_dict=defaultdict(list)
    ST_allgenes_pSgM_dict=defaultdict(list)
    ST_allgenes_pSgR_dict=defaultdict(list)
    ST_allgenes_pRgR_dict=defaultdict(list)
    ST_allgenes_pRgM_dict=defaultdict(list)
    ST_allgenes_pRgS_dict=defaultdict(list)
    ST_allgenes_pIgR_dict=defaultdict(list)
    ST_allgenes_pIgM_dict=defaultdict(list)
    ST_allgenes_pIgS_dict=defaultdict(list)

    make_simpledict(ST_allgenes_pSgS_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoSgenoS_index)
    make_simpledict(ST_allgenes_pSgM_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoSgenoM_index)
    make_simpledict(ST_allgenes_pSgR_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoSgenoR_index)
    make_simpledict(ST_allgenes_pRgR_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoRgenoR_index)
    make_simpledict(ST_allgenes_pRgM_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoRgenoM_index)
    make_simpledict(ST_allgenes_pRgS_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoRgenoS_index)
    make_simpledict(ST_allgenes_pIgR_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoIgenoR_index)
    make_simpledict(ST_allgenes_pIgM_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoIgenoM_index)
    make_simpledict(ST_allgenes_pIgS_dict, 'Serovar_ST', 'all_genes', infile_csvfile, phenoIgenoS_index)
    return(ST_allgenes_pSgS_dict,ST_allgenes_pSgM_dict,ST_allgenes_pSgR_dict,ST_allgenes_pRgR_dict,ST_allgenes_pRgM_dict,ST_allgenes_pRgS_dict,ST_allgenes_pIgR_dict,ST_allgenes_pIgM_dict,ST_allgenes_pIgS_dict)


##########################dictionary for allgene and serovars
def allgene_ST_count():
    allgenes_ST_pSgS_dict=defaultdict(list)
    allgenes_ST_pSgM_dict=defaultdict(list)
    allgenes_ST_pSgR_dict=defaultdict(list)
    allgenes_ST_pRgR_dict=defaultdict(list)
    allgenes_ST_pRgM_dict=defaultdict(list)
    allgenes_ST_pRgS_dict=defaultdict(list)
    allgenes_ST_pIgR_dict=defaultdict(list)
    allgenes_ST_pIgM_dict=defaultdict(list)
    allgenes_ST_pIgS_dict=defaultdict(list)

    reverse_dict_simple(ST_allgenes_pSgS_dict,allgenes_ST_pSgS_dict)
    reverse_dict_simple(ST_allgenes_pSgM_dict,allgenes_ST_pSgM_dict)
    reverse_dict_simple(ST_allgenes_pSgR_dict,allgenes_ST_pSgR_dict)
    reverse_dict_simple(ST_allgenes_pRgR_dict,allgenes_ST_pRgR_dict)
    reverse_dict_simple(ST_allgenes_pRgM_dict,allgenes_ST_pRgM_dict)
    reverse_dict_simple(ST_allgenes_pRgS_dict,allgenes_ST_pRgS_dict)
    reverse_dict_simple(ST_allgenes_pIgR_dict,allgenes_ST_pIgR_dict)
    reverse_dict_simple(ST_allgenes_pIgM_dict,allgenes_ST_pIgM_dict)
    reverse_dict_simple(ST_allgenes_pIgS_dict,allgenes_ST_pIgS_dict)
    return(allgenes_ST_pSgS_dict,allgenes_ST_pSgM_dict,allgenes_ST_pSgR_dict,allgenes_ST_pRgR_dict,allgenes_ST_pRgM_dict,allgenes_ST_pRgS_dict,allgenes_ST_pIgR_dict,allgenes_ST_pIgM_dict,allgenes_ST_pIgS_dict)

##########################dictionary for Serovar and single genes
def ST_gene_count():
    ST_gene_pSgM_dict=defaultdict(list)
    ST_gene_pSgR_dict=defaultdict(list)
    ST_gene_pRgR_dict=defaultdict(list)
    ST_gene_pRgM_dict=defaultdict(list)
    ST_gene_pIgR_dict=defaultdict(list)
    ST_gene_pIgM_dict=defaultdict(list)

    dict_singlegene(phenoSgenoM_index,ST_gene_pSgM_dict,"Serovar_ST", "maybe", infile_csvfile)
    dict_singlegene(phenoSgenoR_index,ST_gene_pSgR_dict,"Serovar_ST", "yes", infile_csvfile)
    dict_singlegene(phenoRgenoR_index,ST_gene_pRgR_dict,"Serovar_ST", "yes", infile_csvfile)
    dict_singlegene(phenoRgenoM_index,ST_gene_pRgM_dict,"Serovar_ST", "maybe", infile_csvfile)
    dict_singlegene(phenoIgenoR_index,ST_gene_pIgR_dict,"Serovar_ST", "yes", infile_csvfile)
    dict_singlegene(phenoIgenoM_index,ST_gene_pIgM_dict,"Serovar_ST", "maybe", infile_csvfile)

    return(ST_gene_pSgM_dict,ST_gene_pSgR_dict,ST_gene_pRgR_dict,ST_gene_pRgM_dict,ST_gene_pIgR_dict,ST_gene_pIgM_dict)

##########################dictionary for single gene and serovar
def gene_ST_count():
    gene_ST_pSgM_dict=defaultdict(list)
    gene_ST_pSgR_dict=defaultdict(list)
    gene_ST_pRgR_dict=defaultdict(list)
    gene_ST_pRgM_dict=defaultdict(list)
    gene_ST_pIgR_dict=defaultdict(list)
    gene_ST_pIgM_dict=defaultdict(list)

    reverse_dict_simple(ST_gene_pSgM_dict,gene_ST_pSgM_dict)
    reverse_dict_simple(ST_gene_pSgR_dict,gene_ST_pSgR_dict)
    reverse_dict_simple(ST_gene_pRgR_dict,gene_ST_pRgR_dict)
    reverse_dict_simple(ST_gene_pRgM_dict,gene_ST_pRgM_dict)
    reverse_dict_simple(ST_gene_pIgR_dict,gene_ST_pIgR_dict)
    reverse_dict_simple(ST_gene_pIgM_dict,gene_ST_pIgM_dict)

    return(gene_ST_pSgM_dict,gene_ST_pSgR_dict,gene_ST_pRgR_dict,gene_ST_pRgM_dict,gene_ST_pIgR_dict,gene_ST_pIgM_dict)

def options():
    if call_a:
        read_infile()
        ST_count()

        #print serovar count
        outfile.write("********************** Pheno: S Geno: S **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pSgS_dict)
        outfile.write("********************** Pheno: S Geno: M **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pSgM_dict)
        outfile.write("********************** Pheno: S Geno: R **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pSgR_dict)
        outfile.write("********************** Pheno: R Geno: R **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pRgR_dict)
        outfile.write("********************** Pheno: R Geno: M **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pRgM_dict)
        outfile.write("********************** Pheno: R Geno: S **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pRgS_dict)
        outfile.write("********************** Pheno: I Geno: R **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pIgR_dict)
        outfile.write("********************** Pheno: I Geno: M **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pIgM_dict)
        outfile.write("********************** Pheno: I Geno: S **********************\n")
        outfile.write("ST\tCount\n")
        dict_loop(ST_ID_pIgS_dict)
        outfile.write("********************** ********************** ********************** \n")


    if call_b:
        read_infile()
        ST_gpattern_count()
        
        #print serovar count
        outfile.write("********************** Pheno: S Geno: M **********************\n")
        outfile.write("ST\tgene_pattern\tCount\n")
        counting_values(ST_gpattern_pSgM_dict)
        outfile.write("********************** Pheno: S Geno: R **********************\n")
        outfile.write("ST\tgene_pattern\tCount\n")
        counting_values(ST_gpattern_pSgR_dict)
        outfile.write("********************** Pheno: R Geno: R **********************\n")
        outfile.write("ST\tgene_pattern\tCount\n")
        counting_values(ST_gpattern_pRgR_dict)
        outfile.write("********************** Pheno: R Geno: M **********************\n")
        outfile.write("ST\tgene_pattern\tCount\n")
        counting_values(ST_gpattern_pRgM_dict)
        outfile.write("********************** Pheno: I Geno: R **********************\n")
        outfile.write("ST\tgene_pattern\tCount\n")
        counting_values(ST_gpattern_pIgR_dict)
        outfile.write("********************** Pheno: I Geno: M **********************\n")
        outfile.write("ST\tgene_pattern\tCount\n")
        counting_values(ST_gpattern_pIgM_dict)
        outfile.write("********************** ********************** ********************** \n")

        
    if call_c:
        read_infile()
        gpattern_ST_count()
        
        #outfile.write serovar count
        outfile.write("********************** Pheno: S Geno: M **********************\n")
        outfile.write("gene_pattern\tST\tCount\n")
        counting_values(gpattern_ST_pSgM_dict)
        outfile.write("********************** Pheno: S Geno: R **********************\n")
        outfile.write("gene_pattern\tST\tCount\n")
        counting_values(gpattern_ST_pSgR_dict)
        outfile.write("********************** Pheno: R Geno: R **********************\n")
        outfile.write("gene_pattern\tST\tCount\n")
        counting_values(gpattern_ST_pRgR_dict)
        outfile.write("********************** Pheno: R Geno: M **********************\n")
        outfile.write("gene_pattern\tST\tCount\n")
        counting_values(gpattern_ST_pRgM_dict)
        outfile.write("********************** Pheno: I Geno: R **********************\n")
        outfile.write("gene_pattern\tST\tCount\n")
        counting_values(gpattern_ST_pIgR_dict)
        outfile.write("********************** Pheno: I Geno: M **********************\n")
        outfile.write("gene_pattern\tST\tCount\n")
        counting_values(gpattern_ST_pIgM_dict)
        outfile.write("********************** ********************** **********************\n")


    if call_d:
        read_infile()
        ST_allgenes_count()

        #print  count
        outfile.write("********************** Pheno: S Geno: S **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pSgS_dict)
        outfile.write("********************** Pheno: S Geno: M **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pSgM_dict)
        outfile.write("********************** Pheno: S Geno: R **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pSgR_dict)
        outfile.write("********************** Pheno: R Geno: R **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pRgR_dict)
        outfile.write("********************** Pheno: R Geno: M **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pRgM_dict)
        outfile.write("********************** Pheno: R Geno: S **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pRgS_dict)
        outfile.write("********************** Pheno: I Geno: R **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pIgR_dict)
        outfile.write("********************** Pheno: I Geno: M **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pIgM_dict)
        outfile.write("********************** Pheno: I Geno: S **********************\n")
        outfile.write("ST\tallgenes\tCount\n")
        counting_values(ST_allgenes_pIgS_dict)
        outfile.write("********************** ********************** ********************** \n")


    if call_e:
        read_infile()
        allgenes_ST_count()

        #print  count
        outfile.write("********************** Pheno: S Geno: S **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pSgS_dict)
        outfile.write("********************** Pheno: S Geno: M **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pSgM_dict)
        outfile.write("********************** Pheno: S Geno: R **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pSgR_dict)
        outfile.write("********************** Pheno: R Geno: R **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pRgR_dict)
        outfile.write("********************** Pheno: R Geno: M **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pRgM_dict)
        outfile.write("********************** Pheno: R Geno: S **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pRgS_dict)
        outfile.write("********************** Pheno: I Geno: R **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pIgR_dict)
        outfile.write("********************** Pheno: I Geno: M **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pIgM_dict)
        outfile.write("********************** Pheno: I Geno: S **********************\n")
        outfile.write("allgenes\tST\tCount\n")
        counting_values(allgenes_ST_pIgS_dict)
        outfile.write("********************** ********************** ********************** \n")


    if call_f:
        read_infile()
        ST_gene_count()
        
        outfile.write("********************** Pheno: S Geno: M **********************\n")
        outfile.write("ST\tgene\tCount\n")
        counting_values(ST_gene_pSgM_dict)
        outfile.write("********************** Pheno: S Geno: R **********************\n")
        outfile.write("ST\tgene\tCount\n")
        counting_values(ST_gene_pSgR_dict)
        outfile.write("********************** Pheno: R Geno: R **********************\n")
        outfile.write("ST\tgene\tCount\n")
        counting_values(ST_gene_pRgR_dict)
        outfile.write("********************** Pheno: R Geno: M **********************\n")
        outfile.write("ST\tgene\tCount\n")
        counting_values(ST_gene_pRgM_dict)
        outfile.write("********************** Pheno: I Geno: R **********************\n")
        outfile.write("ST\tgene\tCount\n")
        counting_values(ST_gene_pIgR_dict)
        outfile.write("********************** Pheno: I Geno: M **********************\n")
        outfile.write("ST\tgene\tCount\n")
        counting_values(ST_gene_pIgM_dict)
        outfile.write("********************** ********************** ********************** \n")

        
    if call_g:
        read_infile()
        gene_ST_count()
        
        outfile.write("********************** Pheno: S Geno: M **********************\n")
        outfile.write("Gene\tST\tCount\n")
        counting_values(gene_ST_pSgM_dict)
        outfile.write("********************** Pheno: S Geno: R **********************\n")
        outfile.write("Gene\tST\tCount\n")
        counting_values(gene_ST_pSgR_dict)
        outfile.write("********************** Pheno: R Geno: R **********************\n")
        outfile.write("Gene\tST\tCount\n")
        counting_values(gene_ST_pRgR_dict)
        outfile.write("********************** Pheno: R Geno: M **********************\n")
        outfile.write("Gene\tST\tCount\n")
        counting_values(gene_ST_pRgM_dict)
        outfile.write("********************** Pheno: I Geno: R **********************\n")
        outfile.write("Gene\tST\tCount\n")
        counting_values(gene_ST_pIgR_dict)
        outfile.write("********************** Pheno: I Geno: M **********************\n")
        outfile.write("Gene\tST\tCount\n")
        counting_values(gene_ST_pIgM_dict)
        outfile.write("********************** ********************** ********************** \n")





#################################################################################
###################################MAIN SCRIPT###################################
#################################################################################

if __name__ == "__main__":
    if not len(sys.argv[1:]):
        usage()
    try:
        opts,args = getopt.getopt(sys.argv[1:],"hi:o:abcdefg",["help", "infile=", "output="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()

    infile = None
    infile_present=False
    outfile = None
    outfile_present=False
    call_a = False
    call_b = False
    call_c = False
    call_d = False
    call_e = False
    call_f = False
    call_g = False
    
    for opt,arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(2)
        elif opt in ("-i","--infile"):
            infile_present=True
            infile = open(arg, 'r')
        elif opt in ("-o","--output"):
            outfile_present=True
            outfile = open(arg, 'w')
        elif opt in ("-a"):
            call_a=True
        elif opt in ("-b"):
            call_b=True
        elif opt in ("-c"):
            call_c=True
        elif opt in ("-d"):
            call_d=True
        elif opt in ("-e"):
            call_e=True
        elif opt in ("-f"):
            call_f=True
        elif opt in ("-g"):
            call_g=True
        else:
            assert False,"Unhandled Option"

    if not infile_present:
        print ("-i was not given")
        usage()
        sys.exit(2)
    elif not outfile_present:
        print ("-o was not given")
        usage()
        sys.exit(2)
        
    if infile_present:
        #display all rows and columns in pandas dataframe
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)

        # read csv file "infile"
        infile_csvfile = pd.read_csv(infile)

        #combine columns Serovar and ST and insert it after column index
        infile_csvfile.insert(infile_csvfile.columns.get_loc('ST'),'Serovar_ST',infile_csvfile['serovar']+" ST"+ infile_csvfile["ST"])

        #drop columns Serovar and ST inplace of current dataframe infile_csvfile
        infile_csvfile.drop(['serovar','ST'], axis=1, inplace=True)

        #counts the number of NaNs and sorts the dataframe in order of the ascending order of most number of null results - this will help avoid removing the wrong duplicates
        infile_csvfile = infile_csvfile.iloc[infile_csvfile.isnull().sum(1).sort_values(ascending=True).index]

        #remove duplicate rows based on column ID, keeping the last entry
        infile_csvfile.drop_duplicates("ID",keep="first", inplace=True)

        #remove indexes that were dropped previously
        infile_csvfile.reset_index(drop=True,inplace=True)
    
    
    #assigning variables
    infile_csvfile,phenoSgenoS_index,phenoSgenoM_index,phenoSgenoR_index,phenoRgenoR_index,phenoRgenoM_index,phenoRgenoS_index,phenoIgenoR_index,phenoIgenoM_index,phenoIgenoS_index=read_infile()
    ID_ST_pSgS_dict, ID_ST_pSgM_dict, ID_ST_pSgR_dict, ID_ST_pRgR_dict, ID_ST_pRgM_dict, ID_ST_pRgS_dict, ID_ST_pIgR_dict, ID_ST_pIgM_dict, ID_ST_pIgS_dict=main_counts()
    ST_ID_pSgS_dict,ST_ID_pSgM_dict,ST_ID_pSgR_dict,ST_ID_pRgR_dict,ST_ID_pRgM_dict,ST_ID_pRgS_dict,ST_ID_pIgR_dict,ST_ID_pIgM_dict,ST_ID_pIgS_dict=ST_count()
    ST_gpattern_pSgM_dict,ST_gpattern_pSgR_dict,ST_gpattern_pRgR_dict,ST_gpattern_pRgM_dict,ST_gpattern_pIgR_dict,ST_gpattern_pIgM_dict=ST_gpattern_count()
    gpattern_ST_pSgM_dict,gpattern_ST_pSgR_dict,gpattern_ST_pRgR_dict,gpattern_ST_pRgM_dict,gpattern_ST_pIgR_dict,gpattern_ST_pIgM_dict=gpattern_ST_count()
    ST_allgenes_pSgS_dict,ST_allgenes_pSgM_dict,ST_allgenes_pSgR_dict,ST_allgenes_pRgR_dict,ST_allgenes_pRgM_dict,ST_allgenes_pRgS_dict,ST_allgenes_pIgR_dict,ST_allgenes_pIgM_dict,ST_allgenes_pIgS_dict=ST_allgenes_count()
    allgenes_ST_pSgS_dict,allgenes_ST_pSgM_dict,allgenes_ST_pSgR_dict,allgenes_ST_pRgR_dict,allgenes_ST_pRgM_dict,allgenes_ST_pRgS_dict,allgenes_ST_pIgR_dict,allgenes_ST_pIgM_dict,allgenes_ST_pIgS_dict=allgene_ST_count()
    ST_gene_pSgM_dict,ST_gene_pSgR_dict,ST_gene_pRgR_dict,ST_gene_pRgM_dict,ST_gene_pIgR_dict,ST_gene_pIgM_dict=ST_gene_count()
    gene_ST_pSgM_dict,gene_ST_pSgR_dict,gene_ST_pRgR_dict,gene_ST_pRgM_dict,gene_ST_pIgR_dict,gene_ST_pIgM_dict=gene_ST_count()


    options()


infile.close()
outfile.close()
#################################################################################
###############################END OF SCRIPT#####################################
#################################################################################
















##############################################################################
##########################EXTRAS (not used)##########################
##############################################################################

##########################dictionary for lab ID and gene_pattern
ID_gpattern_pSgM_dict=defaultdict(list)
ID_gpattern_pSgR_dict=defaultdict(list)
ID_gpattern_pRgR_dict=defaultdict(list)
ID_gpattern_pRgM_dict=defaultdict(list)
ID_gpattern_pIgR_dict=defaultdict(list)
ID_gpattern_pIgM_dict=defaultdict(list)

dict_gpattern(phenoSgenoM_index,ID_gpattern_pSgM_dict,"ID", "maybe", infile_csvfile)
dict_gpattern(phenoSgenoR_index,ID_gpattern_pSgR_dict,"ID", "yes", infile_csvfile)
dict_gpattern(phenoRgenoR_index,ID_gpattern_pRgR_dict,"ID", "yes", infile_csvfile)
dict_gpattern(phenoRgenoM_index,ID_gpattern_pRgM_dict,"ID", "maybe", infile_csvfile)
dict_gpattern(phenoIgenoR_index,ID_gpattern_pIgR_dict,"ID", "yes", infile_csvfile)
dict_gpattern(phenoIgenoM_index,ID_gpattern_pIgM_dict,"ID", "maybe", infile_csvfile)

##########################dictionary for gene_pattern and lab IDs
gpattern_ID_pSgM_dict=defaultdict(list)
gpattern_ID_pSgR_dict=defaultdict(list)
gpattern_ID_pRgR_dict=defaultdict(list)
gpattern_ID_pRgM_dict=defaultdict(list)
gpattern_ID_pIgR_dict=defaultdict(list)
gpattern_ID_pIgM_dict=defaultdict(list)

reverse_dict(ID_gpattern_pSgM_dict,gpattern_ID_pSgM_dict)
reverse_dict(ID_gpattern_pSgR_dict,gpattern_ID_pSgR_dict)
reverse_dict(ID_gpattern_pRgR_dict,gpattern_ID_pRgR_dict)
reverse_dict(ID_gpattern_pRgM_dict,gpattern_ID_pRgM_dict)
reverse_dict(ID_gpattern_pIgR_dict,gpattern_ID_pIgR_dict)
reverse_dict(ID_gpattern_pIgM_dict,gpattern_ID_pIgM_dict)
