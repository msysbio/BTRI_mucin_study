#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:21:06 2022
This script downloads the KEGG and COGG annotation for BT via API. For RI, the annotation has been compiled from EGGNOGG database.
@author: u0145079
"""
####################################################
######Downloading COG and KEGG annotation###########
####################################################

import requests
import json
import pandas as pd
import itertools 
import os
import numpy as np

#Bacteroides_thetaiotaomicron

#URL for the API:ie. https://www.ncbi.nlm.nih.gov/research/cog/api/cog/?organism=Nitrosopumilus_maritimus_SCM1
#NCBI-COG database

url = "https://www.ncbi.nlm.nih.gov/research/cog/api/cog/?organism=Bacteroides_thetaiotaomicron_VPI-5482"

# Grab the search results
response = requests.get(url)
resp_dict = response.json()

#start with the empty list
df=pd.DataFrame()

# Store the first page of results

counter = 0
for i in range(len(resp_dict["results"])): 
    df.loc[counter,"cogs"] = resp_dict["results"][i]["cog"]["cogid"] 
    df.loc[counter, "ids"] = resp_dict["results"][i]["gene_tag"]
    df.loc[counter, "description"] = resp_dict["results"][i]["cog"]["funcats"][0].get('description')
    counter+=1

# While data['next'] isn't empty,vdownload the next page:
while resp_dict['next'] is not None:
    print("Next page found, downloading", resp_dict['next'])
    response = requests.get(resp_dict['next'])
    resp_dict = response.json()
    # Store the current page of results
    for i in range(len(resp_dict["results"])):
        df.loc[counter,"cogs"] = resp_dict["results"][i]["cog"]["cogid"] 
        df.loc[counter, "ids"] = resp_dict["results"][i]["gene_tag"]
        df.loc[counter, "description"] = resp_dict["results"][i]["cog"]["funcats"][0].get('description')
        counter+=1
       
##KEGG KO's extraction
kegg = pd.read_table("/Users/u0145079/Desktop/PlasticDaphnia/OLD/Scripts/Bash/bth_ko.txt",delimiter="\t",header=None)
ec_numbers = pd.read_table("/Users/u0145079/Desktop/PlasticDaphnia/OLD/Scripts/Bash/bth_ec_numbers.txt",delimiter="\t",header=None)

##EC KO's extraction
kegg.rename(columns = {0:'gene', 1:'ko'}, inplace = True)
ec_numbers.rename(columns = {0:'gene', 1:'ec'}, inplace = True)

#extract gene_ids
file1 = open('/Users/u0145079/Desktop/Miscellanious/Bin/Ensemble_indexes/BT/Bacteroides_thetaiotaomicron_vpi_5482_gca_000011065.ASM1106v1.49.gff3', 'r')
Lines = file1.readlines()

#Extract gene_ids and 
gene_ids = []
gene_desp = []
for i in range(len(Lines)):
    if "gene_id" in Lines[i]:   
        id=Lines[i].split("gene_id")[1].split(";")[0]
        desc=Lines[i].split("description")[1].split(";")[0]
        if id not in gene_ids:
            gene_ids.append(id)
            gene_desp.append(desc)

entries = list(zip(gene_ids,gene_desp))

#match the gene_ids with the lists of KEGG, COGG and EC_annotation
bt_annot = pd.DataFrame(entries, columns = ['gene_id','description'])

###make gene_ids in the tables identical before merge
#remove quotation marks in bt_annot table
for i in range(len(bt_annot["gene_id"])):
    bt_annot["gene_id"][i] = bt_annot["gene_id"][i].replace('=','')
    bt_annot["description"][i] = bt_annot["description"][i].replace('=','')
    #bt_annot["gene_id"][i] = bt_annot["gene_id"][i].replace('"','')
    #bt_annot["gene_id"][i] = bt_annot["gene_id"][i].replace(' ','')

#remove bth: and ko: prefixes from kegg table
for i in range(len(kegg["gene"])):
    kegg["gene"][i] = kegg["gene"][i].replace('bth:','')
    kegg["ko"][i] = kegg["ko"][i].replace("ko:","")

#remove bth: and ko: prefixes from kegg table
for i in range(len(ec_numbers["gene"])):
    ec_numbers["gene"][i] = ec_numbers["gene"][i].replace('bth:','')
    ec_numbers["ec"][i] = ec_numbers["ec"][i].replace("ec:","")

#merge all tables into final annotation table
##ko numbers
bt_annot["ko"] = "NA"
for i in range(len(bt_annot["gene_id"])):
    for j in range(len(kegg["gene"])):
        if kegg["gene"][j] == bt_annot["gene_id"][i]:
            bt_annot["ko"][i] = kegg["ko"][j]
##ec numbers
bt_annot["EC_number"] = "NA"
for i in range(len(bt_annot["gene_id"])):
    for j in range(len(ec_numbers["gene"])):
        if ec_numbers["gene"][j] == bt_annot["gene_id"][i]:
            bt_annot["EC_number"][i] = ec_numbers["ec"][j]

#COGs
bt_annot["cogs"] = "NA"
bt_annot["cogs_description"] = "NA"
for i in range(len(bt_annot["gene_id"])):
    for j in range(len(df["ids"])):
        if df["ids"][j] == bt_annot["gene_id"][i]:
            bt_annot["cogs"][i] = df["cogs"][j]
            bt_annot["cogs_description"][i] = df["description"][j]

bt_annot.to_csv('/Users/u0145079/Desktop/Miscellanious/Bin/Bacteroides_thetaiotaomicron/Bacteroides_thetaiotaomicron/Bacteroides_thetaiotaomicron_annotation.csv')


#Roseburia_intestinalis
#concatenate files from the same batch

folder = '/Users/u0145079/Desktop/Miscellanious/Bin/Roseburia_Intestalis/Roseburia_intestinalis/BTRI_batch_mucin'
samples = sorted(os.listdir(folder))
samples.remove('.DS_Store')
data = {}

for samp in samples:
    with open(os.path.join(folder, samp, 'quant.sf' )) as f:
        f.readline()
        for line in f:
            a = line.strip().split('\t')
            if a[0] not in data:
                data[a[0]] = [float(a[4])]
            else:
                data[a[0]].append(float(a[4]))


df = pd.DataFrame(data.items())
split_df = pd.DataFrame(df[1].tolist(), columns= samples)
split_df.index = data.keys()

#delete extra words
split_df["original_tags"] = split_df.index
split_df.index = split_df.index.to_series().str.split('|').str[2]

##Upload PATRIC ffn file:
folder = '/Users/u0145079/Desktop/Miscellanious/Bin/Roseburia_Intestalis/Roseburia_intestinalis/536231.75.PATRIC.ffn'

tags = list()

with open(folder) as f:
    for line in f:
        if line.startswith(">fig"):
            if line not in tags:
                tags.append(line)

for item in tags:
    for i in range(len(split_df["original_tags"])):
        if item.split(" ")[0].replace(">","") == split_df["original_tags"][i]:
            split_df["original_tags"][i] = item.replace(">","")
    
split_df.to_csv('/Users/u0145079/Desktop/Miscellanious/Bin/Roseburia_Intestalis/Roseburia_intestinalis/BTRI_batch_mucin.csv')
            

####Annotation
#extract gene_ids
file1 = open('/Users/u0145079/Desktop/Miscellanious/Bin/Roseburia_Intestalis/Roseburia_intestinalis/genome_files/536231.75.PATRIC.gff', 'r')
Lines = file1.readlines()

#Extract gene_ids and gene description
gene_ids = []
gene_desp = []
for i in range(len(Lines)):
    if "product" in Lines[i]:
        desc=Lines[i].split("product")[1].split(";")[0]
    if "locus_tag" in Lines[i]:
        id=Lines[i].split("locus_tag")[1].split(";")[0]
        if id not in gene_ids:
            gene_ids.append(id)
            gene_desp.append(desc)

entries = list(zip(gene_ids,gene_desp))

#make into dataframe
ri_annot = pd.DataFrame(entries, columns = ['gene_ids','description'])

#remove some extra characters
ri_annot= ri_annot.replace('=', '', regex=True)

#load egg_nog .csv annotation table
eggnogg = pd.read_table("/Users/u0145079/Desktop/Miscellanious/Bin/Roseburia_Intestalis/Roseburia_intestinalis/MM_uvp4lx8d.emapper.annotations.tsv", sep="\t")

for i in range(len(eggnogg["#query"])):
    eggnogg["eggNOG_OGs"][i]= eggnogg["eggNOG_OGs"][i].split("@")[0]
    try:
        if eggnogg["#query"][i].split("|")[2]:
            eggnogg["#query"][i] = eggnogg["#query"][i].split("|")[2]
    except:
        pass

#merge all tables into final annotation table
##ko numbers
ri_annot["ko"] = "NA"
ri_annot["EC_number"] = "NA"
ri_annot["cogs"] = "NA"
ri_annot["cogs_description"] = "NA"
ri_annot["cogs_category"] = "NA"

for i in range(len(ri_annot["gene_ids"])):
    for j in range(len(eggnogg["#query"])):
        if ri_annot["gene_ids"][i] == eggnogg["#query"][j]:
            ri_annot["ko"][i] = eggnogg["KEGG_ko"][j]
            ri_annot["EC_number"][i] = eggnogg["EC"][j]
            ri_annot["cogs"][i] = eggnogg["eggNOG_OGs"][j]
            ri_annot["cogs_description"][i] = eggnogg["Description"][j].replace(r'\n','') 
            ri_annot["cogs_category"][i] = eggnogg["COG_category"][j]
            
ri_annot.to_csv('/Users/u0145079/Desktop/Miscellanious/Bin/Roseburia_Intestalis/Roseburia_intestinalis/Roseburia_intestinalis_annotation_EGGNOGG.csv')





