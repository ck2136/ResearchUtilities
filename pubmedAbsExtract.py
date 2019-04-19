#!/usr/bin/env python
# - - - - - - - - - - - - - - - - - - - - - # 
# Filename      : C.K.
# Purpose       : Script for Downloading Abstracts from PubMed
#                 for the Systematic Review
# Date created  :    
# Created by    :    
# Last modified : Fri 19 Apr 2019 02:05:40 AM MDT
# Modified by   : ck
# - - - - - - - - - - - - - - - - - - - - - # 


# - - - - - - - - - - - - - - - - - - - - - # 
# Import libraries 
# - - - - - - - - - - - - - - - - - - - - - # 
from Bio import Entrez
import time
import sys
import os
import pandas as pd
import xml.etree.ElementTree as ET


try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2
# - - - - - - - - - - - - - - - - - - - - - # 
# Print Python Version
# - - - - - - - - - - - - - - - - - - - - - # 
print(sys.version)

# - - - - - - - - - - - - - - - - - - - - - # 
# Set path to save files
# - - - - - - - - - - - - - - - - - - - - - # 
path = "/home/ck/Downloads/abstracts/"

# - - - - - - - - - - - - - - - - - - - - - # 
# Register to NCBI
# - - - - - - - - - - - - - - - - - - - - - # 
Entrez.email = "chong.kim@ucdenver.edu" # LET ENTREZ KNOW

# - - - - - - - - - - - - - - - - - - - - - # 
# Search pubmed for abstracts
# - - - - - - - - - - - - - - - - - - - - - # 

search_results = Entrez.read(Entrez.esearch(db="pubmed",
    term="Asthma AND (Cost-Effectiveness OR Cost-Benefit OR Economic-Evaluation)",
    reldate=100, # Set relative date (reldate) to however long by days
    datetype="pdat", 
    usehistory="y"))

count = int(search_results["Count"])
print("Found %i results" % count)

# - - - - - - - - - - - - - - - - - - - - - # 
# Download Abstract (using XML format)
# - - - - - - - - - - - - - - - - - - - - - # 

batch_size = 10
batch = 1
for start in range(0,count,batch_size):
    end = min(count, start+batch_size)
    print("Going to download record %i to %i" % (start+1, end))
    attempt = 1
    while attempt <= 3: 
        try: 
            fetch_handle = Entrez.efetch(db="pubmed", 
                    #rettype="abstract", # can change to 'abstract'
                    retmode="xml",retstart=start,
                    retmax=batch_size,
                    webenv=search_results["WebEnv"],
                    query_key=search_results["QueryKey"])
            attempt = 4
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                attempt += 1
                time.sleep(15)
            else:
                raise

    with open(path + "recent_asthmacea_papers_batch"+str(batch)+".xml", "w+") as f:
        print("Saving...")
        data = fetch_handle.read()
        fetch_handle.close()
        f.write(data)

    batch += 1

# - - - - - - - - - - - - - - - - - - - - - # 
# Parse XML files for Title and Abstract Content
# - - - - - - - - - - - - - - - - - - - - - # 

# Create files to iterate over
files = []
for r,d,f in os.walk(path):
   for  file in f:
      if ".xml" in file:
         files.append(os.path.join(r, file))

# Create dataframe to save all relevant information
# 1. PMID
# 2. Title
# 3. Abstract
# 4. Date

df = pd.DataFrame(columns=["pmid","title","abstract","date"])

# Iterate through files and save to df
for f in files:
   print("Current file is %s" % (f))
   tree = ET.parse(f)
   root = tree.getroot()

   for pubmedart in root.iter("PubmedArticle"):
      try:
         pmid = pubmedart.find("MedlineCitation/PMID").text 
         title = pubmedart.find("MedlineCitation/Article/ArticleTitle").text
         abstract = pubmedart.find("MedlineCitation/Article/Abstract/AbstractText").text
         date = pubmedart.find("PubmedData/History/PubMedPubDate/Year").text + "/" + pubmedart.find("PubmedData/History/PubMedPubDate/Month").text + "/" + pubmedart.find("PubmedData/History/PubMedPubDate/Day").text
         df = df.append(pd.Series([pmid, title, abstract, date], index=df.columns), ignore_index=True)
      except:
         print("Exception: ", sys.exc_info()[0])

   print("Save complete for %s" % (f))

# - - - - - - - - - - - - - - - - - - - - - # 
# Save DF to CSV file for viewing
# - - - - - - - - - - - - - - - - - - - - - # 

df.to_csv(path+"abstractdf.csv")
