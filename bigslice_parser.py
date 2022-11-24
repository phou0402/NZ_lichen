'''
This python scirpt is adapted from https://bioinformatics.nl/~kauts001/ltr/bigslice/paper_data/scripts/ and used for parsing BiGSLiCE output sql file.

Usage: python antismash_parser.py --input antismash.json -T mibig_smiliarity_threshold
Output: an csv file containing bigslice_class, bigslice_subclass, gcf, gcf_value, on_contig_edge, len_nt information

If you find our script useful? Please cite!
Feel free to drop us an e-mail (phou@connect.ust.hk, jeremy.owen@vuw.ac.nz) if you have any question!
'''

#!/usr/bin/envpython
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import os, argparse,sys
import csv, sqlite3

db_filename = sys.argv[1]

con = sqlite3.connect(db_filename) # change to 'sqlite:///your_filename.db'
cur = con.cursor()
cur.execute("DROP TABLE IF EXISTS chem_class;")
cur.execute("DROP TABLE IF EXISTS chem_subclass;")
cur.execute("DROP TABLE IF EXISTS chem_subclass_map;")
cur.execute("CREATE TABLE IF NOT EXISTS chem_class (id, name);")
cur.execute("CREATE TABLE IF NOT EXISTS chem_subclass (id,class_id,name);")
cur.execute("CREATE TABLE IF NOT EXISTS chem_subclass_map (class_source,type_source,subclass_id);")# use your column names here
# cur.execute("CREATE TABLE IF NOT EXISTS hmm (id, accession, name, db_id, model_length);")

with open('chem_class.csv','r') as fin: # `with` statement available in 2.5+
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(int(i['id']), i['name']) for i in dr]
cur.executemany("INSERT INTO chem_class (id, name) VALUES (?, ?);", to_db)

with open('chem_subclass.csv','r') as fin: # `with` statement available in 2.5+
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(int(i['id']),int(i['class_id']),i['name']) for i in dr]
cur.executemany("INSERT INTO chem_subclass (id,class_id,name) VALUES (?, ?, ?);", to_db)

with open('chem_subclass_map.csv','r') as fin: # `with` statement available in 2.5+
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(i['class_source'], i['type_source'],int(i['subclass_id'])) for i in dr]
cur.executemany("INSERT INTO chem_subclass_map (class_source,type_source,subclass_id) VALUES (?, ?, ?);", to_db)

con.commit()
con.close()

########################################################################################################################

root_name = os.path.splitext(db_filename)[0]
stat_filename = root_name + "_metadata.csv"
stat_filename2 = root_name + "_metadata-complete.csv"
output_handle = open(stat_filename, "w",newline='')

with sqlite3.connect(db_filename) as con:
    cur = con.cursor()
    print("loading clustergbk metadata..")
    file_names, bgc_ids,gcf_ids,gcf_values,contig_edges, length_nts,folder = list(zip(*cur.execute(
                "select bgc.orig_filename, bgc.id,gcf_membership.gcf_id,gcf_membership.membership_value,bgc.on_contig_edge, bgc.length_nt,bgc.orig_folder"
                " from bgc,gcf_membership"
                " where bgc.id=gcf_membership.bgc_id"
                " order by gcf_membership.gcf_id"
            ).fetchall()))

    print("loading class information..")
    class_titles = sorted(set([
        "{}:{}".format(class_name, subclass_name) \
        for class_name, subclass_name in cur.execute(
            "select chem_class.name, chem_subclass.name from chem_subclass, chem_class"
            " where chem_class.id=chem_subclass.class_id"
        ).fetchall()]))
    class_presences = {}
    for class_title in class_titles:
        class_name, subclass_name = class_title.split(":")
        subclass_presences = pd.Series(
            np.full(len(bgc_ids), False),
            index=bgc_ids
        )

        try:
            subclass_bgc_ids, = list(zip(*cur.execute((
                "select distinct bgc_class.bgc_id from bgc_class, chem_subclass, chem_class "
                "where chem_subclass.class_id=chem_class.id "
                "and bgc_class.chem_subclass_id=chem_subclass.id "
                "and chem_class.name like ? and chem_subclass.name like ?"
            ), (class_name, subclass_name)).fetchall()))
            subclass_presences[subclass_bgc_ids] = True
        except:
            pass
        class_presences[class_title] = subclass_presences

    print("merging datasets...")
    temp_metadata = pd.DataFrame({
        "bgc": file_names,
        "bgc_id":bgc_ids,
        "gcf": gcf_ids,
        "gcf_value":gcf_values,
        "contig_edge": contig_edges,
        "len_nt": length_nts,
        "folder":folder
    }, index=gcf_ids)

    for class_title in sorted(class_titles):
        temp_metadata["class-" + class_title] = class_presences[class_title].values

temp_metadata.to_csv("temp_query-metadata.csv", sep=",",index=False)

class_dict={}
dict_from_csv = pd.read_csv('temp_query-metadata.csv', header=0, index_col=0).T.to_dict()

for k_mibig in dict_from_csv.keys():
    for k,v in dict_from_csv[k_mibig].items():
        if 'class' in k and v==True:
            class_dict[k_mibig]=k
temp_class= pd.DataFrame.from_dict(class_dict, orient="index").reset_index()
temp_class.columns = ['bgc','class']
merged = pd.merge(temp_class[['bgc','class']],temp_metadata[['bgc','bgc_id','gcf','gcf_value','contig_edge','len_nt','folder']],how='left')
merged[['class', 'subclass']] = merged['class'].str.split(':', expand=True)
merged.to_csv(output_handle,index=False)
os.remove("temp_query-metadata.csv")

output_handle.close()
