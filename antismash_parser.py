'''
This python scirpt is used for parsing antiSMASH output json file.
Usage: python antismash_parser.py --input antismash.json -T mibig_smiliarity_threshold
Output: an csv file containing antismash_chem_classes,antismash_subchem_classes,on_the_contig_edge,mibig_similarity information

If you find our script useful? Please cite!
Feel free to drop us an e-mail (phou@connect.ust.hk, jeremy.owen@vuw.ac.nz) if you have any question!
'''

#!/usr/bin/envpython
# -*- coding: utf-8 -*-
import json
import sys, os, argparse
import re
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input')
#parser.add_argument('--loc')
parser.add_argument('--T',type=float, default=0.2)
args = parser.parse_args()

json_filename = args.input
root_name = os.path.splitext(json_filename)[0]
stat_filename= root_name + "_metanti.csv"


with open(json_filename) as json_file:
    data = json.load(json_file)
    index=len(data['records'])

similarity_dict = {}
chemclass_dict = {}
chemsubclass_dict = {}
bgclen_dict = {}
edge_dict = {}

# for j in range(index):
similarity_list = []
for j in range(index):
    if 'antismash.modules.clusterblast' in data['records'][j]['modules'].keys():
        for m in data['records'][j]['modules']['antismash.modules.clusterblast']['knowncluster']['results']:
            node1 = data['records'][j]['modules']['antismash.modules.clusterblast']['knowncluster']['record_id']
            region1 = m['region_number']
            region2 = f'{region1:03}'
            # print(region2)
            dict_keys1 = root_name + "." + str(node1) + "." + 'region' + str(region2)
            rank = m['ranking']
            if len(rank) == 0:
                item = (dict_keys1, "", "")
                similarity_list.append(item)
            else:
                mibig = m['ranking'][0][0]['accession']
                genes = len(m['ranking'][0][0]['tags'])
                mibig_types = m['ranking'][0][0]["description"] + "," + m['ranking'][0][0]["cluster_type"]
                temp_l1 = []
                for item in m['ranking'][0][1]['pairings']:
                    temp_l1.append(item[2]['name'])
                similarity = len(list(set(temp_l1))) / genes
                if similarity < args.T or m['ranking'][0][1]["core_gene_hits"] == 0:
                    item = (dict_keys1, "", "")
                else:
                    item = (dict_keys1, mibig + ".0", similarity)
                similarity_list.append(item)

    node2 = data['records'][j]['id']
    for k in data['records'][j]['features']:
        if k['type'] == 'region':
            # subnode=node2.split("_")[0].split('0')[-1]
            region3 = k['qualifiers']['region_number'][0]
            # region1 = m['region_number']
            region4 = f'{int(region3):03}'
            dict_keys2 = root_name + '.' + node2 + '.' + "region" + region4
            location = re.findall(r"\d+\.?\d*", k['location'])
            bgclen = abs(int(location[0]) - int(location[1]))
            bgclen_dict[dict_keys2] = bgclen
            chemclass_dict[dict_keys2] = k['qualifiers']['product']
            edge_dict[dict_keys2] = k['qualifiers']["contig_edge"]


# print(similarity_list)

n = ["cdps", "nrps", "nrps-like", "thioamide-nrp"]
p=["hgle-ks","pks-like","ppys-ks","t1pks","t2pks","t3pks","transat-pks",'transat-pks-like']
r=["bottromycin","cyanobactin","fungal-ripp","glycocin","lap","lantipeptide class i","lantipeptide class ii","lantipeptide class iii","lantipeptide class iv","lantipeptide class v","lassopeptide","linaridin","lipolanthine","microviridin","proteusin","ranthipeptide","ras-ripp","ripp-like","rre-containing","sactipeptide","thioamitides","thiopeptide","bacteriocin","head_to_tail","lanthidin","lanthipeptide","tfua-related","microcin","lanthipeptide"]
s=["amglyccycl","oligosaccharide","saccharide","cf_saccharide"]
t=["terpene"]
for keys,values in chemclass_dict.items():
    if len(values)==1:
        chemsubclass_dict[keys]=values[0]
        for i in values:
            if i.lower() in n:
                chemclass_dict[keys]="NRP"
            elif i.lower() in p:
                chemclass_dict[keys]="Polyketide"
            elif i.lower() in r:
                chemclass_dict[keys]="RiPP"
            elif i.lower() in s:
                chemclass_dict[keys]="Saccharide"
            elif i.lower() in t:
                chemclass_dict[keys]="Terpene"
            else:
                chemclass_dict[keys]="Other"            
    if len(values)>1:
        chemsubclass_dict[keys]="Hybrid"
        chemclass_dict[keys]="Hybrid"
    

for keys,values in edge_dict.items():
    if values == ['True']:
        edge_dict[keys]='True'
    if values == ['False']:
        edge_dict[keys]='False'

# print(chemclass_dict,edge_dict,bgclen_dict)

pd_sim=pd.DataFrame(similarity_list)
pd_sim.columns=['query','mibig','similarity']
# print(pd_sim)
pd_edge=pd.DataFrame.from_dict(edge_dict, orient='index').reset_index()
pd_edge.columns=['query','on_contig_edge']
pd_chem=pd.DataFrame.from_dict(chemclass_dict, orient='index').reset_index()
pd_chem.columns=['query','product']
pd_subchem=pd.DataFrame.from_dict(chemsubclass_dict, orient='index').reset_index()
pd_subchem.columns=['query','subclass']
meta_pds=[pd_chem,pd_subchem,pd_edge,pd_sim]

import functools as ft
meta = ft.reduce(lambda left, right: pd.merge(left, right, on='query'), meta_pds)
meta['genome'] = root_name
meta.to_csv(stat_filename, header=False,index=None)
print(root_name,"Done!")

