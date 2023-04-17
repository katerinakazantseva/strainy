import random
import pandas as pd
from strainy.params import *


def trim(cl, data):
    clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA', 'Cluster'].values))
    for cluster in clusters:
        reads_to_trim = []
        st_pos = []
        end_pos = []
        for i in cl.loc[cl['Cluster'] == cluster, 'ReadName'].values:
            st_pos.append(data[i]['Start'])
            end_pos.append(data[i]['Stop'])
        #st_pos=sorted(st_pos)
        #end_pos=sorted(end_pos,reverse=True)
        for i in cl.loc[cl['Cluster'] == cluster, 'ReadName'].values:
            if len([x for x in st_pos if x <= data[i]['Start']])<=2 and data[i]['Start']!<start_end_gap : #replace for param
                reads_to_trim.append(i)
            if len([x for x in end_pos if x >= data[i]['Stop']])<=2 and data[i]['Stop']!=max(end_pos): #replace for param
                reads_to_trim.append(i)
        cl=cl[~cl['ReadName'].isin(reads_to_trim)]
        if len(reads_to_trim)>0:
            cl.loc[cl['Cluster'] == cluster, 'Cluster'] = '0%s' % cluster
    return cl


