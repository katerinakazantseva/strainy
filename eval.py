import csv
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mt
from params import *


eval = open('output/eval_cl.txt', 'a')
#eval.write("Edge" + "\t" + "Cluster" + "\t" + "St" + "\n")


import re
from collections import Counter

for edge in edges:
    try:
        cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF), keep_default_na=False)
        clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
        for cluster in clusters:
            res=[]
            for read in cl.loc[cl['Cluster'] == cluster]['ReadName'].values:
                res.append(re.sub(".*\._","",read))
            stat=int(Counter(res).most_common()[0][1])/len(res)
            strain=Counter(res).most_common()[0][0]
            eval.write(edge + "\t" + str(cluster) + "\t" + str(stat) + "\t" + str(strain)+ "\n")


    except(FileNotFoundError): continue
eval.close()