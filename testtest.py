edges=[1,2,3,4,5,6,7,8,9]
import multiprocessing as mp


import multiprocessing
import time


def do_something(i):
   print(edges[i])
   time.sleep(2)
   print(edges[i])

def parallel(edges):
    pool = multiprocessing.Pool(2)
    pool.map(do_something, range(0, len(edges)))
    pool.close()

edge="edge_4"
import pandas as pd
R=1
I=1000
clipp=100
min_mapping_quality=20
min_base_quality=0
min_al_len=1000
de_max=0.05
AF=0.1
cl = pd.read_csv("output/clusters/clusters_%s_%s_%s.csv" % (edge, I, AF),keep_default_na=False)
clusters = sorted(set(cl.loc[cl['Cluster'] != 'NA']['Cluster'].values))
#print(clusters)
clusters = sorted(set(cl['Cluster'].values))
print(clusters)