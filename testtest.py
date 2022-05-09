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


print(set(range(3,5)))
#.intersection(set(range(1,4))))




