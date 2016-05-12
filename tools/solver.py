import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from timeit import timeit

with open('hsim', 'r') as f:
    full_map = []
    lenght = -1
    for line in f:
        break
    for line in f:

        full_map.append(map(float, line.strip().split(' ')))

        assert(lenght == -1 or len(full_map[-1]) == lenght)
        lenght = len(full_map[-1])
    #print full_map
    A = csr_matrix(full_map)

with open('ms', 'r') as f:
    full_map = []
    lenght = -1
    for line in f:
        break
    for line in f:
        full_map.append(map(float, line.strip().split(' ')))
        assert(lenght == -1 or len(full_map[-1]) == lenght)
        lenght = len(full_map[-1])
    B = csr_matrix(full_map)


print "start"
import time
start_time = time.time()
vals, vecs = eigsh(A, 10, B, which='SA')
print("--- %s seconds ---" % (time.time() - start_time))
print vals