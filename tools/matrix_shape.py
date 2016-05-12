import matplotlib.pyplot as plt
import sys
import re

max_ncero = 0
id_max_row  = -1
with open(sys.argv[1], "r") as f:
    line = f.readline().strip()
    rows = int(line.split(' ')[0])
    idrow = 0
    xs = []
    ys = []
    for _ in range(rows):
        # ncero = 0
        idrow += 1
        row = map(float, re.sub(" +", " ", f.readline().strip()).split(" "))
        idcol = 0
        ncero = 10000
        for e in row:
            idcol += 1

            if e != 0:
                ncero = min(ncero, idcol)
                xs.append(idcol)
                ys.append(idrow)

        ncero = idrow - ncero + 1
        if max_ncero < ncero:
            id_max_row = idrow
            max_ncero = ncero


print max_ncero, id_max_row
plt.plot(xs, ys, 'ro')
plt.axis([-1, rows, -1, rows])
plt.gca().invert_yaxis()
plt.show()



"""
365 
KORD = 8
L_Interval = 50



1 .. 8
56 .. 62
110 .. 115
163 .. 167
215 .. 218
266 .. 268
316 .. 317
365 .. 365

"""