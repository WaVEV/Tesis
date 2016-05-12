import matplotlib.pyplot as plt
import sys
import re

xs = []
ys = []
y_max = 0
x_max = 0
with open(sys.argv[1], "r") as f:
    for line in f:
        x, y = map(int, line.split(' '))
        xs.append(x)
        ys.append(y)
        y_max = max(y_max, y)
        x_max = max(x_max, x)


plt.plot(xs, ys, 'ro')
plt.axis([-1, x_max, -1, y_max])
plt.gca().invert_yaxis()
plt.show()