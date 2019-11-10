from scipy.stats import unitary_group
import sys

num = int(sys.argv[1])
filename = sys.argv[2]

with open(filename,'w') as f:
    for i in range(num):
        x = unitary_group.rvs(4)
        xr = x.real
        xi = x.imag
        for j in range(4):
            f.write("(%1.10f, %1.10f);(%1.10f, %1.10f);(%1.10f, %1.10f);(%1.10f, %1.10f);\n" % (xr[j][0],xi[j][0],xr[j][1],xi[j][1],xr[j][2],xi[j][2],xr[j][3],xi[j][3]))
