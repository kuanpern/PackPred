import sys
data = []
fin = open(sys.argv[1])
for line in fin:
	bufferline = line.split()
	data.append(bufferline)
# end for
x = [  int(data[i][1]) for i in range(1, len(data))]
y = [float(data[i][4]) for i in range(1, len(data))]
e = [float(data[i][4]) for i in range(1, len(data))]

import matplotlib.pyplot as plt
plt.grid(True)
plt.title("Difference in Residue Depth upon binding")
plt.xlabel("residue number")
plt.ylabel("delta(depth) (Angstrom)")
plt.errorbar(x,y,e)
outfile = sys.argv[2].replace(".","_")
print "save plot to file:", outfile+".png"
plt.savefig(outfile)
