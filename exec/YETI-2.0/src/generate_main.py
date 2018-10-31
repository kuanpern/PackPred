import commands
import sys

fin = open(sys.argv[1])
wd = '"'+sys.argv[2]+'"'
fout = open(sys.argv[3], 'w')

for line in fin:
	fout.writelines(line.replace('$wd', wd))
# end for
fin.close()
fout.close()

