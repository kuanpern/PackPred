#!/usr/bin/python
# parse html table from Packpred result to tsv file

import sys
from bs4 import *
import lxml
from shorthand import *

if len(sys.argv) != 3:
	print 'python' , sys.argv[0], 'input.html', 'output.tsv'
	sys.exit(0)
# end if

infname  = sys.argv[1]
outfname = sys.argv[2]

# read and parse input htmle
html = open(infname).read()
soup = BeautifulSoup(html, 'lxml')
table = soup.find('table', {'id':'result_table_full'})

# cast into list format
import lxml
data = []

trs = table.findAll('tr')
for tr in trs:
    tds = tr.findAll('td')
    data.append([td.text for td in tds])
# end for    

# print output
print_table(data, outfname)
