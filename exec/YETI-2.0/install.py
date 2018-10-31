#! /usr/bin/python
import os
import sys
import commands

curdir = os.getcwd()
# generate main file
with open('src/main_template.py') as fin:
	content = fin.read()
# end with
with open('main.py', 'w') as fout:
	content = content.replace('$wd', curdir)
	fout.writelines(content)
# end with

pyscripts = ['distance_selector.py', 'PDB_selector.py', 'select_atom.py', 'xyz_generator.py']
for p in pyscripts:
	cmd = "cp src/"+p+" exec/"
	print cmd
	commands.getoutput(cmd)
# end for

cpp_scripts = {'clique.cpp':'build_clique', 'distmap_sort.cpp':'sort_distmatrix', 'linkage.cpp':'linkage', 'g_distmatrix.cpp':'g_distmatrix'}
for p in cpp_scripts.keys():
	cmd = "g++ src/"+p+" -o exec/"+cpp_scripts[p]
	print cmd
	commands.getoutput(cmd)
# end for

