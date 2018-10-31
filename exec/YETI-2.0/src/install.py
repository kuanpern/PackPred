#! /usr/bin/python
import commands

wd = commands.getoutput('pwd')
cmd = 'python src/generate_main.py src/main_template.py '+wd+' main.py'
print cmd
commands.getoutput(cmd)
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

