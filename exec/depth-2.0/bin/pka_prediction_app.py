# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.

# this script compute pKa and write the output into a formatted html file
from apps_backend import *
from numpy import *
import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import commands

def predict_pKa(in_pKa_parameters, out_pKa_parameters, job_parameters):
	# get input
	gpdb_name, workdir, exedir, pka_exedir, cavity_exedir, jmol_dir, home_dir, home_name = job_parameters
	pdb_name = gpdb_name.split('/')[-1]
	#  input parameters
	in_pKa_cycle, in_pKa_solN, in_pKa_EE_cutoff, in_pKa_asa_resolution, in_pKa_probe_radius = in_pKa_parameters
	# output parameters
	out_pKa_ASP, out_pKa_GLU, out_pKa_LYS, out_pKa_HIS = out_pKa_parameters

	# run pKa program
	inputs = ['cd', workdir, ';', 'python', pka_exedir+'pKa_26Dec.py', pdb_name, workdir, pka_exedir, exedir+'/DEPTH', in_pKa_cycle, in_pKa_solN, in_pKa_EE_cutoff, in_pKa_asa_resolution, in_pKa_probe_radius]
	cmd = ' '.join(inputs)

	# output file names

	output_fnames = commands.getoutput(cmd).split('\n')
	output_figure_fnames = []
	pKa_txt_output = "error.txt"
	fout = open(workdir+'/'+"error.txt", 'w')
	fout.writelines(cmd+'\n')
#	fout.writelines(output_fnames+'\n')
	fout.close()


	for fname in output_fnames:
		if fname.startswith('@'):
			pKa_txt_output = fname.split()[1]
		elif fname.startswith('!'):
			output_figure_fnames.append(fname.split()[1])
		# end if
	# end for	

	# Generate output
	# generate navigation list for requested outputs
	B = bookmarks()
	Content = []
	link_top = 'top' # hard-coding for top bookmark here !
	link1 = pKa_txt_output

	img_ASP, img_GLU, img_LYS, img_HIS = ["", "", "", ""]
	for figfname in output_figure_fnames:
		if 'ASP' in figfname:
			img_ASP = figfname
		# end if
		if 'GLU' in figfname:
			img_GLU = figfname
		# end if
		if 'LYS' in figfname:
			img_LYS = figfname
		# end if
		if 'HIS' in figfname:
			img_HIS = figfname
		# end if
	# end for


	if out_pKa_ASP == True:
		B.add('out_pKa_ASP',       'pKa Prediction for ASP')
		Content.append(pKa_output_I('pKa Prediction for ASP', 'out_pKa_ASP', link_top, link1, img_ASP, [home_dir, home_name]))
	if out_pKa_GLU == True:
		B.add('out_pKa_GLU',       'pKa Prediction for GLU')
		Content.append(pKa_output_I('pKa Prediction for GLU', 'out_pKa_GLU', link_top, link1, img_GLU, [home_dir, home_name]))
	if out_pKa_LYS == True:
		B.add('out_pKa_LYS',       'pKa Prediction for LYS')
		Content.append(pKa_output_I('pKa Prediction for LYS', 'out_pKa_LYS', link_top, link1, img_LYS, [home_dir, home_name]))
	if out_pKa_HIS == True:
		B.add('out_pKa_HIS',       'pKa Prediction for HIS')
		Content.append(pKa_output_I('pKa Prediction for HIS', 'out_pKa_HIS', link_top, link1, img_HIS, [home_dir, home_name]))


	# generate content and write to html file
	Navigator = '\n'.join(B.generate())+'\n\n'
	Content  = '\n\n'.join(Content)

	pKa_content_outfile = gpdb_name+'.pKa-content.html'
	fout = open(pKa_content_outfile, 'w')
	fout.writelines(Content)
	fout.close()

	pKa_navigation_outfile = gpdb_name+'.pKa-navigation.html'
	fout = open(pKa_navigation_outfile, 'w')
	fout.writelines(Navigator)
	fout.close()

	pKa_outfiles = [pKa_content_outfile, pKa_navigation_outfile]

	return pKa_outfiles

# end def

def pKa_output_I(title, bookmark, link_top, link1, img1, dir_des): # html output table for other pKa

	home_dir, home_name = dir_des
	link1 = link1.replace(home_dir, home_name)
	img1  =  img1.replace(home_dir, home_name)

	out = []
	out.append('<table> <a name="'+bookmark+'"></a>')
	out.append('	<h4>'+title+'</h4>')
	out.append('	<tbody>	<tr><td> download <a href="'+link1+'">output file</a> in tab-delimited format ')

	out.append('	</td></tr> <tr><td>')
	out.append('	<table>')
	out.append('		<tr><td> <img src ="'+img1+'"></img> </td></tr>')
	out.append('	</table> ')
	out.append('	</td></tr> <tr><td align="right">')
	out.append('	<a href=#'+link_top+'>[back to top]</a>')
	out.append('	</td></tr>')
	out.append('	</tbody>')
	out.append('</table>')

	return '\n'.join(out)

# end def





