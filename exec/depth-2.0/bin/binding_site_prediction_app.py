# This is part of DEPTH.
# DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
# Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan
# 
# DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.


# this script compute site and write the output into a formatted html file
from apps_backend import *
from numpy import *
import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import commands
import sys

def predict_binding_site(in_site_parameters, out_site_parameters, job_parameters):
	# get input
	gpdb_name, workdir, exedir, pka_exedir, cavity_exedir, jmol_dir, home_dir, home_name = job_parameters
	pdb_name = gpdb_name.split('/')[-1]

	depth_exe	  = exedir+'/DEPTH'
	ASA_exe		= exedir+'/ASA'
	hotspot_exe	= exedir+'/predict-binding-site'
	distmatrix_exe = exedir+'/distmatrix_self'
	bin_dir		= exedir

	# output file names
	cavity_pred_out = pdb_name + '.pred'
	cavity_pred_pdb = pdb_name + '.pred.pdb'
	cavity_site_pdb = pdb_name + '.binary.pdb'

	#  input parameters
	in_site_threshold, in_site_cavity_size, in_site_iteration_sol, in_site_dist_sol, in_site_iteration_depth, in_site_solN, in_site_resolution, in_site_probe_radius, use_msa = in_site_parameters

	# output parameters
	out_site_jmol, out_site_depth_asa = out_site_parameters

	# run predict binding site program
	use_msa = str(int(use_msa))
	inputs = ['cd', workdir, ';', 'python', exedir+ 'predict_site.py', pdb_name, depth_exe, ASA_exe, hotspot_exe, bin_dir, distmatrix_exe, in_site_threshold, in_site_cavity_size, cavity_pred_pdb, 
cavity_pred_out, bin_dir, in_site_iteration_sol, in_site_dist_sol, in_site_iteration_depth, in_site_solN, in_site_resolution, in_site_probe_radius, use_msa]

	cmd = ' '.join(inputs)
	fout = open(gpdb_name+'.err', 'w')
	fout.writelines(cmd)
	commands.getoutput(cmd)

	# Generate output
	# generate navigation list for requested outputs
	B = bookmarks()
	Content = []
	link_top = 'top' # hard-coding for top bookmark here !
	link1, link2, link3 = [cavity_pred_pdb, cavity_pred_out, cavity_site_pdb]

	inputs = ['cd', workdir, ';', 'python', exedir+'get_binding_residues.py', cavity_site_pdb]
	cmd = ' '.join(inputs)
	fout.writelines('\n'+cmd+'\n')
	binding_residues = commands.getoutput(cmd)
	fout.writelines(binding_residues+'\n')
	fout.close()

#	recommended_threshold = {2:0.36, 3:0.315, 4:0.31, 5:0.33}
#	recommended_threshold = {2:0.36, 3:0.315, 4:0.50, 5:0.45}
#	thresholds = [recommended_threshold[int(float(in_site_solN))], float(in_site_threshold)]
#	labels	 = ['Recommended Threshold', 'User Defined Threshold']
#	hlines = [thresholds, labels]
	img1 = plot_site_prediction(workdir+'/'+cavity_pred_out)#, hlines)

	# generate jmol script
	pdbfname = cavity_site_pdb.split('/')[-1]
	script1 = gen_site_jmol_script(pdbfname, "../../jmol-12.0.22/")

	B.add('out_site_result',	   'BINDING SITE PREDICTION')
	if out_site_jmol == True:
		Content.append(site_output_II('BINDING SITE PREDICTION', 'out_site_result', binding_residues, link_top, link1, link2, link3, img1, script1, [home_dir, home_name]))
	else:
		Content.append(site_output_I ('BINDING SITE PREDICTION', 'out_site_result', binding_residues, link_top, link1, link2, link3, img1, [home_dir, home_name]))
	# end if

	# generate content and write to html file
	Navigator = '\n'.join(B.generate())+'\n\n'
	Content  = '\n\n'.join(Content)

	site_content_outfile = gpdb_name+'.site-content.html'
	fout = open(site_content_outfile, 'w')
	fout.writelines(Content)
	fout.close()

	site_navigation_outfile = gpdb_name+'.site-navigation.html'
	fout = open(site_navigation_outfile, 'w')
	fout.writelines(Navigator)
	fout.close()

	site_outfiles = [site_content_outfile, site_navigation_outfile]

	return site_outfiles

# end def

def plot_site_prediction(cavity_pred_out, hlines = None):
	# read data
	residue_label, P = zip(*read_table(cavity_pred_out)[1:])
	# quick fix for bug
	P_new = []
	for t in P:
		if isnan(t):
			P_new.append(0)
		else:
			P_new.append(t)
		# end if
	# end for
	
	P = P_new
	

	# output file name
	prediction_plotname = cavity_pred_out+ '.site_prediction.png'
	# generate plots according to requests
	xnames = residue_label
	xlabel, ylabel = ['Residue No.', r'Probability ($\AA$)']

	title = "Probability of Residue Forming a Binding Site"
	plot_type_II(P, xnames, hlines, title, xlabel, ylabel, prediction_plotname, xlim='auto', ylim='auto')

	return prediction_plotname
# end def


def gen_site_jmol_script(pdbfname, jmol_dir):
	lines = []
	lines.append('<script type="text/javascript" language="JavaScript" src="'+jmol_dir+'Jmol.js"></script>')
	lines.append('<script>')
	lines.append('	jmolInitialize("'+jmol_dir+'", "JmolApplet.jar")')
	lines.append("	jmolSetAppletColor('white')")
	lines.append("	jmolApplet(500, 'load "+pdbfname+"; cpk off; frame all; cpk off; wireframe off; spacefill off; cartoon off ; isosurface surf molecular colorscheme bwr property temperature ', 'cavity')")
	lines.append('	jmolBr();')
	lines.append('</script>')
	return '\n'.join(lines)
# end def

def site_output_I(title, bookmark, binding_residues, link_top, link1, link2, link3, img1, dir_des): # html output table for other site
	home_dir, home_name = dir_des
	link1 = link1.replace(home_dir, home_name)
	link2 = link2.replace(home_dir, home_name)
	link3 = link3.replace(home_dir, home_name)
	img1  =  img1.replace(home_dir, home_name)

	out = []
	out.append('<table> <a name="'+bookmark+'"></a>')
	out.append('	<h4>'+title+'</h4>')
	out.append('	<tbody>	<tr><td> download prediction output in <a href="'+link3+'"> PDB </a> format </td></tr>')
	out.append('			<tr><td> download probability output in <a href="'+link1+'"> PDB </a> or <a href="'+link2+'">tab-delimited</a> format ')
	out.append('	</td></tr> <tr><td>')
	out.append('	<table>')
	out.append('		<tr><td> <table border="1"> <tr><td> <img src ="'+img1+'"></img> </td></tr></table> </td></tr>')
	out.append('	</table> ')
	out.append('   	<tr><td> <div> <h4> Predicted Binding Residues </h4><br>  ' + binding_residues +' </div></td></tr>')
	out.append('	</td></tr> <tr><td align="right">')
	out.append('	<a href=#'+link_top+'>[back to top]</a>')
	out.append('	</td></tr>')
	out.append('	</tbody>')
	out.append('</table>')

	return '\n'.join(out)
# end def

def site_output_II(title, bookmark, binding_residues, link_top, link1, link2, link3, img1, script1, dir_des): # html output table for other site
	home_dir, home_name = dir_des
	link1 = link1.replace(home_dir, home_name)
	link2 = link2.replace(home_dir, home_name)
	link3 =	link3.replace(home_dir, home_name)
	img1  =  img1.replace(home_dir, home_name)

	out = []
	out.append('<table> <a name="'+bookmark+'"></a>')
	out.append('	<h4>'+title+'</h4>')
	out.append('	<tbody> <tr><td> download prediction output in <a href="'+link3+'"> PDB </a> format </td></tr>')
	out.append('			<tr><td> download probability output in <a href="'+link1+'"> PDB </a> or <a href="'+link2+'">tab-delimited</a> format ')
	out.append('	</td></tr> <tr><td>')
	out.append('		<tr><td> <div style="float:left; display:block"><table border="1"> <tr><td> '+script1+' </td></tr>')
	out.append('			<tr><td>protein surface: <font color="blue"> blue </font> <br> Residues lining binding cavities: <font color="red">Red</font></td></tr>')
	out.append('	</table></div>')
	out.append('	<div style="float:left; display:block"><img src ="'+img1+'"></img></div>')
	out.append('	<tr><td> <div style="border-style:solid"> <p><b> Predicted Binding Residues </b></p> <p>' + binding_residues +' </p> </div></td></tr>')
	out.append('	</table> ')
	out.append('	</td></tr> <tr><td align="right">')
	out.append('	<a href=#'+link_top+'>[back to top]</a>')
	out.append('	</td></tr>')
	out.append('	</tbody>')
	out.append('</table>')

	return '\n'.join(out)
# end def

