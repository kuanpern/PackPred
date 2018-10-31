# generate output html
import sys
import json
from copy import deepcopy
from collections import OrderedDict
from numpy import *
from bokeh.plotting import ColumnDataSource, figure, show, output_file
from bokeh.models import HoverTool
from bokeh.embed import components
import pandas as pd

if len(sys.argv) != 5:
	print 'python', sys.argv[0], 'input.tsv template.html colorbar.txt output.html'
	sys.exit(1)
# end if

infname        = sys.argv[1]
template_fname = sys.argv[2]
colorbar_fname = sys.argv[3]
outfname       = sys.argv[4]
sel_n = 20

# read data
combined_scores = json.load(open(infname))
mutation_keys = combined_scores.keys()

# sort the data
sorter_keys = []
for key in mutation_keys:
	chain = key.split(':')[0]
	num = key.split(':')[1][1:-1]
	if num[-1] not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']: # fix for iCode
		num = num[:-1]
	# end if
	sorter_keys.append([chain, int(num)])
# end for
mutation_keys = zip(*list(sorted(zip(sorter_keys, mutation_keys))))[1] # sort by (chain, resnum)

# parse data into pandas table
data = {'Mutation': [], 'Score': [], 'Effect': []}
for key in mutation_keys:
	score = combined_scores[key]* 86.5 - 1.15 # <- more readable (and actually corresponding to ddG)
	if score >= -2.0: # <- hard-coded parameter, bad
		eff = 'Neutral'
	else:
		eff = 'Destablizing'
	# end if

	data['Mutation'].append(key)
	data['Score']   .append('%6.4f' % score)
	data['Effect']  .append(eff)
# end for

data_table = zip(*[[key] + data[key] for key in ['Mutation', 'Score', 'Effect']])

data = pd.DataFrame(data)
data = data.set_index('Mutation')

# -- PART A: PLOTTING RESULT -- #

# define color
# RdYlGn style
color_palettes = [['#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#a6d96a', '#66bd63', '#1a9850'],  '#ffffff']
thresholds = [-3.5 , -2.75, -2.0, -1.25, -0.5, 0.25, 1.0]
def color_map(score, color_palettes = color_palettes, thresholds = thresholds):
	colors, null_color = color_palettes
	if isnan(score):
		return null_color
	for i in range(len(thresholds)):
		if score <= thresholds[i]:
			return colors[i]
	return colors[-1]
# end def

#aa_types = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
aa_types = ['G', 'A', 'V', 'L', 'M', 'I', 'S', 'T', 'C', 'P', 'N', 'Q', 'F', 'Y', 'W', 'K', 'R', 'H', 'D', 'E']

# get positions
positions = [dat[:-1] for dat in data.index.values]
positions = [dat.replace(':', ';') for dat in positions] # bokeh does not like ":"
swap = [positions[0]]
for i in range(1, len(positions)):
	if positions[i] not in swap:
		swap.append(positions[i])
positions = swap

# get data
position, aa_list, color, score, effect = [[], [], [], [], []]
for y in positions:
	for x in aa_types:
		aa_list .append(x)
		position.append(y)
		mut = y.replace(';', ':') + x

		if mut.split(':')[1][0] == x:
			mut_score = nan
			eff = 'Native'
		else:
			try:
				mut_score = data['Score'] [mut]
				eff       = data['Effect'][mut]
			except KeyError:
				mut_score = nan
				eff       = 'N/A'
			# end try
		# end if

		mut_score = float(mut_score)
		score.append(mut_score)
		effect.append(eff)
		color.append(color_map(mut_score))
	# end for
# end for

plot_data = {'aa_list': aa_list, 'position': position, 'color': color, 'score': score, 'effect': effect}
source = ColumnDataSource(plot_data)

# define the figure
p = figure(title = '',
	x_range = positions, y_range = list(reversed(aa_types)),
	x_axis_location = 'above', toolbar_location = 'left',
	tools = 'hover',
	plot_width = 15*len(positions), plot_height = 500,)

p.logo = None
p.toolbar_location = None

p.rect('position', 'aa_list', 1, 1, source = source,
	color = 'color', line_color = None)


p.grid.grid_line_color = None
p.axis.axis_line_color = None
p.axis.major_tick_line_color = None
p.axis.major_label_text_font_style = 'bold'
p.axis.major_label_text_font_size = '10pt'
p.axis.major_label_text_font = 'Courier'
p.axis.major_label_standoff = 0
p.xaxis.major_label_orientation = pi/3

hover = p.select(dict(type = HoverTool))
hover.tooltips = OrderedDict([
	('mutation', '@position @aa_list'),
	('score', '@score'),
	('effect', '@effect'),
])

# generate plot's htmls
plot_script, plot_div = components(p)


# PART B - DATATABLE

# generate selected tables
data_f = deepcopy(data_table)
header, data_f = [data_f[0], data_f[1:]]
sorter = [[float(dat[1]), dat] for dat in data_f]
sorter.sort()
data_d = [dat[1] for dat in sorter[:sel_n ]]
data_s = [dat[1] for dat in sorter[-sel_n:]]
sorter = [[abs(float(dat[1])), dat] for dat in data_f]
sorter.sort()
data_n = [dat[1] for dat in sorter[:sel_n ]]

# append header information
data_d, data_s, data_n = [ [header] + dat for dat in [data_d, data_s, data_n]]

def data2htmltable(data, table_id):
	header_html = '<thead><th>' + '</th><th>'.join(header) + '</th></thead>'
	body_html = '<tbody>\n'
	for i in range(1, len(data)):
		dat = data[i]
		body_html += '<tr><td>' + '</td><td>'.join(dat) + '</td></tr>'
	# end for
	body_html += '</tbody>'
	content_html = header_html + body_html

	output_html = '<table id="'+table_id+'">' + content_html + '</table>'
	return output_html
# end def

# generate tables' htmls
table_full_html         = data2htmltable(data_f, 'result_table_full')
table_stablizing_html   = data2htmltable(data_s, 'result_table_stablizing')
table_destablizing_html = data2htmltable(data_d, 'result_table_destablizing')
table_native_html       = data2htmltable(data_n, 'result_table_native')

colorbar_html = '<img alt="colorbar" src="data:image/png;base64,'+open(colorbar_fname).read()+'" style="width:450px" />'

# write output html
html_txt = open(template_fname).read()
html_txt =  html_txt.replace('@plot_result_script@', plot_script) \
					.replace('@plot_html@',          plot_div)    \
					.replace('@colorbar@',       colorbar_html)   \
					.replace('@result_table_full_html@',         table_full_html)         \
					.replace('@result_table_stablizing_html@',   table_stablizing_html)   \
					.replace('@result_table_destablizing_html@', table_destablizing_html) \
					.replace('@result_table_native_html@',       table_native_html)
open(outfname, 'w').writelines(html_txt)









