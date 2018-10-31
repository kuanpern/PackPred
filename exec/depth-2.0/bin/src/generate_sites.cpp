/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/

# include <fstream>
# include <iostream>
# include <set>
# include "easyoption.h"
# include "easystring.h"
# include <map>
# include "shorthand.h"
# include <math.h>
# include "install_dir.h" // generated on compilation to specify install directory
using namespace std;

set <string> add_elements(vector <string> first , vector <string> second){
	set <string> output;
	for (unsigned int i = 0; i < first.size(); ++i){
		output.insert(first[i]);
	} // end for
	for (unsigned int i = 0; i < second.size(); ++i){
		output.insert(second[i]);
	} // end for
	return output;
} // end add_elements

int main(int argc,char* argv[]){
	if (argc == 1){
		cout << "Usage: generate predicted cavity sites" << endl;
		cout << "Syntax: ./program-name -i data -o output -p input.pdb [-s] atom-type [-d] cutoff_distance [-T] threshold" << endl;
		cout << "Default:\n\t-s SC\n\t-d 6.2\n\t-T 2.0" << endl;
		exit(1);
	} // end if

	option parser; // declare command line parser
	parser.set_arg("-T","threshold","2");
	parser.set_arg("-d","dist","6.2");
	parser.set_arg("-i","data","");
	parser.set_arg("-p","pdb_name","");
	parser.set_arg("-s","switch","SC");
	parser.set_arg("-o","outfile","");

	parser.parse(argc, argv);

	// get inputs	
	string data_file = parser.value("data");
	string outfile = parser.value("outfile");
	string pdb_name = parser.value("pdb_name");
	if ( (data_file == "") || (outfile == "") || (pdb_name == "")){
		cout << "please provide input and output filenames" << endl;
		exit(1);
	} // end if
	float threshold = atof(parser.value("threshold").c_str());
	float dist = atof(parser.value("dist").c_str());

	// select atom type to compute prediction. Default to side chain (2)
	map <string, int> selector;
	selector["ALL"] = 1; selector["MC"] = 2; selector["SC"] = 3; selector["SCP"] = 4; selector["SCNP"] = 5;
	int select = selector[parser.value("switch")];

	// cache memory
	string line;
	vector <string> bufferline;

	vector <string> contact_res; // read residue that above threshold
	ifstream fdata; fdata.open(data_file.c_str());
	while (!fdata.eof()){
		getline(fdata, line);
		if (line.size() == 0){
			continue;
		} // end if

		bufferline = split(line, " ");
		if (fabs(atof(bufferline[select].c_str())) > threshold){
			contact_res.push_back(bufferline[0]);
		} // end if
	} // end while
	fdata.close(); fdata.clear();

	ifstream fpdb; fpdb.open(pdb_name.c_str()); // write to pseudo pdb file and get distances
	ofstream ftmp_pdb; ftmp_pdb.open((pdb_name + ".tmp").c_str());
	string resSeq, cmd;
	while (!fpdb.eof()){
		getline(fpdb, line);
		if (line.size() == 0){
			continue;
		} // end if
		if (strip(line.substr(0,6), " ") == "ATOM"){
			resSeq = strip(line.substr(21, 22-21))+":"+strip(line.substr(22, 26-22));
			if (is_element(resSeq, contact_res)){
				ftmp_pdb << line << endl;
			} // end if
		} // end if
	} // end while
	fpdb.close(); fpdb.clear();
	ftmp_pdb.close(); ftmp_pdb.clear();

	cmd = EXE_DIR+"/distmatrix " + pdb_name+".tmp " + pdb_name+".tmp.distmap "+ num2string(dist); // get distances
	system(cmd.c_str());

	// read distmap and cluster residues into sites
	ifstream fdist; fdist.open((pdb_name+".tmp.distmap").c_str());
	set <string> core;
	map <string, vector <string> > sites;
	while (!fdist.eof()){
		getline(fdist, line);
		if (line.size() == 0){
			continue;
		} // end if

		bufferline = split(line, " ");
		if (atof(bufferline[10].c_str()) <= dist){
			sites[bufferline[0]].push_back(bufferline[1]);
			sites[bufferline[1]].push_back(bufferline[0]);
			core.insert(bufferline[0]); core.insert(bufferline[1]);
		} // end if
	} // end while
	fdist.close(); fdist.clear();
	// clear temporary file
	remove((pdb_name+".tmp").c_str());
	remove((pdb_name+".tmp.distmap").c_str());

	// non-redundant
	for (set <string>::iterator iter = core.begin(); iter != core.end(); ++iter){
		sites[*iter] = non_redundant(sites[*iter]);
	} // end for

	vector < vector <string> > sorted_sites, output;
	string winner; unsigned int size;

	set <string> tmp_set; 
	unsigned n, site_n = -1;
	bool included;
	int iter = 0;

	while (1) {
		iter = iter + 1;
		cout << "iteration " << iter << endl;
		while (core.size() != 0){ // sort according to size of sites
			winner = *(core.begin()); size = sites[winner].size();
			for (set <string>::iterator iter = core.begin(); iter != core.end(); ++iter){
				if (sites[*iter].size() < size){
					winner = *iter;
					size = sites[winner].size();
				} // end if
			} // end for
			sorted_sites.push_back(sites[winner]);
			core.erase(winner);
		} // end while

		for (unsigned int i = 0; i < sorted_sites.size(); ++i){
			for (unsigned int j = 0; j < sorted_sites[i].size(); ++j){
				cout << sorted_sites[i][j] << " ";
			} // end for
			cout << endl;
		} // end for

		for (unsigned int i = 0; i < sorted_sites.size(); ++i){ // combine sites if
			included = 0;
			n = int(sorted_sites[i].size() / 2 + 1); cout << "*i= " << i << " *n= " << n << endl;;
			for (unsigned int j = i+1; j < sorted_sites.size(); ++j){
				tmp_set = add_elements(sorted_sites[i], sorted_sites[j]);
//				cout << sorted_sites[j].size() << " " << sorted_sites[i].size() << " " << tmp_set.size() << endl;
				if ( (sorted_sites[j].size() + sorted_sites[i].size()) - tmp_set.size() >= n){ // more than half
					sorted_sites[j] = set2vec(tmp_set);
					included = 1;
				} // end if
				tmp_set.clear();
			} // end for
			if (included == 0){
				output.push_back(sorted_sites[i]); // get output
			} // end if
		} // end for

		// prepare for next iteration of sites combination
		core.clear(); sorted_sites.clear(); sites.clear();
		for (unsigned int i = 0; i < output.size(); ++i){
			sites[num2string(i)] = output[i];
			core.insert(num2string(i));
		} // end for

		// convergence = no more site-combining events
		if (output.size() == site_n){
			break;
		} else {
			site_n = output.size();
			output.clear();
		} // end if
	} // end while

	int counter = 0;
	// write output
	ofstream fout; fout.open(outfile.c_str());
	for (unsigned int i = 0; i < output.size(); ++i){
		if (output[i].size() > 1){
			counter = counter + 1;
			fout << "site " << counter << ": ";
			for (unsigned int j = 0; j < output[i].size(); ++j){
				fout << output[i][j] << " ";
			} // end for
			fout << endl;
		} // end if
	} // end for
	fout.close(); fout.clear();

	return 0;
} // end main
