/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/

# include <iostream>
# include <fstream>
# include <stdlib.h>
# include <stdio.h>
# include <vector>
# include <set>
# include <string>
# include "easystring.h"
# include "easyoption.h"
using namespace std;

int main(int argc,char* argv[]){
	// check inputs
	if (argc == 1){
		cout << "Separate different chains from input PDB file into individual PDB chain files" << endl;
		cout << "OPTION\n-i input.pdb " << endl;
		exit(1);
	} // end if

	// define inputs
	option parser;
	parser.set_arg("-i","input","");

	// parse the command line
	parser.parse(argc,argv);

	// get input filename
	string infilename = parser.value("input");
	if (infilename == ""){
		cerr << "please provide input file (pdb file)" << endl;
		exit(1);
	} // end if

	// deconvolute complex into different chains
	string line, newfilename;
	char chainID_new;
	set <char> chains;
	ifstream fin; fin.open(infilename.c_str());
	while (!fin.eof()){
		getline(fin, line);
		if (line.size() != 0){
			if (line.substr(0,4) == "ATOM"){
				chainID_new = line.substr(21,1)[0];
				chains.insert(chainID_new);
			} // end if
		} // end if
	} // end while
	fin.close(); fin.clear();

	vector <string> chain_fnames;
	ofstream fout[chains.size()];
	map <char, int> chain_map;
	unsigned int tmp = -1;
	for (set <char>::iterator iter = chains.begin(); iter != chains.end(); ++iter){
		tmp = tmp + 1;
		chain_map[*iter] = tmp;
		fout[tmp].open((infilename+(*iter)).c_str());
		chain_fnames.push_back((infilename+(*iter)).c_str());
	} // end for

	ofstream fcomplex; fcomplex.open((infilename+"-complex").c_str());
	fin.open(infilename.c_str());
	while (!fin.eof()){
		getline(fin, line);
		if (line.size() != 0){
			if (line.substr(0,4) == "ATOM"){
				chainID_new = line.substr(21,1)[0];
				fout[chain_map[chainID_new]] << line << endl;
				fcomplex << line << endl;
			} // end if
		} // end if
	} // end while
	fin.close(); fin.clear();
	for (unsigned int i = 0; i < chains.size(); ++i){
		fout[i].close(); fout[i].clear();
	} // end for
	fcomplex.close(); fcomplex.clear();
	return 0;
} // end main
