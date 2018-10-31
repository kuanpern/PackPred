/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/

// c++ version of distance matrix script
// Author: Tan Kuan Pern

# include <iostream>
# include <string>
# include <map>
# include <math.h>
# include <fstream>
# include "self_distmatrix.h"
using namespace std;


int main(int argc, char* argv[]){

	if (argc != 4){
		cout << argc << endl;
		cerr << "syntax: program_name input_pdb_file output_distance_matrix cut-off" << endl;
		exit(1);
	} // end if

	string input_fname  = argv[1];
	string output_fname = argv[2];
	float cut_off       = atof(argv[3]);

	self_distmatrix C = self_distmatrix();
	C.compute(input_fname, cut_off, output_fname);

	return 0;
} // end int
