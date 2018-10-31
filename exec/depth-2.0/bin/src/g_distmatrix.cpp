/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/

# include "cell_list.h"

int main(int argc, char* argv[]){
	if (argc != 5){
		cout << "./g_distmatrix input_1.xyz input_2.xyz cut_off output.distmap" << endl;
		exit(1);
	} // end if

	cell_list C;
	float cutoff = atof(argv[3]);
	C.compute_files(cutoff, string(argv[1]), string(argv[2]), string(argv[4]));
	return 0;
} // end main
