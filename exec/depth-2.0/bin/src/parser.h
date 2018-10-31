/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/

# ifndef PARSER
# define PARSER

# include "easystring.h"
# include <string>
# include <fstream>
# include <vector>

unsigned int read_3dvector(string filename, string FS, string label[], float X[], float Y[], float Z[]){ // parser, read labelled 3d vector from a formatted file to memory. Return number of vectors read.
	int n = -1;
	string line;
	vector <string> bufferline;
	ifstream fin; fin.open(filename.c_str());
	while (! fin.eof()){
		getline(fin, line);
		if (line.size() == 0){
			continue;
		} // end if
		n = n + 1;
		bufferline = split(line, FS);
		label[n] = bufferline[0];
		X[n] = atof(bufferline[1].c_str());
		Y[n] = atof(bufferline[2].c_str());
		Z[n] = atof(bufferline[3].c_str());
	} // end while
	fin.close(); fin.clear();
} // end read_vector


# endif
