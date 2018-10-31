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
