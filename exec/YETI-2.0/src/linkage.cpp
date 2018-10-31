# include <set>
# include <iostream>
# include <string>
# include <map>
# include <math.h>
# include <fstream>
# include <stdlib.h>
# include <vector>
# include <string>
using namespace std;

// user-defined library
# include "easystring.h"

// hold information of a residue pair
class linkage{
	private:
	string res1, res2;
	string restype1, restype2;
	set <int> side_chain_atoms_1, side_chain_atoms_2;
	set <string> atoms_1, atoms_2;
	float current_size;
	bool linked;
	int density_1, density_2;
	string atom_1, atom_2;

	public:
	linkage(){
		current_size = 0;
		linked = 0;
	} // end linkage
	~linkage(){}

	void set_id(string res1_input, string res2_input){ // identity of residue-pair
		res1 = res1_input; res2 = res2_input;
	} // end set_id

	void set_type(string restype1_input, string restype2_input){ // a.a. types of residue-pair
		restype1 = restype1_input; restype2 = restype2_input;
	} // end set_type

	void load_density(int density_1_input, int density_2_input){ // threshold density for a.a. types
		density_1 = density_1_input; density_2 = density_2_input;
	} // end load_density

	void insert_data(float d_input, string atom1_input, string atom2_input){ // inter-atomic distance pairs
		atom_1 = atom1_input;
		atom_2 = atom2_input;
		cout << "--" << atom_1 << "-- | --" << atom_2 << "--" << endl;
		if (linked == 0){ // update only if not-linked, ignore otherwise
			current_size = d_input;
			atoms_1.insert(atom_1); atoms_2.insert(atom_2);
			if ((atoms_1.size() >= density_1) && (atoms_2.size() >= density_2)){ // if both residues have enough atoms
				linked = 1; // linked
			} // end if
		} // end if
		cout << atoms_1.size() << "----" << atoms_2.size() << "--" << linked << endl;
	} // end insert_data

	void print(ofstream &fout){ // print output
		if (linked == 1){ // if linked
			cout << "should not reach here" << endl;
			if (res1 < res2){
				fout << res1 << "\t" << res2 << "\t" << current_size << endl;
			} else { // should not be invoked, just to make sure
				fout << res2 << "\t" << res1 << "\t" << current_size << endl;
			} // end if
		} // end if
	} // end print

}; // end linkage


int main(int argc, char* argv[]){

	if (argc != 4){
		cout << "./linkage input.distmap input.density output.link" << endl;
		exit(1);
	} // end if	
	string distmatrix_file = argv[1]; // distance matrix
	string density_file = argv[2]; // density
	string outfilename = argv[3]; // output filename

	// load density file (linkage definition for each amino acid residue type)
	string line;
	vector <string> bufferline;
	map <string, int> density_map;
	ifstream fdensity; fdensity.open(density_file.c_str()); 
	while (!fdensity.eof()){
		getline(fdensity, line);
		if (line.size() == 0){
			continue;
		} // end if
		bufferline = split(line, "\t");
		density_map[bufferline[0]] = atoi(bufferline[1].c_str());
	} // end while
	fdensity.close(); fdensity.clear();

	// read and group distance matrix file line
	string res1, res2, restype1, restype2, atom1, atom2;
	float d;
	vector <string> dist_element;
	ifstream fin; fin.open(distmatrix_file.c_str());
	pair <string, string> inter_res;
	map < pair<string, string>, linkage > inter_res_record;
	while (!fin.eof()){
		getline(fin, line);
		if (line.size() == 0){
			continue;
		} // end if
		bufferline = split(line, "\t");
		// decode information from distance record
		dist_element = split(bufferline[0], ":");
		res1 = dist_element[0]; restype1 = dist_element[1]; atom1 = dist_element[2];
		dist_element = split(bufferline[1], ":");
		res2 = dist_element[0]; restype2 = dist_element[1]; atom2 = dist_element[2];
		d = atof(bufferline[2].c_str());

		// consolidate information into linkage
		inter_res = make_pair(res1, res2);
		if (inter_res_record.find(inter_res) == inter_res_record.end()){ // new inter-residue pair if not found in record
			// declaration and initialize with information
			inter_res_record[inter_res] = linkage();
			inter_res_record[inter_res].set_id(res1, res2);
			inter_res_record[inter_res].set_type(restype1, restype2);
			inter_res_record[inter_res].load_density(density_map[restype1], density_map[restype2]);
		} // end if
		inter_res_record[inter_res].insert_data(d, atom1, atom2); // feed with inter-atomic distance
	} // end while
	fin.close(); fin.clear();

	// print output to file
	ofstream fout; fout.open(outfilename.c_str());
	for (map < pair<string,string>, linkage >::iterator iter = inter_res_record.begin(); iter != inter_res_record.end(); ++ iter){
		inter_res_record[(*iter).first].print(fout);
	} // end for
	fout.close(); fout.clear();

	return 0;
} // end main
