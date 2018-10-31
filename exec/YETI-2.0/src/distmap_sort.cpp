# include <iostream>
# include <algorithm>
# include <string>
# include <map>
# include <math.h>
# include <fstream>
# include <stdlib.h>
# include <map>
# include <vector>
# include <set>
using namespace std;

// user-defined library
# include "easystring.h"

class inter_residue{
	private:
	string res1, res2;
	string restype1, restype2;
	float dist;
	map < pair<float, int>, pair<string, string> > mapper;
	vector < pair<float, int> > keytosort;
	pair < string, string >  this_pair;
	set < pair<string, string> > used_pair;
	string atom1, atom2;
		
	public:
	inter_residue(){} // end inter_residue
	~inter_residue(){}

	void set_id(pair<string, string> inter_res_id){
		res1 = inter_res_id.first; res2 = inter_res_id.second;
	} // end set_id

	void set_type(pair<string, string> inter_res_type){
		restype1 = inter_res_type.first; restype2 = inter_res_type.second;
	} // end set_id
	
	void feed(pair<float, int> key1, pair<string, string> key2){
		mapper[key1] = key2; // distance match to atom indices
		keytosort.push_back(key1); // going to sort according to distance
	} // end feed

	void sort_n_print(ofstream &fout){
		sort(keytosort.begin(), keytosort.end()); // sort according to distance
		for (vector < pair<float, int> >::iterator iter = keytosort.begin(); iter != keytosort.end(); ++iter){
			atom1 = mapper[*iter].first; atom2 = mapper[*iter].second; // find the atom indices for this distance
			this_pair = make_pair(atom1, atom2);
			if (used_pair.find(this_pair) == used_pair.end()){ // check for duplicate
				// print to file (follow original format)
				fout << res1+":"+restype1+":"+atom1 << "\t";
				fout << res2+":"+restype2+":"+atom2 << "\t";
				fout << (*iter).first << endl;
				used_pair.insert(this_pair);
			} // end if				
		} // end for
	} // end sort_n_print
}; // end class


int main(int argc, char* argv[]){

	if (argc != 3){
		cout << "./sort_distmatrix input.distmap output.distmap" << endl;
		exit(1);
	} // end if

	// variables
	string line; vector <string> bufferline, dist_element;
	string res1, res2, restype1, restype2, atom1, atom2;
	float d;
	unsigned int n = 0;
	vector < pair<string, string> > inter_res_ids;
	pair <string, string> key2, inter_res_id, inter_res_type;
	pair <float, int> key1;
	map < pair<string,string>, inter_residue > inter_res_record;

	// read distances record
	ifstream fin; fin.open(argv[1]);
	while (! fin.eof()){
		getline(fin, line);
		if (line.size() == 0){
			continue;
		} // end if

		n = n + 1;
		bufferline = split(line, "\t");
		// decode information from distance record
		dist_element = split(bufferline[0], ":");
		res1 = dist_element[0]; restype1 = dist_element[1]; atom1 = dist_element[2];
		dist_element = split(bufferline[1], ":");
		res2 = dist_element[0]; restype2 = dist_element[1]; atom2 = dist_element[2];
		d = atof(bufferline[2].c_str());

		// pack information for memory
		key1 = make_pair(d, n);
		if (res1 < res2){
			inter_res_id = make_pair(res1, res2);
			key2 = make_pair(atom1, atom2);
			inter_res_type = make_pair(restype1, restype2);
		} else {
			inter_res_id = make_pair(res2, res1);
			key2 = make_pair(atom2, atom1);
			inter_res_type = make_pair(restype2, restype1);
		} // end if

		// update memory
		if (inter_res_record.find(inter_res_id) == inter_res_record.end()){ // not found before
			inter_res_ids.push_back(inter_res_id);
			inter_res_record[inter_res_id] = inter_residue(); // build one
			inter_res_record[inter_res_id].set_id(inter_res_id); // set identity
			inter_res_record[inter_res_id].set_type(inter_res_type); // set identity
		} // end if
		inter_res_record[inter_res_id].feed(key1, key2); // distance and atom indexes

	} // end while

	// print output
	sort(inter_res_ids.begin(), inter_res_ids.end()); // sort according to residue label
	ofstream fout; fout.open(argv[2]);
	for (vector < pair <string, string> > ::iterator iter = inter_res_ids.begin(); iter != inter_res_ids.end(); ++iter){
		inter_res_record[*iter].sort_n_print(fout); // sort according to distance (in memory) and print to file
	} // end for
	fout.close(); fout.clear();
	return 0;
} // end for

