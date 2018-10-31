// format of clique file = residues split by "\t", last column = size
# include <algorithm>
# include <set>
# include <stdlib.h>
# include <iostream>
# include <string>
# include <vector>
# include <map>
# include <fstream>
# include <sstream>
using namespace std;

// user-defined library
# include "easystring.h"

class clique{
	private:
	vector <string> clique_residues;
	float initial_size, new_size;
	map < pair<string, string>, float >::iterator iter;
	bool to_insert;
	vector <string> output_residue_list;
	int clique_order;

	public:
	clique(vector <string> clique_residues_input, float size){
		clique_residues = clique_residues_input;
		initial_size = size;
		clique_order = clique_residues.size();
	} // end clique
	~clique(){}

	pair < string, string > residue_pair;
	// generate higher order cliques based on current clique and linkage matrix
	vector <clique> derivative(map < pair<string, string>, float > &linkage, vector <string> &residue_list){
		vector <clique> output; // output = a set of cliques
		string residue_new;

		for (unsigned int i = 0; i < residue_list.size(); ++i){
			residue_new = residue_list[i];
			if (residue_new <= clique_residues[clique_residues.size() - 1]){ // this residue considered already
				continue; // skip
			} // end if

			new_size = initial_size;
			// test linkage for every residue already in clique
			to_insert = 1;
			for (int j = 0; j < clique_residues.size(); ++j){
				residue_pair = make_pair(clique_residues[j], residue_new);
				iter = linkage.find(residue_pair);
				if (iter == linkage.end()){ // new residue not linked to any of residues in current clique
					to_insert = 0;
					break;
				} else {
					new_size = max(new_size, (*iter).second); // update clique size
				} // end if
			} // end for
			if (to_insert == 1){ // if pass all residue tests
				output_residue_list.clear();
				for (int j = 0; j < clique_residues.size(); ++j){
					output_residue_list.push_back(clique_residues[j]);
				} // end for
				output_residue_list.push_back(residue_new); // insert
				output.push_back(clique(output_residue_list, new_size)); // append to result
			} // end if
		} // end for
		return output;
	} // end derivative

	int order(){
		return clique_order;
	} // end if

	void print(ofstream &fout){
		for (int i = 0; i < clique_residues.size(); ++i){
			fout << clique_residues[i] << "\t";
		} // end for
		fout << initial_size << endl;
	} // end print
}; // end clique

int main(int argc, char* argv[]){
	// buffer variables
	string line; vector <string> bufferline;
	float clique_size; 	int order;
	ofstream fout, forphan;

	if (argc != 4) {
		cerr << "Usage: program_name linkage_matrix-file proto-clique outfilename " << endl;
		exit(1);
	} // end if

	cout << "step 0: read input parameters" << endl;
	string linkage_file = argv[1];
	string protoclique_file = argv[2];
	string outfilename = argv[3];

	// declare memories
	set <string> residue_list_set; // residue list
	map < pair<string, string>, float > linkage; // linkage
	vector <clique> clique_record; // cliques
	vector <clique> extended_clique_record;
	vector <clique> extended_cliques;
	vector <clique> orphan_cliques;

	cout << "step 1: read linkage file" << endl;
	ifstream flink; flink.open(linkage_file.c_str());
	while (!flink.eof()){
		getline(flink, line);
		if (line.size() == 0){
			continue;
		} // end if
		bufferline = split(line, "\t");
		// update to linkage record (and residue list)
		linkage[make_pair(bufferline[0], bufferline[1])] = atof(bufferline[2].c_str()); // read into linkage
		residue_list_set.insert(bufferline[0]); residue_list_set.insert(bufferline[1]); // read element into residue list
	} // end while
	flink.close(); flink.clear();

	cout << "step 1.2: sorting residue list" << endl;
	vector <string> residue_list;
	for (set <string>::iterator iter = residue_list_set.begin(); iter != residue_list_set.end(); ++iter){
		residue_list.push_back(*iter);
	} // end for
	sort (residue_list.begin(), residue_list.end());

	cout << "step 2: read protoclique file" << endl;
	ifstream fclique; fclique.open(protoclique_file.c_str());
	while (!fclique.eof()){
		getline(fclique, line);
		if (line.size() == 0){
			continue;
		} // end if
		bufferline = split(line, "\t");

		clique_size = atof(bufferline[bufferline.size() - 1].c_str());
		bufferline.pop_back(); // (remove last element = size)
		clique_record.push_back(clique(bufferline, clique_size)); // read into protoclique
	} // end while
	flink.close(); flink.clear();
	if (clique_record.size() == 0){
		cout << "empty proto-clique file" << endl;
		exit(1);
	} // end if

	cout << "step 3: generate derivative of cliques" << endl;
	for (unsigned int i = 0; i < clique_record.size(); ++i){ // for every clique
		extended_cliques.clear();
		extended_cliques = clique_record[i].derivative(linkage, residue_list); // generate derivative clique

		// analyse derivative cliques
		if (extended_cliques.size() == 0){ // record cliques that do not have derivatives
			orphan_cliques.push_back(clique_record[i]);
		} else {
			for (unsigned int j = 0; j < extended_cliques.size(); ++j){
				extended_clique_record.push_back(extended_cliques[j]); // record derivative cliques
			} // end for
		} // end for
	} // end for

	cout << "step 4: print output to file" << endl;
	cout << "step 4.1: print orphan clique" << endl;
	order = orphan_cliques[0].order();
	forphan.open((outfilename+"-orphan-"+num2string(order)).c_str());
	for (int i = 0; i < orphan_cliques.size(); ++i){
		orphan_cliques[i].print(forphan);
	} // end for
	forphan.close(); forphan.clear();

	cout << "step 4.2: print derivative clique" << endl;
	if (extended_clique_record.size() == 0){ // no derivatives
		order = clique_record[0].order();
	} else {
		order = extended_clique_record[0].order();
	} // end if
	fout.open((outfilename+num2string(order)).c_str());
	for (int i = 0; i < extended_clique_record.size(); ++i){
		extended_clique_record[i].print(fout);
	} // end for
	fout.close(); fout.clear();

	return 0;
} // end 

