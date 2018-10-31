/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/

# include <cmath>
# include <string>
# include <map>
# include <vector>
# include <fstream>
# include <cstdlib>
# include "shorthand.h"
# include "easyoption.h"
# include "easystring.h"
# include "modeller_psa.h"
# include "Depth_parser.h"
# include "PDB.h"
# include "install_dir.h" // generated on compilation to specify install directory
using namespace std;

float binning(float s, float bin_size = 0.25){
	return int(s/bin_size)*bin_size;
} // end def

// regression line takes in the form of Y = A*exp(-B*X) + C, with stdev D
class regression_data{
	private:
	float bin, predicted_depth, dev, max_bin, min_bin;
	map <float, float> _mean, _stdev;
	set <float> keys;
	string name;

	public:
	regression_data(){
		min_bin = 999; max_bin = -999;
	} // end regression_data
	~regression_data(){}

	// set & get functions
	void set_name(string s){
		name = s;
	} // end set_name
	void set_mean(float a, float b){
		keys.insert(a);
		if (a > max_bin){
			max_bin = a;
		} // end if
		if (a < min_bin){
			min_bin = a;
		} // end if
		_mean[a] = b;
	} // end set_mean
	void set_stdev(float a, float b){
		keys.insert(a);
		_stdev[a] = b;
	} // end set_mean

	// return deviation for an inquiry (note: binning(discrete) used here)
	float deviation(float depth_in, float asa_in){
		if (depth_in == -99){ // error code
			return -99;
		} // end if
		bin = binning(depth_in);
		if (bin < min_bin){ // for instances not included in database
			cout << "warning: residue depth " << depth_in << " < " << min_bin << " which is smallest instance in database. Default it to " << min_bin << endl;
			bin = min_bin;
		} else if (bin > max_bin){
			bin = max_bin;
			cout << "warning: residue depth " << depth_in << " > " << max_bin << " which is largest instance in database. Default it to " << max_bin << endl;
		} // end if
		predicted_depth = _mean[bin];
		dev = (predicted_depth - asa_in)/_stdev[bin];
		return dev;
	} // end deviation
}; // end class

int main(int argc,char* argv[]){
	if (argc ==1){
		cout << "Usage: find out the deviation of depth and accessibility indicators" << endl;
		cout << "Option:\n-d depth-profile \n-a asa-profile \n-p parameter-file \n-o output \n-c pdb-file (Default NULL) -e output.depth-asa (Default NULL)" << endl;
	} // end if

	option parser;
	parser.set_arg("-d","depth_profile","");
	parser.set_arg("-a","asa_profile","");
	parser.set_arg("-p","parameter_file",LIB_DIR+"/par/survive_2.stats");
	parser.set_arg("-o","output_file","");
	parser.set_arg("-c","in_pdb", "");
	parser.set_arg("-e","out_dev", "");

	parser.parse(argc, argv);

	string depth_profile = parser.value("depth_profile");
	string asa_profile = parser.value("asa_profile");
	string parameter_file = parser.value("parameter_file");
	string outfile = parser.value("output_file");
	string in_pdb = parser.value("in_pdb");
	string out_depth_asa = parser.value("out_dev");

	if ((depth_profile == "") || (asa_profile == "") || (parameter_file == "") || (outfile == "")){
		cerr << "error: please check input files" << endl;
		exit(1);
	} // end if

	string line, resname;
	vector <string> bufferline;
	float d, a;
	// initialize regression lines
	map <string, regression_data> predictor_line;
	predictor_line["ALL"]  = regression_data(); predictor_line["ALL"].set_name("ALL");
	predictor_line["MC"]   = regression_data(); predictor_line["MC"].set_name("MC");
	predictor_line["SC"]   = regression_data(); predictor_line["SC"].set_name("SC");
	predictor_line["SCP"]  = regression_data(); predictor_line["SCP"].set_name("SCP");
	predictor_line["SCNP"] = regression_data(); predictor_line["SCNP"].set_name("SCNP");

	// feed parameters into each line
	ifstream fparameter; fparameter.open(parameter_file.c_str());
	while (!fparameter.eof()){
		getline(fparameter, line);
		if (line.size() == 0){
			continue;
		} // end if
		
		bufferline = split(line, " ");
		predictor_line[bufferline[0]].set_mean(atof(bufferline[1].c_str()), atof(bufferline[2].c_str()));
		predictor_line[bufferline[0]].set_stdev(atof(bufferline[1].c_str()), atof(bufferline[3].c_str()));
	} // end while
	fparameter.close(); fparameter.clear();

	// read residue-depth data
	Depth D = Depth(depth_profile.c_str());
	ASA A = ASA(asa_profile.c_str());

	// get all residue names (according to depth)
	map <string, vector <string> > prediction;
	vector <string> all_residues = D.all_residues();
	map <string, int> finder;

	for (unsigned int i = 0; i < all_residues.size(); ++i){
		resname = all_residues[i];
		finder[resname] = i;
		// all-atom
		d = D.residue(resname).depth_all(); a = A.residue(resname).all_per();
		prediction["ALL"].push_back(strf(predictor_line["ALL"].deviation(d, a),2));
		// MC-atom
		d = D.residue(resname).depth_MC(); a = A.residue(resname).main_per();
		prediction["MC"].push_back(strf(predictor_line["MC"].deviation(d, a),2));
		// SC-atom
		d = D.residue(resname).depth_SC(); a = A.residue(resname).side_per();
		prediction["SC"].push_back(strf(predictor_line["SC"].deviation(d, a),2));
		// SCP-atom
		d = D.residue(resname).depth_SCP(); a = A.residue(resname).polar_per();
		prediction["SCP"].push_back(strf(predictor_line["SCP"].deviation(d, a),2));
		// SCNP-atom
		d = D.residue(resname).depth_SCNP(); a = A.residue(resname).nonP_per();
		prediction["SCNP"].push_back(strf(predictor_line["SCNP"].deviation(d, a),2));
	} // end for


	ofstream fout; fout.open(outfile.c_str());
	fout << "# RES | ALL | MC | SC | SCP | SCNP" << endl;
	for (unsigned int i = 0; i < all_residues.size(); ++i){
		fout << all_residues[i] << " " << prediction["ALL"][i] << " " << prediction["MC"][i] << " " << prediction["SC"][i] << " " << prediction["SCP"][i] << " " << prediction["SCNP"][i] << endl;
	} // end for
	fout.close(); fout.clear();

	// color pdb if requested
	if (in_pdb != ""){
		string resSeq;
		float value;
		PDB mdl = PDB(in_pdb, 0, 1);
		for (int i = 0; i < mdl.size(); ++i){
			resSeq = mdl.chainID(i)+":"+mdl.resSeq(i);
			value = atof(prediction["ALL"][finder[resSeq]].c_str());
			mdl.set_T(i, -value);
		} // end for
		mdl.write(in_pdb+".deviation");
	} // end if

	if (out_depth_asa != ""){
		ofstream fdev; fdev.open(out_depth_asa.c_str());
		fdev << "# RES | ALL | MC | SC | SCP | SCNP" << endl;
		for (unsigned int i = 0; i < all_residues.size(); ++i){
			resname = all_residues[i];
			fdev << resname;
			fdev << " | " << D.residue(resname).depth_all() << " " << A.residue(resname).all_per();
			fdev << " | " << D.residue(resname).depth_MC() << " " << A.residue(resname).main_per();
			fdev << " | " << D.residue(resname).depth_SC() << " " << A.residue(resname).side_per();
			fdev << " | " << D.residue(resname).depth_SCP() << " " << A.residue(resname).polar_per();
			fdev << " | " << D.residue(resname).depth_SCNP() << " " << A.residue(resname).nonP_per();
			fdev << endl;
		} // end for
		fdev.close(); fdev.clear();
	} // end if

	return 0;
} // end main
