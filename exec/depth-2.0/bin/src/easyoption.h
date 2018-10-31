/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/

# ifndef OPTION_H
# define OPTION_H

# include <stdlib.h>
# include <iostream>
# include <map>
# include <string>
using namespace std;

class option{
	// work for string:string only, string:bool (for now)
	map <string, string> flag2arg;
	map <string, string> flag2flag;
	map <string, string> arg_list;
	map <string, bool> flag_list;

	public:
	option(){ // constructor
	} // end atomic depth
	~option(){ //destructor
	} // end ~atomic depth

	// set functions
	void set_arg(string flag, string argname, string default_value){
		flag2arg[flag] = argname;
		arg_list[argname] = default_value;
	} // end insert

	void set_flag(string flag, string flagname, bool default_value){
		flag2flag[flag] = flagname;
		flag_list[flagname] = default_value;
	} // end insert

	void parse(int argc, char* argv[]){
		int i = 0;
		while (i < argc - 1){
			i = i + 1;
			if (argv[i][0] == '-') { // is flag
				// determine is arg or flag
				if (flag2flag.find(string(argv[i])) != flag2flag.end()){ // is flag
					// update flag (reverse default value)
					if (flag_list[flag2flag[string(argv[i])]] == 0){ 
						flag_list[flag2flag[string(argv[i])]] = 1;
					} else {
						flag_list[flag2flag[string(argv[i])]] = 0;
					} // end if
				} else if (flag2arg.find(string(argv[i])) != flag2arg.end()){ // is argument
					// update value
					if (i+1 <= argc){
						arg_list[flag2arg[string(argv[i])]] = string(argv[i+1]);
						i = i + 1; // skip for argument
					} else {
						cerr << "no argument for " << argv[i] << endl;
						cerr << "exiting..." << endl;
						exit(1);
					} // end if
				} else {
					// unknown, exit
					cerr << "unknown argument: " << argv[i] << endl;
					cerr << "exiting..." << endl;
					exit(1);
				} // end if
			} // end if
		} // end while
	} // end parse

	// get functions
	string value(string name){
		if (arg_list.find(name) != arg_list.end()){ // is argument
			return arg_list[name];
		} else {
			cerr << name << ": argument not set" << endl;
			return "";
		} // end if
	} // end argument

	bool ticked(string name){
		if (flag_list.find(name) != flag_list.end()){ // is flag
			return flag_list[name];	
		} else {
			cerr << name << ": flag not set" << endl;
			return 0;
		} // end if
	} // end ticked

}; // end option

# endif
