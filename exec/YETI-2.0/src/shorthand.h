# ifndef SHORTHAND
# define SHORTHAND

# include <cstdlib>
# include <iostream>
# include <string>
# include <vector>
# include <map>
# include <set>
# include <algorithm>
# include <fstream>
# include <sstream>
# include <stdio.h>
# include <math.h>
# include "easystring.h"
using namespace std;

float array_max(unsigned int size, float x[]){ // find maximum number from an array
	float output = x[0];
	for (unsigned int i = 0; i < size; ++i){
		output = max(output, x[i]);
	} // end for
	return output;
} // end array_max

float array_min(unsigned int size, float x[]){ // find minimum number from an array
	float output = x[0];
	for (unsigned int i = 0; i < size; ++i){
		output = min(output, x[i]);
	} // end for
	return output;
} // end array_max

unsigned int line_count(string filename){ // count non-blank line of a file
	unsigned int output = 0;
	string line;
	ifstream fin; fin.open(filename.c_str());
	while (! fin.eof()){
		getline(fin, line);
		if (line.size() == 0){
			continue;
		} // end if
		output = output + 1;
	} // end while
	fin.close(); fin.clear();
	return output;
} // end line_count

float line(float x1, float y1, float x2, float y2, float x){
	return y2 - ((y2-y1) / (x2-x1)) * (x2-x);
} // end line


double round_double(double x, int dec){
	if (dec <= 0){
		return double(int(x));
	} else {
		int w = int(x);
		int p = 10;
		for (int i = 1; i < dec; ++i){
			p = p*10;
		} // end for
		double decimal = p*(x - w);
		double decimal_new;
		if ( decimal - int(decimal) > 0.5){
			decimal_new = int(decimal + 1);
		} else {
			decimal_new = int(decimal);
		} // end if
		decimal_new = decimal_new/p;
		return w + decimal_new;
	} // end if
} // end round_double

vector <int> set2vec(set <int> S){
	vector <int> output;
	for (set <int>::iterator iter = S.begin(); iter != S.end(); ++iter){
		output.push_back(*iter);
	} // end for
	return output;
} // end set2vec

vector <string> set2vec (set <string> s){
	vector <string> output;
	for (set <string>::iterator iter = s.begin(); iter != s.end(); ++iter){
		output.push_back(*iter);
	} // end for
	return output;
} // end set2vec

vector <pair <string, string> > set2vec(set <pair <string, string> > s){
	vector <pair <string, string> > output;
	for (set <pair <string, string> >::iterator iter = s.begin(); iter != s.end(); ++iter){
		output.push_back(*iter);
	} // end for
	return output;
} // end set2vec

set <int> vec2set(vector <int> input){
	set <int> non_redun_list;
	for (unsigned int i = 0; i < input.size(); ++i){
		non_redun_list.insert(input[i]);
	} // end for
	return non_redun_list;
} // end set2vec



vector <int> non_redundant(vector <int> input){
	sort(input.begin(), input.end());
	vector<int>::iterator it;

	// using default comparison:
	it = unique(input.begin(), input.end());
	input.resize(it - input.begin() );

	return input;
} // end for

vector <string> non_redundant(vector <string> input){
	sort(input.begin(), input.end());
	vector<string>::iterator it;

	// using default comparison:
	it = unique(input.begin(), input.end());
	input.resize(it - input.begin() );

	return input;
} // end for


bool is_element(int i, set <int> box){
	if (box.find(i) != box.end()){
		return 1;
	} else {
		return 0;
	} // end if
} // end is_element


bool is_element(int p, vector <int> box){
	for (unsigned int i = 0; i < box.size(); ++i){
		if (box[i] == p){
			return 1;
		} // end if
	} // end for
	return 0;
} // end is_element

bool is_element(string p, vector <string> box){
	for (unsigned int i = 0; i < box.size(); ++i){
		if (box[i] == p){
			return 1;
		} // end if
	} // end for
	return 0;
} // end is_element

bool is_element(float i, set <float> box){
	if (box.find(i) != box.end()){
		return 1;
	} else {
		return 0;
	} // end if
} // end is_element

// round to nearest integer
int near_round(double s){
	double dec = s - int(s);
	if (dec >= 0.5){
		return int(s) + 1;
	} else {
		return int(s);
	} // end if
} // end near_round

// decimal round
double dec_round(double s, int dec){
	double up = 10;
	for (int i = 0; i < dec - 1; ++i){
		up = up*10;
	} // end for
	int m = near_round(s*up);
	return m / up;
} // end dec_round


bool between(float x, float low, float high){
	if ( (low <= x) && (x <= high) ) {
		return 1;
	} else {
		return 0;
	} // end if
} // end between 

# endif
