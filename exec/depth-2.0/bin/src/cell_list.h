/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/
// two layer distance matrix header file
// last revised: Mar 23 2011

// common library
# include <cstdlib>
# include <iostream>
# include <string>
# include <vector>
# include <set>
# include <fstream>
# include <stdio.h>
# include <math.h>
using namespace std;

// user defined libraries
# include "parser.h"
# include "shorthand.h"

float dist(float x1, float y1, float z1, float x2, float y2, float z2){ // Eucliend distance
	float dx = x1 - x2;
	float dy = y1 - y2;
	float dz = z1 - z2;
	return sqrt( dx*dx + dy*dy + dz*dz );
} // end dist

class cell_list{
	private:
	float x_max, y_max, z_max, x_min, y_min, z_min, cutoff;
	unsigned int x_box, y_box, z_box, box_no, labelA_size, labelB_size;
	unsigned int A_size, B_size;
	unsigned int* box_A_size;
	unsigned int* box_B_size;
	float* xA; float* yA; float* zA; float* xB; float* yB; float* zB;
	unsigned int* neighbour_size;
	unsigned int** neighbour_valid;
	unsigned int** box_A; unsigned int** box_B;
	string* labelA; string* labelB;
	set < std::pair <unsigned int, unsigned int> > inter_cells;
	pair <unsigned int, unsigned int>* record_A;
	pair <unsigned int, unsigned int>* record_B;


	// private functions //
	void feed_box_A(){
		unsigned int box_i, current_size;
		// step 1: build box <-> point record
		for (unsigned int i = 0; i < A_size; ++i){
			box_i = find_cell(xA[i], yA[i], zA[i]);
			record_A[i] = make_pair(box_i, i);
		} // end for
		// step 2: parse the record and know size of each box
		for (unsigned int i = 0; i < A_size; ++i){
			box_A_size[record_A[i].first] = box_A_size[record_A[i].first] + 1;
		} // end for
		// step 3: declare memory for each box according to step 2
		for (unsigned int i = 0; i < box_no + 1; ++i){
			box_A[i] = new unsigned int[box_A_size[i]];
		} // end for
		// step 3.5: clear box_A_size record
		for (unsigned int i = 0; i < box_no + 1; ++i){
			box_A_size[i] = 0;
		} // end for
		// step 4: parse the record again, this time assign point to box
		for (unsigned int i = 0; i < A_size; ++i){
			box_i = record_A[i].first; // find the box
			current_size = box_A_size[box_i]; // find the current size of the box
			box_A[box_i][current_size] = record_A[i].second; // update point to the box
			box_A_size[box_i] = box_A_size[box_i] + 1; // update size of the box
		} // end for
	} // end feed_box_A

	void feed_box_B(){
		unsigned int box_i, current_size;
		// step 1: build box <-> point record
		for (unsigned int i = 0; i < B_size; ++i){
			box_i = find_cell(xB[i], yB[i], zB[i]);
			record_B[i] = make_pair(box_i, i);
		} // end for
		// step 2: parse the record and know size of each box
		for (unsigned int i = 0; i < B_size; ++i){
			box_B_size[record_B[i].first] = box_B_size[record_B[i].first] + 1;
		} // end for
		// step 3: declare memory for each box according to step 2
		for (unsigned int i = 0; i < box_no + 1; ++i){
			box_B[i] = new unsigned int[box_B_size[i]];
		} // end for
		// step 3.5: clear box_B_size record
		for (unsigned int i = 0; i < box_no + 1; ++i){
			box_B_size[i] = 0;
		} // end for
		// step 4: parse the record again, this time assign point to box
		for (unsigned int i = 0; i < B_size; ++i){
			box_i = record_B[i].first; // find the box
			current_size = box_B_size[box_i]; // find the current size of the box
			box_B[box_i][current_size] = record_B[i].second; // update point to the box
			box_B_size[box_i] = box_B_size[box_i] + 1; // update size of the box
		} // end for
	} // end feed_box_B


	unsigned int build_neighbour_valid(unsigned int cell_no, unsigned int neighbour_valid_i[]){
	 // immediate neighbour box of a box in 3D
		unsigned int x_i, y_i, z_i;
		unsigned int n = 0;
		revert2_3d(cell_no, x_i, y_i, z_i);
		for (int i = x_i - 1; i <= x_i + 1; ++i){
			if (between(i, 1, x_box)){
				for (int j = y_i - 1; j <= y_i + 1; ++j){
					if (between(j, 1, y_box)){
						for (int k = z_i - 1; k <= z_i + 1; ++k){
							if (between(k, 1, z_box)){
								neighbour_valid_i[n] = find_cell_index(i, j, k);
								n = n + 1;
							} // end if
						} // end for
					} // end if
				} // end for
			} // end if
		} // end for


		return n;
	} // end build_neighbour_valid


	public:
	cell_list(){}
	~cell_list(){}

	void compute_files(float x, string input_fname1, string input_fname2, string outfilename){
		set_cutoff(x);
		read_files(input_fname1, input_fname2);
		build_box();
		assign_points();
		print(outfilename);
	} // end compute

	void compute_vectors(float x, unsigned int size_A, string labelA[], float xA[], float yA[], float zA[], unsigned int size_B, float xB[], float yB[], float zB[], string outfilename){
		set_cutoff(x);
		read_vecs(size_A, labelA, xA, yA, zA, size_B, xB, yB, zB);
		build_box();
		assign_points();
		print(outfilename);
	} // end compute

	void set_cutoff(float x){ // set cutoff threshold for distance computation
		cutoff = x;
	} // end set_cutoff

	unsigned int find_cell(float x, float y, float z){ // find the cell a point will be assigned to
		int x_index = int(( x - x_min) / cutoff) + 1;
		int y_index = int(( y - y_min) / cutoff) + 1;
		int z_index = int(( z - z_min) / cutoff) + 1;
		int cell_index = find_cell_index(x_index, y_index, z_index);
		unsigned int t1, t2, t3;
		revert2_3d(cell_index, t1, t2, t3);
		return cell_index;
	} // end find_cell

	unsigned int find_cell_index(unsigned int x_index, unsigned int y_index, unsigned int z_index){
		return x_index + (y_index-1)*(x_box) + (z_index-1)*(x_box*y_box);
	} // end find_cell_index

	void revert2_3d(unsigned int cell_index, unsigned int &x_index, unsigned int &y_index, unsigned int &z_index){
		unsigned int xy_box = x_box*y_box;
		z_index = 1 + int(cell_index / (xy_box + 0.0));
		cell_index = cell_index - (z_index-1)*(x_box*y_box);
		if (cell_index == 0){
			z_index = z_index - 1;
			cell_index = x_box*y_box;
		} // end if
		y_index = 1 + int(cell_index / (x_box + 0.0));
		cell_index = cell_index - (y_index-1)*(x_box);
		if (cell_index == 0){
			y_index = y_index - 1;
			cell_index = x_box;
		} // end if
		x_index = cell_index;
	} // end revert2_3d

	void read_files(string input_fname1, string input_fname2){ // read 3d vector input file into memory
		// declare some memory
		A_size = line_count(input_fname1);
		B_size = line_count(input_fname2);
		labelA = new string [A_size]; 
		labelB = new string [B_size];
		xA = new float[A_size]; xB = new float[B_size];
		yA = new float[A_size]; yB = new float[B_size];
		zA = new float[A_size]; zB = new float[B_size];
		record_A = new pair <unsigned int, unsigned int> [A_size];
		record_B = new pair <unsigned int, unsigned int> [B_size];

		// read input file into memory
		read_3dvector(input_fname1, "\t", labelA, xA, yA, zA);
		read_3dvector(input_fname2, "\t", labelB, xB, yB, zB);
	} // end initialize


	void read_vecs(unsigned int size_A, string labelA_in[], float xA_in[], float yA_in[], float zA_in[], unsigned int size_B, float xB_in[], float yB_in[], float zB_in[]){
		// declare some memory
		A_size = size_A;
		B_size = size_B;
		string* labelA = new string [A_size]; 
		string* labelB = new string [B_size];
		xA = new float[A_size]; xB = new float[B_size];
		yA = new float[A_size]; yB = new float[B_size];
		zA = new float[A_size]; zB = new float[B_size];
		xA = xA_in; yA = yA_in; zA = zA_in;
		xB = xB_in; yB = yB_in; zB = zB_in;
	} // end vecs


	void build_box(){ // build cell-list

		// step 1: find maximum and minimum
		x_max = max(array_max(A_size, xA), array_max(B_size, xB));
		y_max = max(array_max(A_size, yA), array_max(B_size, yB));
		z_max = max(array_max(A_size, zA), array_max(B_size, zB));
		x_min = min(array_min(A_size, xA), array_min(B_size, xB));
		y_min = min(array_min(A_size, yA), array_min(B_size, yB));
		z_min = min(array_min(A_size, zA), array_min(B_size, zB));

		// step 2: define other box parameters
		x_box = int((x_max - x_min)/cutoff) + 1;
		y_box = int((y_max - y_min)/cutoff) + 1;
		z_box = int((z_max - z_min)/cutoff) + 1;
		box_no = x_box*y_box*z_box;

		// step 3: declare some memory to record neighbouring cells
		neighbour_valid = new unsigned int*[box_no+1]; // record of neighbouring cells
		neighbour_size = new unsigned int[box_no+1]; // number of neighbour of individual cell
		for (unsigned int i = 0; i < box_no + 1; ++i){ // initialize number of neighbour of individual cell to 0s
			neighbour_size[i] = 0;
		} // end for

		unsigned int x_j, y_j, z_j; // !!!
		// step 4: build neighbour list for every box
		for (int i = 0; i < box_no + 1; ++i) {
			neighbour_valid[i] = new unsigned int[27];
			neighbour_size[i] = build_neighbour_valid(i, neighbour_valid[i]);
		} // end for

		// step 5: build non-redundant neighbouring-box pair
		inter_cells.clear();
		for (int i = 1; i < box_no + 1; ++i){
			for (unsigned int n = 0; n < neighbour_size[i]; ++n){
				if (i <= neighbour_valid[i][n]){
					inter_cells.insert(make_pair(i,neighbour_valid[i][n])); // choose inter-cells, prevent double-counting
				} else {
					inter_cells.insert(make_pair(neighbour_valid[i][n],i));
				} // end if
			} // end for
		} // end for

		// allocate size and memory for boxes
		box_A = new unsigned int* [box_no+1]; box_B = new unsigned int* [box_no+1]; // individual cells
		box_A_size = new unsigned int [box_no+1]; box_B_size = new unsigned int [box_no+1]; // sizes of individual cells
		for (unsigned int i = 0; i < box_no + 1; ++i){ // initialize member size of individual cell to 0s
			box_A_size[i] = 0; box_B_size[i] = 0;
		} // end for
	} // end build_box

	void assign_points(){ // assign point to boxes
		for (unsigned int i = 0; i < box_no + 1; ++i){
			box_A_size[i] = 0;
		} // end for 
		feed_box_A();
		for (unsigned int i = 0; i < box_no + 1; ++i){
			box_B_size[i] = 0;
		} // end for 
		feed_box_B();
	} // end assign_point

	void print(string output_fname){ // compute distances and print output to file
		ofstream fout; fout.open(output_fname.c_str()); // open file to write 
		int this_neighbour, index_1, index_2;
		float x1, y1, z1, x2, y2, z2, d;

		for (int i = 1; i < box_no + 1; ++i) { // for every box
			if (box_A_size[i] == 0){ // that is not empty
				continue;
			} // end iff 
			for (unsigned int p = 0; p < box_A_size[i]; ++p){ // for 
				index_1 = box_A[i][p];
				x1 = xA[index_1]; y1 = yA[index_1]; z1 = zA[index_1]; //every contained solute atom
				for (unsigned int n = 0; n < neighbour_size[i]; ++n){ // for every neighbouring 
					this_neighbour = neighbour_valid[i][n]; // solvent cell
					for (unsigned int q = 0; q < box_B_size[this_neighbour]; ++q){ // every its
						index_2 = box_B[this_neighbour][q]; // solvent atom
						x2 = xB[index_2]; y2 = yB[index_2]; z2 = zB[index_2]; //every contained solute atom
						d = dist(x1, y1, z1, x2, y2, z2); // compute their distance = a depth value
						if (d <= cutoff){
							fout << labelA[index_1] << "\t" << labelB[index_2] << "\t" << d << endl;
						} // end if
					} // end for
				} // end for
			} // end for
		} // end for
		fout.close(); fout.clear();
	} // end print

	void print_dimension(){ // print dimension of the cell-list
		cout << "cell dimension: " << x_box << " " << y_box << " " << z_box << endl;
	} // end print_dimension
}; // end class
