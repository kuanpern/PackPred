/*
This is part of DEPTH.
DEPTH (Version: 2.0) computes the closest distance of a residue/atom to bulk solvent and predicts small molecule binding site of a protein. 
Copyright (C) 2013, Kuan Pern Tan, Nguyen Thanh Binh, Raghavan Varadarajan and M.S. Madhusudhan

DEPTH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
DEPTH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with DEPTH.  If not, see <http://www.gnu.org/licenses/>.
*/
std::pair <float,float> rotate(float x, float y, float theta){ //rotation
	float xr, yr;
	xr = x*cos(theta) - y*sin(theta);
	yr = x*sin(theta) + y*cos(theta);
	return std::make_pair(xr,yr);
} // end rotate

std::vector <int> build_neighbour_valid(unsigned int i, int x_box, int y_box, int box_no){
	vector <int> output;
	int xy_box = x_box*y_box;
	int neighbour_buffer[27] = {i, i+1, i+x_box-1, i+x_box, i+x_box+1, i-1, i-x_box-1, i-x_box, i-x_box+1, i-xy_box-x_box-1, i-xy_box-x_box, i-xy_box-x_box+1, i-xy_box-1, i-xy_box, i-xy_box+1, i-xy_box+x_box-1, i-xy_box+x_box, i-xy_box+x_box+1, i+xy_box-x_box-1, i+xy_box-x_box, i+xy_box-x_box+1, i+xy_box-1, i+xy_box, i+xy_box+1, i+xy_box+x_box-1, i+xy_box+x_box, i+xy_box+x_box+1};
	for (int j = 0; j < 27; ++j){
		if (neighbour_buffer[j] >= 0){
			if (neighbour_buffer[j] < box_no+1){
				output.push_back(neighbour_buffer[j]);
			} // end if
		} // end if
	} // end for
	return output;
} // end neighbour

std::vector <int> build_neighbour_2(unsigned int i, int x_box, int y_box, int box_no){
	vector <int> output;
	int xy_box = x_box*y_box; 
	int neighbour_buffer[98] = {i+1+2*x_box, i+1+2*x_box+2*xy_box, i+1+2*x_box+xy_box, i+1+2*x_box-2*xy_box, i+1+2*x_box-xy_box, i+1+2*xy_box, i+1+x_box+2*xy_box, i+1+x_box-2*xy_box, i+1-2*x_box, i+1-2*x_box+2*xy_box, i+1-2*x_box+xy_box, i+1-2*x_box-2*xy_box, i+1-2*x_box-xy_box, i+1-2*xy_box, i+1-x_box+2*xy_box, i+1-x_box-2*xy_box, i+2, i+2*x_box, i+2*x_box+2*xy_box, i+2*x_box+xy_box, i+2*x_box-2*xy_box, i+2*x_box-xy_box, i+2*xy_box, i+2+2*x_box, i+2+2*x_box+2*xy_box, i+2+2*x_box+xy_box, i+2+2*x_box-2*xy_box, i+2+2*x_box-xy_box, i+2+2*xy_box, i+2+x_box, i+2+x_box+2*xy_box, i+2+x_box+xy_box, i+2+x_box-2*xy_box, i+2+x_box-xy_box, i+2+xy_box, i+2-2*x_box, i+2-2*x_box+2*xy_box, i+2-2*x_box+xy_box, i+2-2*x_box-2*xy_box, i+2-2*x_box-xy_box, i+2-2*xy_box, i+2-x_box, i+2-x_box+2*xy_box, i+2-x_box+xy_box, i+2-x_box-2*xy_box, i+2-x_box-xy_box, i+2-xy_box, i+x_box+2*xy_box, i+x_box-2*xy_box, i-1+2*x_box, i-1+2*x_box+2*xy_box, i-1+2*x_box+xy_box, i-1+2*x_box-2*xy_box, i-1+2*x_box-xy_box, i-1+2*xy_box, i-1+x_box+2*xy_box, i-1+x_box-2*xy_box, i-1-2*x_box, i-1-2*x_box+2*xy_box, i-1-2*x_box+xy_box, i-1-2*x_box-2*xy_box, i-1-2*x_box-xy_box, i-1-2*xy_box, i-1-x_box+2*xy_box, i-1-x_box-2*xy_box, i-2, i-2*x_box, i-2*x_box+2*xy_box, i-2*x_box+xy_box, i-2*x_box-2*xy_box, i-2*x_box-xy_box, i-2*xy_box, i-2+2*x_box, i-2+2*x_box+2*xy_box, i-2+2*x_box+xy_box, i-2+2*x_box-2*xy_box, i-2+2*x_box-xy_box, i-2+2*xy_box, i-2+x_box, i-2+x_box+2*xy_box, i-2+x_box+xy_box, i-2+x_box-2*xy_box, i-2+x_box-xy_box, i-2+xy_box, i-2-2*x_box, i-2-2*x_box+2*xy_box, i-2-2*x_box+xy_box, i-2-2*x_box-2*xy_box, i-2-2*x_box-xy_box, i-2-2*xy_box, i-2-x_box, i-2-x_box+2*xy_box, i-2-x_box+xy_box, i-2-x_box-2*xy_box, i-2-x_box-xy_box, i-2-xy_box, i-x_box+2*xy_box, i-x_box-2*xy_box};
	for (int j = 0; j < 98; ++j){
		if (neighbour_buffer[j] >= 0){
			if (neighbour_buffer[j] < box_no+1){
				output.push_back(neighbour_buffer[j]);
			} // end if
		} // end if
	} // end for
	return output;
} // end neighbour
