# include "cell_list.h"

int main(int argc, char* argv[]){
	if (argc != 5){
		cout << "./g_distmatrix input_1.xyz input_2.xyz cut_off output.distmap" << endl;
		exit(1);
	} // end if

	cell_list C;
	float cutoff = atof(argv[3]);
	C.compute_files(cutoff, string(argv[1]), string(argv[2]), string(argv[4]));
	return 0;
} // end main
