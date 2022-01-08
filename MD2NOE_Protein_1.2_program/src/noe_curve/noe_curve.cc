
// Gordon Chalmers

// inputs -
//
//  argv[1] input file
//  argv[2] maximum mixing time
//  argv[3] inter-ms of mixing time
//  argv[4] jpeg,eps,pdf,epslatex
//  argv[5] output directory

//  argv[6]... in spins
//            out spins



#include "glylib.h"
#include <string>
#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
using std::ios;
using namespace std;
#include <cstring>
#include <math.h>
#include <list>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include "package_inverse_calculation.h"
#include "pthread.h"
#define GNUPLOT "gnuplot -persist"




double string_to_double(const std::string& s)
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
	return 0;
	return x;
}


int string_to_int(const std::string& s)
{
	std::istringstream i(s);
	int x;
	if (!(i >> x))
	return 0;
	return x;
}

string int_to_string(int n)
{
	ostringstream strm;
	strm << n;
	return strm.str ();
}




int main(int argc, char *argv[]) {


	string in_file=argv[1];
	double max_mixing_time=string_to_double(argv[2]);
	double mixing_time_step=string_to_double(argv[3]);
	string font=argv[4];
	string out_directory=argv[5];
	string in_spin=argv[6];
	string out_spin=argv[7];

 	int spinpair=0;
	double r_double=0;
	double r_inv_1_6_double=0;

	// index information about complete relaxation rate matrix

	string text;

	string complete_relaxation_rate_matrix = in_file;

	char *completerelaxationratematrix = new char[complete_relaxation_rate_matrix.length() + 1];
	strcpy(completerelaxationratematrix, complete_relaxation_rate_matrix.c_str());
	ifstream crrm(completerelaxationratematrix);

		string noe_information_file = out_directory;
  	noe_information_file=noe_information_file.append("/noe_information.csv");

		char *noeinformationfile = new char[noe_information_file.length() + 1];
		strcpy(noeinformationfile, noe_information_file.c_str());
		ofstream nif(noeinformationfile);

//    nif << "inputs are :" << endl;
//    nif << "spin system " << in_file << endl;
//    nif << "max mixing time (ms) " << max_mixing_time << endl;
//    nif << "mixing time step (ms) " << mixing_time_step << endl;
//    nif << endl << endl;

//    nif << "spin pair,ResName,ResNumber,Proton,ResName,ResNumber,Proton,r,(<1/r^6>)^(-1/6),source file" << endl << endl;


	for(int i=0; i<3; i++){
		getline(crrm,text);
    cout << text << endl;
	}
  string internal_correlation_directory = text;

	getline(crrm,text);
	cout << text << endl;
  double bfieldz=string_to_double(text);

	for(int i=0; i<12; i++){
		getline(crrm,text);
		cout << text << endl;
	}
	int total_protons=string_to_int(text);

	double **complete_rr_matrix;
	complete_rr_matrix = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	complete_rr_matrix[i] = new double[total_protons];

	double * abbr_complete_rr_matrix = new double[total_protons*total_protons];

	// inverse of complete relaxation rate matrix

	double **inverse_complete_rr_matrix;
	inverse_complete_rr_matrix = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	inverse_complete_rr_matrix[i] = new double[total_protons];

	// initialize

	for(int i=0; i<total_protons*total_protons; i++){
		abbr_complete_rr_matrix[i]=0;
	}

	for(int i1=0; i1<total_protons; i1++){
		for(int i2=0; i2<total_protons; i2++){
			complete_rr_matrix[i1][i2]=0;
			inverse_complete_rr_matrix[i1][i2]=0;
		}
	}


    internal_correlation_directory=internal_correlation_directory.append("/file_distance.txt");

		char *internalcorrelationdirectory = new char[internal_correlation_directory.length() + 1];
		strcpy(internalcorrelationdirectory, internal_correlation_directory.c_str());
		ifstream icd(internalcorrelationdirectory);


  crrm.seekg(0);

	while(text.compare("complete relaxation matrix")!=0){
		getline(crrm,text);
    cout << text << endl;
	}
	getline(crrm,text);

	size_t first=0;
	size_t second=0;
	size_t third=0;
//	size_t fourth=0;
//	size_t fifth=0;

	int proton_spin1[10000];
	int proton_spin2[10000];

	for(int i=0; i<10000; i++){
		proton_spin1[i]=0;
		proton_spin2[i]=0;
	}

  int total=0;

  cout << "test" << endl;

	while(text.compare("empty")!=0){
		getline(crrm,text);

    if(text.compare("empty")!=0){

		first=(size_t) text.find(" ");
		second=(size_t) text.find(" ",first+1);

		proton_spin1[total]=string_to_int(text.substr(0,first));
		proton_spin2[total]=string_to_int(text.substr(first+1,second));
		complete_rr_matrix[proton_spin1[total]][proton_spin2[total]]=string_to_double(text.substr(second+1));

    cout << proton_spin1[total] << " " << proton_spin2[total] << " " << complete_rr_matrix[proton_spin1[total]][proton_spin2[total]] << endl << endl;

    total=total+1;

    }
	}

 total=total-1;

	crrm.seekg(0);

	while(text.substr(0,15).compare("list of protons")!=0){
		getline(crrm,text);
	}

	getline(crrm,text);


	string atom_name_noe[total_protons];
	string residue_name_noe[total_protons];
	string residue_number_noe[total_protons];

	int ith=0;

	for(int i=0; i<total_protons; i++){

		getline(crrm,text);

    cout << text << endl;

		first=text.find(" ");
		second=text.find(" ",first+1);
		ith=string_to_int(text.substr(0,first));
		getline(crrm,text);

    cout << text << endl;

		first=text.find(" ");
		second=text.find(" ",first+1);
		third=text.find(" ",second+1);

		cout << ith << endl;

		atom_name_noe[ith]=text.substr(0,first);
		residue_name_noe[ith]=text.substr(first+1,second-first-1);
		residue_number_noe[ith]=text.substr(second+1,third-second);

		cout << atom_name_noe[ith] << endl;
		cout << residue_name_noe[ith] << endl;

// cout << total << endl;

		cout << residue_number_noe[ith] << endl << endl;

	}


	// noe as a function of time
	// X=-{d/dt} R X
	//  X v + w  initial and final bdry conditions

	// eigenvectors/eigenvalues of complete relaxation rate matrix

	for (int ix = 0; ix < total_protons; ix++){
		for (int iy = 0; iy < total_protons; iy++){
			abbr_complete_rr_matrix[ix*total_protons + iy] = complete_rr_matrix[ix][iy];
		}
	}

	double Z[total_protons*total_protons];
	double eigenvalue[total_protons];

	// function

	inverse_of_matrix(total_protons, (double *)abbr_complete_rr_matrix, (double *)eigenvalue, 10, (double *)Z);

	//    Input, int N, the order of the matrix.
	//    Input, double A[N*N], the real symmetric matrix.
	//    Input, int MATZ, is zero if only eigenvalues are desired,
	// 	   and nonzero if both eigenvalues and eigenvectors are desired.
	//    Output, double W[N], the eigenvalues in ascending order.
	//    Output, double Z[N*N], contains the eigenvectors, if MATZ
	//    	is nonzero.

	// for inspection - eigenvalues


	ofstream temp_file("temp_file_eigenvalues_eigenvectors.txt");

	for (int ix = 0; ix < total_protons; ix++){
		temp_file << eigenvalue[ix] << endl;
		for (int iy = 0; iy < total_protons; iy++){
			temp_file << Z[ix*total_protons + iy] << endl;
		}
	}
	temp_file.close();

	ifstream in_file_temp("temp_file_eigenvalues_eigenvectors.txt");

	double eigenvector[total_protons][total_protons];

	for (int i = 0; i < total_protons; i++){
		getline(in_file_temp, text);
		for (int i2 = 0; i2 < total_protons; i2++){
			getline(in_file_temp, text);
			eigenvector[i2][i] = string_to_double(text);
		}
	}


	// eigenvectors

	double rate[total_protons*total_protons];

	for (int ix = 0; ix < total_protons; ix++){
		for (int iy = 0; iy < total_protons; iy++){
			rate[total_protons*ix + iy] = eigenvector[ix][iy];
		}
	}


	// calculate inverse

	double inverse_matrix[total_protons][total_protons];
	double AInverse[total_protons*total_protons];
	int n = total_protons;
	int i, j, iPass, imx, icol, irow;
	double det, temp, pivot, factor;
	double* ac = (double*)calloc(n*n, sizeof(double));
	det = 1;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			AInverse[n*i + j] = 0;
			ac[n*i + j] = rate[n*i + j];
		}
		AInverse[n*i + i] = 1;
	}

	// The current pivot row is iPass.
	// For each pass, first find the maximum element in the //pivot column.

	for (iPass = 0; iPass < n; iPass++)
	{
		imx = iPass;
		for (irow = iPass; irow < n; irow++)
		{
			if (fabs(rate[n*irow + iPass]) > fabs(rate[n*imx + iPass])) imx = irow;
		}
		// Interchange the elements of row iPass and row imx in both A and AInverse.
		if (imx != iPass)
		{
			for (icol = 0; icol < n; icol++)
			{
				temp = AInverse[n*iPass + icol];
				AInverse[n*iPass + icol] = AInverse[n*imx + icol];
				AInverse[n*imx + icol] = temp;

				if (icol >= iPass)
				{
					temp = rate[n*iPass + icol];
					rate[n*iPass + icol] = rate[n*imx + icol];
					rate[n*imx + icol] = temp;
				}
			}
		}

		// The current pivot is now rate[iPass][iPass].
		// The determinant is the product of the pivot elements.
		pivot = rate[n*iPass + iPass];
		det = det * pivot;
		//		if (det == 0)
		//		{
		//			free(ac);
		//			return 0;
		//		}

		for (icol = 0; icol < n; icol++)
		{
			// Normalize the pivot row by dividing by the pivot element.
			AInverse[n*iPass + icol] = AInverse[n*iPass + icol] / pivot;
			if (icol >= iPass) rate[n*iPass + icol] = rate[n*iPass + icol] / pivot;
		}

		for (irow = 0; irow < n; irow++)
		// Add a multiple of the pivot row to each row.  The multiple factor
		// is chosen so that the element of A on the pivot column is 0.
		{
			if (irow != iPass) factor = rate[n*irow + iPass];
			for (icol = 0; icol < n; icol++)
			{
				if (irow != iPass)
				{
					AInverse[n*irow + icol] -= factor * AInverse[n*iPass + icol];
					rate[n*irow + icol] -= factor * rate[n*iPass + icol];
				}
			}
		}
	}
	free(ac);

	for (int ix = 0; ix < total_protons; ix++){
		for (int iy = 0; iy < total_protons; iy++){
			inverse_matrix[ix][iy] = AInverse[total_protons*ix + iy];
		}
	}

	// inverse completed

	// different inputs
	//  inverted  output  -  inverted-output noe curve
	//  selective selective : input atom  output atom
	//  selective all
	//  all selective
	//  all all


	double **initial_state;
	initial_state = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	initial_state[i] = new double[total_protons];

	double **final_state;
	final_state = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	final_state[i] = new double[total_protons];

	double **noe;
	noe = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	noe[i] = new double[total_protons];

	double **relaxation_rate_mixed;
	relaxation_rate_mixed = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	relaxation_rate_mixed[i] = new double[total_protons];

	double **relaxation_rate;
	relaxation_rate = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	relaxation_rate[i] = new double[total_protons];

	double test=0;

	double **t0;
	t0 = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	t0[i] = new double[total_protons];

	double **tinf;
	tinf = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	tinf[i] = new double[total_protons];

	for(int i=0; i<total_protons; i++){
		for(int i2=0; i2<total_protons; i2++){

			initial_state[i][i2]=0;
			final_state[i][i2]=0;
			noe[i][i2]=0;

			relaxation_rate_mixed[i][i2]=0;
			relaxation_rate[i][i2]=0;

			t0[i][i2]=0;
			tinf[i][i2]=0;
		}
	}

	//	string in_directory=argv[1];
	//	double max_mixing_time=string_to_int(argv[2]);
	//	double mixing_time_step=string_to_int(argv[3]);
	//	string font=argv[4];
	//	string out_directory=argv[5];
	//	string in_spin=argv[6]);
	//	string out_spin=argv[7]);


	for(int atom=0; atom<total_protons; atom++){

if(atom_name_noe[atom].compare("H")==0){

		for(int secondatom=0; secondatom<total_protons; secondatom++){

			if(secondatom!=atom){

				for(int i=0; i<total_protons; i++){
					for(int i2=0; i2<total_protons; i2++){
						t0[i][i2]=0;
						tinf[i][i2]=0;
					}
				}

				for (int ix = 0; ix < total_protons; ix++){
					for (int iy = 0; iy < total_protons; iy++){
						tinf[ix][ix] = 1;
					}
				}

				t0[secondatom][secondatom] = -2;
				tinf[secondatom][secondatom] =1;

				// the second proton spin is inverted

				char *atomName1 = new char[atom_name_noe[secondatom].length() + 1];
				char *atomName2 = new char[atom_name_noe[atom].length() + 1];
				char *atomResidue1 = new char[residue_name_noe[secondatom].length() + 1];
				char *atomResidue2 = new char[residue_name_noe[atom].length() + 1];
				char *residue1 = new char[residue_number_noe[secondatom].length() + 1];
				char *residue2 = new char[residue_number_noe[atom].length() + 1];

				strcpy(atomName1, atom_name_noe[secondatom].c_str());
				strcpy(atomName2, atom_name_noe[atom].c_str());
				strcpy(atomResidue1, residue_name_noe[secondatom].c_str());
				strcpy(atomResidue2, residue_name_noe[atom].c_str());
				strcpy(residue1, residue_number_noe[secondatom].c_str());
				strcpy(residue2, residue_number_noe[atom].c_str());

				string file_txt_of_atom_noe_file = out_directory;
				file_txt_of_atom_noe_file.append("/noecurve_");
				file_txt_of_atom_noe_file.append(residue_name_noe[atom]);
				file_txt_of_atom_noe_file.append("_");
				file_txt_of_atom_noe_file.append(residue_number_noe[atom]);
				file_txt_of_atom_noe_file.append("_");
				file_txt_of_atom_noe_file.append(atom_name_noe[atom]);
				file_txt_of_atom_noe_file.append("-");
				file_txt_of_atom_noe_file.append(residue_name_noe[secondatom]);
				file_txt_of_atom_noe_file.append("_");
				file_txt_of_atom_noe_file.append(residue_number_noe[secondatom]);
				file_txt_of_atom_noe_file.append("_");
				file_txt_of_atom_noe_file.append(atom_name_noe[secondatom]);
				file_txt_of_atom_noe_file.append(".txt");

				char *filetxtofatomnoefile = new char[file_txt_of_atom_noe_file.length() + 1];
				strcpy(filetxtofatomnoefile, file_txt_of_atom_noe_file.c_str());

				string fig = out_directory;
				fig.append("/noecurve_");
				fig.append(residue_name_noe[atom]);
				fig.append("_");
				fig.append(residue_number_noe[atom]);
				fig.append("_");
				fig.append(atom_name_noe[atom]);
				fig.append("_");
				fig.append(residue_name_noe[secondatom]);
				fig.append("_");
				fig.append(residue_number_noe[secondatom]);
				fig.append("_");
				fig.append(atom_name_noe[secondatom]);
        fig.append("_");
        fig.append(int_to_string(atom));
        fig.append("_");
        fig.append(int_to_string(secondatom));

				char *noe_fig = new char[fig.length() + 1];
				strcpy(noe_fig, fig.c_str());

				ofstream pair_of_atom(filetxtofatomnoefile);

				double inttotal=0;

				// noe curve

				for (int total = 0; total < ((double) max_mixing_time / mixing_time_step); total++){

					for (int i = 0; i < total_protons; i++){
						for (int i100 = 0; i100 < total_protons; i100++){
							relaxation_rate_mixed[i][i100] = 0;  relaxation_rate[i][i100] = 0;
						}
					}

					inttotal=(double) total* mixing_time_step;
					inttotal = (double)((double) inttotal / (double) 1000.000000);

					for (int ix = 0; ix < total_protons; ix++){
						test = (double)((double)eigenvalue[ix] * (double)inttotal);
						for (int iy = 0; iy < total_protons; iy++){
							relaxation_rate[ix][iy] = (double) relaxation_rate[ix][iy] + (double) inverse_matrix[ix][iy] * exp(-test);

						}
					}

					for (int iz = 0; iz < total_protons; iz++){
						relaxation_rate_mixed[atom][secondatom] = relaxation_rate_mixed[atom][secondatom] + eigenvector[atom][iz] * relaxation_rate[iz][secondatom];
					}

					noe[atom][secondatom] = relaxation_rate_mixed[atom][secondatom];

					double temp=0;
					temp = (double) noe[atom][secondatom] * t0[secondatom][secondatom] + tinf[atom][secondatom];

					pair_of_atom << (double) total*mixing_time_step << " " << temp << endl;

				}


				FILE *gpnoe;
				gpnoe = popen(GNUPLOT, "w");
				if (gpnoe == NULL){
					return 0;
				}

				if (font.compare("jpeg") == 0){
					fprintf(gpnoe, " set terminal jpeg \n unset key \n set output '%s.jpeg' \n  set title 'atom %s %s:%s - atom %s %s:%s, %.2f T' \n set  xlabel 'mixing time, %.0f ms units \n set ylabel 'intensity' \n plot '%s' \n ", noe_fig, atomResidue2, residue2, atomName2, atomResidue1, residue1, atomName1, bfieldz, mixing_time_step, filetxtofatomnoefile);
				}
				if (font.compare("epslatex") == 0){
					fprintf(gpnoe, " set terminal epslatex \n unset key \n set output '%s.jpeg' \n  set title 'atom %s %s:%s - atom %s %s:%s, %.2f T' \n set  xlabel 'mixing tune, %.0f ms units \n set ylabel 'intensity' \n plot '%s' \n ", noe_fig, atomResidue2, residue2, atomName2, atomResidue1, residue1, atomName1, bfieldz, mixing_time_step, filetxtofatomnoefile);
				}
				if (font.compare("pdf") == 0){
					fprintf(gpnoe, " set terminal pdf \n unset key \n set output '%s.jpeg' \n  set title 'atom %s %s:%s - atom %s %s:%s, %.2f T' \n set  xlabel 'mixing time, %.0f ms units \n set ylabel 'intensity' \n plot '%s' \n ", noe_fig, atomResidue2, residue2, atomName2, atomResidue1, residue1, atomName1, bfieldz, mixing_time_step, filetxtofatomnoefile);
				}
				if (font.compare("eps") == 0){
					fprintf(gpnoe, " set terminal eps \n unset key \n set output '%s.jpeg' \n  set title 'atom %s %s:%s - atom %s %s:%s, %.2f T' \n set  xlabel 'mixing tune, %.0f ms units \n set ylabel 'intensity' \n plot '%s' \n ", noe_fig, atomResidue2, residue2, atomName2, atomResidue1, residue1, atomName1, bfieldz, mixing_time_step, filetxtofatomnoefile);
				}

				fclose(gpnoe);

				pair_of_atom.close();


        icd.seekg(0);

				string text_spin_pair = residue_name_noe[atom];
				text_spin_pair.append(",");
				text_spin_pair.append(residue_number_noe[atom]);
				text_spin_pair.append(",");
				text_spin_pair.append(atom_name_noe[atom]);
				text_spin_pair.append(",");
				text_spin_pair.append(residue_name_noe[secondatom]);
				text_spin_pair.append(",");
				text_spin_pair.append(residue_number_noe[secondatom]);
				text_spin_pair.append(",");
				text_spin_pair.append(atom_name_noe[secondatom]);

text="0";

string text2;

while (text==text_spin_pair){

    getline(icd,text);
  	getline(icd,text2);

		first=(size_t) text2.find(" ");
		second=(size_t) text2.find(" ",first+1);

		r_double=string_to_double(text2.substr(0,first));
		r_inv_1_6_double=string_to_double(text2.substr(first+1,second));

}

nif << spinpair << "," << residue_name_noe[atom] << "," << residue_number_noe[atom] << "," << atom_name_noe[atom] << "," << residue_name_noe[secondatom] << "," << residue_number_noe[secondatom] << "," << atom_name_noe[secondatom] << "," << r_double << "," << r_inv_1_6_double << "," << file_txt_of_atom_noe_file << endl;

// cout << spinpair << "," << residue_name_noe[atom] << "," << residue_number_noe[atom] << "," << atom_name_noe[atom] << "," << residue_name_noe[secondatom] << "," << residue_number_noe[secondatom] << "," << atom_name_noe[secondatom] << "," << r_double << "," << r_inv_1_6_double << "," << file_txt_of_atom_noe_file << endl;


spinpair=spinpair+1;

			}

		}

}

	}

 icd.close();

 nif.close();


	return 0;

}
