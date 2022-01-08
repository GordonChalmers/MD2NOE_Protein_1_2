
// Gordon Chalmers

// inputs -
//     
//	argv[1] input directory
//      argv[2] magnetic field in tesla
//      argv[3] maximum frequency of output Fourier transform in MHz
//	argv[4] jpeg,eps,pdf,epslatex
//	argv[5] output directory
//      argv[6] tumbling


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
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include "package_inverse_calculation.h"
# include "pthread.h"
# define GNUPLOT "gnuplot -persist"




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


	string in_directory=argv[1];
	double magnetic_field=string_to_double(argv[2]);
	double max_frequency=string_to_double(argv[3]);
	string font=argv[4];
	string out_directory=argv[5];
	string tumbling=argv[6];

	string text, text_total_protons;

	
	// correlation function information file

	string correlation_function_file_information= in_directory;
	correlation_function_file_information.append("/internal_correlation_function_");
	correlation_function_file_information.append(in_directory);
	correlation_function_file_information.append(".txt");
	char *correlationfunctionfileinformation = new char[correlation_function_file_information.length() + 1];
	strcpy(correlationfunctionfileinformation, correlation_function_file_information.c_str());
	ifstream cffi(correlationfunctionfileinformation);

	while(text.substr(0,16)!="total protons : "){
		getline(cffi,text);
	}
	int total_protons=string_to_int(text.substr(16));	
	
	
//	int proton_sphere[10000];
//	int proton[10000];
//	string atomResidue[10000];
//	string residue[10000];
//	string atomName[10000];

//	string crrm_index_information[10000];
	
//	cffi.seekg(0);
	
//	while(text.substr(0,4).compare("list")!=0){
//		getline(cffi,text);	
//	}
	
//	size_t first=0,second=0,third=0,fourth=0,fifth=0;
	
//	for(int i=0; i<total_protons; i++){
//		getline(cffi,text);
		
//		crrm_index_information[i]=text.substr(text.find(" ")+1);
		
//		first=(size_t) text.find(" ");
//		second=(size_t) text.find(" ",first+1);
//		third=(size_t) text.find(" ",second+1);
//		fourth=(size_t) text.find(" ",third+1);
//		fifth=(size_t) text.find(" ",fourth+1);
		
//		proton_sphere[i]=string_to_int(text.substr(0,first-1));
//		proton[i]=string_to_int(text.substr(first+1,second-1));
//		atomResidue[i]=text.substr(second+1,third-1);
//		residue[i]=text.substr(third+1,fourth-1);
//		atomName[i]=text.substr(fourth+1);
//	}
	
	// input from correlation function information file - how made
	
	cffi.seekg(0);
	
	for(int i=0; i<5; i++){ 
		getline(cffi,text);
	}
	
	//        argv[3] the total time of the simulation in ns
	//        argv[4] the display time 
	//        argv[5] the sample_frequency 
	//        argv[6] maximum distance of spin pair
	//        argv[7] jpeg,eps,pdf,epslatex
	//        argv[8] all/notall
	//        argv[9] output directory 
	//        argv[10] residue name
	//        argv[11] residue number
	//        argv[12] atom name
	//        argv[13] centerRadius

	getline(cffi,text);
	int total_time=string_to_int(text);
	
	getline(cffi,text);
	int display_time=string_to_int(text);
	
	getline(cffi,text);
	int sample_frequency=string_to_int(text);
	
	for(int i=0; i<8; i++){
		getline(cffi,text);
	}

	getline(cffi,text);
	int intfile=string_to_int(text);


	// correlation file names 

	string correlation_file[10000];
	int proton1[10000];
	int proton2[10000];

	for(int i=0; i<10000; i++){
 	   proton1[i]=0;
	   proton2[i]=0;
	}
	
	while(text.substr(0,19).compare("spin pair in sphere")!=0){
		getline(cffi,text);
	}
	for(int i=0; i<3; i++){
		getline(cffi,text); 
	}

//	cout << "correlation files" << endl << endl;
	
	int file_count=0; 
	int first=0;
	
	while(text.substr(0,4).compare("list")!=0){
		getline(cffi,text);
		correlation_file[file_count]=text;	

		for(int i=0; i<3; i++){
			getline(cffi,text);
		}

//cout << text << endl;

		first=(size_t) text.find(" ");
                proton1[file_count]=string_to_int(text.substr(0,first));
		proton2[file_count]=string_to_int(text.substr(first+1));

//cout << first << endl; 

//cout << proton1[file_count] << " " << proton2[file_count] << endl;

		file_count=file_count+1;

		for(int i=0; i<4; i++){
		   getline(cffi,text); 
		}

	}


	double *fourier_input=new double[(int) ((double) intfile/sample_frequency*display_time/total_time)];
	double *fourier_transform=new double[(int) ((double) intfile/sample_frequency*display_time/total_time)];
				
// check

//	double interval = (double) 1/(total_time*.000000001)/1000000;

                        double interval = (double) 1/total_time; 

// cout << interval << endl;

                        interval = (double) interval*intfile; 

// cout << interval << endl;

                        interval = (double) interval/sample_frequency; 

// cout << interval << endl;

                        interval = (double) interval/ .000000001;

// cout << interval << endl;

                        interval = (double) interval/intfile/1000000;
	
// cout << interval << endl;

// cout << total_time << " " << intfile << " " << sample_frequency << endl;


	// output files listed in ccrm
			
	string temp_out_correlation[10000];
	int total_file=0;
	
	string temp_in_correlation_file;
	string temp_out_correlation_file;


	// complete relaxation rate matrix calculation

	double **complete_rr_matrix;
	complete_rr_matrix = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	complete_rr_matrix[i] = new double[total_protons];
	
	double **temp_matrix;
	temp_matrix = new double* [total_protons];
	for (int i=0; i < total_protons; i++)
	temp_matrix[i] = new double[total_protons];

	for (int i = 0; i < total_protons; i++){
		for (int i2 = 0; i2 < total_protons; i2++){
			complete_rr_matrix[i][i2] = 0;
			temp_matrix[i][i2] = 0;
		}
	}
 
	double frequency = magnetic_field * 1000000.0 * 42 * 3.14;
	double distance_normalize = (double) 1.897 * 100000 * 1.897 * 100000*(42.576 / 10.705)*(42.576 / 10.705)/4;
	int frequency1001 = (double) 1/1000000*frequency/interval/2/3.1415;
	int frequency1002 = (double) 1/1000000*2* frequency/interval/2/3.1415;
	double nanosecond = 1000000000.0;
	
	string temp_out_fourier_transform_file;

				FILE *gp;
				gp = popen(GNUPLOT, "w");
				if (gp == NULL){
					return 0;
				}


	// frame every 2 picoseconds
	double picosecond=(double) total_time/intfile*sample_frequency;			
	double angle_coefficient=(double) 2*3.14159265/intfile*sample_frequency;
	

	// output files

	for(int file=0; file<file_count; file++){

		temp_in_correlation_file=in_directory;
		temp_in_correlation_file.append("/internal_c-function_");
		temp_in_correlation_file.append(correlation_file[file]);
		temp_in_correlation_file.append(".txt");
		
		char *tempincorrelationfile = new char[temp_in_correlation_file.length() + 1];
		strcpy(tempincorrelationfile, temp_in_correlation_file.c_str());
		ifstream ticf(tempincorrelationfile);
		

		temp_out_correlation_file=out_directory;
		temp_out_correlation_file.append("/total_c-function_");
		temp_out_correlation_file.append(correlation_file[file]);
		temp_out_correlation_file.append("_");
		temp_out_correlation_file.append(tumbling);
		temp_out_correlation_file.append(".txt");
		
		temp_out_correlation[total_file]=temp_out_correlation_file;
		
		char *tempoutcorrelationfile = new char[temp_out_correlation_file.length() + 1];
		strcpy(tempoutcorrelationfile, temp_out_correlation_file.c_str());
		ofstream tocf(tempoutcorrelationfile);
	
		
		temp_out_fourier_transform_file = out_directory;
		temp_out_fourier_transform_file.append("/total_fourier_");
		temp_out_fourier_transform_file.append(correlation_file[file]);
		temp_out_fourier_transform_file.append("_");
		temp_out_fourier_transform_file.append(tumbling);
		temp_out_fourier_transform_file.append(".txt");
		
		temp_out_correlation[total_file+1]=temp_out_fourier_transform_file;
		
		char *tempoutfouriertransformfile = new char[temp_out_fourier_transform_file.length() + 1];
		strcpy(tempoutfouriertransformfile, temp_out_fourier_transform_file.c_str());
		ofstream toftf(tempoutfouriertransformfile);
		
		total_file=total_file+2;

		double correlation=string_to_double(tumbling);
		
		for(int intfile_i=0; intfile_i<(double) intfile/sample_frequency*display_time/total_time; intfile_i++){	

			// parse the second column of correlation file
			
			getline(ticf,text);
			size_t found = text.find(" ");	
			fourier_input[intfile_i]=string_to_double(text.substr(found+1));

			// fourier input array
			// artificial tumbling correlation output 
			
			fourier_input[intfile_i] = fourier_input[intfile_i]*exp(-(double) 1/correlation*intfile_i/intfile*display_time/sample_frequency);		
			tocf << (double) intfile_i/intfile*display_time/sample_frequency << " " << fourier_input[intfile_i] << endl;
			
		}
		tocf.close();
	
	ticf.close();

				fprintf(gp, " set terminal jpeg \n unset key \n set output '%s.jpeg' \n set xlabel 'nanoseconds' \n set title 'correlation function' \n plot '%s' \n ", tempoutcorrelationfile,tempoutcorrelationfile);

		// fourier transform  
		
		double angle=0;
		double fourier_total=0;
		
		for (int i0 = 0; i0 < (double) max_frequency/interval; i0++){
			fourier_total = 0;
			for (int i2 = 0; i2 < (double) intfile/sample_frequency*display_time/total_time/2; i2++){
				angle = i0*i2*angle_coefficient;
				fourier_total = fourier_total + (double) fourier_input[i2] * cosf(angle);
			}
			fourier_transform[i0]=(double) 2*fourier_total*picosecond;		

			toftf << (double) i0*interval << " " << fourier_transform[i0] << endl;			
		}
		toftf.close();	

//					if (font.compare("jpeg") == 0){
						fprintf(gp, " set terminal jpeg \n unset key \n set output '%s.jpeg' \n set xlabel 'MHz'\n  set title 'fourier transform' \n plot '%s' \n ",tempoutfouriertransformfile,tempoutfouriertransformfile);
//					}
/*
					if (font.compare("epslatex") == 0){
						fprintf(gp, " set terminal epslatex \n unset key \n set output '%s_fourier.epslatex' \n set xlabel 'MHz'\n  set title 'fourier transform' \n plot '%s' \n ",tempoutfouriertransformfile, tempoutfouriertransformfile);
					}
					if (font.compare("pdf") == 0){
						fprintf(gp, " set terminal pdf \n unset key \n set output '%s_fourier.pdf' \n set xlabel 'MHz'\n  set title 'fourier transform' \n plot '%s' \n ",tempoutfouriertransformfile, tempoutfouriertransformfile);
					}
					if (font.compare("eps") == 0){
						fprintf(gp, " set terminal eps \n unset key \n set output '%s_fourier.eps' \n set xlabel 'MHz'\n  set title 'fourier transform' \n plot '%s' \n ",tempoutfouriertransformfile, tempoutfouriertransformfile);
					}
*/


		temp_matrix[proton1[file]][proton2[file]] = (double) 2/5/nanosecond* distance_normalize* (6 * fourier_transform[(int) frequency1002] + 3 * fourier_transform[(int) frequency1001] + fourier_transform[0]);
		temp_matrix[proton2[file]][proton1[file]] = temp_matrix[proton1[file]][proton2[file]];
		complete_rr_matrix[proton1[file]][proton2[file]] = (double) 2/5/nanosecond* distance_normalize * (6 * fourier_transform[(int) frequency1002] - fourier_transform[0]);
		complete_rr_matrix[proton2[file]][proton1[file]] = complete_rr_matrix[proton1[file]][proton2[file]];

// cout << proton1[file] << " " << proton2[file] << endl;

// cout << temp_out_correlation[2*file] << endl;

	}
	
// cout << "test 402" << endl;


	double *diagonal=new double[total_protons];

	for (int i2 = 0; i2 < total_protons; i2++){
		diagonal[i2]=0;
	}
	for (int i20 = 0; i20 < total_protons; i20++){
            for (int i0=0; i0<total_protons; i0++){
                if(i20!=i0){
                  diagonal[i20]=diagonal[i20]+temp_matrix[i20][i0];
                }
            }
        }
	for (int i22 = 0; i22 < total_protons; i22++){
		complete_rr_matrix[i22][i22] = diagonal[i22];
	}


	// complete relaxation rate matrix file
	
	string complete_relaxation_rate_matrix = out_directory;
	complete_relaxation_rate_matrix.append("/complete_relaxation_rate_matrix_");
	complete_relaxation_rate_matrix.append(argv[2]);
	complete_relaxation_rate_matrix.append("_");
	complete_relaxation_rate_matrix.append(argv[6]);
	complete_relaxation_rate_matrix.append(".txt");
	
	char *completerelaxationratematrix = new char[complete_relaxation_rate_matrix.length() + 1];
	strcpy(completerelaxationratematrix, complete_relaxation_rate_matrix.c_str());
	ofstream crrm(completerelaxationratematrix);
		
	crrm << "total correlation function calculation inputs : " << endl << endl; 


	crrm << argv[1] << endl;	
	crrm << argv[2] << endl;
	crrm << argv[3] << endl;
	crrm << argv[4] << endl;
	crrm << argv[5] << endl;
	crrm << argv[6] << endl;
	crrm << endl;

	crrm << "  - from inputs of internal correlation function calculation : " << endl << endl;
	crrm << total_time << endl;
	crrm << display_time << endl;
	crrm << sample_frequency << endl;

	crrm << intfile << endl;
	crrm << total_protons << endl << endl;
	
	crrm << "information about indexing of complete relaxation rate matrix" << endl << endl;

	for(int file=0; file<file_count; file++){
		crrm << temp_out_correlation[2*file] << endl;
		crrm << temp_out_correlation[2*file+1] << endl;
                crrm << proton1[file] << " " << proton2[file] << endl;
		crrm << endl;
	}
	crrm << " " << endl;
	
	crrm << "complete relaxation matrix" << endl << endl;


	for (int ii = 0; ii < total_protons; ii++){
	    for (int ii2=0; ii2 < total_protons; ii2++){
                crrm << ii << " " << ii2 << " " << complete_rr_matrix[ii][ii2] << endl;
	    }
	}


	cffi.seekg(0);

	while(text.substr(0,7).compare("list of")!=0){
	     getline(cffi,text);
	}
	getline(cffi,text);

	crrm << "empty" << endl << endl;

        crrm << "list of protons in the complete relaxation rate matrix" << endl << endl;

	for(int i=0; i<2*total_protons; i++){
           getline(cffi,text);
           crrm << text << endl; 
        }


	cffi.close();

	crrm.close();	
	
	return 0;
}













