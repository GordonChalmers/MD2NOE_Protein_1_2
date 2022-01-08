

// Gordon Chalmers

// inputs -
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
//        argv[14] test fraction



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

double *totalinitial_thread = new double[100000];

#define NUM_THREADS     1

struct thread_data{
	
	double *unit_3_x;
	double *unit_3_y;
	double *unit_3_z;
	double *r_minus_3; 
	double *totalinitial;
	int intfile;
	int sample_frequency;
	int thread_id;
	int total_time;
	int display_time;

};

struct thread_data thread_data_array[NUM_THREADS];


void *correlation_function_thread(void *threadarg){
	
	struct thread_data *my_data;

	my_data = (struct thread_data *) threadarg;
	
	double *unit_3_x_thread = my_data -> unit_3_x;
	double *unit_3_y_thread = my_data -> unit_3_y;
	double *unit_3_z_thread = my_data -> unit_3_z;
	double *r_minus_3_thread = my_data -> r_minus_3;
	
	int thread_id_thread = my_data -> thread_id;
	int intfile_thread = my_data -> intfile;
	int sample_frequency_thread = my_data -> sample_frequency;
	
	int total_time_thread = my_data -> total_time;
	int display_time_thread = my_data -> display_time;
	
	double correlation=0;
	int total = 0;
	double cos_angle=0;	
	
	int file_trajectory_id_thread1=(double) ((double) display_time_thread/total_time_thread)*((double) thread_id_thread/NUM_THREADS)*intfile_thread;
	int file_trajectory_id_thread2=(double) ((double) display_time_thread/total_time_thread)*(thread_id_thread+1)/NUM_THREADS*intfile_thread;
	int file_trajectory_t0_id_thread=(double) ((double) display_time_thread/total_time_thread)*intfile_thread;
	
	cout << "thread " << thread_id_thread << " : frames " << intfile_thread << " : " << file_trajectory_id_thread1 << " " << file_trajectory_id_thread2 <<  endl;
	
	for(int file_trajectory_thread=file_trajectory_id_thread1; file_trajectory_thread<file_trajectory_id_thread2; file_trajectory_thread=file_trajectory_thread+sample_frequency_thread){

		correlation=0;
		total = 0;
		
		for (int file_trajectory_t0 = 0; file_trajectory_t0 < (intfile_thread - file_trajectory_thread); file_trajectory_t0 = file_trajectory_t0+sample_frequency_thread){			
			if(file_trajectory_t0 < file_trajectory_t0_id_thread){
				
				total++;
				cos_angle=(double) (unit_3_x_thread[file_trajectory_t0]*unit_3_x_thread[file_trajectory_thread+file_trajectory_t0]+unit_3_y_thread[file_trajectory_t0]*unit_3_y_thread[file_trajectory_thread+file_trajectory_t0]+unit_3_z_thread[file_trajectory_t0]*unit_3_z_thread[file_trajectory_thread+file_trajectory_t0]);
				correlation = correlation + (double) cos_angle*cos_angle - (double) r_minus_3_thread[file_trajectory_t0]*r_minus_3_thread[file_trajectory_thread+file_trajectory_t0];
				
			}
		}//file_trajectory_t0

		correlation = (double) correlation / total /2;
		
		// totalinitial is the c-function
		totalinitial_thread[(int) (double) file_trajectory_thread/sample_frequency_thread] = correlation;
		
	}
	pthread_exit(NULL);
}



double coordinate=0;
long double ad=0;
double display_time=0;
double total_time=0;
double inv_distance=0;
string text;
double total=0;
double distanceatom = 0;
double distanceatom_test=0;
int sample_frequency=0;
double nanosecond = 1000000000.0;
int mi=0, ri=0, ai=0, bi=0, natom=0;

int centerAtom=0;
long int intfile=0;
double centerRadius=0;

int thread_number=0;
int total_protons=0;

int unique_proton1=0;
int unique_proton2=0;

int proton_set[1000];


int main(int argc, char *argv[]) { /* standard main function */


	double total_time = string_to_double(argv[3]);
	double display_time = string_to_double(argv[4]);
	double sample_frequency = string_to_int(argv[5]);
	double max_distance = string_to_double(argv[6]); 
	string font = argv[7];
	string directory_out = argv[9];

	string center_Residue=argv[10];
	string center_residue=argv[11];
	string center_Atom=argv[12];
	double centerRadius=string_to_double(argv[13]);
	double test_fraction=string_to_double(argv[14]);
	
	pthread_t threads[NUM_THREADS];
	pthread_attr_t attr;

	/* Initialize and set thread detached attribute */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//	void *status;
	
	//	FILE *gp;
	fileset IF;

	// topology file

	/* Input and output file info
	contains file name and pointer
	this simplifies file error reporting */
	assembly A; 		/* an empty, unititialized assembly */
	//	molindex moli;		/* a molecule index ... m,r,a indices */
	molecule *mtmp;		/* Pointer to a molecule (for convenience) */
	residue *rtmp;		/* Pointer to a residue (for convenience) */
	atom *atmp;
	//      atom *a2tmp; 	/* Pointers to atom (for convenience) */
	/* It is always good to check the command line sanity */
	/* glylib contains many convenience functions, for example: */
	//if (argc<2){mywhineusage("Need input prmtop file name.");}

	/* copy the input file name to a convenient location */
	/* this is convenient, not efficient or necessary */
	IF.N = strdup(argv[1]); 	/* Copy name of input file */

	/* Now, load in the contents of the amber prmtop file */
	/* load_amber_prmtop requires the file not be open already */
	A = load_amber_prmtop(IF); 	/* initializes assembly, reads
				file, allocates memory,
				parses file, etc. */

	amber_prmtop *MyAP;

	// initialize the amber prmtop info
	MyAP = (amber_prmtop *)calloc(1, sizeof(amber_prmtop));
	amber_prmtop_init(MyAP);

	// load info as-is into sections
	read_amber_prmtop_asis(IF, MyAP);

	string * atomName=new string[10000];
	string * atomType=new string[10000];
	string * atomResidue=new string[10000];
	string * residue=new string[10000];
	
	// natom is # atoms

	natom = 0;

	if (A.nm < 1){ mywhine("No molecules found in prmtop file!"); }

	//	for (mi = 0; mi < A.nm; mi++) dprint_molecule(A.m[mi],100);

	//printf("nm is %d\n", A.nm);
	for (mi = 0; mi < A.nm; mi++){ /* for each molecule */
		mtmp = &A.m[mi][0]; /* set pointer for code-reading ease */
		//printf("\tnr is %d\n", A.m[mi][0].nr);
		for (ri = 0; ri < mtmp[0].nr; ri++){ /* each residue */
			//printf("\t\tri  is %d\n", ri);
			//printf("\t\tna from A is %d\n", A.m[mi][0].r[ri].na);
			//printf("\t\tna from mtmp is %d\n", mtmp[0].r[ri].na);
			rtmp = &mtmp[0].r[ri]; /* note ampersand*/

			//printf("\t\tna from rtmp is %d\n", rtmp[0].na);
			//printf("\t\tn from rtmp is %d\n", rtmp[0].n);
			
			for (ai = 0; ai < rtmp[0].na; ai++){ /* atom */
				atmp = &rtmp[0].a[ai];
				residue[natom] = int_to_string(rtmp[0].n);
				//residue[natom]=ri;
				//atomMass[natom] = (atmp[0]).m;
				atomName[natom] = (atmp[0]).N;
				atomType[natom] = (atmp[0]).T[0];
				//atomTypeNumber[natom] = (atmp[0]).t;
				//atomNumber[natom] = (atmp[0]).n;
				//atomCharge[natom] = (atmp[0]).ch[0];
				atomResidue[natom] = A.m[atmp[0].moli.m][0].r[atmp[0].moli.r].N;
				natom++;

			}
		}
	}

	/**
	if(argc<14){
		cout << "Not enough inputs, please add proton conditions" << endl;
		exit(1);
	}

	if(ceil(string_to_double(argv[3]))!=string_to_double(argv[3])){
		cout << "simulation time is not integer" << endl;
		exit(1);
	}

	if(ceil(string_to_double(argv[4]))!=string_to_double(argv[4])){
		cout << "Display time is not integer" << endl;
		exit(1);
	}

	if(ceil(string_to_double(argv[5]))!=string_to_double(argv[5])){
		cout << "Sample frequency is not integer" << endl;
		exit(1);
	}

	//if(ceil(string_to_double(argv[8]))!=string_to_double(argv[8])){
	//  cout << "Time of noe curve is not integer" << endl;
	//  exit(1);
	//}

	if(ceil(string_to_double(argv[9]))!=string_to_double(argv[9])){
		cout << "Units of noe curve is not integer" << endl;
		exit(1);
	}

	//if(argv[12]!="all"||argv[12]="notall"){ 
	//  cout << "All/notall of the calculation" << endl;
	//  exit(1);
	//}
**/
	
	// calculates # frames in trajectory
	
	// have to close and then re-open

	intfile=0;
	ifstream fileofcoordinate(argv[2]);
	getline(fileofcoordinate, text);
	while (getline(fileofcoordinate, text) > 0){
		intfile++;
	}
	intfile = floor((double) intfile / ceil((double)natom * 3 / 10));
	fileofcoordinate.close();
	
	cout << intfile << " frames in trajectory " << endl;
	
	int intfile_1_100=(int) intfile*test_fraction;
	
	cout << "test : " << intfile_1_100 << " frames" << endl;
	
	
	double * totalinitial = new double[intfile];
	double * unit_3_x = new double[intfile];
	double * unit_3_y = new double[intfile];
	double * unit_3_z = new double[intfile];
	double * r_minus_3 = new double[intfile];
	
	//	double totalinitial[100000];
	//	double unit_3_x[100000];
	//	double unit_3_y[100000];
	//	double unit_3_z[100000];
	//	double r_minus_3[100000];
	
	// initialize
	
	for(int i=0; i<intfile; i++){
		unit_3_x[i]=0;
		unit_3_y[i]=0;
		unit_3_z[i]=0;
		r_minus_3[i]=0;	
		totalinitial[i]=0;
	}
	
	
	// distance information

	string file_of_distance = argv[9];
	file_of_distance.append("/file_distance.txt");
	char *fileofdistance = new char[file_of_distance.length() + 1];
	strcpy(fileofdistance, file_of_distance.c_str());
	ofstream fod(fileofdistance);


	// correlation function information file 

	string file_of_cross_rate_spin_pairs = argv[9];
	file_of_cross_rate_spin_pairs.append("/internal_correlation_function_");
	file_of_cross_rate_spin_pairs.append(argv[9]);
	file_of_cross_rate_spin_pairs.append(".txt");
	char *fileofcrossratespinpairs = new char[file_of_cross_rate_spin_pairs.length() + 1];
	strcpy(fileofcrossratespinpairs, file_of_cross_rate_spin_pairs.c_str());
	ofstream fcf(fileofcrossratespinpairs);

	// parameters of calculation 
	
	fcf << "internal correlation function calculations" << endl; 
	fcf << "  - from the inputs of correlation_function.cc : " << endl << endl;

	fcf << argv[1] << endl;	
	fcf << argv[2] << endl;
	fcf << argv[3] << endl;
	fcf << argv[4] << endl;
	fcf << argv[5] << endl;
	fcf << argv[6] << endl;
	fcf << argv[7] << endl;
	fcf << argv[8] << endl;
	fcf << argv[9] << endl;
	fcf << argv[10] << endl;
	fcf << argv[11] << endl;
	fcf << argv[12] << endl;
	fcf << argv[13] << endl;
	fcf << intfile << endl;
	fcf << intfile_1_100 << endl;
	fcf << endl;

	
	// find center atom_number 
	
	for(int i=0; i<natom; i++){
		if(atomResidue[i].compare(center_Residue)==0){
			if(residue[i].compare(center_residue)==0){
				if(atomName[i].compare(center_Atom)==0){
					centerAtom=i;
				}
			}
		}
	}
	
	fileofcoordinate.open(argv[2]);

	
	// average center coordinates

	int total=0;

	double center_coordinate_x=0, temp_center_x=0, temp_proton_x=0;
	double center_coordinate_y=0, temp_center_y=0, temp_proton_y=0;
	double center_coordinate_z=0, temp_center_z=0, temp_proton_z=0;

	int initial_point=(3 * natom) % 10;
	long int initial_point_10=81 * ceil((double)natom * 3 / 10) - (10 - (3 * natom) % 10) * 8;
	long int initial_point_0=81 * ceil((double)natom * 3 / 10);

	int center_point_add_16=(3 * centerAtom * 8) % 80 + 16;
	int center_point_add_24=(3 * centerAtom * 8) % 80 + 24;

	double center_point_10=(double) (3 * centerAtom - (3 * centerAtom) % 10) / 10;
	int center_point_0=(3 * centerAtom * 8) % 80;
	int center_point_8=(3 * centerAtom * 8) % 80+8;
	int center_point_16=(3 * centerAtom * 8) % 80+16;

	int sphere_proton_add_16=0;
	int sphere_proton_add_24=0;

	double sphere_proton_10=0;
	int sphere_proton_0=0;
	int sphere_proton_8=0;
	int sphere_proton_16=0;

	double temp=0;
	total=0;
	
	for(int file_trajectory=0; file_trajectory<intfile; file_trajectory=file_trajectory+intfile_1_100){
		total++;
		fileofcoordinate.seekg(0);
		getline(fileofcoordinate, text);
		if (initial_point != 0){
			fileofcoordinate.seekg(file_trajectory*(initial_point_10-1), ios_base::cur);
		}
		if (initial_point == 0){
			fileofcoordinate.seekg(file_trajectory*(initial_point_0-1), ios_base::cur);
		}
		fileofcoordinate.seekg(file_trajectory, ios_base::cur);
		for (int unused = 0; unused < center_point_10; unused++){
			getline(fileofcoordinate, text);
		}
		getline(fileofcoordinate, text);
		temp_center_x = string_to_double(text.substr((size_t) center_point_0, 8))+temp_center_x;
		if (center_point_add_16 != 88){
			temp_center_y = string_to_double(text.substr((size_t) center_point_8, 8))+temp_center_y;
			if (center_point_add_24 <= 80){
				temp_center_z = string_to_double(text.substr((size_t) center_point_16, 8))+temp_center_z;	
			}
			else{
				getline(fileofcoordinate, text);
				temp_center_z= string_to_double(text.substr((size_t) 0, 8))+temp_center_z;
			}
		}
		if (center_point_add_16 == 88){
			getline(fileofcoordinate, text);
			temp_center_y = string_to_double(text.substr((size_t)0, 8))+temp_center_y;
			temp_center_z = string_to_double(text.substr((size_t)8, 8))+temp_center_z;
		}
	}
	center_coordinate_x=temp_center_x/total;
	center_coordinate_y=temp_center_y/total;
	center_coordinate_z=temp_center_z/total;

	cout << "center coordinates : " << center_coordinate_x << " " << center_coordinate_y << " " << center_coordinate_z << endl;
	cout << "center proton : " << centerAtom << " " << atomResidue[centerAtom] << " " << residue[centerAtom] << " " << atomName[centerAtom] << endl;
	
	fcf << "center proton : " << centerAtom << " " << atomResidue[centerAtom] << " " << residue[centerAtom] << " " << atomName[centerAtom] << endl;
	fcf << "center coordinates : " << center_coordinate_x << " " << center_coordinate_y << " " << center_coordinate_z << endl;
	
	
	// find all protons in sphere - using centerRadius - test

	temp_center_x=0;
	temp_center_y=0;
	temp_center_z=0;

	temp_proton_x=0;
	temp_proton_y=0;
	temp_proton_z=0;

	int total_protons_in_sphere_test=0;
	
	int * proton_in_sphere_test=new int[natom];
	int * proton=new int[natom];

	for(int i=0; i<natom; i++){
		proton_in_sphere_test[i]=0;
		proton[i]=0;
	}

	for (int sphere_proton = 0; sphere_proton < natom; sphere_proton++){	
		if (atomType[sphere_proton].compare("H") == 0){
			
//			cout << atomType[sphere_proton] << endl;
			
			sphere_proton_add_16=(3 * sphere_proton * 8) % 80 + 16;
			sphere_proton_add_24=(3 * sphere_proton * 8) % 80 + 24;

			sphere_proton_10=(double) (3 * sphere_proton - (3 * sphere_proton) % 10) / 10;
			sphere_proton_0=(3 * sphere_proton * 8) % 80;
			sphere_proton_8=(3 * sphere_proton * 8) % 80+8;
			sphere_proton_16=(3 * sphere_proton * 8) % 80+16;

			temp=0;	
			total=0;

			for(int file_trajectory=0; file_trajectory<intfile; file_trajectory=file_trajectory+intfile_1_100){
				
				// center coordinate at file_trajectory frame
				
				total++;
				fileofcoordinate.seekg(0);
				getline(fileofcoordinate, text);
				if (initial_point != 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_10-1), ios_base::cur);
				}
				if (initial_point == 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_0-1), ios_base::cur);
				}
				fileofcoordinate.seekg(file_trajectory, ios_base::cur);
				for (int unused = 0; unused < center_point_10; unused++){
					getline(fileofcoordinate, text);
				}
				getline(fileofcoordinate, text);
				temp_center_x = string_to_double(text.substr((size_t) center_point_0, 8));
				if (center_point_add_16 != 88){
					temp_center_y = string_to_double(text.substr((size_t) center_point_8, 8));
					if (center_point_add_24 <= 80){
						temp_center_z = string_to_double(text.substr((size_t) center_point_16, 8));
					}
					else{
						getline(fileofcoordinate, text);
						temp_center_z= string_to_double(text.substr((size_t) 0, 8));
					}
				}
				if (center_point_add_16 == 88){
					getline(fileofcoordinate, text);
					temp_center_y = string_to_double(text.substr((size_t)0, 8));
					temp_center_z = string_to_double(text.substr((size_t)8, 8));
				}			
				
				// sphere proton at file_trajectory frame
				
				fileofcoordinate.seekg(0);
				getline(fileofcoordinate, text);
				if (initial_point != 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_10-1), ios_base::cur);
				}
				if (initial_point == 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_0-1), ios_base::cur);
				}
				fileofcoordinate.seekg(file_trajectory, ios_base::cur);
				for (int unused = 0; unused < sphere_proton_10; unused++){
					getline(fileofcoordinate, text);
				}
				getline(fileofcoordinate, text);
				temp_proton_x = string_to_double(text.substr((size_t) sphere_proton_0, 8));
				if (sphere_proton_add_16 != 88){
					temp_proton_y = string_to_double(text.substr((size_t) sphere_proton_8, 8));
					if (sphere_proton_add_24 <= 80){
						temp_proton_z = string_to_double(text.substr((size_t) sphere_proton_16, 8));
					}
					else{
						getline(fileofcoordinate, text);
						temp_proton_z= string_to_double(text.substr((size_t)0, 8));
					}
				}
				if (sphere_proton_add_16 == 88){
					getline(fileofcoordinate, text);
					temp_proton_y = string_to_double(text.substr((size_t)0, 8));
					temp_proton_z = string_to_double(text.substr((size_t)8, 8));
				}
				temp=temp+(double) sqrt((temp_center_x-temp_proton_x)*(temp_center_x-temp_proton_x)+(temp_center_y-temp_proton_y)*(temp_center_y-temp_proton_y)+(temp_center_z-temp_proton_z)*(temp_center_z-temp_proton_z));
			}
			temp=temp/total;

			// sphere_proton is in sphere of centerRadius about centerAtom
			// sphere_proton indexed by topology file
			
			if(temp<=centerRadius){
				proton_in_sphere_test[sphere_proton]=temp;
				proton[total_protons_in_sphere_test]=sphere_proton;
//				cout << atomType[sphere_proton] << " " << atomType[proton[total_protons_in_sphere_test]] << endl;
				total_protons_in_sphere_test++;
			}
		}
	}
	
	cout << total_protons_in_sphere_test << " protons in the sphere of radius test of " << centerRadius << endl << endl;

	fcf << "sphere of radius : " << centerRadius << endl;	
	fcf << total_protons_in_sphere_test << " protons in the sphere of radius test of " << centerRadius << endl << endl;
	
		
	// spin pairs < max_distance in sphere - test
	
	fcf << "spin pair - test - in sphere, <r> <=: " << max_distance << endl << endl;

	int proton1=0;
	int proton2=0;

	int proton1_add_16=0;
	int proton1_add_24=0;

	int proton1_10=0;
	int proton1_0=0;
	int proton1_8=0;
	int proton1_16=0;

	int proton2_add_16=0;
	int proton2_add_24=0;

	int proton2_10=0;
	int proton2_0=0;
	int proton2_8=0;
	int proton2_16=0;

	double temp_proton1_x=0; 
	double temp_proton1_y=0;
	double temp_proton1_z=0;
	double temp_proton2_x=0; 
	double temp_proton2_y=0;
	double temp_proton2_z=0;

	int total_possible_spin_pairs=0;
	
	//	int * spin_pair_proton1=new int[total_protons_in_sphere];
	//	int * spin_pair_proton2=new int[total_protons_in_sphere];
	//	double * test_spin_pair=new double[total_protons_in_sphere*total_protons_in_sphere];

	double test_spin_pair[10000];
	
	int possible_spin_pair_proton1[10000];
	int possible_spin_pair_proton2[10000];	
	
	int possible_proton1_sphere_spin_pair[10000];
	int possible_proton2_sphere_spin_pair[10000];
	
//	int spin_pair_proton1[10000];
//	int spin_pair_proton2[10000];	
	
	int proton1_sphere_spin_pair[10000];
	int proton2_sphere_spin_pair[10000];
	
		for(int i=0; i<10000; i++){
//			spin_pair_proton1[i]=0;
//			spin_pair_proton2[i]=0;
			proton1_sphere_spin_pair[i]=0;
			proton2_sphere_spin_pair[i]=0;
			
			possible_spin_pair_proton1[i]=0;
			possible_spin_pair_proton2[i]=0;
			possible_proton1_sphere_spin_pair[i]=0;
			possible_proton2_sphere_spin_pair[i]=0;
			
			test_spin_pair[i]=0;
		}

	for (proton1 = 0; proton1 < total_protons_in_sphere_test; proton1++){		
		for (proton2 = proton1+1; proton2 < total_protons_in_sphere_test; proton2++){
		
			proton1_add_16=(3 * proton[proton1] * 8) % 80 + 16;
			proton1_add_24=(3 * proton[proton1] * 8) % 80 + 24;

			proton1_10=(double) (3 * proton[proton1] - (3 * proton[proton1]) % 10) / 10;
			proton1_0=(3 * proton[proton1] * 8) % 80;
			proton1_8=(3 * proton[proton1] * 8) % 80+8;
			proton1_16=(3 * proton[proton1] * 8) % 80+16;

			proton2_add_16=(3 * proton[proton2] * 8) % 80 + 16;
			proton2_add_24=(3 * proton[proton2] * 8) % 80 + 24;

			proton2_10=(double) (3 * proton[proton2] - (3 * proton[proton2]) % 10) / 10;
			proton2_0=(3 * proton[proton2] * 8) % 80;
			proton2_8=(3 * proton[proton2] * 8) % 80+8;
			proton2_16=(3 * proton[proton2] * 8) % 80+16;

			temp=0;	
			total=0;
			
			for(int file_trajectory=0; file_trajectory<intfile; file_trajectory=file_trajectory+intfile_1_100){

				// proton1 coordinate at frame file_trajectory
				
				total++;
				fileofcoordinate.seekg(0);
				getline(fileofcoordinate, text);
				if (initial_point != 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_10-1), ios_base::cur);
				}
				if (initial_point == 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_0-1), ios_base::cur);
				}
				fileofcoordinate.seekg(file_trajectory, ios_base::cur);
				for (int unused = 0; unused < proton1_10; unused++){
					getline(fileofcoordinate, text);
				}
				getline(fileofcoordinate, text);
				temp_proton1_x = string_to_double(text.substr((size_t) proton1_0, 8));
				if (proton1_add_16 != 88){
					temp_proton1_y = string_to_double(text.substr((size_t) proton1_8, 8));
					if (proton1_add_24 <= 80){
						temp_proton1_z = string_to_double(text.substr((size_t) proton1_16, 8));
					}
					else{
						getline(fileofcoordinate, text);
						temp_proton1_z= string_to_double(text.substr((size_t)0, 8));
					}
				}
				if (proton1_add_16 == 88){
					getline(fileofcoordinate, text);
					temp_proton1_y = string_to_double(text.substr((size_t)0, 8));
					temp_proton1_z = string_to_double(text.substr((size_t)8, 8));
				}	
				
				// proton2 coordinate at frame file_trajectory

				fileofcoordinate.seekg(0);
				getline(fileofcoordinate, text);
				if (initial_point != 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_10-1), ios_base::cur);
				}
				if (initial_point == 0){
					fileofcoordinate.seekg(file_trajectory*(initial_point_0-1), ios_base::cur);
				}
				fileofcoordinate.seekg(file_trajectory, ios_base::cur);
				for (int unused = 0; unused < proton2_10; unused++){
					getline(fileofcoordinate, text);
				}
				getline(fileofcoordinate, text);
				temp_proton2_x = string_to_double(text.substr((size_t) proton2_0, 8));
				if (proton2_add_16 != 88){
					temp_proton2_y = string_to_double(text.substr((size_t) proton2_8, 8));
					if (proton2_add_24 <= 80){
						temp_proton2_z = string_to_double(text.substr((size_t) proton2_16, 8));
					}
					else{
						getline(fileofcoordinate, text);
						temp_proton2_z= string_to_double(text.substr((size_t)0, 8));
					}
				}
				if (proton2_add_16 == 88){
					getline(fileofcoordinate, text);
					temp_proton2_y = string_to_double(text.substr((size_t)0, 8));
					temp_proton2_z = string_to_double(text.substr((size_t)8, 8));
				}							
				temp=temp+sqrt((temp_proton1_x-temp_proton2_x)*(temp_proton1_x-temp_proton2_x)+(temp_proton1_y-temp_proton2_y)*(temp_proton1_y-temp_proton2_y)+(temp_proton1_z-temp_proton2_z)*(temp_proton1_z-temp_proton2_z));
			}
			temp=temp/total;
			
			if(temp<=max_distance){
				
				possible_spin_pair_proton1[total_possible_spin_pairs]=proton[proton1];
				possible_spin_pair_proton2[total_possible_spin_pairs]=proton[proton2];
				
				possible_proton1_sphere_spin_pair[total_possible_spin_pairs]=proton1;
				possible_proton2_sphere_spin_pair[total_possible_spin_pairs]=proton2;

				test_spin_pair[total_possible_spin_pairs]=temp;

				cout << proton1 << " " << proton2 << endl;
				cout << total_possible_spin_pairs << " " << possible_spin_pair_proton1[total_possible_spin_pairs] << " " << possible_spin_pair_proton2[total_possible_spin_pairs] << endl;
                // cout << atomType[possible_spin_pair_proton1[total_possible_spin_pairs]] << " " << atomType[possible_spin_pair_proton2[total_possible_spin_pairs]] << endl;				
				cout << test_spin_pair[total_possible_spin_pairs] << endl << endl;

				fcf << "protons : " << proton1 << " " << proton2 << endl;
				fcf << total_possible_spin_pairs << " " << possible_spin_pair_proton1[total_possible_spin_pairs] << " " << possible_spin_pair_proton2[total_possible_spin_pairs] << endl;
				fcf << "test <r> distance : " << test_spin_pair[total_possible_spin_pairs] << endl << endl;
				
				total_possible_spin_pairs++;
			}
		}
	}
	
	cout << total_possible_spin_pairs << " spin pairs in the sphere - test" << endl;
	cout <<	"   - from the " << total_protons_in_sphere_test << " in sphere" << endl << endl << endl;

	fcf << "possible spin pairs in the sphere - test - : " << total_possible_spin_pairs << endl; 
	fcf << "  - from the " << total_protons_in_sphere_test << " in sphere : <r> <= " << max_distance << endl << endl;


	
	
	fcf << "spin pair in sphere from test if - all frames" << endl;
	fcf << "  - <1/r^6>^(-1/6) <= : " << max_distance << endl << endl;

	int total_spin_pairs=0;
	
	// correlation function calculation
	
	for (int possible_spin_pair=0; possible_spin_pair<total_possible_spin_pairs; possible_spin_pair++){

		proton1=possible_spin_pair_proton1[possible_spin_pair];
		proton2=possible_spin_pair_proton2[possible_spin_pair];

		proton1_add_16=(3 * proton1 * 8) % 80 + 16;
		proton1_add_24=(3 * proton1 * 8) % 80 + 24;

		proton1_10=(double) (3 * proton1 - (3 * proton1) % 10) / 10;
		proton1_0=(3 * proton1 * 8) % 80;
		proton1_8=(3 * proton1 * 8) % 80+8;
		proton1_16=(3 * proton1 * 8) % 80+16;

		proton2_add_16=(3 * proton2 * 8) % 80 + 16;
		proton2_add_24=(3 * proton2 * 8) % 80 + 24;

		proton2_10=(double) (3 * proton2 - (3 * proton2) % 10) / 10;
		proton2_0=(3 * proton2 * 8) % 80;
		proton2_8=(3 * proton2 * 8) % 80+8;
		proton2_16=(3 * proton2 * 8) % 80+16;

		total=0;
		double r=0;
		double r_1_6=0;

		temp=0;
		distanceatom=0;
		inv_distance=0;
		int rc=0;
		
		for(int file_trajectory=0; file_trajectory<intfile; file_trajectory=file_trajectory+sample_frequency){
			
			// proton1 coordinate at frame file_trajectory
			
			total++;
			fileofcoordinate.seekg(0);
			getline(fileofcoordinate, text);
			if (initial_point != 0){
				fileofcoordinate.seekg(file_trajectory*(initial_point_10-1), ios_base::cur);
			}
			if (initial_point == 0){
				fileofcoordinate.seekg(file_trajectory*(initial_point_0-1), ios_base::cur);
			}
			fileofcoordinate.seekg(file_trajectory, ios_base::cur);
			for (int unused = 0; unused < proton1_10; unused++){
				getline(fileofcoordinate, text);
			}
			getline(fileofcoordinate, text);
			temp_proton1_x = string_to_double(text.substr((size_t) proton1_0, 8));
			if (proton1_add_16 != 88){
				temp_proton1_y = string_to_double(text.substr((size_t) proton1_8, 8));
				if (proton1_add_24 <= 80){
					temp_proton1_z = string_to_double(text.substr((size_t) proton1_16, 8));
				}
				else{
					getline(fileofcoordinate, text);
					temp_proton1_z= string_to_double(text.substr((size_t)0, 8));
				}
			}
			if (proton1_add_16 == 88){
				getline(fileofcoordinate, text);
				temp_proton1_y = string_to_double(text.substr((size_t)0, 8));
				temp_proton1_z = string_to_double(text.substr((size_t)8, 8));
			}
			
			// proton2 coordinate at frame file_trajectory

			fileofcoordinate.seekg(0);
			getline(fileofcoordinate, text);
			if (initial_point != 0){
				fileofcoordinate.seekg(file_trajectory*(initial_point_10-1), ios_base::cur);
			}
			if (initial_point == 0){
				fileofcoordinate.seekg(file_trajectory*(initial_point_0-1), ios_base::cur);
			}
			fileofcoordinate.seekg(file_trajectory, ios_base::cur);
			for (int unused = 0; unused < proton2_10; unused++){
				getline(fileofcoordinate, text);
			}
			getline(fileofcoordinate, text);
			temp_proton2_x = string_to_double(text.substr((size_t) proton2_0, 8));
			if (proton2_add_16 != 88){
				temp_proton2_y = string_to_double(text.substr((size_t) proton2_8, 8));
				if (proton2_add_24 <= 80){
					temp_proton2_z = string_to_double(text.substr((size_t) proton2_16, 8));
				}
				else{
					getline(fileofcoordinate, text);
					temp_proton2_z= string_to_double(text.substr((size_t)0, 8));
				}
			}
			if (proton2_add_16 == 88){
				getline(fileofcoordinate, text);
				temp_proton2_y = string_to_double(text.substr((size_t)0, 8));
				temp_proton2_z = string_to_double(text.substr((size_t)8, 8));
			}	
			
			r=(double) pow((double) (temp_proton2_x-temp_proton1_x)*(temp_proton2_x-temp_proton1_x)+(temp_proton2_y-temp_proton1_y)*(temp_proton2_y-temp_proton1_y)+(temp_proton2_z-temp_proton1_z)*(temp_proton2_z-temp_proton1_z),.5);
			r_1_6=(double) pow((double) (temp_proton2_x-temp_proton1_x)*(temp_proton2_x-temp_proton1_x)+(temp_proton2_y-temp_proton1_y)*(temp_proton2_y-temp_proton1_y)+(temp_proton2_z-temp_proton1_z)*(temp_proton2_z-temp_proton1_z),-3.0);

			unit_3_x[file_trajectory]=(double) pow(3,.25)*(temp_proton2_x/pow(r,2.5) - temp_proton1_x/pow(r,2.5));		
			unit_3_y[file_trajectory]=(double) pow(3,.25)*(temp_proton2_y/pow(r,2.5) - temp_proton1_y/pow(r,2.5));
			unit_3_z[file_trajectory]=(double) pow(3,.25)*(temp_proton2_z/pow(r,2.5) - temp_proton1_z/pow(r,2.5));
			r_minus_3[file_trajectory]=(double) pow(r,-3.0);
			
			distanceatom=distanceatom+r;		   
			inv_distance = (double) inv_distance + r_1_6;
		}//file_trajectory		
		
		distanceatom = distanceatom / total;
		inv_distance = inv_distance / total;
		
		
		// (1/r^6)^(-1/6) < max_distance - the correlation function is calculated

		if (pow(inv_distance, (double)-1 / 6.00001) < max_distance){
			
//			spin_pair_proton1[total_spin_pairs]=proton1;
//			spin_pair_proton2[total_spin_pairs]=proton2;

			proton1_sphere_spin_pair[total_spin_pairs]=possible_proton1_sphere_spin_pair[possible_spin_pair];
			proton2_sphere_spin_pair[total_spin_pairs]=possible_proton2_sphere_spin_pair[possible_spin_pair];
			
			// noefile - the c-function 
			
			string noe_file = directory_out;
			noe_file.append("/internal_c-function_");

			noe_file.append(atomResidue[proton1]);
			noe_file.append("_");
			noe_file.append(residue[proton1]);
			noe_file.append("_");
			noe_file.append(atomName[proton1]);
			noe_file.append("-");
			noe_file.append(atomResidue[proton2]);
			noe_file.append("_");
			noe_file.append(residue[proton2]);
			noe_file.append("_");
			noe_file.append(atomName[proton2]);
			noe_file.append(".txt");

			char *noefile = new char[noe_file.length() + 1];
			strcpy(noefile, noe_file.c_str());
			ofstream file(noefile);

				
		cout << "spin pair in sphere : " << total_spin_pairs << endl;	
		cout << atomResidue[proton1] << ":" << residue[proton1] << ":" << atomName[proton1] << endl;					
		cout << atomResidue[proton2] << ":" << residue[proton2] << ":" << atomName[proton2] << endl;
	
		cout << "test : " << test_spin_pair[possible_spin_pair] << flush << endl; 
		cout << proton1_sphere_spin_pair[total_spin_pairs] << " " << proton1 << " : " << proton2_sphere_spin_pair[total_spin_pairs] << " " << proton2 << endl;
		cout << "all frames with sampling: <r> " << distanceatom << " " << "<1/r^6>^(-1/6) " << pow(inv_distance, (double)-1 / 6.00001) << endl;
		cout << endl;
			

proton_set[proton1_sphere_spin_pair[total_spin_pairs]]=proton1;
proton_set[proton2_sphere_spin_pair[total_spin_pairs]]=proton2;
		
if(proton1_sphere_spin_pair[total_spin_pairs]>total_protons){
  total_protons=proton1_sphere_spin_pair[total_spin_pairs];
}

if(proton2_sphere_spin_pair[total_spin_pairs]>total_protons){
  total_protons=proton2_sphere_spin_pair[total_spin_pairs];
}

cout << "total protons " << total_protons << endl;

	
			fcf << "spin pair in sphere " << total_spin_pairs << endl;
			fcf << atomResidue[proton1] << "_" << residue[proton1] << "_" << atomName[proton1] << "-" << atomResidue[proton2] << "_" << residue[proton2] << "_" << atomName[proton2] << endl;
			fcf << "internal_c-function_" << atomResidue[proton1] << "_" << residue[proton1] << "_" << atomName[proton1] << "-" << atomResidue[proton2] << "_" << residue[proton2] << "_" << atomName[proton2] << ".txt" << endl;
			fcf << "test : " << proton1 << " " << proton2 << " " << test_spin_pair[possible_spin_pair] << endl; 
			fcf << proton1_sphere_spin_pair[total_spin_pairs] << " " << proton2_sphere_spin_pair[total_spin_pairs]  << endl;
			fcf << proton1 << " " << proton2 << endl;
			fcf << "all frames with sampling: <r> " << distanceatom << " " << "<1/r^6>^(-1/6) " << pow(inv_distance, (double)-1 / 6.00001) << endl;
			fcf << endl;

			fod << atomResidue[proton1] << "," << residue[proton1] << "," << atomName[proton1] << "," << atomResidue[proton2] << "," << residue[proton2] << "," << atomName[proton2] << endl;
			fod << distanceatom << "," << pow(inv_distance,(double) -1/6.0) << endl;

			total_spin_pairs++;
			
			// threads

			//int intfile;
			//int  thread_id;
			//double *unit_3_x;
			//double *unit_3_y;
			//double *unit_3_z;
			//double *r_minus_3;

			// this for all arrays;	
			
			for(int tda=0; tda<NUM_THREADS; tda++){
				thread_data_array[tda].intfile=intfile;
				thread_data_array[tda].total_time=total_time;
				thread_data_array[tda].display_time=display_time;
				thread_data_array[tda].sample_frequency=sample_frequency;
				thread_data_array[tda].unit_3_x=unit_3_x;
				thread_data_array[tda].unit_3_y=unit_3_y;
				thread_data_array[tda].unit_3_z=unit_3_z;
				thread_data_array[tda].r_minus_3=r_minus_3;
			}
			thread_number=0;
			rc=0;

			for(int thread_number=0; thread_number<NUM_THREADS; thread_number++){	
				thread_data_array[thread_number].thread_id = thread_number;
				rc = pthread_create(&threads[thread_number], NULL, correlation_function_thread, (void *)&thread_data_array[thread_number]);
				if (rc){
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
				
				cout << "thread = " << thread_number << endl;

			}

			pthread_attr_destroy(&attr);
			for(int thread_number=0; thread_number<NUM_THREADS; thread_number++){
				rc = pthread_join(threads[thread_number],NULL);
			}
			
			cout << "done with internal correlation function calculation" << endl << endl;
			
			for (int ii = 0; ii < (display_time/total_time)*(intfile/sample_frequency); ii++){
				file << ii*total_time/intfile*sample_frequency << " " << (double) totalinitial_thread[ii] << endl;
			}
			file.close();
			
			
			// figures
			
			char *atomName1 = new char[atomName[proton1].length() + 1];
			char *atomName2 = new char[atomName[proton2].length() + 1];
			char *atomResidue1 = new char[atomResidue[proton1].length() + 1];
			char *atomResidue2 = new char[atomResidue[proton2].length() + 1];
			char *residue1 = new char[residue[proton1].length()+1];
			char *residue2 = new char[residue[proton2].length()+1];
			
			strcpy(atomName1, atomName[proton1].c_str());
			strcpy(atomName2, atomName[proton2].c_str());
			strcpy(atomResidue1, atomResidue[proton1].c_str());
			strcpy(atomResidue2, atomResidue[proton2].c_str());
			strcpy(residue1, residue[proton1].c_str());
			strcpy(residue2, residue[proton2].c_str());
			
			inv_distance=pow(inv_distance,-1/6);
			
			FILE *gp;
			gp = popen(GNUPLOT, "w");
			if (gp == NULL){
				return 0;
			}
			fprintf(gp, " set terminal jpeg \n unset key \n set output '%s/internal_correlation_%s_%s_%s-%s_%s_%s.jpeg' \n set xlabel 'nanoseconds' \n set title 'correlation function - atom %s %s:%s  atom %s %s:%s  %.2f %.2f \n plot '%s' \n ", argv[9], atomResidue1, residue1, atomName1, atomResidue2, residue2, atomName2, atomResidue1, residue1, atomName1, atomResidue2, residue2, atomName2, distanceatom, inv_distance, noefile);
			fclose(gp);
		}
	}
	
	
	fcf << "list of all protons in complete relaxation rate matrix:" << endl << endl;
	
total_protons++;
	
	for(int i=0; i<total_protons; i++){
		fcf << i << " " << proton_set[i] << endl;
		fcf << atomName[proton_set[i]] << " " << atomResidue[proton_set[i]] << " " << residue[proton_set[i]] << endl;
	}
	fcf << endl << "total protons : " << total_protons << endl;
	
	fileofcoordinate.close();
	fcf.close();
	fod.close();
	
	cout << endl << endl << "done" << endl;
	
	return 0;
}
