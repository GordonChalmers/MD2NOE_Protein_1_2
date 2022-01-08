// Gordon Chalmers

// inputs -

//        argv[3] the samplefrequency 
//        argv[4] out directory

//        argv[5] Residue_name_1
//        argv[6] residue1
//        argv[7] atom1
//        argv[8] Residue_name_2 
//        argv[9] residue2 
//        argv[10] atom2


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
#define GNUPLOT "gnuplot -persist"



double string_to_double(const std::string& s)
{
	std::istringstream i(s);
	long double x;
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



string text;
int intfile=0;
int ifirst=0;
int isecond=0;
int natom=0;
long double coordinate=0;
long int ad=0;
long double coordinate_spin_pair[2];
long double coordinateAtom[2];
long double coordinate_dx=0;
long double coordinate_dy=0;
long double coordinate_dz=0;
long double coordinate_ave_dx=0;
long double coordinate_ave_dy=0;
long double coordinate_ave_dz=0;
long double coordinate_spin_1_x=0;
long double coordinate_spin_1_y=0;
long double coordinate_spin_1_z=0;
long double coordinate_spin_2_x=0;
long double coordinate_spin_2_y=0;
long double coordinate_spin_2_z=0;
long double distance_normalized_spin_pair=0;
long double total=0;
long double cosine_angle=0;
long double cosine_angle_x=0;
long double cosine_angle_y=0;
long double cosine_angle_z=0;
long double cosine_angle_x_B=0;
long double cosine_angle_y_B=0;
long double cosine_angle_z_B=0;
long double average_distance=0;
long double bond_distance=0;

long double coordinate_diff_x=0;
long double coordinate_diff_y=0;
long double coordinate_diff_z=0;


long double total_Sxx=0;
long double total_Sxy=0;
long double total_Sxz=0;
long double total_Syz=0;
long double total_Syy=0;

long double total_Sxx_B=0;
long double total_Sxy_B=0;
long double total_Sxz_B=0;
long double total_Syz_B=0;
long double total_Syy_B=0;

int mi, ri, ai, bi;

double d;


int main(int argc, char *argv[]) { /* standard main function */


	int sample_frequency = string_to_int(argv[3]);
	string directory_out=argv[4];

	fileset IF;

	/* Input and output file info
contains file name and pointer
this simplifies file error reporting */
	assembly A; 		/* an empty, unititialized assembly */
//	molindex moli;		/* a molecule index ... m,r,a indices */
	molecule *mtmp;		/* Pointer to a molecule (for convenience) */
	residue *rtmp;		/* Pointer to a residue (for convenience) */
	atom *atmp; 	/* Pointers to atom (for convenience) */
	/* It is always good to check the command line sanity */
	/* glylib contains many convenience functions, for example: */
	//if(argc<2){mywhineusage("Need input prmtop file name.");}


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
	MyAP=(amber_prmtop *)calloc(1,sizeof(amber_prmtop));
	amber_prmtop_init(MyAP); 

	// load info as-is into sections
	read_amber_prmtop_asis(IF,MyAP);

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
				//				atomMass[natom] = (atmp[0]).m;
				atomName[natom] = (atmp[0]).N;
				atomType[natom] = (atmp[0]).T[0];
				//				atomTypeNumber[natom] = (atmp[0]).t;
				//				atomNumber[natom] = (atmp[0]).n;
				//				atomCharge[natom] = (atmp[0]).ch[0];
				atomResidue[natom] = A.m[atmp[0].moli.m][0].r[atmp[0].moli.r].N;
				natom++;

			}
		}
	}


	// calculates # frames in trajectory

	intfile=0;
	ifstream fileofcoordinate(argv[2]);
	getline(fileofcoordinate, text);
	while (getline(fileofcoordinate, text) > 0){
		intfile++;
	}
	intfile = intfile / ceil((double)natom * 3 / 10);;
	cout << intfile << " frames in trajectory " << endl;
	fileofcoordinate.close();

//	cout << intfile/sample_frequency << " used" << endl;

	long double * coordinateAtomx=new long double[intfile];
	long double * coordinateAtomy=new long double[intfile];
	long double * coordinateAtomz=new long double[intfile];

	long double * coordinateAtomx2=new long double[intfile];
	long double * coordinateAtomy2=new long double[intfile];
	long double * coordinateAtomz2=new long double[intfile];
	
	for (int i=0; i<intfile; i++){
		coordinateAtomx[i]=0;
		coordinateAtomy[i]=0;
		coordinateAtomz[i]=0;
		
		coordinateAtomx2[i]=0;
		coordinateAtomy2[i]=0;
		coordinateAtomz2[i]=0;
	}
	
	
	// parameter file
	
	string file_of_order_parameter=argv[4];
	file_of_order_parameter.append("/order_parameter_");
	file_of_order_parameter.append(argv[5]);
	file_of_order_parameter.append("_");
	file_of_order_parameter.append(argv[6]);
	file_of_order_parameter.append("_");
	file_of_order_parameter.append(argv[7]);
	file_of_order_parameter.append("_");
	file_of_order_parameter.append(argv[8]);
	file_of_order_parameter.append("_");
	file_of_order_parameter.append(argv[9]);
	file_of_order_parameter.append("_");
	file_of_order_parameter.append(argv[10]);
	file_of_order_parameter.append(".txt");

	char *file_order_parameter=new char [file_of_order_parameter.length()+1];
	strcpy(file_order_parameter, file_of_order_parameter.c_str());

	ofstream order_parameter(file_order_parameter);


	string file_of_bond_length=argv[4];
	file_of_bond_length.append("/bond_length_");
	file_of_bond_length.append(argv[5]);
	file_of_bond_length.append("_");
	file_of_bond_length.append(argv[6]);
	file_of_bond_length.append("_");
	file_of_bond_length.append(argv[7]);
	file_of_bond_length.append("_");
	file_of_bond_length.append(argv[8]);
	file_of_bond_length.append("_");
	file_of_bond_length.append(argv[9]);
	file_of_bond_length.append("_");
	file_of_bond_length.append(argv[10]);
	file_of_bond_length.append(".txt");

	char *file_bond_length=new char [file_of_bond_length.length()+1];
	strcpy(file_bond_length, file_of_bond_length.c_str());

	ofstream bond_length(file_bond_length);


	order_parameter << "Info for " << argv[5] << ":" << argv[6] << ":" << argv[7] << " " << argv[8] << ":" << argv[9] << ":" << argv[10] << endl << endl;

// find atom numbers

	for (int i = 0; i < natom; i++){
		
		for (int i2 = 0; i2 < natom; i2++){

			if(atomResidue[i].compare(argv[5])==0){
				if(residue[i].compare(argv[6])==0){
					if(atomName[i].compare(argv[7])==0){
						if(atomResidue[i2].compare(argv[8])==0){
							if(residue[i2].compare(argv[9])==0){
								if(atomName[i2].compare(argv[10])==0){
									
									ifirst=i;
									
									isecond=i2;
									
								}								
							}							
						}						
					}					
				}				
			}			
			
		}
		
	}

 cout << ifirst << " " << isecond << endl;

 cout << atomResidue[ifirst] << " " << residue[ifirst] << " " << atomName[ifirst] << endl;

 cout << atomResidue[isecond] << " " << residue[isecond] << " " << atomName[isecond] << endl;

order_parameter << "topology/trajectory " << argv[1] << endl;

order_parameter << argv[2] << endl << endl;
	
order_parameter << "atoms : " << ifirst << " " << isecond << endl << endl;
	
	
// coordinates of atom- coordinateAtomx, coordinateAtomy, coordinateAtomz

				fileofcoordinate.open(argv[2]);
				for (int file_trajectory = 0; file_trajectory < intfile; file_trajectory++){

					fileofcoordinate.seekg(0);
					getline(fileofcoordinate, text);

					if ((3 * natom) % 10 != 0){
						ad = file_trajectory*(81 * ceil((long double)natom * 3 / 10) - (10 - (3 * natom) % 10) * 8);
						fileofcoordinate.seekg(-file_trajectory + ad, ios_base::cur);
					}
					if ((3 * natom) % 10 == 0){
						ad = file_trajectory*(81 * ceil((long double)natom * 3 / 10));
						fileofcoordinate.seekg(-file_trajectory + ad, ios_base::cur);
					}

					fileofcoordinate.seekg(file_trajectory, ios_base::cur);

					for (int unused = 0; unused < ((long double)(3 * ifirst - (3 * ifirst) % 10) / 10); unused++){
						getline(fileofcoordinate, text);
					}
					getline(fileofcoordinate, text);

					coordinate = string_to_double(text.substr((size_t)((3 * ifirst * 8) % 80), 8));
					coordinateAtomx[file_trajectory] = coordinate;

					if (((3 * ifirst * 8) % 80 + 16) != 88){
						coordinate = string_to_double(text.substr((size_t)abs((3 * ifirst * 8) % 80 + 8), 8));
						coordinateAtomy[file_trajectory] = coordinate;

						if (abs((3 * (ifirst)* 8) % 80 + 24) <= 80){
							coordinate = string_to_double(text.substr((size_t)abs((3 * ifirst * 8) % 80 + 16), 8));
							coordinateAtomz[file_trajectory] = coordinate;
						}
						else{
							getline(fileofcoordinate, text);
							coordinate = string_to_double(text.substr((size_t)(0), 8));
							coordinateAtomz[file_trajectory] = coordinate;
						}
					}
					if (abs((3 * ifirst * 8) % 80 + 16) == 88){
						getline(fileofcoordinate, text);
						coordinate = string_to_double(text.substr((size_t)0, 8));
						coordinateAtomy[file_trajectory] = coordinate;
						coordinate = string_to_double(text.substr((size_t)8, 8));
						coordinateAtomz[file_trajectory] = coordinate;
					}
				}
				fileofcoordinate.close();

				// coordinates of atom- coordinateAtomx2, coordinateAtomy2, coordinateAtomz2



				fileofcoordinate.open(argv[2]);
				for (int file_trajectory = 0; file_trajectory < intfile; file_trajectory++){

					fileofcoordinate.seekg(0);
					getline(fileofcoordinate, text);

					if ((3 * natom) % 10 != 0){
						ad = file_trajectory*(81 * ceil((long double)natom * 3 / 10) - (10 - (3 * natom) % 10) * 8);
						fileofcoordinate.seekg(-file_trajectory + ad, ios_base::cur);
					}
					if ((3 * natom) % 10 == 0){
						ad = file_trajectory*(81 * ceil((long double)natom * 3 / 10));
						fileofcoordinate.seekg(-file_trajectory + ad, ios_base::cur);
					}
					fileofcoordinate.seekg(file_trajectory, ios_base::cur);

					for (int unused = 0; unused < ((long double)(3 * isecond - (3 * isecond) % 10) / 10); unused++){
						getline(fileofcoordinate, text);
					}
					getline(fileofcoordinate, text);

					coordinate = string_to_double(text.substr((size_t)((3 * isecond * 8) % 80), 8));
					coordinateAtomx2[file_trajectory] = coordinate;

					if (((3 * isecond * 8) % 80 + 16) != 88){
						coordinate = string_to_double(text.substr((size_t)abs((3 * isecond * 8) % 80 + 8), 8));
						coordinateAtomy2[file_trajectory] = coordinate;

						if (abs((3 * (isecond)* 8) % 80 + 24) <= 80){
							coordinate = string_to_double(text.substr((size_t)abs((3 * isecond * 8) % 80 + 16), 8));
							coordinateAtomz2[file_trajectory] = coordinate;
						}
						else{
							getline(fileofcoordinate, text);
							coordinate = string_to_double(text.substr((size_t)(0), 8));
							coordinateAtomz2[file_trajectory] = coordinate;
						}
					}

					if (abs((3 * isecond * 8) % 80 + 16) == 88){
						getline(fileofcoordinate, text);
						coordinate = string_to_double(text.substr((size_t)0, 8));
						coordinateAtomy2[file_trajectory] = coordinate;
						coordinate = string_to_double(text.substr((size_t)8, 8));
						coordinateAtomz2[file_trajectory] = coordinate;
					}
		}
				fileofcoordinate.close();

	// calculation of average x2-x, ..., spin pair

	coordinateAtom[0]=0;
	coordinateAtom[1]=0;
	coordinateAtom[2]=0;

	coordinate_spin_pair[0]=0;
	coordinate_spin_pair[1]=0;
	coordinate_spin_pair[2]=0;
	
// reusable

	distance_normalized_spin_pair=0;

// difference in coordinate1-coordinate2

        coordinate_diff_x=0;
        coordinate_diff_y=0;
        coordinate_diff_z=0;

        coordinate_spin_1_x=0;
        coordinate_spin_1_y=0;
        coordinate_spin_1_z=0;

        coordinate_spin_2_x=0;
        coordinate_spin_2_y=0;
        coordinate_spin_2_z=0;

        average_distance=0;

// average

	for (int file_trajectory = 0; file_trajectory < intfile; file_trajectory=file_trajectory+sample_frequency){

		coordinate_diff_x = (long double) (coordinateAtomx[file_trajectory]-coordinateAtomx2[file_trajectory]);
		coordinate_diff_y = (long double) (coordinateAtomy[file_trajectory]-coordinateAtomy2[file_trajectory]);
		coordinate_diff_z = (long double) (coordinateAtomz[file_trajectory]-coordinateAtomz2[file_trajectory]);

         distance_normalized_spin_pair=(long double) (pow(coordinate_diff_x*coordinate_diff_x+coordinate_diff_y*coordinate_diff_y+coordinate_diff_z*coordinate_diff_z,(long double) 1/2));

        coordinate_spin_pair[0]=(long double) (coordinate_diff_x+coordinate_spin_pair[0]);
        coordinate_spin_pair[1]=(long double) (coordinate_diff_y+coordinate_spin_pair[1]);
        coordinate_spin_pair[2]=(long double) (coordinate_diff_z+coordinate_spin_pair[2]);

// first atom -  average

        coordinate_spin_1_x=(long double) (coordinateAtomx[file_trajectory]+coordinate_spin_1_x);
        coordinate_spin_1_y=(long double) (coordinateAtomy[file_trajectory]+coordinate_spin_1_y);
        coordinate_spin_1_z=(long double) (coordinateAtomz[file_trajectory]+coordinate_spin_1_z);

// second atom -  average

        coordinate_spin_2_x=(long double) (coordinateAtomx2[file_trajectory]+coordinate_spin_2_x);
        coordinate_spin_2_y=(long double) (coordinateAtomy2[file_trajectory]+coordinate_spin_2_y);
        coordinate_spin_2_z=(long double) (coordinateAtomz2[file_trajectory]+coordinate_spin_2_z);

// average distance 

        bond_distance=(long double) (pow((coordinateAtomx2[file_trajectory]-coordinateAtomx[file_trajectory])*(coordinateAtomx2[file_trajectory]-coordinateAtomx[file_trajectory])+(coordinateAtomy2[file_trajectory]-coordinateAtomy[file_trajectory])*(coordinateAtomy2[file_trajectory]-coordinateAtomy[file_trajectory])+(coordinateAtomz2[file_trajectory]-coordinateAtomz[file_trajectory])*(coordinateAtomz2[file_trajectory]-coordinateAtomz[file_trajectory]),(long double) 1/2));

        bond_length << file_trajectory << " " << bond_distance << endl;

        average_distance=(long double) (pow((coordinateAtomx2[file_trajectory]-coordinateAtomx[file_trajectory])*(coordinateAtomx2[file_trajectory]-coordinateAtomx[file_trajectory])+(coordinateAtomy2[file_trajectory]-coordinateAtomy[file_trajectory])*(coordinateAtomy2[file_trajectory]-coordinateAtomy[file_trajectory])+(coordinateAtomz2[file_trajectory]-coordinateAtomz[file_trajectory])*(coordinateAtomz2[file_trajectory]-coordinateAtomz[file_trajectory]),(long double) 1/2)+average_distance);
	
	}

        bond_length.close();

        average_distance=(long double) average_distance/intfile*sample_frequency;

// average coordinate - not in output - not normalized

        coordinate_spin_1_x=(long double) (coordinate_spin_1_x/intfile*sample_frequency);
        coordinate_spin_1_y=(long double) (coordinate_spin_1_y/intfile*sample_frequency);
        coordinate_spin_1_z=(long double) (coordinate_spin_1_z/intfile*sample_frequency);

        coordinate_spin_2_x=(long double) (coordinate_spin_2_x/intfile*sample_frequency);
        coordinate_spin_2_y=(long double) (coordinate_spin_2_y/intfile*sample_frequency);
        coordinate_spin_2_z=(long double) (coordinate_spin_2_z/intfile*sample_frequency);

// <r_1>-<r_2>

        coordinate_spin_pair[0]=(long double) (coordinate_spin_pair[0]/intfile*sample_frequency);
        coordinate_spin_pair[1]=(long double) (coordinate_spin_pair[1]/intfile*sample_frequency);
        coordinate_spin_pair[2]=(long double) (coordinate_spin_pair[2]/intfile*sample_frequency);

// unit vector - norm of average, average of norm (bond length)

        d=pow(pow(coordinate_spin_pair[0],2)+pow(coordinate_spin_pair[1],2)+pow(coordinate_spin_pair[2],2),.5);

        coordinate_spin_pair[0]=coordinate_spin_pair[0]/d;
		coordinate_spin_pair[1]=coordinate_spin_pair[1]/d;
		coordinate_spin_pair[2]=coordinate_spin_pair[2]/d;
		
		
// average coordinates

order_parameter << "average <x1>,<y1>,<z1> - over frames" << endl;
order_parameter << coordinate_spin_1_x << " " << coordinate_spin_1_y << " " << coordinate_spin_1_z << endl << endl; 
order_parameter << "average <x2>,<y2>,<z2> - over frames" << endl;
order_parameter << coordinate_spin_2_x << " " << coordinate_spin_2_y << " " << coordinate_spin_2_z << endl << endl; 


// 'bond length' normalized coordinates 

//  (\vec x, \vec y)/\sqrt (ave \vec x - \vec y)^2  \avg(distance) 

d = pow( pow((coordinate_spin_1_x-coordinate_spin_2_x),2) + pow((coordinate_spin_1_y-coordinate_spin_2_y),2) + pow((coordinate_spin_1_z-coordinate_spin_2_z),2),.5);

double x1 = (double) coordinate_spin_1_x/d*average_distance;
double y1 = (double) coordinate_spin_1_y/d*average_distance;
double z1 = (double) coordinate_spin_1_z/d*average_distance;

double x2 = (double) coordinate_spin_2_x/d*average_distance;
double y2 = (double) coordinate_spin_2_y/d*average_distance;
double z2 = (double) coordinate_spin_2_z/d*average_distance;

order_parameter << "'bond length' normalized <x1>,<y1>,<z1>" << endl; 
order_parameter << x1 << " " << y1 << " " << z1 << endl << endl;
order_parameter << "'bond length' normalized <x2>,<y2>,<z2>" << endl << endl;
order_parameter << x2 << " " << y2 << " " << z2 << endl << endl;
 

// average distance 

order_parameter << "average distance - average over frames" << endl; 
order_parameter << average_distance << endl << endl;

	
// average of (3cos^2-1)/2 about average x,y,z - over frames

// sampling

	total=0;

        cosine_angle=0;

        cosine_angle_x=0;
        cosine_angle_y=0;
        cosine_angle_z=0;

	for (int file_trajectory = 0; file_trajectory < intfile; file_trajectory = file_trajectory + sample_frequency){
		
// innner product of two unit vectors

        coordinate_diff_x = (long double) (coordinateAtomx[file_trajectory]-coordinateAtomx2[file_trajectory]);
		coordinate_diff_y = (long double) (coordinateAtomy[file_trajectory]-coordinateAtomy2[file_trajectory]);
		coordinate_diff_z = (long double) (coordinateAtomz[file_trajectory]-coordinateAtomz2[file_trajectory]);

         distance_normalized_spin_pair=(long double) (pow(coordinate_diff_x*coordinate_diff_x+coordinate_diff_y*coordinate_diff_y+coordinate_diff_z*coordinate_diff_z,(long double) 1/2));
	
//		cosine_angle_x_B=(long double) (coordinate_diff_x/distance_normalized_spin_pair);
//		cosine_angle_y_B=(long double) (coordinate_diff_y/distance_normalized_spin_pair);
//		cosine_angle_z_B=(long double) (coordinate_diff_z/distance_normalized_spin_pair);

// unit vec

// (r_1-r_2) / |r_1-r_2|

		coordinate_diff_x = (long double) (coordinate_diff_x/distance_normalized_spin_pair);
		coordinate_diff_y = (long double) (coordinate_diff_y/distance_normalized_spin_pair);
		coordinate_diff_z = (long double) (coordinate_diff_z/distance_normalized_spin_pair);

// unit(r_1-r_2) \cdot unit(average (r_1-r_2))

		cosine_angle_x=(long double) (coordinate_diff_x*coordinate_spin_pair[0]);	
		cosine_angle_y=(long double) (coordinate_diff_y*coordinate_spin_pair[1]);
		cosine_angle_z=(long double) (coordinate_diff_z*coordinate_spin_pair[2]);

         cosine_angle=cosine_angle_x+cosine_angle_y+cosine_angle_z;

		total=total+(long double) (3*cosine_angle*cosine_angle-1)/2;	

		total_Sxx=total_Sxx+(long double) (3*cosine_angle_x*cosine_angle_x-1)/2;
		total_Sxy=total_Sxy+(long double) (3*cosine_angle_x*cosine_angle_y)/2;
		total_Sxz=total_Sxz+(long double) (3*cosine_angle_x*cosine_angle_z)/2;
		total_Syz=total_Syz+(long double) (3*cosine_angle_y*cosine_angle_z)/2;
		total_Syy=total_Syy+(long double) (3*cosine_angle_y*cosine_angle_y-1)/2;

		total_Sxx_B=total_Sxx_B+(long double) (3*cosine_angle_x_B*cosine_angle_x_B-1)/2;
		total_Sxy_B=total_Sxy_B+(long double) (3*cosine_angle_x_B*cosine_angle_y_B)/2;
		total_Sxz_B=total_Sxz_B+(long double) (3*cosine_angle_x_B*cosine_angle_z_B)/2;
		total_Syz_B=total_Syz_B+(long double) (3*cosine_angle_y_B*cosine_angle_z_B)/2;
		total_Syy_B=total_Syy_B+(long double) (3*cosine_angle_y_B*cosine_angle_y_B-1)/2;

	}
	
// average

	total=(long double) total/intfile*sample_frequency;

	total_Sxy=(long double) total_Sxy/intfile*sample_frequency;
	total_Sxz=(long double) total_Sxz/intfile*sample_frequency;
	total_Syz=(long double) total_Syz/intfile*sample_frequency;
	total_Sxx=(long double) total_Sxx/intfile*sample_frequency;
	total_Syy=(long double) total_Syy/intfile*sample_frequency;

	total_Sxy_B=(long double) total_Sxy_B/intfile*sample_frequency;
	total_Sxz_B=(long double) total_Sxz_B/intfile*sample_frequency;
	total_Syz_B=(long double) total_Syz_B/intfile*sample_frequency;
	total_Sxx_B=(long double) total_Sxx_B/intfile*sample_frequency;
	total_Syy_B=(long double) total_Syy_B/intfile*sample_frequency;

	
    cout << atomResidue[ifirst] << " " << residue[ifirst] << " " << atomName[ifirst] << "  " << atomResidue[isecond] << " " << residue[isecond] << " " << atomName[isecond] << endl << endl;
	
    order_parameter << "order parameter <(3cos^2-1)/2>" << endl;
	order_parameter << total << endl << endl;

    order_parameter << "average fluctuation: Sxx,Syy,Sxy,Sxz,Syz" << endl;

    order_parameter << total_Sxx << endl;
    order_parameter << total_Syy << endl;
    order_parameter << total_Sxy << endl;
    order_parameter << total_Sxz << endl;
    order_parameter << total_Syz << endl;

    order_parameter << (double) 2/3*pow(total_Sxx*total_Sxx+total_Syy*total_Syy+(total_Sxx+total_Syy)*(total_Sxx+total_Syy)+2*total_Sxy*total_Sxy+2*total_Sxz*total_Sxz+2*total_Syz*total_Syz,(double) 1/2);
    order_parameter << endl << endl;

    order_parameter << "about B-direction: Sxx,Syy,Sxy,Sxz,Syz" << endl;

    order_parameter << total_Sxx_B << endl;
    order_parameter << total_Syy_B << endl;
    order_parameter << total_Sxy_B << endl;
    order_parameter << total_Sxz_B << endl;
    order_parameter << total_Syz_B << endl;

    order_parameter << (double) 2/3*pow(total_Sxx_B*total_Sxx_B+total_Syy_B*total_Syy_B+(total_Sxx_B+total_Syy_B)*(total_Sxx_B+total_Syy_B)+2*total_Sxy_B*total_Sxy_B+2*total_Sxz_B*total_Sxz_B+2*total_Syz_B*total_Syz_B,(double) 1/2);

	order_parameter.close();


// bond figure

	char *atomName1 = new char[atomName[ifirst].length() + 1];
    char *atomName2 = new char[atomName[isecond].length() + 1];
	char *atomResidue1 = new char[atomResidue[ifirst].length() + 1];
	char *atomResidue2 = new char[atomResidue[isecond].length() + 1];
	char *residue1 = new char[residue[ifirst].length()+1];
	char *residue2 = new char[residue[isecond].length()+1];
				
	strcpy(atomName1, atomName[ifirst].c_str());
	strcpy(atomName2, atomName[isecond].c_str());
	strcpy(atomResidue1, atomResidue[ifirst].c_str());
	strcpy(atomResidue2, atomResidue[isecond].c_str());
	strcpy(residue1, residue[ifirst].c_str());
	strcpy(residue2, residue[isecond].c_str());

    FILE *gp;
	gp = popen(GNUPLOT, "w");
	if (gp == NULL){
	   return 0;
        }

	fprintf(gp, " set terminal jpeg \n unset key \n set output '%s/%s_%s_%s-%s_%s_%s_bond_length.jpeg' \n set xlabel 'frame number' \n set title 'bond length in Angstroms - atom %s:%s:%s  atom %s:%s:%s \n plot '%s' \n ", argv[4], atomResidue1, residue1, atomName1, atomResidue2, residue2, atomName2, atomResidue1, residue1, atomName1, atomResidue2, residue2, atomName2,file_bond_length);

	return 0; /* indicate successful completion */

}
