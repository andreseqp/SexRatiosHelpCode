// Gil's ovicide model

//#include "stdafx.h" // delete in g++ version

#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <algorithm>
#include <array>
#include <iostream>
#include <fstream>
#include <functional>
#include "random.h"


double phi;  // efficiency recycling males
double b;  // extra offspring (per queen offspring) added by a helper
const double Sf         = 0.8;  // survival spring females to breed again in summer
const double Sm         = 0.8;  // survival summer males (to mate again in autumn)
const int F             = 5;   // # offspring without help

const int N				= 2500;		// # nests in spring
const double mu			= 0.001;    // mutation rate per allele per generation
const double sdmu		= 0.2;		// standard deviation mutation size

const int NumGen		= 25000;	// # generations
const int Skip			= 10;		// interval for writing output

const int IntroHelp     = 1000;       // generation after which mutation on help locus (and o and z3) is allowed


struct Female // diploid genotype
{
	double z1[2]; // spring sex ratio (proportion sons)
	double z2[2]; // summer sex ratio without helpers
	double z3[2]; // summer sex ratio with helpers
	double h[2];  // helping tendency of spring daughters
	double o[2];  // probability to kill male egg
	// NB: regardless of # helpers, all male eggs subjected to this ovicide probability
	// This may be a bit unrealistic
};

struct Male // haploid genotype
{
	double z1;
	double z2; 
	double z3; //
	double h;
	double o;
};

struct Nest // or just a mated female
{
	Female Q;
	Male K;  // king
	int H; 
};


Nest F1[N]; // spring nests (mated females emerging from hibernation)  
Nest F1_0[100*N]; // mated females that will form next generation's spring nests 
Nest F2[N+N*F];  // summer nests (mated females + count of # helpers
                 // note that helper genotypes are not stored
Male M1[100*N]; // autumn males (summer sons + surviving spring sons)
Male M2[N*F];  // summer sons

int NF1,NF2,NM1,NM2; // counters for arrays

int Generation;          // generation counter
int Seed; // RNG seed

double means[5]; // z1 z2 z3 h o
double sds[5]; // corresponding standard deviations
double Q1[5];
double Q3[5];

FILE *EvolFile;  // output file evolutionary dynamics
FILE *DistFile;  // output file phenotype distributions

// rename random number generators of random.h
int rn(int n){ return RandomNumber(n);}  // random integer
double ru(){ return Uniform();}          // random double between 0 and 1

// initialize stuff
void Init()
{
	int i;
	//Seed=(unsigned)time(NULL);  // "seed" random number generator
	SetSeed(Seed);

	// create first spring nests of mated females without helpers
	for (i=0;i<N;i++)
	{
		F1[i].Q.z1[0]=0.5;
		F1[i].Q.z1[1]=0.5;
		F1[i].Q.z2[0]=0.5;
		F1[i].Q.z2[1]=0.5;
		F1[i].Q.z3[0]=0.5;
		F1[i].Q.z3[1]=0.5;
		F1[i].Q.h[0]=0;
		F1[i].Q.h[1]=0;
		F1[i].Q.o[0]=0;
		F1[i].Q.o[1]=0;

		F1[i].K.z1=0.5;
		F1[i].K.z2=0.5;
		F1[i].K.z3=0.5;
		F1[i].K.h=0;
		F1[i].K.o=0;

		F1[i].H=0;
	}

}



// write data to screen and file
void WriteData()
{
	printf("%8i%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
	Generation,means[0],means[1],means[2],means[3],means[4],sds[0],sds[1],sds[2],sds[3],sds[4]);

	fprintf(EvolFile,"%8i%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
	Generation,means[0],means[1],means[2],means[3],means[4],sds[0],sds[1],sds[2],sds[3],sds[4], Q1[0], Q1[1], Q1[2], Q1[3], Q1[4], Q3[0], Q3[1], Q3[2], Q3[3], Q3[4]);
}


void WriteDist()  // phenotype distribution among spring foundresses
{
	for (int i=0;i<N;i++)
	{
		fprintf(DistFile,"%10.3f%10.3f%10.3f%10.3f%10.3f\n",
		0.5*(F1[i].Q.z1[0]+F1[i].Q.z1[1]),
		0.5*(F1[i].Q.z2[0]+F1[i].Q.z2[1]),
		0.5*(F1[i].Q.z3[0]+F1[i].Q.z3[1]),
		0.5*(F1[i].Q.h[0]+F1[i].Q.h[1]),
		0.5*(F1[i].Q.o[0]+F1[i].Q.o[1]));
	}
}


// mutate an allele
void Mutate(double &g)
{
	if (ru()<mu) g+=Normal(0,sdmu);
	if (g<0) g=0; else if (g>1) g=1;
}

void Statistics()
{
	int i;
	double z1,z2,z3,h,o;
	double ss[5];
	std::array<std::array<double, N>, 5> dist;

	for (i=0;i<5;i++) { means[i]=0; ss[i]=0;}

	for (i=0;i<N;i++)
	{
		z1=0.5*(F1[i].Q.z1[0]+F1[i].Q.z1[1]);
		z2=0.5*(F1[i].Q.z2[0]+F1[i].Q.z2[1]);
		z3=0.5*(F1[i].Q.z3[0]+F1[i].Q.z3[1]);
		h=0.5*(F1[i].Q.h[0]+F1[i].Q.h[1]);
		o=0.5*(F1[i].Q.o[0]+F1[i].Q.o[1]);
		means[0]+=z1;
		dist[0][i] = z1;
		means[1]+=z2;
		dist[1][i] = z2;
		means[2]+=z3;
		dist[2][i] = z3;
		means[3]+=h;
		dist[3][i] = h;
		means[4]+=o;
		dist[4][i] = o;
		ss[0]+=z1*z1;
		ss[1]+=z2*z2;
		ss[2]+=z3*z3;
		ss[3]+=h*h;
		ss[4]+=o*o;
	}

	for (i = 0; i < 5; i++)
	{
		means[i] /= N;
		std::sort(dist[i].begin(), dist[i].end());
	}
	for (i = 0; i < 5; i++)
	{
		sds[i] = sqrt(ss[i] / N - means[i] * means[i]);
		Q1[i] = 0.5*(dist[i][N*0.25]+ dist[i][N*0.25-1]);
		Q3[i] = 0.5*(dist[i][N*0.75] + dist[i][N*0.75 - 1]);
	}
}


void NextGeneration()
{
	int i,j,k;
	double z1, z2, z3, h;
	Male Son;
	Female Daughter;
	int RF;
	double EF;
	int FH;

	NF2=0; // set # spring daughters to zero
	NM2=0; // ditto spring sons

	// reproduction spring
	for (i=0;i<N;i++)
	{
		for (j=0;j<F;j++) // F is the clutch size
		{
			z1=0.5*(F1[i].Q.z1[0]+F1[i].Q.z1[1]);
			if (ru()<z1) // it's a son
			{
				Son.z1=F1[i].Q.z1[rn(2)];
				Son.z2=F1[i].Q.z2[rn(2)];
				Son.z3=F1[i].Q.z3[rn(2)];
				Son.h=F1[i].Q.h[rn(2)];
				Son.o=F1[i].Q.o[rn(2)];

				Mutate(Son.z1);
				Mutate(Son.z2);
				if (Generation>IntroHelp) 
				{
					Mutate(Son.z3);
					Mutate(Son.h);
					Mutate(Son.o);
				}
					else
				{
					Son.z3 = Son.z2;
				}

				
				M2[NM2]=Son; // add son to array
				NM2++;
			}
			else // it's a daughter
			{
				Daughter.z1[0]=F1[i].Q.z1[rn(2)];
				Daughter.z1[1]=F1[i].K.z1;
				Daughter.z2[0]=F1[i].Q.z2[rn(2)];
				Daughter.z2[1]=F1[i].K.z2;
				Daughter.z3[0]=F1[i].Q.z3[rn(2)];
				Daughter.z3[1]=F1[i].K.z3;
				Daughter.h[0]=F1[i].Q.h[rn(2)];
				Daughter.h[1]=F1[i].K.h;
				Daughter.o[0]=F1[i].Q.o[rn(2)];
				Daughter.o[1]=F1[i].K.o;

				for (k=0;k<2;k++)
				{
					Mutate(Daughter.z1[k]);
					Mutate(Daughter.z2[k]);
					if (Generation>IntroHelp) 
					{
						Mutate(Daughter.z3[k]);
						Mutate(Daughter.h[k]);
						Mutate(Daughter.o[k]);
					}
						else
					{
						Daughter.z3[k] = Daughter.z2[k];
					}

				}
				
				// will the daughter become a helper at her mother's nest?
				h=0.5*(Daughter.h[0]+Daughter.h[1]);
				if (ru()<h) F1[i].H++; // yes, she will. Increase # helpers by 1
				else // no, she doesn't. Add her to the non-helping spring daughters
				{
					F2[NF2].Q=Daughter;
					F2[NF2].H=0; // indepdently breeding daughters have no helpers
				    NF2++;
				}
			}
		}
	}

	// mating of daughters with random spring son
	for (i=0;i<NF2;i++) F2[i].K=M2[rn(NM2)];

	// mortality spring foundresses
	// survivors will breed again and may have helpers
	// they are added to the array of non-helping spring daughters
	for (i=0;i<N;i++)
	{
		if (ru()<Sf)
		{
			F2[NF2]=F1[i];
			NF2++;
		}
	}

	// mortality spring sons
	NM1=0;  // set autumn males to zero. 
	// these males will include both summer sons and spring sons that survive to mate again
	for (i=0;i<NM2;i++)
	{
		if (ru()<Sm)
		{
			M1[NM1]=M2[i];
			NM1++;
		}
	}

	// summer reproduction
	NF1=0;
	double EO; // expected ovicide
	double U,Y,Z;

	for (i=0;i<NF2;i++)
	{
		// if # helpers>0 allow ovicide (avg of queen's and king's genotype)
		if (F2[i].H>0) EO=0.25*(F2[i].Q.o[0]+F2[i].Q.o[1])+0.5*F2[i].K.o; else EO=0;
		EF=F*(1+b*F2[i].H); // mean clutch size, dependent on # helpers
		if (ru()<EF-int(EF)) FH=int(EF+1); else FH=int(EF); // round clutch size

		for (j=0;j<FH;j++)
		{
	
			if (F2[i].H==0) // no helpers. Use z2 as sex ratio
			{
				Y=0.5*(F2[i].Q.z2[0]+F2[i].Q.z2[1]);
				Z=1;
			}
			else // helpers; use z3
			{
				z3=0.5*(F2[i].Q.z3[0]+F2[i].Q.z3[1]);
				Y=z3*(1-EO+z3*phi*EO);
				Z=1-z3*EO*(1-phi);
				// note that Y+Z<1 if EO>0
				// that means that the # offspring is smaller when EO>0
			}
			
			U=ru();
			if (U<Y) // a son
			{
				Son.z1=F2[i].Q.z1[rn(2)];
				Son.z2=F2[i].Q.z2[rn(2)];
				Son.z3=F2[i].Q.z3[rn(2)];
				Son.h=F2[i].Q.h[rn(2)];
				Son.o=F2[i].Q.o[rn(2)];

				Mutate(Son.z1);
				Mutate(Son.z2);
				if (Generation>IntroHelp) 
				{
					Mutate(Son.z3);
					Mutate(Son.h);
					Mutate(Son.o);
				}
				else
				{
					Son.z3 = Son.z2;
				}
				
				M1[NM1]=Son;
				NM1++;
			}
			else if (U<Z) // a daughter, or else no offspring due to ovicide
			{
				Daughter.z1[0]=F2[i].Q.z1[rn(2)];
				Daughter.z1[1]=F2[i].K.z1;
				Daughter.z2[0]=F2[i].Q.z2[rn(2)];
				Daughter.z2[1]=F2[i].K.z2;
				Daughter.z3[0]=F2[i].Q.z3[rn(2)];
				Daughter.z3[1]=F2[i].K.z3;
				Daughter.h[0]=F2[i].Q.h[rn(2)];
				Daughter.h[1]=F2[i].K.h;
				Daughter.o[0]=F2[i].Q.o[rn(2)];
				Daughter.o[1]=F2[i].K.o;

				for (k=0;k<2;k++)
				{
					Mutate(Daughter.z1[k]);
					Mutate(Daughter.z2[k]);
					if (Generation>IntroHelp)
					{
						Mutate(Daughter.z3[k]);
						Mutate(Daughter.h[k]);
						Mutate(Daughter.o[k]);
					}
					else
					{
						Daughter.z3[k] = Daughter.z2[k];
					}
				}
				
				F1_0[NF1].Q=Daughter;
			    NF1++;
			}
		}
	} // summer reprod

	// mating summer daughters with summer sons and surviving spring sons
	for (i=0;i<NF1;i++) F1_0[i].K=M1[rn(NM1)];

	// cull excess daughters back to size N
	for (i=0;i<N;i++)
	{
		RF=rn(NF1);
		F1[i]=F1_0[RF];
		F1[i].H=0;
		NF1--;
		F1_0[RF]=F1_0[NF1];
	}

}


int main(int argc, char* argv[])
{
    //char EvolFileName[40] = "";//"/data/p264680/Ovicide/paper/evol_";
	//char DistFileName[40] = "";//"/data/p264680/Ovicide/paper/dist_";
	phi = 0.7;//atof(argv[1]) / 10;
	b   = 0.9;//atof(argv[2]) / 10;

	for (size_t i = 0; i < 10; i++)
	{

		Seed = i;//atoi(argv[3]);
		char EvolFileName[100] = "D:\\quinonesa\\Dropbox\\Haplodiploidy\\Gil\\Andres\\IBD\\evol";
		char DistFileName[100] = "D:\\quinonesa\\Dropbox\\Haplodiploidy\\Gil\\Andres\\IBD\\dist";
		char buffer[40];
		sprintf(buffer, "%1i_seed_%02.0f_phi_%02.0f_b%", Seed, phi*10, b*10);
		strcat(EvolFileName, buffer);
		strcat(EvolFileName, ".txt");
		EvolFile = fopen(EvolFileName, "w");
		//sprintf(Buffer, "%1i_seed_%02.0f_phi_%02.0f_b%", Seed, phi * 10, b * 10);
		strcat(DistFileName, buffer);
		strcat(DistFileName, ".txt");
		DistFile = fopen(DistFileName, "w");
		/*strcat(EvolFileName,".txt");
		strcat(DistFileName,".txt");*/
		
		

		Init();
		Statistics();

		Generation = 0;
		WriteData();


		for (Generation = 1; Generation <= NumGen; Generation++)
		{
			NextGeneration();
			if (Generation%Skip == 0)
			{
				Statistics();
				WriteData();
			}
		}

		WriteDist();

		fclose(EvolFile);
		fclose(DistFile);
	}
	
	;
	return 0;
}
