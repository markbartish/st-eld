#ifndef __PHYLOCOM_H
#define __PHTLOCOM_H



                  /*  **************************  */
                  /*  **************************  */
                  /*  *******            *******  */
                  /*  ******   phylocom   ******  */
                  /*  *******            *******  */
                  /*  **************************  */
                  /*  **************************  */

// Analysis of community phylogenetic structure and associated traits

// Cam Webb, Arnold Arboretum, Harvard University, CTFS
//           cwebb@oeb.harvard.edu
//
// David Ackerly, Dept. of Integrative Biology, UC Berkeley
//           dackerly@berkeley.edu
//
// Steven Kembel, Dept. of Integrative Biology, UC Berkeley
//			 skembel@berkeley.edu

//COPYRIGHT
//Phylocom Version 4.0 software
//Copyright (C) 2008 Cam Webb, David Ackerly, Steven Kembel
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation.
//
//http://www.gnu.org/copyleft/gpl.html
//
//Please acknowledge us if you use code from this software package.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

// Programming notes:                                                      
//                                             
//   <0> indicates counter starts at 0; <1> ... at 1                   
//   Randomization routines differ between flavors of C and UNIX         
//	((int) (((float) X * random()) / (RAND_MAX+1.0))) + 1            
//      produces a random integer between 1 and X (inclusive).         
//      Substitute your routine that does the same thing               

// INCLUDE HEADERS -------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// DEFINITIONS FOR MAIN PROGRAM ------------------------------------------

//#define RAND_MAX ( 2147483647 )
//#define RANDOM_MAX ( 2147483647 )
#define INFILEP "phylo"
#define INFILEN "phylo.new"
#define INFILET "sample"
#define INFILEM "means"
#define INFILEC "traits"
#define INFILEA "ages"
#define INFILEO "lttrsplits"
#define VERSION "4.0 (SVN $Revision: 138 $)\n"

#define MAXNODES 15000 // Higher than the largest no of nodes in phylo
#define MAXTAXA  13000 // Higher than the highest expected code for taxon
#define MAXRUNS 999   // DDA comment: Changed to 999 and significance counters 
                      // initialized at one
#define MAXLEVEL 1.0
#define MAXTAXONLENGTH 100 // Number of chars for taxon name
#define MAXBLLENGTH 15 // Number of chars for BL
#define MAXPLOTLENGTH 100 // Number of chars for plot name
#define MAXNOTELENGTH 100
#define MAXTRAITLINE 1000

// DDA:
#define MAXTRAITS 4   // Maximum number of traits in trait files
// #define MAXPOLYTAXA 100 // maximum polytomy size for independent contrasts

// SWK:
#define MAXSWAPS 1000 // Default number of swaps for each Independent Swap
#define TRUE 1
#define FALSE 0

// ENUMS:
enum OUTPUT_MODE { OUT_ANR_E = 1, OUT_LOGNOBS_E, OUT_LOGNRND_E, OUT_STDEVLOGNRND_E, OUT_DIFFCOEFFS_E, OUT_LASTMODE };

/* FUNCTION DECLARATION -------------------------------------------------- */

void NodeSig();
void Means();
void VMeans();
void Clust();
void ClustInt();  // internal calculation of means
void ReadData();
void SortDistrib();
int  Rel(int A, int B);
void NRI();
float NR();
void NTI();
float NT();
void Slide();
float SlidingN();
void AppendNote();
// void PrintHeader();
void Showlevels();
void Reshuffle();
void Sort();
void Randomize();
void RandomizeB();
void PrintWelcome();
void FormatHelp();
float Relatedness();
void SimpleDist();
void PhyloVarCovar();
void ComDist();
void ComDistNN();
void Randomspp();
void Ltt();
void LttR();
void LttRv2();
struct phylo New2fy();
void Fy2new();
struct sample ReadSample();
struct phylo ReadPhylogeny();
struct means ReadMeans();
struct lttr_intvs ReadLttrIntv();
void AttachSampleToPhylo();
void AttachSampleToTraits();
void AttachTraitsToPhylo();
void AttachPhyloToTraits();
void DistMatrix();
void DistMatrixNN();
float DistToRootNode();
int FindMRCA();
void NewickToNexus();
void WriteNexus();
void NAF();
void AgeNodes();
struct traits ReadTraits();
void PD();
void License();
void Bladj();
int CleanPhy();
int LineOfSight();
void SortAction();
void Adjust();
//void Polytom();
//void ReadDataBladj();
void ComTraitMetric();
// For reading line endings:
char *myfgets();
int whatnewline();
void IComDist();
void IComDistNN();
void VComDist();
void VComDistNN();
void SampleDeepCopy();

//Comnode
void Comnode();

//Ecovolve
struct phylo Prune();
void RandPrune();
void SamplePrune();

// New recursive Newick-writing functions
void Fy2newRec();
struct phylo SetNodePointers();
char *downPar();

// DDA:
void AOT();
void NodeCharF();
void TipStats();
float *TraitsAtNode();
void SigCount();
void PIC();
void binPIC();
void aot_outfile();
void aot_outscreen();
void RandArray();

// to be deleted from aot
float *summaryStats();
float correlation();

// in traits.c - to be deleted from traits.c when aot finished
void PSig();
void PSigRun();
void RandArrayT();

// DA additions to io.c
void MakeUpPassOrder();
void AssignNodeLists();

// SWK:
void ComStruct(); //SWK
void IndependentSwap(); //SWK
void OutputSwappedMatrix(); //SWK
void PhylogenySampleTaxaShuffle(); //SWK
void PhylogenyAttachShuffle(); //SWK
void TraitsAttachShuffle(); //SWK
void RandomizeSampleTaxaShuffle(); //SWK
double MeanDistance(); //SWK
double MeanMinimumDistance(); //SWK
void traitMetric(); //SWK
void CommunityDistance();
void PhyloDiversity();
void RaoDiversity();


/* GLOBAL VARIABLES ------------------------------------------------------ */
// Capital first letter for global variables (generally)
// small first letters for internal vars

FILE *Fp;   // pointer to phylo
FILE *Fn;   // pointer to phylo.new
FILE *Ft;   // pointer to sample
FILE *Fm;   // pointer to means
FILE *Fc;   // pointer to traits
FILE *Fa;   // pointer to age file
FILE *Fo;   // pointer to nodeage file

char PhyloFile[50]; // default name
char SampleFile[50]; // default name
char TraitFile[50]; // default name
char LttrIntvsFile[50]; //default name
//int UseFy; // switch for using .fy format input
int NoBL; // switch for ignoring branch lenghts

int RUNS, TRAITS, SWAPS, XVAR, AOTOUT, SWAPMETHOD, OUTMODE, RNDPRUNEN, RNDPRUNET, MAKENODENAMES;
int Debug;
int Verbose;
float HILEVEL;
int LowSig; // global switch for low vs. high one-tailed sig testing
int UseAbund; //use abundance data when available?
int OutputLttrRawFlag;
int OutputVariancesFlag;

typedef struct nodes {
	float ***tCh; // tip based character stats
	float ***tChLSig; // one tailed low p vals for tip stats
	float ***tChHSig; // one tailed high p vals for tips stats
	float ***nCh; // node based stats at each node
	float ***nChLSig; // one tailed low p vals for node stats
	float ***nChHSig; // one tailed high p vals for node stats
	float **iCon; // independent contrasts
	float *cSt; // st dev of i.contrast
	int *ordTrt;
	int **rndArr; // vector for random tip sorts
} nodes;

typedef struct phylo {
  char phyname[MAXPLOTLENGTH];
  int nnodes; // equals highest node number plus one, because of 0 node)
  int *up; //up[node]
  int **down;
  int *ldown;
  int *rsister;
  int *noat;
  int *depth; //depth[node]
  float *bl; //bl[node]
  float *tbl;
  float *age; //age[node]
  char **taxon; // name of named node - taxon[node][]
  int ntaxa;  // number of named nodes = total number of names
  int termtaxa; // number of terminal taxa
  char **taxalist; //names of terminal taxa - taxalist[0 to termtaxa-1][]
  int *t2n;  //vector of node #s indexed by 0 to termtaxa-1, as taxalist
  float **dist; // matrix of all node-to-node distances dist[node1][node2]
  int arenotes; // 0 | 1
  char **notes;
  // DA additions
  int maxDepth;
  int *ntip;
  int *nint;
  int **tiplist;
  int **intlist;
  int *upo; // up pass order
} phylo;

//TODO change all algs to just use taxon, checking for terminal status (?)

typedef struct sample { // really need to switch this into a S x T matrix! 
  int     nsamples; // # of samples
  int     nrec;     // total # recs in file
  int     maxrec;   // max # of taxa in a sample - use to dim id, abund;
  char  **pname;   // sample names
  int    *srec;     // species per sample srec[sample]
  int	*irec;		// individuals per sample irec[sample]
  int   **id;       // taxon codes (as in taxa) id[sample][rec]
  int   **abund;    // species abundance  abund[sample][rec]
  float **pabund;	// proportional species abundance[sample][rec]
  unsigned long *sppabund;	//species total abundance[id]
  float *psppabund;	//proportional species total abundance[id]
  unsigned long *sppfreq;	//species occurrence frequency[id]
  float *psppfreq;	//proportional species occurrence frequency[id]
  long totabund;	//total abundance across all species
  int     ntaxa;    // number of unique taxa
  char  **taxa;     // vactor of char strings for taxa names
} sample;

typedef struct traits {
  int     ntraits; // # of traits
  char  **trname;
  int     ntaxa;   // # n taxa for which there are traits
  char  **taxon;   // taxon names
  int    *type;    // trait type 0 binary, 3 cont
  float **tr;      // data[taxon#<0>][trait#<0>]
} traits;

typedef struct means {
  float  *mpd;
  float  *sdpd;
  float  *mnnd;
  float  *sdnnd;
} means;

typedef struct lttr_intvs {
  float   *age_bounds;
  int     nintvs;
} lttr_intvs;

// Sorting params, passed from RemoveDups to Shuffle
int TaxaOrder[MAXTAXA+1]; // sp code<1> of the nth ordered taxon<1>
int FirstSingleton, LastSingleton; // the rank<1> of the first/last single

// counters used in main, passed to functions
int Select;       // subset to use
int Sample;       //  The plot (or Sample) counter
int PlotsUsed;    //  The number of plots used, that exceed MINIMUM
char Method[20]; // type of analysis
int TreeView; //TODO fix this and related NodeSig function

// DDA:
// trait Values
float **Char;
int *CharType; //0 = binary; 1 = multistate; 2 = ordered multistate; 3 = continuous // CW changed 9apr04
// trait conservatism
//int **RndArr;
#endif //__PHYLOCOM_H
