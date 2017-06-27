                         /*  ecovolve.c  */

/* THIS SOFTWARE IS COVERED BY A GNU GPL <www.gnu.org> - it must not find */
/* its way into commercial software, and its offspring must contain       */
/* acknowledgment of the original author.                                 */

/* Cam Webb <cam_webb@yahoo.com> */

/* $Header: /cvsroot/phylocom/main/ecovolve.h,v 3.38 2005/09/10 00:53:18 ctenolophon Exp $ */

// Routine:
// new time
//   {
//      determine competetive environ
//      (begin real time)
//      evolve chars on extant lineages
//      (end real time)
//      speciate; new spp get old state; make mid-time marker in t+1 on.
//      extinct; active off (newly speciated clades cannot go extinct)
//        else: pass mid-time lineage markers to future 
//   }

// I.e. for time x lineage is before node




/* INCLUDE HEADERS ------------------------------------------------------- */

#include <stdio.h>
#include <math.h>  // for sqrt()
#include <stdlib.h>
#include "phylocom.h"

/* DEFINITIONS FOR MAIN PROGRAM ------------------------------------------ */

#define OUTMODE 3  // 0 = debug, 1 = .phy out, 2 = lineage-thr-time
#define P_SPECIATE 0.05
#define P_EXTINCT 0.01
//#define COMPETE 0
//#define MAXNODES 50000 // now defined in phylocom
#define MAXTIME 100


// for fy2new
//#define MAXTAXONLENGTH 100
//#define MAXBLLENGTH 15
//#define MAXNOTELENGTH 100

/* FUNCTION DECLARATION -------------------------------------------------- */

void Speciate();
void Extinct();
void CharChange();
void Output();
void Compete();
void MakePhylo();
void WriteTraits();
void DummySample();
struct shift MakeChange();
//struct phylo Prune();  // maybe move to phylocom.h?
float Balance();

/* GLOBAL VARIABLES ------------------------------------------------------ */

int OUT_MODE;
int MAX_TIME;
int PRUNED;
int NCHAR;
int TAPER = 0;
int TAPERFACT = 1;
int COMPETE = 0;

int Lineage; // incremental counter of higher Lineage no
int Node;    // ditto for node
int Name;    // ditto for name

// tree tracing:
int Active_l[MAXNODES]; // active marker during tree gen
int Living_lt[MAXNODES][MAXTIME+1]; // tracer for tree history
int Censor_lt[MAXNODES][MAXTIME+1];     // tracer for cenored tree
int EcoDist_l[MAXNODES];

// lineage states
int LineageUp_l[MAXNODES];          
int NodeUp_l[MAXNODES];
float ***Char_ltc;
int Extant_l[MAXNODES];

// node states:
int Time_n[MAXNODES];
int NodeUp_n[MAXNODES];
int BlUp_n[MAXNODES];
int Name_n[MAXNODES];
float **Char_nc;
int Terminal_n[MAXNODES];

float *CharMin, *CharMax;
float *CharMean, *CharSumXSqr;
float *CharCMin, *CharCMax;
float *CharCMean, *CharCSumXSqr;
int LttZ;
int LttC;

phylo OutTree;
// phylo PrunedTree;

typedef struct shift {
  int n;
  int *x;
} shift;

shift charch;

/* PROGRAM END ---------------------------------------------------------- */

// old log messages:
/*
 * Revision 1.4  2005/08/18 03:37:55  cam
 * NAF added
 *
 * Revision 1.3  2005/05/02 13:06:11  cam
 * *** empty log message ***
 *
 * Revision 1.2  2005/05/01 17:17:43  cam
 * *** empty log message ***
 *
 * Revision 1.1  2005/05/01 16:10:22  cam
 * Initial revision
 *
 * Revision 1.3  2003/05/09 06:03:17  cam
 * end of first version.
 *
 * Revision 1.2  2003/04/23 17:53:13  cam
 * Good working version, including speciate, extinct, evolve
 *
 * Revision 1.1  2003/04/22 15:25:47  cam
 * Initial revision
 *
 * */



