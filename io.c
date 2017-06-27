// io.c Reading data, helpfiles
// Owner: Cam, Steve, David

// includes: ReadData, PrintWelcome, FormatHelp, AppendNote

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "phylocom.h"
#include "nrutil.h"



// ----- Node-as-factor--------------------------------------

void NAF(phylo P[], sample S, traits T)
{

  int i, j, k, x, y, z, found;
  float data;

  printf("Tx\\Dx");
  z = 0;
  for (i = 0; i < T.ntraits; i++)
    {
      printf("\tT:%d:%s", T.type[i], T.trname[i]);
    }
  for (i = 0; i < S.nsamples; i++)
    {
      printf("\tS:%s",S.pname[i]);
    }
  for (i = 0; i < P[0].nnodes; i++)
    {
      if (P[0].noat[i] != 0)
      	{
	  printf("\tN:%s",P[0].taxon[i]);
	}
    }

  // Two header:
/*   printf("\nname:"); */
/*   for (i = 0; i < T.ntraits; i++) */
/*     { */
/*       printf("\t%s", T.trname[i]); */
/*     } */
/*   for (i = 0; i < S.nsamples; i++) */
/*     { */
/*       printf("\t%s", S.pname[i]); */
/*     } */
/*   for (i = 0; (i < P[0].nnodes); i++) */
/*     { */
/*       if (P[0].noat[i] != 0) */
/*       	{ */
/* 	  printf("\t%s", P[0].taxon[i]); */
/* 	} */
/*     } */

   printf("\n");

  for (i = 0; i < P[0].nnodes; i++)
    {
      if (P[0].noat[i] == 0)
	{
	  printf("%s" , P[0].taxon[i]);
	  
	  for (j = 0; j < T.ntraits; j++)
	    {
	      data = 0.0;
	      for (k = 0; k < T.ntaxa; k++)
		{
		  if (strcmp(T.taxon[k], P[0].taxon[i]) == 0) data = T.tr[k][j];
		}
	      printf("\t%f", data);
	    }

	  for (j = 0; j < S.nsamples; j++)
	    {
	      data = 0.0;
	      for (k = 0; k < S.srec[j]; k++)
		{
		  if (strcmp(S.taxa[S.id[j][k]], P[0].taxon[i]) == 0) data = S.abund[j][k];
		}
	      printf("\t%f", data);
	    }

	  for (j = 0; (j < P[0].nnodes); j++)
	    {
	      if (P[0].noat[j] != 0)
		{
		  x = i;
		  found = 0;
		  while (x != -1)
		    {
		      y = x;
		      x = P[0].up[x];
		      if (x == j) 
			{
			  found = 1; 
			  z = y; 
			  break;
			}
		    }
		  if (found == 0) printf("\t.");
		  else printf("\t%s", P[0].taxon[z]);
		}
	    }
	  
	  printf("\n");

	}
    }

}

// ---------- READ MEANS -----------------------------------------------------

means ReadMeans(phylo Intree, char meansfile[50])
{

  char line[500];
  int linesread = 0;
  float meanI, sdI, nnmeanI, nnsdI; 
  int samplesizeI;

  means InM;
  InM.mpd = vector(0, Intree.termtaxa); // using actual # terminal taxa as index
  InM.sdpd = vector(0, Intree.termtaxa);
  InM.mnnd = vector(0, Intree.termtaxa);
  InM.sdnnd = vector(0, Intree.termtaxa);

  if ( (Fm = fopen(meansfile, "r")) == NULL)
    {
      printf("Please generate the `means' file first\n");
      exit(EXIT_FAILURE);
    }
  
  while (fgets(line, 500, Fm) != NULL)
    {
      sscanf(line, "%d\t%10f\t%10f\t%10f\t%10f", \
	     &samplesizeI, &meanI, &sdI, &nnmeanI, &nnsdI); 
      // float width is excess
      linesread++;
      
      InM.mpd[samplesizeI] = meanI;
      InM.sdpd[samplesizeI] = sdI;      
      InM.mnnd[samplesizeI] = nnmeanI;
      InM.sdnnd[samplesizeI] = nnsdI;
    }
  
  fclose(Fm);
  return InM;
}

// =================================

void AgeNodes(phylo P)
{
  int i;
  int maxone = 0;
  int xnode;
  float longest = 0.0;
  float *sum;
  sum = vector(0,P.nnodes);
  
  for (i = 0; i < P.nnodes; i++)
    {
      xnode = i;
      sum[i] = 0.0;
      while (xnode != 0)
	{
	  sum[i] += P.bl[xnode];
	  xnode = P.up[xnode];
	}
      if (sum[i] > longest) 
	{
	  longest = sum[i];
	  maxone = i;
	}
    }
  for (i = 0; i < P.nnodes; i++) P.age[i] = sum[maxone] - sum[i];
  free_vector(sum, 0,P.nnodes);
}

// -------------- WriteNexus ------------------
void WriteNexus(phylo P[], int ntree, sample S, int nsamp, traits T, int ntrf)
{
  // Mesquite style!
  time_t rawtime;
  int i, j, q, k, x, pass, present;
  int makedisc, makecont;
  float abnd;
  int nterm = 0;
  phylo WN[ntree];
  char tmp[MAXTAXONLENGTH+10];

  for (i = 0; i < ntree; i++)
    {
      WN[i] = P[i]; // inefficient to make copy so much, but need to 
      // to create third dimension of taxon array

      // reassign the pointer to a new space - free this!
      WN[i].taxon = cmatrix(0, P[0].nnodes-1, 0, MAXTAXONLENGTH+10);
    }

  // determine number of terminal taxa - assume all trees contain same taxa
  for (i = 0; i < P[0].nnodes; i++)
    {
      if (P[0].noat[i] == 0) nterm++;
    }

  time ( &rawtime );
  strncpy(tmp , ctime(&rawtime), 24);

  printf("#NEXUS\n[output from phylocom, written %s]\n\n", tmp );
  printf("BEGIN TAXA;\n");
  if (TreeView == 0) printf("TITLE Phylocom_Phylogeny_Taxa;\n"); // Needed for correct Mesquite grammar, but V1.1 busted! Will not read interior names correctly.
  printf("\tDIMENSIONS NTAX=%d;\n\tTAXLABELS\n\t", nterm);

  for (i = 0; i < P[0].nnodes; i++)
    {
      if (P[0].noat[i] == 0) printf("   %s", P[0].taxon[i]);
    }
  printf(";\nEND;\n\n");

  if (nsamp > 0)
    {
      // Characters
      printf("BEGIN CHARACTERS;\n\tTITLE  Phylocom_Presence_in_Sample;\n\tDIMENSIONS NCHAR=%d;\n\tFORMAT DATATYPE = STANDARD GAP = - MISSING = ? SYMBOLS = \"  0 1\";\n", S.nsamples);

      printf("\tCHARSTATELABELS\n\t\t");
      printf("%d %s", 1, S.pname[0]);
      for (i = 1; i < S.nsamples; i++)
	{
	  printf(", %d %s", i+1, S.pname[i]);
	}
      printf(";\n\tMATRIX\n");
      for (i = 0; i < P[0].nnodes; i++)
	{
	  if (P[0].noat[i] == 0)
	    {
	      printf("\t%s\t" , P[0].taxon[i]);
	      for (j = 0; j < S.nsamples; j++)
		{
		  present = 0;
		  for (k = 0; k < S.srec[j]; k++)
		    {
		      if (strcmp(S.taxa[S.id[j][k]], P[0].taxon[i]) == 0) present = 1;
		      
		    }
		  printf("%d", present);
		}
	      printf("\n");
	    }
	}
      printf(";\nEND;\n\n");

      // Abundances as continuous
      printf("BEGIN CHARACTERS;\n\tTITLE  Phylocom_Abundance_in_Sample;\n\tDIMENSIONS NCHAR=%d;\n\tFORMAT DATATYPE = CONTINUOUS GAP = - MISSING = ?;\n", S.nsamples);

      printf("\tCHARSTATELABELS\n\t\t");
      printf("%d %s", 1, S.pname[0]);
      for (i = 1; i < S.nsamples; i++)
	{
	  printf(", %d %s", i+1, S.pname[i]);
	}
      printf(";\n\tMATRIX\n");
      for (i = 0; i < P[0].nnodes; i++)
	{
	  if (P[0].noat[i] == 0)
	    {
	      printf("\t%s\t" , P[0].taxon[i]);
	      for (j = 0; j < S.nsamples; j++)
		{
		  abnd = 0.0;
		  for (k = 0; k < S.srec[j]; k++)
		    {
		      if (strcmp(S.taxa[S.id[j][k]], P[0].taxon[i]) == 0) abnd = (float) S.abund[j][k];
		      
		    }
		  printf("  %f", abnd);
		}
	      printf("\n");
	    }
	}
      printf(";\nEND;\n\n");
    }

  if (ntrf > 0)
    {
      makedisc = 0;
      makecont = 0;
      pass = 0;

      for (i = 0; i < T.ntraits; i++)
	{
	  if ((T.type[i] == 0) || (T.type[i] == 1) || (T.type[i] == 2)) makedisc++;
	   
	  if (T.type[i] == 3) makecont++;
	}
      if (makedisc > 0)
	{
	  // Discrete Traits
	  printf("BEGIN CHARACTERS;\n\tTITLE  Phylocom_Discrete_Traits;\n\tDIMENSIONS NCHAR=%d;\n\tFORMAT DATATYPE = STANDARD GAP = - MISSING = ? SYMBOLS = \"  0 1 2 3 4 5 6 7 8 9\";\n", makedisc);
	  x = 1;
	  printf("\tCHARSTATELABELS\n\t\t");
	  // first one
	  for (i = 0; i < T.ntraits; i++)
	    {
	      if ((T.type[i] == 0) || (T.type[i] == 1) || (T.type[i] == 2))
		{
		  printf("%d %s", x, T.trname[i]);
		  x++;
		  pass = i;
		  break;
		}
	    }
	  for (i = pass+1; i < T.ntraits; i++)
	    {
	      if ((T.type[i] == 0) || (T.type[i] == 1) || (T.type[i] == 2))
		{
		  printf(", %d %s", x, T.trname[i]);
		  x++;
		}
	    }
	  
	  printf(";\n\tMATRIX\n");
	  for (i = 0; i < T.ntaxa; i++)
	    {
	      printf("\t%s\t" , T.taxon[i]);
	      for (j = 0; j < T.ntraits; j++)
		{
		  if ((T.type[j] == 0) || (T.type[j] == 1) || (T.type[j] == 2))
		    {
		      printf("%d", (int) T.tr[i][j]);
		    }
		}
	      printf("\n");
	    }
	  printf(";\nEND;\n\n");
	}


      if (makecont > 0)
	{
	  // Continous Traits
	  printf("BEGIN CHARACTERS;\n\tTITLE  Phylocom_Continuous_Traits;\n\tDIMENSIONS NCHAR=%d;\n\tFORMAT DATATYPE = CONTINUOUS GAP = - MISSING = ?;\n", makecont);
	  x=1;
	  printf("\tCHARSTATELABELS\n\t\t");
	  // first one
	  for (i = 0; i < T.ntraits; i++)
	    {
	      if (T.type[i] == 3)
		{
		  printf("%d %s", x, T.trname[i]);
		  x++;
		  pass = i;
		  break;
		}
	    }
	  for (i = pass+1; i < T.ntraits; i++)
	    {
	      if (T.type[i] == 3)
		{
		  printf(", %d %s", x, T.trname[i]);
		  x++;
		}
	    }
	  
	  printf(";\n\tMATRIX\n");
	  for (i = 0; i < T.ntaxa; i++)
	    {
	      printf("\t%s\t" , T.taxon[i]);
	      for (j = 0; j < T.ntraits; j++)
		{
		  if (T.type[j] == 3)
		    {
		      printf("  %f", T.tr[i][j]);
		    }
		}
	      printf("\n");
	    }
	  printf(";\nEND;\n\n");
	}
    }

  printf("BEGIN TREES;\n");
  if (TreeView == 0) printf("\tTITLE Phylocom_Phylogenies;\n\tLINK Taxa = Phylocom_Phylogeny_Taxa;\n"); // Ditto!
  printf("\tTRANSLATE\n\t");


  for (q = 0; q < ntree; q++)
    {
      j = 0;      
      for (i = 0; i < P[0].nnodes; i++)
	{
	  if (P[0].noat[i] == 0) 
	    {
	      j++;
	      if (q == 0)
		{
		  if (i == P[0].nnodes-1) printf(" %d %s;\n", j, P[0].taxon[i]);
		  else printf(" %d %s,", j, P[0].taxon[i]);
		}
	      sprintf(tmp, "%d", j);
	      strcpy(WN[q].taxon[i], tmp);
	    }
	  else if ((strcmp(P[q].taxon[i], "") != 0) && \
		   (strcmp(P[q].taxon[i], ".") != 0))
	    {
	      strcpy(WN[q].taxon[i], "'");
	      strcat(WN[q].taxon[i], P[q].taxon[i]);
	      strcat(WN[q].taxon[i], "'");
	    }
	  else strcpy(WN[q].taxon[i], "");
	  // test if (strcmp(WN[q].notes[i], "") != 0) printf("%s\n", WN[q].notes[i]);
	}
    }
  for (q = 0; q < ntree; q++)
    {

      printf("\tTREE %s = ", WN[q].phyname);
      Fy2newRec(WN[q]);
      free_cmatrix(WN[q].taxon, 0, P[0].nnodes-1, 0, MAXTAXONLENGTH+10);
    }
  printf("END;\n");

  printf("\nBEGIN PHYLOCOM;\n\tTITLE Phylocom_Main;\n\tDATA\n");
  for (i = 0; i < S.nsamples; i++)
    {
      for (j = 0; j < S.srec[i]; j++)
	{
	  printf("%s\t%d\t%s\n", S.pname[i], S.abund[i][j], S.taxa[S.id[i][j]]);
	}
    }
  printf(";\nEND;\n");

  //free_cmatrix(WN.taxon, 0, P.nnodes-1, 0, MAXTAXONLENGTH+10);

}

// -------------- NewickToNexus ------------------
void NewickToNexus(phylo P)
{

  //SWK - A Nexus wrapper takes a Newick file and makes it Mesquite readable
  time_t rawtime;
  int i;
  char tmp[24];

  time ( &rawtime );
  strncpy(tmp , ctime(&rawtime), 24);

  printf("#NEXUS\n[output from phylocom, written %s]\n\n", tmp );
  printf("BEGIN TAXA;\n");
  printf("TITLE Phylocom_Phylogeny_Taxa;\n"); // Needed for correct Mesquite grammar, but V1.1 busted! Will not read interior names correctly.
  printf("\tDIMENSIONS NTAX=%d;\n\tTAXLABELS\n\t", P.termtaxa);

  for (i = 0; i < P.nnodes; i++) {
    if (P.noat[i] == 0) printf("   %s", P.taxon[i]);
  }
 
  printf(";\nEND;\n\n");

  printf("BEGIN TREES;\n");
  printf("\tTITLE Phylocom_Phylogenies;\n\tLINK Taxa = Phylocom_Phylogeny_Taxa;\n"); // Ditto!
  //printf("\tTRANSLATE\n\t");
  //for (i = 0; i < P.ntaxa; i++) printf(" %d %s;\n", i, P.taxon[P.t2n[i]]);

  printf("\tTREE %s = ", P.phyname);
  //Is there a way to exclude node labels? Might be useful... see fy2new.c
  Fy2newRec(P);
  printf("END;\n");

}


// ---------- PRINT WELCOME --------------------------------------------------
void PrintWelcome()
{
  printf("\n  ==========================================================\n\n");
  printf("                       P H Y L O C O M\n");
  printf("        PHYLOGENETIC ANALYSIS OF ECOLOGICAL COMMUNITIES\n");
  printf("          AND SPECIES TRAITS, WITH PHYLOGENETIC TOOLS\n");
  printf("                Cam Webb <cwebb@oeb.harvard.edu>,\n");
  printf("                David Ackerly <dackerly@berkeley.edu>,\n");
  printf("                Steve Kembel <skembel@berkeley.edu>\n");
  printf("\n  ==========================================================\n\n");

  printf("  Copyright (c) 2008 Webb, Ackerly and Kembel\n");
  printf("  This program comes with ABSOLUTELY NO WARRANTY.  This is free software,\n"); 
  printf("  and you are welcome to redistribute it under certain conditions;\n");
  printf("  type `phylocom license' for warranty and license information.\n\n");

  printf("Usage: phylocom method [ options ]\n\n");
  printf("Phylogenetic community structure (intrasample):\n");
  printf("  comstruct: Calculates mean phylogenetic distance and mean minimum\n");
  printf("            phylogenetic distance for each plot and all plots, and\n");
  printf("            compares them to [%d] random runs with matrices created by\n",RUNS);
  printf("            the null model specified by the -m switch.\n");
  printf("            Alter # runs with a second argument.\n");
  printf("  swap:     Swaps the input matrix %d times using the Independent Swap\n",SWAPS);
  printf("            Algorithm (independent swap) and outputs swapped matrix console.\n");
  printf("  pd:       Faith's Phylogenetic diversity; proportion of total branch length.\n");
  printf("  ltt:      Prints the number of lineages at fixed divs of total time\n");
  printf("            proportional to no. of tip taxa.  Phylo must be\n");
  printf("            ultrametric. [ltteld compares with random samples.]\n");
  printf("  ltteld:   Calculates Epoch-specific Lineage Diversity\n");
  printf("  nodesig:  Tests each node for overabundance of terminal taxa.\n");
  printf("\n");
  printf("Phylogenetic community structure (intersample):\n");
  printf("  comdist:  Outputs the pairwise distance matrix between plots, based on mean\n");
  printf("            phylogenetic distance of all possible pairs of taxa in one plot\n");
  printf("            to the taxa in the other. [comdistnn uses nearest neighbor method.]\n");
  printf("  icomdist: Gives the mean distance from each taxon to taxa in other plots.\n");
  printf("  rao     : Calculates Rao's entropy measure\n");
  printf("\n");
  printf("Community trait structure:\n");
  printf("  comtrait: Measures of trait value clustering and evenness within samples\n");
  printf("\n");
  printf("Phylogenetic trait analysis:\n");
  printf("  aot:      Trait analysis algorithms of David Ackerly\n");
  printf("  aotn:     Trait analysis algorithms (Nexus output)\n");
  printf("  aotf:     Trait analysis algorithms (full output)\n");
  printf("\n");
  printf("Misc Phylogeny Tools:\n");  
  printf("  phydist:  Calculates the simple pairwise matrix of phylogenetic distances.\n");
  printf("  phyvar:   Output phylogenetic variance-covariance matrix.\n");
  printf("  agenode:  Output the ages of each node, calculated from BLs.\n");
  printf("  ageterm:  Output the stem age of terminal taxon.\n");
  printf("  makenex:  Convert all data files into a Mesquite-readable Nexus\n");
  printf("  naf:      Convert all data files into a node-as-factor table\n");
  printf("  new2nex:  Convert the phylo file into a Mesquite-readable Nexus\n");
  printf("  bladj:    Branch length adjust alg.  Needs an 'ages' file.\n");
  printf("  comnode:  Find common nodes between two trees of different size.\n");
  printf("  sampleprune: Prunes the phylo for each list of taxa in sample.\n");
  printf("  rndprune: Prunes the phylo [-r n] times, using [-p n] terminals.\n");
  printf("  cleanphy: Prunes `empty' internal nodes (e.g., output from phylomatic)\n");
  printf("  version:  Print software version.\n");
  printf("  help:     This file\n");
  printf("\n");
  printf("Options:\n");
  printf("  -d       Debug mode.\n");
  printf("  -r INT    Number of randomizations to use [%d].\n",RUNS);
  printf("  -m INT    Swap method to use with comstruct/swap: [%d].\n",SWAPMETHOD);
  printf("              0 = Shuffle phylogeny taxa labels.\n");
  printf("              1 = Sample taxa become random draws from sample pool.\n");
  printf("              2 = Sample taxa become random draws from phylogeny pool.\n");	
  printf("              3 = Shuffle sample using independent swap.\n");
  printf("  -w INT    Number of swaps to use per run (for independent swap) [%d].\n",SWAPS);
  printf("  -a INT    Output only specific result lines\n");
  printf("              %d = Output all\n", 0);
  printf("              %d = Output only the stELD line per sample\n", OUT_ANR_E);
  printf("              %d = Output only the LOG Nobs line per sample\n", OUT_LOGNOBS_E);
  printf("              %d = Output only the LOG Nrnd line per sample\n", OUT_LOGNRND_E);
  printf("              %d = Output only the STDEV LOG Nrnd line per sample\n", OUT_STDEVLOGNRND_E);
  printf("              %d = Output all ELD parameters\n", OUT_DIFFCOEFFS_E);
  printf("  -f filename   Use this file as the phylogeny file [phylo].\n");
  printf("  -s filename   Use this file as the sample file [sample].\n");
  printf("  -t filename   Use this file as the traits file [traits].\n");
  printf("  -z filename   Use this file as the ltteld points file [lttrsplits].\n");
  printf("  -p INT    Number of taxa to include in pruned tree [5].\n");
  printf("  -n        Add default internal node names.\n");
  printf("  -e        Suppress output of branch lengths.\n");  
  printf("  -v        Verbose output of raw MPD/MNTD values (comstruct).\n");
  printf("  -ol       Output raw LTTELD data to the file ltt_raw_<samplename>\n");
  printf("  -i        Output variances insttead of standard deviations in ltteld\n");
  printf("\n");
}

// ---------- Append Note -----------------------------------------------------
//void AppendNote(phylo Intree, char **note)
//{
//  int i;
//  // note = cmatrix(0, Intree.nnodes-1, 0, 100);
//  // printf("%s\n\n\n", note[0]);
//
//  for (i = 0; i < Intree.nnodes; i++)
//    {
//      printf("%d\t%d\t%d\t%s\t%d\t%f\t%s%s\n", 
//	     i, Intree.up[i], Intree.noat[i], Intree.downbr[i], 
//	     Intree.depth[i], Intree.bl[i], note[i], Intree.taxon[i]);
//    }
//  printf("\n\n");
//  // free_cmatrix(note, 0, Intree.nnodes-1, 0, 100);
//
//}



// ---------------- ReadLttrIntv ---------------------------
lttr_intvs ReadLttrIntv(char lttr_intvfile[50]){
  
  struct lttr_intvs bounds;
  int lineending;
  char line[500];
  float fBound = 0.0;
  float fMaxBoundSoFar = 0.0;
  int nscanned = 0;
  int iLineNo = 0;

  // pre-pre-read
  lineending = whatnewline(lttr_intvfile);

  // preread
  //READ PROPERLY

  if ((Fo = fopen(lttr_intvfile, "r")) == NULL){
    printf("Cannot open LTTELD datapoints file.");
    exit(0);
  }

  iLineNo = 0;
  while (myfgets(line, 500, Fo, lineending) != NULL){
    if (Debug)
      {
        printf("line = '%s'", line);
      }
    nscanned = sscanf(line, "%f", &fBound);
    if (nscanned < 1){
      printf("Error (LTTELD) in intervals file %s, line %d. Not a floating point number.\n", lttr_intvfile, iLineNo+1);
      exit(0);
      
    }
    iLineNo++;
  }
  fclose(Fo);

  bounds.nintvs = iLineNo;
  bounds.age_bounds = NULL;

  //Create dynamic arrays

  bounds.age_bounds = vector(0, bounds.nintvs-1);

  //Read
  
  if ((Fo = fopen(lttr_intvfile, "r")) == NULL){
    printf("Cannot open LTTELD datapoints file.");
    exit(0);
  }


  iLineNo = 0; //reuse counter
  fMaxBoundSoFar = 0;
  while (myfgets(line, 500, Fo, lineending) != NULL){
    sscanf(line, "%f", &fBound);
    if (fBound < fMaxBoundSoFar){
      //error
      printf("Error (LTTELD): the value %f on line %d in file %s is lower than the highest read so far\n", fBound, iLineNo, lttr_intvfile);
      printf("The values must be in strict ascendant order.\n");
      exit(0);
    }
    bounds.age_bounds[iLineNo] = fBound;
    fMaxBoundSoFar = fBound;
    iLineNo++;
  }
  fclose(Fo);
  

  return bounds;
}



// ---------------- ReadSample ---------------------------
sample ReadSample(char samplefile[50])
{
  struct sample InSample;

  char line[500];
  int abundI;
  char preAbund[50];
  int i, j, z, found;
  char plotname[MAXPLOTLENGTH+1];
  char lastplot[MAXPLOTLENGTH+1];
  char taxaI[MAXTAXONLENGTH+1];
  int lineending;

  strcpy(lastplot, "");
  i = 0;
  InSample.nrec = 0;
  InSample.nsamples = 0;
  InSample.maxrec = 0;
  InSample.ntaxa = 0;
  InSample.totabund = 0;

  // pre-pre-read
  lineending = whatnewline(samplefile);

  // preread:
  // READ PROPERLY
  if ((Ft = fopen(samplefile, "r")) == NULL)
    {
      printf("Cannot open sample file\n"); 
      exit(0);
    }
	
  while (myfgets(line, 500, Ft, lineending) != NULL) 
    // OK to lv length undynamic
    {
      sscanf(line, "%s %s %s", plotname, preAbund, taxaI); // string
      // Test for non-integer abundances:
      z = 0;
      while (preAbund[z] != '\0')
	{
	  // printf("%d %d\n", z, preAbund[z]);
	  if (!isdigit(preAbund[z]))
	    {
	      fprintf(stderr,"Error: sample abundances must be integers\n");
	      exit(1);
	    }
	  z++;
	}
      //if (strlen(taxaI) == 0){
		//printf("sample record number %d empty string", InSample.nrec);
	  //}
      InSample.nrec++;

      // if a new plot:
      if (strcmp(lastplot, plotname) != 0)
	{
	  InSample.nsamples ++;
	  i = 0;
	  strcpy(lastplot, plotname);
	}
      i++;
      if (InSample.maxrec < i) InSample.maxrec = i;
	  
    }
  fclose(Ft);

  //printf("nrec = %d; nsamples = %d; maxrec = %d\n", InSample.nrec, InSample.nsamples, InSample.maxrec);

  // Dimension structure

  InSample.srec = ivector(0, InSample.nsamples - 1);
  InSample.irec = ivector(0, InSample.nsamples - 1);
  InSample.id =   imatrix(0, InSample.nsamples - 1, 0, InSample.maxrec);
  InSample.abund =   imatrix(0, InSample.nsamples - 1, 0, InSample.maxrec);
  InSample.pabund =   matrix(0, InSample.nsamples - 1, 0, InSample.maxrec);
  InSample.taxa  =   cmatrix(0, InSample.nrec-1, 0, MAXTAXONLENGTH); 
  InSample.pname  = cmatrix(0, InSample.nsamples - 1, 0, MAXPLOTLENGTH);

  // Read file for dynamic structure
  strcpy(lastplot, "");
  i = 0;
  InSample.nsamples = 0; // recycling as counter

  Ft = fopen(SampleFile, "r");
  while (myfgets(line, 500, Ft, lineending) != NULL)
    {
      sscanf(line, "%s %d %s", plotname, &abundI, taxaI);

      // do this on the first line of each new plot
      if (strcmp(lastplot, plotname) != 0)
	{
	  strcpy(InSample.pname[InSample.nsamples], plotname);
	  InSample.nsamples ++;
	  i = 0;
	  strcpy(lastplot, plotname);
	  InSample.srec[InSample.nsamples-1]=0;
	}
      // recs per sample counter
      InSample.srec[InSample.nsamples-1]++;

      // set abundance
      InSample.abund[InSample.nsamples-1][i] = abundI;
      
      // set taxon
      // all other cases than first
      found = 0;
      for (j = 0; j < InSample.ntaxa; j++)
	{
	  if (strcmp(InSample.taxa[j], taxaI) == 0)
	    {
	      InSample.id[InSample.nsamples-1][i] = j;
	      found = 1;
	      break;
	    }
	}
      // new taxon found - executed on first line 
      if(found == 0)
	{
	  strcpy(InSample.taxa[InSample.ntaxa], taxaI);
	  InSample.id[InSample.nsamples-1][i] = InSample.ntaxa;
	  InSample.ntaxa++;
	}

      i++;
    }
  fclose(Ft);
  
  //Calculate sample and species total abundances and frequency
  InSample.sppabund = lvector(0, InSample.ntaxa - 1);
  InSample.psppabund = vector(0, InSample.ntaxa -1 );
  InSample.sppfreq = lvector(0, InSample.ntaxa - 1);
  InSample.psppfreq = vector(0, InSample.ntaxa -1 );
  
  // clear by sample, insample
  for (i = 0; i < InSample.nsamples; i++)
    {
      InSample.irec[i] = 0;
      for (j = 0; j < InSample.srec[i];j++)
	{
	  InSample.irec[i] += InSample.abund[i][j];
	  InSample.sppabund[InSample.id[i][j]] = 0;
	}
    }
  
  // clear by taxonNo in all samples
  for (i = 0; i < InSample.ntaxa; i++)
    {
      InSample.sppfreq[i] = 0;
    }

  // calculate
  for (i = 0; i < InSample.nsamples; i++)
    {
      for (j = 0; j < InSample.srec[i];j++)
	{
	  InSample.pabund[i][j] = (float) InSample.abund[i][j] / (float) InSample.irec[i];
	  InSample.sppabund[InSample.id[i][j]] += InSample.abund[i][j];
	  InSample.totabund += InSample.abund[i][j];
	  InSample.sppfreq[InSample.id[i][j]]++;
	  // printf("%d %d %d %s %d\n", i, j, InSample.id[i][j], InSample.taxa[InSample.id[i][j]], (int) InSample.sppabund[InSample.id[i][j]]);
	}
    }
  
  for (i = 0; i < InSample.ntaxa; i++)
  {
	  InSample.psppabund[i] = (float) InSample.sppabund[i] / (float) InSample.totabund;
	  InSample.psppfreq[i] = (float) InSample.sppfreq[i] / (float) InSample.nsamples;
  }
  return InSample;
}

// ----------- ReadTraits -------------
traits ReadTraits(char traitfile[50])
{
  struct traits C;
  FILE *Cread;
  char line[MAXTRAITLINE];
  int i, j, waitingforspace;
  int nline, words;
  int extra = 0;
  char word[(int) (MAXTRAITLINE / 2)][(int) (MAXTRAITLINE / 2)];
  char tmp[MAXTAXONLENGTH + 6];
  int lineending;

  C.ntaxa = 0;
  C.ntraits = 0;

  // pre-pre-read
  lineending = whatnewline(traitfile);
  
  // preread:
  if ((Cread = fopen(traitfile, "r")) == NULL)
    {
      printf("Cannot open traits file\n"); 
      exit(0);
    }
  nline = 0;
  words= 0;

  while (myfgets(line, 1000, Cread, lineending) != NULL)
    {
      // parse the first line
      // should be able to do this with iterative scanf, but I couldn't
      // work it out - this is clunky!
      if (nline < 2)
	{
	  waitingforspace = 1;
	  strcpy(word[0], "");
	  words= 0;

	  for (i = 0; i < strlen(line); i++)
	    {
	      if (isspace(line[i]) && (waitingforspace))
		{
		  words++;
		  strcpy(word[words], "");
		  waitingforspace = 0;
		}
	      else if (isspace(line[i])) // 2nd+ space in a series of white 
		{
		}
	      else
		{
		  strncat(word[words] , &line[i], 1);
		  waitingforspace = 1;
		}
	    }

	  if(nline==1)
	    {
	      // to see if there is a names line
	      if (strcmp(word[0], "name") ==0) extra = 1;
	    }
	}
      
      nline++;
    }
  fclose(Cread);

  C.ntraits = words-1;
  C.type = ivector(0, C.ntraits-1);
  C.ntaxa = nline -1 -extra;
  C.taxon = cmatrix(0 , C.ntaxa-1, 0 , MAXTAXONLENGTH);
  C.tr = matrix(0, C.ntaxa-1, 0, C.ntraits-1);
  C.trname = cmatrix(0, C.ntraits-1, 0, MAXTAXONLENGTH);
  strcpy(C.trname[0], "nonames");

  // printf("traits: %d\ttaxa: %d\n", C.ntraits, C.ntaxa);

  // real read
  Cread = fopen(traitfile, "r");
  nline = 0;

  while (myfgets(line, 1000, Cread, lineending) != NULL)
    {
      words= 0;
      strcpy(word[0], "");
      waitingforspace = 1;

      // parse into words, every line
      for (i = 0; i < strlen(line); i++)
	{
	  if (isspace(line[i]) && (waitingforspace))
	    {
	      words++;
	      strcpy(word[words], "");
	      waitingforspace = 0;
	    }
	  else if (isspace(line[i]))
	    {
	    }
	  else
	    {
	      strncat(word[words] , &line[i], 1);
	      waitingforspace = 1;
	    }
	}

      if(nline==0)
	{
	  // now we know how many traits
	  
	  for (j = 1; j< words; j++)
	    {
	      C.type[j-1] = atoi(word[j]);
	    }
	}

      else // lines 2+
	{
	  // check (again) for trait names
	  if ((nline == 1) && (extra == 1))
	    {
	      for (j = 1; j < words; j++)
		{
		  strcpy(C.trname[j-1], word[j]);
		}
	    }
	  else // either line 3+ or line 2 with no trait names 
	    {
	      // assign taxon name
	      strcpy(C.taxon[(nline-1)-extra], word[0]);
	      for (j = 1; j < words; j++)
		{
		  C.tr[nline-1-extra][j-1] = (float) atof(word[j]);
		  // printf("%s\t%d\t%f\n", C.taxon[nline-1], j-1, (float) atof(word[j]));
		}
	    }
	}

      nline++;
    }
  
  fclose(Cread);

  // assign taxon names if there are none
  if (extra != 1) // if no trait names given
    {
      for (j = 0; j < C.ntraits; j++)
	{
	  sprintf(tmp, "trait_%d", j+1);
	  strcpy(C.trname[j], tmp);
	}
    }

  for (i = 0; i < C.ntraits; i++)
    {
      if (C.type[i] == 1)
	{
	  printf("Char type = 1; multistate characters not implemented yet\n");
	  exit(EXIT_FAILURE);
	}
    }
  
  return C;

}


// --------------- ReadPhylogeny -----------------------
phylo ReadPhylogeny(char PhyloFile[50])
{

      //Read the phylogeny in Newick format
      return New2fy(PhyloFile);
}

// -------------- MakeUpPass Order -----------------
void MakeUpPassOrder(phylo P)
{
  int i, order, depth;

  //P.upo = ivector(0, P.nnodes-1);
  //P.ntip = ivector(0, P.nnodes-1);
  //P.nint = ivector(0, P.nnodes-1);
  order = 0;
  P.maxDepth = 0;

  for (i = 0; i < P.nnodes; i++)
    {
      if (P.maxDepth <= P.depth[i]) P.maxDepth = P.depth[i];
      P.ntip[i] = 0;
      P.nint[i] = 0;
    }

  for (depth = P.maxDepth; depth >= 0; depth--)
    { 
      for (i = 0; i < P.nnodes; i++)
	{
	  // Terminals: Pass up tree
	  if ((P.depth[i] == depth) && (P.noat[i] == 0))
	    {
	      P.ntip[i]++;
	      P.ntip[P.up[i]]++;
	      P.upo[order] = i;
	      order++;
	    }
	  // Internals
	  if ((P.depth[i] == depth) && (P.noat[i] != 0))
	    {
	      P.nint[i]++;
	      P.upo[order] = i;
	      order++;
	      if (i != 0)
		{
		  P.nint[P.up[i]] += P.nint[i];
		  P.ntip[P.up[i]] += P.ntip[i];
		}
	    }
	}
    }

  if (FALSE)
    {
      for (i = 0; i < P.nnodes; i++)
	{
	  printf("%d\t%d\n", i, P.upo[i]);
	}
    }
}

// -------------- AssignNodeLists -----------------
void AssignNodeLists(phylo P)
{
  // Assigns a list of terminal and internal node numbers to vectors
  // associated with each node. These can then be used to slice off
  // trait and contrasts at each node for running standard stats
  int i,j,n,node;
  int *tipin;
  int *intin;


  //P.tiplist = imatrix(0,P.nnodes-1,0,P.nnodes-1);
  //P.intlist = imatrix(0,P.nnodes-1,0,P.nnodes-1);

  tipin = ivector(0, P.nnodes-1);
  intin = ivector(0,P.nnodes-1);

  for (node = 0; node < P.nnodes; node++)
    {
      tipin[node] = 0;
      intin[node] = 0;
    } 

  for (n = 0; n < P.nnodes; n++)
    {
      node = P.upo[n];
      // printf("%d\t%d\n", n, node);
      // Terminals: Pass up tree
      if (P.noat[node] == 0)
	{
	  P.tiplist[node][0] = node;
	  tipin[node]++;
	  P.tiplist[P.up[node]][tipin[P.up[node]]] = node;
	  tipin[P.up[node]]++;
	}

      // Internals
      if (P.noat[node] != 0)
	{
	  P.intlist[node][intin[node]] = node;
	  intin[node]++;
	  if (node != 0)
	    {
	      for (i = tipin[P.up[node]]; 
		   i < tipin[P.up[node]] + P.ntip[node]; i++)
		{
		  j = i - tipin[P.up[node]];
		  P.tiplist[P.up[node]][i] = 
		    P.tiplist[node][j];
		}
	      tipin[P.up[node]] += P.ntip[node];
		    
	      for (i = intin[P.up[node]]; 
		   i < intin[P.up[node]] + P.nint[node]; i++)
		{
		  j = i - intin[P.up[node]];
		  P.intlist[P.up[node]][i] = 
		    P.intlist[node][j];
		}
	      intin[P.up[node]] += P.nint[node];		   
	    }	
	}
    }

  if (FALSE) { //Debug printing
    for (node = 0; node < P.nnodes; node++)
      {
	printf("%d\t%d\t%d\n", node, P.ntip[node], P.nint[node]);
	printf("ints");
	for (i = 0; i < intin[node]; i++) printf("\t%d", P.intlist[node][i]);
	printf("\n");
	printf("tips");
        for (i = 0; i < tipin[node]; i++) printf("\t%d", P.tiplist[node][i]);
	printf("\n");
      } 
  }
  free_ivector(tipin, 0, P.nnodes-1);
  free_ivector(intin, 0, P.nnodes-1);

}

/* myfgets:  get at most n chars from iop */
//TODO: Fix the bug with non-UNIX endlines
char *myfgets(char *s, int n, FILE *iop, int le)
{
  register int c = 0;
  register char *cs;
  cs = s;
 
  // le == 1 => Mac Case
  if (le == 1)
    {
      while (--n > 0 && (c = getc(iop)) != EOF)
	if ((*cs++ = c) == '\r')
          {
            break;
          }
      if (*(cs-1) == '\r')
	{
	  *--cs = '\n'; // replace the \r with \n
	  *++cs = '\0';
	}
      else
	{
	  *cs = '\0';
	}
    }
  
  // le == 2 => Windows Case
  else if (le == 2)
    {
      while (--n > 0 && (c = getc(iop)) != EOF)
	if ((*cs++ = c) == '\r')
	  {
	    getc(iop); // read one extra char first
	    break;
	  }
      if (*(cs-1) == '\r')
	{
	  *--cs = '\n'; // replace the \r with \n
	  *++cs = '\0';
	}
      else
	{
	  *cs = '\0';
	}
    }
  
  // le == 0 => Unix Case
  else if (le == 0)
    {
      while (--n > 0 && (c = getc(iop)) != EOF)
	if ((*cs++ = c) == '\n')
          {
            break;
          }
      *cs = '\0'; 
    }

  return (c == EOF && cs == s) ? NULL : s;
}

int whatnewline(char *filename)
{
  FILE *fp;
  int c = 0;
  int w = 0;
  int m = 0;
  int u = 0;
  int rlast = 0;
  int nl;  // unix = 0; mac = 1; windows = 2

  // Unix:    \n   ASCII = 10
  // Mac:     \r   ASCII = 13
  // Windows: \r\n

  if ((fp = fopen(filename, "r")) == NULL) 
    {
      fprintf(stderr, "testeol: Can't open file: %s\n", filename);
      exit(1);
    } 
  else 
    { 
      while ((c = getc(fp)) != EOF)
	{
	  // if char is a \n ... 
	  if (c == 10)
	    {
	      // ... and the last was an \r : windows delimiter. \r\n\r\n OK
	      if (rlast) {w++; rlast=0;}
	      // ... if not, a unix delim.  \n\n\n\n OK
	      else {u++;}
	    }
	  // if char is an \r ...
	  else if (c == 13) 
	    { 
	      // if one of several in a row, add a mac
	      if (rlast) { m++; }
	      // if first one, store it, could be mac or windows
	      else { rlast = 1 ; }
	    }
	  // when we get to a non-delimeter, after a \r
	  else if (rlast) {m++; rlast = 0;}
	  
	  // putc(c, stdout); // just for show
	}

      fclose(fp);
      
      // in case the very last char was an \r 
      if (rlast) {m++; rlast = 0;}
      
      // printf("w:%d m:%d u:%d\n", w, m, u); 
      if ((w != 0) && (m == 0) && (u == 0)) nl = 2;
      else if ((w == 0) && (m != 0) && (u == 0)) nl = 1;
      else if ((w == 0) && (m == 0) && (u != 0)) nl = 0;
      else 
	{
	  fprintf(stderr, "Mixed delimiters in file: %s\n", filename); 
	  exit(1);
	}
      
      if (nl == 0) fprintf(stderr, "File `%s' has unix newlines\n", filename);
      if (nl == 1) fprintf(stderr, "File `%s' has mac newlines\n", filename);
      if (nl == 2) fprintf(stderr, "File `%s' has windows newlines\n", filename);

      return nl;
    }
}

void License()
{

printf("\n       === License and warranty information for phylocom ===\n");
printf("\n");
printf(" Copyright (c) 2008, Campbell O Webb, David D Ackerly, Steven W Kembel\n");
printf(" All rights reserved.\n");
printf("\n");
printf(" Redistribution and use in source and binary forms, with or without\n");
printf(" modification, are permitted provided that the following conditions\n");
printf(" are met:\n");
printf("\n");
printf("    * Redistributions of source code must retain the above copyright\n");
printf("      notice, this list of conditions and the following disclaimer.\n");
printf("\n");
printf("    * Redistributions in binary form must reproduce the above\n");
printf("      copyright notice, this list of conditions and the following\n");
printf("      disclaimer in the documentation and/or other materials provided\n");
printf("      with the distribution.\n");
printf("\n");
printf("    * The names of the authors may not be used to endorse or promote\n");
printf("      products derived from this software without specific prior\n");
printf("      written permission.\n");
printf("\n");
printf(" THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n");
printf(" \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n");
printf(" LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n");
printf(" A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n");
printf(" OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n");
printf(" SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n");
printf(" LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n");
printf(" DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n");
printf(" THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n");
printf(" (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n");
printf(" OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n");
printf("\n");
printf(" Please cite phylocom thus:\n");
printf("\n");
printf("   Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom:\n");
printf("   software for the analysis of phylogenetic community structure\n");
printf("   and trait evolution. Version 4.0.\n");
printf("   <http://www.phylodiversity.net/phylocom/>\n");
printf("\n");

}

// ========== LOG =============================================================

// $Log: io.c,v $
// Revision 3.43  2005/11/19 04:57:46  skembel
// phyvar: a phylogenetic variance-covariance function
//
// Revision 3.42  2005/09/28 06:05:46  skembel
// Phylogeny-Trait attachment moved out of trait file loading function
//
// Revision 3.41  2005/09/19 19:21:06  skembel
// -v (verbose) switch outputs raw MPD/MNTD values from comstruct
//
// Revision 3.40  2005/09/18 02:29:22  dackerly
// aot tips stats working, preliminary output function
//
// Revision 3.39  2005/09/15 07:30:31  dackerly
// aot development update
//
// Revision 3.38  2005/09/14 22:29:09  dackerly
// Functions MakeUpPassOrder and AssignNodeLists added to io.c
//
// Revision 3.37  2005/09/09 19:06:25  ctenolophon
// updated main revision number
//
// Revision 3.36  2005/09/09 19:04:31  ctenolophon
// updated main revision numberphylocom.h
//
// Revision 1.1.1.1  2005/09/09 18:35:24  ctenolophon
// Initial import.
//
// Revision 3.26  2005/08/18 03:37:55  cam
// NAF added
//
// Revision 3.25  2005/06/15 08:26:59  cam
// Changes associated with phylomatic
//
// Revision 3.24  2005/06/14 19:35:20  cam
// Added new, recursive Fy2new - yeah!  Works with big files
//
// Revision 3.23  2005/06/12 01:41:14  cam
// Fixed confusion over ntaxa or termtaxa.
//
// Revision 3.22  2005/06/11 18:45:59  cam
// *** empty log message ***
//
// Revision 3.21  2005/06/11 16:01:15  cam
// Big joining with Steve's code.
//
// Revision 3.20  2005/01/05 23:35:52  cam
// Added clust without means!
//
// Revision 3.19  2005/01/05 22:50:26  cam
// *** empty log message ***
//
// Revision 3.18  2004/08/28 14:59:57  cam
// minor
//
// Revision 3.17  2004/07/01 20:52:42  cam
// David's edits
//
// Revision 3.16  2004/05/24 18:03:28  cam
// added PD
//
// Revision 3.15  2004/05/17 21:48:46  cam
// *** empty log message ***
//
// Revision 3.14  2004/04/25 04:46:27  cam
// minor
//
// Revision 3.13  2004/04/12 15:01:19  cam
// minor fixes
//
// Revision 3.12  2004/04/09 08:05:47  cam
// Dynamic traits
//
// Revision 3.11  2004/04/05 16:27:44  cam
// Add dynamic sample and string plot names to all
//
// Revision 3.10  2004/04/04 18:02:50  cam
// Nexus out
//
// Revision 3.9  2004/04/01 23:01:11  cam
// minor fix in fy2new.c
//
// Revision 3.9  2004/03/31 22:17:10  cam
// dynamic sample!
//
// Revision 3.8  2004/03/30 03:26:56  cam
// dynamic version of comdist
//
// Revision 3.7  2004/03/29 15:10:52  cam
// minor changes
//
// Revision 3.6  2004/03/27 00:07:56  cam
// combined new2fy
//
// Revision 3.5  2004/02/01 22:00:14  cam
// More DDA fixes - moving to dynamic for all...
//
// Revision 3.4  2004/01/30 01:01:43  cam
// pre- big fix
//
// Revision 3.3  2003/12/18 18:17:00  cam
// New distributable
//
// Revision 3.2  2003/12/11 22:05:32  cam
// working
//
// Revision 3.1  2003/12/09 01:11:42  cam
// Major overhaul
//
// Revision 3.0  2003/12/05 20:03:40  cam
// Split into multiple source files.  Seems to all be working!
//

// ========== EOF =============================================================
