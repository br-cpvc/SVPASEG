// ******************************************************
// Functions for mixture model optimization using genetic algorithms
// (C) Jussi Tohka, 2005 - 2009 
// Laboratory of NeuroImaging, Department of Neurology, UCLA Medical School, CA, USA
// Institute of Signal Processing, Tampere University of Technology, Finland  
// Version 2.0.
// e-mail jussi.tohka@tut.fi
// *****************************************************
// ******************************************************
// Permission to use, copy, modify, and distribute this software 
// for any purpose and without fee is hereby
// granted, provided that the above copyright notice appears in all
// copies.  The author, Tampere University of Technology and University of
// California, Los Angeles make no representations
// about the suitability of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
// *****************************************************
// This file belongs to SVPASEG package. The SVPASEG implements 
// a flexible pipeline to tissue classify MR-images. The most important 
// references are:
// [1]  J. Tohka , E. Krestyannikov, I.D. Dinov , A. MacKenzie-Graham,
// D.W. Shattuck , U. Ruotsalainen, and A.W. Toga. 
// Genetic algorithms for finite mixture model based voxel classification
// in neuroimaging.  
// IEEE Transactions on Medical Imaging  , 26(5):696 - 711, 2007.
//
// [2] J. Tohka , I.D. Dinov, D.W. Shattuck, and A.W. Toga. Brain MRI 
// Tissue Classification Based on Local Markov Random Fields, 
// Magnetic Resonance Imaging , in press, 2009. 
// ********************************************************************

// -------------------------------------
// At the moment this module cannot handle correlation matrices
// At the moment only for 1-d distributions

#ifndef GAMIXTURE_H
#define GAMIXTURE_H

#include "analyze.h"
#include "atlasspec.h"
#include "nhmrf.h"
#include "parzen.h"
using namespace std;


// Default parameter values
#define DEFAULT_ALPHA           0.5
#define DEFAULT_POPSIZE         100
#define DEFAULT_XOVERRATE         1
#define DEFAULT_MAXGENERATIONS  500      
#define DEFAULT_TERMINATIONTHR  0.0005
#define DEFAULT_SORTPOP         1
#define DEFAULT_PARZENN         101
#define DEFAULT_PARZENSIGMA     1
#define DEFAULT_EQUALVAR        0
#define DEFAULT_RESTARTS        10
#define DEFAULT_DETERMINISTIC   1
#define DEFAULT_RANDOM          1


struct SortStruct {  // I know this is clumsy, but couldn't come up with a better alternative in ten minutes or so.... 
  float energy;
  int index;
};

struct Population {
  int size;
  int dim;            //at the moment this must be equal to 1
  int numberOfLabels; // label 0 is reserved for background, this should be counted
  int numberOfPveLabels;
  LabelType* labelTypes; // first labels should be non PVE and then the PVE labels follow
  float* lowLimit;       // gives the lower limit on the probability of the label
  float* upLimit;        // gives the upper limit on the probability of the label
  float** mixtures; 
  float* energies; 
};

struct Parameters {
  float alpha;
  int size;
  float terminationThr;
  float xoverRate; 
  int maxGenerations;
  int sortPop;
  int parzenN;  
  float parzenSigma;
  bool equalVar;
  int restarts;
  int deterministic;
  int random;
};



inline int cmpStruct(const void* a, const void* b)
{
  //  a = (SortStruct*) a;
  // b = (SortStruct*) b;
  if( (((SortStruct*) a)->energy) > (((SortStruct*) b)->energy))
    return(1);
  else
    return(-1);
}

inline void bubbleSort(float* individual , int nofMixtures, int componentSize)
{
  int i,j,k;
  float tmp;  

  for(i = 0;i < (nofMixtures - 2);i++) {
    for(j = 1; j < (nofMixtures - 1 - i);j++) {
      if(individual[(j + 1)*componentSize] < individual[j*componentSize]) {
        for(k = 0;k < componentSize;k++) {
          tmp = individual[j*componentSize + k];
          individual[j*componentSize + k] = individual[(j + 1)*componentSize + k];
          individual[(j + 1)*componentSize + k] = tmp;
	}  
      }
    }
  }
}


// this is only 1-d version
inline float gaGetProb(Population* pop, int indNumber, int classNumber) 
{ 
  int num = 3;  
  if(classNumber > (pop->numberOfLabels - pop->numberOfPveLabels - 1)) {
    return(pop->mixtures[indNumber][(pop->numberOfLabels - pop->numberOfPveLabels)*num 
				    + classNumber - (pop->numberOfLabels - pop->numberOfPveLabels)]);
  }
  else {
    return pop->mixtures[indNumber][classNumber*num + num - 1];
  }
}

// this is only 1-d version
inline float gaGetMu(Population* pop, int indNumber, int classNumber) 
{
  int num = 3;
  return pop->mixtures[indNumber][classNumber*num ];
}

// this is only 1-d version
inline float gaGetSigma2(Population* pop, int indNumber, int classNumber) 
{
  int num = 3;
  return pop->mixtures[indNumber][classNumber*num + 1];
}
// this is only 1-d version
inline void gaSetProb(Population* pop, int indNumber, int classNumber, float val) 
{ 
  int num = 3;  
  if(classNumber > (pop->numberOfLabels - pop->numberOfPveLabels - 1)) {
    pop->mixtures[indNumber][(pop->numberOfLabels - pop->numberOfPveLabels)*num + classNumber
        - (pop->numberOfLabels - pop->numberOfPveLabels)] = val;
  }
  else {
    pop->mixtures[indNumber][classNumber*num + num - 1] = val;
  }
}

// this is only 1-d version
inline void gaSetMu(Population* pop, int indNumber, int classNumber,float val) 
{
  int num = 3;
  pop->mixtures[indNumber][classNumber*num] = val;
}

// this is only 1-d version
inline void gaSetSigma2(Population* pop, int indNumber, int classNumber, float val) 
{
  int num = 3;
  pop->mixtures[indNumber][classNumber*num + 1] = val;
}

// copies "essential" parts of the population

inline void copyPartialPopulation(Population* sourcePop, Population* targetPop)
{
  int i,j,mixtureSize;

  targetPop->size = sourcePop->size;
  targetPop->dim  = sourcePop->dim;
  targetPop->numberOfLabels = sourcePop->numberOfLabels;
  targetPop->numberOfPveLabels = sourcePop->numberOfPveLabels;
  targetPop->labelTypes = NULL;
  targetPop->upLimit = NULL;
  targetPop->lowLimit = NULL;

  // compute the size of each mixture
  if (sourcePop->dim == 1) { 
    mixtureSize = 2*(sourcePop->dim)*(sourcePop->numberOfLabels - sourcePop->numberOfPveLabels) + sourcePop->numberOfLabels;
  }
  else {
    mixtureSize = (2*(sourcePop->dim) + (sourcePop->dim)*(sourcePop->dim))
                   *(sourcePop->numberOfLabels - sourcePop->numberOfPveLabels) + sourcePop->numberOfLabels;
  }
  // Allocate
  targetPop->mixtures = new float* [sourcePop->size];
  for(i = 0;i < sourcePop->size;i++) {
    targetPop->mixtures[i] = new float[mixtureSize];
  }
  targetPop->energies = new float[sourcePop->size];
  for(i = 0;i < sourcePop->size;i++) {
    targetPop->energies[i] = sourcePop->energies[i];
    for(j = 0;j < mixtureSize;j++) {
      targetPop->mixtures[i][j] = sourcePop->mixtures[i][j];
    }
  }
}


// frees memory of the population created by the partial copying

inline void freePartialPopulation(Population* pop) 
{
  int i;
  for(i = 0;i < pop->size;i++) {
    delete[] pop->mixtures[i];
  }
  delete[] pop->mixtures;
  delete[] pop->energies;
}

inline void freePopulation(Population* pop)
{
  int i;
  for(i = 0;i < pop->size;i++) {
    delete[] pop->mixtures[i];
  }
  delete[] pop->mixtures;
  delete[] pop->energies;
  delete[] pop->upLimit;
  delete[] pop->lowLimit;
  delete[] pop->labelTypes;
}

int parseParams(Parameters* param,int n,char** arguments); 
void displayParameterHelp();
int gaInitializePopulation(Population* pop,int size, int dim, int numberOfLabels,int numberOfPveLabels,
                           LabelType* labelTypes, float* lowLimit, float* upLimit,  float* randomFloats, bool equalVar = false);
void gaEvaluate(Population* pop,PdfEstimate* hatf);
void gaTournamentSelection(Population* selPop,Population* pop,int elitism, unsigned int* randomUnsignedInts);
void gaReorder(Population* pop); // Re-order individuals in the population according the fitness
void gaBLX(Population* pop,Population* selPop,float xoverRate, int elitism, float alpha, float* lowLimit, float* upLimit, unsigned int* randomUnsignedInts, float* randomFloats, bool equalVar = false);
void gaSortPopulation(Population* pop,int sortDim);
bool gaTerminate(Population* pop, float thr);

#endif
