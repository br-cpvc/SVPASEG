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



#include "gamixture.h"
#include <time.h>

// This function parses upperLimit and lowerLimit from the command line 
// input to the program. 

int parseParams(Parameters* param,int n,char** arguments) 
{
  // Initialize with the default values
  int i;  

  param->alpha = DEFAULT_ALPHA;
  param->size  = DEFAULT_POPSIZE;
  param->terminationThr = DEFAULT_TERMINATIONTHR;
  param->xoverRate = DEFAULT_XOVERRATE;
  param->maxGenerations = DEFAULT_MAXGENERATIONS;
  param->sortPop   = DEFAULT_SORTPOP;
  param->parzenN  = DEFAULT_PARZENN;
  param->parzenSigma = DEFAULT_PARZENSIGMA;
  param->equalVar = DEFAULT_EQUALVAR;
  param->restarts = DEFAULT_RESTARTS;
  param->deterministic = DEFAULT_DETERMINISTIC;
  param->random = DEFAULT_RANDOM;

  if( ((n - 1) % 2) > 0 ) return(1);
  for(i = 1;i < (n - 1)/2 + 1;i++) {
    if(!strcmp(arguments[2*i - 1],"-alpha"))
      param->alpha = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-size"))
      param->size = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-terminationthr"))
      param->terminationThr = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-xoverrate"))
      param->xoverRate = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-maxgen"))
      param->maxGenerations = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-sortpop"))
      param->sortPop = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-parzenn"))
      param->parzenN = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-parzensigma"))
      param->parzenSigma = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-equalvar"))
      param->equalVar = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-restarts"))  // Version 1.1
      param->restarts = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-deterministic"))
      param->deterministic = atoi(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-random"))
      param->random = atoi(arguments[2*i]);
  }
  return(0);
}

void displayParameterHelp()
{

}

int gaInitializePopulation(Population* pop,int size, int dim, int numberOfLabels,int numberOfPveLabels,
                           LabelType* labelTypes, float* lowLimit, float* upLimit, float* randomFloats, bool equalVar)
{
  int i,j;
  int mixtureSize;
  int bgVar;
  float tmpRand;  
  float probsum;

  pop->size = size;
  pop->dim = dim;
  pop->numberOfLabels = numberOfLabels;
  pop->numberOfPveLabels = numberOfPveLabels;
  pop->labelTypes = labelTypes;
  pop->lowLimit = lowLimit;
  pop->upLimit = upLimit;

  
  // compute the size of each mixture
  if (pop->dim == 1) { 
    mixtureSize = 2*(pop->dim)*(numberOfLabels - numberOfPveLabels) + numberOfLabels;
  }
  else {
    mixtureSize = (2*(pop->dim) + (pop->dim)*(pop->dim))*(numberOfLabels - numberOfPveLabels) + numberOfLabels;
  }
  // compute the number of background variables
  if (pop->dim == 1) {
    bgVar = 3;
  }
  else {
    bgVar = 2*dim + 1 + dim*dim;
  } 
  // Allocate
  pop->mixtures = new float* [pop->size];
  for(i = 0;i < size;i++) {
    pop->mixtures[i] = new float[mixtureSize];
  }
  pop->energies = new float[pop->size];
  
  for(i = 0;i < pop->size;i++) {
    gaSetMu(pop,i,0,0.0);    // background
    gaSetProb(pop,i,0,0.0);  // background
    for(j = bgVar; j < mixtureSize;j++) {
      tmpRand = randomFloats[i*mixtureSize+j];
      pop->mixtures[i][j] = tmpRand*(upLimit[j] - lowLimit[j]) + lowLimit[j];
    }
    gaSetSigma2(pop,i,0,gaGetSigma2(pop,i,1));  // background variance
    if(equalVar) {
      for(j = 2;j < (pop->numberOfLabels - pop->numberOfPveLabels);j++) {
        gaSetSigma2(pop,i,j,gaGetSigma2(pop,i,1));
      }
    }
    probsum = 0.0;
    for(j = 1;j < pop->numberOfLabels;j++) {
      probsum = probsum + gaGetProb(pop,i,j);
    }
    for(j = 1;j < pop->numberOfLabels;j++) {
      gaSetProb(pop,i,j,gaGetProb(pop,i,j) / probsum);
    }
  }
  return(0);
}


// computes the Kullback Leibler divergence between the mixture and the 
// given data density
// ASSUMES THAT X-AXIS OF THE DENSITY IS EQUALLY SPACED

void gaEvaluate(Population* pop,PdfEstimate* hatf)
{
  int i,j,k;
  float interval;
  float klAddEnergy;
  float mixtureVal;
  float tmp;
  int pveIntervals = 10;

  klAddEnergy = 0.0;
  interval = hatf->x[1] - hatf->x[0];
  for(i = 0;i < hatf->n ; i++) {
    if ( hatf->y[i] > 0 ) {
      klAddEnergy = klAddEnergy + hatf->y[i]*log( hatf->y[i] );
    }
  }
  klAddEnergy = interval*klAddEnergy;
  
  for(i = 0;i < pop->size;i++) {
    pop->energies[i] = 0.0;
    for(j = 0;j < hatf->n;j++) {
      mixtureVal = 0.0;
      for(k = 1;k < (pop->numberOfLabels - pop->numberOfPveLabels);k++) {
        mixtureVal = mixtureVal + gaGetProb(pop,i,k)*computeGaussianLikelihood(hatf->x[j],gaGetMu(pop,i,k),gaGetSigma2(pop,i,k));
       
      }
      
      for(k = (pop->numberOfLabels - pop->numberOfPveLabels);k < pop->numberOfLabels;k++) {
        mixtureVal = mixtureVal + gaGetProb(pop,i,k)*computeMarginalizedLikelihood(hatf->x[j],gaGetMu(pop,i,pop->labelTypes[k].mixed[0]), 
                                                                                   gaGetMu(pop,i,pop->labelTypes[k].mixed[1]), 
                                                                                   gaGetSigma2(pop,i,pop->labelTypes[k].mixed[0]), 
                                                                                   gaGetSigma2(pop,i,pop->labelTypes[k].mixed[1]),10 ,DEFAULT_PVESTART, DEFAULT_PVEEND);
      }
      if(mixtureVal == 0) {
        mixtureVal = 1/(interval*1000000);
      }
      pop->energies[i] = pop->energies[i] + hatf->y[j]* log(mixtureVal);
    }
    pop->energies[i] = pop->energies[i]*interval; 
    pop->energies[i] = klAddEnergy - pop->energies[i];
  }
  
}


// selects pop->size - elitism from the population by means of the tournament selection
// The tournament size is 2 in the current implementation.
// Assumes that the population has been sorted.

 void gaTournamentSelection(Population* selPop,Population* pop,int elitism, unsigned int* randomUnsignedInts)
{
  int i,j,mixtureSize,r,r1,r2;
 

  // compute the size of each mixture
  if (pop->dim == 1) { 
    mixtureSize = 2*(pop->dim)*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  else {
    mixtureSize = (2*(pop->dim) + (pop->dim)*(pop->dim))*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  // copyPartialPopulation(pop,selPop);
  for(i = 0;i < elitism;i++) {
    for(j = 0;j < mixtureSize;j++) {
      selPop->mixtures[i][j] = pop->mixtures[i][j];  
    }                                                
    selPop->energies[i] = pop->energies[i];
  }    
  for(i = elitism;i < pop->size;i++) {
    r1 = randomUnsignedInts[i*2];
    r2 = randomUnsignedInts[i*2+1];
    r = MIN(r1,r2);
    for(j = 0;j < mixtureSize;j++) {
      selPop->mixtures[i][j] = pop->mixtures[r][j];  // this works because population is 
    }                                                // ordered according to fitness
    selPop->energies[i] = pop->energies[r];
  } 
}

void gaReorder(Population* pop)
{
  int i,j,mixtureSize;
  SortStruct* ss;
  float** mixTmp;
 
  if (pop->dim == 1) { 
    mixtureSize = 2*(pop->dim)*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  else {
    mixtureSize = (2*(pop->dim) + (pop->dim)*(pop->dim))*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
 
  ss = new SortStruct [pop->size];
  for(i = 0;i < pop->size;i++) {
    ss[i].energy = pop->energies[i];
    ss[i].index = i;
  }
  mixTmp = new float*[pop->size];
  for(i = 0;i < pop->size;i++) {
    mixTmp[i] = new float[mixtureSize];
  }
  qsort(ss, pop->size, sizeof(SortStruct), cmpStruct);
  for(i = 0;i < pop->size;i++) {
    for(j = 0;j < mixtureSize;j++) {
      mixTmp[i][j] = pop->mixtures[i][j];
    }
  }
  for(i = 0;i < pop->size;i++) {
    pop->energies[i] = ss[i].energy;
    for(j = 0;j < mixtureSize;j++) {
      pop->mixtures[i][j] = mixTmp[ss[i].index][j];
    }
  }
  delete[] ss;
  for(i = 0;i < pop->size;i++) {
    delete[] mixTmp[i];
  }
  delete[] mixTmp;

}

// Implements the blended crossover
// Assumes that population pop is sorted and 
// that the selPop is generated by the tournament selection

void gaBLX(Population* pop,Population* selPop,float xoverRate, int elitism, float alpha, float* lowLimit, float* upLimit, unsigned int* randomUnsignedInts, float* randomFloats, bool equalVar)
{
  int i,j;
  int leaveAlone;
  int mixtureSize;
  int bgVar;
  int par1,par2;
  float b,probsum;

  // first skip the indivuals based on the elistism
  // i.e elistism first (that is fittest) indivuals survive 
  // automatically and unchanged to the next generation

  // then if xoverRate < 1.0 copy necessary number of indivuals from selPop
  // to pop. this is not the typical way to implement crossover but it nevertheless 
  // works because the selPop is in random order
  if (pop->dim == 1) { 
    mixtureSize = 2*(pop->dim)*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  else {
    mixtureSize = (2*(pop->dim) + (pop->dim)*(pop->dim))*(pop->numberOfLabels - pop->numberOfPveLabels) + pop->numberOfLabels;
  }
  if (pop->dim == 1) {
    bgVar = 3;
  }
  else {
    bgVar = 2*pop->dim + 1 + (pop->dim)*(pop->dim);
  } 
  leaveAlone = (int) floor( (1 - xoverRate) * (pop->size));
  for(i = elitism;i < elitism + leaveAlone;i++) {
    for(j = 0;j < mixtureSize;j++) {
      pop->mixtures[i][j] = selPop->mixtures[i][j];
    }
  }
  for(i = elitism + leaveAlone;i < pop->size;i++) {
    par1 = randomUnsignedInts[i*2];
    par2 = randomUnsignedInts[i*2+1];
    for(j = bgVar;j < mixtureSize;j++) {
      b = randomFloats[i*mixtureSize+j];
      b = (1 + 2*alpha)*b - alpha;
      if((b < (-0.5)) || (b > 1.5)) cout << "b" << b << endl; 
      pop->mixtures[i][j] = b*(selPop->mixtures[par1][j]) + (1 - b)*(selPop->mixtures[par2][j]);
      pop->mixtures[i][j] = MIN(pop->mixtures[i][j],upLimit[j]);
      pop->mixtures[i][j] = MAX(pop->mixtures[i][j],lowLimit[j]);
    }
    probsum = 0.0;
    for(j = 1;j < pop->numberOfLabels;j++) {
      probsum = probsum + gaGetProb(pop,i,j);
    }
    for(j = 1;j < pop->numberOfLabels;j++) {
      gaSetProb(pop,i,j,gaGetProb(pop,i,j) / probsum);
    }
    // take care of the background variance
    pop->mixtures[i][1] = pop->mixtures[i][4];
    if(equalVar) {
      for(j = 2;j < (pop->numberOfLabels - pop->numberOfPveLabels);j++) {
        gaSetSigma2(pop,i,j,gaGetSigma2(pop,i,1));
      }
    }
  }
  
  
}


// Implements the sort operator. Uses the trivial Bubble Sort to achieve its goal.
// The Bubble Sort should be handy because the number of floats that need to be sorted is small.

void gaSortPopulation(Population* pop,int sortDim)
{
  int i;
  int componentSize,pureLabels;
  
  pureLabels = pop->numberOfLabels - pop->numberOfPveLabels;
   if (pop->dim == 1) { 
    componentSize = 3;
  }
  else {
    componentSize = 2*(pop->dim) + (pop->dim)*(pop->dim) + 1;
  }  
  for(i = 0;i < pop->size;i++) {
    bubbleSort(pop->mixtures[i],pureLabels,componentSize);
  }
}

// whether to terminate
// Assumes sorted population
// return 'true' if termination condition is fulfilled

bool gaTerminate(Population* pop, float thr) 
{
  float mean,best;
  int i;
  
  mean = 0.0;
  for(i = 0;i < pop->size;i++) {
    mean = mean + pop->energies[i];
  }
  mean = mean / pop->size;
  best = pop->energies[0];
  return( (mean - best) < thr);

}




