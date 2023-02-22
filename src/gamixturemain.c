// ******************************************************
// Mixture model optimization using genetic algorithms
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

// Version 1.0: The first version of the software.
// Version 1.1: Added the possibility of restarts, and speeded up the Parzen window computation. JT 17th Jul 2007.
// Version 2.0: The first public release of the SVPASEG. No changes to gamixture 

#include "gamixture.h"


int main(int argc,char** argv)
{
  AnalyzeImage img;
  AnalyzeImage mask;
  AnalyzeLabelImage mask2;
  AnalyzeImage** atlasImages;
  MixtureSpec mixture;
  AtlasSpec atlas; 
  PdfEstimate hatf;
  Parameters params;  
  Population popRuns;

  float maxMu,minMu,minVar,maxVar,mean;
  float* lowLimit;
  float* upLimit;  
  //  float* fitness;

  bool boolstatus;
  bool changed_n = false;
  bool signedData = true;

  int intstatus,rr; 
  int pveLabels,pureLabels;  

  if(argc > 2) {
    if(!(strcmp(argv[1],"-unsigned"))) {
      signedData = false;
      argc--;
      for(int i = 1;i < argc;i++) {
        argv[i] = argv[i + 1];
      }
    }
  }

  if(argc < 5) {
    cout << "Genetic algorithm for mixture model optimization" << endl;
    cout << "Usage: gamixture [-unsigned] imagefile maskfile atlasfile fmmparamfile [parameters]" << endl;
    cout << "Optional parameters are :" << endl;
    cout << "-alpha:             parameter for blended crossover, defaults to " << DEFAULT_ALPHA << endl;
    cout << "-size:              population size, defaults to " << DEFAULT_POPSIZE << endl;
    cout << "-terminationthr:    threshold for terminating the algorithm, defaults to "<< DEFAULT_TERMINATIONTHR << endl; 
    cout << "-xoverrate:         crossover rate, defaults to " << DEFAULT_XOVERRATE << endl;
    cout << "-maxgenerations:    maximum number of generations, defaults to " << DEFAULT_MAXGENERATIONS << endl;
    cout << "-sortpop            whether to use the permutation operator, defaults to " << DEFAULT_SORTPOP << endl;
    cout << "-parzenn            number of points for Parzen estimate (for approximate ML), defaults to " << DEFAULT_PARZENN << endl;
    cout << "-parzensigma        window width parameter for the Parzen estimate, defaults to " << DEFAULT_PARZENSIGMA << endl;         
    cout << "-equalvar           whether component densities should have equal variances, defaults to " << DEFAULT_EQUALVAR << endl;  
    cout << "-restarts           the number of individual runs of the GA, defaults to " << DEFAULT_RESTARTS << endl;  
   return(1);
  }
  intstatus = readImage(argv[1],&img,signedData);
  if(intstatus != 0) {
    cout << "Could not read image file " << argv[1] << " " << intstatus << endl;
    return(2);
  }  
  cout << "Image dimensions: " << img.header.x_dim << " " << img.header.y_dim << " " << img.header.z_dim << endl;
  if(!(strcmp(argv[2],"default"))) {
    
    copyImage(&img,&mask);
    thresholdImage(&mask,VERY_SMALL);
    boolstatus = binaryMask(&img,&mask2,VERY_SMALL);
    
  }
  else {
    intstatus = readImage(argv[2],&mask);
    if(intstatus != 0) {
      cout << "Could not read image (brainmask) file " << argv[2] << " " << intstatus << endl;
      return(3);
    }  
    thresholdImage(&mask,0.0001);
    binaryMask(&mask,&mask2,0.5);
  } 
   
  boolstatus = readAtlasSpec(&atlas,argv[3]);
  if(boolstatus == false) {
    cout << "Could not read atlas file " << argv[3] << endl;
    return(4);
  }   
 
  intstatus = parseParams(&params,argc - 4,&(argv[4]));
  if(intstatus != 0) {  
    cout << "Incorrect parameter value input" << endl;
    return(6);
  }
  if(atlas.n > 0) {
    atlasImages = new AnalyzeImage*[atlas.n];
    for(int i = 0;i < atlas.n;i++) {
      atlasImages[i] = new AnalyzeImage;
    }
    intstatus = readAtlasImages(&atlas,atlasImages);
    if(intstatus != 0) {
      cout << "Could not read probabilistic atlas. Error: " << intstatus << endl;
    return(6);
    }  
    cout << "The atlas files have been read" << endl;
    if(maskAtlas(&atlas,atlasImages,&mask) == false) {
      cout << "Could not mask atlas" << endl;
      return(7);
    }  
    cout << "The atlas has been masked" << endl;
  }
  else {
   // taking care for the case where atlas->n == 0
   // then we assume that the atlas is defined by the mask
    atlas.n = 1;
    changed_n = true;
    atlasImages = new AnalyzeImage*[atlas.n];
    atlasImages[0] = new AnalyzeImage;
    boolstatus = copyImage(&mask,atlasImages[0]);
  } 
  allocateMixtureSpec(&atlas,&mixture);
  
  // compute the number of pve labels and pure labels.
  // remember that pve labels have to have indeces greater than pure labels
  pureLabels = 0;
  pveLabels = 0;
  for(int i = 0;i < atlas.numberOfLabels;i++) {
    if(atlas.labelTypes[i].pureLabel == true) pureLabels++;
    else pveLabels++;     
  }
  cout << "purelabels:  " << pureLabels << " Pvelabels: " << pveLabels << endl; 
  
  
  // initializing the Parzen windows
  hatf.n = params.parzenN;
  computeX(&hatf,&img,&mask2,0.999);
  setSigma(&hatf,params.parzenSigma); 
  
  // Setting the limits for GA
  maxMu = hatf.x[hatf.n - 1];
  minMu = hatf.x[0];
  minVar = hatf.sigma * hatf.sigma;
  maxVar = 0.0;  // this is re-set later on  

  lowLimit = new float[3*pureLabels + pveLabels]; 
  upLimit = new float[3*pureLabels + pveLabels]; 
  // initialize by defults, then adapt the necessary values
  for(int i = 0; i < pureLabels;i++) {
    upLimit[3*i] = maxMu;
    lowLimit[3*i] = minMu;
    upLimit[3*i + 1] = maxVar;  // this is re-set later on
    lowLimit[3*i + 1] = minVar;
    if(changed_n)  {
      upLimit[3*i + 2] = 1.0;  
    }
    else { 
      upLimit[3*i + 2] = 0.0; 
    } // this is changed afterwards
    lowLimit[3*i + 2] = 0.0;
  }
  upLimit[0] = 0.0;  // adjustment for the background label
  lowLimit[0] = 0.0;
  upLimit[2] = 0.0;
   
  // pve classes
  for(int i = 0; i <pveLabels;i++) {
    if(changed_n) {
      upLimit[pureLabels*3 + i] = 1.0;
    }
    else {
      upLimit[pureLabels*3 + i] = 0.0;
    }
    lowLimit[pureLabels*3 + i] = 0.0;
  }  
  for(int i = 0; i < atlas.n;i++) {
    rr = i + 1;
    cout << "Region " << rr << ":computing parzen estimate..." << endl;;  
    computeY(&hatf,&img,atlasImages[i]);
    //  if(i == 5) writeEstimate("tmp.parzen",&hatf,true);
  
    cout << "Region " << rr << ":Genetic algorithm..." << endl;
    // set the region wise limits
    maxVar = 0.0;
    mean = 0.0;
    for(int j = 0;j < hatf.n; j++) {
      mean = mean +  (hatf.x[1] - hatf.x[0])*(hatf.y[j])*(hatf.x[j]);
    }
    for(int j = 0;j < hatf.n; j++) {
      maxVar = maxVar + (hatf.x[1] - hatf.x[0])*hatf.y[j]*(hatf.x[j] - mean)*(hatf.x[j] - mean);
    }
    maxVar = maxVar/pureLabels;
    for(int j = 0; j < pureLabels;j++) {
       upLimit[3*j + 1] = maxVar;
    }
    if(!atlas.regionLowProb.empty()) {
      for(int j = 1; j < pureLabels;j++) {
        lowLimit[j*3 + 2] = atlas.regionLowProb[i][j];
        upLimit[j*3 + 2] = atlas.regionUpProb[i][j];
      }
            for(int j = 0; j < pveLabels;j++) {
        lowLimit[pureLabels*3 + j] = atlas.regionLowProb[i][pureLabels + j];
        upLimit[pureLabels*3 + j] = atlas.regionUpProb[i][pureLabels + j];
      }
    }
    else {
      if(!changed_n) {
        for(int j = 1; j < pureLabels;j++) { // skip background
          upLimit[3*j + 2] = 0.0; 
        }
        for(int j = 0; j < pveLabels;j++) {
          upLimit[pureLabels*3 + j] = 0.0;
        }
        for(int j = 0;j < atlas.permittedLabels[i].len;j++) {
	  //   cout << atlas.permittedLabels[i].list[j] << " ";  
          if(atlas.permittedLabels[i].list[j] > (pureLabels -1)) {
            upLimit[3*pureLabels + atlas.permittedLabels[i].list[j] - pureLabels] = 1.0;
          }
          else {
            upLimit[atlas.permittedLabels[i].list[j]*3 + 2] = 1.0;
          }
        }
      }
    } 

    RFRandom kInitRNG;
    kInitRNG.Seed(i+1);
    gaInitializePopulation(&popRuns,params.restarts,1,pureLabels + pveLabels, pveLabels,
		atlas.labelTypes.data(),lowLimit,upLimit, kInitRNG, params.equalVar);

    #pragma omp parallel for
    for(int n = 0;n < params.restarts;n++) {
      Population pop, selPop;
      RFRandom kRNG;
      kRNG.Seed((n+1)*117*(i+1));
      gaInitializePopulation(&pop,params.size,1,pureLabels + pveLabels,pveLabels,
			atlas.labelTypes.data(),lowLimit,upLimit, kRNG, params.equalVar);
      gaSortPopulation(&pop,1);
      gaEvaluate(&pop,&hatf);
      gaReorder(&pop);
      copyPartialPopulation(&pop,&selPop);
      bool terminate = false;
      int itercount = 0;
      while((!terminate) && (itercount < params.maxGenerations)) {
			gaTournamentSelection(&selPop,&pop,1,kRNG);
			gaBLX(&pop,&selPop,params.xoverRate,1,params.alpha,lowLimit,upLimit,kRNG,params.equalVar); 
        if(params.sortPop) {
          gaSortPopulation(&pop,1);
        }
        gaEvaluate(&pop,&hatf);
        gaReorder(&pop);
        terminate = gaTerminate(&pop,params.terminationThr);
        itercount++;
      }
      if(!(params.sortPop)) gaSortPopulation(&pop,1); 
      cout << "GA converged after " << itercount << " iterations." << endl;
      cout << "the KL distance is " << pop.energies[0] << "." << endl;
      // put the best invividual into popRuns
     for(int j = 0;j < pureLabels;j++) {
        gaSetMu(&popRuns,n,j,gaGetMu(&pop,0,j));
        gaSetSigma2(&popRuns,n,j,gaGetSigma2(&pop,0,j));
        gaSetProb(&popRuns,n,j,gaGetProb(&pop,0,j));     
      }
      for(int j = 0;j < pveLabels;j++) {
        gaSetProb(&popRuns,n,j + pureLabels,gaGetProb(&pop,0,j + pureLabels));
      }
      popRuns.energies[n] = pop.energies[0];
    }
    if( params.restarts > 1) gaReorder(&popRuns); 
    // convert the best individual to mixtureSpec
    for(int j = 0;j < pureLabels;j++) {
      putMu(&mixture,i,j,gaGetMu(&popRuns,0,j));
      putSigma2(&mixture,i,j,gaGetSigma2(&popRuns,0,j));
      putProb(&mixture,i,j,gaGetProb(&popRuns,0,j));     
    }
    for(int j = 0;j < pveLabels;j++) {
      putProb(&mixture,i,j + pureLabels,gaGetProb(&popRuns,0,j + pureLabels));
    }
  }
  printMixture(&mixture);
  intstatus = writeMixtureParameters(argv[4],&atlas,&mixture, true);
  if(intstatus != 0) {
    cout << "Error in writing the mixture parameters file" << argv[4] << endl;
    return(7);
  }
  freeAtlasImages(&atlas,atlasImages);
  return(0); 
}
