// ******************************************************
// SVPASEG main function
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


#include "analyze.h"
#include "atlasspec.h"
#include "nhmrf.h"

#include <vector>

int main(int argc, char** argv) 
{
  AnalyzeImage img;
  AnalyzeImage mask;
  AnalyzeLabelImage labelImg;
  AnalyzeLabelImage pveLabelImg;
  std::vector<AnalyzeImage> atlasImages;
  std::vector<AnalyzeImage> labelLikelihoods;
  AnalyzeImage** TPMImages;
  MixtureSpec mixture;
  AtlasSpec atlas;
  SvpasegParameters params;
  
  int intstatus,i,itercount,pureLabels; 
  
  bool boolstatus,clobber,markov;
  bool signedData = true;  
  bool useTPM = false;
  bool onlyPureLabels = false;
 
  markov = false; // XXXXXXX was undefined

  float beta1,beta2;

  // the call should be :    svpaseg imagefilename brainmask_filename
  //                             atlas_def_filename mixture_params_filename 
  //                             labelfilename [optional parameters]
  
  if(argc > 2) {
    if(!(strcmp(argv[1],"-unsigned"))) {
      signedData = false;
      argc--;
      for(i = 1;i < argc;i++) {
        argv[i] = argv[i + 1];
      }
    }
  }
 
  if(argc < 6) {
     cout << "Tissue classification based on local Markov Random Fields" << endl;
    cout << "Usage: svpaseg [-unsigned] imagefile maskfile atlasfile fmmparamfile labelfile [parameters]" << endl;
     cout <<  "If you say 'default' for maskfile, it is assumed that all non zero intensities in" << endl;
      cout <<  "input image are within the brain mask" << endl;  
    cout << "Optional parameters are :" << endl;
    cout << "-beta1:             regularization parameter first first-order interactions," << endl;
    cout << "                    defaults depend on particular model setting." << endl;
    cout << "                    It is NOT recommended to change the defaults. " << endl;
    cout << "-beta2:             regularization parameter for second order interactions" << endl; 
    cout << "                    i.e. controls the strenght of the spatial regularization component. " << endl; 
    cout << "                    defaults to "<< DEFAULT_BETA2 << endl;
    cout << "-pvelabels:         controls whether intermediate PVE classification labels are written," << endl;
    cout << "                     not relevant except for p-type atlases.  " << endl; 
   
  
    return(1); 
  }

  clobber = true;  // remember to check writeLabelImage
  
  intstatus = readImage(argv[1],&img,signedData);
  if(intstatus != 0) {
    cout << "Could not read image file " << argv[1] << " " << intstatus << endl;
    return(2);
  }  
  if(!(strcmp(argv[2],"default"))) {
    copyImage(&img,&mask);
    thresholdImage(&mask,VERY_SMALL);
  }
  else {
    intstatus = readImage(argv[2],&mask);
    if(intstatus != 0) {
      cout << "Could not read image (brainmask) file " << argv[2] << " " << intstatus << endl;
      return(3);
    }  
    thresholdImage(&mask,0.0001);
  }
  boolstatus = readAtlasSpec(&atlas,argv[3]);
  if(boolstatus == false) {
    cout << "Could not read atlas file " << argv[3] << endl;
    return(4);
  }  

  intstatus = readMixtureParameters(argv[4],&atlas,&mixture);
  if(intstatus != 0) {
    cout << "Could not read mixture file " << argv[4] << " " << intstatus << endl;
    return(5);
  }  
  intstatus = parseParamsSvpaseg(&params,argc - 5,&(argv[5]),&atlas);
  if(intstatus != 0) {  
    cout << "Incorrect parameter value input" << endl;
    return(6);
  }

 

  if(atlas.n > 0) {
	atlasImages.resize(atlas.n);
    intstatus = readAtlasImages(&atlas,atlasImages);
    if(intstatus != 0) {
      cout << "Could not read probabilistic atlas. Error: " << intstatus << endl;
    return(7);
    }  
    cout << "The atlas files have been read" << endl;
    if(maskAtlas(&atlas,atlasImages,&mask) == false) {
      cout << "Could not mask atlas" << endl;
      return(8);
    }  
    cout << "The atlas has been masked" << endl;
  }
  else {
   // taking care for the case where atlas->n == 0
   // then we assume that the atlas is defined by the mask
    atlas.n = 1;
    mixture.patlas->n = 1;
   
    atlasImages.resize(1);
    boolstatus = copyImage(&mask,&atlasImages[0]);
  }  
  

  if(!(newLabelImage(&pveLabelImg,&img))) {
    cout << "Failed to create pveLabel image" << endl;
    return(9);
  }  

  if(!(newLabelImage(&labelImg,&img))) {
    cout << "Failed to create pveLabel image" << endl;
    return(9);
  }  
  
  if(atlas.onlyPureLabels) {
    pureLabels = countPureLabels(&atlas);
    if(pureLabels == 0) {
      cout << "Wrong label order: Pure labels should be first" << endl;
      return(10);
    }
  }
/*  if(atlas.useTPM) {
    TPMImages = new AnalyzeImage*[pureLabels];
    for(i = 0;i < pureLabels;i++) {
      TPMImages[i] = new AnalyzeImage;
    }
    intstatus = readTPMimages(&atlas,TPMImages,&mask,pureLabels);
    if(intstatus != 0) {
      cout << "Could not read TPMs. Error: " << intstatus << endl;
      return(11);
    }  
    //    intstatus = processTPMimages(&atlas, TPMimages,&mask,pureLabels);
    // if(intstatus != 0) {
    //  cout << "Could not process TPMs. Error: " << intstatus << endl;
    //  return(11);   
    //  }
    // cout << "TPMs have been read and processed" << endl;
  }*/

  cout << "Computing the ML classification" << endl; 
  //labelLikelihoods = new AnalyzeImage*[atlas.numberOfLabels - 1];
  labelLikelihoods.resize(atlas.numberOfLabels - 1);
  for(i = 0;i < (atlas.numberOfLabels - 1);i++) {
    if(!(newImage(&labelLikelihoods[i],&img))) {
      cout << "Failed to create likelihood image" << endl;
      return(12);
    }  
  }  
  // maximum likelihood labeling.
  intstatus = computeVoxelLikelihood(&mixture,&img,&mask,atlasImages,labelLikelihoods);
  if(intstatus != 0) {
    cout << "something wrong with the ML classification " << intstatus << endl;
    return(13);
  }
  printLabelInfo(&atlas);
  cout << "Iterated Conditional modes..." << endl;
/*  if(markov) {
    itercount = computeMRF(&pveLabelImg,&mixture,&mask,labelLikelihoods,atlasImages, params.beta1, params.beta2,50,true);
  }
  else {
    if(atlas.useTPM) {
      itercount = computeGibbsAtlas(&pveLabelImg,&mixture,&mask,labelLikelihoods,atlasImages, TPMImages, params.beta1, params.beta2,50,true);
    }
    else { 
      if(atlas.onlyPureLabels) {
        itercount = computeGibbsPure(&pveLabelImg,&mixture,&mask,labelLikelihoods,atlasImages, params.beta1, params.beta2,50,true);
     } 
      else {*/
        itercount = computeGibbs(&pveLabelImg,&mixture,&mask,labelLikelihoods,atlasImages, params.beta1, params.beta2,50,true);
     /* }
    }
  }*/
  cout << "Iterations: " << itercount << endl;
  intstatus = convertPVElabels(&labelImg,&pveLabelImg,&img,atlasImages,&mixture);
  if(intstatus != 0) {
    cout << "Conversion to pure labels did not succeed" << endl;
    return(14);
  }  
  cout << argv[5] << endl;
  intstatus = writeLabelImage(argv[5],&labelImg,clobber);
  if(intstatus != 0) {
    cout << "Could not write labeled image. Error: " << intstatus << endl;
    return(15);
  }  

  //TODO: this return statement was put in because of a segmentation fault. dig into this
  exit(0);

  if(params.writePveLabelImage) {
    intstatus = writeLabelImage(params.pveLabelImage,&pveLabelImg,clobber);
    if(intstatus != 0) {
      cout << "Could not write pve labeled image. Error: " << intstatus << endl;
      return(16);
    }
  }

  freeImage(&img);
  freeImage(&mask);
  freeLabelImage(&pveLabelImg);
  freeLabelImage(&labelImg);  
 
 
  return(0);
  // think about freeing the rest of the images as well

}


