// ******************************************************
// Functions for SVPA based non homogeneous MRFs 
//
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

#include "nhmrf.h"

// reads the parameters for the mixture density fron a file.
// the exact format of the file is decided by the given MixtureSpec,
// more precisely AtlasSpec attached to it.
// the format is 
// region 1 params 
// region 2 params
//  ....
// region n params

// the classes that are not allowed for some region should have 
// the probability parameter value 0.0 for these regions.

int readMixtureParameters(char* filename,AtlasSpec* atlas,MixtureSpec* mixture)
{
  int r,l;
  float val;
  bool changed_n = false;

  ifstream ifile;
  ifile.open(filename);
  if(!ifile.is_open()) return(1); 
  mixture->patlas = atlas;
  if(mixture->patlas->n == 0) {
    mixture->patlas->n = 1;
    changed_n = true;
  }

  // reserving the memory
  allocateMixtureSpec(mixture->patlas,mixture);
  
  for(r = 0;r < mixture->patlas->n;r++) {
    for(l = 1;l < mixture->patlas->numberOfLabels;l++) {
      if(mixture->patlas->labelTypes[l].pureLabel) {
	ifile >> val;
        putMu(mixture,r,l,val);       
        ifile >> val;
        putSigma2(mixture,r,l,val);
        ifile >> val;
        putProb(mixture,r,l,val);
      }
      else {
        ifile >> val;
        putProb(mixture,r,l,val);
      }
    }
  }
  // place values for the background class
  for(r = 0;r < mixture->patlas->n;r++) {
    putMu(mixture,r,0,0.0); 
    putSigma2(mixture,r,0,getSigma2(mixture,r,1));
    putProb(mixture,r,0,0.0);
  }
  
  if(changed_n) mixture->patlas->n = 0;
  ifile.close();
  return(0);
}

int writeMixtureParameters(char* filename,AtlasSpec* atlas,MixtureSpec* mixture,bool overwrite)
{
  int r,l;
  float val;
  bool changed_n = false;

  if(! overwrite) {
    ifstream ifile;
    ifile.open(filename);
    if(ifile.is_open()) {
      ifile.close();
      return(1); 
    }
  }
  ofstream ofile;
  ofile.open(filename);
  if(!ofile.is_open()) return(2);
   if(mixture->patlas->n == 0) {
    mixture->patlas->n = 1;
    changed_n = true;
  }
  for(r = 0;r < mixture->patlas->n;r++) {
    for(l = 1;l < mixture->patlas->numberOfLabels;l++) {
      if(mixture->patlas->labelTypes[l].pureLabel) {
	ofile << getMu(mixture,r,l);
        ofile << ' ';       
        ofile << getSigma2(mixture,r,l);
        ofile << ' ';
        ofile << getProb(mixture,r,l);
        ofile << ' ';
      }
      else {
        ofile << getProb(mixture,r,l);
        ofile << ' ';
      }
    }
    ofile << endl;
  }  
  if(changed_n) mixture->patlas->n = 0;
  ofile.close();
  return(0);
}

int parseParamsSvpaseg(SvpasegParameters* param,int n,char** arguments,AtlasSpec* atlas) 
{
  // Initialize with the default values
  int i;  

  param->beta1 = DEFAULT_BETA1;
  param->beta2  = DEFAULT_BETA2;
  param->writePveLabelImage = false;
  param->pveLabelImage = NULL;
  param->markov = false;  // this is the only setting supported at this moment

  // check that beta1 matches with the choice of the method
  if(atlas->useTPM) param->beta1 = 1.0;
  if(!(atlas->onlyPureLabels)) param->beta1 = 1.0;
  if((!(atlas->useTPM)) && (atlas->onlyPureLabels)) param->beta1 = 0.0; 
    

  if( ((n - 1) % 2) > 0 ) return(1);
  for(i = 1;i < (n - 1)/2 + 1;i++) {
    if(!strcmp(arguments[2*i - 1],"-beta1"))
      param->beta1 = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-beta2"))
      param->beta2 = atof(arguments[2*i]);
    if(!strcmp(arguments[2*i - 1],"-pvelabels")) {
      param->pveLabelImage = arguments[2*i];
      param->writePveLabelImage = true;
    }
  }
  return(0);
}


// Computes the value of the likelihood term of an intensity value in an image for each label
// AnalyzeImage(s) pointed by labelLikelihoods should be allocated.  
// Atlas images should be in the same space than the image
// The mask is assumed to be thresholded
// A value 0 is returned if everything is ok.


int computeVoxelLikelihood(MixtureSpec* mixture,AnalyzeImage* img,AnalyzeImage* mask,std::vector<AnalyzeImage> & atlasImages,std::vector<AnalyzeImage> & labelLikelihoods)
{
  int dimx,dimy,dimz;
  // First test that the mask size and image size and atlas size match

  //  printMixture(mixture);

  dimx = img->header.x_dim;
  dimy = img->header.y_dim; 
  dimz = img->header.z_dim;

  if(( dimx != mask->header.x_dim ) || ( dimy != mask->header.y_dim ) 
                                    || ( dimz != mask->header.z_dim ))
    return(2);

  if(( dimx != atlasImages[0].header.x_dim ) || ( dimy != atlasImages[0].header.y_dim ) 
     || ( dimz != atlasImages[0].header.z_dim )) {
    if(( dimx > atlasImages[0].header.x_dim ) || ( dimy > atlasImages[0].header.y_dim ) 
       || ( dimz > atlasImages[0].header.z_dim )) {
      return(1);
    }
    else {
      cout << "Warning: image dimensions do not match to the atlas dimensions->continuing but check the results" << endl;
    }
  }
  
  // remember: the label 0 is reserved for background 
  if((mixture->patlas->n) == 0) {
    int lvalue;
    int label1,label2;
    int j = 0;
    int i = 0;
    if(mask->data[i] > 0.5) {
      for(int k = 1;k < (mixture->patlas->numberOfLabels);k++) {
        labelLikelihoods[k - 1].data[i] = 0.0;
      }
      for(int k = 1;k < (mixture->patlas->numberOfLabels);k++) {
        if(fabs(getProb(mixture,j,k)) > 0.0001) {
          if(mixture->patlas->labelTypes[k].pureLabel) {
            lvalue = computeGaussianLikelihood(img->data[i],getMu(mixture,j,k),
                                                   getSigma2(mixture,j,k));
          }
          else {
            label1 = mixture->patlas->labelTypes[k].mixed[0];
            label2 = mixture->patlas->labelTypes[k].mixed[1];
            lvalue = computeMarginalizedLikelihood(img->data[i],getMu(mixture,j,label1),
                                                       getMu(mixture,j,label2),
                                                       getSigma2(mixture,j,label1),
						       getSigma2(mixture,j,label2),10,0.0,1.0);
          }
        }
        else {
          lvalue = 0.0;
        }
        labelLikelihoods[k - 1].data[i] = labelLikelihoods[k - 1].data[i] 
                                         + lvalue;      
      }
    }
  }
  else {
    #pragma omp parallel for
    for(int z = 0;z < dimz;z++) {
      int label1,label2;
      float lvalue,tmpval;

      for(int y = 0;y < dimy;y++) {
        for(int x = 0;x < dimx;x++) {
          if(getVoxelValue(mask,x,y,z) > 0.5) {
            for(int k = 1;k < (mixture->patlas->numberOfLabels);k++) {
              putVoxelValue(&labelLikelihoods[k - 1],x,y,z,0.0);
            }
            for(int j = 0; j < mixture->patlas->n; j++) {
              if(getVoxelValue(&atlasImages[j],x,y,z) > 0.000001) {
                for(int k = 1;k < (mixture->patlas->numberOfLabels);k++) {
                  if(fabs(getProb(mixture,j,k)) > 0.0001) {
                    if(mixture->patlas->labelTypes[k].pureLabel) {
                      lvalue = computeGaussianLikelihood(getVoxelValue(img,x,y,z),getMu(mixture,j,k),
                                                   getSigma2(mixture,j,k));
                    }
                    else {
                      label1 = mixture->patlas->labelTypes[k].mixed[0];
                      label2 = mixture->patlas->labelTypes[k].mixed[1];
                      lvalue = computeMarginalizedLikelihood(getVoxelValue(img,x,y,z),getMu(mixture,j,label1),
                                                       getMu(mixture,j,label2),
                                                       getSigma2(mixture,j,label1),
							   getSigma2(mixture,j,label2),10,DEFAULT_PVESTART,DEFAULT_PVEEND);
                    }
                  }
		  else {
                    lvalue = 0.0;
                  }
		  // labelLikelihoods[k - 1]->data[i] = labelLikelihoods[k - 1]->data[i] 
		  //                     + (atlasImages[j]->data[i]) * lvalue; 
		  tmpval = getVoxelValue(&labelLikelihoods[k - 1],x,y,z) + getVoxelValue(&atlasImages[j],x,y,z) * lvalue; 
                  putVoxelValue(&labelLikelihoods[k - 1],x,y,z,tmpval);
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(0);
}
 
 
// minimizes the posterior by the ICM algorithm. returns the number of iterations required
// labels that are forbidden for each region should receive the probability 0 in the mixture. 
// no explicit checking for forbidden labels is done

/*
int computeMRF(AnalyzeLabelImage* labels,MixtureSpec* mixture,AnalyzeImage* mask,AnalyzeImage** labelLikelihoods, AnalyzeImage** atlasImages, float beta1, float beta2, int maxIterations,bool verbose)
{
  int x,y,z;
  int i,j,k,r;
  char l;
  int dimx,dimy,dimz;
  int numberOfLabels,numberOfRegions;
  float* voxelProb;
  float* gibbsProb;
  float* posteriorProb;
  bool changed = true;
  int iteration = 0;
  float distLookup[27];
  float sliceWidth[3], sliceWidthMin;
  char newLabel;
  AnalyzeImage tmpimg;

  //  bool testverbose;
  // int count = 0;

  if(verbose) cout << "Initializing... " << endl;
  // save few commonly needed values  
  dimx = mask->header.x_dim;
  dimy = mask->header.y_dim;  
  dimz = mask->header.z_dim;
  numberOfLabels = mixture->patlas->numberOfLabels;
  numberOfRegions =  mixture->patlas->n;
  sliceWidthMin = MIN(fabs(labels->header.x_size),fabs(labels->header.y_size));
  sliceWidthMin = MIN(sliceWidthMin,fabs(labels->header.z_size));

  sliceWidth[0] = (labels->header.x_size) / sliceWidthMin;
  sliceWidth[1] = (labels->header.y_size) / sliceWidthMin;
  sliceWidth[2] = (labels->header.z_size) / sliceWidthMin;
  cout << sliceWidth[0] << " " << sliceWidth[1] << " " << sliceWidth[2] << " " << endl;
  // allocate
  voxelProb = new float[numberOfLabels - 1];
  gibbsProb = new float[numberOfLabels - 1];  
  posteriorProb = new float[numberOfLabels - 1];

  

  // fill the look up table
  for(i = (-1);i < 2;i++) {
    for(j = (-1);j < 2;j++) {
      for(k = (-1);k < 2;k++) {
        distLookup[ (i + 1) * 9 + (j + 1) * 3 + k + 1 ] = sqrt(pow(sliceWidth[0] * abs(i),2) +
                                                    pow(sliceWidth[1] * abs(j),2) +
                                                    pow(sliceWidth[2] * abs(k),2));
      }
    }
  }
  // start by initializing the ICM
  for(x = 0; x < dimx ; x++) {
    for(y = 0; y < dimy; y++) {
      for(z = 0;z < dimz; z++) {
	//         if((x == 160) && (y == 160) && (z == 190))  testverbose = true;
        // else testverbose = false;  
        if(getVoxelValue(mask,x,y,z) > 0.5) {
          collectValuesFromImagePP(labelLikelihoods,voxelProb,x,y,z,numberOfLabels - 1);
          for(l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = 0.0;
          }
          for(r = 0;r < numberOfRegions;r++) {
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = posteriorProb[l]  +  getVoxelValue(atlasImages[r],x,y,z) * (exp ( beta1 * log (getProb(mixture,r,(l + 1)) + 0.0001)));
            }
          }
          for(l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = posteriorProb[l] * voxelProb[l];
          }
          putLabelValue(labels,x,y,z,maxArg(posteriorProb,numberOfLabels - 1) + 1);
          
	  //  count++;
	  //  if( testverbose ) cout << "test " <<  (int) (maxArg(voxelProb,numberOfLabels - 1) + 1) << endl;
        }
        else {
          putLabelValue(labels,x,y,z,0);
          // if( testverbose ) cout << "test2 " << endl; 
        }
      }
    }
  }
  

  // start the ICM
  while(changed && (iteration < maxIterations)) {
    changed = false;
    iteration++;
    if( verbose ) cout << "iteration " << iteration << endl;
    for(x = 0; x < dimx ; x++) {
      for(y = 0; y < dimy; y++) {
        for(z = 0;z < dimz; z++) {
          if(getVoxelValue(mask,x,y,z) > 0.5) {
            
            //  compute the second order term in the prior, non-normalized
            // this is the same for every region.
            for(l = 0;l < (numberOfLabels - 1);l++) {
              gibbsProb[l] = secondOrderGibbs(l + 1,labels,mixture->patlas->mrfConstants[l + 1],
                                              x,y,z,distLookup,beta2);
	      //  if(testverbose) cout << " " << gibbsProb[l];
            }
	    //            if(testverbose) cout << endl;
            // then start to calculate the posterior probabilities region by region
            // initialize label probabilities to 0
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = 0.0;
            }
            // then compute region wise prior
            for(r = 0;r < numberOfRegions;r++) {
              for(l = 0;l < (numberOfLabels - 1);l++) {
                voxelProb[l] = ( gibbsProb[l] ) * (exp ( beta1 * log (getProb(mixture,r,(l + 1)) + 0.0001)));
              }
              normalize(voxelProb,numberOfLabels - 1);
              for(l = 0;l < (numberOfLabels - 1);l++) {
                posteriorProb[l] = posteriorProb[l] 
                                 + getVoxelValue(atlasImages[r],x,y,z) * voxelProb[l];
              }    
            }
	    //  if(testverbose) {
            //  for(l = 0;l < (numberOfLabels - 1);l++) {
	    //   cout << " " <<  posteriorProb[l];
            //  }
            //  cout << endl;
            // } 
            // now the prior probability is computed and it remains to 
            // multiply it with the likelihood term to get the posterior
            for(l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] =  posteriorProb[l] * getVoxelValue(labelLikelihoods[l],x,y,z);
            }
            // then just find the minimum posterior and update the label
            newLabel = maxArg(posteriorProb,(numberOfLabels - 1)) + 1;
            if(newLabel != getLabelValue(labels,x,y,z)) {
              changed = true;
              putLabelValue(labels,x,y,z,newLabel);  
            }
            // if(em) {
	    // cout << "Re-estimating parameters" << endl;
            //  if(mixture->patlas->labelTypes[newLabel].pureLabel) {
                

          } // endif
        }   // end for z
      }     // end for y
    }       // end for x
  }         // endwhile

  // free memory
  delete[] voxelProb; 
  delete[] gibbsProb; 
  delete[] posteriorProb;
 
  return(iteration);
}
*/
int computeGibbs(AnalyzeLabelImage* labels,MixtureSpec* mixture,AnalyzeImage* mask,std::vector<AnalyzeImage> & labelLikelihoods, std::vector<AnalyzeImage> & atlasImages, float beta1, float beta2, int maxIterations,bool verbose)
{
  int dimx,dimy,dimz;

  bool changed = true;
  int iteration = 0;
  float distLookup[27];
  float sliceWidth[3], sliceWidthMin;

  bool testverbose;

  // int count = 0;
  if(verbose) cout << "Initializing... " << endl;
  
  // save few commonly needed values  
  dimx = mask->header.x_dim;
  dimy = mask->header.y_dim;  
  dimz = mask->header.z_dim;
  const int numberOfLabels = mixture->patlas->numberOfLabels;
  const int numberOfRegions =  mixture->patlas->n;
  sliceWidthMin = MIN(fabs(labels->header.x_size),fabs(labels->header.y_size));
  sliceWidthMin = MIN(sliceWidthMin,fabs(labels->header.z_size));

  sliceWidth[0] = (labels->header.x_size) / sliceWidthMin;
  sliceWidth[1] = (labels->header.y_size) / sliceWidthMin;
  sliceWidth[2] = (labels->header.z_size) / sliceWidthMin;
  cout << "voxel dimensions:" << sliceWidth[0] << " x " << sliceWidth[1] << " x " << sliceWidth[2] << " x " << endl;
  
 // unsigned char* new_label_values = new unsigned char[dimx*dimy*dimz];

  // fill the look up table
  for(int i = (-1);i < 2;i++) {
    for(int j = (-1);j < 2;j++) {
      for(int k = (-1);k < 2;k++) {
        distLookup[ (i + 1) * 9 + (j + 1) * 3 + k + 1 ] = sqrt(pow(sliceWidth[0] * abs(i),2) +
                                                    pow(sliceWidth[1] * abs(j),2) +
                                                    pow(sliceWidth[2] * abs(k),2));
      }
    }
  }
  // start by initializing the ICM
  #pragma omp parallel for
  for(int z = 0; z < dimz ; z++) {
    // allocate
    float *posteriorProb = new float[numberOfLabels - 1];

    for(int y = 0; y < dimy; y++) {
      for(int x = 0;x < dimx; x++) {
        if(getVoxelValue(mask,x,y,z) > 0.5) {
          // collectValuesFromImagePP(labelLikelihoods,voxelProb,x,y,z,numberOfLabels - 1);
          for(char l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = 0.0;
          }
          for(int r = 0;r < numberOfRegions;r++) {
            for(char l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = posteriorProb[l]  +  getVoxelValue(&atlasImages[r],x,y,z) * (exp ( beta1 * log (getProb(mixture,r,(l + 1)) + 0.0001)));
            }
          }
         
          for(char l = 0;l < (numberOfLabels - 1);l++) {
            posteriorProb[l] = posteriorProb[l] * getVoxelValue(&labelLikelihoods[l],x,y,z);
            // posteriorProb[l] = posteriorProb[l] * voxelProb[l]; // MAP init
	    //  posteriorProb[l] = voxelProb[l]; // mlinit 
          }
          char newLabel = maxArg(posteriorProb,numberOfLabels - 1) + 1;
          putLabelValue(labels,x,y,z,newLabel);
          //new_label_values[dimx*dimy*z+dimx*y+x] = newLabel;
	 
        }
        else {
          putLabelValue(labels,x,y,z,0);
          //new_label_values[dimx*dimy*z+dimx*y+x] = 0;
        }
      }
    }
    // free memory
    delete[] posteriorProb;
  }
  
 
  // start the ICM
  while(changed && (iteration < maxIterations)) {
    changed = false;
    iteration++;
    if( verbose ) cout << "iteration " << iteration << endl;

    for(int odd = 0; odd < 2; odd++) {
    #pragma omp parallel for
    for(int z = odd; z < dimz ; z += 2) {
      // allocate
      float voxelProb[numberOfLabels - 1];
      float gibbsProb[numberOfLabels - 1];
      float posteriorProb[numberOfLabels - 1];

      for(int y = 0; y < dimy; y++) {
        for(int x = 0;x < dimx; x++) {
          if(getVoxelValue(mask,x,y,z) > 0.5) {
            //  compute the second order term in the prior, non-normalized
            // this is the same for every region.
            for(char l = 0;l < (numberOfLabels - 1);l++) {
              gibbsProb[l] = secondOrderGibbs(l + 1,labels,mixture->patlas->mrfConstants[l + 1].data(),
                                              x,y,z,distLookup,beta2);
	     
            }
	   
            // then start to calculate the posterior probabilities region by region
            // initialize label probabilities to 0
            for(char l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] = 0.0;
            }
            // then compute region wise prior
            for(int r = 0;r < numberOfRegions;r++) {
              for(char l = 0;l < (numberOfLabels - 1);l++) {
                voxelProb[l] = ( gibbsProb[l] ) * (exp ( beta1 * log (getProb(mixture,r,(l + 1)) + 0.0001)));
              }
	      //  normalize(voxelProb,numberOfLabels - 1); this is the difference between Markov and Gibbs
              for(char l = 0;l < (numberOfLabels - 1);l++) {
                posteriorProb[l] = posteriorProb[l] 
                                 + getVoxelValue(&atlasImages[r],x,y,z) * voxelProb[l];
              }    
            }
            
	   
            // now the prior probability is computed and it remains to 
            // multiply it with the likelihood term to get the posterior
            for(char l = 0;l < (numberOfLabels - 1);l++) {
              posteriorProb[l] =  posteriorProb[l] * getVoxelValue(&labelLikelihoods[l],x,y,z);
            }
            // then just find the minimum posterior and update the label
            unsigned char newLabel = maxArg(posteriorProb,(numberOfLabels - 1)) + 1;
            #pragma omp critical
            {
            if(newLabel != getLabelValue(labels,x,y,z)) {
              changed = true;
              putLabelValue(labels,x,y,z,newLabel);  
            }
            }
            //new_label_values[dimx*dimy*z+dimx*y+x] = maxArg(posteriorProb,(numberOfLabels - 1)) + 1;
          } // end if inside mask
        }   // end for x
      }     // end for y
    }       // end for z
    }       // end for odd
  }         // endwhile
 
  // delete[] new_label_values;
  return(iteration);
}

// Computes the ICM algorithm when the Atlas priors are used. Assumes that 
// only the pure labels can exist. 
int computeGibbsAtlas(AnalyzeLabelImage* labels,MixtureSpec* mixture,AnalyzeImage* mask,std::vector<AnalyzeImage> & labelLikelihoods, std::vector<AnalyzeImage> & atlasImages, std::vector<AnalyzeImage> & tissueProbMaps, float beta1, float beta2, int maxIterations, bool verbose)
{
  int x,y,z;
  int i,j,k,r;
  char l;
  int dimx,dimy,dimz;
  int numberOfLabels,numberOfRegions,pureLabels;
  float* voxelProb;
  float* gibbsProb;
  float* posteriorProb;
  bool changed = true;
  int iteration = 0;
  float distLookup[27];
  float sliceWidth[3], sliceWidthMin;
  char newLabel;
  AnalyzeImage tmpimg;

  bool testverbose;

 
  if(verbose) cout << "Initializing... " << endl;
  
  // save few commonly needed values  
  dimx = mask->header.x_dim;
  dimy = mask->header.y_dim;  
  dimz = mask->header.z_dim;
  numberOfLabels = mixture->patlas->numberOfLabels;
  numberOfRegions =  mixture->patlas->n;
  sliceWidthMin = MIN(fabs(labels->header.x_size),fabs(labels->header.y_size));
  sliceWidthMin = MIN(sliceWidthMin,fabs(labels->header.z_size));

  sliceWidth[0] = (labels->header.x_size) / sliceWidthMin;
  sliceWidth[1] = (labels->header.y_size) / sliceWidthMin;
  sliceWidth[2] = (labels->header.z_size) / sliceWidthMin;

  cout << "voxel dimensions:" << sliceWidth[0] << " x " << sliceWidth[1] << " x " << sliceWidth[2] << " x " << endl;
  pureLabels = countPureLabels(mixture->patlas);
 
  // allocate
  voxelProb = new float[pureLabels];
  gibbsProb = new float[pureLabels];  
  posteriorProb = new float[pureLabels];

  

  // fill the look up table
  for(i = (-1);i < 2;i++) {
    for(j = (-1);j < 2;j++) {
      for(k = (-1);k < 2;k++) {
        distLookup[ (i + 1) * 9 + (j + 1) * 3 + k + 1 ] = sqrt(pow(sliceWidth[0] * abs(i),2) +
                                                    pow(sliceWidth[1] * abs(j),2) +
                                                    pow(sliceWidth[2] * abs(k),2));
      }
    }
  }
  // start by initializing the ICM
  for(x = 0; x < dimx ; x++) {
    for(y = 0; y < dimy; y++) {
      for(z = 0;z < dimz; z++) {
	
        if(getVoxelValue(mask,x,y,z) > 0.5) {
          // collectValuesFromImagePP(labelLikelihoods,voxelProb,x,y,z,pureLabels);
          for(l = 0;l < (pureLabels);l++) {
            posteriorProb[l] = 0.0;
          }
/*
          // TODO: This does not seem to be used!!
          for(r = 0;r < numberOfRegions;r++) {
            for(l = 0;l < (pureLabels);l++) {
              posteriorProb[l] = posteriorProb[l]  +  getVoxelValue(&atlasImages[r],x,y,z) * (exp ( beta1 * log (getVoxelValue(&tissueProbMaps[l],x,y,z) + 0.0001)));
            }
          }
*/
	 
          for(l = 0;l < (pureLabels);l++) {
	    // posteriorProb[l] = posteriorProb[l] * voxelProb[l]; // MAP init
	    // posteriorProb[l] = voxelProb[l]; // mlinit 
	    posteriorProb[l] = getVoxelValue(&labelLikelihoods[l],x,y,z); // mlinit
	  }
          putLabelValue(labels,x,y,z,maxArg(posteriorProb,pureLabels) + 1);
	  
        }
        else {
          putLabelValue(labels,x,y,z,0);
         
        }
      }
    }
  }
  
 
  while(changed && (iteration < maxIterations)) {
    changed = false;
    iteration++;
    if( verbose ) cout << "iteration " << iteration << endl;
    for(x = 0; x < dimx ; x++) {
      for(y = 0; y < dimy; y++) {
        for(z = 0;z < dimz; z++) {
          if(getVoxelValue(mask,x,y,z) > 0.5) {
            
            //  compute the second order term in the prior, non-normalized
            // this is the same for every region.
            for(l = 0;l < (pureLabels);l++) {
              gibbsProb[l] = secondOrderGibbs(l + 1,labels,mixture->patlas->mrfConstants[l + 1].data(),
                                              x,y,z,distLookup,beta2);
	      
            }
	   
            // then start to calculate the posterior probabilities region by region
            // initialize label probabilities to 0
            for(l = 0;l < (pureLabels);l++) {
              posteriorProb[l] = 0.0;
            }
            // then compute region wise prior
            for(r = 0;r < numberOfRegions;r++) {
              for(l = 0;l < (pureLabels);l++) {
                voxelProb[l] = ( gibbsProb[l] ) * (exp ( beta1 * log (getVoxelValue(&tissueProbMaps[l],x,y,z) + 0.0001)));
              }
	      //  normalize(voxelProb,numberOfLabels - 1); this is the difference between Markov and Gibbs
              for(l = 0;l < (pureLabels);l++) {
                posteriorProb[l] = posteriorProb[l] 
                                 + getVoxelValue(&atlasImages[r],x,y,z) * voxelProb[l];
              }    
            }
            
	    
            // now the prior probability is computed and it remains to 
            // multiply it with the likelihood term to get the posterior
            for(l = 0;l < (pureLabels);l++) {
              posteriorProb[l] =  posteriorProb[l] * getVoxelValue(&labelLikelihoods[l],x,y,z);
            }
            // then just find the minimum posterior and update the label
            newLabel = maxArg(posteriorProb,(pureLabels)) + 1;
            if(newLabel != getLabelValue(labels,x,y,z)) {
              changed = true;
              putLabelValue(labels,x,y,z,newLabel);  
            }
          } // endif
        }   // end for z
      }     // end for y
    }       // end for x
  }         // endwhile

  // free memory
  delete[] voxelProb; 
  delete[] gibbsProb; 
  delete[] posteriorProb;
 
  return(iteration);
}
// Computes the ICM algorithm when the 
// only the pure labels can exist. 
/*
int computeGibbsPure(AnalyzeLabelImage* labels,MixtureSpec* mixture,AnalyzeImage* mask,AnalyzeImage** labelLikelihoods, AnalyzeImage** atlasImages, float beta1, float beta2, int maxIterations,bool verbose)
{
  int x,y,z;
  int i,j,k,r;
  char l;
  int dimx,dimy,dimz;
  int numberOfLabels,numberOfRegions,pureLabels;
  float* voxelProb;
  float* gibbsProb;
  float* posteriorProb;
  bool changed = true;
  int iteration = 0;
  float distLookup[27];
  float sliceWidth[3], sliceWidthMin;
  char newLabel;
  AnalyzeImage tmpimg;

  bool testverbose;

  // int count = 0;
  if(verbose) cout << "Initializing... " << endl;
  
  // save few commonly needed values  
  dimx = mask->header.x_dim;
  dimy = mask->header.y_dim;  
  dimz = mask->header.z_dim;
  numberOfLabels = mixture->patlas->numberOfLabels;
  numberOfRegions =  mixture->patlas->n;
  sliceWidthMin = MIN(fabs(labels->header.x_size),fabs(labels->header.y_size));
  sliceWidthMin = MIN(sliceWidthMin,fabs(labels->header.z_size));

  sliceWidth[0] = (labels->header.x_size) / sliceWidthMin;
  sliceWidth[1] = (labels->header.y_size) / sliceWidthMin;
  sliceWidth[2] = (labels->header.z_size) / sliceWidthMin;
  cout << sliceWidth[0] << " " << sliceWidth[1] << " " << sliceWidth[2] << " " << endl;
 
  pureLabels = countPureLabels(mixture->patlas);
 
  // allocate
  voxelProb = new float[pureLabels];
  gibbsProb = new float[pureLabels];  
  posteriorProb = new float[pureLabels];

  

  // fill the look up table
  for(i = (-1);i < 2;i++) {
    for(j = (-1);j < 2;j++) {
      for(k = (-1);k < 2;k++) {
        distLookup[ (i + 1) * 9 + (j + 1) * 3 + k + 1 ] = sqrt(pow(sliceWidth[0] * abs(i),2) +
                                                    pow(sliceWidth[1] * abs(j),2) +
                                                    pow(sliceWidth[2] * abs(k),2));
      }
    }
  }
  // start by initializing the ICM
  for(x = 0; x < dimx ; x++) {
    for(y = 0; y < dimy; y++) {
      for(z = 0;z < dimz; z++) {

        if(getVoxelValue(mask,x,y,z) > 0.5) {
          collectValuesFromImagePP(labelLikelihoods,voxelProb,x,y,z,pureLabels);
          for(l = 0;l < (pureLabels);l++) {
            posteriorProb[l] = 0.0;
          }
          for(r = 0;r < numberOfRegions;r++) {
            for(l = 0;l < (pureLabels);l++) {
               posteriorProb[l] = posteriorProb[l]  +  getVoxelValue(atlasImages[r],x,y,z) * (exp ( beta1 * log (getProb(mixture,r,(l + 1)) + 0.0001)));
            
            }
          }
	 
          for(l = 0;l < (pureLabels);l++) {
	    posteriorProb[l] = posteriorProb[l] * voxelProb[l]; // MAP init
	    //  posteriorProb[l] = voxelProb[l]; // mlinit 
	  }
          putLabelValue(labels,x,y,z,maxArg(posteriorProb,pureLabels) + 1);
	  
        }
        else {
          putLabelValue(labels,x,y,z,0);
         
        }
      }
    }
  }
  
 
  while(changed && (iteration < maxIterations)) {
    changed = false;
    iteration++;
    if( verbose ) cout << "iteration " << iteration << endl;
    for(x = 0; x < dimx ; x++) {
      for(y = 0; y < dimy; y++) {
        for(z = 0;z < dimz; z++) {
          if(getVoxelValue(mask,x,y,z) > 0.5) {
            
            //  compute the second order term in the prior, non-normalized
            // this is the same for every region.
            for(l = 0;l < (pureLabels);l++) {
              gibbsProb[l] = secondOrderGibbs(l + 1,labels,mixture->patlas->mrfConstants[l + 1],
                                              x,y,z,distLookup,beta2);
	      
            }
	   
            // then start to calculate the posterior probabilities region by region
            // initialize label probabilities to 0


            for(l = 0;l < (pureLabels);l++) {
              posteriorProb[l] = 0.0;
            }
            // then compute region wise prior
            for(r = 0;r < numberOfRegions;r++) {
              for(l = 0;l < (pureLabels);l++) {
                voxelProb[l] = ( gibbsProb[l] ) * (exp ( beta1 * log (getProb(mixture,r,(l + 1)) + 0.0001)));
              }
	      //  normalize(voxelProb,numberOfLabels - 1); this is the difference between Markov and Gibbs
              for(l = 0;l < (pureLabels);l++) {
                posteriorProb[l] = posteriorProb[l] 
                                 + getVoxelValue(atlasImages[r],x,y,z) * voxelProb[l];
              }    
            }            
	    
            // now the prior probability is computed and it remains to 
            // multiply it with the likelihood term to get the posterior
            for(l = 0;l < (pureLabels);l++) {
              posteriorProb[l] =  posteriorProb[l] * getVoxelValue(labelLikelihoods[l],x,y,z);
            }
            // then just find the minimum posterior and update the label
            newLabel = maxArg(posteriorProb,(pureLabels)) + 1;
            if(newLabel != getLabelValue(labels,x,y,z)) {
              changed = true;
              putLabelValue(labels,x,y,z,newLabel);  
            }
          } // endif
        }   // end for z
      }     // end for y
    }       // end for x
  }         // endwhile

  // free memory
  delete[] voxelProb; 
  delete[] gibbsProb; 
  delete[] posteriorProb;
 
  return(iteration);
}
*/


int convertPVElabels(AnalyzeLabelImage* crispLabels, AnalyzeLabelImage* pveLabels, AnalyzeImage* img, std::vector<AnalyzeImage> & atlasImages, MixtureSpec* mixture)
{
 

  int n,r;
  int i,j,k;
  char label,pureLabel1,pureLabel2;;
  float likelihoodValue[INTERVALS + 1];
  float rprob,mean1,mean2,var1,var2;;
  float maxindex,maxvalue;
  float t;
  
 

  for(i = 0;i < pveLabels->header.x_dim;i++) {
    for(j = 0;j < pveLabels->header.y_dim;j++) {
      for(k = 0;k < pveLabels->header.z_dim;k++) {
        label = getLabelValue(pveLabels,i,j,k);
        if(mixture->patlas->labelTypes[label].pureLabel) {
          putLabelValue(crispLabels,i,j,k,label);
        }
        else {
          pureLabel1 = mixture->patlas->labelTypes[label].mixed[0]; 
          pureLabel2 = mixture->patlas->labelTypes[label].mixed[1];
          if(pureLabel1 == 0) {
             putLabelValue(crispLabels,i,j,k,pureLabel2);
          }
          else if(pureLabel2 == 0) {
            putLabelValue(crispLabels,i,j,k,pureLabel1);
          }
          else {
            for(n = 0;n < (INTERVALS + 1);n++) {
              likelihoodValue[n] = 0.0;
            }
            for(r = 0;r < mixture->patlas->n;r++) {
              rprob = getVoxelValue(&atlasImages[r],i,j,k);
              if(rprob > 0.001) {
	        mean1 = getMu(mixture,r,pureLabel1);
	        mean2 = getMu(mixture,r,pureLabel2);          
                var1 = getSigma2(mixture,r,pureLabel1);
                var2 = getSigma2(mixture,r,pureLabel2);
                for(n = 0;n < (INTERVALS + 1);n++) {
		  t = n * (float) 1/ INTERVALS;
                  likelihoodValue[n] = likelihoodValue[n] + rprob * 
		                       computeGaussianLikelihood(getVoxelValue(img,i,j,k),
                                                               t * mean1 + ( 1 - t) * mean2,
                                                               t * t * var1 + (1 - t) * ( 1-  t) *var2);
                }
              }
            }
            maxvalue = likelihoodValue[0];
            maxindex = 0;
            for(n = 0;n < (INTERVALS + 1) ;n++) {
	      if(likelihoodValue[n] > maxvalue) {
                maxvalue = likelihoodValue[n];
                maxindex = n;
              }
	    } 
            if(maxindex > (INTERVALS/2)) {
              putLabelValue(crispLabels,i,j,k,pureLabel1);   
            }   
            else {
              putLabelValue(crispLabels,i,j,k,pureLabel2);
            }
          } // end else
        } // end else
      }   // end for k
    }     // end for j
  }       // end for i
  return(0);
}
