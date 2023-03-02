// ******************************************************
// Functions for SVPA based non homogeneous MRFs 
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

#ifndef NHMRF_H
#define NHMRF_H

#include "analyze.h"
#include "atlasspec.h"

#define VERY_SMALL 1E-16
#define INTERVALS 100

#define DEFAULT_PVESTART 0.0
#define DEFAULT_PVEEND 1.0
#define DEFAULT_BETA1 0.0
#define DEFAULT_BETA2 0.05

struct SvpasegParameters {
  float beta1;
  float beta2;
  bool writePveLabelImage;
  char* pveLabelImage;
  bool markov;
};

struct MixtureSpec { 
  AtlasSpec* patlas;  // remember that forbidden labels will mean that the corresponding probability is zero
  float* prob;        //  prior probabilities
  float* mu;          //  means, 
  float* sigma2;      //   variances
};

// then a few inline helper functions
// **************************************************
// allocates a mixtureSpec and ties it to a particular atlas

inline void allocateMixtureSpec(AtlasSpec* atlas, MixtureSpec* mixture)
{
  mixture->patlas = atlas;
  mixture->prob = new float[(atlas->n) * (atlas -> numberOfLabels)];
  mixture->mu = new float[(atlas->n) * (atlas -> numberOfLabels)];
  mixture->sigma2 = new float[(atlas->n) * (atlas -> numberOfLabels)];   
  
};

inline void copyMixtureSpec(MixtureSpec* source,MixtureSpec* target) 
{
  int i;
  for(i = 0;i < ((source->patlas->n) * (source->patlas -> numberOfLabels)); i++) {
    target->prob[i] = source->prob[i];
    target->mu[i] = source->mu[i];
    target->sigma2[i] = source->sigma2[i]; 
  }
}

inline float getProb(MixtureSpec* mixture,int region,int label)
{
  return(mixture->prob[(mixture->patlas->numberOfLabels)*region + label]);  
};

inline void putProb(MixtureSpec* mixture,int region,int label, float val)
{
  mixture->prob[(mixture->patlas->numberOfLabels)*region + label] = val;  
};

inline float getMu(MixtureSpec* mixture,int region,int label)
{
  return(mixture->mu[(mixture->patlas->numberOfLabels)*region + label]);  
};

inline void putMu(MixtureSpec* mixture,int region,int label, float val)
{
  mixture->mu[(mixture->patlas->numberOfLabels)*region + label] = val;  
};

inline float getSigma2(MixtureSpec* mixture,int region,int label)
{
  return(mixture->sigma2[(mixture->patlas->numberOfLabels)*region + label]);  
};

inline void putSigma2(MixtureSpec* mixture,int region,int label, float val)
{
  mixture->sigma2[(mixture->patlas->numberOfLabels)*region + label] = val;  
};

inline void printMixture(MixtureSpec* mixture) 
{
  int i,j;
  for(i = 0;i < mixture->patlas->n;i++) {
    if(i > 1) { 
      cout << mixture->patlas->regionnames[i] << endl;
    }
    for(j = 0;j < mixture->patlas->numberOfLabels;j++) {
      if(mixture->patlas->labelTypes[j].pureLabel)
        cout << getMu(mixture,i,j) << " " << getSigma2(mixture,i,j) << " " << getProb(mixture,i,j) << endl;
      else
        cout << getProb(mixture,i,j) << endl;
    }
    cout << "----------------------------------------" << endl;
  }
}

// Normalizes n probalities to sum to one

inline void normalize(float* pval, char n)

{ float sca = 0.0;
  char i;

  for(i = 0;i < n;i++) {
    sca = sca + pval[i];
  }
  if(fabs(sca) >  VERY_SMALL) {       // To avoid divisions by zero 
    for(i = 0;i < n;i++) {
      pval[i] = pval[i]/sca;
    }
  }
}

// Finds maximum argument out of the n possibilities
// moved to analyze.h
// inline char maxArg(float* pval,int n)
// {
//  float maximum;
//  int i,index;
//  
//  maximum = pval[0];
//  index = 0;
//  for(i = 1;i < n;i++) {
//    if(pval[i] > maximum) {
//      index = i;
//      maximum = pval[i];
//    }
//  }
//  return( (char) index);
// }

/*
inline void collectValuesFromImagePP(std::vector<AnalyzeImage> & imagePP,float* collectHere,
                                     int x, int y, int z,int n) 
{
  int i;
  for(i = 0;i < n;i++) {
    collectHere[i] =getVoxelValue(&imagePP[i],x,y,z);
  } 
}
*/

// Computes likelihood of value given parameters mean and variance. 
// Returns the likelihood.

inline float computeGaussianLikelihood(float value, float mean , float var)

{ 
  return(exp(-((value - mean) *(value - mean))/(2 * var))/(sqrt(2 * M_PI * var)));

}

// Computes the likelihoods for the mixed classes. Returns the likelihood.
// var1,var2 are the variances of pdfs representing pure classes.
// So the model for the variable y (representing the 
// intensity value) that is composed of t * tissue1 and (1 - t)* tissue2 becomes :

// y = t*x1 + (1 - t)*x2 ,
// x1 ~ N(mean1,var1) , x2 ~ N(mean2,var2) 

inline float computeMarginalizedLikelihood(float value, float mean1 , float mean2, 
                                       float var1, float var2, 
					   unsigned int nof_intervals, float start, float end)

{ 
  float lh, tmean , tvar, t, interval_len;
  int i;  
  
  interval_len = (float) (end - start) / nof_intervals;
  lh = 0;
  for(i = 0; i < nof_intervals; i++) {
    t = (i + 0.5) * interval_len + start;
    tmean = t * mean1 + ( 1 - t ) * mean2;
    tvar = pow(t,2) * var1 + pow((1 - t),2) * var2;
    lh = lh + computeGaussianLikelihood(value, tmean,tvar) / nof_intervals;
  }
  return(lh);
}

inline float secondOrderGibbs(char testLabel,AnalyzeLabelImage* labels,float* mrfConstants,
                              int x, int y,int z,float* distanceLookup, float beta)
{
  int i,j,k;
  float exponent = 0.0;
  float tmp;

  for(i = (-1);i < 2;i++) {
    for(j = (-1);j < 2;j++) {
      for(k = (-1);k < 2;k++) {
        if(! ((i == 0) && (j == 0) && (k == 0)) ) {
          tmp = mrfConstants[getSafeLabelValue(labels,x + i,y + j, z + k)]/distanceLookup[ (i + 1) * 9  + (j + 1) * 3 + (k + 1) ];
          exponent = exponent + tmp;
	  
        }
      }
    }
  }  
  return(exp(- (beta * exponent)));
}


int readMixtureParameters(char* filename,AtlasSpec* atlas,MixtureSpec* mixture);
int writeMixtureParameters(char* filename,AtlasSpec* atlas,MixtureSpec* mixture,bool overwrite);
int parseParamsSvpaseg(SvpasegParameters* param,int n,char** arguments,AtlasSpec* atlas);

int computeVoxelLikelihood(MixtureSpec* mixture,AnalyzeImage* img,AnalyzeImage* mask,std::vector<AnalyzeImage> & atlasImages,std::vector<AnalyzeImage> & labelLikelihoods);

//int computeMRF(AnalyzeLabelImage* labels,MixtureSpec* mixture,AnalyzeImage* mask,AnalyzeImage** labelLikelihoods, AnalyzeImage** atlasImages, float beta1, float beta2,int maxIterations, bool verbose);

// the next three functions do the same under 3 slightly different settings and 
// they should be combined

int computeGibbs(AnalyzeLabelImage* labels,MixtureSpec* mixture, AnalyzeImage* mask, std::vector<AnalyzeImage> & labelLikelihoods, std::vector<AnalyzeImage> & atlasImages, float beta1,float beta2,int maxIterations, bool verbose );

int computeGibbsAtlas(AnalyzeLabelImage* labels,MixtureSpec* mixture,AnalyzeImage* mask,std::vector<AnalyzeImage> & labelLikelihoods, std::vector<AnalyzeImage> & atlasImages, std::vector<AnalyzeImage> & tissueProbMaps, float beta1, float beta2, int maxIterations, bool verbose);

//int computeGibbsPure(AnalyzeLabelImage* labels,MixtureSpec* mixture, AnalyzeImage* mask,AnalyzeImage** labelLikelihoods, AnalyzeImage** atlasImages, float beta1,float beta2,int maxIterations, bool verbose );
//{
// };
 
int convertPVElabels(AnalyzeLabelImage* crispLabels, AnalyzeLabelImage* pveLabels, AnalyzeImage* img, std::vector<AnalyzeImage> & atlasImages, MixtureSpec* mixture);


#endif
