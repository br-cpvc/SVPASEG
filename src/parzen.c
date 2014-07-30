#include "parzen.h"

// ******************************************************
// Parzen estimates for GAMIXTURE 
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

// computes the x coordinates for the estimate

void computeX(PdfEstimate* hatf,AnalyzeImage* img, AnalyzeLabelImage* mask,float percentage)
{
  float maximum,minimum,step;
  int i;  

  if(fabs(percentage - 1) < 0.00001) { 
    maximum = imageMax(img,mask);
  }
  else {
    maximum = imageKth(img,mask,percentage);
  }
  minimum = imageMin(img,mask);
  hatf->x = new float[hatf->n];
  step = (maximum - minimum)/((float) (hatf->n - 1));
  for(i = 0;i < hatf->n;i++) {
    hatf->x[i] = minimum +(i * step);
  }
}

void computeY(PdfEstimate* hatf,AnalyzeImage* img, AnalyzeImage* mask)
{
  int x,y,z,imgsize,j,enditer,startiter;
  float masksum;  
  float c;
  float beta;
  float step,ysum;

  hatf->y = new float[hatf->n];
 
  for(j = 0;j < hatf->n;j++) {
    hatf->y[j] = 0;
  }   

  step = hatf->x[1] - hatf->x[0]; // assumes equal x-axis spacing 
  beta = 2 * (hatf->sigma) * (hatf->sigma);
  c = sqrt(2 * M_PI) * (hatf->sigma);
  masksum = 0;
  ysum = 0;  

  // imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim);
  for(x = 0;x < img->header.x_dim;x++) {
    for(y = 0;y < img->header.y_dim;y++) {
      for(z = 0;z < img->header.z_dim;z++) {
       
        if(getVoxelValue(mask,x,y,z) > 0.001) {
	  masksum = masksum + getVoxelValue(mask,x,y,z);
          startiter = (int) floor((getVoxelValue(img,x,y,z) - 5*(hatf->sigma) - hatf->x[0])/step); // speeding up. Assumes equal
          enditer = (int) ceil((getVoxelValue(img,x,y,z) + 5*(hatf->sigma) - hatf->x[0])/step) + 1;   // x-axis spacing 
          startiter = MAX(startiter,0);
          enditer = MIN(enditer,hatf->n); 
          for(j = startiter;j < enditer;j++) {
            hatf->y[j] = hatf->y[j] +getVoxelValue(mask,x,y,z)*exp(-(pow((hatf->x[j] - getVoxelValue(img,x,y,z)),2)/beta));
	  }
	}
      }
    }
  }
  
  for(j = 0;j < hatf->n;j++) {
    hatf->y[j] = hatf->y[j]/(c*masksum);
    ysum = ysum + hatf->y[j];
  } 
  ysum = ysum*step;
  cout << "Parzen integrand " << ysum << endl;  
  // take care that the pdf intergrates to 1. Note that we study intergration range from x[0] - step/2 to x[n] + step/2 
  // This precaution is because the speed up
  for(j = 0;j < hatf->n;j++) {
    hatf->y[j] = hatf->y[j]/(ysum);
  } 
 
}

// writes pdf estimate to a text file
// format of the text file:
// 1st line : float sigma int n 
// 2nd - nth + 1 line: x[i] y[i] 

bool writeEstimate(char* filename,PdfEstimate* hatf, bool overwrite) 
{
  int i;
  
  if(! overwrite) {
    ifstream ifile(filename);
    if(ifile) {
      ifile.close();
      return(false);
    }
  }
  ofstream ofile(filename);
  if(!ofile) return(false);
  ofile << hatf->sigma << " " << hatf->n << "\n";
  for(i = 0;i < hatf->n;i++) {
    ofile << hatf->x[i] << " " << hatf->y[i] << "\n"; 
  }
  ofile.close();
  return(true);
}

bool readEstimate(char* filename,PdfEstimate* hatf)
{
  int i;

  ifstream ifile(filename);
  if(!ifile) return(false);
  ifile >> hatf->sigma;
  ifile >> hatf->n;
  hatf->x = new float[hatf->n];
  hatf->y = new float[hatf->n];
  for(i = 0;i < hatf->n;i++) {
    ifile >> hatf->x[i];
    ifile >> hatf->y[i];
  }
  return(true);
}
