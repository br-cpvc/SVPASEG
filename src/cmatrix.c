// 
// ******************************************************
// Confusion matrices and related stuff for Analyze images
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

int main(int argc,char** argv)
{
  AnalyzeImage tmp;
  AnalyzeLabelImage img1,img2;
  int i,j,k;
  bool jaccardStatus,completeStatus,mcrStatus;
  int dimx,dimy,dimz;
  unsigned char max1,max2;
  int tmpint,imgsize;
  int intstatus;
  int** confMatrix;
  float jaccard,mcr;
  int a,b,c; 
  unsigned char x,y;

  if(argc < 2) {
    cout << "Usage: cmatrix img1 img2 [jaccard/mcr]" << endl;
    return(1);
  }
  
  intstatus = readLabelImage(argv[1],&img1);
  if(intstatus != 0) {
    cout << "File " << argv[1] <<  "could not be opened " << intstatus << endl;
    return(2);
  }
  intstatus = readLabelImage(argv[2],&img2);
  if((intstatus != 0) && (intstatus < 4)) {
    cout << "File " << argv[2] << " could not be opened " << intstatus << endl;
    return(2);
  }
  if(intstatus == 4) {
    intstatus = readImage(argv[2],&tmp);
    if(intstatus != 0) {
      cout << "File " << argv[2] <<  "could not be opened " << intstatus << endl;
      return(2);
    }    
    newLabelImage(&img2,&tmp);
    imgsize = tmp.header.x_dim*tmp.header.y_dim*tmp.header.z_dim;
    for(i = 0;i < imgsize;i++) {
      img2.data[i] = (char) (ceil((tmp.data[i]) - 0.5));
    }
  }

  // check that the dimensions match
 
  dimx = img1.header.x_dim;
  dimy = img1.header.y_dim;
  dimz = img1.header.z_dim;

  if(!( (dimx == img2.header.x_dim) && (dimy == img2.header.y_dim) && (dimz == img2.header.z_dim))) {
    cout << "Dimensions did not match" << endl;
    return(3);
  }
  
  jaccardStatus = true;
  mcrStatus = true;
  completeStatus = true;

  if(argc > 3) {
    if(!strcmp("jaccard",argv[3])) {
      mcrStatus = false;
      completeStatus = false;
    }
    if(!strcmp("mcr",argv[3])) {
      jaccardStatus = false;
      completeStatus = false;
    }
  }
 
  max1 = labelMax(&img1);
  max2 = labelMax(&img2);
 

  if(max1 != max2) {
    if(max1 == 255) {
      for(i = 0;i < dimx*dimy*dimz;i++) {
        tmpint = img1.data[i] * max2;
        img1.data[i] = tmpint / 255;
      }
      max1 = labelMax(&img1);
    }
    if(max2 == 255) {
      for(i = 0;i < dimx*dimy*dimz;i++) {
        tmpint = img2.data[i] * max1;
        img2.data[i] = tmpint / 255;
      }
      max2 = labelMax(&img2);
    }
  }
  max1 = MAX(max1,max2);
  max2 = MAX(max1,max2);
  

  confMatrix = new int*[max1 + 1];
  for(i = 0;i < (max1 + 1);i++) {
    confMatrix[i] = new int[max2 + 1];
  }
  for(i = 0;i < (max1 + 1);i++) {
    for(j = 0;j < (max2 + 1);j++) {
      confMatrix[i][j] = 0;
    }
  }
 
  for(i = 0;i < dimx;i++) {
    for(j = 0;j < dimy;j++) {
      for(k = 0;k < dimz;k++) {
        x = getLabelValue(&img1,i,j,k);
	y = getLabelValue(&img2,i,j,k);
        confMatrix[x][y] = confMatrix[x][y] + 1;
      }
    }
  }
 
  if(completeStatus) {
    cout << "Confusion matrix" << endl;
    for(i = 0;i < (max1 + 1);i++) {
      for(j = 0;j < (max2 + 1);j++) {
        cout << confMatrix[i][j] << " ";
      }
      cout << endl; 
    }
    cout << "Jaccard" << endl;
  }
  if(jaccardStatus) {
    for(i = 0;i <  ( max1 + 1);i++) {
      a =  confMatrix[i][i];
      b = -a;
      for(j = 0;j <  ( max1 + 1);j++) {
        b = b + confMatrix[j][i];
      }
      c = -a;
      for(j = 0;j <  ( max2 + 1);j++) {
        c = c + confMatrix[i][j];
      } 
      jaccard = (float) a/( (float) (a + b + c));
      cout << jaccard;
      cout << " ";  
    }
    cout << endl;
  }
  if(mcrStatus) {
    a = 0;
    for(i = 1;i <  ( max1 + 1);i++) {
      a =  a + confMatrix[i][i];
    }
    b = 0;
    for(i = 1;i <  ( max1 + 1);i++) {
      for(j = 1;j <  ( max1 + 1);j++) {
        b =  b + confMatrix[i][j];
      }
    }
    mcr = (float) a / (float) b;
    cout << mcr;
    cout << endl;
  }
  return(0);

}
