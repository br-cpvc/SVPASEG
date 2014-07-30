// 
// ******************************************************
// Spatial translation of Analyze images
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
  AnalyzeImage img,timg;
  int intstatus;
  int i,j,k;
  int xtrans,ytrans,ztrans,xstart,ystart,zstart,xend,yend,zend;
  bool status;

  if( argc < 6) {
    cout << "Usage: translate_analyze inputfile outputfile x_trans y_trans z_trans" << endl;
    cout << "Images are assumed to be transaxial." << endl;
    cout << "Translations are given relative to voxels and not in millimeters" << endl;   
  }
  intstatus = readImage(argv[1],&img);
  if(intstatus != 0) {
    cout << "Could not read the file" << argv[1] << endl;
    return(2);
  }
  xtrans = atoi(argv[3]);
  ytrans = atoi(argv[4]);
  ztrans = atoi(argv[5]);

  if(xtrans < 0) {
    xstart = -xtrans;
    xend = img.header.x_dim;
  }
  else {  
    xstart = 0;
    xend = img.header.x_dim - xtrans;
  }   
  if(ytrans < 0) {
    ystart = -ytrans;
    yend = img.header.y_dim;
  }
  else {  
    ystart = 0;
    yend = img.header.y_dim - ytrans;
  }   
  if(ztrans < 0) {
    zstart = -ztrans;
    zend = img.header.z_dim;
  }
  else {  
    zstart = 0;
    zend = img.header.z_dim - ztrans;
  }   
 
  status = newImage(&timg,&img);
  for(i = xstart;i < xend;i++) {  
    for(j = ystart;j < yend;j++) {
      for(k = zstart;k < zend;k++) {
        putVoxelValue(&timg,i + xtrans ,j + ytrans,k + ztrans,getVoxelValue(&img,i,j,k));
      }
    }
  }
  intstatus = writeImage(argv[2],&timg,true);
  if(intstatus != 0) {
    cout << "Could not write the file" << argv[2] << endl;
    return(3);
  }
  return(0);
}
