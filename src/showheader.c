#include "analyze.h"

// Displays some information Contained in the header of the analyze image.
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

int main(int argc,char** argv) 
{
  AnalyzeHeader hdr;
  bool bs;

  if(argc < 2) {
    cout << "Usage: showheader analyzefile" << endl;
    return(1);
  } 
  if(readHeader(argv[1],hdr) == false) {
    cout << "File does not exist" << endl; 
    return(1);
  } 
  bs = byteSwapNecessary(&hdr);
  if(bs) swapHeader(hdr);
  cout << "datatype: " << hdr.datatype<< endl; 
  cout << "dimensions: " << hdr.x_dim << " " << hdr.y_dim << " " << hdr.z_dim << endl;
  cout << "voxel size: " << hdr.x_size << " " << hdr.y_size << " " << hdr.z_size << endl;
  cout << "origin: " << hdr.x_orig << " " << hdr.y_orig << " " << hdr.z_orig << endl;
  return(0); 

}
