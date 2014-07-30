// 
// ******************************************************
// Relabeling a labeled Analyze image
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
  AnalyzeImage imgin,imgout;
  int x,y,z,labelsToChange,i;
  int intstatus;  
  int inlabel,outlabel;
  bool clobber = true;
  bool status;
  div_t divtmp;

  if(argc < 5) {
    cout << "Usage: relabel_analyze inputfile outputfile inlabel1 outlabel1" << endl;
    return(1);
  }
  intstatus = readImage(argv[1],&imgin);
  if(intstatus != 0) {
    cout << "Could not read the file" << argv[1] << endl;
    return(2);
  }
  status = copyImage(&imgin,&imgout);
  divtmp = div(argc - 3,2);
  if(divtmp.rem == 0) {
    labelsToChange = divtmp.quot;
    
  }
  else {
    cout << " wrong number of inputs" << endl;
    return(3);
  }
  for(i = 0;i < labelsToChange;i++) {
    inlabel = atoi(argv[2*i + 3]);
    outlabel = atoi(argv[2*i + 4]);
    cout << "Changing" << inlabel << "->" << outlabel << endl;
    for(x = 0;x < imgin.header.x_dim;x++) {
      for(y = 0;y < imgin.header.y_dim;y++) {
        for(z = 0;z < imgin.header.z_dim;z++) {
          if(fabs(getVoxelValue(&imgin,x,y,z) - inlabel) < 0.0001)
	    putVoxelValue(&imgout,x,y,z,outlabel);
        }
      }
    }
  }
  intstatus = writeImage(argv[2],&imgout,clobber);
  if(intstatus != 0) {
    cout << "Could not write the file" << argv[2] << endl;
    return(3);
  }
  return(0);



}
