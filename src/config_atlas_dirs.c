// ******************************************************
// A simple configuration function for atlas specification files 
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

#include "atlasspec.h"

int main(int argc, char** argv) 
{
  
  bool boolstatus;
  AtlasSpec atlas;
  int intstatus;
  int i,purelabels;
  char atlasType;
  char filename[256];
  char pathname[256];
  char fullpath[256];


  if(argc < 3) {
    cout << "Usage: conf_atlas_dirs atlas_spec_file atlas_dir" << endl;
   
    return(0);
  }   

  boolstatus = readAtlasSpec(&atlas,argv[1]);
  if(boolstatus == false) {
    cout << "Could not read atlas file " << argv[1] << endl;
    return(1);
  }  
  if(atlas.onlyPureLabels) {
    if(atlas.useTPM) atlasType = 't';
    else atlasType = 'r';
  }
  else { 
    atlasType = 'p';
  } 
  strcpy(fullpath,argv[2]);
  if(fullpath[strlen(fullpath) - 1] != '/') strcat(fullpath,"/");

  for(i = 0;i < atlas.n;i++) {
    extractpath(filename,pathname,atlas.filenames[i]);
    strcpy(pathname,fullpath);
    strcat(pathname,filename);
    strcpy(atlas.filenames[i],pathname);
  }
  if(atlas.useTPM) {
    purelabels = countPureLabels(&atlas);
    for(i = 0;i < purelabels;i++) {
      extractpath(filename,pathname,atlas.TPMfilenames[i]);
      strcpy(pathname,fullpath);
      strcat(pathname,filename);
      strcpy(atlas.TPMfilenames[i],pathname);
    }
  }
 
  intstatus = writeAtlasSpec(&atlas,argv[1],atlasType,true);
  if(intstatus != 0) {
    cout << "Could not write atlas file " << argv[1] << " " << intstatus << endl;
    return(2);
  }  
  
  return(0);
}  
