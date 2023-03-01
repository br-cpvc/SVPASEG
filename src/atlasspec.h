// ******************************************************
// Functions for atlas specification files
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



#ifndef ATLASSPEC_H
#define ATLASSPEC_H

#include "analyze.h"
#include <vector>
#include <string>

#define MAX_NAME_LEN 256
#define MAX_LABELS 20 
#define PURELABEL  0
#define PVELABEL 1

struct LabelType {
  bool pureLabel;
  int mixed[2];  
};

struct LabelList {
  int len;
  int list[MAX_LABELS];
};

struct AtlasSpec {
  int n; // Gives the number of different fuzzy regions
  int numberOfLabels; // Gives the number of possible labels, background included 
   
  std::vector<std::string> labelnames;    // Gives the names of labels; 0 is reserved for background
  std::vector<std::string> filenames;    // Gives the names of the files containg fuzzy masks
  std::vector<std::string> regionnames;  // Gives the description of each region (optional)
  std::vector<LabelList> permittedLabels; // Labels permitted for each region 
                       
  std::vector<std::vector<float> > mrfConstants; // pairwise interactions in the mrf
  std::vector<LabelType> labelTypes;  // types of labels (pve or gaussian)
  std::vector<std::vector<float> > regionLowProb;
  std::vector<std::vector<float> > regionUpProb;  
  std::vector<std::string> TPMfilenames; // Gives the tissue probability maps for each region 
  bool onlyPureLabels;
  bool useTPM;                
};




bool readAtlasSpec(AtlasSpec* atlas,char* filename);
int writeAtlasSpec(AtlasSpec* atlas,char* filename, char atlasType, bool overwrite);
void freeAtlas(AtlasSpec* atlas);
int readAtlasImages(AtlasSpec* atlas,AnalyzeImage** atlasImages); 
int readAtlasImages(AtlasSpec* atlas,std::vector<AnalyzeImage> & atlasImages); 
void freeAtlasImages(AtlasSpec* atlas,AnalyzeImage** atlasImages);
bool maskAtlas(AtlasSpec* atlas, std::vector<AnalyzeImage> & atlasImages,AnalyzeImage* mask);
bool maskAtlas(AtlasSpec* atlas, AnalyzeImage** atlasImages,AnalyzeImage* mask);
int readTPMimages(AtlasSpec* atlas,std::vector<AnalyzeImage>& TPMImages, AnalyzeImage* mask ,int pureLabels);
void freeTPMimages(AtlasSpec* atlas,AnalyzeImage** TPMImages, int pureLabels);
// int processTPMimages(AtlasSpec* atlas, AnalyzeImage** atlasImages,AnalyzeImage* mask, int pureLabels);


inline void constructLabelList(LabelList* labels,int maxLabel,int* remThese, int howManyToRemove) 
{
  int i,j,k;
  bool addLabel;

  k = 0;
  for(i = 1;i < (maxLabel + 1);i++) {
    addLabel = true;
    for(j = 0;j < howManyToRemove;j++) {
      if(remThese[j] == i) addLabel = false;
    }
    if(addLabel) {
      labels->list[k] = i;
      k++;
    }
  }
  labels->len = k;
} 

// inline void invertLabelList(LabelList* newlabels,LabelList* labels,int maxLabel, bool inclBG)
// {}

inline void printLabelInfo(AtlasSpec* atlas) 
{
  int i,j;
  cout << "Use TPMs :" << atlas->useTPM << endl;
  cout << "Only pure labels:" << atlas->onlyPureLabels << endl;
  cout << "Number of regions:" << atlas->n << " Number of Labels:" << atlas->numberOfLabels << endl;
  cout << "Label:   Name   Pure   compos " << endl;
  for(i = 0; i < atlas->numberOfLabels;i++) {
    cout << i << "        " << atlas->labelnames[i] << "   " << atlas->labelTypes[i].pureLabel << "   ";
    if(atlas->labelTypes[i].pureLabel) cout << endl;
    else cout <<  atlas->labelTypes[i].mixed[0] << " " <<  atlas->labelTypes[i].mixed[1] << endl;   
  }
  if(atlas->useTPM) {
    cout << "TPM filenames:" << endl;
    for(i = 0;i < 3;i++) {
      cout << atlas->TPMfilenames[i] << endl;
    }
  }
}

// counts the number of pure labels. returns 0 if labels are not in the correct order
// Background label is not counted

inline int countPureLabels(AtlasSpec* atlas)
{
  int i,pureLabels,pve;

  pureLabels = 0;
  pve = 0;

  for(i = 1; i < atlas->numberOfLabels;i++) {
    if(atlas->labelTypes[i].pureLabel) {
      pureLabels++;
    }
  }
  // Check that the labels are in correct order 
  for(i = 1; i < atlas->numberOfLabels;i++) {
    if(!atlas->labelTypes[i].pureLabel) {
      pve = 1;
    }
    if(pve) {
      if(atlas->labelTypes[i].pureLabel) {
        pureLabels = 0;
      }
    }
  }
  if(!pureLabels) {
    cout << "Warning: Label order is not correct. There are PVE labels after pure labels" << endl;
  }

  return(pureLabels); 
}

#endif

