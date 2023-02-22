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



#include "atlasspec.h"


bool readAtlasSpec(AtlasSpec* atlas,char* filename) 
{
  char str[256];
  char* pstr = str;
  int i,j;
  char c;
  int tmp[256];
  bool probLimits;
  bool useTPM;
  bool onlyPureLabels;
  
  ifstream ifile;
  ifile.open(filename);
  if(!ifile.is_open()) return(false);
  ifile >> pstr;
  if(!(strcmp(pstr,"p"))) {
    probLimits = true;
    useTPM = false;
    onlyPureLabels = false;
  }
  else { 
    if(!(strcmp(pstr,"t"))) {
      probLimits = true;
      useTPM = true;
      onlyPureLabels = true;
    }
    else { 
      if(!(strcmp(pstr,"r"))) {
        probLimits = true;
        useTPM = false;
        onlyPureLabels = true;
      }
      else {
        probLimits = false;   
        useTPM = false;
        onlyPureLabels = false;
      }
    }
  }
  atlas->useTPM = useTPM;
  atlas->onlyPureLabels = onlyPureLabels;
  ifile >> atlas->n;
  ifile >> atlas->numberOfLabels;
 
  if(!probLimits) {
    atlas->regionLowProb.clear();
    atlas->regionUpProb.clear();
  }
  else {
    if(atlas->n > 0) {
      atlas->regionLowProb.resize(atlas->n);
      atlas->regionUpProb.resize(atlas->n);
      for(i = 0;i < atlas->n;i++) {
        atlas->regionLowProb[i].resize(atlas->numberOfLabels);
        atlas->regionUpProb[i].resize(atlas->numberOfLabels);
      }
    }
    else { 
      atlas->regionLowProb.resize(1);
      atlas->regionUpProb.resize(1);
      atlas->regionLowProb[0].resize(atlas->numberOfLabels);
      atlas->regionUpProb[0].resize(atlas->numberOfLabels);
    }
  }   

  // reserving the memory 
  if(atlas->n > 0) {
    atlas->filenames.resize(atlas->n);
    atlas->regionnames.resize(atlas->n); 
    atlas->permittedLabels.resize(atlas->n); 
  }
 
  atlas->labelnames.resize(atlas->n);
  atlas->labelTypes.resize(atlas->numberOfLabels);
  atlas->mrfConstants.resize(atlas->numberOfLabels);
  atlas->TPMfilenames.resize(atlas->n); 

  /*for(i = 0;i < atlas->n;i++) {
    atlas->regionnames[i] = new char[MAX_NAME_LEN];
    atlas->filenames[i] = new char[MAX_NAME_LEN];
  } */
 
  for(i = 0;i < atlas->numberOfLabels;i++) {
    //atlas->labelnames[i] = new char[MAX_NAME_LEN];
    atlas->mrfConstants[i].resize(atlas->numberOfLabels);
    //atlas->TPMfilenames[i] =  new char[MAX_NAME_LEN];
   }

  // set the background label info
  atlas->labelnames[0] = std::string("background");
  atlas->labelTypes[0].pureLabel = true; // pure label;  

  if(probLimits && atlas->n == 0 ) {
    atlas->regionLowProb[0][0] = 0.0;
    atlas->regionUpProb[0][0] = 0.0;
    for(j = 0; j < (atlas->numberOfLabels - 1);j++) {
      ifile >> atlas->regionLowProb[0][j + 1];
      if(atlas->regionLowProb[0][j + 1] < 0.0) atlas->regionLowProb[0][j + 1] = 0.0;
      ifile >> atlas->regionUpProb[0][j + 1];
      if( atlas->regionUpProb[0][j + 1] > 1.0) atlas->regionUpProb[0][j + 1] = 1.0;
    }
  }
  // read the region info
  for(i = 0;i < atlas->n;i++) {
    ifile >> pstr;
    atlas->regionnames[i] = std::string(pstr);
    c = ifile.get();
    ifile >> pstr;
    atlas->filenames[i] = std::string(pstr); 
    if(probLimits) {
      atlas->regionLowProb[i][0] = 0.0;
      atlas->regionUpProb[i][0] = 0.0;
      for(j = 0; j < (atlas->numberOfLabels - 1);j++) {
	ifile >> atlas->regionLowProb[i][j + 1];
        if(atlas->regionLowProb[i][j + 1] < 0.0) atlas->regionLowProb[i][j + 1] = 0.0;
        ifile >> atlas->regionUpProb[i][j + 1];
        if( atlas->regionUpProb[i][j + 1] > 1.0) atlas->regionUpProb[i][j + 1] = 1.0;
      }
      constructLabelList(&(atlas->permittedLabels[i]),atlas->numberOfLabels - 1,NULL,0); 
    }
    else {
      c = ifile.get();
      j = 0;
      if(c == ' ') {     
        while(c == ' ') {
          ifile >> tmp[j];
          c = ifile.get();   
          j++;
        }
      }      
      constructLabelList(&(atlas->permittedLabels[i]),atlas->numberOfLabels - 1,tmp,j);
    }
  } 
  // read the label info
  for(i = 1;i< atlas->numberOfLabels;i++) {
    ifile >> pstr;
    atlas->labelnames[i] = std::string(pstr);
    ifile >> atlas->labelTypes[i].pureLabel;
    if(ifile.get() != '\n') {
      ifile >> atlas->labelTypes[i].mixed[0];
      ifile >> atlas->labelTypes[i].mixed[1];
    }
    if(useTPM && (atlas->labelTypes[i].pureLabel)) {
      ifile >> pstr;
      atlas->TPMfilenames[i - 1] = std::string(pstr);
    }
    else {
      atlas->TPMfilenames[i - 1] = std::string("");
    }
  }
  
  // read the mrf matrix
  for(i = 0;i < atlas->numberOfLabels;i++) {
    for(j = 0;j < atlas->numberOfLabels;j++) {
      ifile >> atlas->mrfConstants[i][j];
    }
  }
  //  cout << "TPM1 " << atlas->TPMfilenames[0] << endl; 
  ifile.close();
  return true;
}

// This function works only for atlases of the type 'p', 't', or 'r'.
// not much checking is done, so use with care.

int writeAtlasSpec(AtlasSpec* atlas,char* filename, char atlasType, bool overwrite)
{
  
  int i,j;
  char c;
  
  ofstream ofile;
  ifstream ifile;

 // first few checks :

  switch(atlasType)
  {
  case 'p': 
   
    if(atlas->useTPM != false) return(1);
    if(atlas->onlyPureLabels != false) return(1);
    break;
  case 't': 
   
    if(atlas->useTPM != true) return(1);
    if(atlas->onlyPureLabels != true) return(1);
    break;
  case 'r': 
    
    if(atlas->useTPM != false) return(1);
    if(atlas->onlyPureLabels != true) return(1);
    break;
  default: return(2); 
  }
  if(!overwrite) { 
    // do not overwrite files
    ifstream ifile(filename);
    if(ifile.is_open()) {
      ifile.close();
      return(3);
    }
  } 

  ofile.open(filename);
  if(!ofile.is_open()) return(4);

  ofile << atlasType;
  ofile << ' ';
 
  ofile << atlas->n << ' ' << atlas->numberOfLabels << endl;
 
  // first n == 0 case
  if(atlas->n == 0 ) {
    for(j = 0; j < (atlas->numberOfLabels - 1);j++) {
      ofile << atlas->regionLowProb[0][j + 1] << ' ' << atlas->regionUpProb[0][j + 1] << endl;
    }
  }
  // write the region info when n > 0
  for(i = 0;i < atlas->n;i++) {
    ofile << atlas->regionnames[i];
    ofile << ' ';
    ofile << atlas->filenames[i] << endl; 
    for(j = 0; j < (atlas->numberOfLabels - 1);j++) {
      ofile << atlas->regionLowProb[i][j + 1] << ' ' << atlas->regionUpProb[i][j + 1] << endl;
    }
  } 
  // write the label info
  for(i = 1;i< atlas->numberOfLabels;i++) {
    ofile << atlas->labelnames[i];
    ofile << ' ';
    ofile << atlas->labelTypes[i].pureLabel;
    if(atlas->labelTypes[i].pureLabel) {
      ofile << " 0 0" << endl;
    }
    else {   
      ofile << ' ' << atlas->labelTypes[i].mixed[0] << ' ' << atlas->labelTypes[i].mixed[1] << endl;
    }
    if((atlas->useTPM) && (atlas->labelTypes[i].pureLabel)) {
      ofile << atlas->TPMfilenames[i - 1] << endl;;
    }
  }
  
  // write the mrf matrix
  for(i = 0;i < atlas->numberOfLabels;i++) {
    for(j = 0;j < (atlas->numberOfLabels - 1);j++) {
      ofile << atlas->mrfConstants[i][j] << ' ';
    }
    ofile << atlas->mrfConstants[i][atlas->numberOfLabels - 1] << endl;
  }
  //  cout << "TPM1 " << atlas->TPMfilenames[0] << endl; 
  ofile.close();
  return(0);
}


// void freeAtlas(AtlasSpec* atlas) 
// {
// }

// reads the images in the atlas and ensures that each
// voxelsums  are equal to the unity
// returns 0 if ok

int readAtlasImages(AtlasSpec* atlas,AnalyzeImage** atlasImages)
{
  int i,j,imgsize;
  int intstatus;
  int dimx,dimy,dimz;
  float probsum;
  AnalyzeImage psum;
  bool filter = false;

  // reserve memory
 // atlasImages = new AnalyzeImage*[atlas->n];

  // read images
  for(i = 0;i < atlas->n;i++) {
   // atlasImages[i] = new AnalyzeImage;
    intstatus = readImage(atlas->filenames[i].c_str(),atlasImages[i]);
    if(intstatus != 0) return(intstatus + i*100);
  }
  // check that the dimensions of the atlas match
  dimx = atlasImages[0]->header.x_dim;
  dimy = atlasImages[0]->header.y_dim; 
  dimz = atlasImages[0]->header.z_dim;
  for(i = 1;i < atlas->n;i++) {
    if(dimx != atlasImages[i]->header.x_dim) return(5 + i*10);
    if(dimy != atlasImages[i]->header.y_dim) return(5 + i*10);
    if(dimz != atlasImages[i]->header.z_dim) return(5 + i*10);
  }
  // make sure that every voxel sums to the unity in prob atlas
  imgsize = dimx*dimy*dimz;
  copyImage(atlasImages[0],&psum);

  for(i = 0;i < imgsize;i++) {
    probsum = 0.0;
    for(j = 0;j < atlas->n;j++) {
      probsum = probsum + atlasImages[j]->data[i];
    } 
    psum.data[i] = probsum;
    if(fabs(probsum) > 0.0001) {
      for(j = 0;j < atlas->n;j++) {
        atlasImages[j]->data[i] =  atlasImages[j]->data[i]/probsum;
      }
    }       
  }
  if(filter) {
   for(i = 0;i < atlas->n;i++) {
      averageFilterZeros(atlasImages[i],&psum,0.0001,1);
    }
    for(i = 0;i < imgsize;i++) {
      probsum = 0.0;
      for(j = 0;j < atlas->n;j++) {
        probsum = probsum + atlasImages[j]->data[i];
      } 
   
      if(fabs(probsum) > 0.0001) {
        for(j = 0;j < atlas->n;j++) {
          atlasImages[j]->data[i] =  atlasImages[j]->data[i]/probsum;
        }
      }       
    } 
  }   
  freeImage(&psum);
  return(0);
}

int readAtlasImages(AtlasSpec* atlas,std::vector<AnalyzeImage> & atlasImages)
{
  int i,j,imgsize;
  int intstatus;
  int dimx,dimy,dimz;
  float probsum;
  AnalyzeImage psum;
  bool filter = false;

  // reserve memory
 // atlasImages = new AnalyzeImage*[atlas->n];

  // read images
  for(i = 0;i < atlas->n;i++) {
   // atlasImages[i] = new AnalyzeImage;
    intstatus = readImage(atlas->filenames[i].c_str(),&atlasImages[i]);
    if(intstatus != 0) return(intstatus + i*100);
  }
  // check that the dimensions of the atlas match
  dimx = atlasImages[0].header.x_dim;
  dimy = atlasImages[0].header.y_dim; 
  dimz = atlasImages[0].header.z_dim;
  for(i = 1;i < atlas->n;i++) {
    if(dimx != atlasImages[i].header.x_dim) return(5 + i*10);
    if(dimy != atlasImages[i].header.y_dim) return(5 + i*10);
    if(dimz != atlasImages[i].header.z_dim) return(5 + i*10);
  }
  // make sure that every voxel sums to the unity in prob atlas
  imgsize = dimx*dimy*dimz;
  copyImage(&atlasImages[0],&psum);

  for(i = 0;i < imgsize;i++) {
    probsum = 0.0;
    for(j = 0;j < atlas->n;j++) {
      probsum = probsum + atlasImages[j].data[i];
    } 
    psum.data[i] = probsum;
    if(fabs(probsum) > 0.0001) {
      for(j = 0;j < atlas->n;j++) {
        atlasImages[j].data[i] =  atlasImages[j].data[i]/probsum;
      }
    }       
  }
  if(filter) {
   for(i = 0;i < atlas->n;i++) {
      averageFilterZeros(&atlasImages[i],&psum,0.0001,1);
    }
    for(i = 0;i < imgsize;i++) {
      probsum = 0.0;
      for(j = 0;j < atlas->n;j++) {
        probsum = probsum + atlasImages[j].data[i];
      } 
   
      if(fabs(probsum) > 0.0001) {
        for(j = 0;j < atlas->n;j++) {
          atlasImages[j].data[i] =  atlasImages[j].data[i]/probsum;
        }
      }       
    } 
  }   
  freeImage(&psum);
  return(0);
}

// releases memory allocated by readAtlasImages

void freeAtlasImages(AtlasSpec* atlas,AnalyzeImage** atlasImages)
{
  int i;

  for(i = 0;i < atlas->n;i++) {
    delete[] atlasImages[i]->data;
  }
}

// takes a mask and masks the atlas (makes atlas values out of mask 0)
// in addition tries to ensure that each voxel within the mask is represented in the atlas
// mask should have value 0 for background and 1 for foreground (i.e. brain) 
// atlas images should be normalized to sum to the unity.

 
bool maskAtlas(AtlasSpec* atlas, std::vector<AnalyzeImage> & atlasImages,AnalyzeImage* mask)
{
    AnalyzeImage probmask;
    int i,j,k,l,ci,cj,ck,imgsize;
    
    
   // check that atlas images correspond to the mask image
   if((atlasImages[0].header.x_dim != mask->header.x_dim) ||(atlasImages[0].header.y_dim != mask->header.y_dim) 
       || (atlasImages[0].header.z_dim != mask->header.z_dim)) {
     if((atlasImages[0].header.x_dim < mask->header.x_dim) ||(atlasImages[0].header.y_dim < mask->header.y_dim) 
	|| (atlasImages[0].header.z_dim <  mask->header.z_dim)) {
       cout << "ERROR: Atlas dimensions do not match to the maskfile dimensions" << endl;
       return(false);
     }
     else {
       cout << "Warning:  Atlas dimensions do not match to the maskfile dimensions -> continuing but check the results" << endl;
     }
   }
   newImage(&probmask,&atlasImages[0]);
   cout << "probmask created ok" << endl;
   imgsize = (atlasImages[0].header.x_dim)*(atlasImages[0].header.y_dim)*(atlasImages[0].header.z_dim); 
   for(i = 0;i < imgsize;i++) {
     probmask.data[i] = 0;
     for(j = 0;j < atlas->n;j++) {
       probmask.data[i] = probmask.data[i] + atlasImages[j].data[i];
     }
   }
  

   cout << "probmask generated" << endl;
   for(i = 0;i < mask->header.x_dim;i++ ) {
     for(j = 0;j < mask->header.y_dim;j++ ) {
       for(k = 0;k < mask->header.z_dim;k++ ) {
         // first check if the atlas spatially cover the mask
         if((getVoxelValue(mask,i,j,k) > 0.5) && (getVoxelValue(&probmask,i,j,k) < 0.5) ) { 
	   //  cout << "*";
           if(findClosestNonZero(&probmask,i,j,k,&ci,&cj,&ck,0.5) == false) {
	     cout << "here" << i  << " " << j << " " << k<< endl;
             return(false);
           }
           for(l = 0;l < atlas->n;l++) {
              putVoxelValue(&atlasImages[l],i,j,k,getVoxelValue(&atlasImages[l],ci,cj,ck));
           }
         }
         
         // then remove those parts of the atlas that are not necessary
         if(getVoxelValue(mask,i,j,k) < 0.5) {
           for(l = 0;l < atlas->n;l++) {
              putVoxelValue(&atlasImages[l],i,j,k,0.0);
           }
         }
       }
     }
   }
   // cout << "everything ok" << endl;
   freeImage(&probmask);
   return(true);
}

bool maskAtlas(AtlasSpec* atlas, AnalyzeImage** atlasImages, AnalyzeImage* mask)
{
	AnalyzeImage probmask;
	int i,j,k,l,ci,cj,ck,imgsize;


	// check that atlas images correspond to the mask image
	if((atlasImages[0]->header.x_dim != mask->header.x_dim) ||(atlasImages[0]->header.y_dim != mask->header.y_dim) 
		|| (atlasImages[0]->header.z_dim != mask->header.z_dim)) {
			if((atlasImages[0]->header.x_dim < mask->header.x_dim) ||(atlasImages[0]->header.y_dim < mask->header.y_dim) 
				|| (atlasImages[0]->header.z_dim <  mask->header.z_dim)) {
					cout << "ERROR: Atlas dimensions do not match to the maskfile dimensions" << endl;
					return(false);
			}
			else {
				cout << "Warning:  Atlas dimensions do not match to the maskfile dimensions -> continuing but check the results" << endl;
			}
	}
	newImage(&probmask,atlasImages[0]);
	cout << "probmask created ok" << endl;
	imgsize = (atlasImages[0]->header.x_dim)*(atlasImages[0]->header.y_dim)*(atlasImages[0]->header.z_dim); 
	for(i = 0;i < imgsize;i++) {
		probmask.data[i] = 0;
		for(j = 0;j < atlas->n;j++) {
			probmask.data[i] = probmask.data[i] + atlasImages[j]->data[i];
		}
	}


	cout << "probmask generated" << endl;
	for(i = 0;i < mask->header.x_dim;i++ ) {
		for(j = 0;j < mask->header.y_dim;j++ ) {
			for(k = 0;k < mask->header.z_dim;k++ ) {
				// first check if the atlas spatially cover the mask
				if((getVoxelValue(mask,i,j,k) > 0.5) && (getVoxelValue(&probmask,i,j,k) < 0.5) ) { 
					//  cout << "*";
					if(findClosestNonZero(&probmask,i,j,k,&ci,&cj,&ck,0.5) == false) {
						cout << "here: " << i  << " " << j << " " << k<< endl;
						return(false);
					}
					for(l = 0;l < atlas->n;l++) {
						putVoxelValue(atlasImages[l],i,j,k,getVoxelValue(atlasImages[l],ci,cj,ck));
					}
				}

				// then remove those parts of the atlas that are not necessary
				if(getVoxelValue(mask,i,j,k) < 0.5) {
					for(l = 0;l < atlas->n;l++) {
						putVoxelValue(atlasImages[l],i,j,k,0.0);
					}
				}
			}
		}
	}
	// cout << "everything ok" << endl;
	freeImage(&probmask);
	return(true);
}
/*
int readTPMimages(AtlasSpec* atlas,AnalyzeImage** TPMImages, AnalyzeImage* mask, int pureLabels)
{

  int i,j,k,l,ci,cj,ck,imgsize;
  int intstatus;
  int dimx,dimy,dimz;
  float probsum;
  AnalyzeImage psum;
  AnalyzeImage probmask;
   

  // reserve memory
 // atlasImages = new AnalyzeImage*[atlas->n];
  cout << "Entering readTPMimages" << endl;
  // read images
  for(i = 0;i < pureLabels;i++) {
   // atlasImages[i] = new AnalyzeImage;
    cout << i << endl;
    intstatus = readImage(atlas->TPMfilenames[i].c_str(),TPMImages[i]);
    if(intstatus != 0) return(intstatus + i*100);
  }
  cout << "TPM imagess read" << endl;
  // check that the dimensions of the atlas match
  dimx = TPMImages[0]->header.x_dim;
  dimy = TPMImages[0]->header.y_dim; 
  dimz = TPMImages[0]->header.z_dim;
  for(i = 1;i < pureLabels;i++) {
    if(dimx != TPMImages[i]->header.x_dim) return(5 + i*10);
    if(dimy != TPMImages[i]->header.y_dim) return(5 + i*10);
    if(dimz != TPMImages[i]->header.z_dim) return(5 + i*10);
  }
  // make sure that every voxel sums to the unity in prob atlas
  imgsize = dimx*dimy*dimz;
  copyImage(TPMImages[0],&psum);

  for(i = 0;i < imgsize;i++) {
    probsum = 0.0;
    for(j = 0;j < pureLabels;j++) {
      probsum = probsum + TPMImages[j]->data[i];
    } 
    psum.data[i] = probsum;
    if(fabs(probsum) > 0.0001) {
      for(j = 0;j < pureLabels;j++) {
        TPMImages[j]->data[i] =  TPMImages[j]->data[i]/probsum;
      }
    }       
  }
 
  freeImage(&psum);

   if((TPMImages[0]->header.x_dim != mask->header.x_dim) ||(TPMImages[0]->header.y_dim != mask->header.y_dim) 
       || (TPMImages[0]->header.z_dim != mask->header.z_dim)) {
     if((TPMImages[0]->header.x_dim < mask->header.x_dim) ||(TPMImages[0]->header.y_dim < mask->header.y_dim) 
	|| (TPMImages[0]->header.z_dim <  mask->header.z_dim)) {
       cout << "ERROR: TPM dimensions do not match to the maskfile dimensions" << endl;
       return(28);
     }
     else {
       cout << "Warning:  TPM dimensions do not match to the maskfile dimensions -> continuing but check the results" << endl;
     }
   }
   newImage(&probmask,TPMImages[0]);
   
 
   for(i = 0;i < imgsize;i++) {
     probmask.data[i] = 0;
     for(j = 0;j < pureLabels;j++) {
       probmask.data[i] = probmask.data[i] + TPMImages[j]->data[i];
     }
   }
  

   cout << "TPM generated" << endl;
   for(i = 0;i < mask->header.x_dim;i++ ) {
     for(j = 0;j < mask->header.y_dim;j++ ) {
       for(k = 0;k < mask->header.z_dim;k++ ) {
         // first check if the atlas spatially cover the mask
         if((getVoxelValue(mask,i,j,k) > 0.5) && (getVoxelValue(&probmask,i,j,k) < 0.5) ) { 
	   //  cout << "*";
           if(findClosestNonZero(&probmask,i,j,k,&ci,&cj,&ck,0.5) == false) {
             for(l = 0;l < pureLabels;l++) {
               putVoxelValue(TPMImages[l],i,j,k,0.333333333);
             }
           }	     
	   else {
             for(l = 0;l < pureLabels;l++) {
                putVoxelValue(TPMImages[l],i,j,k,getVoxelValue(TPMImages[l],ci,cj,ck));
             }
           }
         }
         
         // then remove those parts of the atlas that are not necessary
         if(getVoxelValue(mask,i,j,k) < 0.5) {
           for(l = 0;l < pureLabels;l++) {
              putVoxelValue(TPMImages[l],i,j,k,0.0);
           }
         }
       }
     }
   }
   // cout << "everything ok" << endl;
   freeImage(&probmask);


  return(0);
}

void freeTPMimages(AtlasSpec* atlas,AnalyzeImage** TPMImages, int pureLabels)
{
  int i;

  for(i = 0;i < pureLabels;i++) {
    delete[] TPMImages[i]->data;
  }
}
*/