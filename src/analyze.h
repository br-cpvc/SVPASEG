// ******************************************************
// Functions to manipulate Analyze 7.5 images 
// (C) Jussi Tohka, 2005 - 2009 
// Laboratory of NeuroImaging, Department of Neurology, UCLA Medical School, CA, USA
// Department of Signal Processing, Tampere University of Technology, Finland  
// Version 2.1:
// Modified to include support for NIFTI format 14 Oct 2010 
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
// Change log and notes
// 2.0 -> 2.1: Updates to handle Nifti-files: This is based on the API provided
// in http://afni.nimh.nih.gov/pub/dist/src/nifti/niftilib
// Files  nifti1.h, nifti1_io.h, nifti1_io.c, znzlib.c znzlib.h are required.
// Note that nifti_io.h is modified in order to avoid the need to translate znzlib 
// into a library. Nifti support means that nifti-files can be read and written. 
// However, reading analyze in and writing nifti is not supported. 
// gzipped files are not supported, but the support is likely to be easy to add.
// The modification is a hack, and could be done better. The idea is to add the
// Nifti-header information to the AnalyzeImage type. 


#ifndef ANALYZE_H
#define ANALYZE_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <values.h>
#include <math.h>
#include "nifti1.h"
#include "nifti1_io.h"

using namespace std;


#define DT_NONE                 0 
#define DT_UNKNOWN              0      /*Unknown data type*/ 
#define DT_BINARY               1      /*Binary (1 bit per voxel)*/ 
#define DT_UNSIGNED_CHAR				2      /*Unsigned character (8 bits per voxel)*/ 
#define DT_SIGNED_SHORT					4      /*Signed short (16 bits per voxel)*/ 
#define DT_SIGNED_INT           8      /*Signed integer (32 bits per voxel)*/ 
#define DT_FLOAT                16     /*Floating point (32 bits per voxel)*/ 
#define DT_COMPLEX              32     /*Complex (64 bits per voxel; 2 floating point numbers) */
#define DT_DOUBLE               64     /*Double precision (64 bits per voxel)*/ 
#define DT_RGB                  128    /* */
#define DT_ALL                  255    /* */ 

#define F_UNKNOWN              0   /* Image format is unknown 14 OCT 10 */
#define F_ANA                  1   /* Image format is analyze 14 OCT 10 */
#define F_NII                  2   /* Image format is Nifti   14 OCT 10 */ 

#ifdef   MAX
#undef   MAX
#endif
#define  MAX( x, y )  ( ((x) >= (y)) ? (x) : (y) )

#ifdef   MIN
#undef   MIN
#endif
#define  MIN( x, y )  ( ((x) <= (y)) ? (x) : (y) )


#define ROUND(x) (ceil((x) - 0.5))

#define ARE_EQUALF(x,y) (fabs((x) - (y)) < 0.0000001)   

#define FLOAT_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

// some definitions borrowed from fslio.h
// these are for convertBufferToScaledFloat to work 
// Added 14 OCT 10

typedef unsigned char   THIS_UINT8; 
typedef char            THIS_INT8;
typedef unsigned short  THIS_UINT16;
typedef short           THIS_INT16;
typedef unsigned int    THIS_UINT32;
typedef int             THIS_INT32;
typedef unsigned long   THIS_UINT64;
typedef long            THIS_INT64;
typedef float           THIS_FLOAT32;
typedef double          THIS_FLOAT64;

 

struct AnalyzeHeader {
	AnalyzeHeader()
	{
		memset((char *)this,0,sizeof(*this));
		regular = 'r';
		sizeof_hdr=sizeof(*this);
	}
	int		sizeof_hdr;	/* For ANALYZE compatibility only */
	char		pad1[28];
	int		extents;	/* For ANALYZE compatibility only */
	char		pad2[2];
	char		regular;	/* For ANALYZE compatibility only */
	char		pad3;
	short int	dims;		/* For ANALYZE compatibility only */
	short int	x_dim;		/* AIR */
	short int	y_dim;		/* AIR */
	short int	z_dim;		/* AIR */
	short int	t_dim;		/* For ANALYZE compatibility only */
	char		pad4[20];
	short int	datatype;	/* For ANALYZE compatibility only */
	short int	bits;		/* AIR */
	char		pad5[6];
	float		x_size;		/* AIR */
	float		y_size;		/* AIR */
	float		z_size;		/* AIR */
	float		x_orig;		/* dws */
	float		y_orig;		/* dws */
	float		z_orig;		/* dws */
	char		pad6[36];	// was 48, now 12 for origin codes.
	int		glmax;		/* AIR */
	int		glmin;		/* AIR */
	char		descrip[80];	/* AIR (non-essential) */
	char		pad7[120];
};

struct AnalyzeImage {
  AnalyzeHeader header;
  float* data; // Internal data presented using floats 
  int format;  // Image format ana, nii or unknown 14 OCT 10
  nifti_image nii_header; // Pointer to NIFTI header if the format is nifti 14 OCT 10
} ; 

struct AnalyzeLabelImage {
  AnalyzeHeader header;
  unsigned char* data; // If memory is valuable
  int format;  // Image format ana, nii or unknown 14 OCT 10
  nifti_image nii_header; // Pointer to NIFTI header if the format is nifti 14 OCT 10
} ;

void swapHeader(AnalyzeHeader &h);

inline void headername(char *dst, const char *src)
{
	strcpy(dst,src);
	char *p = strrchr(dst,'.');
	if (!p)
	{
		strcat(dst,".hdr");
	}
	else
	{
		strcpy(p,".hdr");
	}
}

inline void imagename(char *dst, const char *src)
{
	strcpy(dst,src);
	char *p = strrchr(dst,'.');
	if (!p)
	{
		strcat(dst,".img");
	}
	else
	{
		strcpy(p,".img");
	}
}

// gets the filename from full path name

inline void extractpath(char* filename, char* path, const char* fullpath)
{
  strcpy(filename,fullpath);
  strcpy(path,fullpath);
  char *p = strrchr(path,'/');
  if (!p) {
    strcpy(path,"");
  }
  else {
    strcpy(filename,p + 1);
    strcpy(p + 1,"");
  }       
} 

// Finds out the image format based on filename 
// 14 OCT 10

inline int imageFormat(const char *ifname) 
{
   const char *p = ifname + strlen(ifname) - 4; // .img or .hdr or .nii or none of these
        if ('.' != *p)
        {
	   return F_UNKNOWN;
        }
        else
        {
           if (!(strcmp(p + 1,"nii"))) return F_NII;
           else if (!(strcmp(p + 1,"hdr"))) return F_ANA;
           else if (!(strcmp(p + 1,"img"))) return F_ANA;
           else return F_UNKNOWN;
           
        }

}



inline bool readHeader(const char *ifname, AnalyzeHeader& header)
{
	char hfname[256];
	headername(hfname,ifname);
	ifstream ifile(hfname);
	if (!ifile) return false;
	ifile.read((char *)&header,sizeof(AnalyzeHeader));
	return true;
}

inline bool writeHeader(const char *ofname, AnalyzeHeader& header)
{
	char hfname[256];
	headername(hfname,ofname);
	ofstream ofile(hfname);
	if (!ofile) return false;
	ofile.write((char *)&header,sizeof(AnalyzeHeader));
	return true;
}


inline void byteswap(unsigned short &X)
{
	unsigned short a = (unsigned char)(X>>8);
	unsigned char b = (unsigned char)X;
	X = (int(b)<<8)|a;
}

inline void byteswap(short &X)
{
	unsigned char a = (unsigned char)(X>>8);
	unsigned char b = (unsigned char)X;
	X = (int(b)<<8)|a;
}

inline void swapx(unsigned int &X)
{
	unsigned short a = X>>16;
	byteswap(a);
	unsigned short b = (unsigned short)(X);
	byteswap(b);
	X = (int(b)<<16)|a;
}

inline void swapf(unsigned int *X)
{
	unsigned short a = *X>>16;
	byteswap(a);
	unsigned short b = (unsigned short)(*X);
	byteswap(b);
	*X = (int(b)<<16)|a;
}

inline void swapx(int &X)
{
	unsigned short a = X>>16;
	byteswap(a);
	unsigned short b = (unsigned short)(X);
	byteswap(b);
	X = (int(b)<<16)|a;
}

inline void swapHeader(AnalyzeHeader &h)
{
	swapx(h.sizeof_hdr); // int		sizeof_hdr;
	//char		pad1[28];
	swapx(h.extents);		//int		extents;	/* For ANALYZE compatibility only */
	//char		pad2[2];
	//char		regular;	/* For ANALYZE compatibility only */
	//char		pad3;
	byteswap(h.dims); //short int	dims;		/* For ANALYZE compatibility only */
	byteswap(h.x_dim);//short int	x_dim;		/* AIR */
	byteswap(h.y_dim);//short int	y_dim;		/* AIR */
	byteswap(h.z_dim);//short int	z_dim;		/* AIR */
	byteswap(h.t_dim);//short int	t_dim;		/* For ANALYZE compatibility only */
	//char		pad4[20];
	byteswap(h.datatype);//short int	datatype;	/* For ANALYZE compatibility only */
	byteswap(h.bits);//short int	bits;		/* AIR */
	//char		pad5[6];
	swapf((unsigned int *)&h.x_size);//float		x_size;		/* AIR */
	swapf((unsigned int *)&h.y_size);//float		y_size;		/* AIR */
	swapf((unsigned int *)&h.z_size);//float		z_size;		/* AIR */
	swapf((unsigned int *)&h.x_orig);//float		x_size;		/* AIR */
	swapf((unsigned int *)&h.y_orig);//float		y_size;		/* AIR */
	swapf((unsigned int *)&h.z_orig);//float		z_size;		/* AIR */
	//char		pad6[48];
	swapx(h.glmax);//int		glmax;		/* AIR */
	swapx(h.glmin);//int		glmin;		/* AIR */
	//char		descrip[80];	/* AIR (non-essential) */
	//char		pad7[120];
};
 
inline bool byteSwapNecessary(AnalyzeHeader* h) {
  if( h->sizeof_hdr == 348) return false;
  else return true;
};  

int copyNifti2Analyze(nifti_image* nim, AnalyzeImage* img);
// int copyAnalyze2Nifti(AnalyzeImage* img, nifti_image* nim);
int copyNifti2AnalyzeLabel(nifti_image* nim, AnalyzeLabelImage* img);

int readImage(const char* filename,AnalyzeImage* img, bool signedData = true);
int writeImage(const char* filename, AnalyzeImage* img, bool overwrite, bool signedData = true);
int readLabelImage(char* filename,AnalyzeLabelImage* img);
int writeLabelImage(char* filname,AnalyzeLabelImage* img,bool overwrite);

bool copyImage(AnalyzeImage* source,AnalyzeImage* target);
bool copyLabelImage(AnalyzeLabelImage* source,AnalyzeLabelImage* target);

bool newImage(AnalyzeImage* img, AnalyzeImage* alike);  
bool newLabelImage(AnalyzeLabelImage* img, AnalyzeImage* alike);

bool binaryMask(AnalyzeImage* img,AnalyzeLabelImage* mask,float tr); 
bool mask2Image(AnalyzeLabelImage* source, AnalyzeImage* target);


bool findClosestNonZero(AnalyzeImage*,int x,int y, int z,int* cx, int* cy, int* cz,float tr);
void averageFilterZeros(AnalyzeImage* img,AnalyzeImage* refImg,float tr,int times);

inline void thresholdImage(AnalyzeImage* img,float val) {
  int i,imgsize;
  imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim); 
  for(i = 0;i< imgsize;i++) {
    if(img->data[i] > val) img->data[i] = 1;
    else img->data[i] = 0;
  }
}   

inline void freeImage(AnalyzeImage* img) 
{
  delete[] img->data;
};

inline void freeLabelImage(AnalyzeLabelImage *img) 
{
  delete[] img->data;
};

inline float getVoxelValue(AnalyzeImage* img,int x, int y, int z) 
{
  return img->data[x + y*(img->header.x_dim) + z*(img->header.y_dim)*(img->header.x_dim)];
};

inline float getSafeVoxelValue(AnalyzeImage* img,int x, int y, int z)
{
  if((x < 0) || (y < 0) || (z < 0) || (x > (img->header.x_dim - 1)) || 
     (y > (img->header.y_dim - 1)) || (z > (img->header.z_dim - 1))) 
    return(0);
  else  return img->data[x + y*(img->header.x_dim) + z*(img->header.y_dim)*(img->header.x_dim)];
}

inline void putVoxelValue(AnalyzeImage* img,int x, int y, int z, float val) 
{
  img->data[x + y*(img->header.x_dim) + z*(img->header.y_dim)*(img->header.x_dim)] = val;
};


inline unsigned char getLabelValue(AnalyzeLabelImage* img,int x, int y, int z) 
{
  return img->data[x + y*(img->header.x_dim) + z*(img->header.y_dim)*(img->header.x_dim)];
};

// this is the same as the previous one exept returns 0 if we are ouside the image

inline unsigned char getSafeLabelValue(AnalyzeLabelImage* img,int x, int y, int z)
{
  if((x < 0) || (y < 0) || (z < 0) || (x > (img->header.x_dim - 1)) || 
     (y > (img->header.y_dim - 1)) || (z > (img->header.z_dim - 1))) 
    return(0);
  else  return img->data[x + y*(img->header.x_dim) + z*(img->header.y_dim)*(img->header.x_dim)];
}

inline void putLabelValue(AnalyzeLabelImage* img,int x, int y, int z, unsigned char val) 
{
  img->data[x + y*(img->header.x_dim) + z*(img->header.y_dim)*(img->header.x_dim)] = val;
};

// inline double gaussian(double x,double mu,double sigma) {
// return(exp(-(pow((x - mu),2))/(2 * sigma * sigma))/(sqrt(2 * PI) * sigma))) };

inline float imageMax(AnalyzeImage* img,AnalyzeLabelImage* mask) 
{   
  int i,j,k;
  float maximum;
  
  maximum = MINFLOAT;
  for(i = 0;i < img->header.x_dim;i++) {
    for(j = 0;j < img->header.y_dim;j++) {
       for(k = 0;k < img->header.z_dim;k++) {
         if(getLabelValue(mask,i,j,k) > 0) {
           if(getVoxelValue(img,i,j,k) > maximum) maximum = getVoxelValue(img,i,j,k);
         }
       } 
    }
  }
  return(maximum);
};

inline unsigned char labelMax(AnalyzeLabelImage* img) 
{   
  int i,j,k;
  unsigned char maximum;
  
  maximum = 0;
  for(i = 0;i < img->header.x_dim;i++) {
    for(j = 0;j < img->header.y_dim;j++) {
       for(k = 0;k < img->header.z_dim;k++) {
         if(getLabelValue(img,i,j,k) > maximum) maximum = getLabelValue(img,i,j,k);
         
       } 
    }
  }
  return(maximum);
};

  
inline float imageMin(AnalyzeImage* img,AnalyzeLabelImage* mask) 
{   
  int i,j,k;
  float minimum;
  
  minimum = MAXFLOAT;
  for(i = 0;i < img->header.x_dim;i++) {
    for(j = 0;j < img->header.y_dim;j++) {
       for(k = 0;k < img->header.z_dim;k++) {
         if(getLabelValue(mask,i,j,k) > 0) {
           if(getVoxelValue(img,i,j,k) < minimum) minimum = getVoxelValue(img,i,j,k);
         }
       } 
    }
  }
  return(minimum);
};

float imageKth(AnalyzeImage* img, AnalyzeLabelImage* mask, float percentage);


// 3D erosion with 3 x 3 x 3 structuring element

inline void erode3D(AnalyzeLabelImage* img) 
{
  int i,j,k,x,y,z;
  bool boundaryVoxel;
  AnalyzeLabelImage tmp;

  copyLabelImage(img,&tmp);
  for(i = 0;i < img->header.x_dim;i++) {
    for(j = 0;j < img->header.y_dim;j++) {
      for(k = 0;k < img->header.z_dim;k++) {
        if(getLabelValue(img,i,j,k) > 0.5) {
          boundaryVoxel = false;
          for(x = -1; x < 2; x++) {
            for(y = -1; y < 2; y++) {
              for(z = -1; z < 2; z++) {
                if(getSafeLabelValue(&tmp,i + x, j + y, k + z) < 0.5) {
		  boundaryVoxel = true;
                }
	      }
	    }
          }
	  if(boundaryVoxel) {
            putLabelValue(img,i,j,k,0);
          }
	}
      }
    }
  }
  freeLabelImage(&tmp);
}

void medianFilter3D(AnalyzeImage* img, AnalyzeLabelImage* mask, int windowX, int windowY, int windowZ);
void knnFilter(AnalyzeImage* img, AnalyzeLabelImage* mask, int windowLenX, int windowLenY, int windowLenZ, int knn);
  
inline char maxArg(float* pval,int n)
{
  float maximum;
  int i,index;
  
  maximum = pval[0];
  index = 0;
  for(i = 1;i < n;i++) {
    if(pval[i] > maximum) {
      index = i;
      maximum = pval[i];
    }
  }
  return( (char) index);
}


#endif

