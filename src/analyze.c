#include "analyze.h"
#include <math.h>
using namespace std;

// ******************************************************
// Functions to manipulate Analyze 7.5 images 
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




/***************************************************************
 * convertBufferToScaledFloat, modified from the function  
   convertBufferToScaledDouble from the public domain fslio library  
   
 ***************************************************************/
/*! \fn int  convertBufferToScaledFloat(float *outbuf, void *inbuf, long len, 
                                        float slope, float inter, int nifti_datatype )
    \brief allocate a 4D buffer, use 1 contiguous buffer for the data 

        Array is indexed as buf[0..th-1][0..zh-1][0..yh-1][0..xh-1].  
        <br>To access all elements as a vector, use buf[0][0][0][i] where
        i can range from 0 to th*zh*yh*xh - 1.

    \param outbuf pointer to array of float of size len
    \param inbuf void pointer to an array of len items of datatype nifti_datatype
    \param len number of elements in outbuf and inbuf
    \param slope slope term of scaling to be applied
    \param inter intercept term of scaling to be applied:  out = (in*slope)+inter
    \param nifti_datatype NIFTI datatype code for the datatype of the elements in inbuf
    \return error code: 0=OK -1=error
 */
int  convertBufferToScaledFloat(float *outbuf, void *inbuf, int len, float slope, float inter, int nifti_datatype ) 
{

        int i;


    /** fill the buffer */
    for (i=0; i<len; i++)
        switch(nifti_datatype) {
            case NIFTI_TYPE_UINT8:
                outbuf[i] = (float) ( *((THIS_UINT8 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT8:
                outbuf[i] = (float) ( *((THIS_INT8 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_UINT16:
                outbuf[i] = (float) ( *((THIS_UINT16 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT16:
                outbuf[i] = (float) ( *((THIS_INT16 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_UINT64:
                outbuf[i] = (float) ( *((THIS_UINT64 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT64:
                outbuf[i] = (float) ( *((THIS_INT64 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_UINT32:
                outbuf[i] = (float) ( *((THIS_UINT32 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_INT32:
                outbuf[i] = (float) ( *((THIS_INT32 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_FLOAT32:
                outbuf[i] = (float) ( *((THIS_FLOAT32 *)(inbuf)+i) * slope + inter);
                break;
            case NIFTI_TYPE_FLOAT64:
                outbuf[i] = (float) ( *((THIS_FLOAT64 *)(inbuf)+i) * slope + inter);
                break;

            case NIFTI_TYPE_FLOAT128:
            case NIFTI_TYPE_COMPLEX128:
            case NIFTI_TYPE_COMPLEX256:
            case NIFTI_TYPE_COMPLEX64:
            default:
	      //    fprintf(stderr, "\nWarning, cannot support %s yet.\n",nifti_datatype_string(nifti_datatype));
                return(-1);
        }

return(0);
}
/***************************************************************
 * convertBufferToChar, modified from the function  
   convertBufferToScaledDouble from the public domain fslio library  
   
 ***************************************************************/
/*! \fn int  convertBufferToChar(float *outbuf, void *inbuf, long len, 
                                        float slope, float inter, int nifti_datatype )
    \brief allocate a 4D buffer, use 1 contiguous buffer for the data 

        Array is indexed as buf[0..th-1][0..zh-1][0..yh-1][0..xh-1].  
        <br>To access all elements as a vector, use buf[0][0][0][i] where
        i can range from 0 to th*zh*yh*xh - 1.

    \param outbuf pointer to array of float of size len
    \param inbuf void pointer to an array of len items of datatype nifti_datatype
    \param len number of elements in outbuf and inbuf
  
    \return error code: 0=OK -1=error
 */
int  convertBufferToChar(unsigned char *outbuf, void *inbuf, int len) 
{

  int i; 
 /** fill the buffer */
  for (i=0; i<len; i++) {
    outbuf[i] = (unsigned char) ( *((THIS_UINT8 *)(inbuf)+i));
  }

return(0);
}


int copyNifti2Analyze(nifti_image* nim, AnalyzeImage* img)
{
  int code,imgsize;
  float slope,inter;

  img->header.sizeof_hdr = 348;
  img->header.dims = 3;
  img->header.x_dim = nim->nx;
  img->header.y_dim = nim->ny;
  img->header.z_dim = nim->nz;
  img->header.x_size = nim->dx;
  img->header.y_size = nim->dy;
  img->header.z_size = nim->dz;
  img->header.datatype = DT_FLOAT; // This is to match the internal data representation
  img->nii_header = *nim; // copy the Nifti header info
  imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim); 
  img->data = new float[imgsize]; 
  slope = nim->scl_slope;
  inter = nim->scl_inter;
  if(slope == 0) slope = 1;
  if(isnan(slope)) slope = 1;
  if(isnan(inter)) inter = 1; 

  code = convertBufferToScaledFloat(img->data, nim->data,imgsize,slope, inter, nim->datatype);
  nifti_image_free(nim);
  img->nii_header.scl_slope = 1;
  img->nii_header.scl_inter = 0;
  img->nii_header.datatype = NIFTI_TYPE_FLOAT32;
  img->nii_header.data = img->data;
  img->nii_header.nbyper = 4;   
  return(code);

}

int copyNifti2AnalyzeLabel(nifti_image* nim, AnalyzeLabelImage* img)
{
  int code,imgsize;
 
  img->header.sizeof_hdr = 348;
  img->header.dims = 3;
  img->header.x_dim = nim->nx;
  img->header.y_dim = nim->ny;
  img->header.z_dim = nim->nz;
  img->header.x_size = nim->dx;
  img->header.y_size = nim->dy;
  img->header.z_size = nim->dz;
  img->header.datatype = DT_UNSIGNED_CHAR; 
  img->nii_header = *nim; // copy the Nifti header info
  imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim); 
  img->data = new unsigned char[imgsize]; 
  //slope = nim->scl_slope;
  //inter = nim->scl_inter;
  code = convertBufferToChar(img->data, nim->data,imgsize);
  nifti_image_free(nim);
  img->nii_header.scl_slope = 1;
  img->nii_header.scl_inter = 0;
  img->nii_header.datatype =  NIFTI_TYPE_UINT8;
  img->nii_header.data = img->data;
  img->nii_header.nbyper = 1;   
  return(code);

}

// Changes only the data and other relevent parts of the existing
// NIFTi file  
// int copyAnalyze2Nifti(AnalyzeImage* img, nifti_image* nim)
// { 
//  return(0);
// }


// Reads an Analyze image 
// Supports datatypes UNSIGNED_CHAR, SIGNED SHORT, SIGNED INT AND SIGNED FLOAT
// Returns:  0 if succesful
//           1 if the header could not be opened
//           2 if the memory allocation did not work as usual
//           3 if the image file could not be opened
//           4 if the datatype was not supported
//           5 if nii file could not be opened
//           6 if nii file could not be copieed to internal format

int readImage(const char* filename,AnalyzeImage* img, bool signedData)
{
 
  char imgfname[256];
  unsigned char* tmpchar;
  short* tmpshort;
  unsigned short* tmpushort;
  int* tmpint;  
  int i,imgsize;
  bool bs;
  int format;    // 14 OCT 10
  nifti_image * nim=NULL; // 14 OCT 10 

  format = imageFormat(filename);
  if (format == F_NII) {
    img->format = F_NII;
    nim = nifti_image_read(filename, 1);
    if( !nim ) return(5);
    if(copyNifti2Analyze(nim,img)) return(6);            
  }
  else { // reading Analyze file 
    img->format = F_ANA;
    if(readHeader(filename,img->header) == false) return(1);
    bs = byteSwapNecessary(&img->header);
 
    if(bs) swapHeader(img->header);  
 
    imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim); 
    // img->data = (float *) malloc(imgsize*sizeof(float));  
    img->data = new float[imgsize];
    if(img->data == NULL) return(2);
    imagename(imgfname,filename);
    ifstream ifile(imgfname,ios::binary);
    if (!ifile) return(3);
  
    switch(img->header.datatype) {
	case DT_UNSIGNED_CHAR:
	     tmpchar = new unsigned char [imgsize];
             ifile.read((char*) tmpchar,imgsize);
             for(i = 0;i < imgsize;i++) {
               img->data[i] =  (float) tmpchar[i];
             }
             delete[] tmpchar;
             ifile.close();
             break;
        case DT_SIGNED_SHORT: 
	     if(signedData) {
	       tmpshort = new short [imgsize];  
               ifile.read((char*)tmpshort,imgsize*sizeof(short)); // corrected 27 Sep 2010
               for(i = 0;i < imgsize;i++) {
                 if (bs) byteswap(tmpshort[i]);
                 img->data[i] = (float) tmpshort[i];
               }            
               delete[] tmpshort;
               ifile.close();
             }
             else {
               tmpushort = new unsigned short [imgsize];  
               ifile.read((char*)tmpushort,imgsize*sizeof(short)); // corrected 27 Sep 2010
               for(i = 0;i < imgsize;i++) {
                 if (bs) byteswap(tmpushort[i]);
                 img->data[i] = (float) tmpushort[i];
               }            
               delete[] tmpushort;
               ifile.close();
             }
             break;
        case DT_SIGNED_INT: 
	     tmpint = new int [imgsize];  
	     ifile.read((char*) tmpint,imgsize*sizeof(int));
             for(i = 0;i < imgsize;i++) {
               if (bs) swapx(tmpint[i]);
               img->data[i] = (float) tmpint[i];
             }
             ifile.close();
             delete[] tmpint;
             break;  
        case DT_FLOAT:             
	     ifile.read((char*)img->data,imgsize*sizeof(float));
             for(i = 0;i < imgsize;i++) {
               if (bs) swapf((unsigned int *)&img->data[i]);
               
             }
             ifile.close();
             break;
        default: return(4); break;
    }
  }
  return(0);

}
   

// Writes an Analyze Image 
// returns 0 if ok
// returns 1 if the header could not be opened for writing
// returns 2 if the image file could not be opened for writing 
// returns 3 if the datatype was not supported
// returns 4 if the file exists and overwriting is not permitted
// returns 5 if the nifti file does not seem valid  

int writeImage(char* filename,AnalyzeImage* img,bool overwrite, bool signedData) 
{
  char imgfname[256];
 
  int imgsize;
  int i; 
  int format;    // 14 OCT 10
 
 
  imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim); 
   format = imageFormat(filename);
   
  if ((format == F_NII) && (img->format == F_NII)) {
    if(! overwrite) {
      ifstream ifile(filename,ios::binary);
      if(ifile.is_open()) {
        ifile.close();
        return(4);
      }
    }
 
    if(!(nifti_nim_is_valid( &(img->nii_header), 1))) return(5);
    strcpy(img->nii_header.iname,filename);
    strcpy(img->nii_header.fname,filename);
    //   if( nifti_set_filenames(&(img->nii_header), filename, 1, 1) ) return(5);
    img->nii_header.data = img->data;    
    nifti_image_write( &(img->nii_header) );             
    
  }
  else { 
    imagename(imgfname,filename);
    if(! overwrite) {
      ifstream ifile(imgfname,ios::binary);
      if(ifile.is_open()) {
        ifile.close();
        return(4);
      }
    }
    if(writeHeader(filename,img->header) == false) return(1);
    ofstream ofile(imgfname,ios::binary);
    //  else ofstream ofile(imgfname,ios::binary|ios::noreplace); 
    if (!ofile) return(2);
    switch(img->header.datatype) {
	case DT_UNSIGNED_CHAR:
	    for(i = 0;i<imgsize;i++) {
              unsigned char tmp;
	      img->data[i] = ROUND(img->data[i]);
              img->data[i] = MIN(img->data[i],255);
              img->data[i] = MAX(img->data[i],0);
              tmp = (unsigned char) img->data[i];
              ofile.write((char *)&tmp,sizeof(char)); 
            } 
            break;
  	case DT_SIGNED_SHORT:
	    for(i = 0;i<imgsize;i++) {
              if(signedData) {
                short tmp;
	        img->data[i] = ROUND(img->data[i]);
                img->data[i] = MIN(img->data[i],32768);
                img->data[i] = MAX(img->data[i],-32767);
                tmp = (short) img->data[i];
                ofile.write((char *)&tmp,sizeof(short));
              }
              else {
                unsigned short tmp;
                img->data[i] = ROUND(img->data[i]);
                img->data[i] = MIN(img->data[i],65535);
                img->data[i] = MAX(img->data[i],0);
                tmp = (unsigned short) img->data[i];
                ofile.write((char *)&tmp,sizeof(short));
	      }
            } 
            break;
        case DT_SIGNED_INT:
	    for(i = 0;i<imgsize;i++) {
              int tmp;
	      img->data[i] = ROUND(img->data[i]);
              img->data[i] = MIN(img->data[i],2147000000);
              img->data[i] = MAX(img->data[i],-2147000000);
              tmp = (int) img->data[i];
              ofile.write((char *)&tmp,sizeof(int)); 
            }   
            break;
        case DT_FLOAT:
	    for(i = 0;i<imgsize;i++) {
              float tmp;
              tmp =  img->data[i];
              ofile.write((char *)&tmp,sizeof(float)); 
            }
            break; 
       default: return(3); break;
    }
    ofile.close();
  }
     
  return(0);
}

int readLabelImage(char* filename,AnalyzeLabelImage* img)
{
  char imgfname[256];
  unsigned char* tmpchar;
  int i,imgsize;
  bool bs;
  int format;    // 14 OCT 10
  nifti_image * nim=NULL; // 14 OCT 10 

  format = imageFormat(filename); // 14 OCT 10 
  if (format == F_NII) {
    img->format = F_NII;
    nim = nifti_image_read(filename, 1);
    if( !nim ) return(5);
    if(!(nim->datatype ==  NIFTI_TYPE_UINT8)) return(6);
    if(copyNifti2AnalyzeLabel(nim,img)) return(7);            
  }
  else { // reading Analyze file 
    img->format = F_ANA;
    if(readHeader(filename,img->header) == false) return(1);
    bs = byteSwapNecessary(&img->header);
    if(bs) swapHeader(img->header);  

    imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim); 
    // img->data = (float *) malloc(imgsize*sizeof(float));  
    img->data = new unsigned char[imgsize];
    if(img->data == NULL) return(2);
    imagename(imgfname,filename);
    ifstream ifile(imgfname,ios::binary);
    if (!ifile) return(3);
    if(img->header.datatype != DT_UNSIGNED_CHAR) return(4);
    else {
      ifile.read((char*) img->data,imgsize);
      ifile.close();
    }  
  }
  return(0); 
}

int writeLabelImage(char* filename,AnalyzeLabelImage* img,bool overwrite)
{
  char imgfname[256];
 
  int imgsize;
  int i; 
  int format; // 14 OCT 10

  imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim);
  format = imageFormat(filename);
  if ((format == F_NII) && (img->format == F_NII)) {
    if(! overwrite) {
      ifstream ifile(filename,ios::binary);
      if(ifile.is_open()) {
        ifile.close();
        return(4);
      }
    }
    if(!(nifti_nim_is_valid( &(img->nii_header), 1))) return(5);
	int filenamelength = strlen(filename)+2;
	img->nii_header.iname = (char*)malloc(filenamelength*sizeof(char*));
	img->nii_header.fname = (char*)malloc(filenamelength*sizeof(char*));
    strcpy(img->nii_header.iname,filename);
    strcpy(img->nii_header.fname,filename);
       //  if( nifti_set_filenames(&(img->nii_header), filename, 1, 1) ) return(5);
    img->nii_header.data = img->data; 
    nifti_image_write( &(img->nii_header) );             
  }
  else { 
    imagename(imgfname,filename); 
    if(! overwrite) {
      ifstream ifile(imgfname,ios::binary);
      if(ifile.is_open()) {
        ifile.close();
        return(4);
      }
    }
    if(writeHeader(filename,img->header) == false) return(1);  
    ofstream ofile(imgfname,ios::binary);
    if (!ofile) return(2);
    ofile.write((char *)img->data,sizeof(char)*imgsize);
    ofile.close();
  }
  return(0);

}

// bool copyHeader(AnalyzeHeader* source, AnalyzeHeader* target);

bool copyImage(AnalyzeImage* source,AnalyzeImage* target) 
{
  int imgsize,i;

  target->header = source->header;
  target->nii_header = source->nii_header; // 14 OCT 10
  target->format = source->format;  // 14 OCT 10
  imgsize = (source->header.x_dim)*(source->header.y_dim)*(source->header.z_dim);
  target->data = new float [imgsize];
  if(target->data == NULL) return false;
  for(i = 0;i < imgsize;i++) {
    target->data[i] = source->data[i];
  }
  return true;
}

bool copyLabelImage(AnalyzeLabelImage* source,AnalyzeLabelImage* target)
{
  int imgsize,i;

  target->header = source->header;
  target->nii_header = source->nii_header; // 14 OCT 10
  target->format = source->format;             // 14 OCT 10
  imgsize = (source->header.x_dim)*(source->header.y_dim)*(source->header.z_dim);
  target->data = new unsigned char [imgsize];
  if(target->data == NULL) return false;
  for(i = 0;i < imgsize;i++) {
    target->data[i] = source->data[i];
  }
  return true;
}

bool newImage(AnalyzeImage* img, AnalyzeImage* alike) 
{
  int imgsize,i;

  img->header = alike->header;
  img->nii_header = alike->nii_header; // 14 OCT 10
  img->format = alike->format;             // 14 OCT 10
  imgsize = (alike->header.x_dim)*(alike->header.y_dim)*(alike->header.z_dim);
  img->data = new float [imgsize];
  if(img->data == NULL) return false;
  for(i = 0;i < imgsize;i++) {
    img->data[i] = 0;
  }
  return true;
}

bool newLabelImage(AnalyzeLabelImage* img, AnalyzeImage* alike)
{
  int imgsize,i;

  img->header = alike->header;
  img->nii_header = alike->nii_header; // 14 OCT 10
  img->format = alike->format;             // 14 OCT 10
  img->header.datatype = DT_UNSIGNED_CHAR;
  img->nii_header.datatype =  NIFTI_TYPE_UINT8;
  img->nii_header.nbyper = 1;
  img->nii_header.scl_slope = 1;
  img->nii_header.scl_inter = 0;
  imgsize = (alike->header.x_dim)*(alike->header.y_dim)*(alike->header.z_dim);
  img->data = new unsigned char [imgsize];
  if(img->data == NULL) return false;
  for(i = 0;i < imgsize;i++) {
    img->data[i] = 0;
  }
  return true;
}

bool mask2Image(AnalyzeLabelImage* source, AnalyzeImage* target) 
{
{
  int imgsize,i;

  target->header = source->header;
  target->nii_header = source->nii_header; // 14 OCT 10
  target->format = source->format;             // 14 OCT 10
  imgsize = (source->header.x_dim)*(source->header.y_dim)*(source->header.z_dim);
  target->data = new float [imgsize];
  if(target->data == NULL) return false;
  for(i = 0;i < imgsize;i++) {
    target->data[i] = source->data[i];
  }
  return true;
}
}

bool binaryMask(AnalyzeImage* img,AnalyzeLabelImage* mask,float tr) 
{
  int i;
  int imgsize; 

  if(!newLabelImage(mask,img)) return(false);
  imgsize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim); 
  for(i = 0;i < imgsize;i++) {
    if(img->data[i] > tr) mask->data[i] = 1;
    else mask->data[i] = 0;
  }
  return(true);
}

// finds the closest voxel with the value greater than tr to the voxel (x,y,z)
// extremely slow but correct 
// returns false if no voxel value exceeds the threshold

bool findClosestNonZero(AnalyzeImage* img,int x,int y, int z,int* cx, int* cy, int* cz,float tr)
{
  int i,j,k;
  float d;
  bool found; 
  int windowSize = 11;
  // Initialize
  *cx = x;  // closest x coordinate
  *cy = y;  // closest y coordinate
  *cz = z;  // closest z coordinate
  found = false;

int retries = 0;
while (!found) {
  // if value at (x,y,z) > tr then ok
  if(getVoxelValue(img,*cx,*cy,*cz) > tr) return(true);
  d = pow((float) (img->header.x_dim + img->header.y_dim  + img->header.z_dim),2);
  for(i = (-windowSize);i < (windowSize + 1);i++) {
    for(j = (-windowSize);j < (windowSize + 1) ;j++) {
      for(k = (-windowSize);k < (windowSize + 1);k++) {
         if(getSafeVoxelValue(img,x + i,y + j,z + k) > tr) {
           float fDistanceHere = (float)(i*i + j*j + k*k);
           if(fDistanceHere < d) {
             d = fDistanceHere;
             *cx = x + i;
             *cy = y + j;
             *cz = z + k;
             found = true;
        
           }
         }
       }
     }
   }
  if (!found) {
    cout << "WARNING: No atlas region found within a window size of " << windowSize << "! Retrying with twice the window size." << endl;
			windowSize *= 2;
  }
  if (retries>10) {
    break;
  }
}
   return(found);

}
void averageFilterZeros(AnalyzeImage* img,AnalyzeImage* refImg,float tr,int times)
{
  int x,y,z,i;
  int x1,y1,z1;
  float tmpf;  

  for(i = 0;i < times;i++) {
    for(x = 0;x < img->header.x_dim;x++) {
      for(y = 0;y < img->header.y_dim;y++) {
        for(z = 0;z < img->header.z_dim;z++) {
          if(getVoxelValue(refImg,x,y,z) < tr) {
	    tmpf = 0.0;
            for(x1 = -3;x1 < 4;x1++) {
              for(y1 = -3;y1 < 4;y1++) {
                for(z1 = -3;z1 < 4;z1++) {
                  tmpf = tmpf + getSafeVoxelValue(img,x + x1,y + y1,z + z1);
		}
              }
            }
            putVoxelValue(img,x,y,z,tmpf/(7*7*7));
	  }
	}
      }
    }
  } 
}

// finds a value such that exatcly a certain percentage of voxel values (in img) 
// are smaller than it. Considers only voxels within the mask.

float imageKth(AnalyzeImage* img, AnalyzeLabelImage* mask, float percentage)
{
  long int i,j,k,l,m,n;
  int imgSize;
  float x;
  float* a; 
  
  imgSize = (img->header.x_dim)*(img->header.y_dim)*(img->header.z_dim);
  n = 0;
  for(i = 0;i < imgSize;i++) {
    if(mask->data[i] > 0.5) {
      n++;
    }
  }
  k = (long int) floor(percentage*n);
  if(k > n) return(0);
  if(k < 0) return(0);
  a = new float[n];
  n = 0;
  for(i = 0;i < imgSize;i++) {
    if(mask->data[i] > 0.5) {
      a[n] = img->data[i];
      n++;
    }
  }
  
 l=0 ; m=n-1 ;
 while (l<m) {
   x=a[k] ;
   i=l ;
   j=m ;
   do {
     while (a[i]<x) i++ ;
     while (x<a[j]) j-- ;
     if (i<=j) {
       FLOAT_SWAP(a[i],a[j]) ;
       i++ ; j-- ;
     }
   } while (i<=j) ;
   if (j<k) l=i ;
   if (k<i) m=j ;
 }
 x = a[k];
 delete[] a;
 return(x);
}

// Median filters image img; 
// Assumes windowLengths are odd integers
// for example (medianFilter3D(..., ..., 1, 1, 1);
// uses 3 x 3 x 3 mask.

void medianFilter3D(AnalyzeImage* img, AnalyzeLabelImage* mask, int windowLenX, int windowLenY, int windowLenZ)
{
  int windowSize;
  float* window;
  int x,y,z,i,j,k,l,m;
  int n;
  int med;
  float f;

  windowSize = (2*windowLenX + 1)*(2*windowLenY + 1)*(2*windowLenZ + 1);  
  med = (windowSize - 1)/2;
  window = new float[windowSize];


  for(x = 0;x < img->header.x_dim;x++) {
    for(y = 0;y < img->header.y_dim;y++) {
      for(z = 0;z < img->header.z_dim;z++) {
        if(getLabelValue(mask,x,y,z) > 0.5) {
          n = 0;
          for(i = (-windowLenX);i < (windowLenX + 1);i++) {
            for(j = (-windowLenY);j < (windowLenX + 1);j++) {
              for(k = (-windowLenZ);k < (windowLenX + 1);k++) {
                window[n] = getSafeVoxelValue(img,x + i, y + j, z + k);
                n++;
              }
            }
          }
          l=0 ; m=n-1 ;
          while (l<m) {
            f=window[med] ;
            i=l ;
            j=m ;
            do {
              while (window[i]<f) i++ ;
              while (f<window[j]) j-- ;
              if (i<=j) {
                FLOAT_SWAP(window[i],window[j]) ;
                i++ ; j-- ;
              }
            } while (i<=j) ;
            if (j<med) l=i ;
            if (med<i) m=j ;
	  }
	  putVoxelValue(img,x,y,z,window[med]);
        }
      }
    }
  }
  delete[] window;
}

void knnFilter(AnalyzeImage* img, AnalyzeLabelImage* mask, int windowLenX, int windowLenY, int windowLenZ, int knn)
{
  float* window;
  float* distances;
  int x,y,z,i,j,k;
  int n;
  float f,sum,voxelVal;
  char c;

  window = new float[knn]; 
  distances = new float[knn];

  for(x = 0;x < img->header.x_dim;x++) {
    for(y = 0;y < img->header.y_dim;y++) {
      for(z = 0;z < img->header.z_dim;z++) {
        if(getLabelValue(mask,x,y,z) > 0.5) {
          n = 0;
          voxelVal = getVoxelValue(img,x,y,z);
          for(i = (-windowLenX);i < (windowLenX + 1);i++) {
            for(j = (-windowLenY);j < (windowLenX + 1);j++) {
              for(k = (-windowLenZ);k < (windowLenX + 1);k++) { 
		if(n < knn) {
                  window[n] = getSafeVoxelValue(img,x + i, y + j, z + k);
                  distances[n] = fabs(window[n] - voxelVal);
                   n++;
                }
                else {
                  c = maxArg(distances,knn);
                  f = fabs(getSafeVoxelValue(img,x + i, y + j, z + k) - voxelVal);
                  if(f < distances[c]) {
                    distances[c] = f;
                    window[c] = getSafeVoxelValue(img,x + i, y + j, z + k);
		  }
                }
              }
            }
          }
          sum = 0.0;
          for(i = 0;i < knn;i++) {
            sum = sum + window[i];
          }
          putVoxelValue(img,x,y,z,sum/((float) knn));
        }
      }  
    }
  }
  delete[] window;
  delete[] distances;
}
