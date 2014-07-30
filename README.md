SVPASEG
=======

A tool written in C implementing a flexible MRI brain tissue classification framework. Various aspects of the algorithm 
can be controlled by editing a configuration file. Deals with the intensity nonuniformities by using several 
voxel intensity models across the brain. These models are then combined by  novel inhomogeneous MRF. See

J. Tohka , I.D. Dinov, D.W. Shattuck, and A.W. Toga. 
Brain MRI Tissue Classification Based on Local Markov Random Fields, 
Magnetic Resonance Imaging, 28(4): 557 - 573 , 2010.

http://www.cs.tut.fi/~jupeto/atlas_mrf_mri_postprint.pdf

for more details.

The software is meant for obtaining a hard tissue classification that is easily tuned for specific applications.
If, however, you're interested in partial volume estimation in high quality T1-weighted brain MRI, I would suggest to try 

https://github.com/jussitohka/pvemri 

running in Matlab.

INSTALLATION
============

1. unpack 
2. cd src
3. make all
4. cd ../atlases
5. unzip atlas-files, assuming the zip-files (2 of them) have been downloaded to this directory, 
6. cd ..
7. run configure_atlases

for more info check the documentation under doc subdir
