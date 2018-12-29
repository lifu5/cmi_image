# SPECT system matrix with attenuation

(..underconstucted)
## part 1 basic functions

basic structure
 - cmi_image : 
    name: the description of image
    n_row,n_col,n_depth: the size of image (n_depth==1 means 2d image)
    data_type: support four different type (bool, int, float, double)
    data: raw data void* type 
     
basic function
 - allocimage: 
   constructe method, similar to zeros function in MATLAB
 - allocimage2d: 
   construct a 2d image.
 - writerawimage:
   save image data as binary file.
 
 
 ## part 2 system matrix calculation 
