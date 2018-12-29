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
 In this code, we use strip model (overlap area) to compute the weight of pixels.
 - Area of intersection between a strip of finite width and a pixel.
 - Can be computed fast by pre-computing the pixel footprint.
 
 the possible overplap situation can be complicated, we need discuss each pixel belong to which overlap situation.
 
 the main idea of this code is firstly get each pixel's four corners belong to which strips and record it. Then use the results to compute the area of intersection. For pixels along to one ray (strip), we sort it according to its distance to the detector plane in order to compute the attenuation map easiler.
 
 we seperate the whole process into two part 1. attenuation mat A 2. system mat H 
 
 - compute_system_mat: 
   basically the main function of computing system matrix process.
 - assign_grid: 
   pre-compute which stripe points on grid belong to for every angles.
 - pdistlist: 
   a link list to record all pixels along every strips.
 - calcu_weight_parallel: 
   according to different intersection cases to get the area of intersection between a strip and a pixel.
 
