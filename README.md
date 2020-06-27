## Non-Linear Fit to Produce T1-maps from the Modified Look-Locker Inversion Recovery (MOLLI) Method 

# molli_nonlinear.m

This repo contains molli_nonlinear.m which should ask you to select your folder where your raw MOLLI dicoms are stored (you should have 11 dicoms in that folder), will perform a nonlinear fit to estimate the T1 at each voxel and will save a T1 map as a dicom. 

**This script requires:**
- MATLAB's Curve Fitting Toolbox: https://www.mathworks.com/products/curvefitting.html
- DICOM Series reader: https://www.mathworks.com/matlabcentral/fileexchange/71807-dicom-series-reader 
- MATLAB's Image Processing Toolbox https://www.mathworks.com/products/image.html.


Note - the molli_nonlienar_par_absolute.m file stored here might not calculate the T1s correctly - I am not sure, I haven't checked. There were some issueswith previous MOLLIs and the changes made to molli_nonlinear made the fit work well, so perhaps use that instead for now.
