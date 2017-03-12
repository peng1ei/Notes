This set of files contains the Matlab code
for VCA algorithm describe in the paper:

Jos?M. P. Nascimento and Jos?M. B. Dias 
"Vertex Component Analysis: A Fast Algorithm to Unmix Hyperspectral Data"
 submited to IEEE Trans. Geosci. Remote Sensing, vol. .., no. .., pp. .-., 2004

A short version can be found in  
Jos?M. P. Nascimento and Jos?M. B. Dias
"Vertex Component Analysis: A fast algorithm to extract endmembers spectra from Hyperspectral data"
in Pattern Recognition and Image Analysis, ser. Lecture Notes in Computer Science,
F. j. Perales, A.J.C. Campilho, N. P. Blanca and A. Sanfeliu, eds.,
vol. 2652, Springer-Verlag, 2003, pp. 626-635.


Files:
       readme.txt - this file
     demo_vca.m   - Matlab script file to start demo
          vca.m   - Matlab function file with VCA algorithm
dirichlet_rnd.m   - Matlab function file for generation of 
                    abundance fractions（丰度分数） with 
                    dirichlet distributions 
   signatures.mat - Matlab data file with a set of mineral signatures（矿物标志）
                    extracted（提取） from USGS spectral library.
                    This file contains three variables:
                    wavlen (224 x 1) - wavelength
                    A (224 x 21)     - 21 mineral signatures
                    names (21 x 29)  - names of the minerals
      cup_ref.mat - Matlab data file with the Cuprite, Nevada
                    reflectance subimage (250 x 190 pixels)
                    from data set acquired on the AVIRIS flight 19 June 1997
                    the subimage starts at line 1620 and column 420 and it 
                    ends on line 1869 and column 610, noisy channels 
                    {1, 2, 104...113, 148...167, 221...224} were removed.

                    This file contains three variables:
                    wavlen (224 x 1) - wavelength
                    BANDS  (1 x 188) - selected bands
                    Lines            - number of lines of the subimage
                    Columns          - number of columns of the subimage
                    L                - number of channels selected
                    x  (188 x 47750) - subimage (channels x number of pixels)
                    
Getting started:

To run different demos, set the "demo_number" on the first line of the demo_vca.m file:
demo_number=1
demo_number=2
demo_number=3

this will set different simulation scenarios（模拟场景）:
demo 1:  simulated mixing of  3 endmembers. The miximig in piecewise 
         (strips in an image) constant;
demo 2:   simulated mixing of  3 endmembers. Abundances follow a Dirichlet distributions;
demo 3:  real data from Cuprite (AVIRIS).

Run the demo_vca.m script Matlab file


Programming with VCA:

Modify demo_vca.m to personalize your demo:
Parameters
p (number of endmembers);
signatures (A) (from 1 to 21);
dimension of the image (Lines and Columns);
the signal-to-noise ration (SNR) in dB;
and the illumination perturbation (q). 
