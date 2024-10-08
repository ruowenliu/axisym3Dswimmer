# axisym3Dswimmer

## This repository stores the programs for:

_Shape optimization of slip-driven axisymmetric microswimmers_

## This repository contains the following parts
1. Problem solver (kernals, ect.)

   Folder **quadrature_and_kernal** contains the programs for solving the forward (and adjoint) problem. 
   
2. Functions for shape representation are stored in the folder **shape_classes**

3. Shape sensitivity verifications

   folder **verify_sensitivity**
   
   run **sens_verify_main.m** 
    
   It reproduces the main results in Section 5.1. The results are stored in **sens_verify_result.txt**.

5. Numerical results

  5.1. Calculation of maximum efficiency on prolate spheroids with various reduced volumes

   folder **prolate_spheroids**
  
   run **runProlateSpheroid.m**
   
  5.2. Shape optimization from an arbitrary shape

   folder **maxE_from_peanut**
  
   run **max_efficiency_main.m** directly 
   
   The initial shape is set as a peanut-alike swimmer, then the shape is optimized to obtain the maximum swimming efficiency.
 
  5.3. Library of optimal shapes for minimum drag force, with various reduced volume

   folder **min_drag_1** and folder **min_drag_2**
   
   run **min_drag_main_program.m**
   
   It reads the preset data from input_parameters.csv automatically. The process and results are stored in .txt file, labeled by the value of nu times 100, for example, file "..._nu_060.txt" corresponds to nu=0.60. Additional, the design vectors are also storee in .txt file.

  5.4. Library of optimal shapes for maximum swimming efficiency, with various reduced volume

   folder **maxE_nu_050etc**, run **maxE_main_program.m** for nu = 0.50, 0.525, 0.55, 0.575
   
   folder **maxE_nu_060**, run **maxE_060.m** for nu = 0.60. Same for the rest nu up to 0.95.
   
   Note: Reuse the results of nu=0.70 from the previous peanut-alike example for best E.
   
6. Generate figures (velocity fields)

   folder **plot_figures** stores the codes to generate figures

  _brewermap_ (c) 2014-2022 Stephen Cobeldick is used to control colorbars.
  
7. (Optional) Test for Convergence of the Bspline representation

   folder **verify_Bspline**
   
   There are 3 programs can be used for checking the well-posedness of Bspline for sharp poles, the test is based on the iteration 17 in the peanut-like optimization. 
   
   program #1: **plot_a_single_shape.m** generates figures of the slip profile and geometry R, Z 

   program #2: **change_NL_testGeometry.m** generates figures of convergence test for R, Z for various NL
      
   program #3: **change_NL_testE.m** generates figures of errors in R, Z, E, and drag force for various NL
   
