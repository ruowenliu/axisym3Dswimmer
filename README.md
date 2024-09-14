# axisym3Dswimmer

## This repository stores the programs for:

_Shape optimization of slip-driven axisymmetric microswimmers_

## This repository contains the following parts
1. Problem solver (kernals, ect.)

  The folder **quadrature_and_kernal** contains the programs for solving the forward (and adjoint) problem. 
   
  _brewermap_ (c) 2014 Stephen Cobeldick is used to control colorbars.
   
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

   **max_eff_060 is kind of done, run "test_solo.m". to-do: vary sigma for different nu, polish test_solo to get dimvec from "final_designvec" txt file, etc.**

