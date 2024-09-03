# axisym3Dswimmer

## This repository stores the programs for the preprint:

_Shape optimization of slip-driven axisymmetric microswimmers_ [arxiv](https://arxiv.org/abs/2405.00656/)

## This repository contains the following parts
1. Problem solver (kernals, ect.)

   The _solver_ folder contains the programs for solving the forward (and adjoint) problem. 

   Some despcription of the programs are listed below.
   
   _brewermap_ (c) 2014 Stephen Cobeldick is used to control colorbars.

2. Calculation of maximum efficiency on prolate spheroids with various reduced volumes

   **To run the program, use "runProlateSpheroid.m"**

3. Parametrization of shapes

   **To run the program, use xxxx.m**

4. Shape sensitivity verifications

   **To run the program, go to the folder "verify_sensitivity", then run "sens_verify_main.m". It reproduces the main results in Section 5.1.**

5. Numerical results

   **In the folder "example_max_efficiency", run the main program "max_efficiency_main.m" directly can repeat the iterations in Section 5.2. The initial shape is set as a peanut-like swimmer, then the shape is optimized to obtain the maximum swimming efficiency.**
   
   **In the folder "min_drag_force", the process and results are recorded in .txt file, labeled by the value of nu times 100, for example, file "..._nu_060.txt" corresponds to nu=0.60. You can run "mindragforce_main_program.m" to reproduce all the results, reading preset data from input_parameters.csv automatically. Additional, the design vectors of the process are also storee in .txt file.**

   **max_eff_060 is kind of done, run "test_solo.m". to-do: vary sigma for different nu, polish test_solo to get dimvec from "final_designvec" txt file, etc.**

