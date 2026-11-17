# ihmt-push
MT-weighted data and PUSH pulses adapted for universal pulse design. This repo contains some of the scripts that were adapted from the original PUSH code (https://github.com/mriphysics/PUSH/) for UP pulse design. If you find it useful, please cite the publication once available.

## If used, please reference the following publication:
1. Leitão D, Tomi-Tricot R, Bridgen P, et al. Parallel transmit pulse design for saturation homogeneity (PUSH) for magnetization transfer imaging at 7T. Magn Reson Med. 2022; 88: 180-194. doi:10.1002/mrm.29199
2. ihMT PUSH publication once available

## Requirements:
 - Derived from code: https://github.com/mriphysics/PUSH/
 - Requires proprietary Siemens code, which may be made available upon request, subject to an agreement with the vendor

## INI folder
- contains RF waveforms,which have been scaled to saturation flip angle (FA) of 20 deg (Nominal FA)
- The target saturation flip angle is the in filename, therefore the RF waveform in the INI file can be scaled by targetFA/NominalFA, i.e. 191/20 to achieve target FA
- Can scaled to even higher or lower FA depending on SAR constraints or if higher FA is desired! 

## DATA folder
- contains pre-processed MT-weighted images for Subjects 1,2,3
- ### Subject X:
  - ### 2 kHz 
      - ants_CP_Spos ( pre-processed single-sided CP saturation MTw image )
      - ants_UP_Spos ( pre-processed single-sided UP saturation MTw image )
      - ants_SS_Spos ( pre-processed single-sided SS saturation MTw image )
      - ants_Sref    ( pre-processed no saturation reference image )
  - ### 7 kHz 
      -   ants_CP_S_4d
      -   ants_UP_S_4d
      -   ants_SS_S_4d
      #### These 4D datasets are organized as the following:
            - Sref_reg = S_4d[:,:,:,0]  ( pre-processed no saturation 3D reference image )
            - Spos_reg = S_4d[:,:,:,1]  ( pre-processed single-sided saturation MTw image (+7 kHz) )
            - Sneg_reg = S_4d[:,:,:,2]  ( pre-processed single-sided saturation MTw image (-7 kHz) )
            - Sdual_reg = S_4d[:,:,:,3] ( pre-processed dual-sided saturation MTw image (+/-7 kHz) )
   
## PUSH_UP folder
- To be used with the original PUSH scripts found here: https://github.com/mriphysics/PUSH/
- Requires proprietary Siemens code, which may be made available upon request, subject to an agreement with the vendor
- ### run_PUSH_UP_optimization:
        - This script modifies PUSH_optimisation.m to take in multiple B1 field maps and the MTboost saturation scheme
        - Adj.S: Channel-wise B1 maps 







