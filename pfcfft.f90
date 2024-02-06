!############## HOW  TO COMPILE THIS CODE ##################
! type "make" in terminal without the quotation marks

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     PROGRAM Test_Nucleation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
USE modNoisePFC
USE modCalcPFC
!###########################################################
!Purpose:
!Liquid-Soild system
!Growth of a crystalline phase from a supercooled liquid
!Reference: K. Elder and M. Grant
!            Phys. Rev. E. 70, 051605 (2004)
!Title: Modeling elastic & plastic deformation in non eqm 
!       processing using phase field crystals
!###########################################################



IMPLICIT NONE
include 'fftw3.f'

CALL SetInputs
CALL InitConfig
!CALL InitConfigRndSeed
CALL EOM

END PROGRAM Test_Nucleation

