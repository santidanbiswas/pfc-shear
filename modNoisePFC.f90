!******************************************************************************
                         MODULE modNoisePFC
!******************************************************************************

USE modGlobalPFC

IMPLICIT NONE ; SAVE

CONTAINS

!===============================================================================
                   SUBROUTINE GasDev(sdNoise,rf)
!To generate the random noise
!===============================================================================
#ifdef __INTEL_COMPILER
use IFPORT
#endif

IMPLICIT NONE
INTEGER,PARAMETER :: DP = KIND(1.D0)
       
REAL(DP),INTENT(OUT):: rf
REAL(DP),INTENT(IN)::sdNoise
REAL(DP)::rsq,v1,v2,fac,deviate,gset!,rand
INTEGER::iset
SAVE iset,gset
DATA iset/0/
!WRITE(*,*)'sdNoiseSUBrtn1:',sdNoise
!DO j = 1,3
!  Returns a normally distributed deviate with zero mean and variance sigma
   IF( iset == 0) THEN
!891   v1 = 2.0d0*ranF()-1.0d0
!      v2 = 2.0d0*ranF()-1.0d0
891   v1 = 2.0d0*rand()-1.0d0
      v2 = 2.0d0*rand()-1.0d0

      rsq = v1**2.0d0+v2**2.0d0
      IF (rsq >= 1.0d0 .OR. rsq == 0.0d0) GOTO 891
      fac=DSQRT(-2.0d0*(DLOG(rsq))/rsq)
      gset=v1*fac
      deviate=v2*fac
      iset=1
   ELSE
      deviate=gset
      iset=0
   ENDIF
   rf=deviate*sdNoise
!ENDDO
!WRITE(*,*)'sdNoiseSUBrtn2:',sdNoise

RETURN
!===============================================================================
                       END SUBROUTINE GasDev
!===============================================================================

!===============================================================================
                   SUBROUTINE GetNoiseTable
!To generate the conserved noise table i.e noise value at each site
!===============================================================================
IMPLICIT NONE

INTEGER::i,j
REAL(DP) :: noise1

DO j = 1, ySide
  DO i = 1, xSide
    CALL GasDev(sdNoiseV,noise1)
    nu1(i,j) = noise1
    CALL GasDev(sdNoiseV,noise1)
    nu2(i,j) = noise1
  ENDDO
ENDDO

RETURN
!===============================================================================
                  END SUBROUTINE GetNoiseTable
!===============================================================================

!===============================================================================
                 SUBROUTINE GetRandNoise
!===============================================================================
IMPLICIT NONE

INTEGER :: i,j
INTEGER :: xplus,yPlus

DO j = 1, ySide
  DO i = 1, xSide
! Implementing Periodic Boundary condition
  xPlus = i + 1
  IF (i==xSide) xPlus = 1
  yPlus = j + 1
  IF (j==ySide) yPlus = 1
  
  zeta(i,j) = dx_Inv*(nu1(xplus,j) - nu1(i,j) + nu2(i,yPlus) &
                     - nu2(i,j))
  
  ENDDO
ENDDO

RETURN
!===============================================================================
                      END SUBROUTINE GetRandNoise
!===============================================================================

END MODULE modNoisePFC

