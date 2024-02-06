!*******************************************************************************
                       MODULE modGlobalPFC
!*******************************************************************************
IMPLICIT NONE

!GLOBAL VARIABLES
INTEGER, PARAMETER :: DP = KIND(1.D0) !Declaring Double Precision
!Variable used while opening files:
INTEGER             :: error 
CHARACTER(LEN=50)   :: fileIn
CHARACTER(LEN=50)   :: fileOut

!Initializing input file units(10-19):
INTEGER, PARAMETER  ::  fNoInp = 10                     

!Initializing output file units(20-29):
INTEGER, PARAMETER  ::fNoOut = 20

INTEGER  :: xSide, ySide !dimensions of the lattice x,y
INTEGER  :: xHalfp1, yHalfp1 !x half  plus 1, y half  plus 1
REAL(DP) :: r !parameter used in Eqn of motion

!Constants used in the code
REAL(DP), PARAMETER :: piby4 = DATAN(1.D0)
REAL(DP), PARAMETER :: pi = 4.D0*DATAN(1.D0)
REAL(DP), PARAMETER :: v0 = 0.6D0
REAL(DP), PARAMETER :: l  = 64.D0
!Dynamical array allocation of variables used in the code

!psi(x,y), FreeNrg(x,y) represnts the density field & the free energy
REAL(DP), ALLOCATABLE :: psi(:,:), psiDot(:,:),freeNrg(:,:)
REAL(DP), ALLOCATABLE :: psiK(:,:), psiKDot(:,:)
!REAL(DP), ALLOCATABLE :: psiold(:,:), psiDotold(:,:)

!Shear profile variable v
REAL(DP), ALLOCATABLE :: v(:)

!nu1, nu2 are two noise generated to get the final conserved noise zeta 
!in real space
REAL(DP), ALLOCATABLE :: nu1(:,:),nu2(:,:),zeta(:,:)
REAL(DP) :: dx,dx_Inv,dxSq_Inv,dy,dy_Inv,dysq_Inv,boxSize_Inv 
REAL(DP) :: t,dt,dt_SqRoot,dtSqRt_Inv,tFinal
REAL(DP) :: sdNoise,sdNoiseV,sdNoise2
REAL(DP) :: psiAvg 
INTEGER  :: iter
CONTAINS

!===============================================================================
                   SUBROUTINE SetInputs
!===============================================================================
IMPLICIT NONE
CHARACTER(LEN=50)::fileIn

!fileIn = 'inputVPFC'
fileIn = 'inputNigel'
!fileIn = 'input'
  OPEN(1,FILE= fileIn)
    READ(1,*) xSide
    READ(1,*) ySide
    READ(1,*) dx
    READ(1,*) dy
    READ(1,*) dt
    READ(1,*) tFinal
    READ(1,*) psiAvg
    READ(1,*) r
    READ(1,*) sdNoiseV
  CLOSE(1)
!dx = piby4
!dy = piby4
ALLOCATE(psi(xSide,yside))
ALLOCATE(psiDot(xSide,yside))
ALLOCATE(psiK(xSide,yside))
ALLOCATE(psiKDot(xSide,yside))
ALLOCATE(zeta(xSide, ySide))
ALLOCATE(nu1(xSide, ySide))
ALLOCATE(nu2(xSide, ySide))
ALLOCATE(v(ySide))

!initialize specific parameters 
iter=0
t=0.0d0
!nu1=0.D0
!nu2=0.D0
!zeta=0.D0
psi = 0.D0
psiDot = 0.D0
psiK = 0.D0
psiKDot = 0.D0

  !define inverse square mesh size
  dx_Inv = 1.D0/dx
  dxSq_Inv = dx_Inv*dx_Inv 
  dt_SqRoot = DSQRT(dt)
  dtSqRt_Inv = 1.D0/dt_SqRoot
  dy_Inv = 1.D0/dy
  dySq_Inv = dy_Inv*dy_Inv 
  xHalfp1 = xSide/2 + 1
  yHalfp1 = ySide/2 + 1
  boxSize_Inv = 1.D0/(DFLOAT(xSide*ySide))
 WRITE(*,*)'Global data:'
 WRITE(*,*)'_________________________________'
 WRITE(*,10)  xSide,ySide,dx,dy,dt,tFinal,psiAvg,r,sdNoiseV
 10 FORMAT(' xSide =',I4,/ ' ySide =', I4, /  &
             ' dx =',F6.4, / ' dy =',F6.4, / ' dt =',F6.4,/ &
             ' tFinal =',F10.4,/  ' psiAvg =',F8.4,/ ' r =',F10.4,/ &
              ' sdNoiseV =',F6.4,/)
           
 WRITE(*,*)'_________________________________'


!===============================================================================
                   END SUBROUTINE SetInputs
!===============================================================================

!===============================================================================
                   SUBROUTINE FileError(fileName,error,modul)
!===============================================================================
IMPLICIT NONE

!Shared
INTEGER::error
CHARACTER(LEN=*)::fileName
CHARACTER(LEN=*),OPTIONAL :: modul

IF(error/=0)THEN
  WRITE(*,1)
  WRITE(*,'(" ===>> FATAL ERROR :: file not found : ")')fileName
  IF(PRESENT(modul))WRITE(*,'(1x,"           Module :: ",A20,/)')modul
  WRITE(*,1)
  STOP
ENDIF

1 FORMAT(1x,'*==================================================',&
            '==================*') 
!===============================================================================
                   END SUBROUTINE FileError
!===============================================================================

!===============================================================================
                        SUBROUTINE Fatal(routine,string,modul)
!===============================================================================
IMPLICIT NONE

!Shared
CHARACTER(LEN=*)  :: routine, string
CHARACTER(LEN=*),OPTIONAL :: modul


WRITE(*,*)
WRITE(*,1)
WRITE(*,'(1x, "PROGRAM FAILED")')

WRITE(*,*)string
WRITE(*,'(1x,"===>> FATAL ERROR in the routine :",A30)')routine
IF(PRESENT(modul))WRITE(*,'(5x,"       MODULE :",A30,/)')modul
WRITE(*,1)

STOP

1 FORMAT(1x,' ==================================================',&
             '================== ') 
!===============================================================================
                           END SUBROUTINE Fatal
!===============================================================================

!*******************************************************************************
                          END MODULE modGlobalPFC
!*******************************************************************************
