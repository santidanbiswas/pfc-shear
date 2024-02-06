!*******************************************************************************
                          MODULE modCalcPFC
!*******************************************************************************
USE modGlobalPFC
USE modNoisePFC
CONTAINS
!===============================================================================
                        SUBROUTINE InitConfig
!===============================================================================
IMPLICIT NONE
INTEGER :: i,j,k,ii,jj
REAL(DP):: sdNoise1,noise1
!OPEN(3,FILE='new.txt')
 
  DO j = 1, ySide
    DO i = 1, xSide
      !READ(3,*)ii,jj,psi(i,j)
      psi(i,j) = psiAvg 
    ENDDO
      !IF(i==xSide) READ(3,*) 
  ENDDO
!CLOSE(3)

sdNoise1 = 1.D0
!   only for a single nucleus at centre
 DO j = ySide/2 -5,ySide/2 + 5
  DO i =  xSide/2 -5,xSide/2 + 5
   CALL GasDev(sdNoise1,noise1)
 !psi(i,j) = psiAvg
    psi(i,j) = psi(i,j) + 0.1D0*noise1
  ENDDO
 ENDDO
OPEN(1,FILE='Output/outtemp.dat',STATUS='unknown')
CALL OutputPsi(psi,1)
CLOSE(1)

!############################ Addition of psidot initial condition#######
CALL GetNoiseTable
CALL GetRandNoise
  DO j = 1, ySide
    DO i = 1, xSide
      psiDot(i,j) = 0.D0 !zeta(i,j)*dtSqRt_Inv 
      !WRITE(*,*) i,j,psiDot(i,j)
    ENDDO
  ENDDO
!OPEN(1,FILE='Output/out0.dat',STATUS='unknown')
!     CALL OutputPsi(psi,1)
!CLOSE(1)
!##############################Shearing Profile###########################
DO j = 1, ySide/2
    v(j) = v0*DEXP(-j/l) 
ENDDO
DO j = yHalfp1,ySide
    v(j) = -1.D0*v0*DEXP(-(ySide-j)/l) 
ENDDO

RETURN
!==============================================================================
                        END SUBROUTINE InitConfig
!===============================================================================
!===============================================================================
                        SUBROUTINE InitConfigRndSeed
!===============================================================================
IMPLICIT NONE
INTEGER :: i,j,k,ii,jj
REAL(DP):: d1,d2,rad,mid_sys1,mid_sys2,qo,theta,xr,yr
REAL(DP):: sdNoise1,noise1
REAL(DP):: d,seedXpos,seedYpos,theta1,c1
INTEGER, ALLOCATABLE :: flagSite(:,:)

ALLOCATE(flagSite(xSide,ySide))
flagSite = 0
psi = psiAvg
sdNoise1 = 1.D0
!mid_sys1 = xSide/4.D0
!mid_sys2 = 3.D0*xSide/4.D0
qo = 0.8660D0
rad = 15.D0
theta1 = 2.D0*pi/3.D0
!OPEN(3,FILE='out15000.dat')
c1 = 0.D0 
DO ii = 1,3
  seedXpos = 160*ii
  DO jj = 1,6
    seedYpos = 160*jj
    IF(MOD(ii,2)==1) THEN
      c1 = c1 + 1.D0
      theta = theta1/c1
    ENDIF
    IF(MOD(ii,2)==0) THEN
      c1 = c1 - 1.D0
      theta = theta1/(c1+1.D0)
    ENDIF
  DO j = 1, ySide
    DO i = 1, xSide
      !READ(3,*)ii,jj,psi(i,j)
      !psi(i,j) = psiAvg
      !d1=sqrt((i*dx-mid_sys1*dx)**2+(j*dy-mid_sys1*dy)**2)-rad*dx
      

      d=DSQRT((i*dx-seedXpos*dx)**2+(j*dy-seedYPos*dy)**2)-rad*dx
      !if(i.gt.mid_sys1-Rad.and.i.lt.mid_sys1+Rad)then !case 1: two long strips of solid
      !IF(d1 <= 0.0D0) THEN  !case 2: two circular seeds
                            ! small sinusoidal seed at left
       !       psi(i,j) = 0.5D0*( DCOS(2.D0*qo*(j-1)*dx/DSQRT(3.0D0))/2.0D0 &
       !     - DCOS(qo*(i-1)*dx)*DCOS(qo*(j-1)*dx/DSQRT(3.0D0))) + psi(i,j)
            !else if(i.gt.mid_sys2-Rad.and.i.lt.mid_sys2+Rad)then !case 1
      !ELSE
      IF(d <= 0.0d0)THEN !case 2:
                        ! small sinusoidal seed at right, misoreinted by theta
            xr=  (i-1)*dx*DCOS(theta)+(j-1)*dx*DSIN(theta)
            yr= -(i-1)*dx*DSIN(theta)+(j-1)*dx*DCOS(theta)
            psi(i,j) = 0.5*(DCOS(2*qo*yr/DSQRT(3.0D0))/2.D0 &
                       - DCOS(qo*xr)*DCOS(qo*yr/DSQRT(3.0D0))) + psi(i,j)
            flagSite(i,j) = 1
      !ELSE    !  everywhere else
      !     CALL GasDev(sdNoise1,noise1)
           !psi(i,j) = psi(i,j) + 0.0001D0*noise1
      ENDIF
    ENDDO
  ENDDO
ENDDO
ENDDO
!CLOSE(3)

!sdNoise1 = 1.D0
!   only for a single nucleus at centre
! DO j = ySide/2 -5,ySide/2 + 5
!  DO i =  xSide/2 -5,xSide/2 + 5
  !psi(i,j) = psiAvg
!    psi(i,j) = psi(i,j) + 0.1D0*noise1
!  ENDDO
! ENDDO
OPEN(1,FILE='Output/outtemp.dat',STATUS='unknown')
CALL OutputPsi(psi,1)
CLOSE(1)

!############################ Addition of psidot initial condition#######
CALL GetNoiseTable
CALL GetRandNoise
  DO j = 1, ySide
    DO i = 1, xSide
      psiDot(i,j) = 0.0D0 !zeta(i,j)*dtSqRt_Inv 
      !WRITE(*,*) i,j,psiDot(i,j)
    ENDDO
  ENDDO
!OPEN(1,FILE='Output/out0.dat',STATUS='unknown')
!     CALL OutputPsi(psi,1)
!CLOSE(1)
!##############################Shearing Profile###########################
DO j = 1, ySide/2
    v(j) = v0*DEXP(-j/l) 
ENDDO
DO j = yHalfp1,ySide
    v(j) = -1.D0*v0*DEXP(-(ySide-j)/l) 
ENDDO
DEALLOCATE(flagSite)
RETURN
!==============================================================================
                        END SUBROUTINE InitConfigRndSeed
!===============================================================================
!===============================================================================
                             SUBROUTINE EOM
!===============================================================================
IMPLICIT NONE

INTEGER :: i,j,il
CHARACTER(LEN=9) :: cn
REAL(DP) :: psi3(xSide,ySide)  
REAL(DP) :: tShear(xSide,ySide)
!OPEN(1,FILE='Output/outtemp.dat',STATUS='unknown')
!CALL OutputPsi(psi,1)
!CLOSE(1)

DO 
  t = t + dt
  iter = iter + 1  
  IF(t > tFinal) THEN
    WRITE(*,*)t,tFinal
    EXIT
  ENDIF
!  DO j=1, ySide
!    DO i=1,xSide
!      WRITE(*,*)i,j,k,psi(i,j,k)
!    ENDDO 
!  ENDDO

     
!#################### Shear Term ###############################
tShear = 0.D0
IF (t>1000.D0) THEN
DO j = 1,ySide
  DO i =1,xSide
    tShear(i,j) = v(j)*psi(i,j)
  ENDDO
ENDDO
ENDIF
!##################################################

     CALL Main(psi,psiDot,tShear)
      psi = psi*boxSize_Inv
      psiDot = psiDot*boxSize_Inv

!     OPEN(1,FILE='Output/outtemp.dat',STATUS='unknown')
!      CALL OutputPsi(psi,1)
!     CLOSE(1)

  !initialize to zero average for this time step
  psiAvg=0.0d0
    DO j=1, ySide
      DO i=1,xSide
        psiAvg = psiAvg + psi(i,j)
      ENDDO 
    ENDDO 
  psiAvg = psiAvg*boxSize_Inv

  WRITE(2,*) psiAvg
  IF (MOD(iter,10000)+1 == 1) THEN
     CALL chari(iter,cn,il)
     OPEN(1,FILE='Output/out'//cn(1:il)//'.dat',STATUS='unknown')
       CALL OutputPsi(psi,1)
     CLOSE(1)
  ENDIF  
ENDDO
CLOSE(2)  !Close the file opened to store average psi
RETURN
!===============================================================================
                        END SUBROUTINE EOM
!===============================================================================


!===============================================================================
                       SUBROUTINE Main(dummyF,dummyFdot,tShear)
!===============================================================================
IMPLICIT NONE
include 'fftw3.f'
REAL(DP), INTENT(INOUT) ::dummyF(xSide,ySide) 
REAL(DP), INTENT(INOUT) ::dummyFdot(xSide,ySide)  
REAL(DP), INTENT(IN)    ::tShear(xSide,ySide)
REAL(DP), ALLOCATABLE :: in1(:,:),in2(:,:),in3(:,:),in4(:,:)!,in2(:,:)
DOUBLE COMPLEX,ALLOCATABLE :: out1(:,:),out2(:,:),out3(:,:),out4(:,:)
INTEGER :: i,j
INTEGER*8 ::plan

ALLOCATE(in1(xSide,ySide))
ALLOCATE(in2(xSide,ySide))
ALLOCATE(in3(xSide,ySide))
ALLOCATE(in4(xSide,ySide))
ALLOCATE(out1(xhalfp1, ySide))
ALLOCATE(out2(xhalfp1, ySide))
ALLOCATE(out3(xhalfp1, ySide))
ALLOCATE(out4(xhalfp1, ySide))

in1 = dummyF
in2 = dummyF*dummyF*dummyF
in3 = dummyFdot
in4 = tShear
! Forward Fourier transform
call dfftw_plan_dft_r2c_2d(plan,xSide,ySide,in1,out1,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, in1, out1)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan,xSide,ySide,in2,out2,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, in2, out2)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan,xSide,ySide,in3,out3,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, in3, out3)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan,xSide,ySide,in4,out4,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, in4, out4)
call dfftw_destroy_plan(plan)

CALL CalcEq(out1,out2,out3,out4)

call dfftw_plan_dft_c2r_2d(plan,xSide,ySide,out1,in1,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, out1, in1)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_2d(plan,xSide,ySide,out3,in3,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, out3, in3)
call dfftw_destroy_plan(plan)

dummyF = in1
dummyFdot = in3
DEALLOCATE(in1,in2,in3,out1,out2,out3)
RETURN
!===============================================================================
                       END SUBROUTINE Main
!===============================================================================

!==============================================================================
SUBROUTINE CalcEqOld(out1,psik3)
!           (xSide,ySide,delT,delT_SqRoot,sdNoise,laplc3,gradSq,psi)
!==============================================================================
IMPLICIT NONE
COMPLEX(DP),INTENT(INOUT):: out1(xhalfp1,ySide)
COMPLEX(DP),INTENT(IN)   :: psik3(xhalfp1,ySide) 
INTEGER :: i,j,i1,j1
!REAL(DP),PARAMETER :: alpha2=225.D0, beta=0.9D0
!REAL(DP),PARAMETER :: alpha2=1.D0, beta=1.D00
!REAL(DP),PARAMETER :: r = -0.25D0
REAL(DP) :: kx, ky, k2, k4, delKX, delKY 
REAL(DP) :: term1, term2,term3 
REAL(DP):: sdNoise1,noise1
REAL(DP) ::Lx, Ly
COMPLEX(DP), ALLOCATABLE:: out1new(:,:)


ALLOCATE(out1new(xhalfp1,ySide))
sdNoise1 = 1.E-04_DP
!sdNoise1 = 10.0D0

Lx = xSide*dx
Ly = ySide*dy

delKX = (2.D0*pi)/(Lx)
delKY = (2.D0*pi)/(Ly)


!DO j = 1, ySide
!    DO i = 1, xSide
!      WRITE(100,*) i,j,dummyF(i,j)
!      IF(i==xSide) WRITE(100,*)''
!    ENDDO
!  ENDDO
!STOP
!WRITE(13,*) "iter", iter
  DO j = 1, ySide
    j1 = j - 1
    IF(j > (yHalfp1)) j1 = j - (ySide + 1)
    ky = j1*delKY
    DO i = 1,xHalfp1
      i1 = i -1
      kx = i1*delKX     
      k2 = kx*kx + ky*ky
      k4 = k2*k2
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!     Integrating factor method
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      term1 = -1.D0*(r + ((1.D0 - k2)*(1.D0 - k2)))*k2*dt 
      term2 = DEXP(term1) !+ (-1.D0*k2*psik3(i,j))
      CALL GasDev(sdNoise1,noise1)
      out1new(i,j) = term2*(out1(i,j) + (-1.D0*k2*psik3(i,j)*dt) + &
                            (-1.D0*DSQRT(k2)*noise1*DSQRT(dt)))
    ENDDO
  ENDDO
!  DO j = ySide/2+1 , ySide
!    jj = ySide - (j + 1)
!    DO i = 1, xSide
!      dummyF(i,j) =  dummyF(i,jj)
!    ENDDO
!  ENDDO
out1 = out1new

DEALLOCATE(out1new)
RETURN
!==============================================================================
                     END SUBROUTINE CalcEqOld
!==============================================================================
!==============================================================================
SUBROUTINE CalcEq(out1,psik3,out3,out4)
!           (xSide,ySide,delT,delT_SqRoot,sdNoise,laplc3,gradSq,psi)
!==============================================================================
IMPLICIT NONE
COMPLEX(DP),INTENT(INOUT):: out1(xhalfp1,ySide),out3(xhalfp1,ySide)
COMPLEX(DP),INTENT(INOUT):: out4(xhalfp1,ySide)
COMPLEX(DP),INTENT(IN)   :: psik3(xhalfp1,ySide) 
INTEGER :: i,j,i1,j1
!REAL(DP),PARAMETER :: alpha2=225.D0, beta=0.9D0
!REAL(DP),PARAMETER :: alpha2=1.D0, beta=1.D00
!REAL(DP),PARAMETER :: r = -0.25D0
REAL(DP) :: kx, ky, k2, k4, delKX, delKY 
REAL(DP) :: term1, term2,term3 
REAL(DP):: sdNoise1,noise1
REAL(DP) ::Lx, Ly, xr, xi
REAL(DP) :: alphaSq, beta
COMPLEX(DP), ALLOCATABLE:: out1new(:,:),out3new(:,:),out5(:,:)


ALLOCATE(out1new(xhalfp1,ySide))
ALLOCATE(out3new(xhalfp1,ySide))
ALLOCATE(out5(xhalfp1,ySide))
!sdNoise1 = 1.E-04_DP
sdNoise1  = 100.D0
!sdNoise1 = 1.225D0 !1.E-04_DP
beta = 1.D0 ! 1.D0
alphaSq = 1.D0 !1.D0
Lx = xSide*dx
Ly = ySide*dy

delKX = (2.D0*pi)/(Lx)
delKY = (2.D0*pi)/(Ly)

!out5 = (AIMAG(out4),DBLE(out4))
DO j = 1, ySide
  DO i = 1, xHalfp1
      xr =  AIMAG(out4(i,j))
      xi =  -1.D0*DBLE(out4(i,j))
      out5(i,j) = DCMPLX(xr,xi)
      !PRINT*, out4(i,j),out5(i,j)
ENDDO
ENDDO
  DO j = 1, ySide
    j1 = j - 1
    IF(j > (yHalfp1)) j1 = j - (ySide + 1)
    ky = j1*delKY
    DO i = 1,xHalfp1
      i1 = i -1
      kx = i1*delKX     
      k2 = kx*kx + ky*ky
      k4 = k2*k2
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!     Integrating factor method
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      term1 = DEXP(-1.D0*beta*dt)
      term2 = -1.D0*(r + ((1.D0 - k2)*(1.D0 - k2)))*k2*alphaSq 
      CALL GasDev(sdNoise1,noise1)
      out3new(i,j) = term1*(out3(i,j) + &
                   (((term2*out1(i,j)) +(-1.D0*k2*psik3(i,j)))*dt) + &
                   (1.D0*kx*out5(i,j)*dt)  + &
                   (-1.D0*DSQRT(k2)*noise1*DSQRT(dt)))
    ENDDO
  ENDDO
  out1new = out1 + out3new*dt
!  DO j = ySide/2+1 , ySide
!    jj = ySide - (j + 1)
!    DO i = 1, xSide
!      dummyF(i,j) =  dummyF(i,jj)
!    ENDDO
!  ENDDO
out1 = out1new
out3 = out3new
DEALLOCATE(out1new,out3new,out5)
RETURN
!==============================================================================
                     END SUBROUTINE CalcEq
!==============================================================================

!==============================================================================
                SUBROUTINE chari(i,ci,il)
!converts integer "i" into a character string, referenced as "cn(1:il)" commands
!==============================================================================
    integer ::i,ii,i3,kk,k,j,j1,il
        real    ::ri
        character(len=9) ci
        character(len=10) str
        ii=i
4       if(ii.gt.999999999) then
                !ri=ii! This IF part is used to rescale iteration
                !ii=nint(ri/10) !instead I chose to stop the code
                WRITE(*,*) "DEVISE A DIFFERENT METHOD TO NAME OUTPUT FILES"
                STOP
                !goto 4
        end if
        i3=ii
        str='0123456789' !String 
        do 11 k=1,9   !To find the number of digits
                j=10**k
                j1=10**(k-1)
                if((i3.ge.j1).and.(i3.lt.j)) il=k
11      continue
        do 22 k=il,1,-1 !Creating the string as same as ii
                kk=mod(ii,10)+1
                ci(k:k)=str(kk:kk)
                ii=ii/10
22      continue
        return
!==============================================================================
                 END SUBROUTINE chari
!==============================================================================

!==============================================================================
SUBROUTINE OutputPsi(dummyF,fileNo)
!==============================================================================
IMPLICIT NONE
REAL(DP),INTENT(IN):: dummyF(xSide,ySide)
INTEGER, INTENT(IN):: fileNo
INTEGER :: i,j


  DO j = 1, ySide
    DO i = 1, xSide
      WRITE(fileNo,*) i,j,dummyF(i,j)
      !WRITE(fileNo) (i,j,psi(i,j),i=1, xSide) !binary formatting
      !IF(j==ySide) WRITE(fileNo,*)''
      IF(i==(xSide)) WRITE(fileNo,*)''
    ENDDO
  ENDDO


RETURN
!==============================================================================
                        END SUBROUTINE OutputPsi
!==============================================================================

!*******************************************************************************
                           END MODULE modCalcPFC
!*******************************************************************************
