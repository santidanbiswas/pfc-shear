! test of fft4g2d.f
! ooura's code
!
!******************************************************************************
                       MODULE modFFT2d
!******************************************************************************

CONTAINS
!==============================================================================
       SUBROUTINE PassData2dFFT(nmax,nmaxsqrt,xSide,ySide,intSign,dummyF)
!==============================================================================

!
      INTEGER, INTENT(IN) :: nmax,nmaxsqrt, intSign
!      integer  nmax, nmaxsqrt
!      parameter (nmax = 1024)
!      parameter (nmaxsqrt = 32)
      INTEGER, INTENT(IN) :: xSide, ySide 
      integer ip(0 : nmaxsqrt + 1), n1, n2, i
      REAL*8,INTENT(INOUT) :: dummyF(0:xSide-1,0:ySide-1)
      real*8 a(0 : nmax - 1, 0 : nmax - 1), t(0 : 2 * nmax - 1), &
         w(0 : nmax * 3 / 2 - 1), err, errorcheck2d
      !REAL*8 :: scalefac
     !
      !write (*, *) 'data length n1=? (n1 = power of 2) '
      !read (*, *) n1
      !write (*, *) 'data length n2=? (n2 = power of 2) '
      !read (*, *) n2
      n1 = xSide
      n2 = ySide
     !scalefac = 2.D0/(n1*n2)
      ip(0) = 0
      pi = 4.D0*DATAN(1.D0)
!   check of RDFT
!      call putdata2d(nmax, n1, n2, a) !enter input
      !WRITE(40,*) "FT"
      !DO ii = 0,n1 -1
      !  WRITE(40,'(64(F6.2,1x))')(a(ii,jj),jj=0,n2-1)
        !WRITE(40,*)(a(ii,jj),jj=0,n2-1)
      !ENDDO

      DO j2 = 0,(n2-1)
        DO j1 = 0, (n1-1)
          a(j1,j2) = dummyF(j1,j2)
          !WRITE(40,*) j1, j2, a(j1,j2)
          !IF(j1==n2-1) WRITE(40,*)''
        ENDDO
      ENDDO

      call rdft2d(nmax, n1, n2,intSign, a, t, ip, w) !fourier transform
!       DO j2 = 0,(n2/2-1)
       !DO j2 = 0,(n2-1)
!        DO j1 = 0, (n1 - 1)
          !IF(MOD(j1,2)==0) WRITE(30,'(2(I3,2x),F6.2)') &
          !  j1,2*j2,a(j1,j2)
          !IF(MOD(j1,2)==0) WRITE(20,'(2(F6.2,2x),F6.2)') &
          !  j1*pi/N1,2.D0*j2*pi/N2,a(j1,j2)
!     c       j1,2*j2,a(j1,j2)
!          IF(MOD(j1,2)==1) WRITE(30,*)
          !IF((j1/2)*2==j1) WRITE(30,*)
!     c    pi*(j1-2)/N1,2.D0*pi*(j2-2)/N2,a(j1,j2)
           !IF(j1==n2-1) WRITE(20,*)''
!          IF(j1==n2-1) WRITE(30,*)''
!        ENDDO
!      ENDDO
!      DO j2 = 1,(n2-1)
!        j22 = j2 - 1
!        DO j1 = 1, (n1 - 1)
!          j11 = j1 - 1
!          dummyF(j1,j2) = a(j11,j22)
          !WRITE(40,*) j1, j2, a(j1,j2)
          !IF(j1==n2-1) WRITE(40,*)''
!        ENDDO
!      ENDDO

!      do j2 = 0, n2 - 1
!        do j1 = 0,n1-1
!          WRITE(50,*) j1,j2,a(j1,j2)*2.D0/(n1*n2)
!          IF(j1==n1-1) WRITE(50,*)''
!        enddo
!      enddo
      DO j2 = 0,(n2-1)
        DO j1 = 0, (n1 - 1)
          dummyF(j1,j2) = a(j1,j2)!*scalefac
          IF(MOD(j1,2)==0) WRITE(12,*) j1,j2, a(j1,j2),dummyF(j1,j2)
        ENDDO
      ENDDO
STOP
RETURN
!==============================================================================
                          END SUBROUTINE PassData2dFFT
!==============================================================================


!
! -------- FOR 2D Real DFT / Inverse of Real DFT --------
!     [definition]
!         <case1> RDFT
!             R(k1,k2) = sum_j1=0^n1-1 sum_j2=0^n2-1 a(j1,j2) * 
!                            cos(2*pi*j1*k1/n1 + 2*pi*j2*k2/n2), 
!                            0<=k1<n1, 0<=k2<n2
!             I(k1,k2) = sum_j1=0^n1-1 sum_j2=0^n2-1 a(j1,j2) * 
!                            sin(2*pi*j1*k1/n1 + 2*pi*j2*k2/n2), 
!                            0<=k1<n1, 0<=k2<n2
!         <case2> IRDFT (excluding scale)
!             a(k1,k2) = (1/2) * sum_j1=0^n1-1 sum_j2=0^n2-1
!                            (R(j1,j2) * 
!                            cos(2*pi*j1*k1/n1 + 2*pi*j2*k2/n2) + 
!                            I(j1,j2) * 
!                            sin(2*pi*j1*k1/n1 + 2*pi*j2*k2/n2)), 
!                            0<=k1<n1, 0<=k2<n2
!         (notes: R(n1-k1,n2-k2) = R(k1,k2), 
!                 I(n1-k1,n2-k2) = -I(k1,k2), 
!                 R(n1-k1,0) = R(k1,0), 
!                 I(n1-k1,0) = -I(k1,0), 
!                 R(0,n2-k2) = R(0,k2), 
!                 I(0,n2-k2) = -I(0,k2), 
!                 0<k1<n1, 0<k2<n2)
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call rdft2d(n1max, n1, n2, 1, a, t, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call rdft2d(n1max, n1, n2, -1, a, t, ip, w)
!     [parameters]
!         n1max  :row size of the 2D array (integer)
!         n1     :data length (integer)
!                 n1 >= 2, n1 = power of 2
!         n2     :data length (integer)
!                 n2 >= 2, n2 = power of 2
!         a(0:n1-1,0:n2-1)
!                :input/output data (real*8)
!                 <case1>
!                     output data
!                         a(2*k1,k2) = R(k1,k2) = R(n1-k1,n2-k2), 
!                         a(2*k1+1,k2) = I(k1,k2) = -I(n1-k1,n2-k2), 
!                             0<k1<n1/2, 0<k2<n2, 
!                         a(2*k1,0) = R(k1,0) = R(n1-k1,0), 
!                         a(2*k1+1,0) = I(k1,0) = -I(n1-k1,0), 
!                             0<k1<n1/2, 
!                         a(0,k2) = R(0,k2) = R(0,n2-k2), 
!                         a(1,k2) = I(0,k2) = -I(0,n2-k2), 
!                         a(1,n2-k2) = R(n1/2,k2) = R(n1/2,n2-k2), 
!                         a(0,n2-k2) = -I(n1/2,k2) = I(n1/2,n2-k2), 
!                             0<k2<n2/2, 
!                         a(0,0) = R(0,0), 
!                         a(1,0) = R(n1/2,0), 
!                         a(0,n2/2) = R(0,n2/2), 
!                         a(1,n2/2) = R(n1/2,n2/2)
!                 <case2>
!                     input data
!                         a(2*j1,j2) = R(j1,j2) = R(n1-j1,n2-j2), 
!                         a(2*j1+1,j2) = I(j1,j2) = -I(n1-j1,n2-j2), 
!                             0<j1<n1/2, 0<j2<n2, 
!                         a(2*j1,0) = R(j1,0) = R(n1-j1,0), 
!                         a(2*j1+1,0) = I(j1,0) = -I(n1-j1,0), 
!                             0<j1<n1/2, 
!                         a(0,j2) = R(0,j2) = R(0,n2-j2), 
!                         a(1,j2) = I(0,j2) = -I(0,n2-j2), 
!                         a(1,n2-j2) = R(n1/2,j2) = R(n1/2,n2-j2), 
!                         a(0,n2-j2) = -I(n1/2,j2) = I(n1/2,n2-j2), 
!                             0<j2<n2/2, 
!                         a(0,0) = R(0,0), 
!                         a(1,0) = R(n1/2,0), 
!                         a(0,n2/2) = R(0,n2/2), 
!                         a(1,n2/2) = R(n1/2,n2/2)
!         t(0:2*n2-1)
!                :work area (real*8)
!         ip(0:*):work area for bit reversal (integer)
!                 length of ip >= 2+sqrt(n)  ; if mod(n,4) = 0
!                                 2+sqrt(n/2); otherwise
!                 (n = max(n1/2, n2))
!                 ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:*) :cos/sin table (real*8)
!                 length of w >= max(n1/4, n2/2) + n1/4
!                 w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call rdft2d(n1max, n1, n2, 1, a, t, ip, w)
!         is 
!             call rdft2d(n1max, n1, n2, -1, a, t, ip, w)
!             do j2 = 0, n2 - 1
!                 do j1 = 0, n1 - 1
!                     a(j1, j2) = a(j1, j2) * (2.0d0 / (n1 * n2))
!                 end do
!             end do
!         .
!
!

!******************************************************************************
                 SUBROUTINE rdft2d(n1max, n1, n2, isgn, a, t, ip, w)
!Real discrete fourier transform in 2-dimensions
!******************************************************************************
      integer n1max, n1, n2, isgn, ip(0 : *), n, nw, nc, &
         n2h, i, j, j2
      real*8 a(0 : n1max - 1, 0 : n2 - 1), t(0 : 2 * n2 - 1), &
         w(0 : *), xi
      n = max(n1, 2 * n2)
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n1 .gt. 4 * nc) then
          nc = n1 / 4
          call makect(nc, ip, w(nw))
      end if
      n2h = n2 / 2
      if (isgn .lt. 0) then
          do i = 1, n2h - 1
              j = n2 - i
              xi = a(0, i) - a(0, j)
              a(0, i) = a(0, i) + a(0, j)
              a(0, j) = xi
              xi = a(1, j) - a(1, i)
              a(1, i) = a(1, i) + a(1, j)
              a(1, j) = xi
          end do
          do i = 0, n1 - 2, 2
              do j = 0, n2 - 1
                  j2 = 2 * j
                  t(j2) = a(i, j)
                  t(j2 + 1) = a(i + 1, j)
              end do
              call cdft(2 * n2, isgn, t, ip, w)
              do j = 0, n2 - 1
                  j2 = 2 * j
                  a(i, j) = t(j2)
                  a(i + 1, j) = t(j2 + 1)
              end do
          end do
          do j = 0, n2 - 1
              call rdft(n1, isgn, a(0, j), ip, w)
          end do
      else
          do j = 0, n2 - 1
              !WRITE(*,*) j,a(0,j)
              !WRITE(*,*) "j =", j, "a(:,:)"
              !DO ii = 0,n1 -1
              !   WRITE(*,'(4(F6.2,1x))')(a(ii,jj),jj=0,n2-1)
              !ENDDO

              call rdft(n1, isgn, a(0, j), ip, w)
              !WRITE(*,*) j,a(0,0),a(1,0),a(2,0),a(3,0)
          end do
          ! DO ii = 0,n1 -1
          !    WRITE(*,'(4(F6.2,1x))')(a(ii,jj),jj=0,n2-1)
          ! ENDDO

          !STOP
          do i = 0, n1 - 2, 2
              do j = 0, n2 - 1
                  j2 = 2 * j
                  t(j2) = a(i, j)
                  !WRITE(*,*) "i",i,"j",j,"j2",j2,a(i,j),a(i + 1, j)
                  t(j2 + 1) = a(i + 1, j)
                  !WRITE(*,*) 
              end do
              !WRITE(*,*) "t : =>"
              !WRITE(*,'(16(F6.2,1X))')(t(ii),ii=0,15)
              call cdft(2 * n2, isgn, t, ip, w)
              !WRITE(*,*) "t : =>"
              !WRITE(*,'(16(F6.2,1X))')(t(ii),ii=0,15)
              do j = 0, n2 - 1
                  j2 = 2 * j
                  a(i, j) = t(j2)
                  a(i + 1, j) = t(j2 + 1)
              end do
          end do
        !  DO ii = 0,n1 -1
        !    WRITE(*,'(4(F6.2,1x))')(a(ii,jj),jj=0,n2-1)
        !  ENDDO

          do i = 1, n2h - 1
              j = n2 - i
              a(0, j) = 0.5d0 * (a(0, i) - a(0, j))
              a(0, i) = a(0, i) - a(0, j)
              a(1, j) = 0.5d0 * (a(1, i) + a(1, j))
              a(1, i) = a(1, i) - a(1, j)
          end do
         ! DO ii = 0,n1 -1
         !   WRITE(*,'(4(F6.2,1x))')(a(ii,jj),jj=0,n2-1)
         ! ENDDO
      end if

!******************************************************************************
                      END SUBROUTINE rdft2d
!******************************************************************************

! -------- Complex DFT (Discrete Fourier Transform) --------
!     [definition]
!         <case1>
!             X(k) = sum_j=0^n-1 x(j)*exp(2*pi*i*j*k/n), 0<=k<n
!         <case2>
!             X(k) = sum_j=0^n-1 x(j)*exp(-2*pi*i*j*k/n), 0<=k<n
!         (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, -1, a, ip, w)
!     [parameters]
!         2*n          :data length (integer)
!                       n >= 1, n = power of 2
!         a(0:2*n-1)   :input/output data (real*8)
!                       input data
!                           a(2*j) = Re(x(j)), 
!                           a(2*j+1) = Im(x(j)), 0<=j<n
!                       output data
!                           a(2*k) = Re(X(k)), 
!                           a(2*k+1) = Im(X(k)), 0<=k<n
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n)  ; if mod(n,4) = 0
!                                       2+sqrt(n/2); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call cdft(2*n, -1, a, ip, w)
!         is 
!             call cdft(2*n, 1, a, ip, w)
!             do j = 0, 2 * n - 1
!                 a(j) = a(j) / n
!             end do
!         .
!
!
      subroutine cdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *)
      real*8 a(0 : n - 1), w(0 : *)
      if (n .gt. 4 * ip(0)) then
          call makewt(n / 4, ip, w)
      end if
          !WRITE(*,*)"bitrv2 called from cdft"
      if (n .gt. 4) call bitrv2(n, ip(2), a)
      if (isgn .lt. 0) then
          call cftfsub(n, a, w)
      else
          call cftbsub(n, a, w)
      end if
      end subroutine cdft

! -------- Real DFT / Inverse of Real DFT --------
!     [definition]
!         <case1> RDFT
!             R(k) = sum_j=0^n-1 a(j)*cos(2*pi*j*k/n), 0<=k<=n/2
!             I(k) = sum_j=0^n-1 a(j)*sin(2*pi*j*k/n), 0<k<n/2
!         <case2> IRDFT (excluding scale)
!             a(k) = R(0)/2 + R(n/2)/2 + 
!                    sum_j=1^n/2-1 R(j)*cos(2*pi*j*k/n) + 
!                    sum_j=1^n/2-1 I(j)*sin(2*pi*j*k/n), 0<=k<n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call rdft(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call rdft(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       <case1>
!                           output data
!                               a(2*k) = R(k), 0<=k<n/2
!                               a(2*k+1) = I(k), 0<k<n/2
!                               a(1) = R(n/2)
!                       <case2>
!                           input data
!                               a(2*j) = R(j), 0<=j<n/2
!                               a(2*j+1) = I(j), 0<j<n/2
!                               a(1) = R(n/2)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2); if mod(n,4) = 2
!                                       2+sqrt(n/4); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call rdft(n, 1, a, ip, w)
!         is 
!             call rdft(n, -1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         
!
      subroutine rdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), nw, nc
      real*8 a(0 : n - 1), w(0 : *), xi
      !WRITE(*,*)(a(ii),ii=0,n-1)
      !STOP
      nw = ip(0)
      if (n .gt. 4 * nw) then
          !WRITE(*,*)" 1leaves rdft"
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 4 * nc) then
          !WRITE(*,*)" 2leaves rdft"
          nc = n / 4
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          !WRITE(*,*)" 3leaves rdft"
          a(1) = 0.5d0 * (a(0) - a(1))
          a(0) = a(0) - a(1)
          if (n .gt. 4) then
              call rftfsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
          end if
          call cftfsub(n, a, w)
      else
          !WRITE(*,*)"  rdft", a(0),a(1)
          if (n .gt. 4) call bitrv2(n, ip(2), a)
          call cftbsub(n, a, w)
          !WRITE(*,*)"  rdft2", a(0),a(1),a(2),a(3)
          if (n .gt. 4) call rftbsub(n, a, nc, w(nw))
          xi = a(0) - a(1)
          a(0) = a(0) + a(1)
          a(1) = xi
          !WRITE(*,*)" leaves rdft", a(0),a(1)
      end if
      end subroutine rdft

!
! -------- initializing routines --------
!
!******************************************************************************
                    SUBROUTINE makewt(nw, ip, w)
!******************************************************************************
      integer nw, ip(0 : *), nwh, j
      real*8 w(0 : nw - 1), delta, x, y

      ip(0) = nw
      ip(1) = 1
      if (nw .gt. 2) then
          !WRITE(*,*)"makewt 1"
          nwh = nw / 2
          delta = atan(1.0d0) / nwh
          w(0) = 1
          w(1) = 0
          w(nwh) = cos(delta * nwh)
          w(nwh + 1) = w(nwh)
          do j = 2, nwh - 2, 2
              x = cos(delta * j)
              y = sin(delta * j)
              w(j) = x
              w(j + 1) = y
              w(nw - j) = y
              w(nw - j + 1) = x
          end do
          !WRITE(*,*)"bitrv2 called from makewt"
          call bitrv2(nw, ip(2), w)
      end if
!******************************************************************************
                     END SUBROUTINE makewt
!******************************************************************************
!
      subroutine makect(nc, ip, c)
      integer nc, ip(0 : *), nch, j
      real*8 c(0 : nc - 1), delta
      ip(1) = nc
      if (nc .gt. 1) then
          nch = nc / 2
          delta = atan(1.0d0) / nch
          c(0) = 0.5d0
          c(nch) = 0.5d0 * cos(delta * nch)
          do j = 1, nch - 1
              c(j) = 0.5d0 * cos(delta * j)
              c(nc - j) = 0.5d0 * sin(delta * j)
          end do
      end if
      end subroutine makect
!
! -------- child routines --------
!
      subroutine bitrv2(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, m2
      real*8 a(0 : n - 1), xr, xi
      !WRITE(*,*) "Subroutine bitrv2 entered"
      ip(0) = 0
      l = n
      m = 1
      do while (4 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      if (4 * m .gt. l) then
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  a(j1) = a(k1)
                  a(j1 + 1) = a(k1 + 1)
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
          end do
      else
          m2 = 2 * m
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  a(j1) = a(k1)
                  a(j1 + 1) = a(k1 + 1)
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  a(j1) = a(k1)
                  a(j1 + 1) = a(k1 + 1)
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
          end do
      end if
      end subroutine bitrv2
!
      subroutine cftbsub(n, a, w)
      integer n, j, j1, j2, j3, k, k1, ks, l, m
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      l = 2

      !WRITE(*,'("cftbsub",3(F6.2,1x))') (a(ii),ii=0,n-1)
      do while (2 * l .lt. n)
          m = 4 * l
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              !WRITE(*,'("cftbsub",3(F6.2,1x))') (a(ii),ii=0,n-1)
              !STOP
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i - x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i + x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i - x3r
          end do
          if (m .lt. n) then
              wk1r = w(2)
              do j = m, l + m - 2, 2
                  j1 = j + l
                  j2 = j1 + l
                  j3 = j2 + l
                  x0r = a(j) + a(j1)
                  x0i = a(j + 1) + a(j1 + 1)
                  x1r = a(j) - a(j1)
                  x1i = a(j + 1) - a(j1 + 1)
                  x2r = a(j2) + a(j3)
                  x2i = a(j2 + 1) + a(j3 + 1)
                  x3r = a(j2) - a(j3)
                  x3i = a(j2 + 1) - a(j3 + 1)
                  a(j) = x0r + x2r
                  a(j + 1) = x0i + x2i
                  a(j2) = x2i - x0i
                  a(j2 + 1) = x0r - x2r
                  x0r = x1r - x3i
                  x0i = x1i + x3r
                  a(j1) = wk1r * (x0r - x0i)
                  a(j1 + 1) = wk1r * (x0r + x0i)
                  x0r = x3i + x1r
                  x0i = x3r - x1i
                  a(j3) = wk1r * (x0i - x0r)
                  a(j3 + 1) = wk1r * (x0i + x0r)
              end do
              k1 = 1
              ks = -1
              do k = 2 * m, n - m, m
                  k1 = k1 + 1
                  ks = -ks
                  wk1r = w(2 * k1)
                  wk1i = w(2 * k1 + 1)
                  wk2r = ks * w(k1)
                  wk2i = w(k1 + ks)
                  wk3r = wk1r - 2 * wk2i * wk1i
                  wk3i = 2 * wk2i * wk1r - wk1i
                  do j = k, l + k - 2, 2
                      j1 = j + l
                      j2 = j1 + l
                      j3 = j2 + l
                      x0r = a(j) + a(j1)
                      x0i = a(j + 1) + a(j1 + 1)
                      x1r = a(j) - a(j1)
                      x1i = a(j + 1) - a(j1 + 1)
                      x2r = a(j2) + a(j3)
                      x2i = a(j2 + 1) + a(j3 + 1)
                      x3r = a(j2) - a(j3)
                      x3i = a(j2 + 1) - a(j3 + 1)
                      a(j) = x0r + x2r
                      a(j + 1) = x0i + x2i
                      x0r = x0r - x2r
                      x0i = x0i - x2i
                      a(j2) = wk2r * x0r - wk2i * x0i
                      a(j2 + 1) = wk2r * x0i + wk2i * x0r
                      x0r = x1r - x3i
                      x0i = x1i + x3r
                      a(j1) = wk1r * x0r - wk1i * x0i
                      a(j1 + 1) = wk1r * x0i + wk1i * x0r
                      x0r = x1r + x3i
                      x0i = x1i - x3r
                      a(j3) = wk3r * x0r - wk3i * x0i
                      a(j3 + 1) = wk3r * x0i + wk3i * x0r
                  end do
              end do
          end if
          l = m
      end do
      if (l .lt. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = a(j + 1) - a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = a(j + 1) + a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end subroutine cftbsub
!
      subroutine cftfsub(n, a, w)
      integer n, j, j1, j2, j3, k, k1, ks, l, m
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      l = 2
      do while (2 * l .lt. n)
          m = 4 * l
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i - x2i
              a(j1) = x1r + x3i
              a(j1 + 1) = x1i - x3r
              a(j3) = x1r - x3i
              a(j3 + 1) = x1i + x3r
          end do
          if (m .lt. n) then
              wk1r = w(2)
              do j = m, l + m - 2, 2
                  j1 = j + l
                  j2 = j1 + l
                  j3 = j2 + l
                  x0r = a(j) + a(j1)
                  x0i = a(j + 1) + a(j1 + 1)
                  x1r = a(j) - a(j1)
                  x1i = a(j + 1) - a(j1 + 1)
                  x2r = a(j2) + a(j3)
                  x2i = a(j2 + 1) + a(j3 + 1)
                  x3r = a(j2) - a(j3)
                  x3i = a(j2 + 1) - a(j3 + 1)
                  a(j) = x0r + x2r
                  a(j + 1) = x0i + x2i
                  a(j2) = x0i - x2i
                  a(j2 + 1) = x2r - x0r
                  x0r = x1r + x3i
                  x0i = x1i - x3r
                  a(j1) = wk1r * (x0i + x0r)
                  a(j1 + 1) = wk1r * (x0i - x0r)
                  x0r = x3i - x1r
                  x0i = x3r + x1i
                  a(j3) = wk1r * (x0r + x0i)
                  a(j3 + 1) = wk1r * (x0r - x0i)
              end do
              k1 = 1
              ks = -1
              do k = 2 * m, n - m, m
                  k1 = k1 + 1
                  ks = -ks
                  wk1r = w(2 * k1)
                  wk1i = w(2 * k1 + 1)
                  wk2r = ks * w(k1)
                  wk2i = w(k1 + ks)
                  wk3r = wk1r - 2 * wk2i * wk1i
                  wk3i = 2 * wk2i * wk1r - wk1i
                  do j = k, l + k - 2, 2
                      j1 = j + l
                      j2 = j1 + l
                      j3 = j2 + l
                      x0r = a(j) + a(j1)
                      x0i = a(j + 1) + a(j1 + 1)
                      x1r = a(j) - a(j1)
                      x1i = a(j + 1) - a(j1 + 1)
                      x2r = a(j2) + a(j3)
                      x2i = a(j2 + 1) + a(j3 + 1)
                      x3r = a(j2) - a(j3)
                      x3i = a(j2 + 1) - a(j3 + 1)
                      a(j) = x0r + x2r
                      a(j + 1) = x0i + x2i
                      x0r = x0r - x2r
                      x0i = x0i - x2i
                      a(j2) = wk2r * x0r + wk2i * x0i
                      a(j2 + 1) = wk2r * x0i - wk2i * x0r
                      x0r = x1r + x3i
                      x0i = x1i - x3r
                      a(j1) = wk1r * x0r + wk1i * x0i
                      a(j1 + 1) = wk1r * x0i - wk1i * x0r
                      x0r = x1r - x3i
                      x0i = x1i + x3r
                      a(j3) = wk3r * x0r + wk3i * x0i
                      a(j3 + 1) = wk3r * x0i - wk3i * x0r
                  end do
              end do
          end if
          l = m
      end do
      if (l .lt. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = a(j + 1) - a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = a(j + 1) + a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end subroutine cftfsub
!
      subroutine rftbsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      ks = 4 * nc / n
      kk = 0
      do k = n / 2 - 2, 2, -2
          j = n - k
          kk = kk + ks
          wkr = 0.5d0 - c(kk)
          wki = c(nc - kk)
          xr = a(k) - a(j)
          xi = a(k + 1) + a(j + 1)
          yr = wkr * xr - wki * xi
          yi = wkr * xi + wki * xr
          a(k) = a(k) - yr
          a(k + 1) = a(k + 1) - yi
          a(j) = a(j) + yr
          a(j + 1) = a(j + 1) - yi
      end do
      end subroutine rftbsub
!
      subroutine rftfsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      ks = 4 * nc / n
      kk = 0
      do k = n / 2 - 2, 2, -2
          j = n - k
          kk = kk + ks
          wkr = 0.5d0 - c(kk)
          wki = c(nc - kk)
          xr = a(k) - a(j)
          xi = a(k + 1) + a(j + 1)
          yr = wkr * xr + wki * xi
          yi = wkr * xi - wki * xr
          a(k) = a(k) - yr
          a(k + 1) = a(k + 1) - yi
          a(j) = a(j) + yr
          a(j + 1) = a(j + 1) - yi
      end do
      end  subroutine rftfsub
!
!******************************************************************************
                     END  MODULE modFFT2d
!******************************************************************************

