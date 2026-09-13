C     Standalone reference-value generator for the D2 dipole.
C
C     Calls the REAL, unmodified GEANT3 routine src/mitray_dipole.f (the
C     DRAGON simulation's own MIT-RAYTRACE dipole field code) at a grid of
C     test points, using the D2 dipole's actual parameters from
C     dat/dragon_2014_DSSSD.dat, and writes the resulting B-field to a
C     text file. This is the ground truth that the C++ port
C     (MitrayDipoleField) is checked against.
C
      PROGRAM mitray_dipole_probe_d2
C
      IMPLICIT none
C
      include 'gcflag.inc'
      include 'gcunit.inc'
      include 'mitray_diag.inc'
      include 'diagnostic.inc'
C
      REAL*8 DATA(75), XPOS(3), BFLD(3)
      INTEGER i, ix, iy, iz, nx, ny, nz
      REAL*8 x, y, z, xlo, xhi, ylo, yhi, zlo, zhi
C
      lout = 6
      ldiag = .false.
C
      DO i = 1, 75
        DATA(i) = 0.D0
      ENDDO
C
C     'DIPO' 'D2  ' card, dat/dragon_2014_DSSSD.dat
C
      DATA(1)  = 1.D0
      DATA(2)  = 1.D0
      DATA(3)  = 1.D0
      DATA(4)  = 1.D0
      DATA(5)  = 3.D0
      DATA(6)  = 0.D0
      DATA(11) = 0.D0
      DATA(12) = 0.D0
      DATA(13) = 12.D0
      DATA(14) = 81.3D0
      DATA(15) = 0.26520546D0
      DATA(16) = 75.D0
      DATA(17) = 29.D0
      DATA(18) = 29.D0
      DATA(19) = 0.D0
      DATA(20) = 0.D0
      DATA(21) = 0.D0
      DATA(22) = 0.D0
      DATA(25) = 36.D0
      DATA(26) = -30.D0
      DATA(27) = -30.D0
      DATA(28) = 36.D0
      DATA(29) = 0.3295D0
      DATA(30) = 3.31886D0
      DATA(31) = -1.2036D0
      DATA(32) = 0.181157D0
      DATA(33) = 0.1103868D0
      DATA(34) = -0.029513D0
      DATA(35) = 0.3295D0
      DATA(36) = 3.31886D0
      DATA(37) = -1.2036D0
      DATA(38) = 0.181157D0
      DATA(39) = 0.1103868D0
      DATA(40) = -0.029513D0
      DATA(41) = 0.D0
      DATA(42) = 0.D0
      DATA(43) = 0.D0
      DATA(44) = 0.D0
      DATA(45) = 0.D0
      DATA(46) = 0.D0
      DATA(47) = 0.D0
      DATA(48) = 0.D0
      DATA(49) = 80.D0
      DATA(50) = 80.D0
C     DATA(51..64) already zero -> flat pole faces (NSRF=0)
C
      OPEN(UNIT=30, FILE='d2_fortran_reference.csv', STATUS='UNKNOWN')
      WRITE(30,'(A)') 'x_cm,y_cm,z_cm,bx_T,by_T,bz_T'
C
C     Grid: covers the aperture in x/y (small, since Y is the vertical
C     gap direction), and a wide sweep in z spanning entrance fringe,
C     uniform (bend) region, and exit fringe.
C
      nx = 7
      ny = 5
      nz = 25
C     Kept well inside XBMAX=WDIP1/2=40cm even after the ALPHA=29deg
C     A->B rotation (unlike D1's ALPHA=5.8deg, this is not small), so no
C     grid point falls into the "unknown region" catch-all.
      xlo = -15.0D0
      xhi =  15.0D0
      ylo = -3.0D0
      yhi =  3.0D0
      zlo = -55.0D0
      zhi =  55.0D0
C
      DO iz = 0, nz-1
        z = zlo + (zhi-zlo)*DBLE(iz)/DBLE(nz-1)
        DO ix = 0, nx-1
          x = xlo + (xhi-xlo)*DBLE(ix)/DBLE(nx-1)
          DO iy = 0, ny-1
            y = ylo + (yhi-ylo)*DBLE(iy)/DBLE(ny-1)
C
            XPOS(1) = x
            XPOS(2) = y
            XPOS(3) = z
C
            CALL MITRAY_DIPOLE(DATA, XPOS, BFLD)
C
            WRITE(30,100) x, y, z, BFLD(1), BFLD(2), BFLD(3)
  100       FORMAT(F10.4,',',F10.4,',',F10.4,',',
     +             E16.8,',',E16.8,',',E16.8)
C
          ENDDO
        ENDDO
      ENDDO
C
      CLOSE(30)
C
      WRITE(6,*) 'Wrote d2_fortran_reference.csv'
C
      END
