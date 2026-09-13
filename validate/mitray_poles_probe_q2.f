C     Standalone reference-value generator for the Q2 quadrupole.
C
C     Calls the REAL, unmodified GEANT3 routines src/mitray_zone.f and
C     src/mitray_poles.f (the DRAGON simulation's own MIT-RAYTRACE field
C     code) at a grid of test points, using the Q2 quadrupole's actual
C     parameters from dat/dragon_2014_DSSSD.dat, and writes the resulting
C     B-field to a text file. Unlike Q1, Q2 has a nonzero hexapole (BHX)
C     component, so this also exercises the multipole terms beyond the
C     pure quadrupole one.
C
      PROGRAM mitray_poles_probe_q2
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
C     'POLE' 'Q2  ' card, dat/dragon_2014_DSSSD.dat, laid out per the
C     field indices read in src/mitray_setup.f
C
      DATA(1)  = 3.D0
      DATA(2)  = 3.D0
      DATA(3)  = 3.D0
      DATA(10) = 0.D0
      DATA(11) = 0.D0
      DATA(12) = 33.385D0
      DATA(13) = 7.9375D0
      DATA(14) = 0.08637239D0
      DATA(15) = 0.0045675D0
      DATA(16) = 0.D0
      DATA(17) = 0.D0
      DATA(18) = 0.D0
      DATA(19) = 23.813D0
      DATA(20) = -19.844D0
      DATA(21) = -19.844D0
      DATA(22) = 23.813D0
      DATA(23) = 0.22D0
      DATA(24) = 5.367112D0
      DATA(25) = -1.99912D0
      DATA(26) = 0.911917D0
      DATA(27) = -0.663814D0
      DATA(28) = 0.348883D0
      DATA(29) = 0.22D0
      DATA(30) = 5.367112D0
      DATA(31) = -1.99912D0
      DATA(32) = 0.911917D0
      DATA(33) = -0.663814D0
      DATA(34) = 0.348883D0
      DATA(35) = 0.D0
      DATA(36) = 0.D0
      DATA(37) = 0.D0
      DATA(38) = 0.D0
      DATA(39) = 0.D0
      DATA(40) = 0.D0
      DATA(41) = 0.D0
      DATA(42) = 0.D0
C
      OPEN(UNIT=30, FILE='q2_fortran_reference.csv', STATUS='UNKNOWN')
      WRITE(30,'(A)') 'x_cm,y_cm,z_cm,bx_T,by_T,bz_T'
C
C     Grid: covers the (larger than Q1's) aperture in x/y, and entrance
C     fringe / uniform / exit fringe / far-field zones in z.
C
      nx = 5
      ny = 5
      nz = 21
      xlo = -4.5D0
      xhi =  4.5D0
      ylo = -4.5D0
      yhi =  4.5D0
      zlo = -35.0D0
      zhi =  35.0D0
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
            CALL MITRAY_POLES(DATA, XPOS, BFLD)
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
      WRITE(6,*) 'Wrote q2_fortran_reference.csv'
C
      END
