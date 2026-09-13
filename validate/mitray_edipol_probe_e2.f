C     Standalone reference-value generator for the E2 electrostatic
C     deflector.
C
C     Calls the REAL, unmodified GEANT3 routine src/mitray_edipol.f (the
C     DRAGON simulation's own MIT-RAYTRACE electric-dipole field code) at
C     a grid of test points, using E2's actual parameters from
C     dat/dragon_2014_DSSSD.dat, and writes the resulting E-field to a
C     text file. This is the ground truth that the C++ port
C     (MitrayEdipoleField) is checked against.
C
      PROGRAM mitray_edipol_probe_e2
C
      IMPLICIT none
C
      include 'gcflag.inc'
      include 'gcunit.inc'
      include 'mitray_diag.inc'
C
      REAL*8 DATA(75), XPOS(3), EFLD(3)
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
C     'EDIP' 'E2 ' card, dat/dragon_2014_DSSSD.dat
C
      DATA(1)  = 1.D0
      DATA(2)  = 100.D0
      DATA(3)  = 100.D0
      DATA(4)  = 0.6D0
      DATA(11) = 0.D0
      DATA(12) = 0.D0
      DATA(13) = 10.D0
      DATA(14) = 250.D0
      DATA(15) = -3.77557049D0
      DATA(16) = 35.D0
      DATA(17) = 0.D0
      DATA(18) = 0.D0
      DATA(19) = 50.D0
      DATA(20) = 0.D0
      DATA(25) = 15.1658D0
      DATA(26) = -14.6504D0
      DATA(27) = -14.6504D0
      DATA(28) = 15.1658D0
      DATA(29) = 0.07901D0
      DATA(30) = 3.90918D0
      DATA(31) = -0.65329D0
      DATA(32) = 1.91401D0
      DATA(33) = 0.22838D0
      DATA(34) = -0.80791D0
      DATA(35) = 0.07901D0
      DATA(36) = 3.90918D0
      DATA(37) = -0.65329D0
      DATA(38) = 1.91401D0
      DATA(39) = 0.22838D0
      DATA(40) = -0.80791D0
C
      OPEN(UNIT=30, FILE='e2_fortran_reference.csv', STATUS='UNKNOWN')
      WRITE(30,'(A)') 'x_cm,y_cm,z_cm,ex_kVcm,ey_kVcm,ez_kVcm'
C
      nx = 7
      ny = 5
      nz = 25
      xlo = -4.0D0
      xhi =  4.0D0
      ylo = -2.0D0
      yhi =  2.0D0
      zlo = -30.0D0
      zhi =  30.0D0
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
            CALL MITRAY_EDIPOL(DATA, XPOS, EFLD)
C
            WRITE(30,100) x, y, z, EFLD(1), EFLD(2), EFLD(3)
  100       FORMAT(F10.4,',',F10.4,',',F10.4,',',
     +             E16.8,',',E16.8,',',E16.8)
C
          ENDDO
        ENDDO
      ENDDO
C
      CLOSE(30)
C
      WRITE(6,*) 'Wrote e2_fortran_reference.csv'
C
      END
