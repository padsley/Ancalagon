C     Standalone reference-value generator for Q3 through Q14 (every
C     'POLE' card in dat/dragon_2014_DSSSD.dat besides Q1/Q2, which already
C     have their own probes).
C
C     Calls the REAL, unmodified GEANT3 routines src/mitray_zone.f and
C     src/mitray_poles.f at a small grid of test points per quadrupole,
C     using each card's actual parameters, and writes the resulting
C     B-field to one combined text file (with a name column) -- this is
C     the ground truth the C++ port (MitrayPoleData::Q3()..Q14()) is
C     checked against.
C
      PROGRAM mitray_poles_probe_rest
C
      IMPLICIT none
C
      include 'gcflag.inc'
      include 'gcunit.inc'
      include 'mitray_diag.inc'
      include 'diagnostic.inc'
C
      REAL*8 DATA(75), XPOS(3), BFLD(3)
      INTEGER i, ix, iy, iz, nx, ny, nz, iq
      REAL*8 x, y, z, xlo, xhi, ylo, yhi, zlo, zhi
      CHARACTER*4 name(12)
      REAL*8 lval(12), rad(12), zfar(12), znear(12)
      REAL*8 bqd(12), bhx(12), boc(12), bdc(12), bdd(12)
      REAL*8 c0(12), c1(12), c2(12), c3(12), c4(12), c5(12)
C
      lout = 6
      ldiag = .false.
C
C     Order: Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14
C
      name(1)='Q3  '
      lval(1)=18.75D0
      rad(1)=7.95D0
      zfar(1)=20.D0
      znear(1)=-20.0D0
      bqd(1)=0.D0
      bhx(1)=0.01831480D0
      boc(1)=0.D0
      bdc(1)=0.D0
      bdd(1)=0.D0
      c0(1)=0.D0
      c1(1)=0.D0
      c2(1)=0.D0
      c3(1)=0.D0
      c4(1)=0.D0
      c5(1)=0.D0
C
      name(2)='Q4  '
      lval(2)=33.38D0
      rad(2)=7.9375D0
      zfar(2)=23.813D0
      znear(2)=-19.844D0
      bqd(2)=0.07875767D0
      bhx(2)=0.D0
      boc(2)=0.D0
      bdc(2)=0.D0
      bdd(2)=0.D0
      c0(2)=0.225D0
      c1(2)=6.22466D0
      c2(2)=-2.38148D0
      c3(2)=0.2341D0
      c4(2)=-0.72032D0
      c5(2)=0.72371D0
C
      name(3)='Q5  '
      lval(3)=33.38D0
      rad(3)=7.9375D0
      zfar(3)=23.813D0
      znear(3)=-19.844D0
      bqd(3)=-0.10400018D0
      bhx(3)=0.D0
      boc(3)=0.D0
      bdc(3)=0.D0
      bdd(3)=0.D0
      c0(3)=0.225D0
      c1(3)=6.22466D0
      c2(3)=-2.38148D0
      c3(3)=0.2341D0
      c4(3)=-0.72032D0
      c5(3)=0.72371D0
C
      name(4)='Q6  '
      lval(4)=33.38D0
      rad(4)=7.9375D0
      zfar(4)=23.813D0
      znear(4)=-19.844D0
      bqd(4)=0.05728922D0
      bhx(4)=0.D0
      boc(4)=0.D0
      bdc(4)=0.D0
      bdd(4)=0.D0
      c0(4)=0.225D0
      c1(4)=6.22466D0
      c2(4)=-2.38148D0
      c3(4)=0.2341D0
      c4(4)=-0.72032D0
      c5(4)=0.72371D0
C
      name(5)='Q7  '
      lval(5)=18.75D0
      rad(5)=7.95D0
      zfar(5)=20.D0
      znear(5)=-20.0D0
      bqd(5)=0.D0
      bhx(5)=0.00385606D0
      boc(5)=0.D0
      bdc(5)=0.D0
      bdd(5)=0.D0
      c0(5)=0.D0
      c1(5)=0.D0
      c2(5)=0.D0
      c3(5)=0.D0
      c4(5)=0.D0
      c5(5)=0.D0
C
      name(6)='Q8  '
      lval(6)=25.23D0
      rad(6)=5.3975D0
      zfar(6)=18.89D0
      znear(6)=-13.494D0
      bqd(6)=-0.05091622D0
      bhx(6)=0.D0
      boc(6)=0.D0
      bdc(6)=0.D0
      bdd(6)=0.D0
      c0(6)=0.295D0
      c1(6)=6.30221D0
      c2(6)=-3.51059D0
      c3(6)=0.29528D0
      c4(6)=1.19866D0
      c5(6)=-0.423408D0
C
      name(7)='Q9  '
      lval(7)=33.38D0
      rad(7)=7.9375D0
      zfar(7)=23.813D0
      znear(7)=-19.844D0
      bqd(7)=0.0731129D0
      bhx(7)=0.D0
      boc(7)=0.D0
      bdc(7)=0.D0
      bdd(7)=0.D0
      c0(7)=0.225D0
      c1(7)=6.22466D0
      c2(7)=-2.38148D0
      c3(7)=0.2341D0
      c4(7)=-0.72032D0
      c5(7)=0.72371D0
C
      name(8)='Q10 '
      lval(8)=19.9D0
      rad(8)=8.0D0
      zfar(8)=20.D0
      znear(8)=-20.0D0
      bqd(8)=0.D0
      bhx(8)=0.0020298D0
      boc(8)=0.D0
      bdc(8)=0.0007228D0
      bdd(8)=0.D0
      c0(8)=0.D0
      c1(8)=0.D0
      c2(8)=0.D0
      c3(8)=0.D0
      c4(8)=0.D0
      c5(8)=0.D0
C
      name(9)='Q11 '
      lval(9)=33.38D0
      rad(9)=7.9375D0
      zfar(9)=23.813D0
      znear(9)=-19.844D0
      bqd(9)=0.05420925D0
      bhx(9)=0.D0
      boc(9)=0.D0
      bdc(9)=0.D0
      bdd(9)=0.D0
      c0(9)=0.225D0
      c1(9)=6.22466D0
      c2(9)=-2.38148D0
      c3(9)=0.2341D0
      c4(9)=-0.72032D0
      c5(9)=0.72371D0
C
      name(10)='Q12 '
      lval(10)=19.9D0
      rad(10)=8.0D0
      zfar(10)=20.D0
      znear(10)=-20.0D0
      bqd(10)=0.D0
      bhx(10)=0.015526D0
      boc(10)=0.D0
      bdc(10)=0.005701D0
      bdd(10)=0.D0
      c0(10)=0.D0
      c1(10)=0.D0
      c2(10)=0.D0
      c3(10)=0.D0
      c4(10)=0.D0
      c5(10)=0.D0
C
      name(11)='Q13 '
      lval(11)=46.7D0
      rad(11)=7.5D0
      zfar(11)=20.25D0
      znear(11)=-18.75D0
      bqd(11)=-0.04192126D0
      bhx(11)=0.D0
      boc(11)=0.D0
      bdc(11)=0.D0
      bdd(11)=0.D0
      c0(11)=0.2535D0
      c1(11)=5.840314D0
      c2(11)=-3.40247D0
      c3(11)=1.456423D0
      c4(11)=1.44575D0
      c5(11)=-0.754832D0
C
      name(12)='Q14 '
      lval(12)=46.7D0
      rad(12)=7.5D0
      zfar(12)=20.25D0
      znear(12)=-18.75D0
      bqd(12)=0.04687356D0
      bhx(12)=0.D0
      boc(12)=0.D0
      bdc(12)=0.D0
      bdd(12)=0.D0
      c0(12)=0.2535D0
      c1(12)=5.840314D0
      c2(12)=-3.40247D0
      c3(12)=1.456423D0
      c4(12)=1.44575D0
      c5(12)=-0.754832D0
C
      OPEN(UNIT=30, FILE='q3_q14_fortran_reference.csv',
     +     STATUS='UNKNOWN')
      WRITE(30,'(A)') 'name,x_cm,y_cm,z_cm,bx_T,by_T,bz_T'
C
      nx = 5
      ny = 5
      nz = 11
C
      DO iq = 1, 12
C
        DO i = 1, 75
          DATA(i) = 0.D0
        ENDDO
C
        DATA(1)  = 3.D0
        DATA(2)  = 3.D0
        DATA(3)  = 3.D0
        DATA(10) = 0.D0
        DATA(11) = 0.D0
        DATA(12) = lval(iq)
        DATA(13) = rad(iq)
        DATA(14) = bqd(iq)
        DATA(15) = bhx(iq)
        DATA(16) = boc(iq)
        DATA(17) = bdc(iq)
        DATA(18) = bdd(iq)
        DATA(19) = zfar(iq)
        DATA(20) = znear(iq)
        DATA(21) = znear(iq)
        DATA(22) = zfar(iq)
        DATA(23) = c0(iq)
        DATA(24) = c1(iq)
        DATA(25) = c2(iq)
        DATA(26) = c3(iq)
        DATA(27) = c4(iq)
        DATA(28) = c5(iq)
        DATA(29) = c0(iq)
        DATA(30) = c1(iq)
        DATA(31) = c2(iq)
        DATA(32) = c3(iq)
        DATA(33) = c4(iq)
        DATA(34) = c5(iq)
C
        xlo = -0.5D0*rad(iq)
        xhi =  0.5D0*rad(iq)
        ylo = -0.5D0*rad(iq)
        yhi =  0.5D0*rad(iq)
        zlo = -1.3D0*lval(iq)/2.D0 - 5.D0
        zhi =  1.3D0*lval(iq)/2.D0 + 5.D0
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
              WRITE(30,100) name(iq), x, y, z,
     +                      BFLD(1), BFLD(2), BFLD(3)
  100         FORMAT(A4,',',F10.4,',',F10.4,',',F10.4,',',
     +               E16.8,',',E16.8,',',E16.8)
C
            ENDDO
          ENDDO
        ENDDO
C
      ENDDO
C
      CLOSE(30)
C
      WRITE(6,*) 'Wrote q3_q14_fortran_reference.csv'
C
      END
