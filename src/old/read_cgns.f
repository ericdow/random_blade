C
C     ******************************************************************
C
      SUBROUTINE READ_CGNS (FNAME_IN)
C
C     ******************************************************************
C     *                                                                *
C     *   READ IN CGNS FILE                                            *
C     *                                                                *
C     ******************************************************************
C
      USE MESH_VAR
C
      IMPLICIT NONE
C
      INCLUDE 'cgnslib_f.h'
C
C     ******************************************************************
C
C     LOCAL VARIABLES
C
C     ******************************************************************
C
      INTEGER   :: IER, I_FILE, I_BASE, I_ZONE, I_COORD, I_FLOW, I_FIELD

      INTEGER   :: I, J, K, N, IMAX, JMAX, KMAX, LOC
      INTEGER   :: MAXS, C(3), D(3)
      REAL      :: DX, DY, DZ, S1, S2
      REAL(8), POINTER :: SOL(:,:,:), DIJK(:,:)
      REAL(8), POINTER :: XC(:,:,:), YC(:,:,:), ZC(:,:,:)
      CHARACTER :: ZONENAME*32, SOLNNAME*32, BASENAME*32, COORDNAME*32
      CHARACTER :: FIELDNAME*32, FNAME_IN*32, STATE*32, ANAME*32
C
C     OPEN CGNS FILE TO READ
C
      CALL CG_OPEN_F(FNAME_IN, CG_MODE_READ, I_FILE, IER)

      IF (IER.NE.CG_OK) THEN
           CALL CG_ERROR_EXIT_F
           STOP
      END IF
      
      CALL CG_NBASES_F(I_FILE, N_BASE, IER)

      I_BASE = 1

      CALL CG_BASE_READ_F(I_FILE, I_BASE, BASENAME,
     .                    ICELLDIM, IPHYSDIM, IER)

C
C     DEFINE THE MESH DIMENSION
C
      CALL CG_NZONES_F(I_FILE, I_BASE, N_ZONE, IER)
      ALLOCATE(DIJK(N_ZONE,3))

      MAXS = -1

      DO I_ZONE = 1,N_ZONE

          CALL CG_ZONE_TYPE_F(I_FILE, I_BASE, I_ZONE, ZONETYPE, IER)

          IF (ZONETYPE .NE. Structured) THEN
              PRINT *, 'ZONETYPE of ZONE ', I_ZONE, ' must be ',
     .                 Structured, ', Now ', ZONETYPE
              STOP
          END IF

          CALL CG_ZONE_READ_F(I_FILE, I_BASE, I_ZONE,
     .                        ZONENAME, ISIZE, IER)
C      
C         CHECK THE NUMBER OF FLOW SOLUTIONS AND NUMBER OF FIELDS
C
          CALL CG_NSOLS_F(I_FILE, I_BASE, I_ZONE, N_FLOW, IER)

          CALL CG_NFIELDS_F(I_FILE, I_BASE, I_ZONE, 1, N_FIELD, IER)

          CALL CG_SOL_INFO_F(I_FILE, I_BASE, I_ZONE, 1, SOLNNAME, LOC, 
     .                       IER)

C          IF ((LOC.NE.Vertex).AND.(LOC.NE.CG_Null)) THEN
C              PRINT *, 'GridLocation must be Vertex'
C              STOP
C          END IF

C
C         CHECK THE COORDINATES
C
          CALL CG_NCOORDS_F(I_FILE, I_BASE, I_ZONE, N_COORD, IER)
C
C         ALLOCATE SPACE
C
          ALLOCATE(X(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          ALLOCATE(Y(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          ALLOCATE(Z(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          IRMIN(1:3) = 1
          IRMAX(1:3) = ISIZE(1:3,1)
C
C         READ THE COORDINATES
C
          CALL CG_COORD_INFO_F(I_FILE, I_BASE, I_ZONE, 1,
     .                         DATATYPE, COORDNAME, IER)

          CALL CG_COORD_READ_F(I_FILE, I_BASE, I_ZONE, COORDNAME,
     .         RealDouble, IRMIN, IRMAX, X, IER)

          CALL CG_COORD_INFO_F(I_FILE, I_BASE, I_ZONE, 2,
     .                         DATATYPE, COORDNAME, IER)

          CALL CG_COORD_READ_F(I_FILE, I_BASE, I_ZONE, COORDNAME,
     .         RealDouble, IRMIN, IRMAX, Y, IER)

          CALL CG_COORD_INFO_F(I_FILE, I_BASE, I_ZONE, 3,
     .                         DATATYPE, COORDNAME, IER)

          CALL CG_COORD_READ_F(I_FILE, I_BASE, I_ZONE, COORDNAME,
     .         RealDouble, IRMIN, IRMAX, Z, IER)
C
C         ZONE DIMENSIONS
C
          IMAX = ISIZE(1,1)
          JMAX = ISIZE(2,1)
          KMAX = ISIZE(3,1)

          DIJK(I_ZONE,1) = (X(IMAX,1,1)-X(1,1,1))**2 + 
     .         (Y(IMAX,1,1)-Y(1,1,1))**2 + (Z(IMAX,1,1)-Z(1,1,1))**2
          DIJK(I_ZONE,2) = (X(1,JMAX,1)-X(1,1,1))**2 + 
     .         (Y(1,JMAX,1)-Y(1,1,1))**2 + (Z(1,JMAX,1)-Z(1,1,1))**2
          DIJK(I_ZONE,3) = (X(1,1,KMAX)-X(1,1,1))**2 +
     .         (Y(1,1,KMAX)-Y(1,1,1))**2 + (Z(1,1,KMAX)-Z(1,1,1))**2

          IF (MIN(DIJK(I_ZONE,1),DIJK(I_ZONE,2),DIJK(I_ZONE,3))  
     .                                               < 1.0E-04) THEN
             IF (IMAX+JMAX+KMAX > MAXS) THEN
                O_ZONE = I_ZONE
                MAXS = IMAX+JMAX+KMAX
             END IF
          END IF
          
          PRINT *, ZONENAME
          PRINT *, 'IMAX      JMAX      KMAX'
          PRINT 5, IMAX, JMAX, KMAX      
C      
C         EXTRACT SOLUTION
C         
C          ALLOCATE(SOL(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
C
C          DO N = 1,N_FIELD
C             
C             CALL CG_FIELD_INFO_F(I_FILE, I_BASE, I_ZONE, 1, N, DTYPE, 
C     .                            FIELDNAME, IER)
C             PRINT *, FIELDNAME 
C
C             IF (N.EQ.4) THEN
C                 CALL CG_FIELD_READ_F(I_FILE, I_BASE, I_ZONE, 1, 
C     .                 FIELDNAME, DTYPE, IRMIN, IRMAX, SOL, IER)
C             END IF
C
C          END DO

      END DO

      IF (MAXS.EQ.-1) THEN
         PRINT *, 'O_ZONE NOT FOUND'
         STOP
      END IF

      PRINT *, 'O_ZONE: ', O_ZONE
C
C     STORE THE CORRDINATES FOR THE O-MESH
C
      CALL CG_ZONE_READ_F(I_FILE, I_BASE, O_ZONE,
     .                    ZONENAME, OSIZE, IER)

      ALLOCATE(X(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
      ALLOCATE(Y(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
      ALLOCATE(Z(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
      IRMIN(1:3) = 1
      IRMAX(1:3) = OSIZE(1:3,1)

      CALL CG_COORD_INFO_F(I_FILE, I_BASE, O_ZONE, 1,
     .                     DATATYPE, COORDNAME, IER)

      CALL CG_COORD_READ_F(I_FILE, I_BASE, O_ZONE, COORDNAME,
     .     RealDouble, IRMIN, IRMAX, X, IER)

      CALL CG_COORD_INFO_F(I_FILE, I_BASE, O_ZONE, 2,
     .                     DATATYPE, COORDNAME, IER)

      CALL CG_COORD_READ_F(I_FILE, I_BASE, O_ZONE, COORDNAME,
     .     RealDouble, IRMIN, IRMAX, Y, IER)

      CALL CG_COORD_INFO_F(I_FILE, I_BASE, O_ZONE, 3,
     .                     DATATYPE, COORDNAME, IER)

      CALL CG_COORD_READ_F(I_FILE, I_BASE, O_ZONE, COORDNAME,
     .     RealDouble, IRMIN, IRMAX, Z, IER)

      CALL CG_CLOSE_F(I_FILE, IER)
C
C     DETERMINE THE ORIENTATION OF THE O-MESH
C
      CDIR = MINLOC(DIJK(O_ZONE,1:3),1)
      SDIR = MAXLOC(DIJK(O_ZONE,1:3),1)
      FDIR = 6 - CDIR - SDIR

      S1 = 0
      S2 = 0
      DO I=1,OSIZE(CDIR,1)-1
         C(FDIR) = 1
         C(CDIR) = I+1
         C(SDIR) = 1
         DX = X(C(1),C(2),C(3))
         DY = Y(C(1),C(2),C(3))
         DZ = Z(C(1),C(2),C(3))
         C(CDIR) = I
         DX = DX - X(C(1),C(2),C(3))
         DY = DY - Y(C(1),C(2),C(3))
         DZ = DZ - Z(C(1),C(2),C(3))
         S1 = S1 + SQRT(DX*DX + DY*DY + DZ*DZ)
         C(FDIR) = OSIZE(FDIR,1)
         C(CDIR) = I+1
         DX = X(C(1),C(2),C(3))
         DY = Y(C(1),C(2),C(3))
         DZ = Z(C(1),C(2),C(3))
         C(CDIR) = I
         DX = DX - X(C(1),C(2),C(3))
         DY = DY - Y(C(1),C(2),C(3))
         DZ = DZ - Z(C(1),C(2),C(3))
         S2 = S2 + SQRT(DX*DX + DY*DY + DZ*DZ)
      END DO

      FLIP = 0
      IF (S2 < S1) THEN
         FLIP = 1
         ALLOCATE(XC(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
         ALLOCATE(YC(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
         ALLOCATE(ZC(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
         DO I = 1,OSIZE(CDIR,1)
         DO J = 1,OSIZE(SDIR,1)
         DO K = 1,OSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = OSIZE(FDIR,1)-K+1
            D(CDIR) = I
            D(SDIR) = J
            D(FDIR) = K
            XC(D(1),D(2),D(3)) = X(C(1),C(2),C(3))
            YC(D(1),D(2),D(3)) = Y(C(1),C(2),C(3))
            ZC(D(1),D(2),D(3)) = Z(C(1),C(2),C(3))
         END DO
         END DO
         END DO
         DO I = 1,OSIZE(CDIR,1)
         DO J = 1,OSIZE(SDIR,1)
         DO K = 1,OSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = K
            X(C(1),C(2),C(3)) = XC(C(1),C(2),C(3))
            Y(C(1),C(2),C(3)) = YC(C(1),C(2),C(3))
            Z(C(1),C(2),C(3)) = ZC(C(1),C(2),C(3))
         END DO
         END DO
         END DO
      END IF

C
C     WRITE OUT BLADE SURFACE
C
      OPEN(UNIT=304,FILE='blade_surf.dat')
      WRITE(304,'(A,I10)') 'CDIM:', OSIZE(CDIR,1) 
      WRITE(304,'(A,I10)') 'SDIM:', OSIZE(SDIR,1) 

      DO I=1,OSIZE(CDIR,1)
      DO J=1,OSIZE(SDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = J
         WRITE(304,'(3E20.8)') X(C(1),C(2),C(3)), Y(C(1),C(2),C(3)),
     .                         Z(C(1),C(2),C(3))
      END DO
      END DO

      CLOSE(304)
C
C     ******************************************************************
C

5     FORMAT (1X,I3,6X,I3,6X,I3)
10    FORMAT (I9,I9,I9,I9)
15    FORMAT (E26.19,1X,E26.19)

      RETURN
C
      END
