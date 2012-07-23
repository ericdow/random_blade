C
C     ******************************************************************
C
      PROGRAM SETBC
C
C     ******************************************************************
C     *                                                                *
C     *   CONVERT FINETURBO MESH TO UTCFD FORMAT                       *
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
      INTEGER   :: PREC, OUT_ZONE, OUTSIZE(3,3)
      INTEGER   :: RDIR, TDIR, XDIR, OUTFACE
      INTEGER   :: I, J, K, N, NP, IMAX, JMAX, KMAX, LOC
      INTEGER   :: C(3), D(3)
      INTEGER, DIMENSION(:), ALLOCATABLE     :: IY
      REAL(8)   :: GC,GA,LC,TC,PC,RHOINPAV,UINPAV,CV,CP,TS,P,FOO
      REAL(8)   :: XMAX
      REAL(8)   :: DR(3), DX(3)
      REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: SOLS(:,:,:)
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SOL(:,:,:)
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: XC, YC, ZC
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: R, T
      REAL(8), DIMENSION(:,:),   ALLOCATABLE :: PT, TT, M
      REAL(8), DIMENSION(:),     ALLOCATABLE :: RF, PTF, TTF, MF
      CHARACTER :: ZONENAME*32, SOLNNAME*32, BASENAME*32, COORDNAME*32
      CHARACTER :: FIELDNAME*32, FNAME_IN*32, STR*32
C
C     OBTAIN INPUT AND OUTPUT FILES
C
      N = IARGC()
      IF (N.NE.2) THEN
          PRINT *, 'Usage: setbc [infile] [P]'
          STOP
      END IF
      CALL GETARG(1, FNAME_IN)
      CALL GETARG(2, STR)
      READ(STR,*) P
      PRINT *, P

      PREC = 0
C
C     OPEN CGNS FILE TO MODIFY
C
      CALL CG_OPEN_F(FNAME_IN, CG_MODE_MODIFY, I_FILE, IER)

      IF (IER.NE.CG_OK) THEN
           CALL CG_ERROR_EXIT_F
           STOP
      END IF
C      
      CALL CG_NBASES_F(I_FILE, N_BASE, IER)
C
      XMAX = -10000.0     
C
C     DEFINE THE MESH DIMENSION
C
      CALL CG_NZONES_F(I_FILE, 1, N_ZONE, IER)

      DO I_ZONE = 1,N_ZONE

          CALL CG_ZONE_TYPE_F(I_FILE, 1, I_ZONE, ZONETYPE, IER)

          IF (ZONETYPE .NE. Structured) THEN
              PRINT *, 'ZONETYPE of ZONE ', I_ZONE, ' must be ',
     .                 Structured, ', Now ', ZONETYPE
              STOP
          END IF

          CALL CG_ZONE_READ_F(I_FILE, 1, I_ZONE,
     .                        ZONENAME, ISIZE, IER)
C
C         CHECK THE COORDINATES
C
          CALL CG_NCOORDS_F(I_FILE, 1, I_ZONE, N_COORD, IER)
C
C         ALLOCATE SPACE
C
          IF (ALLOCATED(R)) DEALLOCATE(R,T,X)
          ALLOCATE(R(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          ALLOCATE(T(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          ALLOCATE(X(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          IRMIN(1:3) = 1
          IRMAX(1:3) = ISIZE(1:3,1)
C
C         READ THE COORDINATES
C
          IF (PREC == 0) THEN

             IF (ALLOCATED(XS)) DEALLOCATE(XS,YS,ZS)
             ALLOCATE(XS(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
             ALLOCATE(YS(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
             ALLOCATE(ZS(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))

             CALL CG_COORD_INFO_F(I_FILE, 1, I_ZONE, 1,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, 1, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, XS, IER)

             CALL CG_COORD_INFO_F(I_FILE, 1, I_ZONE, 2,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, 1, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, YS, IER)

             CALL CG_COORD_INFO_F(I_FILE, 1, I_ZONE, 3,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, 1, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, ZS, IER)

             DO I = 1,ISIZE(1,1)
             DO J = 1,ISIZE(2,1)
             DO K = 1,ISIZE(3,1)
                R(I,J,K) = DBLE(XS(I,J,K))
                T(I,J,K) = DBLE(YS(I,J,K))
                X(I,J,K) = DBLE(ZS(I,J,K))
             END DO
             END DO
             END DO

          ELSE

             CALL CG_COORD_INFO_F(I_FILE, 1, I_ZONE, 1,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, 1, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, R, IER)

             CALL CG_COORD_INFO_F(I_FILE, 1, I_ZONE, 2,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, 1, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, T, IER)

             CALL CG_COORD_INFO_F(I_FILE, 1, I_ZONE, 3,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, 1, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, X, IER)

          END IF
C
C         ZONE DIMENSIONS
C
          IMAX = ISIZE(1,1)
          JMAX = ISIZE(2,1)
          KMAX = ISIZE(3,1)

          PRINT *, ZONENAME
          PRINT *, 'IMAX      JMAX      KMAX'
          PRINT 5, IMAX, JMAX, KMAX      
C
C         DETERMINE WHICH ZONE IS THE OUTLET
C          
             DO I = 1,ISIZE(1,1)
             DO J = 1,ISIZE(2,1)
             DO K = 1,ISIZE(3,1)
                IF (X(I,J,K) > XMAX) THEN
                   XMAX = X(I,J,K)
                   OUT_ZONE = I_ZONE
                   OUTSIZE(1:3,1:3) = ISIZE(1:3,1:3)
                END IF
             END DO
             END DO
             END DO

      END DO

      PRINT *, 'OUTLET ZONE #: ',OUT_ZONE
C
C     STORE THE CORRDINATES FOR THE INLET
C
      CALL CG_ZONE_READ_F(I_FILE, 1, OUT_ZONE,
     .                    ZONENAME, OUTSIZE, IER)

      IF (ALLOCATED(R)) DEALLOCATE(R,T,X)
      ALLOCATE(R(OUTSIZE(1,1), OUTSIZE(2,1), OUTSIZE(3,1)))
      ALLOCATE(T(OUTSIZE(1,1), OUTSIZE(2,1), OUTSIZE(3,1)))
      ALLOCATE(X(OUTSIZE(1,1), OUTSIZE(2,1), OUTSIZE(3,1)))
      IRMIN(1:3) = 1
      IRMAX(1:3) = OUTSIZE(1:3,1)

      IF (PREC == 0) THEN

         IF (ALLOCATED(XS)) DEALLOCATE(XS,YS,ZS)
         ALLOCATE(XS(OUTSIZE(1,1), OUTSIZE(2,1), OUTSIZE(3,1)))
         ALLOCATE(YS(OUTSIZE(1,1), OUTSIZE(2,1), OUTSIZE(3,1)))
         ALLOCATE(ZS(OUTSIZE(1,1), OUTSIZE(2,1), OUTSIZE(3,1)))

         CALL CG_COORD_INFO_F(I_FILE, 1, OUT_ZONE, 1,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, 1, OUT_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, XS, IER)

         CALL CG_COORD_INFO_F(I_FILE, 1, OUT_ZONE, 2,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, 1, OUT_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, YS, IER)

         CALL CG_COORD_INFO_F(I_FILE, 1, OUT_ZONE, 3,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, 1, OUT_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, ZS, IER)

         DO I = 1,OUTSIZE(1,1)
         DO J = 1,OUTSIZE(2,1)
         DO K = 1,OUTSIZE(3,1)
            R(I,J,K) = DBLE(XS(I,J,K))
            T(I,J,K) = DBLE(YS(I,J,K))
            X(I,J,K) = DBLE(ZS(I,J,K))
         END DO
         END DO
         END DO

      ELSE

         CALL CG_COORD_INFO_F(I_FILE, 1, OUT_ZONE, 1,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, 1, OUT_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, R, IER)

         CALL CG_COORD_INFO_F(I_FILE, 1, OUT_ZONE, 2,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, 1, OUT_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, T, IER)

         CALL CG_COORD_INFO_F(I_FILE, 1, OUT_ZONE, 3,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, 1, OUT_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, X, IER)

      END IF
C
C     ORIENT THE OUTLET BLOCK
C
      DO I=1,3
         C(1:3) = 1
         C(I) = OUTSIZE(I,1)
         DR(I) = ABS(R(C(1),C(2),C(3)) - R(1,1,1))
         DX(I) = ABS(X(C(1),C(2),C(3)) - X(1,1,1))
      END DO
      RDIR = MAXLOC(DR(1:3),1)
      XDIR = MAXLOC(DX(1:3),1)
      TDIR = 6 - RDIR - XDIR
C
      C(1:3) = 1
      OUTFACE = OUTSIZE(XDIR,1)
      IF (ABS((X(C(1),C(2),C(3)) - XMAX)/XMAX) < 0.01) OUTFACE = 1
C      
C     CHECK THE NUMBER OF FLOW SOLUTIONS AND NUMBER OF FIELDS
C
      CALL CG_NSOLS_F(I_FILE, 1, OUT_ZONE, N_FLOW, IER)

      CALL CG_NFIELDS_F(I_FILE, 1, OUT_ZONE, 1, N_FIELD, IER)

      CALL CG_SOL_INFO_F(I_FILE, 1, OUT_ZONE, 1, SOLNNAME, LOC, IER)
C      
C     EXTRACT AND MODIFY SOLUTION
C         
      ALLOCATE(SOL(OUTSIZE(1,1), OUTSIZE(2,1), OUTSIZE(3,1)))
      IF (PREC==0) ALLOCATE(SOLS(OUTSIZE(1,1), OUTSIZE(2,1), 
     .                           OUTSIZE(3,1)))

      DO N = 1,N_FIELD
         
         CALL CG_FIELD_INFO_F(I_FILE, 1, OUT_ZONE, 1, N, DTYPE, 
     .                        FIELDNAME, IER)
         IRMIN(1:3) = 1
         IRMAX(1:3) = OUTSIZE(1:3,1)
C
C        PRESSURE
C            
         IF (FIELDNAME.EQ.'Pressure') THEN
             IF (PREC.EQ.0) THEN
                CALL CG_FIELD_READ_F(I_FILE, 1, OUT_ZONE, 1, 
     .                  FIELDNAME, DTYPE, IRMIN, IRMAX, SOLS, IER)
                DO I = 1,OUTSIZE(1,1)
                DO J = 1,OUTSIZE(2,1)
                DO K = 1,OUTSIZE(3,1)
                   SOL(I,J,K) = DBLE(SOLS(I,J,K))
                END DO
                END DO
                END DO
             ELSE
                CALL CG_FIELD_READ_F(I_FILE, 1, OUT_ZONE, 1, 
     .                  FIELDNAME, DTYPE, IRMIN, IRMAX, SOL, IER)
             END IF
C
C            MODIFY OUTLET FACE OF SOL
C            
             DO I = 1,OUTSIZE(RDIR,1)
             DO J = 1,OUTSIZE(TDIR,1)
                C(RDIR) = I
                C(TDIR) = J
                C(XDIR) = OUTFACE
                SOL(C(1),C(2),C(3)) = P
             END DO
             END DO
C
C            WRITE BACK TO CGNS FILE
C            
             IF (PREC.EQ.0) THEN
                DO I = 1,OUTSIZE(1,1)
                DO J = 1,OUTSIZE(2,1)
                DO K = 1,OUTSIZE(3,1)
                   SOLS(I,J,K) = SNGL(SOL(I,J,K))
                END DO
                END DO
                END DO
                CALL CG_FIELD_WRITE_F(I_FILE, 1, OUT_ZONE, 1, 
     .                  DTYPE, FIELDNAME, SOLS, N,  IER)
             ELSE
                CALL CG_FIELD_WRITE_F(I_FILE, 1, OUT_ZONE, 1, 
     .                  DTYPE, FIELDNAME, SOL, N,  IER)
             END IF
         END IF

      END DO
C
C     ******************************************************************
C
      CALL CG_CLOSE_F(I_FILE, IER)
C
C     ******************************************************************
C
5     FORMAT (1X,I3,6X,I3,6X,I3)
10    FORMAT (I9,I9,I9,I9)
15    FORMAT (E26.19,1X,E26.19)
C
      END
