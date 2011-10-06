C
C     ******************************************************************
C
      SUBROUTINE WRITE_CGNS (FNAME_OUT, PREC)
C
C     ******************************************************************
C     *                                                                *
C     *   WRITE CGNS FILE                                              *
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
      INTEGER   :: I, J, K, C(3), D(3), PREC
      REAL(8), POINTER :: XC(:,:,:), YC(:,:,:), ZC(:,:,:)

      CHARACTER :: FNAME_OUT*32

      CALL CG_OPEN_F(FNAME_OUT, CG_MODE_MODIFY, I_FILE, IER)
C
C     WRITE OUT THE MODIFIED COORDINATES
C
      IF (FLIP.EQ.1) THEN
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
            XC(D(1),D(2),D(3)) = XMOD(C(1),C(2),C(3))
            YC(D(1),D(2),D(3)) = YMOD(C(1),C(2),C(3))
            ZC(D(1),D(2),D(3)) = ZMOD(C(1),C(2),C(3))
         END DO
         END DO
         END DO
         DO I = 1,OSIZE(CDIR,1)
         DO J = 1,OSIZE(SDIR,1)
         DO K = 1,OSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = K
            XMOD(C(1),C(2),C(3)) = XC(C(1),C(2),C(3))
            YMOD(C(1),C(2),C(3)) = YC(C(1),C(2),C(3))
            ZMOD(C(1),C(2),C(3)) = ZC(C(1),C(2),C(3))
         END DO
         END DO
         END DO
      END IF
C
C     WRITE THE O-MESH COORDINATES
C
      IF (PREC == 0) THEN

         IF (ALLOCATED(XS)) DEALLOCATE(XS,YS,ZS)
         ALLOCATE(XS(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
         ALLOCATE(YS(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
         ALLOCATE(ZS(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
         DO I = 1,OSIZE(1,1)
         DO J = 1,OSIZE(2,1)
         DO K = 1,OSIZE(3,1)
            XS(I,J,K) = REAL(XMOD(I,J,K))
            YS(I,J,K) = REAL(YMOD(I,J,K))
            ZS(I,J,K) = REAL(ZMOD(I,J,K))
         END DO
         END DO
         END DO

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, O_ZONE, DATATYPE,
     .        'CoordinateX', XS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, O_ZONE, DATATYPE,
     .        'CoordinateY', YS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, O_ZONE, DATATYPE,
     .        'CoordinateZ', ZS, I_COORD, IER)


      ELSE

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, O_ZONE, DATATYPE,
     .        'CoordinateX', XMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, O_ZONE, DATATYPE,
     .        'CoordinateY', YMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, O_ZONE, DATATYPE,
     .        'CoordinateZ', ZMOD, I_COORD, IER)

      END IF
C
C     WRITE THE TIP CLEARANCE O-MESH COORDINATES
C
      IF (PREC == 0) THEN

         DEALLOCATE(XS,YS,ZS)
         ALLOCATE(XS(TSIZE(1,1), TSIZE(2,1), TSIZE(3,1)))
         ALLOCATE(YS(TSIZE(1,1), TSIZE(2,1), TSIZE(3,1)))
         ALLOCATE(ZS(TSIZE(1,1), TSIZE(2,1), TSIZE(3,1)))
         DO I = 1,TSIZE(1,1)
         DO J = 1,TSIZE(2,1)
         DO K = 1,TSIZE(3,1)
            XS(I,J,K) = REAL(XTMOD(I,J,K))
            YS(I,J,K) = REAL(YTMOD(I,J,K))
            ZS(I,J,K) = REAL(ZTMOD(I,J,K))
         END DO
         END DO
         END DO

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, T_ZONE, DATATYPE,
     .        'CoordinateX', XS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, T_ZONE, DATATYPE,
     .        'CoordinateY', YS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, T_ZONE, DATATYPE,
     .        'CoordinateZ', ZS, I_COORD, IER)

      ELSE

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, T_ZONE, DATATYPE,
     .        'CoordinateX', XTMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, T_ZONE, DATATYPE,
     .        'CoordinateY', YTMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, T_ZONE, DATATYPE,
     .        'CoordinateZ', ZTMOD, I_COORD, IER)

      END IF
C
C     WRITE THE TIP CLEARANCE H-MESH COORDINATES
C
      IF (PREC == 0) THEN

         DEALLOCATE(XS,YS,ZS)
         ALLOCATE(XS(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
         ALLOCATE(YS(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
         ALLOCATE(ZS(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
         DO I = 1,HSIZE(1,1)
         DO J = 1,HSIZE(2,1)
         DO K = 1,HSIZE(3,1)
            XS(I,J,K) = REAL(XHMOD(I,J,K))
            YS(I,J,K) = REAL(YHMOD(I,J,K))
            ZS(I,J,K) = REAL(ZHMOD(I,J,K))
         END DO
         END DO
         END DO

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateX', XS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateY', YS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateZ', ZS, I_COORD, IER)

      ELSE

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateX', XHMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateY', YHMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateZ', ZHMOD, I_COORD, IER)

      END IF

      RETURN

      END
