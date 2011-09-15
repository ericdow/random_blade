C
C     ******************************************************************
C
      SUBROUTINE WRITE_CGNS (FNAME_OUT)
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
      INTEGER   :: I, J, K, C(3), D(3)
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

      I_BASE = 1

      CALL CG_COORD_WRITE_F(I_FILE, I_BASE, O_ZONE, RealDouble,
     .     'CoordinateX', XMOD, I_COORD, IER)

      CALL CG_COORD_WRITE_F(I_FILE, I_BASE, O_ZONE, RealDouble,
     .     'CoordinateY', YMOD, I_COORD, IER)

      CALL CG_COORD_WRITE_F(I_FILE, I_BASE, O_ZONE, RealDouble,
     .     'CoordinateZ', ZMOD, I_COORD, IER)

      CALL CG_CLOSE_F(I_FILE, IER)

      RETURN

      END
