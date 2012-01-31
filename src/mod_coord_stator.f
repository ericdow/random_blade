C
C     ******************************************************************
C
      SUBROUTINE MOD_COORD_STATOR
C
C     ******************************************************************
C     *                                                                *
C     *   MODIFY COORDINATES                                           *
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
      INTEGER   :: I, J, K, C(3), D(3), JJ, N
      REAL      :: DX, DY, DZ, R, RR, MAG, DENOM, HTIP
      REAL      :: DXI, DYI, DZI, DXO, DYO, DZO
      REAL      :: DXS, DYS, DZS, DXH, DYH, DZH, DXB, DYB, DZB
      REAL(8), DIMENSION(:), ALLOCATABLE     :: DS, S, RO, RI, ZI
      REAL(8), DIMENSION(:,:), ALLOCATABLE   :: DS2, S2, SJK, TJK
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: XSH, YSH, ZSH
      INTEGER, DIMENSION(:), ALLOCATABLE     :: IY
      CHARACTER :: STR*32

      ALLOCATE(XHMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(YHMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(ZHMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))

      DO I = 1,HSIZE(1,1)
      DO J = 1,HSIZE(2,1)
      DO K = 1,HSIZE(3,1)
         XHMOD(I,J,K) = XH(I,J,K)
         YHMOD(I,J,K) = YH(I,J,K)
         ZHMOD(I,J,K) = ZH(I,J,K)
      END DO
      END DO
      END DO
C
C     READ IN BLADE SURFACE
C      
      OPEN(UNIT=304,FILE='blade_surf_mod.dat')

      READ(304,'(A,I10)') STR
      READ(304,'(A,I10)') STR

      DO I = 1,OSIZE(CDIR,1)
      DO J = 1,OSIZE(SDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = J
         READ(304,'(3E20.8)') XMOD(C(1),C(2),C(3)),YMOD(C(1),C(2),C(3)),
     .                        ZMOD(C(1),C(2),C(3))
      END DO
      END DO

      CLOSE(304)
C
C     CONVERT MESH TO CARTESIAN
C
      ALLOCATE(SJK(TSIZE(FDIR,1),TSIZE(SDIR,1)))
      ALLOCATE(TJK(TSIZE(FDIR,1),TSIZE(SDIR,1)))

      ALLOCATE(XSH(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
      ALLOCATE(YSH(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
      ALLOCATE(ZSH(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))

      DO I = 1,OSIZE(1,1)
      DO J = 1,OSIZE(2,1)
      DO K = 1,OSIZE(3,1)
         XSH(I,J,K) = XMOD(I,J,K)
         YSH(I,J,K) = YMOD(I,J,K)
         ZSH(I,J,K) = ZMOD(I,J,K)
      END DO
      END DO
      END DO
C      
C     INTERPOLATE O-MESH COORDINATES ON HUB AND SHROUD SURFACES
C
      ALLOCATE(S(OSIZE(FDIR,1)))
      ALLOCATE(DS(OSIZE(FDIR,1)))

      DO I = 1,OSIZE(CDIR,1)
      DO JJ = 1,2 
         IF (JJ.EQ.1) J = 1 
         IF (JJ.EQ.2) J = OSIZE(SDIR,1)
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1

         DX = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
         DY = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
         DZ = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
         
         S(1) = 0.0
         DS(1) = 0.0
         DO K = 2,OSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = K
            D(CDIR) = I
            D(SDIR) = J
            D(FDIR) = K-1
            DS(K) = SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                   (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                   (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
            S(K) = S(K-1) + DS(K)
         END DO
         DO K = 2,OSIZE(FDIR,1)-1
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = K
            XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) + 
     .                    DX*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1))
            YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) + 
     .                    DY*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1))
            ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) + 
     .                    DZ*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1))
         END DO
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = OSIZE(FDIR,1)
         XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3))
         YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3))
         ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3))
      END DO
      END DO
C
C     PROJECT THE ENDS OF THE MESH ONTO THE HUB/SHROUD
C
      ALLOCATE(RO(OSIZE(CDIR,1)-1))
      ALLOCATE(RI(OSIZE(CDIR,1)-1))
      ALLOCATE(ZI(OSIZE(CDIR,1)-1))
      ALLOCATE(IY(OSIZE(CDIR,1)-1))
C
C     HUB SURFACE
C
      DO I = 1,OSIZE(CDIR,1)-1
         IY(I) = I
         C(CDIR) = I
         C(SDIR) = 1
         C(FDIR) = OSIZE(FDIR,1)
         ZI(I) = Z(C(1),C(2),C(3))
         RO(I) = SQRT(X(C(1),C(2),C(3))**2 + Y(C(1),C(2),C(3))**2)
      END DO

      CALL SSORT(ZI,IY,OSIZE(CDIR,1)-1)

      DO I=1,OSIZE(CDIR,1)-1
         RI(I) = RO(IY(I))
      END DO

      DO I = 1,OSIZE(CDIR,1)
      DO J = 1,OSIZE(FDIR,1)-1
         C(CDIR) = I
         C(SDIR) = 1
         C(FDIR) = J
         R = SQRT(XMOD(C(1),C(2),C(3))**2 + YMOD(C(1),C(2),C(3))**2)
         DO K = 1,OSIZE(CDIR,1)-1
            IF (ZMOD(C(1),C(2),C(3)).GE.ZI(K)) THEN
               RR = RI(K) + (ZMOD(C(1),C(2),C(3)) - ZI(K))*
     .                      (RI(K-1) - RI(K)) / (ZI(K-1) - ZI(K))
               EXIT 
            END IF
         END DO
         XMOD(C(1),C(2),C(3)) = XMOD(C(1),C(2),C(3))*RR/R
         YMOD(C(1),C(2),C(3)) = YMOD(C(1),C(2),C(3))*RR/R
      END DO   
      END DO
C
C     SHROUD SURFACE
C
      DO I = 1,OSIZE(CDIR,1)-1
         IY(I) = I
         C(CDIR) = I
         C(SDIR) = OSIZE(SDIR,1)
         C(FDIR) = OSIZE(FDIR,1)
         ZI(I) = Z(C(1),C(2),C(3))
         RO(I) = SQRT(X(C(1),C(2),C(3))**2 + Y(C(1),C(2),C(3))**2)
      END DO

      CALL SSORT(ZI,IY,OSIZE(CDIR,1)-1)

      DO I=1,OSIZE(CDIR,1)-1
         RI(I) = RO(IY(I))
      END DO

      DO I = 1,OSIZE(CDIR,1)
      DO J = 1,OSIZE(FDIR,1)-1
         C(CDIR) = I
         C(SDIR) = OSIZE(SDIR,1)
         C(FDIR) = J
         R = SQRT(XMOD(C(1),C(2),C(3))**2 + YMOD(C(1),C(2),C(3))**2)
         DO K = 1,OSIZE(CDIR,1)-1
            IF (ZMOD(C(1),C(2),C(3)).GE.ZI(K)) THEN
               RR = RI(K) + (ZMOD(C(1),C(2),C(3)) - ZI(K))*
     .                      (RI(K-1) - RI(K)) / (ZI(K-1) - ZI(K))
               EXIT 
            END IF
         END DO
         XMOD(C(1),C(2),C(3)) = XMOD(C(1),C(2),C(3))*RR/R
         YMOD(C(1),C(2),C(3)) = YMOD(C(1),C(2),C(3))*RR/R
      END DO   
      END DO
C
C     ADJUST THE FIRST N LAYERS OF THE BLADE SURFACE MESH
C
      ALLOCATE(S2(OSIZE(CDIR,1),15))
      ALLOCATE(DS2(OSIZE(CDIR,1),15))
C
C     HUB END
C
      S2(1:OSIZE(CDIR,1),1:15) = 0.0
      DS2(1:OSIZE(CDIR,1),1:15) = 0.0
      DO I=1,OSIZE(CDIR,1)
      DO J=2,15
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = OSIZE(FDIR,1)
         D(CDIR) = I
         D(SDIR) = J-1
         D(FDIR) = OSIZE(FDIR,1)
         DS2(I,J) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
         S2(I,J) = S2(I,J-1) + DS2(I,J)
      END DO
      END DO

      DO I=1,OSIZE(CDIR,1)
      C(CDIR) = I
      C(SDIR) = 1
      C(FDIR) = 1
      D(CDIR) = I
      D(SDIR) = 15
      D(FDIR) = 1
      R = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .         (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .         (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
      DO J=2,15
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J-1
         D(FDIR) = 1
C        DIRECTION OF PROJECTED MESH         
         DX = XSH(C(1),C(2),C(3)) - XSH(D(1),D(2),D(3))
         DY = YSH(C(1),C(2),C(3)) - YSH(D(1),D(2),D(3))
         DZ = ZSH(C(1),C(2),C(3)) - ZSH(D(1),D(2),D(3))
         MAG = SQRT(DX*DX + DY*DY + DZ*DZ)
         DX = DX/MAG
         DY = DY/MAG
         DZ = DZ/MAG
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,J)/S2(I,15)*R*DX
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,J)/S2(I,15)*R*DY
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,J)/S2(I,15)*R*DZ
      END DO
      END DO
C
C     SHROUD END
C
      N = TSIZE(SDIR,1)/2
      DEALLOCATE(S2)
      DEALLOCATE(DS2)
      ALLOCATE(S2(OSIZE(CDIR,1),N))
      ALLOCATE(DS2(OSIZE(CDIR,1),N))
      S2(1:OSIZE(CDIR,1),1:N) = 0.0
      DS2(1:OSIZE(CDIR,1),1:N) = 0.0

      S(1) = 0.0
      DS(1) = 0.0
      DO I=1,OSIZE(CDIR,1)
      DO J=OSIZE(SDIR,1)-1,OSIZE(SDIR,1)-N+1,-1
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = OSIZE(FDIR,1)
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = OSIZE(FDIR,1)
         K = OSIZE(SDIR,1)-J+1
         DS2(I,K) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
         S2(I,K) = S2(I,K-1) + DS2(I,K)
      END DO
      END DO

      DO I=1,OSIZE(CDIR,1)
      C(CDIR) = I
      C(SDIR) = OSIZE(SDIR,1)
      C(FDIR) = 1
      D(CDIR) = I
      D(SDIR) = OSIZE(SDIR,1)-N+1
      D(FDIR) = 1
      R = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .         (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .         (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
      DO J=OSIZE(SDIR,1)-1,OSIZE(SDIR,1)-N+1,-1
C         C(CDIR) = I
C         C(SDIR) = J
C         C(FDIR) = 1
C         D(CDIR) = I
C         D(SDIR) = J+1
C         D(FDIR) = 1
C         DX = XSH(C(1),C(2),C(3)) - XSH(D(1),D(2),D(3))
C         DY = YSH(C(1),C(2),C(3)) - YSH(D(1),D(2),D(3))
C         DZ = ZSH(C(1),C(2),C(3)) - ZSH(D(1),D(2),D(3))
         C(CDIR) = I
         C(SDIR) = OSIZE(SDIR,1)-N
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = OSIZE(SDIR,1)
         D(FDIR) = 1
         DX = XMOD(C(1),C(2),C(3)) - XMOD(D(1),D(2),D(3))
         DY = YMOD(C(1),C(2),C(3)) - YMOD(D(1),D(2),D(3))
         DZ = ZMOD(C(1),C(2),C(3)) - ZMOD(D(1),D(2),D(3))
         MAG = SQRT(DX*DX + DY*DY + DZ*DZ)
         DX = DX/MAG
         DY = DY/MAG
         DZ = DZ/MAG
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = 1
         K = OSIZE(SDIR,1)-J+1
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*R*DX
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*R*DY
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*R*DZ
      END DO
      END DO
C
C     BILINEAR INTERPOLATION OF O-MESH INTERIOR
C
      DEALLOCATE(SJK)
      DEALLOCATE(TJK)
      ALLOCATE(SJK(OSIZE(FDIR,1),OSIZE(SDIR,1)))
      ALLOCATE(TJK(OSIZE(FDIR,1),OSIZE(SDIR,1)))
      DO I = 1,OSIZE(CDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = OSIZE(SDIR,1)

         SJK(1,1:OSIZE(SDIR,1)) = 0.0
         TJK(1:OSIZE(FDIR,1),1) = 0.0
C        FILL SJK
         DO J = 2,OSIZE(FDIR,1)
         DO K = 1,OSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J-1
            D(SDIR) = K
            SJK(J,K) = SJK(J-1,K) + 
     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        FILL TJK
         DO J = 1,OSIZE(FDIR,1)
         DO K = 2,OSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J
            D(SDIR) = K-1
            TJK(J,K) = TJK(J,K-1) + 
     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        SHIFT INTERIOR POINTS USING BILINEAR INTERPOLATION
         DO J = 2,OSIZE(FDIR,1)-1
C           SHIFT VECTOR AT HUB            
            C(FDIR) = J
            C(CDIR) = I
            C(SDIR) = 1
            DXH = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
            DYH = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
            DZH = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C           SHIFT VECTOR AT SHROUD            
            C(SDIR) = OSIZE(SDIR,1)
            DXS = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
            DYS = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
            DZS = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
            DO K = 2,OSIZE(SDIR,1)-1
C              SHIFT VECTOR AT BLADE SURFACE
               C(FDIR) = 1
               C(SDIR) = K
               DXB = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
               DYB = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
               DZB = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C              MOVE THE INTERIOR POINTS
               DENOM = 1.0/SJK(J,K)+1.0/(SJK(OSIZE(FDIR,1),K)-SJK(J,K))
     .               + 1.0/TJK(J,K)+1.0/(TJK(J,OSIZE(SDIR,1))-TJK(J,K))
               DX = (DXH/TJK(J,K) + DXS/(TJK(J,OSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DXB/SJK(J,K))/DENOM
               DY = (DYH/TJK(J,K) + DYS/(TJK(J,OSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DYB/SJK(J,K))/DENOM
               DZ = (DZH/TJK(J,K) + DZS/(TJK(J,OSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DZB/SJK(J,K))/DENOM
               C(FDIR) = J
               C(CDIR) = I
               C(SDIR) = K
               XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) + DX
               YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) + DY
               ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) + DZ
            END DO
         END DO
      END DO
C
C     MOVE THE TIP CLEARANCE O-MESH
C
      DO I=1,TSIZE(CDIR,1)
      DO J=1,TSIZE(SDIR,1)
         C(FDIR) = TSIZE(FDIR,1)
         C(SDIR) = J
         C(CDIR) = I
         D(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = OSIZE(SDIR,1)-TSIZE(SDIR,1)+J
C        MOVE OUTER LAYER
         XTMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3))
         YTMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3))
         ZTMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3))
C        MOVE INNER LAYERS
         DO K=1,TSIZE(FDIR,1)-1
            D(FDIR) = K
            D(CDIR) = I
            D(SDIR) = J
            XTMOD(D(1),D(2),D(3)) = XT(D(1),D(2),D(3)) 
     .                      + XTMOD(C(1),C(2),C(3)) - XT(C(1),C(2),C(3))
            YTMOD(D(1),D(2),D(3)) = YT(D(1),D(2),D(3)) 
     .                      + YTMOD(C(1),C(2),C(3)) - YT(C(1),C(2),C(3))
            ZTMOD(D(1),D(2),D(3)) = ZT(D(1),D(2),D(3)) 
     .                      + ZTMOD(C(1),C(2),C(3)) - ZT(C(1),C(2),C(3))
         END DO
      END DO
      END DO
C
C     PROJECT TIP CLEARANCE O-MESH SHROUD SURFACE 
C
      DO I = 1,TSIZE(CDIR,1)
      DO J = 1,TSIZE(FDIR,1)-1
         C(CDIR) = I
         C(SDIR) = TSIZE(SDIR,1)
         C(FDIR) = J
         R = SQRT(XTMOD(C(1),C(2),C(3))**2 + YTMOD(C(1),C(2),C(3))**2)
         DO K = 1,TSIZE(CDIR,1)-1
            IF (ZTMOD(C(1),C(2),C(3)).GE.ZI(K)) THEN
               RR = RI(K) + (ZTMOD(C(1),C(2),C(3)) - ZI(K))*
     .                      (RI(K-1) - RI(K)) / (ZI(K-1) - ZI(K))
               EXIT 
            END IF
         END DO
         XTMOD(C(1),C(2),C(3)) = XTMOD(C(1),C(2),C(3))*RR/R
         YTMOD(C(1),C(2),C(3)) = YTMOD(C(1),C(2),C(3))*RR/R
      END DO   
      END DO
C
C     SHIFT INTERIOR TIP CLEARANCE O-MESH EDGE POINTS
C
      DEALLOCATE(S2)
      DEALLOCATE(DS2)
      ALLOCATE(S2(TSIZE(CDIR,1),N))
      ALLOCATE(DS2(TSIZE(CDIR,1),N))

      S(1) = 0.0
      DS(1) = 0.0
      S2(1:TSIZE(CDIR,1),1:N) = 0.0
      DS2(1:TSIZE(CDIR,1),1:N) = 0.0
      DO I=1,TSIZE(CDIR,1)
      DO J=TSIZE(SDIR,1)-1,TSIZE(SDIR,1)-N+1,-1
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = TSIZE(FDIR,1)
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = TSIZE(FDIR,1)
         K = TSIZE(SDIR,1)-J+1
         DS2(I,K)=SQRT((XTMOD(C(1),C(2),C(3))-XTMOD(D(1),D(2),D(3)))**2+
     .                 (YTMOD(C(1),C(2),C(3))-YTMOD(D(1),D(2),D(3)))**2+
     .                 (ZTMOD(C(1),C(2),C(3))-ZTMOD(D(1),D(2),D(3)))**2)
         S2(I,K) = S2(I,K-1) + DS2(I,K)
      END DO
      END DO

      DO I=1,TSIZE(CDIR,1)
      S(1) = 0.0
      DS(1) = 0.0
      DO J=TSIZE(SDIR,1)-1,TSIZE(SDIR,1)-N+1,-1
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = 1
         K = TSIZE(SDIR,1)-J+1
         DS(K) = SQRT((XTMOD(C(1),C(2),C(3))-XTMOD(D(1),D(2),D(3)))**2+
     .                (YTMOD(C(1),C(2),C(3))-YTMOD(D(1),D(2),D(3)))**2+
     .                (ZTMOD(C(1),C(2),C(3))-ZTMOD(D(1),D(2),D(3)))**2)
         S(K) = S(K-1) + DS(K)
      END DO
      DO J=TSIZE(SDIR,1)-1,TSIZE(SDIR,1)-N+1,-1
C         C(CDIR) = I
C         C(SDIR) = J
C         C(FDIR) = 1
C         D(CDIR) = I
C         D(SDIR) = J+1
C         D(FDIR) = 1
C         DX = XSH(C(1),C(2),C(3)) - XSH(D(1),D(2),D(3))
C         DY = YSH(C(1),C(2),C(3)) - YSH(D(1),D(2),D(3))
C         DZ = ZSH(C(1),C(2),C(3)) - ZSH(D(1),D(2),D(3))
         C(CDIR) = I
         C(SDIR) = TSIZE(SDIR,1)-N
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = TSIZE(SDIR,1)
         D(FDIR) = 1
         DX = XTMOD(C(1),C(2),C(3)) - XTMOD(D(1),D(2),D(3))
         DY = YTMOD(C(1),C(2),C(3)) - YTMOD(D(1),D(2),D(3))
         DZ = ZTMOD(C(1),C(2),C(3)) - ZTMOD(D(1),D(2),D(3))
         MAG = SQRT(DX*DX + DY*DY + DZ*DZ)
         DX = DX/MAG
         DY = DY/MAG
         DZ = DZ/MAG
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = 1
         K = TSIZE(SDIR,1)-J+1
         XTMOD(C(1),C(2),C(3)) = XTMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*S(N)*DX
         YTMOD(C(1),C(2),C(3)) = YTMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*S(N)*DY
         ZTMOD(C(1),C(2),C(3)) = ZTMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*S(N)*DZ
      END DO
      END DO
C
C     BILINEAR INTERPOLATION OF TIP CLEARANCE O-MESH INTERIOR
C
      DEALLOCATE(SJK)
      DEALLOCATE(TJK)
      ALLOCATE(SJK(TSIZE(FDIR,1),TSIZE(SDIR,1)))
      ALLOCATE(TJK(TSIZE(FDIR,1),TSIZE(SDIR,1)))
      DO I = 1,TSIZE(CDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = TSIZE(SDIR,1)

         SJK(1,1:TSIZE(SDIR,1)) = 0.0
         TJK(1:TSIZE(FDIR,1),1) = 0.0
C        FILL SJK
         DO J = 2,TSIZE(FDIR,1)
         DO K = 1,TSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J-1
            D(SDIR) = K
            SJK(J,K) = SJK(J-1,K) + 
     .           SQRT((XT(C(1),C(2),C(3))-XT(D(1),D(2),D(3)))**2 + 
     .                (YT(C(1),C(2),C(3))-YT(D(1),D(2),D(3)))**2 +
     .                (ZT(C(1),C(2),C(3))-ZT(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        FILL TJK
         DO J = 1,TSIZE(FDIR,1)
         DO K = 2,TSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J
            D(SDIR) = K-1
            TJK(J,K) = TJK(J,K-1) + 
     .           SQRT((XT(C(1),C(2),C(3))-XT(D(1),D(2),D(3)))**2 + 
     .                (YT(C(1),C(2),C(3))-YT(D(1),D(2),D(3)))**2 +
     .                (ZT(C(1),C(2),C(3))-ZT(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        SHIFT INTERIOR POINTS USING BILINEAR INTERPOLATION
         DO J = 2,TSIZE(FDIR,1)-1
C           SHIFT VECTOR AT HUB            
            C(FDIR) = J
            C(CDIR) = I
            C(SDIR) = 1
            DXH = XTMOD(C(1),C(2),C(3)) - XT(C(1),C(2),C(3)) 
            DYH = YTMOD(C(1),C(2),C(3)) - YT(C(1),C(2),C(3)) 
            DZH = ZTMOD(C(1),C(2),C(3)) - ZT(C(1),C(2),C(3))
C           SHIFT VECTOR AT SHROUD            
            C(SDIR) = TSIZE(SDIR,1)
            DXS = XTMOD(C(1),C(2),C(3)) - XT(C(1),C(2),C(3)) 
            DYS = YTMOD(C(1),C(2),C(3)) - YT(C(1),C(2),C(3)) 
            DZS = ZTMOD(C(1),C(2),C(3)) - ZT(C(1),C(2),C(3))
            DO K = 2,TSIZE(SDIR,1)-1
C              SHIFT VECTOR AT BLADE SURFACE
               C(FDIR) = 1
               C(SDIR) = K
               DXB = XTMOD(C(1),C(2),C(3)) - XT(C(1),C(2),C(3)) 
               DYB = YTMOD(C(1),C(2),C(3)) - YT(C(1),C(2),C(3)) 
               DZB = ZTMOD(C(1),C(2),C(3)) - ZT(C(1),C(2),C(3))
C              SHIFT VECTOR AT OUTER SURFACE
               C(FDIR) = TSIZE(FDIR,1)
               C(SDIR) = K
               DXO = XTMOD(C(1),C(2),C(3)) - XT(C(1),C(2),C(3)) 
               DYO = YTMOD(C(1),C(2),C(3)) - YT(C(1),C(2),C(3)) 
               DZO = ZTMOD(C(1),C(2),C(3)) - ZT(C(1),C(2),C(3))
C              MOVE THE INTERIOR POINTS
               DENOM = 1.0/SJK(J,K)+1.0/(SJK(TSIZE(FDIR,1),K)-SJK(J,K))
     .               + 1.0/TJK(J,K)+1.0/(TJK(J,TSIZE(SDIR,1))-TJK(J,K))
               DX = (DXH/TJK(J,K) + DXS/(TJK(J,TSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DXB/SJK(J,K) + DXO/(SJK(TSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               DY = (DYH/TJK(J,K) + DYS/(TJK(J,TSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DYB/SJK(J,K) + DYO/(SJK(TSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               DZ = (DZH/TJK(J,K) + DZS/(TJK(J,TSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DZB/SJK(J,K) + DZO/(SJK(TSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               C(FDIR) = J
               C(CDIR) = I
               C(SDIR) = K
               XTMOD(C(1),C(2),C(3)) = XT(C(1),C(2),C(3)) + DX
               YTMOD(C(1),C(2),C(3)) = YT(C(1),C(2),C(3)) + DY
               ZTMOD(C(1),C(2),C(3)) = ZT(C(1),C(2),C(3)) + DZ
            END DO
         END DO
      END DO
C
C     MOVE THE OUTER FACES OF THE TIP CLEARANCE H-MESH
C
      DO I=1,HSIZE(CDIR,1)
      DO J=1,HSIZE(SDIR,1)
         C(FDIR) = 1
         C(SDIR) = J
         C(CDIR) = I
         D(FDIR) = 1
         D(SDIR) = J
         D(CDIR) = HOL(I,1)
         XHMOD(C(1),C(2),C(3)) = XTMOD(D(1),D(2),D(3))
         YHMOD(C(1),C(2),C(3)) = YTMOD(D(1),D(2),D(3))
         ZHMOD(C(1),C(2),C(3)) = ZTMOD(D(1),D(2),D(3))
         C(FDIR) = HSIZE(FDIR,1)
         D(CDIR) = HOL(I,2)
         XHMOD(C(1),C(2),C(3)) = XTMOD(D(1),D(2),D(3))
         YHMOD(C(1),C(2),C(3)) = YTMOD(D(1),D(2),D(3))
         ZHMOD(C(1),C(2),C(3)) = ZTMOD(D(1),D(2),D(3))
      END DO
      END DO
      DO I=1,HSIZE(FDIR,1)
      DO J=1,HSIZE(SDIR,1)
         C(FDIR) = I
         C(SDIR) = J
         C(CDIR) = 1
         D(FDIR) = 1
         D(SDIR) = J
         D(CDIR) = HOW(I,1)
         XHMOD(C(1),C(2),C(3)) = XTMOD(D(1),D(2),D(3))
         YHMOD(C(1),C(2),C(3)) = YTMOD(D(1),D(2),D(3))
         ZHMOD(C(1),C(2),C(3)) = ZTMOD(D(1),D(2),D(3))
         C(CDIR) = HSIZE(CDIR,1)
         D(CDIR) = HOW(I,2)
         XHMOD(C(1),C(2),C(3)) = XTMOD(D(1),D(2),D(3))
         YHMOD(C(1),C(2),C(3)) = YTMOD(D(1),D(2),D(3))
         ZHMOD(C(1),C(2),C(3)) = ZTMOD(D(1),D(2),D(3))
      END DO
      END DO
C      
C     INTERPOLATE THE UPPER AND LOWER FACES OF THE H-MESH
C
      DEALLOCATE(S,DS)
      ALLOCATE(S(HSIZE(FDIR,1)))
      ALLOCATE(DS(HSIZE(FDIR,1)))
      DO I=2,HSIZE(CDIR,1)-1
C        BLADE FACE
         S(1) = 0.0
         DS(1) = 0.0
         DO J = 2,HSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = 1
            C(FDIR) = J
            D(CDIR) = I
            D(SDIR) = 1
            D(FDIR) = J-1
            DS(J) = SQRT((XH(C(1),C(2),C(3))-XH(D(1),D(2),D(3)))**2 + 
     .                   (YH(C(1),C(2),C(3))-YH(D(1),D(2),D(3)))**2 +
     .                   (ZH(C(1),C(2),C(3))-ZH(D(1),D(2),D(3)))**2)
            S(J) = S(J-1) + DS(J)
         END DO
         C(CDIR) = I
         C(SDIR) = 1
         C(FDIR) = 1
         DXB = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3))
         DYB = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3))
         DZB = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
         C(FDIR) = HSIZE(FDIR,1)
         DXO = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3))
         DYO = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3))
         DZO = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
         DO J = 2,HSIZE(FDIR,1)-1
            DX = DXB + (DXO-DXB)*S(J)/S(HSIZE(FDIR,1))
            DY = DYB + (DYO-DYB)*S(J)/S(HSIZE(FDIR,1))
            DZ = DZB + (DZO-DZB)*S(J)/S(HSIZE(FDIR,1))
            C(CDIR) = I
            C(SDIR) = 1
            C(FDIR) = J
            XHMOD(C(1),C(2),C(3)) = XH(C(1),C(2),C(3)) + DX
            YHMOD(C(1),C(2),C(3)) = YH(C(1),C(2),C(3)) + DY
            ZHMOD(C(1),C(2),C(3)) = ZH(C(1),C(2),C(3)) + DZ
         END DO
C        SHROUD FACE
         S(1) = 0.0
         DS(1) = 0.0
         DO J = 2,HSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = HSIZE(SDIR,1)
            C(FDIR) = J
            D(CDIR) = I
            D(SDIR) = HSIZE(SDIR,1)
            D(FDIR) = J-1
            DS(J) = SQRT((XH(C(1),C(2),C(3))-XH(D(1),D(2),D(3)))**2 + 
     .                   (YH(C(1),C(2),C(3))-YH(D(1),D(2),D(3)))**2 +
     .                   (ZH(C(1),C(2),C(3))-ZH(D(1),D(2),D(3)))**2)
            S(J) = S(J-1) + DS(J)
         END DO
         C(CDIR) = I
         C(SDIR) = HSIZE(SDIR,1)
         C(FDIR) = 1
         DXB = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3))
         DYB = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3))
         DZB = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
         C(FDIR) = HSIZE(FDIR,1)
         DXO = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3))
         DYO = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3))
         DZO = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
         DO J = 2,HSIZE(FDIR,1)-1
            DX = DXB + (DXO-DXB)*S(J)/S(HSIZE(FDIR,1))
            DY = DYB + (DYO-DYB)*S(J)/S(HSIZE(FDIR,1))
            DZ = DZB + (DZO-DZB)*S(J)/S(HSIZE(FDIR,1))
            C(CDIR) = I
            C(SDIR) = HSIZE(SDIR,1)
            C(FDIR) = J
            XHMOD(C(1),C(2),C(3)) = XH(C(1),C(2),C(3)) + DX
            YHMOD(C(1),C(2),C(3)) = YH(C(1),C(2),C(3)) + DY
            ZHMOD(C(1),C(2),C(3)) = ZH(C(1),C(2),C(3)) + DZ
         END DO
      END DO
C
C     PROJECT SHROUD FACE OF TIP H-MESH
C
      DO I = 2,HSIZE(CDIR,1)-1
      DO J = 2,HSIZE(FDIR,1)-1
         C(CDIR) = I
         C(SDIR) = HSIZE(SDIR,1)
         C(FDIR) = J
         R = SQRT(XHMOD(C(1),C(2),C(3))**2 + YHMOD(C(1),C(2),C(3))**2)
         DO K = 1,OSIZE(CDIR,1)-1
            IF (ZHMOD(C(1),C(2),C(3)).GE.ZI(K)) THEN
               RR = RI(K) + (ZHMOD(C(1),C(2),C(3)) - ZI(K))*
     .                      (RI(K-1) - RI(K)) / (ZI(K-1) - ZI(K))
               EXIT 
            END IF
         END DO
         XHMOD(C(1),C(2),C(3)) = XHMOD(C(1),C(2),C(3))*RR/R
         YHMOD(C(1),C(2),C(3)) = YHMOD(C(1),C(2),C(3))*RR/R
      END DO   
      END DO
C
C     BILINEAR INTERPOLATION OF TIP H-MESH INTERIOR
C
      DEALLOCATE(SJK)
      DEALLOCATE(TJK)
      ALLOCATE(SJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
      ALLOCATE(TJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
      DO I = 2,HSIZE(CDIR,1)-1
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = HSIZE(SDIR,1)

         SJK(1,1:HSIZE(SDIR,1)) = 0.0
         TJK(1:HSIZE(FDIR,1),1) = 0.0
C        FILL SJK
         DO J = 2,HSIZE(FDIR,1)
         DO K = 1,HSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J-1
            D(SDIR) = K
            SJK(J,K) = SJK(J-1,K) + 
     .           SQRT((XH(C(1),C(2),C(3))-XH(D(1),D(2),D(3)))**2 + 
     .                (YH(C(1),C(2),C(3))-YH(D(1),D(2),D(3)))**2 +
     .                (ZH(C(1),C(2),C(3))-ZH(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        FILL TJK
         DO J = 1,HSIZE(FDIR,1)
         DO K = 2,HSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J
            D(SDIR) = K-1
            TJK(J,K) = TJK(J,K-1) + 
     .           SQRT((XH(C(1),C(2),C(3))-XH(D(1),D(2),D(3)))**2 + 
     .                (YH(C(1),C(2),C(3))-YH(D(1),D(2),D(3)))**2 +
     .                (ZH(C(1),C(2),C(3))-ZH(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        SHIFT INTERIOR POINTS USING BILINEAR INTERPOLATION
         DO J = 2,HSIZE(FDIR,1)-1
C           SHIFT VECTOR AT HUB            
            C(FDIR) = J
            C(CDIR) = I
            C(SDIR) = 1
            DXH = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3)) 
            DYH = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3)) 
            DZH = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
C           SHIFT VECTOR AT SHROUD            
            C(SDIR) = HSIZE(SDIR,1)
            DXS = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3)) 
            DYS = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3)) 
            DZS = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
            DO K = 2,HSIZE(SDIR,1)-1
C              SHIFT VECTOR AT BLADE SURFACE
               C(FDIR) = 1
               C(SDIR) = K
               DXB = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3)) 
               DYB = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3)) 
               DZB = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
C              SHIFT VECTOR AT OUTER SURFACE
               C(FDIR) = HSIZE(FDIR,1)
               C(SDIR) = K
               DXO = XHMOD(C(1),C(2),C(3)) - XH(C(1),C(2),C(3)) 
               DYO = YHMOD(C(1),C(2),C(3)) - YH(C(1),C(2),C(3)) 
               DZO = ZHMOD(C(1),C(2),C(3)) - ZH(C(1),C(2),C(3))
C              MOVE THE INTERIOR POINTS
               DENOM = 1.0/SJK(J,K)+1.0/(SJK(HSIZE(FDIR,1),K)-SJK(J,K))
     .               + 1.0/TJK(J,K)+1.0/(TJK(J,HSIZE(SDIR,1))-TJK(J,K))
               DX = (DXH/TJK(J,K) + DXS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DXB/SJK(J,K) + DXO/(SJK(HSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               DY = (DYH/TJK(J,K) + DYS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DYB/SJK(J,K) + DYO/(SJK(HSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               DZ = (DZH/TJK(J,K) + DZS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DZB/SJK(J,K) + DZO/(SJK(HSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               C(FDIR) = J
               C(CDIR) = I
               C(SDIR) = K
               XHMOD(C(1),C(2),C(3)) = XH(C(1),C(2),C(3)) + DX
               YHMOD(C(1),C(2),C(3)) = YH(C(1),C(2),C(3)) + DY
               ZHMOD(C(1),C(2),C(3)) = ZH(C(1),C(2),C(3)) + DZ
            END DO
         END DO
      END DO

      RETURN
C
      END
