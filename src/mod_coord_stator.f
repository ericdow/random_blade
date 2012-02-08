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
      REAL      :: XTMP, YTMP, ZTMP, RTMP, TH1, TH2, DTH, PI
      REAL      :: DXS, DYS, DZS, DXH, DYH, DZH, DXB, DYB, DZB
      REAL(8), DIMENSION(:), ALLOCATABLE     :: DS, S, RO, RI, ZI
      REAL(8), DIMENSION(:,:), ALLOCATABLE   :: DS2, S2, SJK, TJK
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: XSH, YSH, ZSH
      INTEGER, DIMENSION(:), ALLOCATABLE     :: IY
      CHARACTER :: STR*32

      PI = 4.0*ATAN(1.0)

      ALLOCATE(XMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(YMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(ZMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))

      DO I = 1,HSIZE(1,1)
      DO J = 1,HSIZE(2,1)
      DO K = 1,HSIZE(3,1)
         XMOD(I,J,K) = X(I,J,K)
         YMOD(I,J,K) = Y(I,J,K)
         ZMOD(I,J,K) = Z(I,J,K)
      END DO
      END DO
      END DO
C
C     READ IN BLADE SURFACE
C      
      OPEN(UNIT=304,FILE='blade_surf_mod.dat')

      READ(304,'(A,I10)') STR
      READ(304,'(A,I10)') STR

C     PRESSURE SIDE
      DO I=LE_IND,TE_IND
      DO J=1,HSIZE(SDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = J
         READ(304,'(3E20.8)') XMOD(C(1),C(2),C(3)),
     .                        YMOD(C(1),C(2),C(3)),
     .                        ZMOD(C(1),C(2),C(3))   
      END DO
      END DO

C     SUCTION SIDE      
      DO I=TE_IND-1,LE_IND+1,-1
      DO J=1,HSIZE(SDIR,1)
         READ(304,'(3E20.8)') XTMP,YTMP,ZTMP
C        ROTATE SUCTION SIDE FOR 46 BLADES
         C(FDIR) = HSIZE(FDIR,1)
         C(CDIR) = I
         C(SDIR) = J
         D(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J
         TH1 = ATAN2(ZTMP,YTMP)
         TH2 = ATAN2(ZMOD(D(1),D(2),D(3)),YMOD(D(1),D(2),D(3)))
         RTMP = SQRT(YTMP**2 + ZTMP**2)
         IF (TH2 > TH1) THEN
            DTH = 2.0*PI/46.0
         ELSE
            DTH = -2.0*PI/46.0
         ENDIF        
C        TODO HACK HERE (ANGLE DETECT WRONG)          
         DTH = -2.0*PI/46.0
         XMOD(C(1),C(2),C(3)) = XTMP
         YMOD(C(1),C(2),C(3)) = RTMP*COS(TH1+DTH)
         ZMOD(C(1),C(2),C(3)) = RTMP*SIN(TH1+DTH) 
      END DO
      END DO

      CLOSE(304)
C
C     FIX SUCTION SIDE LE/TE
C     
      DO J=1,HSIZE(SDIR,1)
C        LE
         C(FDIR) = HSIZE(FDIR,1)
         C(CDIR) = LE_IND
         C(SDIR) = J
         D(FDIR) = 1
         D(CDIR) = LE_IND
         D(SDIR) = J
         TH1 = ATAN2(ZMOD(C(1),C(2),C(3)),YMOD(C(1),C(2),C(3)))
         TH2 = ATAN2(ZMOD(D(1),D(2),D(3)),YMOD(D(1),D(2),D(3)))
         RTMP = SQRT(YMOD(D(1),D(2),D(3))**2 + ZMOD(D(1),D(2),D(3))**2)
         IF (TH2 > TH1) THEN
            DTH = -2.0*PI/46.0
         ELSE
            DTH = 2.0*PI/46.0
         ENDIF
         YMOD(C(1),C(2),C(3)) = RTMP*COS(TH2+DTH)
         ZMOD(C(1),C(2),C(3)) = RTMP*SIN(TH2+DTH) 
C        TE
         C(CDIR) = TE_IND
         D(CDIR) = TE_IND
         TH1 = ATAN2(ZMOD(C(1),C(2),C(3)),YMOD(C(1),C(2),C(3)))
         TH2 = ATAN2(ZMOD(D(1),D(2),D(3)),YMOD(D(1),D(2),D(3)))
         RTMP = SQRT(YMOD(D(1),D(2),D(3))**2 + ZMOD(D(1),D(2),D(3))**2)
         IF (TH2 > TH1) THEN
            DTH = -2.0*PI/46.0
         ELSE
            DTH = 2.0*PI/46.0
         ENDIF
         YMOD(C(1),C(2),C(3)) = RTMP*COS(TH2+DTH)
         ZMOD(C(1),C(2),C(3)) = RTMP*SIN(TH2+DTH) 
      END DO
      
      ALLOCATE(XSH(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(YSH(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(ZSH(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))

      DO I = 1,HSIZE(1,1)
      DO J = 1,HSIZE(2,1)
      DO K = 1,HSIZE(3,1)
         XSH(I,J,K) = XMOD(I,J,K)
         YSH(I,J,K) = YMOD(I,J,K)
         ZSH(I,J,K) = ZMOD(I,J,K)
      END DO
      END DO
      END DO
C      
C     INTERPOLATE H-MESH COORDINATES ON HUB AND SHROUD SURFACES
C
      ALLOCATE(S(HSIZE(FDIR,1)))
      ALLOCATE(DS(HSIZE(FDIR,1)))

      DO I = LE_IND,TE_IND
      DO JJ = 1,2 
         IF (JJ.EQ.1) J = 1 
         IF (JJ.EQ.2) J = HSIZE(SDIR,1)
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1

         DX = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
         DY = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
         DZ = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
         
         S(1) = 0.0
         DS(1) = 0.0
         DO K = 2,HSIZE(FDIR,1)
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
         DO K = 2,HSIZE(FDIR,1)-1
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = K
            XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) + 
     .                    DX*(S(HSIZE(FDIR,1)) - S(K))/S(HSIZE(FDIR,1))
            YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) + 
     .                    DY*(S(HSIZE(FDIR,1)) - S(K))/S(HSIZE(FDIR,1))
            ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) + 
     .                    DZ*(S(HSIZE(FDIR,1)) - S(K))/S(HSIZE(FDIR,1))
         END DO
      END DO
      END DO
CC
CC     PROJECT THE ENDS OF THE MESH ONTO THE HUB/SHROUD
CC
C      ALLOCATE(RO(HSIZE(CDIR,1)-1))
C      ALLOCATE(RI(HSIZE(CDIR,1)-1))
C      ALLOCATE(ZI(HSIZE(CDIR,1)-1))
C      ALLOCATE(IY(HSIZE(CDIR,1)-1))
CC
CC     HUB SURFACE
CC
C      DO I = 1,HSIZE(CDIR,1)-1
C         IY(I) = I
C         C(CDIR) = I
C         C(SDIR) = 1
C         C(FDIR) = HSIZE(FDIR,1)
C         ZI(I) = Z(C(1),C(2),C(3))
C         RO(I) = SQRT(X(C(1),C(2),C(3))**2 + Y(C(1),C(2),C(3))**2)
C      END DO
C
C      CALL SSORT(ZI,IY,HSIZE(CDIR,1)-1)
C
C      DO I=1,HSIZE(CDIR,1)-1
C         RI(I) = RO(IY(I))
C      END DO
C
C      DO I = LE_IND,TE_IND
C      DO J = 1,HSIZE(FDIR,1)-1
C         C(CDIR) = I
C         C(SDIR) = 1
C         C(FDIR) = J
C         R = SQRT(XMOD(C(1),C(2),C(3))**2 + YMOD(C(1),C(2),C(3))**2)
C         DO K = 1,HSIZE(CDIR,1)-1
C            IF (ZMOD(C(1),C(2),C(3)).GE.ZI(K)) THEN
C               RR = RI(K) + (ZMOD(C(1),C(2),C(3)) - ZI(K))*
C     .                      (RI(K-1) - RI(K)) / (ZI(K-1) - ZI(K))
C               EXIT 
C            END IF
C         END DO
C         XMOD(C(1),C(2),C(3)) = XMOD(C(1),C(2),C(3))*RR/R
C         YMOD(C(1),C(2),C(3)) = YMOD(C(1),C(2),C(3))*RR/R
C      END DO   
C      END DO
CC
CC     SHROUD SURFACE
CC
C      DO I = 1,HSIZE(CDIR,1)-1
C         IY(I) = I
C         C(CDIR) = I
C         C(SDIR) = HSIZE(SDIR,1)
C         C(FDIR) = HSIZE(FDIR,1)
C         ZI(I) = Z(C(1),C(2),C(3))
C         RO(I) = SQRT(X(C(1),C(2),C(3))**2 + Y(C(1),C(2),C(3))**2)
C      END DO
C
C      CALL SSORT(ZI,IY,HSIZE(CDIR,1)-1)
C
C      DO I=1,HSIZE(CDIR,1)-1
C         RI(I) = RO(IY(I))
C      END DO
C
C      DO I = LE_IND,TE_IND
C      DO J = 1,HSIZE(FDIR,1)-1
C         C(CDIR) = I
C         C(SDIR) = HSIZE(SDIR,1)
C         C(FDIR) = J
C         R = SQRT(XMOD(C(1),C(2),C(3))**2 + YMOD(C(1),C(2),C(3))**2)
C         DO K = 1,HSIZE(CDIR,1)-1
C            IF (ZMOD(C(1),C(2),C(3)).GE.ZI(K)) THEN
C               RR = RI(K) + (ZMOD(C(1),C(2),C(3)) - ZI(K))*
C     .                      (RI(K-1) - RI(K)) / (ZI(K-1) - ZI(K))
C               EXIT 
C            END IF
C         END DO
C         XMOD(C(1),C(2),C(3)) = XMOD(C(1),C(2),C(3))*RR/R
C         YMOD(C(1),C(2),C(3)) = YMOD(C(1),C(2),C(3))*RR/R
C      END DO   
C      END DO
CC
CC     ADJUST THE FIRST N LAYERS OF THE BLADE SURFACE MESH
CC
C      N = 15
C      ALLOCATE(S2(HSIZE(CDIR,1),N))
C      ALLOCATE(DS2(HSIZE(CDIR,1),N))
CC
CC     HUB END
CC
C      S2(1:HSIZE(CDIR,1),1:N) = 0.0
C      DS2(1:HSIZE(CDIR,1),1:N) = 0.0
C      DO I=1,HSIZE(CDIR,1)
C      DO J=2,N
C         C(CDIR) = I
C         C(SDIR) = J
C         C(FDIR) = HSIZE(FDIR,1)
C         D(CDIR) = I
C         D(SDIR) = J-1
C         D(FDIR) = HSIZE(FDIR,1)
C         DS2(I,J) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
C     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
C     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
C         S2(I,J) = S2(I,J-1) + DS2(I,J)
C      END DO
C      END DO
C
C      DO I=LE_IND,TE_IND
C      C(CDIR) = I
C      C(SDIR) = 1
C      C(FDIR) = 1
C      D(CDIR) = I
C      D(SDIR) = N
C      D(FDIR) = 1
C      R = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
C     .         (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
C     .         (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
C      DO J=2,N
C         C(CDIR) = I
C         C(SDIR) = J
C         C(FDIR) = 1
C         D(CDIR) = I
C         D(SDIR) = J-1
C         D(FDIR) = 1
CC        DIRECTION OF PROJECTED MESH         
C         DX = XSH(C(1),C(2),C(3)) - XSH(D(1),D(2),D(3))
C         DY = YSH(C(1),C(2),C(3)) - YSH(D(1),D(2),D(3))
C         DZ = ZSH(C(1),C(2),C(3)) - ZSH(D(1),D(2),D(3))
C         MAG = SQRT(DX*DX + DY*DY + DZ*DZ)
C         DX = DX/MAG
C         DY = DY/MAG
C         DZ = DZ/MAG
C         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + 
C     .                          DS2(I,J)/S2(I,15)*R*DX
C         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + 
C     .                          DS2(I,J)/S2(I,15)*R*DY
C         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + 
C     .                          DS2(I,J)/S2(I,15)*R*DZ
C      END DO
C      END DO
CC
CC     SHROUD END
CC
C      DEALLOCATE(S2)
C      DEALLOCATE(DS2)
C      ALLOCATE(S2(HSIZE(CDIR,1),N))
C      ALLOCATE(DS2(HSIZE(CDIR,1),N))
C      S2(1:HSIZE(CDIR,1),1:N) = 0.0
C      DS2(1:HSIZE(CDIR,1),1:N) = 0.0
C
C      S(1) = 0.0
C      DS(1) = 0.0
C      DO I=1,HSIZE(CDIR,1)
C      DO J=HSIZE(SDIR,1)-1,HSIZE(SDIR,1)-N+1,-1
C         C(CDIR) = I
C         C(SDIR) = J
C         C(FDIR) = HSIZE(FDIR,1)
C         D(CDIR) = I
C         D(SDIR) = J+1
C         D(FDIR) = HSIZE(FDIR,1)
C         K = HSIZE(SDIR,1)-J+1
C         DS2(I,K) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
C     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
C     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
C         S2(I,K) = S2(I,K-1) + DS2(I,K)
C      END DO
C      END DO
C
C      DO I=LE_IND,TE_IND
C      C(CDIR) = I
C      C(SDIR) = HSIZE(SDIR,1)
C      C(FDIR) = 1
C      D(CDIR) = I
C      D(SDIR) = HSIZE(SDIR,1)-N+1
C      D(FDIR) = 1
C      R = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
C     .         (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
C     .         (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
C      DO J=HSIZE(SDIR,1)-1,HSIZE(SDIR,1)-N+1,-1
CC         C(CDIR) = I
CC         C(SDIR) = J
CC         C(FDIR) = 1
CC         D(CDIR) = I
CC         D(SDIR) = J+1
CC         D(FDIR) = 1
CC         DX = XSH(C(1),C(2),C(3)) - XSH(D(1),D(2),D(3))
CC         DY = YSH(C(1),C(2),C(3)) - YSH(D(1),D(2),D(3))
CC         DZ = ZSH(C(1),C(2),C(3)) - ZSH(D(1),D(2),D(3))
C         C(CDIR) = I
C         C(SDIR) = HSIZE(SDIR,1)-N
C         C(FDIR) = 1
C         D(CDIR) = I
C         D(SDIR) = HSIZE(SDIR,1)
C         D(FDIR) = 1
C         DX = XMOD(C(1),C(2),C(3)) - XMOD(D(1),D(2),D(3))
C         DY = YMOD(C(1),C(2),C(3)) - YMOD(D(1),D(2),D(3))
C         DZ = ZMOD(C(1),C(2),C(3)) - ZMOD(D(1),D(2),D(3))
C         MAG = SQRT(DX*DX + DY*DY + DZ*DZ)
C         DX = DX/MAG
C         DY = DY/MAG
C         DZ = DZ/MAG
C         C(CDIR) = I
C         C(SDIR) = J
C         C(FDIR) = 1
C         D(CDIR) = I
C         D(SDIR) = J+1
C         D(FDIR) = 1
C         K = HSIZE(SDIR,1)-J+1
C         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + 
C     .                          DS2(I,K)/S2(I,N)*R*DX
C         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + 
C     .                          DS2(I,K)/S2(I,N)*R*DY
C         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + 
C     .                          DS2(I,K)/S2(I,N)*R*DZ
C      END DO
C      END DO
CC
CC     BILINEAR INTERPOLATION OF H-MESH INTERIOR
CC
C      ALLOCATE(SJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
C      ALLOCATE(TJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
C      DO I = LE_IND,TE_IND
C         C(FDIR) = 1
C         C(CDIR) = I
C         C(SDIR) = HSIZE(SDIR,1)
C
C         SJK(1,1:HSIZE(SDIR,1)) = 0.0
C         TJK(1:HSIZE(FDIR,1),1) = 0.0
CC        FILL SJK
C         DO J = 2,HSIZE(FDIR,1)
C         DO K = 1,HSIZE(SDIR,1)
C            C(FDIR) = J
C            C(SDIR) = K
C            D(CDIR) = I
C            D(FDIR) = J-1
C            D(SDIR) = K
C            SJK(J,K) = SJK(J-1,K) + 
C     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
C     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
C     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
C         END DO
C         END DO
CC        FILL TJK
C         DO J = 1,HSIZE(FDIR,1)
C         DO K = 2,HSIZE(SDIR,1)
C            C(FDIR) = J
C            C(SDIR) = K
C            D(CDIR) = I
C            D(FDIR) = J
C            D(SDIR) = K-1
C            TJK(J,K) = TJK(J,K-1) + 
C     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
C     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
C     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
C         END DO
C         END DO
CC        SHIFT INTERIOR POINTS USING BILINEAR INTERPOLATION
C         DO J = 2,HSIZE(FDIR,1)-1
CC           SHIFT VECTOR AT HUB            
C            C(FDIR) = J
C            C(CDIR) = I
C            C(SDIR) = 1
C            DXH = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
C            DYH = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
C            DZH = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
CC           SHIFT VECTOR AT SHROUD            
C            C(SDIR) = HSIZE(SDIR,1)
C            DXS = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
C            DYS = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
C            DZS = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C            DO K = 2,HSIZE(SDIR,1)-1
CC              SHIFT VECTOR AT BLADE SURFACE
C               C(FDIR) = 1
C               C(SDIR) = K
C               DXB = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
C               DYB = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
C               DZB = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
CC              MOVE THE INTERIOR POINTS
C               DENOM = 1.0/SJK(J,K)+1.0/(SJK(HSIZE(FDIR,1),K)-SJK(J,K))
C     .               + 1.0/TJK(J,K)+1.0/(TJK(J,HSIZE(SDIR,1))-TJK(J,K))
C               DX = (DXH/TJK(J,K) + DXS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
C     .            +  DXB/SJK(J,K))/DENOM
C               DY = (DYH/TJK(J,K) + DYS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
C     .            +  DYB/SJK(J,K))/DENOM
C               DZ = (DZH/TJK(J,K) + DZS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
C     .            +  DZB/SJK(J,K))/DENOM
C               C(FDIR) = J
C               C(CDIR) = I
C               C(SDIR) = K
C               XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) + DX
C               YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) + DY
C               ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) + DZ
C            END DO
C         END DO
C      END DO

C
C     TODO (REMOVE THIS) REMOVE FRONT/BACK MESH
C
      DO J=1,HSIZE(SDIR,1)
      DO K=1,HSIZE(FDIR,1)
      DO I=1,LE_IND-1
         C(FDIR) = K
         C(CDIR) = I
         C(SDIR) = J
         D(FDIR) = K
         D(CDIR) = LE_IND
         D(SDIR) = J
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3))
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3))
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3))
      END DO
      DO I=TE_IND+1,HSIZE(CDIR,1)
         C(FDIR) = K
         C(CDIR) = I
         C(SDIR) = J
         D(FDIR) = K
         D(CDIR) = TE_IND
         D(SDIR) = J
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3))
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3))
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3))
      END DO
      END DO
      END DO
      
      RETURN
C
      END
