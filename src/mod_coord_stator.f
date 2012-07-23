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
C
C     PROJECT THE ENDS OF THE MESH ONTO THE HUB/SHROUD
C
      ALLOCATE(RO(HSIZE(CDIR,1)-1))
      ALLOCATE(RI(HSIZE(CDIR,1)-1))
      ALLOCATE(ZI(HSIZE(CDIR,1)-1))
      ALLOCATE(IY(HSIZE(CDIR,1)-1))
C
C     HUB SURFACE
C
      DO I = 1,HSIZE(CDIR,1)-1
         IY(I) = I
         C(CDIR) = I
         C(SDIR) = 1
         C(FDIR) = HSIZE(FDIR,1)
         ZI(I) = Z(C(1),C(2),C(3))
         RO(I) = SQRT(X(C(1),C(2),C(3))**2 + Y(C(1),C(2),C(3))**2)
      END DO

      CALL SSORT(ZI,IY,HSIZE(CDIR,1)-1)

      DO I=1,HSIZE(CDIR,1)-1
         RI(I) = RO(IY(I))
      END DO

      DO I = LE_IND,TE_IND
      DO J = 1,HSIZE(FDIR,1)-1
         C(CDIR) = I
         C(SDIR) = 1
         C(FDIR) = J
         R = SQRT(XMOD(C(1),C(2),C(3))**2 + YMOD(C(1),C(2),C(3))**2)
         DO K = 1,HSIZE(CDIR,1)-1
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
      DO I = 1,HSIZE(CDIR,1)-1
         IY(I) = I
         C(CDIR) = I
         C(SDIR) = HSIZE(SDIR,1)
         C(FDIR) = HSIZE(FDIR,1)
         ZI(I) = Z(C(1),C(2),C(3))
         RO(I) = SQRT(X(C(1),C(2),C(3))**2 + Y(C(1),C(2),C(3))**2)
      END DO

      CALL SSORT(ZI,IY,HSIZE(CDIR,1)-1)

      DO I=1,HSIZE(CDIR,1)-1
         RI(I) = RO(IY(I))
      END DO

      DO I = LE_IND,TE_IND
      DO J = 1,HSIZE(FDIR,1)-1
         C(CDIR) = I
         C(SDIR) = HSIZE(SDIR,1)
         C(FDIR) = J
         R = SQRT(XMOD(C(1),C(2),C(3))**2 + YMOD(C(1),C(2),C(3))**2)
         DO K = 1,HSIZE(CDIR,1)-1
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
      N = 15
      ALLOCATE(S2(HSIZE(CDIR,1),N))
      ALLOCATE(DS2(HSIZE(CDIR,1),N))
C
C     HUB END
C
      S2(1:HSIZE(CDIR,1),1:N) = 0.0
      DS2(1:HSIZE(CDIR,1),1:N) = 0.0
      DO I=1,HSIZE(CDIR,1)
      DO J=2,N
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = HSIZE(FDIR,1)
         D(CDIR) = I
         D(SDIR) = J-1
         D(FDIR) = HSIZE(FDIR,1)
         DS2(I,J) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
         S2(I,J) = S2(I,J-1) + DS2(I,J)
      END DO
      END DO

      DO I=LE_IND,TE_IND
      C(CDIR) = I
      C(SDIR) = 1
      C(FDIR) = 1
      D(CDIR) = I
      D(SDIR) = N
      D(FDIR) = 1
      R = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .         (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .         (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
      DO J=2,N
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
      DEALLOCATE(S2)
      DEALLOCATE(DS2)
      ALLOCATE(S2(HSIZE(CDIR,1),N))
      ALLOCATE(DS2(HSIZE(CDIR,1),N))
      S2(1:HSIZE(CDIR,1),1:N) = 0.0
      DS2(1:HSIZE(CDIR,1),1:N) = 0.0

      S(1) = 0.0
      DS(1) = 0.0
      DO I=1,HSIZE(CDIR,1)
      DO J=HSIZE(SDIR,1)-1,HSIZE(SDIR,1)-N+1,-1
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = HSIZE(FDIR,1)
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = HSIZE(FDIR,1)
         K = HSIZE(SDIR,1)-J+1
         DS2(I,K) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
         S2(I,K) = S2(I,K-1) + DS2(I,K)
      END DO
      END DO

      DO I=LE_IND,TE_IND
      C(CDIR) = I
      C(SDIR) = HSIZE(SDIR,1)
      C(FDIR) = 1
      D(CDIR) = I
      D(SDIR) = HSIZE(SDIR,1)-N+1
      D(FDIR) = 1
      R = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .         (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .         (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
      DO J=HSIZE(SDIR,1)-1,HSIZE(SDIR,1)-N+1,-1
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
         C(SDIR) = HSIZE(SDIR,1)-N
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = HSIZE(SDIR,1)
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
         K = HSIZE(SDIR,1)-J+1
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*R*DX
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*R*DY
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,N)*R*DZ
      END DO
      END DO
C
C     BILINEAR INTERPOLATION OF H-MESH INTERIOR
C
      ALLOCATE(SJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
      ALLOCATE(TJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
      DO I = LE_IND,TE_IND
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
     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
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
     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        SHIFT INTERIOR POINTS USING BILINEAR INTERPOLATION
         DO J = 2,HSIZE(FDIR,1)-1
C           SHIFT VECTOR AT HUB            
            C(FDIR) = J
            C(CDIR) = I
            C(SDIR) = 1
            DXH = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
            DYH = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
            DZH = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C           SHIFT VECTOR AT SHROUD            
            C(SDIR) = HSIZE(SDIR,1)
            DXS = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
            DYS = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
            DZS = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
            DO K = 2,HSIZE(SDIR,1)-1
C              SHIFT VECTOR AT BLADE SURFACE
               C(FDIR) = 1
               C(SDIR) = K
               DXB = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
               DYB = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
               DZB = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C              MOVE THE INTERIOR POINTS
               DENOM = 1.0/SJK(J,K)+1.0/(SJK(HSIZE(FDIR,1),K)-SJK(J,K))
     .               + 1.0/TJK(J,K)+1.0/(TJK(J,HSIZE(SDIR,1))-TJK(J,K))
               DX = (DXH/TJK(J,K) + DXS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DXB/SJK(J,K))/DENOM
               DY = (DYH/TJK(J,K) + DYS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DYB/SJK(J,K))/DENOM
               DZ = (DZH/TJK(J,K) + DZS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
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

CC
CC     TODO (REMOVE THIS) REMOVE FRONT/BACK MESH
CC
C      DO J=1,HSIZE(SDIR,1)
C      DO K=1,HSIZE(FDIR,1)
C      DO I=1,LE_IND-1
C         C(FDIR) = K
C         C(CDIR) = I
C         C(SDIR) = J
C         D(FDIR) = K
C         D(CDIR) = LE_IND
C         D(SDIR) = J
C         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3))
C         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3))
C         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3))
C      END DO
C      DO I=TE_IND+1,HSIZE(CDIR,1)
C         C(FDIR) = K
C         C(CDIR) = I
C         C(SDIR) = J
C         D(FDIR) = K
C         D(CDIR) = TE_IND
C         D(SDIR) = J
C         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3))
C         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3))
C         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3))
C      END DO
C      END DO
C      END DO
      
      RETURN
C
      END
