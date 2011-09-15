C
C     ******************************************************************
C
      SUBROUTINE MOD_COORD
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
      INTEGER   :: I, J, K, C(3), D(3)
      REAL      :: DX, DY, DZ, R, RR, MAG
      REAL      :: DXI, DYI, DZI, DXO, DYO, DZO
      REAL(8), POINTER :: DS(:), S(:), RO(:), RI(:), ZI(:)
      REAL(8), POINTER :: DS2(:,:), S2(:,:)
      REAL(8), POINTER :: XSH(:,:,:), YSH(:,:,:), ZSH(:,:,:)
      INTEGER, POINTER :: IY(:)
      CHARACTER :: STR*32

      ALLOCATE(XMOD(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
      ALLOCATE(YMOD(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))
      ALLOCATE(ZMOD(OSIZE(1,1), OSIZE(2,1), OSIZE(3,1)))

      DO I = 1,OSIZE(1,1)
      DO J = 1,OSIZE(2,1)
      DO K = 1,OSIZE(3,1)
         XMOD(I,J,K) = X(I,J,K)
         YMOD(I,J,K) = Y(I,J,K)
         XMOD(I,J,K) = Z(I,J,K)
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
C     INTERPOLATE O-MESH COORDINATES
C
      ALLOCATE(S(OSIZE(FDIR,1)))
      ALLOCATE(DS(OSIZE(FDIR,1)))

      DO I = 1,OSIZE(CDIR,1)
      DO J = 1,OSIZE(SDIR,1)
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
C     PROJECT THE ENDS OF THE MESH ONTO THE HUB/SHROUD
C
      ALLOCATE(RO(OSIZE(CDIR,1)-1))
      ALLOCATE(RI(OSIZE(CDIR,1)-1))
      ALLOCATE(ZI(OSIZE(CDIR,1)-1))
      ALLOCATE(IY(OSIZE(CDIR,1)-1))
C
C     INNER SURFACE
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
C     OUTER SURFACE
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
C     ADJUST THE ENDS OF THE SURFACE MESH
C
      DEALLOCATE(S)
      DEALLOCATE(DS)
      ALLOCATE(S(15))
      ALLOCATE(DS(15))
      ALLOCATE(S2(OSIZE(CDIR,1),15))
      ALLOCATE(DS2(OSIZE(CDIR,1),15))

      S(1) = 0.0
      DS(1) = 0.0
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
      S(1) = 0.0
      DS(1) = 0.0
      DO J=2,15
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J-1
         D(FDIR) = 1
         DS(J) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
         S(J) = S(J-1) + DS(J)
      END DO
      DO J=2,15
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J-1
         D(FDIR) = 1
         DX = XSH(C(1),C(2),C(3)) - XSH(D(1),D(2),D(3))
         DY = YSH(C(1),C(2),C(3)) - YSH(D(1),D(2),D(3))
         DZ = ZSH(C(1),C(2),C(3)) - ZSH(D(1),D(2),D(3))
         MAG = SQRT(DX*DX + DY*DY + DZ*DZ)
         DX = DX/MAG
         DY = DY/MAG
         DZ = DZ/MAG
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,J)/S2(I,15)*S(15)*DX
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,J)/S2(I,15)*S(15)*DY
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,J)/S2(I,15)*S(15)*DZ
      END DO
      END DO
      
      S(1) = 0.0
      DS(1) = 0.0
      DO I=1,OSIZE(CDIR,1)
      DO J=OSIZE(SDIR,1)-1,OSIZE(SDIR,1)-15+1,-1
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
      S(1) = 0.0
      DS(1) = 0.0
      DO J=OSIZE(SDIR,1)-1,OSIZE(SDIR,1)-15+1,-1
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = 1
         K = OSIZE(SDIR,1)-J+1
         DS(K) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
         S(K) = S(K-1) + DS(K)
      END DO
      DO J=OSIZE(SDIR,1)-1,OSIZE(SDIR,1)-15+1,-1
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = 1
         DX = XSH(C(1),C(2),C(3)) - XSH(D(1),D(2),D(3))
         DY = YSH(C(1),C(2),C(3)) - YSH(D(1),D(2),D(3))
         DZ = ZSH(C(1),C(2),C(3)) - ZSH(D(1),D(2),D(3))
         MAG = SQRT(DX*DX + DY*DY + DZ*DZ)
         DX = DX/MAG
         DY = DY/MAG
         DZ = DZ/MAG
         K = OSIZE(SDIR,1)-J+1
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,15)*S(15)*DX
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,15)*S(15)*DY
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + 
     .                          DS2(I,K)/S2(I,15)*S(15)*DZ
      END DO
      END DO
C
C     ADJUST THE INTERIOR OF THE MESH
C
      DEALLOCATE(S)
      DEALLOCATE(DS)
      ALLOCATE(S(OSIZE(FDIR,1)))
      ALLOCATE(DS(OSIZE(FDIR,1)))

C     LOWER

      DO I=1,OSIZE(CDIR,1)
      DO J=2,15
         S(1) = 0.0
         DS(1) = 0.0
         DO K = 2,OSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = J-1
            C(FDIR) = K
            D(CDIR) = I
            D(SDIR) = J-1
            D(FDIR) = K-1
            DS(K) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
            S(K) = S(K-1) + DS(K)
         END DO
         DO K = 2,OSIZE(FDIR,1)-1
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = 1
            D(CDIR) = I
            D(SDIR) = J-1
            D(FDIR) = 1
            DXI = XMOD(C(1),C(2),C(3)) - XMOD(D(1),D(2),D(3))
            DYI = YMOD(C(1),C(2),C(3)) - YMOD(D(1),D(2),D(3))
            DZI = ZMOD(C(1),C(2),C(3)) - ZMOD(D(1),D(2),D(3))
            C(FDIR) = OSIZE(FDIR,1)
            D(FDIR) = OSIZE(FDIR,1)
            DXO = XMOD(C(1),C(2),C(3)) - XMOD(D(1),D(2),D(3))
            DYO = YMOD(C(1),C(2),C(3)) - YMOD(D(1),D(2),D(3))
            DZO = ZMOD(C(1),C(2),C(3)) - ZMOD(D(1),D(2),D(3))

            C(FDIR) = K
            D(FDIR) = K
            XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) +  
     .           DXI*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1)) + 
     .           DXO*S(K)/S(OSIZE(FDIR,1))
            YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) +  
     .           DYI*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1)) + 
     .           DYO*S(K)/S(OSIZE(FDIR,1))
            ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) +  
     .           DZI*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1)) + 
     .           DZO*S(K)/S(OSIZE(FDIR,1))
         END DO
      END DO
      END DO

C     UPPER

      DO I=1,OSIZE(CDIR,1)
      DO J=OSIZE(SDIR,1)-1,OSIZE(SDIR,1)-15+1,-1
         S(1) = 0.0
         DS(1) = 0.0
         DO K = 2,OSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = J+1
            C(FDIR) = K
            D(CDIR) = I
            D(SDIR) = J+1
            D(FDIR) = K-1
            DS(K) = SQRT((XMOD(C(1),C(2),C(3))-XMOD(D(1),D(2),D(3)))**2+
     .                   (YMOD(C(1),C(2),C(3))-YMOD(D(1),D(2),D(3)))**2+
     .                   (ZMOD(C(1),C(2),C(3))-ZMOD(D(1),D(2),D(3)))**2)
            S(K) = S(K-1) + DS(K)
         END DO

         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1
         D(CDIR) = I
         D(SDIR) = J+1
         D(FDIR) = 1
         DXI = XMOD(C(1),C(2),C(3)) - XMOD(D(1),D(2),D(3))
         DYI = YMOD(C(1),C(2),C(3)) - YMOD(D(1),D(2),D(3))
         DZI = ZMOD(C(1),C(2),C(3)) - ZMOD(D(1),D(2),D(3))
         C(FDIR) = OSIZE(FDIR,1)
         D(FDIR) = OSIZE(FDIR,1)
         DXO = XMOD(C(1),C(2),C(3)) - XMOD(D(1),D(2),D(3))
         DYO = YMOD(C(1),C(2),C(3)) - YMOD(D(1),D(2),D(3))
         DZO = ZMOD(C(1),C(2),C(3)) - ZMOD(D(1),D(2),D(3))

         DO K = 2,OSIZE(FDIR,1)-1

            C(FDIR) = K
            D(FDIR) = K
            XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) +  
     .           DXI*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1)) + 
     .           DXO*S(K)/S(OSIZE(FDIR,1))
            YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) +  
     .           DYI*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1)) + 
     .           DYO*S(K)/S(OSIZE(FDIR,1))
            ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) +  
     .           DZI*(S(OSIZE(FDIR,1)) - S(K))/S(OSIZE(FDIR,1)) + 
     .           DZO*S(K)/S(OSIZE(FDIR,1))
           
         END DO
      END DO
      END DO

C      OPEN(UNIT=304,FILE='blade_surf_mod_proj.dat')
C
C      I = 93
C      DO J=1,OSIZE(SDIR,1)
C      DO K=1,OSIZE(FDIR,1)
C         C(FDIR) = K
C         C(CDIR) = I
C         C(SDIR) = J
C         WRITE(304,'(3E20.8)') XMOD(C(1),C(2),C(3)), 
C     .                         YMOD(C(1),C(2),C(3)),
C     .                         ZMOD(C(1),C(2),C(3))
C      END DO
C      END DO
C
C      CLOSE(304)

      RETURN
C
      END
