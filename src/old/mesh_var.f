C
C     ******************************************************************
C
      MODULE MESH_VAR
C
C     ******************************************************************
C     *                                                                *
C     *   MESH VARIABLES                                               *
C     *                                                                *
C     ******************************************************************
C
      INTEGER   :: N_BASE, N_COORD, N_ZONE, N_FLOW, N_FIELD
      INTEGER   :: ZONETYPE, DATATYPE, LOCATION
      INTEGER   :: ICELLDIM, IPHYSDIM, DTYPE
      INTEGER   :: ISIZE(3,3), IRMIN(3), IRMAX(3)
      INTEGER   :: O_ZONE, OSIZE(3,3)
      INTEGER   :: CDIR, FDIR, SDIR
      INTEGER   :: FLIP

      REAL(8), POINTER :: X(:,:,:), Y(:,:,:), Z(:,:,:)
      REAL(8), POINTER :: XMOD(:,:,:), YMOD(:,:,:), ZMOD(:,:,:)

      END MODULE MESH_VAR
