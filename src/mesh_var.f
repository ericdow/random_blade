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
      INTEGER   :: N_BASE, N_COORD, N_ZONE, N_FLOW, N_FIELD, M_BASE
      INTEGER   :: ZONETYPE, DATATYPE, LOCATION
      INTEGER   :: ICELLDIM, IPHYSDIM, DTYPE
      INTEGER   :: ISIZE(3,3), IRMIN(3), IRMAX(3)
      INTEGER   :: O_ZONE,OSIZE(3,3),T_ZONE,TSIZE(3,3),H_ZONE,HSIZE(3,3)
      INTEGER   :: CDIR, FDIR, SDIR
      INTEGER   :: FLIP

      INTEGER, POINTER :: HOL(:,:), HOW(:,:)
      REAL(8), POINTER :: X(:,:,:), Y(:,:,:), Z(:,:,:)
      REAL(8), POINTER :: XT(:,:,:), YT(:,:,:), ZT(:,:,:)
      REAL(8), POINTER :: XH(:,:,:), YH(:,:,:), ZH(:,:,:)
      REAL(8), POINTER :: XMOD(:,:,:), YMOD(:,:,:), ZMOD(:,:,:)
      REAL(8), POINTER :: XTMOD(:,:,:), YTMOD(:,:,:), ZTMOD(:,:,:)
      REAL(8), POINTER :: XHMOD(:,:,:), YHMOD(:,:,:), ZHMOD(:,:,:)

      END MODULE MESH_VAR
