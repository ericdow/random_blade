C
C     ******************************************************************
C
      PROGRAM RANDOM_BLADE_STATOR
C
C     ******************************************************************
C     *                                                                *
C     *   CREATE RANDOM GEOMETRY                                       *
C     *                                                                *
C     ******************************************************************
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
      INTEGER   :: N, PREC
      CHARACTER :: FNAME_IN*32, FNAME_OUT*32
C
C     OBTAIN INPUT AND OUTPUT FILES
C
      N = IARGC()
      IF (N > 2) THEN
          PRINT *, 'Usage: random_blade [infile] [outfile]'
          STOP
      END IF
      CALL GETARG(1, FNAME_IN)

      PREC = 0 
      CALL READ_CGNS_STATOR(FNAME_IN, PREC)
      
      IF (N.EQ.2) THEN
        
          CALL GETARG(2, FNAME_OUT)

          CALL MOD_COORD_STATOR

          CALL WRITE_CGNS_STATOR(FNAME_OUT, PREC)

      END IF

      END
