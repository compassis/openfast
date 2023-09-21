        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 21 18:54:22 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DISCON__genmod
          INTERFACE 
            SUBROUTINE DISCON(AVRSWAP,FROM_SC_GLOB,FROM_SC,TO_SC,AVIFAIL&
     &,ACCINFILE,AVCOUTNAME,AVCMSG) BIND(C, NAME = 'DISCON')
              REAL(KIND=4), INTENT(INOUT) :: AVRSWAP(*)
              REAL(KIND=4), INTENT(IN) :: FROM_SC_GLOB(*)
              REAL(KIND=4), INTENT(IN) :: FROM_SC(*)
              REAL(KIND=4), INTENT(INOUT) :: TO_SC(*)
              INTEGER(KIND=4), INTENT(INOUT) :: AVIFAIL
              CHARACTER(LEN=1), INTENT(IN) :: ACCINFILE(NINT(AVRSWAP((50&
     &))))
              CHARACTER(LEN=1), INTENT(IN) :: AVCOUTNAME(NINT(AVRSWAP(( &
     &51))))
              CHARACTER(LEN=1), INTENT(INOUT) :: AVCMSG(NINT(AVRSWAP((49&
     &))))
            END SUBROUTINE DISCON
          END INTERFACE 
        END MODULE DISCON__genmod
