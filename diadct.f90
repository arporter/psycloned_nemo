MODULE diadct
  LOGICAL, PUBLIC, PARAMETER :: lk_diadct = .FALSE.
  PUBLIC
  CONTAINS
  SUBROUTINE dia_dct_init
    IMPLICIT NONE
    WRITE(*, FMT = *) 'dia_dct_init: You should not have seen this print! error?'
  END SUBROUTINE dia_dct_init
  SUBROUTINE dia_dct(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('dia_dct', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'dia_dct: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE dia_dct
END MODULE diadct