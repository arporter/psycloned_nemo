MODULE dyncor_c1d
  CONTAINS
  SUBROUTINE cor_c1d
    IMPLICIT NONE
  END SUBROUTINE cor_c1d
  SUBROUTINE dyn_cor_c1d(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('dyn_cor_c1d', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'dyn_cor_c1d: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE dyn_cor_c1d
END MODULE dyncor_c1d