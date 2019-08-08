MODULE step_c1d
  CONTAINS
  SUBROUTINE stp_c1d(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('stp_c1d', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'stp_c1d: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE stp_c1d
END MODULE step_c1d