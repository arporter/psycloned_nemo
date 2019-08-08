MODULE floats
  CONTAINS
  SUBROUTINE flo_stp(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('flo_stp', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'flo_stp: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE flo_stp
  SUBROUTINE flo_init
    IMPLICIT NONE
  END SUBROUTINE flo_init
END MODULE floats