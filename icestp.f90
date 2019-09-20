MODULE icestp
  CONTAINS
  SUBROUTINE ice_stp(kt, ksbc)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kt, ksbc
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('ice_stp', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'ice_stp: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE ice_stp
  SUBROUTINE ice_init
  END SUBROUTINE ice_init
END MODULE icestp