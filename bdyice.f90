MODULE bdyice
  CONTAINS
  SUBROUTINE bdy_ice(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('bdy_ice', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'bdy_ice: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE bdy_ice
END MODULE bdyice