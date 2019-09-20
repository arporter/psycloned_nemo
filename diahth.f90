MODULE diahth
  LOGICAL, PUBLIC, PARAMETER :: lk_diahth = .FALSE.
  CONTAINS
  SUBROUTINE dia_hth(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('dia_hth', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'dia_hth: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE dia_hth
END MODULE diahth