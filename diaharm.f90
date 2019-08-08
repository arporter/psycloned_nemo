MODULE diaharm
  LOGICAL, PUBLIC, PARAMETER :: lk_diaharm = .FALSE.
  CONTAINS
  SUBROUTINE dia_harm(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT( IN ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('dia_harm', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'dia_harm: you should not have seen this print'
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE dia_harm
END MODULE diaharm