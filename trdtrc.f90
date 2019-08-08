MODULE trdtrc
  CONTAINS
  SUBROUTINE trd_trc(ptrtrd, kjn, ktrd, kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER :: kt, kjn, ktrd
    REAL :: ptrtrd(:, :, :)
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('trd_trc', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'trd_trc : You should not have seen this print! error?', ptrtrd(1, 1, 1)
    WRITE(*, FMT = *) '  "      "      : You should not have seen this print! error?', kjn, ktrd, kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE trd_trc
END MODULE trdtrc