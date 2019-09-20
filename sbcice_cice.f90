MODULE sbcice_cice
  CONTAINS
  SUBROUTINE sbc_ice_cice(kt, ksbc)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: kt, ksbc
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('sbc_ice_cice', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'sbc_ice_cice: You should not have seen this print! error?', kt
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE sbc_ice_cice
  SUBROUTINE cice_sbc_init(ksbc)
    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: ksbc
    WRITE(*, FMT = *) 'cice_sbc_init: You should not have seen this print! error?', ksbc
  END SUBROUTINE cice_sbc_init
  SUBROUTINE cice_sbc_final
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    IMPLICIT NONE
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('cice_sbc_final', 'r0', psy_profile0)
    WRITE(*, FMT = *) 'cice_sbc_final: You should not have seen this print! error?'
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE cice_sbc_final
END MODULE sbcice_cice