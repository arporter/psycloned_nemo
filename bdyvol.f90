MODULE bdyvol
  USE oce
  USE bdy_oce
  USE sbc_oce
  USE dom_oce
  USE phycst
  USE sbcisf
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_vol
  CONTAINS
  SUBROUTINE bdy_vol(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jb, jgrd
    INTEGER :: ib_bdy, ii, ij
    REAL(KIND = wp) :: zubtpecor, z_cflxemp, ztranst
    TYPE(OBC_INDEX), POINTER :: idx
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('bdy_vol', 'r0', psy_profile0)
    IF (ln_vol) THEN
      IF (kt == nit000) THEN
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'bdy_vol : Correction of velocities along unstructured OBC'
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
      END IF
      z_cflxemp = SUM((emp(:, :) - rnf(:, :) + fwfisf(:, :)) * bdytmask(:, :) * e1e2t(:, :)) / rau0
      IF (lk_mpp) CALL mpp_sum(z_cflxemp)
      zubtpecor = 0._wp
      DO ib_bdy = 1, nb_bdy
        idx => idx_bdy(ib_bdy)
        jgrd = 2
        DO jb = 1, idx % nblenrim(jgrd)
          DO jk = 1, jpkm1
            ii = idx % nbi(jb, jgrd)
            ij = idx % nbj(jb, jgrd)
            zubtpecor = zubtpecor + idx % flagu(jb, jgrd) * ua(ii, ij, jk) * e2u(ii, ij) * e3u_n(ii, ij, jk)
          END DO
        END DO
        jgrd = 3
        DO jb = 1, idx % nblenrim(jgrd)
          DO jk = 1, jpkm1
            ii = idx % nbi(jb, jgrd)
            ij = idx % nbj(jb, jgrd)
            zubtpecor = zubtpecor + idx % flagv(jb, jgrd) * va(ii, ij, jk) * e1v(ii, ij) * e3v_n(ii, ij, jk)
          END DO
        END DO
      END DO
      IF (lk_mpp) CALL mpp_sum(zubtpecor)
      IF (nn_volctl == 1) THEN
        zubtpecor = (zubtpecor - z_cflxemp) / bdysurftot
      ELSE
        zubtpecor = zubtpecor / bdysurftot
      END IF
      ztranst = 0._wp
      DO ib_bdy = 1, nb_bdy
        idx => idx_bdy(ib_bdy)
        jgrd = 2
        DO jb = 1, idx % nblenrim(jgrd)
          DO jk = 1, jpkm1
            ii = idx % nbi(jb, jgrd)
            ij = idx % nbj(jb, jgrd)
            ua(ii, ij, jk) = ua(ii, ij, jk) - idx % flagu(jb, jgrd) * zubtpecor * umask(ii, ij, jk)
            ztranst = ztranst + idx % flagu(jb, jgrd) * ua(ii, ij, jk) * e2u(ii, ij) * e3u_n(ii, ij, jk)
          END DO
        END DO
        jgrd = 3
        DO jb = 1, idx % nblenrim(jgrd)
          DO jk = 1, jpkm1
            ii = idx % nbi(jb, jgrd)
            ij = idx % nbj(jb, jgrd)
            va(ii, ij, jk) = va(ii, ij, jk) - idx % flagv(jb, jgrd) * zubtpecor * vmask(ii, ij, jk)
            ztranst = ztranst + idx % flagv(jb, jgrd) * va(ii, ij, jk) * e1v(ii, ij) * e3v_n(ii, ij, jk)
          END DO
        END DO
      END DO
      IF (lk_mpp) CALL mpp_sum(ztranst)
      IF (lwp .AND. MOD(kt, nwrite) == 0) THEN
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) 'bdy_vol : time step :', kt
        IF (lwp) WRITE(numout, FMT = *) '~~~~~~~ '
        IF (lwp) WRITE(numout, FMT = *) '          cumulate flux EMP             =', z_cflxemp, ' (m3/s)'
        IF (lwp) WRITE(numout, FMT = *) '          total lateral surface of OBC  =', bdysurftot, '(m2)'
        IF (lwp) WRITE(numout, FMT = *) '          correction velocity zubtpecor =', zubtpecor, '(m/s)'
        IF (lwp) WRITE(numout, FMT = *) '          cumulated transport ztranst   =', ztranst, '(m3/s)'
      END IF
    END IF
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE bdy_vol
END MODULE bdyvol