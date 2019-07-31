MODULE iom_nf90
  USE dom_oce
  USE lbclnk
  USE iom_def
  USE netcdf
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iom_nf90_open, iom_nf90_close, iom_nf90_varid, iom_nf90_get, iom_nf90_gettime, iom_nf90_rstput
  PUBLIC :: iom_nf90_getatt, iom_nf90_putatt
  INTERFACE iom_nf90_get
    MODULE PROCEDURE iom_nf90_g0d, iom_nf90_g123d
  END INTERFACE
  INTERFACE iom_nf90_getatt
    MODULE PROCEDURE iom_nf90_giatt, iom_nf90_gratt, iom_nf90_gcatt
  END INTERFACE
  INTERFACE iom_nf90_putatt
    MODULE PROCEDURE iom_nf90_piatt, iom_nf90_pratt, iom_nf90_pcatt
  END INTERFACE
  INTERFACE iom_nf90_rstput
    MODULE PROCEDURE iom_nf90_rp0123d
  END INTERFACE
  CONTAINS
  SUBROUTINE iom_nf90_open(cdname, kiomid, ldwrt, ldok, kdompar, kdlev)
    CHARACTER(LEN = *), INTENT(INOUT) :: cdname
    INTEGER, INTENT(  OUT) :: kiomid
    LOGICAL, INTENT(IN   ) :: ldwrt
    LOGICAL, INTENT(IN   ) :: ldok
    INTEGER, DIMENSION(2, 5), INTENT(IN   ), OPTIONAL :: kdompar
    INTEGER, INTENT(IN   ), OPTIONAL :: kdlev
    CHARACTER(LEN = 256) :: clinfo
    CHARACTER(LEN = 256) :: cltmp
    INTEGER :: iln
    INTEGER :: istop
    INTEGER :: if90id
    INTEGER :: idmy
    INTEGER :: jl
    INTEGER :: ichunk
    INTEGER :: imode
    INTEGER :: ihdf5
    LOGICAL :: llclobber
    INTEGER :: ilevels
    clinfo = '                    iom_nf90_open ~~~  '
    istop = nstop
    IF (PRESENT(kdlev)) THEN
      ilevels = kdlev
    ELSE
      ilevels = jpk
    END IF
    IF (nn_chunksz > 0) THEN
      ichunk = nn_chunksz
    ELSE
      ichunk = nf90_sizehint_default
    END IF
    llclobber = ldwrt .AND. ln_clobber
    IF (ldok .AND. .NOT. llclobber) THEN
      IF (ldwrt) THEN
        IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' open existing file: ' // TRIM(cdname) // ' in WRITE mode'
        IF (snc4set % luse) THEN
          CALL iom_nf90_check(NF90_OPEN(TRIM(cdname), NF90_WRITE, if90id), clinfo)
        ELSE
          CALL iom_nf90_check(NF90_OPEN(TRIM(cdname), NF90_WRITE, if90id, chunksize = ichunk), clinfo)
        END IF
        CALL iom_nf90_check(NF90_SET_FILL(if90id, NF90_NOFILL, idmy), clinfo)
      ELSE
        IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' open existing file: ' // TRIM(cdname) // ' in READ mode'
        CALL iom_nf90_check(NF90_OPEN(TRIM(cdname), NF90_NOWRITE, if90id, chunksize = ichunk), clinfo)
      END IF
    ELSE
      iln = INDEX(cdname, '.nc')
      IF (ldwrt) THEN
        IF (jpnij > 1) THEN
          WRITE(cltmp, FMT = '(a,a,i4.4,a)') cdname(1 : iln - 1), '_', narea - 1, '.nc'
          cdname = TRIM(cltmp)
        END IF
        IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' create new file: ' // TRIM(cdname) // ' in WRITE mode'
        IF (llclobber) THEN
          imode = IOR(nf90_64bit_offset, nf90_clobber)
        ELSE
          imode = IOR(nf90_64bit_offset, nf90_noclobber)
        END IF
        IF (snc4set % luse) THEN
          IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' creating file: ' // TRIM(cdname) // ' in hdf5 (netcdf4) mode'
          CALL GET_NF90_SYMBOL("NF90_HDF5", ihdf5)
          IF (llclobber) THEN
            imode = IOR(ihdf5, nf90_clobber)
          ELSE
            imode = IOR(ihdf5, nf90_noclobber)
          END IF
          CALL iom_nf90_check(NF90_CREATE(TRIM(cdname), imode, if90id), clinfo)
        ELSE
          CALL iom_nf90_check(NF90_CREATE(TRIM(cdname), imode, if90id, chunksize = ichunk), clinfo)
        END IF
        CALL iom_nf90_check(NF90_SET_FILL(if90id, NF90_NOFILL, idmy), clinfo)
        CALL iom_nf90_check(NF90_DEF_DIM(if90id, 'x', kdompar(1, 1), idmy), clinfo)
        CALL iom_nf90_check(NF90_DEF_DIM(if90id, 'y', kdompar(2, 1), idmy), clinfo)
        IF (PRESENT(kdlev)) THEN
          CALL iom_nf90_check(NF90_DEF_DIM(if90id, 'numcat', ilevels, idmy), clinfo)
        ELSE
          CALL iom_nf90_check(NF90_DEF_DIM(if90id, 'nav_lev', ilevels, idmy), clinfo)
        END IF
        CALL iom_nf90_check(NF90_DEF_DIM(if90id, 'time_counter', NF90_UNLIMITED, idmy), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_number_total', jpnij), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_number', narea - 1), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_dimensions_ids', (/1, 2/)), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_size_global', (/jpiglo, jpjglo/)), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_size_local', kdompar(:, 1)), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_position_first', kdompar(:, 2)), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_position_last', kdompar(:, 3)), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_halo_size_start', kdompar(:, 4)), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_halo_size_end', kdompar(:, 5)), clinfo)
        CALL iom_nf90_check(NF90_PUT_ATT(if90id, NF90_GLOBAL, 'DOMAIN_type', 'BOX'), clinfo)
      ELSE
        CALL ctl_stop(TRIM(clinfo), ' should be impossible case...')
      END IF
    END IF
    IF (istop == nstop) THEN
      kiomid = 0
      DO jl = jpmax_files, 1, - 1
        IF (iom_file(jl) % nfid == 0) kiomid = jl
      END DO
      iom_file(kiomid) % name = TRIM(cdname)
      iom_file(kiomid) % nfid = if90id
      iom_file(kiomid) % iolib = jpnf90
      iom_file(kiomid) % nvars = 0
      iom_file(kiomid) % irec = - 1
      iom_file(kiomid) % nlev = ilevels
      CALL iom_nf90_check(NF90_Inquire(if90id, unlimitedDimId = iom_file(kiomid) % iduld), clinfo)
      IF (iom_file(kiomid) % iduld .GE. 0) THEN
        CALL iom_nf90_check(NF90_Inquire_Dimension(if90id, iom_file(kiomid) % iduld, name = iom_file(kiomid) % uldname, len = iom_file(kiomid) % lenuld), clinfo)
      END IF
      IF (lwp) WRITE(numout, FMT = *) '                   ---> ' // TRIM(cdname) // ' OK'
    ELSE
      kiomid = 0
    END IF
  END SUBROUTINE iom_nf90_open
  SUBROUTINE iom_nf90_close(kiomid)
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = 100) :: clinfo
    clinfo = '      iom_nf90_close    , file: ' // TRIM(iom_file(kiomid) % name)
    CALL iom_nf90_check(NF90_CLOSE(iom_file(kiomid) % nfid), clinfo)
  END SUBROUTINE iom_nf90_close
  FUNCTION iom_nf90_varid(kiomid, cdvar, kiv, kdimsz, kndims)
    INTEGER, INTENT(IN   ) :: kiomid
    CHARACTER(LEN = *), INTENT(IN   ) :: cdvar
    INTEGER, INTENT(IN   ) :: kiv
    INTEGER, DIMENSION(:), INTENT(  OUT), OPTIONAL :: kdimsz
    INTEGER, INTENT(  OUT), OPTIONAL :: kndims
    INTEGER :: iom_nf90_varid
    INTEGER :: if90id
    INTEGER :: ji
    INTEGER :: ivarid
    INTEGER :: i_nvd
    INTEGER, DIMENSION(jpmax_dims) :: idimid
    LOGICAL :: llok
    CHARACTER(LEN = 100) :: clinfo
    clinfo = '          iom_nf90_varid, file: ' // TRIM(iom_file(kiomid) % name) // ', var: ' // TRIM(cdvar)
    iom_nf90_varid = 0
    IF (PRESENT(kdimsz)) kdimsz(:) = 0
    if90id = iom_file(kiomid) % nfid
    llok = NF90_INQ_VARID(if90id, TRIM(cdvar), ivarid) == nf90_noerr
    IF (llok) THEN
      iom_nf90_varid = kiv
      iom_file(kiomid) % nvars = kiv
      iom_file(kiomid) % nvid(kiv) = ivarid
      iom_file(kiomid) % cn_var(kiv) = TRIM(cdvar)
      CALL iom_nf90_check(NF90_Inquire_Variable(if90id, ivarid, ndims = i_nvd), clinfo)
      iom_file(kiomid) % ndims(kiv) = i_nvd
      CALL iom_nf90_check(NF90_Inquire_Variable(if90id, ivarid, dimids = idimid(1 : i_nvd)), clinfo)
      iom_file(kiomid) % luld(kiv) = .FALSE.
      iom_file(kiomid) % dimsz(:, kiv) = 0
      DO ji = 1, i_nvd
        CALL iom_nf90_check(NF90_Inquire_Dimension(if90id, idimid(ji), len = iom_file(kiomid) % dimsz(ji, kiv)), clinfo)
        IF (idimid(ji) == iom_file(kiomid) % iduld) iom_file(kiomid) % luld(kiv) = .TRUE.
      END DO
      llok = NF90_Inquire_attribute(if90id, ivarid, 'scale_factor') == nf90_noerr
      IF (llok) THEN
        CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, 'scale_factor', iom_file(kiomid) % scf(kiv)), clinfo)
      ELSE
        iom_file(kiomid) % scf(kiv) = 1.
      END IF
      llok = NF90_Inquire_attribute(if90id, ivarid, 'add_offset') == nf90_noerr
      IF (llok) THEN
        CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, 'add_offset', iom_file(kiomid) % ofs(kiv)), clinfo)
      ELSE
        iom_file(kiomid) % ofs(kiv) = 0.
      END IF
      IF (PRESENT(kdimsz)) THEN
        IF (i_nvd <= SIZE(kdimsz)) THEN
          kdimsz(1 : i_nvd) = iom_file(kiomid) % dimsz(1 : i_nvd, kiv)
        ELSE
          WRITE(ctmp1, FMT = *) i_nvd, SIZE(kdimsz)
          CALL ctl_stop(TRIM(clinfo), 'error in kdimsz size' // TRIM(ctmp1))
        END IF
      END IF
      IF (PRESENT(kndims)) kndims = iom_file(kiomid) % ndims(kiv)
    ELSE
      iom_nf90_varid = - 1
    END IF
  END FUNCTION iom_nf90_varid
  SUBROUTINE iom_nf90_g0d(kiomid, kvid, pvar, kstart)
    INTEGER, INTENT(IN   ) :: kiomid
    INTEGER, INTENT(IN   ) :: kvid
    REAL(KIND = wp), INTENT(  OUT) :: pvar
    INTEGER, DIMENSION(1), INTENT(IN   ), OPTIONAL :: kstart
    CHARACTER(LEN = 100) :: clinfo
    clinfo = 'iom_nf90_g0d , file: ' // TRIM(iom_file(kiomid) % name) // ', var: ' // TRIM(iom_file(kiomid) % cn_var(kvid))
    CALL iom_nf90_check(NF90_GET_VAR(iom_file(kiomid) % nfid, iom_file(kiomid) % nvid(kvid), pvar, start = kstart), clinfo)
  END SUBROUTINE iom_nf90_g0d
  SUBROUTINE iom_nf90_g123d(kiomid, kvid, knbdim, kstart, kcount, kx1, kx2, ky1, ky2, pv_r1d, pv_r2d, pv_r3d)
    INTEGER, INTENT(IN   ) :: kiomid
    INTEGER, INTENT(IN   ) :: kvid
    INTEGER, INTENT(IN   ) :: knbdim
    INTEGER, DIMENSION(:), INTENT(IN   ) :: kstart
    INTEGER, DIMENSION(:), INTENT(IN   ) :: kcount
    INTEGER, INTENT(IN   ) :: kx1, kx2, ky1, ky2
    REAL(KIND = wp), DIMENSION(:), INTENT(  OUT), OPTIONAL :: pv_r1d
    REAL(KIND = wp), DIMENSION(:, :), INTENT(  OUT), OPTIONAL :: pv_r2d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(  OUT), OPTIONAL :: pv_r3d
    CHARACTER(LEN = 100) :: clinfo
    INTEGER :: if90id
    INTEGER :: ivid
    clinfo = 'iom_nf90_g123d , file: ' // TRIM(iom_file(kiomid) % name) // ', var: ' // TRIM(iom_file(kiomid) % cn_var(kvid))
    if90id = iom_file(kiomid) % nfid
    ivid = iom_file(kiomid) % nvid(kvid)
    IF (PRESENT(pv_r1d)) THEN
      CALL iom_nf90_check(NF90_GET_VAR(if90id, ivid, pv_r1d(:), start = kstart(1 : knbdim), count = kcount(1 : knbdim)), clinfo)
    ELSE IF (PRESENT(pv_r2d)) THEN
      CALL iom_nf90_check(NF90_GET_VAR(if90id, ivid, pv_r2d(kx1 : kx2, ky1 : ky2), start = kstart(1 : knbdim), count = kcount(1 : knbdim)), clinfo)
    ELSE IF (PRESENT(pv_r3d)) THEN
      CALL iom_nf90_check(NF90_GET_VAR(if90id, ivid, pv_r3d(kx1 : kx2, ky1 : ky2, :), start = kstart(1 : knbdim), count = kcount(1 : knbdim)), clinfo)
    END IF
  END SUBROUTINE iom_nf90_g123d
  SUBROUTINE iom_nf90_giatt(kiomid, cdatt, pv_i0d, cdvar)
    INTEGER, INTENT(IN   ) :: kiomid
    CHARACTER(LEN = *), INTENT(IN   ) :: cdatt
    INTEGER, INTENT(  OUT) :: pv_i0d
    CHARACTER(LEN = *), INTENT(IN   ), OPTIONAL :: cdvar
    INTEGER :: if90id
    INTEGER :: ivarid
    LOGICAL :: llok
    CHARACTER(LEN = 100) :: clinfo
    if90id = iom_file(kiomid) % nfid
    IF (PRESENT(cdvar)) THEN
      llok = NF90_INQ_VARID(if90id, TRIM(cdvar), ivarid) == nf90_noerr
      IF (llok) THEN
        llok = NF90_Inquire_attribute(if90id, ivarid, cdatt) == nf90_noerr
      ELSE
        CALL ctl_warn('iom_nf90_getatt: no variable ' // cdvar // ' found')
      END IF
    ELSE
      llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt) == nf90_noerr
      ivarid = NF90_GLOBAL
    END IF
    IF (llok) THEN
      clinfo = 'iom_nf90_getatt, file: ' // TRIM(iom_file(kiomid) % name) // ', giatt: ' // TRIM(cdatt)
      CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values = pv_i0d), clinfo)
    ELSE
      CALL ctl_warn('iom_nf90_getatt: no attribute ' // cdatt // ' found')
      pv_i0d = - 999
    END IF
  END SUBROUTINE iom_nf90_giatt
  SUBROUTINE iom_nf90_gratt(kiomid, cdatt, pv_r0d, cdvar)
    INTEGER, INTENT(IN   ) :: kiomid
    CHARACTER(LEN = *), INTENT(IN   ) :: cdatt
    REAL(KIND = wp), INTENT(  OUT) :: pv_r0d
    CHARACTER(LEN = *), OPTIONAL, INTENT(IN   ) :: cdvar
    INTEGER :: if90id
    INTEGER :: ivarid
    LOGICAL :: llok
    CHARACTER(LEN = 100) :: clinfo
    if90id = iom_file(kiomid) % nfid
    IF (PRESENT(cdvar)) THEN
      llok = NF90_INQ_VARID(if90id, TRIM(cdvar), ivarid) == nf90_noerr
      IF (llok) THEN
        llok = NF90_Inquire_attribute(if90id, ivarid, cdatt) == nf90_noerr
      ELSE
        CALL ctl_warn('iom_nf90_getatt: no variable ' // cdvar // ' found')
      END IF
    ELSE
      llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt) == nf90_noerr
      ivarid = NF90_GLOBAL
    END IF
    IF (llok) THEN
      clinfo = 'iom_nf90_getatt, file: ' // TRIM(iom_file(kiomid) % name) // ', gratt: ' // TRIM(cdatt)
      CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values = pv_r0d), clinfo)
    ELSE
      CALL ctl_warn('iom_nf90_getatt: no attribute ' // cdatt // ' found')
      pv_r0d = - 999._wp
    END IF
  END SUBROUTINE iom_nf90_gratt
  SUBROUTINE iom_nf90_gcatt(kiomid, cdatt, pv_c0d, cdvar)
    INTEGER, INTENT(IN   ) :: kiomid
    CHARACTER(LEN = *), INTENT(IN   ) :: cdatt
    CHARACTER(LEN = *), INTENT(  OUT) :: pv_c0d
    CHARACTER(LEN = *), OPTIONAL, INTENT(IN   ) :: cdvar
    INTEGER :: if90id
    INTEGER :: ivarid
    LOGICAL :: llok
    CHARACTER(LEN = 100) :: clinfo
    if90id = iom_file(kiomid) % nfid
    IF (PRESENT(cdvar)) THEN
      llok = NF90_INQ_VARID(if90id, TRIM(cdvar), ivarid) == nf90_noerr
      IF (llok) THEN
        llok = NF90_Inquire_attribute(if90id, ivarid, cdatt) == nf90_noerr
      ELSE
        CALL ctl_warn('iom_nf90_getatt: no variable ' // cdvar // ' found')
      END IF
    ELSE
      llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt) == nf90_noerr
      ivarid = NF90_GLOBAL
    END IF
    IF (llok) THEN
      clinfo = 'iom_nf90_getatt, file: ' // TRIM(iom_file(kiomid) % name) // ', gcatt: ' // TRIM(cdatt)
      CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values = pv_c0d), clinfo)
    ELSE
      CALL ctl_warn('iom_nf90_getatt: no attribute ' // cdatt // ' found')
      pv_c0d = '!'
    END IF
  END SUBROUTINE iom_nf90_gcatt
  SUBROUTINE iom_nf90_piatt(kiomid, cdatt, pv_i0d, cdvar)
    INTEGER, INTENT(IN   ) :: kiomid
    CHARACTER(LEN = *), INTENT(IN   ) :: cdatt
    INTEGER, INTENT(IN   ) :: pv_i0d
    CHARACTER(LEN = *), INTENT(IN   ), OPTIONAL :: cdvar
    INTEGER :: if90id
    INTEGER :: ivarid
    LOGICAL :: llok
    LOGICAL :: lenddef
    CHARACTER(LEN = 100) :: clinfo
    if90id = iom_file(kiomid) % nfid
    lenddef = .FALSE.
    IF (PRESENT(cdvar)) THEN
      llok = NF90_INQ_VARID(if90id, TRIM(cdvar), ivarid) == nf90_noerr
      IF (.NOT. llok) THEN
        CALL ctl_warn('iom_nf90_putatt: no variable ' // cdvar // ' found')
      END IF
    ELSE
      llok = .TRUE.
      ivarid = NF90_GLOBAL
    END IF
    IF (llok) THEN
      clinfo = 'iom_nf90_putatt, file: ' // TRIM(iom_file(kiomid) % name) // ', piatt: ' // TRIM(cdatt)
      IF (iom_file(kiomid) % irec /= - 1) THEN
        CALL iom_nf90_check(NF90_REDEF(if90id), clinfo)
        lenddef = .TRUE.
      END IF
      CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values = pv_i0d), clinfo)
      IF (lenddef) THEN
        CALL iom_nf90_check(NF90_ENDDEF(if90id), clinfo)
      END IF
    ELSE
      CALL ctl_warn('iom_nf90_putatt: no attribute ' // cdatt // ' written')
    END IF
  END SUBROUTINE iom_nf90_piatt
  SUBROUTINE iom_nf90_pratt(kiomid, cdatt, pv_r0d, cdvar)
    INTEGER, INTENT(IN   ) :: kiomid
    CHARACTER(LEN = *), INTENT(IN   ) :: cdatt
    REAL(KIND = wp), INTENT(IN   ) :: pv_r0d
    CHARACTER(LEN = *), OPTIONAL, INTENT(IN   ) :: cdvar
    INTEGER :: if90id
    INTEGER :: ivarid
    LOGICAL :: llok
    LOGICAL :: lenddef
    CHARACTER(LEN = 100) :: clinfo
    if90id = iom_file(kiomid) % nfid
    lenddef = .FALSE.
    IF (PRESENT(cdvar)) THEN
      llok = NF90_INQ_VARID(if90id, TRIM(cdvar), ivarid) == nf90_noerr
      IF (.NOT. llok) THEN
        CALL ctl_warn('iom_nf90_putatt: no variable ' // cdvar // ' found')
      END IF
    ELSE
      llok = .TRUE.
      ivarid = NF90_GLOBAL
    END IF
    IF (llok) THEN
      clinfo = 'iom_nf90_putatt, file: ' // TRIM(iom_file(kiomid) % name) // ', pratt: ' // TRIM(cdatt)
      IF (iom_file(kiomid) % irec /= - 1) THEN
        CALL iom_nf90_check(NF90_REDEF(if90id), clinfo)
        lenddef = .TRUE.
      END IF
      CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values = pv_r0d), clinfo)
      IF (lenddef) THEN
        CALL iom_nf90_check(NF90_ENDDEF(if90id), clinfo)
      END IF
    ELSE
      CALL ctl_warn('iom_nf90_putatt: no attribute ' // cdatt // ' written')
    END IF
  END SUBROUTINE iom_nf90_pratt
  SUBROUTINE iom_nf90_pcatt(kiomid, cdatt, pv_c0d, cdvar)
    INTEGER, INTENT(IN   ) :: kiomid
    CHARACTER(LEN = *), INTENT(IN   ) :: cdatt
    CHARACTER(LEN = *), INTENT(IN   ) :: pv_c0d
    CHARACTER(LEN = *), OPTIONAL, INTENT(IN   ) :: cdvar
    INTEGER :: if90id
    INTEGER :: ivarid
    LOGICAL :: llok
    LOGICAL :: lenddef
    CHARACTER(LEN = 100) :: clinfo
    if90id = iom_file(kiomid) % nfid
    lenddef = .FALSE.
    IF (PRESENT(cdvar)) THEN
      llok = NF90_INQ_VARID(if90id, TRIM(cdvar), ivarid) == nf90_noerr
      IF (.NOT. llok) THEN
        CALL ctl_warn('iom_nf90_putatt: no variable ' // cdvar // ' found')
      END IF
    ELSE
      llok = .TRUE.
      ivarid = NF90_GLOBAL
    END IF
    IF (llok) THEN
      clinfo = 'iom_nf90_putatt, file: ' // TRIM(iom_file(kiomid) % name) // ', pcatt: ' // TRIM(cdatt)
      IF (iom_file(kiomid) % irec /= - 1) THEN
        CALL iom_nf90_check(NF90_REDEF(if90id), clinfo)
        lenddef = .TRUE.
      END IF
      CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values = pv_c0d), clinfo)
      IF (lenddef) THEN
        CALL iom_nf90_check(NF90_ENDDEF(if90id), clinfo)
      END IF
    ELSE
      CALL ctl_warn('iom_nf90_putatt: no attribute ' // cdatt // ' written')
    END IF
  END SUBROUTINE iom_nf90_pcatt
  SUBROUTINE iom_nf90_gettime(kiomid, kvid, ptime, cdunits, cdcalendar)
    INTEGER, INTENT(IN   ) :: kiomid
    INTEGER, INTENT(IN   ) :: kvid
    REAL(KIND = wp), DIMENSION(:), INTENT(  OUT) :: ptime
    CHARACTER(LEN = *), OPTIONAL, INTENT(  OUT) :: cdunits
    CHARACTER(LEN = *), OPTIONAL, INTENT(  OUT) :: cdcalendar
    CHARACTER(LEN = 100) :: clinfo
    clinfo = 'iom_nf90_gettime, file: ' // TRIM(iom_file(kiomid) % name) // ', var: ' // TRIM(iom_file(kiomid) % cn_var(kvid))
    CALL iom_nf90_check(NF90_GET_VAR(iom_file(kiomid) % nfid, iom_file(kiomid) % nvid(kvid), ptime(:), start = (/1/), count = (/iom_file(kiomid) % dimsz(1, kvid)/)), clinfo)
    IF (PRESENT(cdunits)) THEN
      CALL iom_nf90_check(NF90_GET_ATT(iom_file(kiomid) % nfid, iom_file(kiomid) % nvid(kvid), "units", values = cdunits), clinfo)
    END IF
    IF (PRESENT(cdcalendar)) THEN
      CALL iom_nf90_check(NF90_GET_ATT(iom_file(kiomid) % nfid, iom_file(kiomid) % nvid(kvid), "calendar", values = cdcalendar), clinfo)
    END IF
  END SUBROUTINE iom_nf90_gettime
  SUBROUTINE iom_nf90_rp0123d(kt, kwrite, kiomid, cdvar, kvid, ktype, pv_r0d, pv_r1d, pv_r2d, pv_r3d)
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kwrite
    INTEGER, INTENT(IN) :: kiomid
    CHARACTER(LEN = *), INTENT(IN) :: cdvar
    INTEGER, INTENT(IN) :: kvid
    INTEGER, INTENT(IN), OPTIONAL :: ktype
    REAL(KIND = wp), INTENT(IN), OPTIONAL :: pv_r0d
    REAL(KIND = wp), DIMENSION(:), INTENT(IN), OPTIONAL :: pv_r1d
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN), OPTIONAL :: pv_r2d
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN), OPTIONAL :: pv_r3d
    INTEGER :: idims
    INTEGER :: idvar
    INTEGER :: jd
    INTEGER :: ix1, ix2, iy1, iy2
    INTEGER, DIMENSION(4) :: idimsz
    INTEGER, DIMENSION(4) :: idimid
    CHARACTER(LEN = 256) :: clinfo
    CHARACTER(LEN = 12), DIMENSION(4) :: cltmp
    INTEGER :: if90id
    INTEGER :: idmy
    INTEGER :: itype
    INTEGER, DIMENSION(4) :: ichunksz
    INTEGER :: ichunkalg, ishuffle, ideflate, ideflate_level
    LOGICAL :: lchunk
    INTEGER :: idlv
    clinfo = '          iom_nf90_rp0123d, file: ' // TRIM(iom_file(kiomid) % name) // ', var: ' // TRIM(cdvar)
    if90id = iom_file(kiomid) % nfid
    IF (iom_file(kiomid) % nvars == 0) THEN
      IF (iom_file(kiomid) % irec /= - 1) THEN
        CALL iom_nf90_check(nf90_redef(if90id), clinfo)
        iom_file(kiomid) % irec = - 1
      END IF
      IF (iom_file(kiomid) % nlev == jpk) THEN
        cltmp = (/'nav_lon     ', 'nav_lat     ', 'nav_lev     ', 'time_counter'/)
      ELSE
        cltmp = (/'nav_lon     ', 'nav_lat     ', 'numcat      ', 'time_counter'/)
      END IF
      CALL iom_nf90_check(NF90_DEF_VAR(if90id, TRIM(cltmp(1)), NF90_FLOAT, (/1, 2/), iom_file(kiomid) % nvid(1)), clinfo)
      CALL iom_nf90_check(NF90_DEF_VAR(if90id, TRIM(cltmp(2)), NF90_FLOAT, (/1, 2/), iom_file(kiomid) % nvid(2)), clinfo)
      CALL iom_nf90_check(NF90_DEF_VAR(if90id, TRIM(cltmp(3)), NF90_FLOAT, (/3/), iom_file(kiomid) % nvid(3)), clinfo)
      CALL iom_nf90_check(NF90_DEF_VAR(if90id, TRIM(cltmp(4)), NF90_DOUBLE, (/4/), iom_file(kiomid) % nvid(4)), clinfo)
      iom_file(kiomid) % nvars = 4
      iom_file(kiomid) % luld(1 : 4) = (/.FALSE., .FALSE., .FALSE., .TRUE./)
      iom_file(kiomid) % cn_var(1 : 4) = cltmp
      iom_file(kiomid) % ndims(1 : 4) = (/2, 2, 1, 1/)
      iom_file(kiomid) % dimsz(1, 1) = 0
      IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' define dimension variables done'
    END IF
    IF (kvid <= 0) THEN
      ichunkalg = 0
      ishuffle = 1
      ideflate = 1
      ideflate_level = 1
      idvar = iom_file(kiomid) % nvars + 1
      IF (iom_file(kiomid) % irec /= - 1) THEN
        CALL iom_nf90_check(nf90_redef(if90id), clinfo)
        iom_file(kiomid) % irec = - 1
      END IF
      IF (PRESENT(pv_r0d)) THEN
        idims = 0
      ELSE IF (PRESENT(pv_r1d)) THEN
        idims = 2
        idimid(1 : idims) = (/3, 4/)
      ELSE IF (PRESENT(pv_r2d)) THEN
        idims = 3
        idimid(1 : idims) = (/1, 2, 4/)
      ELSE IF (PRESENT(pv_r3d)) THEN
        idims = 4
        idimid(1 : idims) = (/1, 2, 3, 4/)
      END IF
      IF (PRESENT(ktype)) THEN
        SELECT CASE (ktype)
        CASE (jp_r8)
          itype = nf90_double
        CASE (jp_r4)
          itype = nf90_float
        CASE (jp_i4)
          itype = nf90_int
        CASE (jp_i2)
          itype = nf90_short
        CASE (jp_i1)
          itype = nf90_byte
        CASE DEFAULT
          CALL ctl_stop(TRIM(clinfo) // ' unknown variable type')
        END SELECT
      ELSE
        itype = NF90_DOUBLE
      END IF
      IF (PRESENT(pv_r0d)) THEN
        CALL iom_nf90_check(NF90_DEF_VAR(if90id, TRIM(cdvar), itype, iom_file(kiomid) % nvid(idvar)), clinfo)
      ELSE
        CALL iom_nf90_check(NF90_DEF_VAR(if90id, TRIM(cdvar), itype, idimid(1 : idims), iom_file(kiomid) % nvid(idvar)), clinfo)
      END IF
      lchunk = .FALSE.
      IF (snc4set % luse .AND. idims == 4) lchunk = .TRUE.
      iom_file(kiomid) % nvars = idvar
      iom_file(kiomid) % cn_var(idvar) = TRIM(cdvar)
      iom_file(kiomid) % scf(idvar) = 1.
      iom_file(kiomid) % ofs(idvar) = 0.
      iom_file(kiomid) % ndims(idvar) = idims
      IF (.NOT. PRESENT(pv_r0d)) THEN
        iom_file(kiomid) % luld(idvar) = .TRUE.
      ELSE
        iom_file(kiomid) % luld(idvar) = .FALSE.
      END IF
      DO jd = 1, idims
        CALL iom_nf90_check(NF90_INQUIRE_DIMENSION(if90id, idimid(jd), len = iom_file(kiomid) % dimsz(jd, idvar)), clinfo)
        IF (lchunk) ichunksz(jd) = iom_file(kiomid) % dimsz(jd, idvar)
      END DO
      IF (lchunk) THEN
        ichunksz(1) = MIN(ichunksz(1), MAX((ichunksz(1) - 1) / snc4set % ni + 1, 16))
        ichunksz(2) = MIN(ichunksz(2), MAX((ichunksz(2) - 1) / snc4set % nj + 1, 16))
        ichunksz(3) = MIN(ichunksz(3), MAX((ichunksz(3) - 1) / snc4set % nk + 1, 1))
        ichunksz(4) = 1
        CALL iom_nf90_check(SET_NF90_DEF_VAR_CHUNKING(if90id, idvar, ichunkalg, ichunksz), clinfo)
        CALL iom_nf90_check(SET_NF90_DEF_VAR_DEFLATE(if90id, idvar, ishuffle, ideflate, ideflate_level), clinfo)
        IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' chunked ok. Chunks sizes: ', ichunksz
      END IF
      IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' defined ok'
    ELSE
      idvar = kvid
    END IF
    IF (kt == kwrite) THEN
      IF (iom_file(kiomid) % irec == - 1) THEN
        CALL iom_nf90_check(nf90_enddef(if90id), clinfo)
        iom_file(kiomid) % irec = 0
      END IF
      IF (PRESENT(pv_r2d) .OR. PRESENT(pv_r3d)) THEN
        idimsz(1 : 2) = iom_file(kiomid) % dimsz(1 : 2, idvar)
        IF (idimsz(1) == (nlei - nldi + 1) .AND. idimsz(2) == (nlej - nldj + 1)) THEN
          ix1 = nldi
          ix2 = nlei
          iy1 = nldj
          iy2 = nlej
        ELSE IF (idimsz(1) == nlci .AND. idimsz(2) == nlcj) THEN
          ix1 = 1
          ix2 = nlci
          iy1 = 1
          iy2 = nlcj
        ELSE IF (idimsz(1) == jpi .AND. idimsz(2) == jpj) THEN
          ix1 = 1
          ix2 = jpi
          iy1 = 1
          iy2 = jpj
        ELSE
          CALL ctl_stop('iom_nf90_rp0123d: should have been an impossible case...')
        END IF
        IF (iom_file(kiomid) % dimsz(1, 1) == 0) THEN
          CALL iom_nf90_check(NF90_INQ_VARID(if90id, 'nav_lon', idmy), clinfo)
          CALL iom_nf90_check(NF90_PUT_VAR(if90id, idmy, glamt(ix1 : ix2, iy1 : iy2)), clinfo)
          CALL iom_nf90_check(NF90_INQ_VARID(if90id, 'nav_lat', idmy), clinfo)
          CALL iom_nf90_check(NF90_PUT_VAR(if90id, idmy, gphit(ix1 : ix2, iy1 : iy2)), clinfo)
          IF (iom_file(kiomid) % nlev == jpk) THEN
            CALL iom_nf90_check(NF90_INQ_VARID(if90id, 'nav_lev', idmy), clinfo)
            CALL iom_nf90_check(NF90_PUT_VAR(if90id, idmy, gdept_1d), clinfo)
          ELSE
            CALL iom_nf90_check(NF90_INQ_VARID(if90id, 'numcat', idmy), clinfo)
            CALL iom_nf90_check(NF90_PUT_VAR(if90id, idmy, (/(idlv, idlv = 1, iom_file(kiomid) % nlev)/)), clinfo)
          END IF
          CALL iom_nf90_check(NF90_INQ_VARID(if90id, 'time_counter', idmy), clinfo)
          CALL iom_nf90_check(NF90_PUT_VAR(if90id, idmy, kt), clinfo)
          CALL iom_nf90_check(NF90_INQUIRE_DIMENSION(if90id, 1, len = iom_file(kiomid) % dimsz(1, 1)), clinfo)
          CALL iom_nf90_check(NF90_INQUIRE_DIMENSION(if90id, 2, len = iom_file(kiomid) % dimsz(2, 1)), clinfo)
          iom_file(kiomid) % dimsz(1 : 2, 2) = iom_file(kiomid) % dimsz(1 : 2, 1)
          CALL iom_nf90_check(NF90_INQUIRE_DIMENSION(if90id, 3, len = iom_file(kiomid) % dimsz(1, 3)), clinfo)
          iom_file(kiomid) % dimsz(1, 4) = 1
          IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' write dimension variables done'
        END IF
      END IF
      IF (PRESENT(pv_r0d)) THEN
        CALL iom_nf90_check(NF90_PUT_VAR(if90id, idvar, pv_r0d), clinfo)
      ELSE IF (PRESENT(pv_r1d)) THEN
        CALL iom_nf90_check(NF90_PUT_VAR(if90id, idvar, pv_r1d(:)), clinfo)
      ELSE IF (PRESENT(pv_r2d)) THEN
        CALL iom_nf90_check(NF90_PUT_VAR(if90id, idvar, pv_r2d(ix1 : ix2, iy1 : iy2)), clinfo)
      ELSE IF (PRESENT(pv_r3d)) THEN
        CALL iom_nf90_check(NF90_PUT_VAR(if90id, idvar, pv_r3d(ix1 : ix2, iy1 : iy2, :)), clinfo)
      END IF
      IF (iom_file(kiomid) % luld(idvar)) iom_file(kiomid) % dimsz(iom_file(kiomid) % ndims(idvar), idvar) = iom_file(kiomid) % dimsz(iom_file(kiomid) % ndims(idvar), idvar) + 1
      IF (lwp) WRITE(numout, FMT = *) TRIM(clinfo) // ' written ok'
    END IF
  END SUBROUTINE iom_nf90_rp0123d
  SUBROUTINE iom_nf90_check(kstatus, cdinfo)
    INTEGER, INTENT(IN) :: kstatus
    CHARACTER(LEN = *), INTENT(IN) :: cdinfo
    IF (kstatus /= nf90_noerr) CALL ctl_stop('iom_nf90_check : ' // TRIM(nf90_strerror(kstatus)), TRIM(cdinfo))
  END SUBROUTINE iom_nf90_check
END MODULE iom_nf90