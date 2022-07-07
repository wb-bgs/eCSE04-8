program gd2gc_cmdln
! IN:
!  cltgd geodetic colatitude (deg)
!  alt   altitude above sea-level (km)
!  xgd   north magnetic component in geodetic coords
!  zgd   vertical magnetic component in geodetic coords
!
! OUT:
!  cltgc geocentic colatitude (deg)
!  r     radius from centre of earth (km)
!  xgc   north magnetic component in geocentric coords
!  zgc   vertical magnetic component in geocentric coords
  implicit none
  real(kind=8) :: cltgd, alt, xgd, zgd, cltgc,r,xgc,zgc
  integer :: ioStatus

  write(*,*) "  cltgd      alt       xgd       zgd  |   cltgc        r       xgc       zgc"
  do
    read (*,*,iostat=ioStatus) cltgd, alt, xgd, zgd
    if (ioStatus.ne.0) exit
    call gd2gc (cltgd,alt,xgd,zgd,cltgc,r,xgc,zgc)
    write(*,'(F8.3,A1,F8.3,A1,F9.2,A1,F9.2,A3,F8.3,A1,F8.3,A1,F9.2,A1,F9.2)') &
     cltgd,' ',alt,' ',xgd,' ',zgd,'  |',cltgc, ' ',r,' ',xgc,' ',zgc
  end do

end program gd2gc_cmdln
