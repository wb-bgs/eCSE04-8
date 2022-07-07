ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine dy2mjd
c
c       V. Lesur 16/06/06
c
c tranform a date given in decimal year in Julian days from 00:00 
c on 01/01/2000
c
c subroutine called: IDAYSINYR
c
c input:
c         dy          decimal year (real*8)
c
c output:
c         djd         julian day (real*8)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine dy2mjd(dy,djd)
c
        implicit real*8 (a-h,o-z)
        implicit integer*4 (i-n)
c
        real*8 dy,djd
c
        integer iy,nd
        real*8 ddiy
c
c year and day number in year
        iy=idint(dy)
        ddiy=dble(IDAYSINYR(iy))
c
c number of days from 01/01/2000 at 00:00 to 01/01/iy at 00:00
        nd=0
        if(iy.ge.2000)then
          do i=2000,iy-1
            nd=nd+IDAYSINYR(i) 
          enddo
        else
          do i=iy,1999
             nd=nd+IDAYSINYR(i)
          enddo
          nd=-nd
        endif
c
c add days, hours, minutes and seconds 
        djd=dble(nd)+(dy-dble(iy))*ddiy
c
        return
        end

