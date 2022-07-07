ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine write_comp
c		V. Lesur	07/03/2002
c
c       write a file of nl data line.
c       The line before the first data line will be '##########'
c       previous lines must start with a '#' and may contain comments
c       data lines have nd+3 real records:
c       point position (nd real), data, fit, (data-fit) 
c       The data point printed are selected by type
c
c       input:
c         nom		file name to write
c         nl            number of data line to write
c         npm           number max of data point
c         nd            number of space dimmension
c         nt(*)         data type (dim nl)
c         ts            type selected
c         ddat(nd+1,*)  position and data (dim (nd+1,npm))      
c         xyzf(*)       fit to data (dim nl)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine write_comp(nom,nl,npm,nd,nt,ts,ddat,xyzf)
c
        implicit none
c
	character nom*100,theformat*15
        integer i,j,nl,npm,nd,nt(*),ts,nb,nd1,nd2
        real*8 ddat(*),xyzf(*),adif,sdif
c
	open(10,file=nom,status='unknown')
c
        write(10,'(A)')'##########'
c
        nb=0
        adif=0.0d0
        sdif=0.0d0
        nd1=int((nd+3)/10.)
        nd2=(nd+3)-nd1*10
        theformat='('//char(nd2+48)//'e17.9)'
        if(nd1.ne.0)theformat='('//char(nd1+48)//char(nd2+48)//'e17.9)'
        nd1=nd+1
c
        i=0
        do while (i.lt.nl)
          i=i+1
          if(nt(i).eq.ts)then
            nb=nb+1
            write(10,theformat)(ddat(j+(i-1)*nd1),j=1,nd1)
     >                ,xyzf(i),ddat(i*nd1)-xyzf(i)
            adif=adif+ddat(i*nd1)-xyzf(i)
            sdif=sdif+(ddat(i*nd1)-xyzf(i))**2
          endif
        enddo
c
        close(10)
c
        adif=adif/nb
        sdif=dsqrt((sdif-nb*adif**2)/nb)
c
        write(*,'(2(A,f16.7))')'Residual average : ',adif
     >                        ,' with L2 std : ',sdif
c
        return
        end
