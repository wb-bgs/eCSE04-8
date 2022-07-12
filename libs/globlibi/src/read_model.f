ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine read_model
c		V. Lesur	21/01/2003
c
c       13/07/04:   modified for refyear in Julian days & gotos removed
c
c       read a model file with a maximum of nl model coefficients
c
c       The first lines of the file may be used for comments.
c       These comment lines should all start with '# ' 
c       The line before the first data line starts with '##'.
c       The first line after should contain the reference date of the model
c       after, start the model and on each line:
c          a dummy character (usually g or h)
c          two dummy integer (usually l and m)
c          the model parameter value 
c          a dummy real (usually the parameter variance)
c       
c
c       input:
c         fname		file name to read
c         ncoeffs       Maximum number of line to read
c
c       output:
c         ryg           ref year for model in Julian days
c         ncoeffs       number of data line
c         bc            real array (ncoeffs)
c         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine read_model(fname, ryg, ncoeffs, bc)
c
        implicit none
c
	character fname*100
        integer ryg, ncoeffs
        real*8 bc(*)
c        
        character ligne*80, cd
        integer i, id
        real*8 ry, rd
c
	open(10, file=fname, status='old')
c
c  read header       
        ligne(1:2)='  '
        do while (ligne(1:2) .ne. '##')
          read(10,*) ligne
        enddo
c
c  read year
        read(10,*) ry
        call dy2mjd(ry, ryg)
c
c  read coefficients
        i=0
        do while (i .le. ncoeffs)
          i=i+1
          read(10,*,END=22222) cd,id,id,bc(i),rd
        enddo
c
        write(*,*) 'ERROR READING INPUT FILE (NL > NLMAX)'
        close(10)
        stop
c
22222   ncoeffs=i-1
        close(10)
c
        return
        end