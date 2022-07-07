cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine thread_segmenter.f
c                  M.Rother 2010-04-09          
c
c   modified 12.04.2010 to usual Lesur's format (V.Lesur)
c
c   Given an array of size "nn", define a segmentation in "size" 
c   more or less equal parts. (replace the usual segmentation used 
c   in v.lesur's codes that can fail for small nn)
c
c   input:
c      size INTEGER number of threads available (from MPI)
c      nn   INTEGER array length in question
c
c   output:
c      proc_np(*) INTEGER  block lengths (dim min = size)
c      proc_ip(*) INTEGER  data pointers (dim min = size)
c
c   remarks:
c      1- no checks at 'proc_np' and 'proc_ip', content overwritten,
c         correct dimension etc. is assumed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine thread_segmenter(size,nn,proc_np,proc_ip)
c
        implicit none
c
        integer ,intent(in)    :: size,nn
        integer ,intent(out)   :: proc_np(*),proc_ip(*)
c
        integer                :: i,k,n
c
        proc_np(1:size)=0 
        proc_ip(1:size)=0
c
c Define segment size & remaining
        k=int(nn/size)
        n=mod(nn,size)
c
c Init by identical elements of width k
        proc_np(1:size)=k       
c
c Deal with remaining part
        if (n.gt.0) proc_np(1:n)=proc_np(1:n)+1
c
c Set integer pointers
        proc_ip(1)=1
        do i=2,size
          proc_ip(i)=proc_ip(i-1)+proc_np(i-1) 
        enddo
c
        return
        end subroutine thread_segmenter
