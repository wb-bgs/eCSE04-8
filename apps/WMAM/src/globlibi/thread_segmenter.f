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
c      nranks INTEGER number of threads available (from MPI)
c      nelems INTEGER array length in question
c
c   output:
c      proc_np(*) INTEGER  block lengths (dim min = size)
c      proc_ip(*) INTEGER  data indexes (dim min = size)
c
c   remarks:
c      1- no checks at 'proc_np' and 'proc_ip', content overwritten,
c         correct dimension etc. is assumed.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine thread_segmenter(nranks, nelems, proc_np, proc_ip)
c
        implicit none
c
        integer ,intent(in)    :: nranks, nelems
        integer ,intent(out)   :: proc_np(1:nranks), proc_ip(1:nranks)
c
        integer                :: i,k,n
c
        proc_np(1:nranks) = 0 
        proc_ip(1:nranks) = 0
c
c Define segment size & remaining
        k = int(nelems/nranks)
        n = mod(nelems,nranks)
c
c Init by identical elements of width k
        proc_np(1:nranks) = k       
c
c Deal with remaining part
        if (n.gt.0) proc_np(1:n) = proc_np(1:n)+1
c
c Set indexes
        proc_ip(1) = 1
        do i = 2,nranks
          proc_ip(i) = proc_ip(i-1) + proc_np(i-1) 
        enddo
c
        return
        end subroutine thread_segmenter