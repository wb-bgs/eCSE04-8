cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p2
c		Vincent Lesur 09/02/2005
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals_p2
c
c       input:
c         nd            Space dimension
c         nlocpts       number of data+sampling points local to rank
c         nlocdatpts    number of data points local to rank
c         shdeg         max SH degree value
c         d2a           pre-computed array for mk_lf_dlf()
c         ppos(nd+1,*)  point position in ndD
c         nb            Number of base functions
c         bc            base coefficients
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p2(nd, nlocpts, nlocdatpts, shdeg,
     >                             d2a, ppos, nb, bc, xyzf)
c
        implicit none
c
        include 'mpif.h'
c
        integer :: nd, nb, nlocpts, nlocdatpts, shdeg
        real*8  :: d2a(0:shdeg), ppos(nd+1,nlocpts), bc(nb)
        real*8  :: xyzf(nlocpts)
c
        real*8, parameter :: RAG = 6371.2d0
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c        
        integer i
        real*8 p1, p2, ra
        real*8 bex, bey, bez
c
        real*8 XYZsph_bi0_fun
c
#if defined(OMP_OFFLOAD_CPTP2)
        logical, save :: firstcall = .TRUE.
#endif
c 
c
        integer ierr, nranks, rank 
        integer j, k
        real*8 d2a_sum, bc_sum, ppos_sum
        real*8 xyzf_sum
        call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
c        
#if defined(OMP_OFFLOAD_CPTP2)
!$OMP TARGET DATA if(firstcall)
!$omp& map(to: nb, nd)
!$omp& map(to: nlocpts, nlocdatpts, shdeg)
!$omp& map(to: d2a(0:shdeg))
!$omp& map(to: ppos(1:nd+1,1:nlocpts))
        if (firstcall) then
          firstcall = .FALSE.
        endif
c
!$OMP TARGET DATA
!$omp& map(to: bc(1:nb))
#endif
c
c
#if defined(OMP_OFFLOAD_CPTP2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!$omp& map(from: xyzf(1:nlocdatpts))
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
!$omp& schedule(static)
        do i = 1,nlocdatpts
c
          p1 = ppos(1,i)*D2R
          p2 = ppos(2,i)*D2R
          ra = RAG / ppos(3,i)
c
          bex = ppos(5,i)
          bey = ppos(6,i)
          bez = ppos(7,i)
c
          xyzf(i) = XYZsph_bi0_fun(shdeg, nb,
     >                             d2a, bc,
     >                             p1, p2, ra,
     >                             bex, bey, bez)
c
        enddo
#if defined(OMP_OFFLOAD_CPTP2)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP END PARALLEL DO
#endif
c
c
#if defined(OMP_OFFLOAD_CPTP2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!$omp& map(from: xyzf(nlocdatpts+1:nlocpts))
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
!$omp& schedule(static)
        do i = nlocdatpts+1,nlocpts
c
          p1 = ppos(1,i)*D2R
          p2 = ppos(2,i)*D2R
          ra = RAG / ppos(3,i)
c
          bex = ppos(5,i)
          bey = ppos(6,i)
          bez = ppos(7,i)
c
          call XYZsph_bi0_sample(shdeg, nb,
     >                           d2a, bc,
     >                           p1, p2, ra, 
     >                           bex, bey, bez)
c
          xyzf(i) = XYZsph_bi0_fun(shdeg, nb,
     >                             d2a, bc,
     >                             p1, p2, ra,
     >                             bex, bey, bez)
c
        enddo
#if defined(OMP_OFFLOAD_CPTP2)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
!$OMP END TARGET DATA
!$OMP END TARGET DATA
#else
!$OMP END PARALLEL DO
#endif
c
c
        do i = 0,nranks-1
          if (i .eq. rank) then
            d2a_sum = 0.d0
            do j = 0,shdeg
              d2a_sum = d2a_sum + d2a(j)
            enddo
            bc_sum = 0.d0
            do j = 1,nb
              bc_sum = bc_sum + bc(j)
            enddo
            ppos_sum = 0.d0
            do j = 1,nlocpts
              do k = 1,3
                ppos_sum = ppos_sum + ppos(k,j)
              enddo
              do k = 5,7
                ppos_sum = ppos_sum + ppos(k,j)
              enddo
            enddo 
            xyzf_sum = 0.d0
            do j = 1,nlocpts
              xyzf_sum = xyzf_sum + xyzf(j)
            enddo
            write(*,*) 'cpt_dat_vals_p2: ', rank,
     >          ': ', nd, ' ', nb,
     >          ' ', nlocpts, ' ', nlocdatpts,
     >          ' ', shdeg, ' ', d2a_sum,
     >          ' ', bc_sum, ' ', ppos_sum,
     >          ' ', xyzf_sum
          endif
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        enddo

        return
        end