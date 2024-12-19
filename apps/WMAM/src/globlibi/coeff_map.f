cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	module coeff_map
c
c   The spherical harmonic coefficients are held within a single array (bc) of size l(l+2),
c   where l is the maximum spherical harmonic degree.
c
c   Previously, this cofficient array was accessed non-sequentially.
c
c   The first loop accessed coefficients at positions, 1^2, 2^2, 3^2, ..., l^2.
c   The second loop visited those same positions offset by one.
c   Subsequent loops then proceeded through the following offsets, 3, 5, 7, ..., 2(l-1) with
c   the number of iterations being one less than the previous loop. 
c   
c   This access pattern cannot take advantage of cach locality.
c   And so this module creates a coefficient map, allowing one to work with a scratch coefficient
c   array that can be accessed sequentially by the code in the "XYZsph_bi0.f" source file.
c
c   The coefficient map is used to re-arrange the terms within a coefficient array such that these
c   can be accessed sequentially within the XYZsph_bi0_sample() and XYZsph_bi0_sub() subroutines
c   and XYZsph_bi0_fun() function. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        module coeff_map
c
          implicit none
c
          public
     >    create_coeff_map, destroy_coeff_map,
     >    sequentialise, desequentialise
c
          private
          integer :: shd  = 0
          integer :: ncfs = 0
          integer :: cf_map_len = 0
          integer, allocatable :: cf_map(:)
          real*8, allocatable :: buf(:)
c
c
        contains
c
c
          subroutine create_coeff_map(shdeg)
c
            integer shdeg
            integer i, j, k
c
            if (0 .ne. cf_map_len) then
              call destroy_coeff_map()
            endif
c
            shd  = shdeg
            ncfs = shd*(shd + 2)
            cf_map_len = shd + int((ncfs - shd)/2)
c
            allocate(cf_map(1:cf_map_len))
c
            do i = 1,shd
              cf_map(i) = i**2
            enddo
c  
            i = shd+1
            do j = 1,shd
              do k = j,shd
                cf_map(i) = (k**2 + 2*(j-1)) + 1
                i = i+1
              enddo
            enddo
c
            allocate(buf(1:ncfs))
c
          end subroutine create_coeff_map
c
c
          subroutine destroy_coeff_map()
c
            if (0 .ne. cf_map_len) then
              deallocate(buf)
              deallocate(cf_map)
              cf_map_len = 0
              shd  = 0
              ncfs = 0
            endif
c
          end subroutine destroy_coeff_map
c
c
          subroutine sequentialise(coeffs)
c
            real*8  coeffs(1:ncfs)
            integer i, j, nu
c
            buf(1:ncfs) = coeffs(1:ncfs)
c
            do i = 1,shd
              nu = cf_map(i)
              coeffs(i) = buf(nu)
            enddo
c
            j = shd + 1
            do i = shd+1,cf_map_len
              nu = cf_map(i)
              coeffs(j) = buf(nu)
              coeffs(j+1) = buf(nu+1)
              j = j + 2
            enddo
c
          end subroutine sequentialise
c
c
          subroutine desequentialise(coeffs) 
c
            real*8  coeffs(1:ncfs)
            integer i, j, nu
c
            buf(1:ncfs) = coeffs(1:ncfs)
c
            do i = 1,shd
              nu = cf_map(i)
              coeffs(nu) = buf(i)
            enddo
c
            j = shd + 1
            do i = shd+1,cf_map_len
              nu = cf_map(i)
              coeffs(nu) = buf(j)
              coeffs(nu+1) = buf(j+1)
              j = j + 2
            enddo
c
          end subroutine desequentialise
c
c
        end module coeff_map