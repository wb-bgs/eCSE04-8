# Coefficient Map

Initialising an index map for sequentialising the coefficient array.


```fortran
subroutine create_coeff_map(deg)

  ...
  
  do i = 1,deg
    cf_map(i) = i**2
  enddo
  
  i = deg+1
  do j = 1,deg
    do k = j,deg
      cf_map(i) = (k**2 + 2*(j-1)) + 1
      i = i+1
    enddo
  enddo

  ...

end subroutine create_coeff_map
```
