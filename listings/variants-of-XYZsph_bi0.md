# Variants of `XYZsph_bi0` subroutine

The calls to variant subroutines of `XYZsph_bi0()`.
  

```fortran
subroutine sub_sph_wmam_l(iof, nub, nb, bc, bp, be)
    

  ...
  

  if (sampled_point) then  
    
    call XYZsph_bi0_sample(deg, rag, bp, d2a,
                           bc, dx, dy, dz)
    ...

    bex = bex2
    bey = bey2
    bez = bez2

  else

    ! input point
    bex = bp[5]
    bey = bp[6]
    bez = bp[7]

  endif


  if (iof .eq. 'f') then 

    call XYZsph_bi0_fun(deg, rag, bp, d2a,
                        bc, bedotbc, bex, bey, bez)

  else

    call XYZsph_bi0_sub(deg, rag, bp, d2a,
                        be, bex, bey, bez)

  end if 


end subroutine sub_sph_wmam_l
```


Previously, the `dwx`, `dwy` and `dwz` arrays would be passed into the `XYZsph_bi0` subroutine to be populated and then
for the majority of calls only be used to calculate dot products.

We now call variants of the `XYZsph_bi0` subroutine, trapping the times when only dot products are needed and instead
accumulate these directly within the subroutine, returning them as simple scalar variables, e.g. `dx`, `dy` and `dz`. 


```fortran
subroutine XYZsph_bi0_sample(ilg, rag, pos, d2a, bc,
                             dx, dy, dz)
  
  ...

  real*8 dx,dy,dz

  ...

  do il=1,ilg
    dx = dx + bx * bc(nu) 
    dy = dy + by * bc(nu) 
    dz = dz + bz * bc(nu) 
  enddo

  ...

end subroutine XYZsph_bi0_sample
```