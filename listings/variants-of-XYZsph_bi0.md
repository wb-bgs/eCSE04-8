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
