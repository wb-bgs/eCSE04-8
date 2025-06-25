# Extra operations on `dw` arrays

The additional operations requiring the `dw` arrays. The `be` array assignments are given in an abbreviated form compared to the actual code.


```fortran
subroutine sub_sph_wmam_l(iof, nb, bc, bp, be)

  
  ...

  
  if (sampled_point) then  
  
    ...
  
    call XYZsph_bi0(deg, rag, bp, dwx, dwy, dwz)
  
    ...
  
    be(1:nb) = (dwx(1:nb)*bex
             +  dwy(1:nb)*bey
             +  dwz(1:nb)*bez) / dd            
  
  else
  
    ! input point
  
    ...
  
    call XYZsph_bi0(deg, rag, bp, dwx, dwy, dwz)
 
     be(1:nb) = (dwx(1:nb)*bp(5)
              +  dwy(1:nb)*bp(6)
              +  dwz(1:nb)*bp(7)) / dd         
  
  endif


  if (iof .eq. 'f') then
    bedotbc = dot_product(be(1:nb), bc(1:nb))
  endif     

  
  ...

  
end subroutine sub_sph_wmam_l
```
