# Sampled dot product

The dot products performed for each sampled point following a call to `XYZsph_bi0()`, wherein the `dw` arrays are populated.
  

```fortran
subroutine sub_sph_wmam_l(iof, nb, bc, bp, be)

  ...

  if (sampled_point) then  

    ...

    call XYZsph_bi0(deg, rag, bp, dwx, dwy, dwz)
      
    dx = dot_product(dwx(1:nb), bc(1:nb))
    dy = dot_product(dwy(1:nb), bc(1:nb))
    dz = dot_product(dwz(1:nb), bc(1:nb))

    ...

  endif

  ...

end subroutine sub_sph_wmam_l
```
