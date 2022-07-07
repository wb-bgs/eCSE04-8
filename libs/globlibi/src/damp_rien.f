cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               Subroutine damp_rien
c                       V. Lesur 27/05/2002
c
c             compute the damping matrix "Dt.R.D",
c               scale it by dl(1)  and add it to the normal matrix GG
c               scale ir by dl(2) and substract "Dt.R.D.bc" to BB
c             Assume a time variation of 2 years
c
c      input
c
c      output
c        GG = GG + dl(1) * Dt.R.D
c        BB = BB - dl(2) * Dt.R.D * bc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine damp_rien(dl,nb,GG,BB,bc)
c
         integer nb
         real*8 dl(*),gg(nb,*),bc(*),BB(*)
c
         return
         end
