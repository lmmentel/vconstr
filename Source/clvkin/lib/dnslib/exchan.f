      subroutine exchan(p,norb,px,valmo)
c*******************************************************************
c calculate exchange-hole-density : px(2|1) = p(1,2)p(2,1) / p(1)
c
c density stored in array p in triangle form : p(i*(i-1)/2+j)= pij + pji
c
c !!!!  assumed that alpha-spin-density and beta-spin-density are equal
c
c arguments:
c p  : density matrix in mo basis
c norb  : dimension of density matrix
c x,y,z  : position 1 of reference electron
c px  : density of exchange-hole (returned)
c valmo : array to hold values of mos in ref. pos.
c*******************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension px(norb*(norb+1)/2)
      dimension p(norb*(norb+1)/2),valmo(norb)
      dimension thelp(norb)
c
      px(1:norb*(norb+1)/2)=0.d0
      thelp(1:norb)=0.d0
c
      do i=1,norb
        do j=1,norb
          if (i.ge.j) then
            ij=i*(i-1)/2+j
          else
            ij=j*(j-1)/2+i
          endif
          fac=5.d-1
          if (i.eq.j) fac=1.d0
          thelp(i)=thelp(i)+fac*p(ij)*valmo(j)
        enddo
      enddo
c
      do j=1,norb
        do k=1,j-1
          jk=j*(j-1)/2+k
          px(jk)=thelp(j)*thelp(k)*2.d0
        enddo
        jk=j*(j+1)/2
        px(jk)=thelp(j)*thelp(j)
      enddo
c
      dens = valmat(p,norb,valmo,.false.)
      factor = -1.d0/dens/2.d0
      px(1:norb*(norb+1)/2)=px(1:norb*(norb+1)/2)*factor
c
      return
      end
