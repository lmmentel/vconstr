      subroutine valgrd(xp1,yp1,zp1,pmo,norb,value, 
     +                  valuex,valuey,valuez)
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension pmo(norb*(norb+1)/2),value(norb), 
     +     valuex(norb),valuey(norb),valuez(norb)
c
c*** diagonal
c
      xp1 = 0.d0
      yp1 = 0.d0
      zp1 = 0.d0
      do 20 i=1,norb
        ii=i*(i-1)/2
        xp1=xp1+2.d0*pmo(ii+i)*value(i)*valuex(i)
        yp1=yp1+2.d0*pmo(ii+i)*value(i)*valuey(i)
        zp1=zp1+2.d0*pmo(ii+i)*value(i)*valuez(i)
c
c*** now off-diagonal
c
        do 30 j=1,i-1
          xp1=xp1+pmo(ii+j)*(value(i)*valuex(j)
     +        +value(j)*valuex(i))
          yp1=yp1+pmo(ii+j)*(value(i)*valuey(j)
     +        +value(j)*valuey(i))
          zp1=zp1+pmo(ii+j)*(value(i)*valuez(j)
     +        +value(j)*valuez(i))
 30     continue
 20   continue
      return
      end
c
      function valsgrd(pmo,val,grdx,grdy,grdz,grds,norb)
c***********************************************************************
c
c Calculate the Laplacian of the density matrix, where pmo is given in 
c   HF-MO basis
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension pmo(norb*(norb+1)/2),val(norb),grdx(norb),grdy(norb),
     +  grdz(norb),grds(norb)
c
      sk=0.d0
      do 110 i=1,norb
        ii=i*(i-1)/2
        sk=sk+2*pmo(ii+i)*(grdx(i)**2+grdy(i)**2+grdz(i)**2+
     + val(i)*grds(i))
        do 120 j=1,i-1
          sk=sk+pmo(ii+j)*(2*(grdx(i)*grdx(j)+grdy(i)*grdy(j)+
     +      grdz(i)*grdz(j))+val(i)*grds(j)+grds(i)*val(j))
  120   continue
  110 continue
      valsgrd=sk
c
      return
      end

      subroutine val2grd(dxx,dxy,dxz,dyy,dyz,dzz,pmo,norb,val,valx,
     + valy,valz,valxx,valxy,valxz,valyy,valyz,valzz)
c***********************************************************************
c
c Calculate the Laplacian of the density matrix, where pmo is given in 
c   HF-MO basis
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension pmo(norb*(norb+1)/2),val(norb),valx(norb),valy(norb),
     + valz(norb),valxx(norb),valxy(norb),valxz(norb),valyy(norb),
     + valyz(norb),valzz(norb)
c
      dxx=0.d0
      dxy=0.d0
      dxz=0.d0
      dyy=0.d0
      dyz=0.d0
      dzz=0.d0
      do 110 i=1,norb
        ii=i*(i-1)/2
        dxx=dxx+2*pmo(ii+i)*(valx(i)*valx(i)+val(i)*valxx(i))
        dxy=dxy+2*pmo(ii+i)*(valx(i)*valy(i)+val(i)*valxy(i))
        dxz=dxz+2*pmo(ii+i)*(valx(i)*valz(i)+val(i)*valxz(i))
        dyy=dyy+2*pmo(ii+i)*(valy(i)*valy(i)+val(i)*valyy(i))
        dyz=dyz+2*pmo(ii+i)*(valy(i)*valz(i)+val(i)*valyz(i))
        dzz=dzz+2*pmo(ii+i)*(valz(i)*valz(i)+val(i)*valzz(i))
        do 120 j=1,i-1
          dxx=dxx+pmo(ii+j)*(2*valx(i)*valx(j)+val(i)*valxx(j)+
     + valxx(i)*val(j))
          dxy=dxy+pmo(ii+j)*(valx(i)*valy(j)+valy(i)*valx(j)+
     + val(i)*valxy(j)+valxy(i)*val(j))
          dxz=dxz+pmo(ii+j)*(valx(i)*valz(j)+valz(i)*valx(j)+
     + val(i)*valxz(j)+valxz(i)*val(j))
          dyy=dyy+pmo(ii+j)*(2*valy(i)*valy(j)+val(i)*valyy(j)+
     + valyy(i)*val(j))
          dyz=dyz+pmo(ii+j)*(valy(i)*valz(j)+valz(i)*valy(j)+
     + val(i)*valyz(j)+valyz(i)*val(j))
          dzz=dzz+pmo(ii+j)*(2*valz(i)*valz(j)+val(i)*valzz(j)+
     + valzz(i)*val(j))
  120   continue
  110 continue
c
      return
      end
