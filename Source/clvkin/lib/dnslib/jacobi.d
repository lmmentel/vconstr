      subroutine jacobi(a,v,d,n,nrot)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension a(n,n),d(n),v(n,n)
      dimension b(n),z(n)
c
      v(1:n,1:n)=0.d0
      do ip=1,n
        v(ip,ip)=1.d0
      enddo
c
      do ip=1,n
        b(ip)=a(ip,ip)
      enddo
      d(1:n)=b(1:n)
      z(1:n)=0.d0
c
      nrot=0
      do i=1,50
        sm=0.d0
        do ip=1,n-1
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
          enddo
        enddo
        if (sm.eq.0.d0) return
        if (i.lt.4) then
          thresh=0.2*sm/n**2
        else
          thresh=0.d0
        endif
        do ip=1,n-1
          do iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if ((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     + .and.(abs(d(iq))+g.eq.abs(d(iq)))) then
              a(ip,iq)=0.d0
            elseif (abs(a(ip,iq)).gt.thresh) then
              h=d(iq)-d(ip)
              if (abs(h)+g.eq.abs(h)) then
                t=a(ip,iq)/h
              else
                theta=5.d-1*h/a(ip,iq)
                t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                if (theta.lt.0.d0) t=-t
              endif
              c=1.d0/sqrt(1.d0+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              enddo
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              enddo
              nrot=nrot+1
            endif
          enddo
        enddo
        b(1:n)=b(1:n)+z(1:n)
        d(1:n)=b(1:n)
        z(1:n)=0
      enddo
      write(6,'(''ERROR; too many iterations in jacobi'')')
c
      return
      end
