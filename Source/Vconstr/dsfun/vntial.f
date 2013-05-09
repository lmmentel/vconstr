      subroutine vntial(npnt,dns,ddns,vxc,alpha,beta,gamma)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(xconst = -1.240700981799d0,
     +          restrc =  0.7937005259841d0,
     +          third  = 1d0 / 3d0)
c
      dimension dns(npnt),ddns(npnt,3),vxc(npnt)
c
      cxalph = alpha*xconst*restrc
c
      do ipnt=1,npnt
        rh=dns(ipnt)
        absgrd=sqrt(ddns(ipnt,1)**2+ddns(ipnt,2)**2+
     + ddns(ipnt,3)**2)
        vxc(ipnt) = cxalph*rh**third+
     + beta*epsvwn(rh)+gamma*epsbek(rh,absgrd)
      enddo
c
      return
      end
