      function epsbek(dnsty,grd)
c
c     Becke's nonlocal exchange energy density
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(zero= 0.0d0)
      parameter(one = 1.0d0)
      parameter(two = 2.0d0)
      parameter(three = 3.0d0)
      parameter(six = 6.0d0)
      parameter(frthrd=4.0d0/3.0d0)
      parameter(third =1.0d0/3.0d0)
      parameter(bbecke=0.0042d0)
      parameter(tol = 1d-3)
c
      halfgr = abs(grd/two)
      halfrh = dnsty/two
      rho13  = halfrh**third
      rho43  = halfrh**frthrd
      if (rho43.gt.1.d-320) then
        x = halfgr/rho43
        xsq = x**2
        sqrt1x = sqrt(one+xsq)
        tmp = x+sqrt1x
        if (abs(tmp-one).gt.tol) then
          sinhix = log(tmp)
        else
          tmp = tmp-one
          sinhix = tmp-tmp**2/two+tmp**3/three
        endif
        fbecke = xsq/(one+six*bbecke*x*sinhix)
        epsbek=-bbecke*fbecke*rho13
      else
        epsbek=zero
      endif
c
      return
      end
c
      function vwnpotential(dnsty)
c
c     This function calculates the VWN-potential
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(xp    = -0.10498d0)
      parameter(bp    = 3.72744d0)
      parameter(cp    = 12.9352d0)
      parameter(qp    = 6.15199d0)
      parameter(cap   = 0.0310907d0)
      parameter(cbp   = -0.0311676d0)
      parameter(ccp   = 1.24742d0)
c
      epscp=hvwn(dnsty,xp,bp,cp,qp,cap,cbp,ccp)
      depscp=rhdh(dnsty,xp,bp,cp,qp,cap,cbp,ccp)
      vwn=epscp + depscp
      return
      end
c
      function epsvwn(dnsty)
c
c     VWN-correlation energy density
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(xp    = -0.10498d0)
      parameter(bp    = 3.72744d0)
      parameter(cp    = 12.9352d0)
      parameter(qp    = 6.15199d0)
      parameter(cap   = 0.0310907d0)
      parameter(cbp   = -0.0311676d0)
      parameter(ccp   = 1.24742d0)
c
      epsvwn=hvwn(dnsty,xp,bp,cp,qp,cap,cbp,ccp)
      return
      end
c

      function hvwn(dnsty,x,b,c,q,ca,cb,cc)
c
c     This function is called by the VWN-routines
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(one   =1.0d0)
      parameter(two   =2.0d0)
      parameter(three =3.0d0)
      parameter(four  =4.0d0)
      parameter(third =1.0d0/3.0d0)
      parameter(sixth =1.0d0/6.0d0)
c
      fourpi=four*acos(-one)
      g=(three/fourpi)**third
      sg=sqrt(g)
c
      d3=dnsty**third
      d6=dnsty**sixth
c
      hvwn=ca*(log(g)-two*cb*log(sg-x*d6)+
     +  (cb-one)*log(g+b*sg*d6+c*d3)+
     +  cc*atan(q*d6/(two*sg+b*d6)))
c
      return
      end

      function rhdh(dnsty,x,b,c,q,ca,cb,cc)
c
c     This function is called by the VWN-potential routine
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(one   =1.0d0)
      parameter(two   =2.0d0)
      parameter(three =3.0d0)
      parameter(four  =4.0d0)
      parameter(third =1.0d0/3.0d0)
      parameter(sixth =1.0d0/6.0d0)
c
      fourpi=four*acos(-one)
      g=(three/fourpi)**third
      sg=sqrt(g)
c
      d3=dnsty**third
      d6=dnsty**sixth
c
      rhdh=ca*(third*cb*x*d6/(sg-x*d6) +
     +  (cb-one)*sixth*(b*sg*d6+two*c*d3)/(g+b*sg*d6+c*d3)
     +  + third*cc*q*sg*d6/((two*sg+b*d6)**two+q*q*d3))
c
      return
      end
