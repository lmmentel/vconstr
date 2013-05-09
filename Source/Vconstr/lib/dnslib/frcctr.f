      subroutine frcctr(nq,eorb,qorb,occorb,dqmax,parin,epse,epsq,kkk)
c ======================================================================
c
c purpose:  Charge transfer as described by Averil and Painter
c
c input  :  nq - nr. of levels participating in charge transfer
c           qorb() - occution numbers
c           eorb() - orbital energies
c output :  qorb() - the new occupation numbers
c
c remarks:  See article : F. W. Averil and G. S. Painter,
c                         Phys. Rev. B, vol. 46, 1992, page 2498
c
c *=====================================================================
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter (eps   = 1.0d-8)
c
      dimension eorb(nq),qorb(nq)
c
      logical lmask
      dimension lmask(nq)
      dimension gamma(nq),qorbin(nq)
c
c** sum of orbital energies
c
      esum = sum(eorb)
      esum = esum-nq*eorb(kkk)
c
c** determine set of parameters (gamma)
c
      kmin = kkk
      kmax = kkk
      qorbin(1:nq) = qorb(1:nq)
      do iq = 1,nq
        if (iq.ne.kkk) then
          de = eorb(iq)-eorb(kkk)
          if (de.gt.0.d0) then
            gamma(iq) = qorbin(iq)/de
          elseif (de.lt.0.d0) then
            gamma(iq) = (qorbin(iq)-occorb)/de
          endif
          if ((abs(de).lt.epse).and.
     + (abs(qorbin(iq)-qorbin(kkk)).lt.epsq)) then
            kmin = min(iq,kmin)
            kmax = max(iq,kmax)
          endif
        else
          if (esum.gt.0.d0) then
            gamma(kkk) = (occorb-qorbin(kkk))/esum
          else
            gamma(kkk) = -qorbin(kkk)/esum
          endif
        endif
      enddo
c
c** recalculate sum and gamma(kkk) excluding states with gamma = 0
c** also a possibly degeneracy of accepting (kkk=1) or donating (kkk=nq) 
c** level is taken care of
c
      ndeg = kmax-kmin+1
      gamma(kkk) = gamma(kkk)*esum*ndeg
      esum = 0.d0
      do iq = 1,nq
        if ((iq.lt.kmin).or.(iq.gt.kmax)) then
          if (gamma(iq).gt.eps) then
            esum = esum+eorb(iq)-eorb(kkk)
          endif
        endif
      enddo
      if (abs(esum).lt.eps) then
        write(6,'(''WARNING; No states left in frcctr'')')
        return
      endif
c
      rdum = gamma(kkk)/esum
      gamma(kmin:kmax) = rdum
c
c** Minimal value of gamma
c
      tmp = parin
      do iq = 1,nq
        if (gamma(iq).gt.eps) then
          tmp = min(tmp,gamma(iq))
        endif
      enddo
c
c** Charge transfer
c
      tmp = min(tmp,dqmax/abs(esum))
c
      qorb(1:nq) = qorbin(1:nq)
      do iq = 1,nq
        if (gamma(iq).gt.eps) then
          if ((iq.lt.kmin).or.(iq.gt.kmax)) then
            qorb(iq) = qorbin(iq)-tmp*(eorb(iq)-eorb(kkk))
          endif
        endif
      enddo
c
      do k = kmin,kmax
        qorb(k) = qorbin(k)+tmp*esum/ndeg
      enddo
c
      do iq = 1,nq
        if (qorb(iq).lt.0.d0) then
          qorb(iq)=0.d0
        endif
      enddo
c
c** check for degeneracy
c
      do iq = 1,nq-1
        ndeg = 0
        lmask(1:nq) = .false.
        do jq = iq,nq
          if ((abs(eorb(jq)-eorb(iq)).lt.epse).and.
     + (abs(qorbin(jq)-qorbin(iq)).lt.epsq)) then
            ndeg = ndeg+1
            lmask(jq) = .true.
          endif
        enddo
c
c** divide charge equally
c
        qdeg = sum(qorb,1,lmask)
        do jq = iq,nq
          if (lmask(jq)) qorb(jq) = qdeg/ndeg
        enddo
c
      enddo
c
      return
      end
