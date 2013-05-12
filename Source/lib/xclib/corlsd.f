      subroutine corlsd (np, lnlpot, rs, zet,  nilrho,
     +                   eu, eurs, ep, eprs, alfm, alfrsm,
     +                   ec, ecrs, eczet)
c
c ======================================================================
c
c purpose: uniform-gas correlation of perdew and wang 1991
c
c input  : rs - seitz radius (rs)
c          zet - relative spin polarization
c
c output : ec - correlation energy per electron
c          ecrs, eczet - derivatives of ec wrt rs & zet
c
c scratch: eu, eurs, ep, eprs, alfm
c
c history: 1996.08.02, GtV: error in second call GCOR(.false. rather
c          than LNLPOT), yielding incorrect results for spin-polac
c          ALFRSM
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    lnlpot, nilrho
c
      parameter (zero = 0d0,      one   = 1d0,
     +           two  = 2d0,      four  = 4d0,
     +           thrd = 1d0/3d0,  thrd4 = 4d0/3d0,
     +           tiny = 1d-30)
c
      dimension rs(np), zet(np), ec(np), ecrs(np),
     +          eczet(np), nilrho(np), eu(np), eurs(np), ep(np),
     +          eprs(np), alfm(np), alfrsm(np)
c
      data gam,fzz /0.5198421d0, 1.709921d0/
c
c ======================================================================
c
      glue1 = 0.0310907d0
      glue2 = 0.21370d0
      glue3 = 7.5957d0
      glue4 = 3.5876d0
      glue5 = 1.6382d0
      glue6 = 0.49294d0
      glue7 = 1.00d0
c
      call gcor (lnlpot, np,
     +           glue1, glue2, glue3, glue4, glue5, glue6, glue7,
     +           rs, nilrho, eu, eurs)
c
      glue11 = 0.01554535d0
      glue12 = 0.20548d0
      glue13 = 14.1189d0
      glue14 = 6.1977d0
      glue15 = 3.3662d0
      glue16 = 0.62517d0
      glue17 = 1.00d0
c
      call gcor (lnlpot, np,
     +           glue11, glue12, glue13, glue14, glue15, glue16, glue17,
     +           rs, nilrho ,ep, eprs)
c
      glue21 = 0.0168869d0
      glue22 = 0.11125d0
      glue23 = 10.357d0
      glue24 = 3.6231d0
      glue25 = 0.88026d0
      glue26 = 0.49671d0
      glue27 = 1.00d0
c
      call gcor (lnlpot, np,
     +           glue21, glue22, glue23, glue24, glue25, glue26, glue27,
     +           rs, nilrho, alfm, alfrsm)
c
      do 10 ip = 1, np
        alfm(ip) = alfm(ip) / fzz
   10 continue
      if (lnlpot) then
        do 20 ip = 1, np
          alfrsm(ip) = alfrsm(ip) / fzz
   20   continue
      endif
c
c
      do 100 ip = 1, np
        if (nilrho(ip)) goto 100
c
        tmp1 = zero
        s = one + zet(ip)
        if (s.gt.tiny) tmp1 = s ** thrd4
c
        tmp2 = zero
        s = one - zet(ip)
        if (s.gt.tiny) tmp2 = s ** thrd4
c
        f = (tmp1 + tmp2 - two) / gam
c
c       --------------------------------
c       alfm is minus the spin stiffness
c       --------------------------------
c
        z3 = zet(ip)**3
        z4 = zet(ip) * z3
c
        ec(ip) = eu(ip) - f*alfm(ip) +
     +           f*z4 * (ep(ip)-eu(ip)+alfm(ip))
c
        if (lnlpot) then
c
c         ------------------------------
c         energy done. now the potential
c         ------------------------------
c
          ecrs(ip) = eurs(ip) - f*alfrsm(ip) +
     +               f*z4 * (eprs(ip)-eurs(ip)+alfrsm(ip))
c
          tmp1 = zero
          s = one + zet(ip)
          if (s.gt.tiny) tmp1 = s**thrd
c
          tmp2 = zero
          s = one - zet(ip)
          if (s.gt.tiny) tmp2 = s**thrd
c
          fz = thrd4 * (tmp1 - tmp2) / gam
c
          eczet(ip) = z3 * (four*f - zet(ip)*fz) *
     +                     (ep(ip) - eu(ip) + alfm(ip)) -
     +                fz * alfm(ip)
        endif
c
 100  continue
c
c
      return
      end
