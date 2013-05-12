      block data io
c
      implicit real *8   (a-h,p-w), integer (i-n), logical (o)
      implicit character*8 (z), character*1 (x)
      implicit character*4 (y)
c
      logical oswit,orevis
c
c...  ynam : atmol-file name
c...  fname : unix-file_name (restricted to 44 chars (arbitrary)
c
      common/disc/isel,iselr,iselw,irep,ichek,ipos(16),iblksz(16),
     * ilen(16),nam(16)
      character*44 fname
      character*4 ynam
      common/discc/ynam(16),fname(16)
      common/work/jrec,jump,istrt(40),inumb(40),iwid,nend(80),
     * nline,noline,iwidth,ierr,oswit
      character * 80 ia,ja
      common/workc/ia,ja
      common/scra/i205(340,4)
      common/local/v(3680),ib(3680),
     *i(4,3680),j(4,3680),k(4,3680),l(4,3680)
      common/sector/num3,iblk3,orevis(2),
     * apos(204),maxb,iblkla
      common/iofile/iread,iwr,nav,oprint
      common/timess/timlim,date,time,accno,anam,isec
c
      data ynam/'ed0','ed1','ed2','ed3','ed4','ed5','ed6','ed7',
     1          'mt0','mt1','mt2','mt3','mt4','mt5','mt6','mt7'/
      data fname/'ed0','ed1','ed2','ed3','ed4','ed5','ed6','ed7',
     1           'mt0','mt1','mt2','mt3','mt4','mt5','mt6','mt7'/
c
      data ipos/16*-1/,iblksz/16*1/,isel,iselr,iselw,irep,ichek/5*0/
      data ilen/16*0/
      data nam/61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76/
c
      data maxb,iblkla/99999,1/
      data jrec/-1/
      data iwidth/72/,oswit/.false./
      data ierr/999/
      data nline,noline/0,0/
      data iwid/80/
c
      data num3,iblk3,orevis/4,1,2*.false./
      data apos/204*0.0d0/
c
c       *** nav = # integer/(real*8) word
c
      data iread,iwr/5,6/,oprint/.false./,nav/2/
c
      end
