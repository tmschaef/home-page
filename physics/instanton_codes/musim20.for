      program mutsim
c-------------------------------------------------------------------------c
c     finite temperature and density instanton liquid model. tsim runs    c
c     metropolis algorithm using the full gauge interaction and fermionic c
c     determinant.                                                        c 
c-------------------------------------------------------------------------c
c     version:             2.0                                            c
c     creation date:     09-21-96                                         c
c     last modification: 05-11-98                                         c
c-------------------------------------------------------------------------c
c     based on tsim, version: 4.4, created 11-03-93, modified 07-28-95    c
c-------------------------------------------------------------------------c
c     version 3.2: minimal fix for small T behavior of T_IA               c
c-------------------------------------------------------------------------c
c     version 3.3: eliminated anisotropic part of gluonic IA interaction. c
c-------------------------------------------------------------------------c
c     version 3.4: new alpha version. IMSL subroutines EVLHF, LFTHF re-   c
c     placed by LAPACK counterpart CHEEV. In order to convert back to     c
c     IMSL uncomment RNSET in main, change to RNUNF in rang, change to    c
c     EVLHF/LFTHF/LFDHF in logdet and spect.                              c
c-------------------------------------------------------------------------c
c     version 3.5: save configurations in outfile. Can be used as initial c
c     configuration for tsim, or ensemble for correlators.                c
c-------------------------------------------------------------------------c
c     version 3.6: change definition of relative orientation matrix to be c
c     in agreement with icor etc. Some minor improvements.                c
c-------------------------------------------------------------------------c
c     version 3.7: fix small T behavior of gluonic interaction.           c
c-------------------------------------------------------------------------c
c     version 3.8: replaced random number generator.                      c
c-------------------------------------------------------------------------c
c     version 3.9: core in ratio ansatz, new parameter fcut.              c
c-------------------------------------------------------------------------c
c     version 4.0: same input parameters as in tsimf, also corrected      c
c     normalisation of distribution of orientation angles.                c
c-------------------------------------------------------------------------c
c     version 4.1: calculate susceptibilities, print order par histories. c
c-------------------------------------------------------------------------c
c     version 4.2: changed dimension for simulations up to N=256.         c
c-------------------------------------------------------------------------c
c     version 4.3: additional input parameter for hot/cold start, study   c
c     dependence of susceptibilities on valence mass.                     c
c-------------------------------------------------------------------------c
c     version 4.4: include double precison LAPACK subroutines.            c
c-------------------------------------------------------------------------c
c     version 1.0: first version of musim.                                c
c-------------------------------------------------------------------------c
c     version 1.1: some bug fixes, more statistics for matrix elements.   c
c-------------------------------------------------------------------------c
c     version 1.3: more efficient calculation of determinant. Correct     c
c     phase averaging.                                                    c
c-------------------------------------------------------------------------c
c     version 1.4: remove all the finite temperature stuff.               c
c-------------------------------------------------------------------------c
c     version 1.5: asymmetric (mu!=0) overlap matrix elements.            c
c-------------------------------------------------------------------------c
c     version 1.6: improved parametrization of matrix elements.           c
c-------------------------------------------------------------------------c
c     version 1.7: simplified parametrization, collect T+M->T             c
c-------------------------------------------------------------------------c
c     version 1.8: calculate baryon density, chiral susceptibility. Small c
c     error in phase averaging fixed.                                     c
c-------------------------------------------------------------------------c
c     version 1.9: study collectivity of eigenvectors.                    c
c-------------------------------------------------------------------------c
c     version 2.0: minor additions, claculate P(lambda).                  c
c-------------------------------------------------------------------------c
c     input file: inmusim.dat                                             c
c     nc,nf   number of colors,flavors                                    c
c     nin     number of instantons                                        c
c     icold   cold start (icold=1), hot start (icold=0)                   c
c     rh0     average size                                                c
c     sg      metropolis step in group space                              c
c     dz(4)   metropolis step in position                                 c
c     drh     metropolis step in size                                     c
c     nit     number of ierations                                         c
c     alpha   Shuryaks overlap                                            c
c     fcut    core in ratio ansatz                                        c
c     rmu,rms fermion masses                                              c
c     ieq,kp1 number of ieration to equilibrate, print                    c
c     iread   input from file (iread=1), random start (iread=0)           c
c     steig   bin size for eigenvalue plot                                c
c     strho   same for size distribution                                  c
c     al(4)   box size, al(4) should be equal to tinv                     c
c     tinv    inverse temperature                                         c
c     inew    zero temperature check (inew=0), full interaction (inew=1)  c
c     ipis    pisarski factor included (ipis=1) or not (ipis=0)           c
c     isize   size renormalisation included (isize=1) or not              c
c     itwo    two loop (itwo=1) or one loop (itwo=0) beta function used   c
c     xmu     chemical potential                                          c
c     imu     options (1) quenched, (2) no phase, (3) unquenched          c
c-------------------------------------------------------------------------c
c     input file: infile.dat                                              c
c        determines start configuration if iread=1 is used. may use old   c
c        configuration as stored in outfile.                              c
c-------------------------------------------------------------------------c
c     output file: outmusim.dat                                           c
c        contains full statistics as well as simple histograms.           c
c     output file: outfile.dat                                            c
c        save configurations                                              c
c     output file: outplot.dat                                            c
c        simulation history, quark condensate, phase of determinant, etc. c
c        and data from histograms, size, eigenvalue distributions, etc.   c
c     output file: outmuspec.dat                                          c
c        complex eigenvalues for scatter plots of eigenvalue distribution.c
c-------------------------------------------------------------------------c
c     note: everything in units of lambda_QCD^(-1)                        c
c-------------------------------------------------------------------------c

      parameter(n=3, ni=256, ni2=128, nbin=150)
      complex u, cprod, ceig, u0, v, csom
      complex cqq, cqqt, cqq2, cqq2t, ci
      complex crho, crhop, crhot
      complex cqqp, cqqpt, cqq2p, cqq2pt
      complex eip, eipt, eip2t, eipa    
      complex zwr
      real imzwr, imzw 

      dimension rlp(ni2,ni2), rcos(ni2,ni2), rcii(ni2,ni2)
      dimension rnii(ni2,ni2), rnia(ni2,ni2), mol(ni2), p(ni2)
      dimension r4ii(ni2,ni2), r4ia(ni2,ni2)
      dimension rdii(ni2,ni2), rdia(ni2,ni2)
      dimension rmulp(ni2,ni2)
      dimension zr(ni,5),e1(ni,6),e2(ni,6), wr(ni2)  
      dimension zwr(ni2),rezwr(ni2),imzwr(ni2),zwre(ni2)
      dimension irho(nbin), ieig(nbin), ilap(nbin), icos(nbin)
      dimension inii(nbin), inia(nbin), icii(nbin), iact(nbin)
      dimension i4ia(nbin), i4ii(nbin), idii(nbin), idia(nbin)
      dimension mcos(nbin), mnia(nbin), m4ia(nbin), mlap(nbin)
      dimension mdia(nbin), mxmu(nbin)
      dimension ipar(nbin), ilow(nbin), ihig(nbin)
      dimension ireeig(nbin), iimeig(nbin), ieigre(nbin)
      dimension plam(nbin), plam2(nbin), nlam(nbin)
      dimension plamav(nbin), plamerr(nbin), iplam(nbin)
      dimension prelam(nbin), prelam2(nbin), nrelam(nbin)
      dimension prelamav(nbin), prelamerr(nbin), iprelam(nbin)
      dimension pimlam(nbin), pimlam2(nbin), nimlam(nbin)
      dimension pimlamav(nbin), pimlamerr(nbin), ipimlam(nbin)

      dimension itdeig(nbin,nbin), ictdeig(nbin,nbin)
      complex   ctdeig(nbin,nbin)

c-------------------------------------------------------------------------c
c     valence mass dependence                                             c
c-------------------------------------------------------------------------c

      dimension qqtv(nbin), qq2tv(nbin), qq4tv(nbin)
      dimension chiptv(nbin), chidtv(nbin), chip2tv(nbin)
      dimension chid2tv(nbin)

      common /param/al(4), alpha,rh0,sg,dz(4),drh, nc, nf, rmu, rms
      common /cden/ b,alc, p1, p2
      common /metr/ itot, ifail1, ifail2, ifail3, actt
      common /acti/ acold
      common /pi/   pi, eps
      common /sij/  sij(ni,ni)
      common /pibi/ pibi,tpar1,tpar2,ipis,isize,itwo
      common /cuts/ fcut,tc,delt
      common /tfac/ tfac(10)
      common /inew/ tinv, inew
      common /seed/ iseed
      common /musim/xmu,phase,imu

c------------------------------------------------------------------------c
c     read input file                                                    c
c------------------------------------------------------------------------c

      open(unit=1,file='inmusim.dat',status='unknown')
      open(unit=2,file='outmusim.dat',status='unknown')
      open(unit=7,file='outmuspec.dat',status='unknown')
      open(unit=8,file='outmueig.dat',status='unknown')
      open(unit=14,file='outfile.dat',status='unknown')
      open(unit=15,file='outplot.dat',status='unknown')

      read(1,*) nc, nf, nin, icold, rh0, sg
      read(1,*) (dz(k),k=1,4), drh, nit, alpha
      read(1,*) fcut,tc,delt
      read(1,*) rmu,rms,ieq,kp1
      read(1,*) iread, steig, strho
      read(1,*) (al(k),k=1,4),tinv,inew
      read(1,*) ipis,isize,itwo
      read(1,*) xmu,imu,bndim
      read(1,*) iout
      close (unit=1)

c------------------------------------------------------------------------c
c     echo input                                                         c
c------------------------------------------------------------------------c

      nlast  = (nit/kp1)*kp1
      nwrite = (nlast-ieq)/kp1+1
      call datime(id,it)
      call timed(time)

      write(15,*) ' conf     <qq>        S_av       ',
     &        ' Re(<qq>)    Im(<qq>)    phase'

      if (iout .eq. 1) then
      write(14,204) nc, nf, nin, al, rh0
      write(14,206) sg, (dz(k),k=1,4), drh, alpha
      write(14,207) rmu,rms,nwrite
      endif

      write(2,*)  ' musim version 2.0'
      write(2,998)
      write(2,997) id,it
      write(2,*)
      write(2,501) nc, nf, nin
      write(2,502) al
      write(2,503) rh0, sg, drh
      write(2,504) (dz(k),k=1,4)
      write(2,505) nit, alpha
      write(2,506) fcut,tc,delt
      write(2,507) rmu, rms
      write(2,508) ieq,kp1,iread
      write(2,509) tinv,inew,icold
      write(2,510) ipis,isize,itwo
      write(2,511) steig,strho
      write(2,500) xmu,imu,bndim

  501 format(1x,' N_c = ',i5,5x,' N_f = ',i5,5x,' N_in = ',i5)
  502 format(1x,' a_1 = ',f10.4,' a_2 = ',f10.4,
     2          ' a_3 = ',f10.4,' a_4 = ',f10.4)
  503 format(1x,' rho = ',f10.4,' du  = ',f10.4,
     2          ' drh = ',f10.4)
  504 format(1x,' d_1 = ',f10.4,' d_2 = ',f10.4,
     2          ' d_3 = ',f10.4,' d_4 = ',f10.4)
  505 format(1x,' nit = ',i5,5x,' alp = ',f10.4)
  506 format(1x,' fcut= ',f10.4,' tc  = ',f10.4,' delt = ',f10.4)
  507 format(1x,' m_u = ',f10.4,' m_s = ',f10.4)
  508 format(1x,' ieq = ',i5,5x,' ikp = ',i5,5x,' iread= ',i5)
  509 format(1x,' tin = ',f10.4,' inew= ',i5,5x,' icold= ',i5)
  510 format(1x,' ipis= ',i5,5x,' isiz= ',i5,5x,' itwo = ',i5)
  511 format(1x,' stei= ',f10.4,' strh= ',f10.4)
  500 format(1x,' xmu = ',f10.4,' imu = ',i2,8x,' bndim= ',f10.4)
  997 format(1x,' date : ',i6,' time : ',i4)

c------------------------------------------------------------------------c
c     parameters for size distribution                                   c
c------------------------------------------------------------------------c

      nd   = ni
      pi   = 4*atan(1.0)
      pibi = pi/tinv
      eps  = exp(-40*alog(2.0))
      vol  = al(1)*al(2)*al(3)*al(4)
      b    = 11.0/3.0*nc -2.0/3.0*nf
      bp   = 34.0/3.0*nc*nc - 13.0/3.0*nc*nf +nf/float(nc)
      p1   = 2*nc-bp/2/b
      p2   = bp/2/b
      cnc  = 1.34**nf *4.66*exp(-1.68*nc)/pi/pi/(nc-1)*(b/2)**(bp/2/b)
      alc  = alog(cnc)
      betas = tinv**2
      tpar1 = 2/3.*nc + nf/3.
      tpar2 = 1+nc/6.-nf/6.

c     iseed =-476
      iseed =-9234
c     iseed =-56789

      write(2,*)
      write(2,512) b,bp
      write(2,513) p1,p2
      write(2,514) cnc,iseed

  512 format(1x,' b   = ',f10.4,' bp  = ',f10.4)
  513 format(1x,' p1  = ',f10.4,' p2  = ',f10.4)
  514 format(1x,' Cnc = ',f10.4,' isd = ',i8)

c------------------------------------------------------------------------c
c     clear arrays etc.                                                  c
c------------------------------------------------------------------------c

c     call rnset(iseed)
      dum = rang( )

      call taumat

      call zero(nbin,ieig)
      call zero(nbin,ireeig)
      call zero(nbin,iimeig)
      call zero(nbin,ieigre)
      call zero(nbin,irho)
      call zero(nbin,ilap)
      call zero(nbin,icos)
      call zero(nbin,icii)
      call zero(nbin,inia)
      call zero(nbin,inii)
      call zero(nbin,i4ia)
      call zero(nbin,i4ii)
      call zero(nbin,idia)
      call zero(nbin,idii)

      call zero(nbin,ipar)
      call zero(nbin,ilow)
      call zero(nbin,ihig)

      call zero(nbin,nlam)
      call zero(nbin,nrelam)
      call zero(nbin,nimlam)
      call zerof(nbin,plam)
      call zerof(nbin,plam2)
      call zerof(nbin,prelam)
      call zerof(nbin,prelam2)
      call zerof(nbin,pimlam)
      call zerof(nbin,pimlam2)

      call zero2(nbin,nbin,itdeig)
      call zero2c(nbin,nbin,ctdeig)

      call zero(nbin,mlap)
      call zero(nbin,mcos)
      call zero(nbin,mnia)
      call zero(nbin,m4ia)
      call zero(nbin,mdia)
      call zero(nbin,mxmu)

      call zerof(nbin,qqtv)
      call zerof(nbin,qq2tv)
      call zerof(nbin,qq4tv)
      call zerof(nbin,chiptv)
      call zerof(nbin,chidtv)
      call zerof(nbin,chip2tv)
      call zerof(nbin,chid2tv)

      itott =  0
      actt  = 0.0
      rmdett= 0.0
      rvt   = 0.0
      uint  = 0.0
      qqt   = 0.0
      qq2t  = 0.0
      qq4t  = 0.0 

      cqqt  =(0.0,0.0)
      cqqpt =(0.0,0.0)
      cqqre2t = 0.0
      cqqim2t = 0.0
      cqqpre2t= 0.0
      cqqpim2t= 0.0

      cqq2t =(0.0,0.0)
      cqq2pt=(0.0,0.0)
      cqq2re2t = 0.0 
      cqq2im2t = 0.0
      cqq2pre2t= 0.0
      cqq2pim2t= 0.0

      rhoret  = 0.0
      rhore2t = 0.0
      rhore4t = 0.0
      rhoimt  = 0.0
      rhoim2t = 0.0
      rhoim4t = 0.0

      crhoret  = 0.0
      crhore2t = 0.0
      crhore4t = 0.0
      crhoimt  = 0.0
      crhoim2t = 0.0
      crhoim4t = 0.0
      
      chipt = 0.0
      chidt = 0.0
      chip2t= 0.0
      chid2t= 0.0  
      pht   = 0.0
      ph2t  = 0.0
      eipt  = 0.0
      eipre2t = 0.0
      eipim2t = 0.0

      ifailt1 = 0
      ifailt2 = 0
      ifailt3 = 0
      itel = 0
      nconf= 0
      ci = (0.0,1.0)

c------------------------------------------------------------------------c
c     initial configuration is made or read                              c
c------------------------------------------------------------------------c

      call setup(nc,nin,nd,zr,e1,e2,iread,icold)

      do 5 i = 1, nin
      do 5 j = 1, nin
         sij(j,i) = abs(sign(1.0,j-i+0.5)+sign(1.0,j-i-0.5))/2
    5 continue

c------------------------------------------------------------------------c
c     calculate action and fermionic determinant for first configuration c
c------------------------------------------------------------------------c

      call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,zwr,act,rv,uin,rlp,
     2          rmulp,rcos,rcii,rnii,rnia,r4ii,r4ia,rdii,rdia,mol,p)
      acold = act

      call dens(nc,nin,nd,zr,e1,e2,crho)

      qq = 0.0
      do 7 k=1,nin/2
         qq = qq + 2.0*rmu/(rmu**2+wr(k)**2)
   7  continue
      qq = qq/vol

      cqq =(0.0,0.0)
      do 8 k=1,nin/2
         cqq = cqq + 2.0*rmu/(rmu**2+zwr(k)**2)
   8  continue
      cqq = cqq/vol
      
      cqqp = cqq*exp(ci*phase)

      write(2,*)
      write(2,*)   ' Starting values '
      write(2,998)
      write(2,515) rmdet, rv, uin
      write(2,516) act, qq, cqq
      write(2,517) crho

  515 format(1x,' rmdet= ',f10.4,' rv  = ',f10.4,' uin = ',f10.4)
  516 format(1x,' act  = ',f10.4,' <qq>= ',f10.4,' <qq>= ',2f10.4)
  517 format(1x,' rhob = ',2f10.4)
      write(2,*)
      write(2,104)
      write(2,999)

c------------------------------------------------------------------------c
c     iterations start                                                   c
c------------------------------------------------------------------------c

      do 10 i = 1, nit

         itot = 0
         ifail1 = 0
         ifail2 = 0
         ifail3 = 0

         call iter(nc, nin,nd,zr,e1,e2)

         itott   = itott + itot
         ifailt1 = ifailt1 + ifail1
         ifailt2 = ifailt2 + ifail2
         ifailt3 = ifailt3 + ifail3

         if (i .eq. ieq) then
            rmdett= 0.0
            rvt   = 0.0
            uint  = 0.0
            actt  = 0.0

            qqt   = 0.0
            qq2t  = 0.0
            qq4t  = 0.0

            cqqt  =(0.0,0.0)
            cqqpt =(0.0,0.0)
            cqqre2t = 0.0 
            cqqim2t = 0.0
            cqqpre2t= 0.0
            cqqpim2t= 0.0

            cqq2t =(0.0,0.0)
            cqq2pt=(0.0,0.0)
            cqq2re2t = 0.0 
            cqq2im2t = 0.0
            cqq2pre2t= 0.0
            cqq2pim2t= 0.0

            chipt = 0.0
            chip2t= 0.0
            chidt = 0.0
            chid2t= 0.0   

            pht   = 0.0
            ph2t  = 0.0
            eipt  =(0.0,0.0)
            eipre2t = 0.0
            eipim2t = 0.0  

            rhoret  = 0.0
            rhore2t = 0.0
            rhore4t = 0.0
            rhoimt  = 0.0
            rhoim2t = 0.0
            rhoim4t = 0.0

            crhoret  = 0.0
            crhore2t = 0.0
            crhore4t = 0.0
            crhoimt  = 0.0
            crhoim2t = 0.0
            crhoim4t = 0.0

            call zerof(nbin,qqtv)
            call zerof(nbin,qq2tv)
            call zerof(nbin,qq4tv)
            call zerof(nbin,chiptv)
            call zerof(nbin,chidtv)
            call zerof(nbin,chip2tv)
            call zerof(nbin,chid2tv)
            itel  =  0
         end if
         itel = itel + 1

c------------------------------------------------------------------------c
c     calculate action and fermionic determinant                         c
c------------------------------------------------------------------------c

         call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,zwr,act,rv,uin,rlp,
     2           rmulp,rcos,rcii,rnii,rnia,r4ii,r4ia,rdii,rdia,mol,p)

         call dens(nc,nin,nd,zr,e1,e2,crho)

c------------------------------------------------------------------------c
c     average n(rho), gauge interaction, ferm. det., and total action    c
c------------------------------------------------------------------------c

         rvt   = rvt  + rv
         uint  = uint + uin
         rmdett= rmdett + rmdet
         actt  = actt + act

c------------------------------------------------------------------------c
c     quark condensate, susceptibilities                                 c
c------------------------------------------------------------------------c

         qq   = 0.0
         cqq  =(0.0,0.0)
         chip = 0.0
         chid = 0.0

         do 16 k=1,nin/2
            qq   = qq + 2.0*rmu/(rmu**2+wr(k)**2)
            chip = chip + 2.0/(rmu**2+wr(k)**2)
            chid = chid + 2.*(wr(k)**2-rmu**2)/(rmu**2+wr(k)**2)**2
  16     continue

         nreal = 0

         do 17 k=1,nin/2
            cqq  = cqq + 2.0*rmu/(rmu**2+zwr(k)**2)
            rezwr(k) = real(zwr(k))
            imzwr(k) = aimag(zwr(k))
            if(iout .eq. 1) write(7,*) rezwr(k),imzwr(k)
            if (abs(imzwr(k)) .lt. bndim) then
               nreal = nreal + 1
               zwre(nreal) = rezwr(k)
            endif
  17     continue                               

         qq    = qq/vol  
         cqq   = cqq/vol
         cqqp  = cqqp/vol
         chip  = chip/vol
         chid  = chid/vol

         eip  = exp(ci*phase)
         cqqp = cqq*eip
         crhop= crho*eip
         cqq2 = cqq**2
         cqq2p= cqq2*eip

         qqt   = qqt + qq
         qq2t  = qq2t + qq**2
         qq4t  = qq4t + qq**4

         cqqt  = cqqt + cqq
         cqqpt = cqqpt + cqqp
         cqq2t = cqq2t + cqq**2         
         cqq2pt= cqq2pt+ cqq2p

         rhoret  = rhoret + real(crho)
         rhore2t = rhore2t+ real(crho)**2
         rhore4t = rhore4t+ real(crho)**4
         rhoimt  = rhoimt + aimag(crho)
         rhoim2t = rhoim2t+ aimag(crho)**2
         rhoim4t = rhoim4t+ aimag(crho)**4

         crhoret  = crhoret + real(crhop)
         crhore2t = crhore2t+ real(crhop)**2
         crhore4t = crhore4t+ real(crhop)**4
         crhoimt  = crhoimt + aimag(crhop)
         crhoim2t = crhoim2t+ aimag(crhop)**2
         crhoim4t = crhoim4t+ aimag(crhop)**4
         
         cqqre2t = cqqre2t + real(cqq)**2
         cqqim2t = cqqim2t + aimag(cqq)**2
         cqqpre2t= cqqpre2t+ real(cqqp)**2
         cqqpim2t= cqqpim2t+ aimag(cqqp)**2

         cqq2re2t= cqq2re2t + real(cqq2)**2
         cqq2im2t= cqq2im2t + aimag(cqq2)**2
         cqq2pre2t= cqq2pre2t+ real(cqq2p)**2
         cqq2pim2t= cqq2pim2t+ aimag(cqq2p)**2

         chipt = chipt + chip
         chidt = chidt + chid
         chip2t= chip2t + chip**2
         chid2t= chid2t + chid**2 

         pht   = pht  + phase
         ph2t  = ph2t + phase**2
         eipt  = eipt + eip

         eipre2t = eipre2t + real(eip)**2
         eipim2t = eipim2t + aimag(eip)**2

c------------------------------------------------------------------------c
c     valence mass dependence                                            c
c------------------------------------------------------------------------c

         rmvmin = 0.0001
         rmvmax = 10.000
         nmv = 40
         xlgmax = log(rmvmax)
         xlgmin = log(rmvmin)
         dlog   = (xlgmax-xlgmin)/float(nmv-1)

         do 18 imv=1,nmv
            rmv = exp(xlgmin+(imv-1)*dlog)
            xqqv  = 0.0
            xchipv= 0.0
            xchidv = 0.0
            do 19 k=1,nin/2
               xqqv   = xqqv + 2.0*rmv/(rmv**2+wr(k)**2)
               xchipv = xchipv + 2.0/(rmv**2+wr(k)**2)
               xchidv = xchidv + 2.*(wr(k)**2-rmv**2)/
     1                                (rmv**2+wr(k)**2)**2
  19        continue
            xqqv    = xqqv/vol
            xchipv  = xchipv/vol
            xchidv  = xchidv/vol
            qqtv(imv)   = qqtv(imv) + xqqv
            qq2tv(imv)  = qq2tv(imv) + xqqv**2
            qq4tv(imv)  = qq4tv(imv) + xqqv**4
            chiptv(imv) = chiptv(imv) + xchipv
            chidtv(imv) = chidtv(imv) + xchidv
            chip2tv(imv)= chip2tv(imv) + xchipv**2
            chid2tv(imv)= chid2tv(imv) + xchidv**2
  18     continue

c------------------------------------------------------------------------c
c     histogram action, size, eigenvalues, relative orientation          c
c------------------------------------------------------------------------c

         write(15,993) i,qq,act/nin,cqq,phase

         if (i .ge. ieq) then
            stlap = 0.030
            stxmu = 0.030
            stcos = 0.025
            stdii = 0.040
            sact  = 0.100
            stc4  = 0.025
            stre  = steig
            stim  = 0.025
            nsre  = 40
            nsim  = 10
            stpar = 1
            xlow  = 0.1
            xhig  = 0.8
            call lens(act/nin,-6.0,sact,40,iact)
            do 50 k = 1, nin
               call lens(zr(k,5),0.0,strho,40,irho)
   50       continue
            do 55 k = 1, nreal
               call lens(zwre(k),0.0,steig,40,ieigre)
   55       continue
            do 60 k = 1, nin/2

               call lens(wr(k),0.0,steig,60,ieig)
               call lens(rezwr(k),0.0,steig,60,ireeig)
               call lens(imzwr(k),-30*steig,steig,60,iimeig)

               call tdlens(rezwr(k),imzwr(k),0.0,
     &           -nsim/2*stim,stre,stim,nsre,nsim,itdeig,nbin)     
               call ctdlens(rezwr(k),imzwr(k),eip,0.0,
     &           -nsim/2*stim,stre,stim,nsre,nsim,ctdeig,nbin)     

               call lens(rlp(mol(k),k),0.0,stlap,50,mlap)
               call lens(rmulp(mol(k),k),0.0,stxmu,50,mxmu)
               call lens(rcos(mol(k),k),-stcos,stcos,42,mcos)
               call lens(rnia(mol(k),k),-stcos,stcos,42,mnia)
               call lens(r4ia(mol(k),k),-stcos,stcos,42,m4ia)
               call lens(rdia(mol(k),k),0.0,stdii,50,mdia)

               azwr = abs(zwr(k))
               rezw = abs(real(zwr(k)))
               imzw = abs(aimag(zwr(k)))
               pk = p(k)
               call lens(pk,1.0,stpar,50,ipar)
               if(azwr .lt. xlow) then
                  call lens(pk,1.0,stpar,50,ilow)
               else if(azwr .gt. xhig) then
                  call lens(pk,1.0,stpar,50,ihig)
               endif

               call binadd(pk,azwr,plam,plam2,nlam,0.0,steig,60)
               if(imzw .lt. bndim) then
               call binadd(pk,rezw,prelam,prelam2,nrelam,0.0,steig,60)
               else if (rezw .lt. bndim) then
               call binadd(pk,imzw,pimlam,pimlam2,nimlam,0.0,steig,60)
               endif
                   
               do 60 l = 1, nin/2
                  call lens(rlp(k,l), 0.0, stlap,50, ilap)
                  call lens(rcos(k,l),-stcos,stcos,42,icos)
                  call lens(rnia(k,l),-stcos,stcos,42,inia)
                  call lens(r4ia(k,l),-stcos,stcos,42,i4ia)
                  if (l .ne. k) then
                     call lens(rnii(k,l),-stcos,stcos,42,inii)
                     call lens(rcii(k,l),-stcos,stcos,42,icii)
                     call lens(r4ii(k,l),-stcos,stcos,42,i4ii)
                  end if
   60       continue

         end if

c------------------------------------------------------------------------c
c     every kp1 sweeps do statistics, save configuration                 c
c------------------------------------------------------------------------c

         if ((i/kp1)*kp1 .eq. i) then

            write(6,*) 'iteration number ',i
            write(6,701) (ifail1+ifail2+ifail3)/(3.0*nin)
            write(6,702) act/nin
            write(6,703) qq
            write(6,704) cqq
            write(6,*)
            write(2,103) i,itott,ifailt1,ifailt2,ifailt3,
     1           act/nin,actt/(itel*nin),rmdet,rmdett/itel,
     2           uint/(nin*itel),rvt/(nin*itel),qqt/itel
            ifailt1 = 0
            ifailt2 = 0
            ifailt3 = 0
            itott = 0

            if(i .ge. ieq) then

            nconf = nconf + 1
            if(iout .eq. 1)then
            write(14,*) nconf
            do 25 j = 1, nin
               write(14,102) (zr(j,k),k=1,5)
               write(14,102) (e1(j,k), k=1,6)
               write(14,102) (e2(j,k), k=1,6)
   25       continue
            endif

            endif

         endif

 701     format(1x,'rejection rate   ',f10.4)
 702     format(1x,'average action   ',f10.4)
 703     format(1x,'quark condensate ',f10.4)
 704     format(1x,'quark condensate ',2f10.4)

c------------------------------------------------------------------------c
c     next iteration                                                     c
c------------------------------------------------------------------------c

   10 continue

c------------------------------------------------------------------------c
c     metropolis finished, calculate averages                            c
c------------------------------------------------------------------------c

      ntot = itel
      qqa  = qqt/ntot
      qqe  = disp(ntot,qqt,qq2t)
      qq2a = qq2t/ntot
      qq2e = disp(ntot,qq2t,qq4t)

      chipa= chipt/ntot
      chipe= disp(ntot,chipt,chip2t)
      chida= chidt/ntot
      chide= disp(ntot,chidt,chid2t)
      chica= (qq2a-qqa**2)*vol
      chice= (qq2e+2.0*qqa*qqe)*vol      

c------------------------------------------------------------------------c
c     complex condensate, phase not included                             c
c------------------------------------------------------------------------c

      cqqret = real(cqqt)
      cqqimt = aimag(cqqt)
      cqqrea = cqqret/ntot
      cqqima = cqqimt/ntot
      cqqree = disp(ntot,cqqret,cqqre2t)
      cqqime = disp(ntot,cqqimt,cqqim2t)

c------------------------------------------------------------------------c
c     complex condensate, phase included                                 c
c------------------------------------------------------------------------c

      cqqpret = real(cqqpt)
      cqqpimt = aimag(cqqpt)
      cqqprea = cqqpret/ntot
      cqqpima = cqqpimt/ntot
      cqqpree = disp(ntot,cqqpret,cqqpre2t)
      cqqpime = disp(ntot,cqqpimt,cqqpim2t)

c------------------------------------------------------------------------c
c     baryon density                                                     c
c------------------------------------------------------------------------c

      rhorea = rhoret/ntot
      rhoree = disp(ntot,rhoret,rhore2t)
      rhore2a= rhore2t/ntot
      rhore2e= disp(ntot,rhore2t,rhore4t)
    
      rhoima = rhoimt/ntot
      rhoime = disp(ntot,rhoimt,rhoim2t)
      rhoim2a= rhoim2t/ntot
      rhoim2e= disp(ntot,rhoim2t,rhoim4t)

c------------------------------------------------------------------------c
c     baryon density, phase included                                     c
c------------------------------------------------------------------------c

      crhorea = crhoret/ntot
      crhoree = disp(ntot,crhoret,crhore2t)
      crhore2a= crhore2t/ntot
      crhore2e= disp(ntot,crhore2t,crhore4t)

      crhoima = crhoimt/ntot
      crhoime = disp(ntot,crhoimt,crhoim2t)
      crhoim2a= crhoim2t/ntot
      crhoim2e= disp(ntot,crhoim2t,crhoim4t)

c------------------------------------------------------------------------c
c     phase and exp(i*phase)                                             c
c------------------------------------------------------------------------c

      pha  = pht/ntot
      pha  = mod(pha,pi)
      phe  = disp(ntot,pht,ph2t)

      eipa   = eipt/ntot 
      eipret = real(eipt)
      eipimt = aimag(eipt)
      eiprea = eipret/ntot
      eipima = eipimt/ntot
      eipree = disp(ntot,eipret,eipre2t)
      eipime = disp(ntot,eipimt,eipim2t)
      
c------------------------------------------------------------------------c
c     normalize complex condensate                                       c
c------------------------------------------------------------------------c

      call caverage(cqqprea,cqqpima,cqqpree,cqqpime,eiprea,eipima,
     &  eipree,eipime,cqqnrea,cqqnima,cqqnree,cqqnime)

c------------------------------------------------------------------------c
c     normalize complex baryon density                                   c
c------------------------------------------------------------------------c
 
      call caverage(crhorea,crhoima,crhoree,crhoime,eiprea,eipima,
     &  eipree,eipime,crhonrea,crhonima,crhonree,crhonime)
      call caverage(crhore2a,crhoim2a,crhore2e,crhoim2e,eiprea,eipima,
     &  eipree,eipime,crhonre2a,crhonim2a,crhonre2e,crhonim2e)

c------------------------------------------------------------------------c
c     baryon number susceptibility                                       c
c------------------------------------------------------------------------c

      rhosusa = crhonre2a - crhonrea**2
      rhosuse = crhonre2e - 2.0*crhonrea*crhonree
   
c------------------------------------------------------------------------c
c     complex scalar susceptibility, phase not included                  c
c------------------------------------------------------------------------c

      cqq2ret = real(cqq2t)
      cqq2imt = aimag(cqq2t)
      cqq2rea = cqq2ret/ntot
      cqq2ima = cqq2imt/ntot
      cqq2ree = disp(ntot,cqq2ret,cqq2re2t)
      cqq2ime = disp(ntot,cqq2imt,cqq2im2t)

c------------------------------------------------------------------------c
c     complex scalar susecptibility, phase included                      c
c------------------------------------------------------------------------c

      cqq2pret = real(cqq2pt)
      cqq2pimt = aimag(cqq2pt)
      cqq2prea = cqq2pret/ntot
      cqq2pima = cqq2pimt/ntot
      cqq2pree = disp(ntot,cqq2pret,cqq2pre2t)
      cqq2pime = disp(ntot,cqq2pimt,cqq2pim2t)
      
c------------------------------------------------------------------------c
c     normalize complex scalar susceptibility                            c
c------------------------------------------------------------------------c

      call caverage(cqq2prea,cqq2pima,cqq2pree,cqq2pime,eiprea,eipima,
     &  eipree,eipime,cqq2nrea,cqq2nima,cqq2nree,cqq2nime)
 
c------------------------------------------------------------------------c
c     subtract susceptibility, only real parts                           c
c------------------------------------------------------------------------c

      csusa = (cqq2rea-cqqrea**2)*vol
      csuse = (cqq2ree+2.0*cqqrea*cqqree)*vol      

      cnsusa = (cqq2nrea-cqqnrea**2)*vol
      cnsuse = (cqq2nree+2.0*cqqnrea*cqqnree)*vol      

c------------------------------------------------------------------------c
c     normalize complex eigenvalue distribution                          c
c------------------------------------------------------------------------c

      do 200 i=1,nsre
      do 200 j=1,nsim
         ictdeig(i,j) = real(ctdeig(i,j)/eipa)
         xre = (i-1)*stre
         xim = (j-1-nsim/2)*stim
         write(8,*) xre,xim,itdeig(i,j),ictdeig(i,j)
 200  continue

c------------------------------------------------------------------------c
c     average participation number                                       c
c------------------------------------------------------------------------c

      call binav(plam,plam2,nlam,60,plamav,plamerr)
      call binav(prelam,prelam2,nrelam,60,prelamav,prelamerr)
      call binav(pimlam,pimlam2,nimlam,60,pimlamav,pimlamerr)
      do 220 i=1,60
        iplam(i)   = int(plamav(i))
        iprelam(i) = int(prelamav(i))
        ipimlam(i) = int(pimlamav(i))
 220  continue

c------------------------------------------------------------------------c
c     output                                                             c
c------------------------------------------------------------------------c

      write(2,*)
      write(2,*)   ' susceptibilities'
      write(2,998)
      write(2,801) qqa,qqe
      write(2,802) qq2a,qq2e
      write(2,803) cqqrea,cqqima,cqqree,cqqime 
      write(2,804) cqqprea,cqqpima,cqqpree,cqqpime 
      write(2,805) cqqnrea,cqqnima,cqqnree,cqqnime
      write(2,806) chipa,chipe
      write(2,807) chida,chide
      write(2,808) chica,chice
      write(2,809) csusa,csuse
      write(2,810) cnsusa,cnsuse 
      write(2,811) pha,phe
      write(2,812) eiprea,eipima,eipree,eipime
      write(2,813) xmu**3/(3.0*pi**2)*vol
      write(2,814) rhorea,rhoima,rhoree,rhoime 
      write(2,815) crhorea,crhoima,crhoree,crhoime 
      write(2,816) crhonrea,crhonima,crhonree,crhonime
      write(2,817) rhosusa,rhosuse
      write(2,*)

 801  format(1x,' <QQ>          = ',f8.5,' +/- ',f8.5)
 802  format(1x,' <QQ^2>        = ',f8.5,' +/- ',f8.5)
 803  format(1x,' <qq>_||       =(',f8.5,',',f8.5,')+/-(',
     &                           f8.5,',',f8.5,')')
 804  format(1x,' <qq e^iph>_|| =(',f8.5,',',f8.5,')+/-(',
     &                           f8.5,',',f8.5,')')
 805  format(1x,' <qq e^iph>    =(',f8.5,',',f8.5,')+/-(',
     &                           f8.5,',',f8.5,')')
 806  format(1x,' Chi_pi (mu=0) = ',f8.5,' +/- ',f8.5)
 807  format(1x,' Chi_del(mu=0) = ',f8.5,' +/- ',f8.5)
 808  format(1x,' Chi_dis(mu=0) = ',f8.5,' +/- ',f8.5)   
 809  format(1x,' chi_dis       = ',f8.5,' +/- ',f8.5)
 810  format(1x,' chi_dis e^ip  = ',f8.5,' +/- ',f8.5)
 811  format(1x,' phase         = ',f8.5,' +/- ',f8.5)
 812  format(1x,' e^iph         =(',f8.5,',',f8.5,')+/-(',
     &                           f8.5,',',f8.5,')')
 813  format(1x,' <rho>_pert    = ',f8.5)
 814  format(1x,' <rho>_||      =(',f8.5,',',f8.5,')+/-(',
     &                           f8.5,',',f8.5,')')
 815  format(1x,' <rho e^iph>_||=(',f8.5,',',f8.5,')+/-(',
     &                           f8.5,',',f8.5,')')
 816  format(1x,' <rho e^iph>   =(',f8.5,',',f8.5,')+/-(',
     &                           f8.5,',',f8.5,')')
 817  format(1x,' rho_sus       = ',f8.5,' +/- ',f8.5)

c------------------------------------------------------------------------c
c     valence mass dependence                                            c
c------------------------------------------------------------------------c

      write(2,*)
      write(2,*) ' valence mass dependence'
      write(2,998)
      write(2,*)
      write(2,907)
      write(2,999)

      do 210 imv=1,nmv
         rmv  = exp(xlgmin+(imv-1)*dlog)
         qqa  = qqtv(imv)/ntot
         qqe  = disp(ntot,qqtv(imv),qq2tv(imv))
         qq2a = qq2tv(imv)/ntot
         qq2e = disp(ntot,qq2tv(imv),qq4tv(imv))
         chipa= chiptv(imv)/ntot
         chipe= disp(ntot,chiptv(imv),chip2tv(imv))
         chida= chidtv(imv)/ntot
         chide= disp(ntot,chidtv(imv),chid2tv(imv))
         chica= (qq2a-qqa**2)*vol
         chice= (qq2e+2.0*qq*qqe)*vol
         write(2,908) rmv,qqa,qqe,chida,chide,chica,chice
  210 continue

  907 format(1x,' m_val',7x,'<qq>',17x,'chi_del',14x,'chi_disc')
  908 format(1x,f8.5,1x,3(f11.5,1x,f9.5))

c------------------------------------------------------------------------c
c     participation numbers                                              c
c------------------------------------------------------------------------c

      write(2,*)
      write(2,*) ' participation numbers '
      write(2,998)
      write(2,*)
      write(2,909)
      write(2,999)
    
      do 215 i=1,60
         xlam = i*steig
         write(2,977) xlam,plamav(i),plamerr(i),prelamav(i),
     &   prelamerr(i),pimlamav(i),pimlamerr(i)
 215  continue

 909  format(2x,'lam',5x,'P(lam)',20x,'P|lam real',20x,'P|lam im')

c------------------------------------------------------------------------c
c     distributions, histograms, etc.                                    c
c------------------------------------------------------------------------c

      write(15,204)  nc, nf, nin, al, rh0
      write(15,206)  sg, (dz(k),k=1,4), drh, alpha
      write(15,205)  rmu,rms
  204 format(1x,3i5,5f10.4)
  206 format(1x,7f10.4)
  205 format(1x,2f10.4)
  207 format(1x,2f10.4,1x,i5)
  102 format(6f12.5)

      write(15,*) ' size distribution'
      do 601 k=1,40
         x = k*strho
         write(15,902) x,irho(k)
  601 continue
      write(15,*) ' eigenvalue distribution (mu=0)'
      do 602 k=1,60
         x = k*steig
         write(15,902) x,ieig(k)
  602 continue
      write(15,*) ' overlap distribution'
      do 603 k=1,50
         x = k*stlap
         write(15,902) x,ilap(k)
  603 continue
      write(15,*) ' max TIA overlap distribution'
      do 604 k=1,50
         x = k*stlap
         write(15,902) x,mlap(k)
  604 continue
      write(15,*) ' max MIA overlap distribution'
      do 1604 k=1,50
         x = k*stxmu
         write(15,902) x,mxmu(k)
 1604 continue
c     write(15,*) ' distribution of cos(th)^2 in IA'
c     do 605 k=1,42
c        x = -stcos+k*stcos
c        write(15,902) x,icos(k)
c 605 continue
c     write(15,*) ' distribution of cos(th)^2 in II'
c     do 606 k=1,42
c        x = -stcos+k*stcos
c        write(15,902) x,icii(k)
c 606 continue
c     write(15,*) ' distribution of spatial cos(th)^2 in IA'
c     do 607 k=1,42
c        x = -stcos+k*stcos
c        write(15,902) x,inia(k)
c 607 continue
c     write(15,*) ' distribution of spatial cos(th)^2 in II'
c     do 608 k=1,42
c        x = -stcos+k*stcos
c        write(15,902) x,inii(k)
c 608 continue
c     write(15,*) ' distribution of cos(al)^2 in IA'
c     do 609 k=1,42
c        x = -stcos+k*stcos
c        write(15,902) x,i4ia(k)
c 609 continue
c     write(15,*) ' distribution of cos(al)^2 in II'
c     do 610 k=1,42
c        x = -stcos+k*stcos
c        write(15,902) x,i4ii(k)
c 610 continue
      write(15,*) ' distribution of re(lam)'
      do 705 k=1,60
         x = k*steig
         write(15,902) x,ireeig(k)
  705 continue
      write(15,*) ' distribution of im(lam)'
      do 706 k=1,60
         x = k*steig
         write(15,902) x,iimeig(k)
  706 continue
      write(15,*) ' distribution of lam with im(lam) < ',bndim
      do 707 k=1,60
         x = k*steig
         write(15,902) x,ieigre(k)
  707 continue
      write(15,*) ' distribution of participation numbers'
      do 708 k=1,50
         x = k*stpar
         write(15,902) x,ipar(k)
  708 continue
      write(15,*) ' participation number for |lam| <',xlow
      do 709 k=1,50
         x = k*stpar
         write(15,902) x,ilow(k)
  709 continue
      write(15,*) ' participation number for |lam| >',xhig
      do 710 k=1,50
         x = k*stpar
         write(15,902) x,ihig(k)
  710 continue

  901 format(100(i6,/))
  902 format(1x,f10.5,i8)

c------------------------------------------------------------------------c
c     plot histograms for size, eigenvalues, overlaps, orientation, etc  c
c------------------------------------------------------------------------c

      write (2,113)
  113 format(1x,//,1x, '   Distribution of the sizes ')
      call lev(0.0,strho,40,2,irho)
      write (2,112)
  112 format(1x,//,1x, ' Distribution of the eigenvalues ')
      call lev(0.0,steig,60,2,ieig)
      write (2,114)
  114 format(1x,//,1x, ' Distribution of re(lam) ')
      call lev(0.0,steig,60,2,ireeig)
      write (2,115)
  115 format(1x,//,1x, ' Distribution of im(lam) ')
      call lev(0.0,steig,60,2,iimeig)
      write (2,116) bndim
  116 format(1x,//,1x, ' Distribution of lam with im(lam)<',f11.4)
      call lev(0.0,steig,60,2,ieigre)
c     write (2,117)
c 117 format(1x,' Distribution of the overlap matrix elements ')
c     call lev(0.0,stlap,50,2,ilap)
      write (2,117)
  117 format(1x,' Distribution of max TIA overlap matrix elements ')
      call lev(0.0,stlap,50,2,mlap)      
      write (2,1117)
 1117 format(1x,' Distribution of max MIA overlap matrix elements ')
      call lev(0.0,stxmu,50,2,mxmu)
      write (2,1118)
 1118 format(1x,' Distribution of separation in pairs ')
      call lev(0.0,stlap,50,2,mdia)
      write (2,1119)
 1119 format(1x,' distribution of participation numbers ')
      call lev(0.0,stpar,50,2,ipar)
      write (2,1120) xlow
 1120 format(1x,' participation number for |lam| <',f8.5)
      call lev(0.0,stpar,50,2,ilow)
      write (2,1121) xhig
 1121 format(1x,' participation number for |lam| >',f8.5)
      call lev(0.0,stpar,50,2,ihig)
      write (2,1122)
 1122 format(1x,' participation number P(lam) ')
      call lev(0.0,steig,60,2,iplam)
      write (2,1123)
 1123 format(1x,' participation number P(lam), lam real ')
      call lev(0.0,steig,60,2,iprelam)
      write (2,1124)
 1124 format(1x,' participation number P(lam), lam imag ')
      call lev(0.0,steig,60,2,ipimlam)

c     write (2,118)
c 118 format(1x, ' Distribution of cos^2(th) in IA ')
c     call lev(-stcos,stcos,42,2,icos)
c     write (2,119)
c 119 format(1x, ' Distribution of cos^2(th) in II ')
c     call lev(-stcos,stcos,42,2,icii)

c     write (2,131)
c 131 format(1x, ' Distribution of spatial cos^2(th) in IA ')
c     call lev(-stcos,stcos,42,2,inia)
c     write (2,132)
c 132 format(1x, ' Distribution of spatial cos^2(th) in II ')
c     call lev(-stcos,stcos,42,2,inii)
c     write (2,133)
c 133 format(1x, ' Distribution of cos^2(al) in IA ')
c     call lev(-stcos,stcos,42,2,i4ia)
c     write (2,134)
c 134 format(1x, ' Distribution of cos^2(al) in II ')
c     call lev(-stcos,stcos,42,2,i4ii)

c------------------------------------------------------------------------c
c     normalized distributions                                           c
c------------------------------------------------------------------------c

      call norm(nc,icos,42,0.0,stcos)
      call norm(nc,icii,42,0.0,stcos)
      call norm(2,inia,42,0.0,stcos)
      call norm(2,inii,42,0.0,stcos)
      call norm(nc,i4ia,42,0.0,stcos)
      call norm(nc,i4ii,42,0.0,stcos)

c     write(2,124)
c 124 format(1x, ' normalized Distribution of cos^2(th) ')
c     call lev(-stcos,stcos,42,2,icos)
c     write(2,125)
c 125 format(1x, ' normalized Distribution of sp cos^2(th) ')
c     call lev(-stcos,stcos,42,2,inia)
c     write(2,135)
c 135 format(1x, ' normalized Distribution of cos^2(al) ')
c     call lev(-stcos,stcos,42,2,i4ia)

      write(15,*) ' norm dist of cos(th)^2 in IA'
      do 611 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,icos(k)
  611 continue
      write(15,*) ' norm dist of cos(th)^2 in II'
      do 612 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,icii(k)
  612 continue
      write(15,*) ' norm dist of spatial cos(th)^2 in IA'
      do 613 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,inia(k)
  613 continue
      write(15,*) ' norm dist of spatial cos(th)^2 in II'
      do 614 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,inii(k)
  614 continue
      write(15,*) ' norm dist of cos(al)^2 in IA'
      do 615 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,i4ia(k)
  615 continue
      write(15,*) ' norm dist of cos(al)^2 in II'
      do 616 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,i4ii(k)
  616 continue

c------------------------------------------------------------------------c
c     distribution in molecules                                          c
c------------------------------------------------------------------------c

      write(2,120)
  120 format(1x, ' Distribution of cos^2(th) in pairs ')
      call lev(-stcos,stcos,42,2,mcos)
      write (2,121)
  121 format(1x, ' Distribution of sp. cos^2(th) in pairs ')
      call lev(-stcos,stcos,42,2,mnia)
      write (2,136)
  136 format(1x, ' Distribution of cos^2(al) in pairs ')
      call lev(-stcos,stcos,42,2,m4ia)

      call norm(nc,mcos,42,0.0,stcos)
      call norm(2,mnia,42,0.0,stcos)
      call norm(nc,m4ia,42,0.0,stcos)

      write(2,122)
  122 format(1x, ' normalized distribution of cos^2(th) in pairs ')
      call lev(-stcos,stcos,42,2,mcos)
      write (2,123)
  123 format(1x, ' normalized distribution sp. of cos^2(th) in pairs')
      call lev(-stcos,stcos,42,2,mnia)
      write (2,137)
  137 format(1x, ' normalized distribution of cos^2(al) in pairs')
      call lev(-stcos,stcos,42,2,m4ia)

      write(15,*) ' norm dist of cos(th)^2 in pairs'
      do 617 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,mcos(k)
  617 continue
      write(15,*) ' norm dist of spatial cos(th)^2 in pairs'
      do 618 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,mnia(k)
  618 continue
      write(15,*) ' norm dist of cos(al)^2 in pairs'
      do 619 k=1,42
         x = -stcos+k*stcos
         write(15,902) x,m4ia(k)
  619 continue

c     write (2,219)
c 219 format(1x, ' Distribution of the i-a norm  ')
c     call lev(-strii,strii,42,2,irii)
c     write (2,218)
c 218 format(1x, ' Distribution of the i-i norm  ')
c     call lev(-strii,strii,42,2,iria)
      write (2,318)
  318 format(1x, ' Distribution of the total action ')
      call lev(-6.0, sact,40,2,iact)

      call timed(time)
      write(2,*)
      write(2,996) time

  103 format(1x,5i5,6f8.3,f6.3)
  104 format(1x,'    i itot ifl1 ifl2 ifl3   S       S_av    m',
     2   '       m_av    u_av    w(rh) <qq>')
  992 format(1x,i5,f12.5,f12.5)  
  993 format(1x,i5,f12.5,f12.5,2f12.5,f12.5)
  977 format(1x,f8.5,6f11.5)
  996 format(1x,' execution time : ',f11.2)
  998 format(1x,20('-'))
  999 format(1x,78('-'))

      stop
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      subroutine caverage(orea,oima,oree,oime,eiprea,eipima,
     &  eipree,eipime,onrea,onima,onree,onime)
c------------------------------------------------------------------------c
c     calculate average and uncertainty for observabel normalized to     c
c     complex phase eip.                                                 c
c------------------------------------------------------------------------c
c     orea,oima   real, imaginary part of observable                     c
c     oree,oime   uncertainties in real, imaginary part                  c
c     eiprea,..   same for complex phase e^(i*phi)                       c
c     onrea,..    same for observable normalized to complex phase (outp) c   
c------------------------------------------------------------------------c

      den   = eiprea**2 + eipima**2 
      xnum1 = orea*eiprea + oima*eipima
      xnum2 = oima*eiprea - orea*eipima
      xnum3 = eiprea**2 - eipima**2
      xnum4 = orea*xnum3 + 2.0*oima*eiprea*eipima
      xnum5 =-oima*xnum3 + 2.0*orea*eiprea*eipima
      onrea = xnum1/den
      onima = xnum2/den
      sigre = ((eiprea*oree)**2+(eipima*oime)**2)/den**2
     & + ((xnum4*eipree)**2+(xnum5*eipime)**2)/den**4 
      sigim = ((eipima*oree)**2+(eiprea*oime)**2)/den**2
     & + ((xnum5*eipree)**2+(xnum4*eipime)**2)/den**4 
      onree = sqrt(sigre)
      onime = sqrt(sigim)

      return
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      function disp(n,xtot,x2tot)
c------------------------------------------------------------------------c
c     determine rms fluctuation, disp=del x                              c
c------------------------------------------------------------------------c
c     n     number of events                                             c
c     xtot  sum of x_i                                                   c
c     x2tot sum of x_i^2                                                 c
c------------------------------------------------------------------------c

      xav  = xtot/n
      del2 = x2tot/n**2 - xav**2/n
      del2 = max(del2,0.0)
      disp = sqrt(del2)

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine rsu(n, x, y)
c------------------------------------------------------------------------c
c     generate random su(n) matrix                                       c
c------------------------------------------------------------------------c
c     n     number of colors                                             c
c     x(n)  first row (real,imag) of su(n) matrix                        c
c     y(n)  first row (real,imag) of su(n) matrix                        c
c------------------------------------------------------------------------c

      dimension x(n), y(n), xs(6)

      call r6(x, n)
      if (n.eq.4) return

      xs(1) = x(2)
      xs(2) = - x(1)
      xs(3) = x(4)
      xs(4) = -x(3)
      xs(5) = x(6)
      xs(6) = -x(5)

      call r6(y,n)

      call rdot(n, x, y, xdy)
      call rdot(n, xs, y, xsdy)

      r2 = 0.0
      do 10 i = 1, n
        y(i) = y(i) - xdy*x(i) -xsdy*xs(i)
        r2 = r2 + y(i)*y(i)
   10 continue
      r = sqrt(r2)
      do 20 i = 1, n
        y(i) = y(i)/r
   20 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine r6(x,n)
c------------------------------------------------------------------------c
c     random n dimensional unit vector                                   c
c------------------------------------------------------------------------c

      dimension x(n)

 100  continue
      r2 = 0.0
      do 10 i = 1, n
        x(i) = 2*rang( ) - 1
        r2 = r2+x(i)*x(i)
  10  continue
      if (r2 .gt. 1) goto 100
      r = sqrt(r2)
      do 20 i = 1, n
        x(i) = x(i)/r
   20 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine rdot(n, x, y, xdy)
c------------------------------------------------------------------------c
c     n dimensional vector product xdy=x(n)*y(n)                         c
c------------------------------------------------------------------------c

      dimension x(n), y(n)

      xdy = 0.0
      do 10 i = 1, n
        xdy = xdy + x(i)*y(i)
   10 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine rg(n,e,sg,x)
c------------------------------------------------------------------------c
c     random rotation of n (even) dimensional unit vector e(n). Size of  c
c     rotation controlled by parameter sg.                               c
c------------------------------------------------------------------------c
c     input:   n     dimension of unit vector (=2*nc)                    c
c              e(n)  n dimensional unit vetor                            c
c              sg    parameter for size of mc step                       c
c     output:  x(n)  n dimensional unit vector                           c
c------------------------------------------------------------------------c

      dimension x(n), e(n)
      common /pi/ pi, eps

c------------------------------------------------------------------------c
c     random vector with gaussian distributed components                 c
c------------------------------------------------------------------------c

      do 10 i = 1, n/2
        x1 = rang( )
        x2 = rang( )
        x3 = rang( )
        x4 = rang( )
        ap = sqrt(-2*alog(x3))
        a  = sqrt(-2*alog(x1))
        x(2*i-1) = sg * a * cos(2*pi*x2)
        x(2*i)   = sg * ap* cos(2*pi*x4)
  10  continue

      call rdot(n,x,e,xde)
      r2 = 0.0
      do 30 i = 1, n
c       x(i) = e(i) + x(i) - xde*e(i)
        x(i) = e(i) + x(i)
        r2 = r2 + x(i)*x(i)
   30 continue

      r = sqrt(r2)
      do 20 i = 1, n
        x(i) = x(i)/r
   20 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine rgu(n, e1, e2, sg)
c------------------------------------------------------------------------c
c     random rotation of su(n/2) matrix. group element characterized by  c
c     first two rows e1(n),e2(n). sg controlls mc step size.             c
c------------------------------------------------------------------------c
c     n      dimension of unit vectors (=2*nc)                           c
c     e1(n)  first row of su(n/2) matrix (real,imag)                     c
c     e2(n)  second row                                                  c
c     sg     average size of mc step                                     c
c------------------------------------------------------------------------c

      complex zp
      dimension x(6), y(6), z(6), xs(6), e1(n), e2(n)
      common /pi/ pi, eps

      call rg(n, e1, sg, x)
      do 30 k = 1, n
        e1(k) = x(k)
  30  continue
      if (n.eq.4) return

      xs(1) = x(2)
      xs(2) =-x(1)
      xs(3) = x(4)
      xs(4) =-x(3)
      xs(5) = x(6)
      xs(6) =-x(5)

      call rg(n,e2,sg,y)

      call rdot(n, x, y, xdy)
      call rdot(n, xs, y, xsdy)

      r2 = 0.0

      do 10 i = 1, n
        y(i) = y(i) - xdy*x(i) - xsdy*xs(i)
        r2   =  r2  + y(i)*y(i)
   10 continue

      r = sqrt(r2)

      do 20 i = 1, n
        e2(i) = y(i)/r
        e1(i) = x(i)
   20 continue

      return
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      function rang()
c------------------------------------------------------------------------c
c     various random number generators                                   c
c------------------------------------------------------------------------c
      common /seed/ iseed
c     rang = rnunf()
      rang = ran2(iseed)
c     rang = ran(iseed)

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      function ran2(idum)
c----------------------------------------------------------------------c
c     numerical recipes random number generator ran2 (revised version) c
c     copr. 1986-92 numerical recipes software. Reinitialize with idum c
c     negative, then do not later idum between successive calls.       c
c----------------------------------------------------------------------c

      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     *ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      function ran2old(idum)
c---------------------------------------------------------------------------c
c     numerical recipes random number generator, set iseed to a negative    c
c     value to initialize sequence. on first call, sequence is initialized  c
c     automatically.                                                        c
c---------------------------------------------------------------------------c
c     note: cray version requires save statements not present in nr !       c
c---------------------------------------------------------------------------c
      parameter (m=714025,ia=1366,ic=150889,rm=1.4005112E-6)
      dimension ir(97)
      save ir,iy
      data iff /0/

      if (idum .lt. 0 .or. iff .eq. 0) then
        iff = 1
        idum= mod(ic-idum,m)
        do 11 j=1,97
          idum  = mod(ia*idum+ic,m)
          ir(j) = idum
11      continue
        idum = mod(ia*idum+ic,m)
        iy   = idum
      endif

      j = 1+(97*iy)/m
      if ( j .gt. 97 .or. J .lt. 1) pause
      iy   = ir(j)
      ran2old = iy*rm
      idum = mod(ia*idum+ic,m)
      ir(j)= idum

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine overl(n,nin,nd,nd2,zr,e1,e2,clp,clm,cl4,bfc,rhofc,
     &               uintt,il,iu,iop)
c------------------------------------------------------------------------c
c     calculate overlap matrix elements, gauge field action and inst.    c
c     density at finite temperature. can be used in two different modes: c
c     evaluate all matrix elements (iop=0), or update all matrix elem.   c
c     affected by changing the collective coordinates of a single inst.  c
c------------------------------------------------------------------------c
c     input:  n            number of colors                              c
c             nin          number of instantons                          c
c             nd           max number of instantons                      c
c             nd2          nd/2                                          c
c             zr(i,5)      position, size of instanton i                 c
c             e1(i,6)      first row of orientation matrix of inst i     c
c             e2(i,6)      second row                                    c
c     output: clp(i,j)     n/2*n/2 matrix of overlap matrix elements     c
c             clm(i,j)     n*n overlap matrix, chem potential included   c
c             cl4(i,j)     n*n matrix of gamma_4 matrix elements         c
c             bfc          instanton measure                             c
c             rhofc        Pisarski factor                               c
c             uintt        gauge field action                            c
c     input:  il           index of first I to be updated                c
c             iu           same for last I                               c
c             iop          mode of operation (see below)                 c
c------------------------------------------------------------------------c
c     iop=0   everything calculated                                      c
c     iop=1   rho update for I=il,..,iu                                  c
c     iop=2   position update for I=il,..,iu                             c
c     iop=3   orientation update for I=il,..,iu                          c
c------------------------------------------------------------------------c
c     common /param/  box size                                           c
c     common /cden/   parameters for two loop beta function              c
c     common /d/      store overlaps for update (iop=1,2,3)              c
c     common /tfac/   parametrization of overlap matrix elements         c
c     common /inew/   inverse temperature, inew=0 for T=0 check          c
c     common /pibi/   pi/T (no input), parameters for Pisarski           c
c------------------------------------------------------------------------c

      parameter(ni=128,n2=256)
      complex u1, u2, tr, cl, clp, clm, ctr, uc, utr, udotr
      complex u1c, u2c, u3c, u4c, ofac, ci, cusr, cusrt
      complex c4ij, ctij, c4ji, ctji, cl4
      dimension zr(nd,5), u1(n2,3,3),u2(n2,3,3),e1(nd,6), e2(nd,6)
      dimension tr(n2,2,2), clp(nd2,nd2), cl(ni,ni)
      dimension clm(nd,nd), cl4(nd,nd)
      dimension z1(5), z2(5), uc(n2,2,2)
      dimension ra(n2), rai(n2), si(n2)
      dimension ffr(n2)
      dimension ahl(4)
      dimension du1(n2), du2(n2)
      dimension drfc(n2), rfc(n2), zr2(n2), bet(n2), betp(n2)

      common /param/al(4),alpha,rh0,sg,dz(4),drh, nc, nf, rmu, rms
      common /cden/ b,alc,p1,p2
      common /pi/   pi, eps
      common /sij/  sij(n2,n2)
      common /tfac/ tfac(10)
      common /inew/ tinv, inew
      common /pibi/ pibi, tpar1, tpar2, ipis, isize, itwo
      common /cuts/ fcut, tc, delt
      common /d/    uint(n2,n2), dis(n2,n2,4), pcl1(n2,n2),
     &       pcl2(n2,n2), r2s(n2,n2), r2t(n2,n2), ofac(n2,n2),
     &       u1c(n2,n2), u2c(n2,n2), u3c(n2,n2), u4c(n2,n2),
     &       rcl1(n2,n2), rcl2(n2,n2), pcl1p(n2,n2), pcl2p(n2,n2),
     &       rcl1p(n2,n2), rcl2p(n2,n2)
      common /musim/xmu, phase, imu

c------------------------------------------------------------------------c
c     check : inew=0 --> zero temperature                                c
c------------------------------------------------------------------------c

      if (inew .eq. 0) then
         bett = 100.0
      else if (inew.eq. 1) then
         bett = tinv
      end if

      ci = cmplx(0.0, 1.0)
      nih = nin/2
      ahl(1) = al(1)/2
      ahl(2) = al(2)/2
      ahl(3) = al(3)/2
      ahl(4) = al(4)/2

c------------------------------------------------------------------------c
c     si = -/+ for instanton/antiinstanton                               c
c------------------------------------------------------------------------c

      do 5 j = 1, nin
        si(j) = sign(1.0,j-nih-0.5)
        zr2(j) = zr(j,5)*zr(j,5)
  5   continue

c------------------------------------------------------------------------c
c     reconstruct orientation U1, U2=U1^+                                c
c------------------------------------------------------------------------c

      call su3(n,nin,nd,e1, e2,u1)

      do 30 i1 = 1, n
      do 30 i2 = 1, n
         do 31 i = 1, nin
            u2(i,i1,i2) = conjg(u1(i,i2,i1))
 31      continue
 30   continue

      rhofc = 0.0
      uintt = 0.0
      do 216 i = 1, nin
         rfc(i) = 1.0
 216  continue

c------------------------------------------------------------------------c
c     for full calculation or position update recalculate distances      c
c------------------------------------------------------------------------c

      if(iop.eq.0 .or. iop.eq.2) then

      do 121 i = il, iu
         do 21 j = 1, nin
            r2s(j,i) = 0.0
            r2t(j,i) = 0.0
  21    continue
  121 continue

      do 212 i = il, iu
         do 20 ir = 1, 3
         do 22 j = 1, nin
            ds = zr(i,ir) - zr(j,ir)
            ads = abs(ds)
            asg = 1+sign(1.0,ads-ahl(ir))
            dis(j,i,ir) = ds - ahl(ir)*asg*sign(1.0,ds)
            r2s(j,i) = r2s(j,i) + dis(j,i,ir)**2
   22    continue
   20    continue

c------------------------------------------------------------------------c
c     symmetrize 4-direction ?                                           c
c------------------------------------------------------------------------c

         ir = 4
         do 922 j = 1, nin
            ds = zr(i,ir) - zr(j,ir)
            ads = abs(ds)
            asg = 1+sign(1.0,ads-ahl(ir))
            dis(j,i,ir) = ds - ahl(ir)*asg*sign(1.0,ds)
            r2t(j,i) = r2t(j,i) + dis(j,i,ir)**2
  922    continue
  212 continue

c------------------------------------------------------------------------c
c     for position update, also have to correct other ordering           c
c------------------------------------------------------------------------c

      do 215 i = il, iu
         do 214 ir = 1, 4
         do 214 j = 1, nin
            dis(i,j,ir) = -dis(j,i,ir)
  214    continue
         do 215 j = 1, nin
            r2s(i,j) = r2s(j,i)
            r2t(i,j) = r2t(j,i)
  215 continue

c------------------------------------------------------------------------c
c     end of position update                                             c
c------------------------------------------------------------------------c

      end if

c------------------------------------------------------------------------c
c     calculate orientation vectors (except for rho update)              c
c------------------------------------------------------------------------c

      if(iop.eq.0 .or. iop.eq.2 .or. iop.eq.3) then

      do 150 i = il, iu

c------------------------------------------------------------------------c
c     relative orientation matrix UC=U_I^+*U_J (DP conventions)          c
c------------------------------------------------------------------------c

         do 35 k = 1, 2
         do 35 l = 1, 2
         do 95 j = 1, nin
            uc(j,k,l) = cmplx(0.0, 0.0)
   95    continue

         do 35 m = 1, n
            do 97 j = 1, nin
               uc(j,k,l) = uc(j,k,l) + u2(i,k,m)*u1(j,m,l)
   97       continue
   35    continue

c------------------------------------------------------------------------c
c     orientation vector u1c_JI = tr(UC_IJ \tau_1^(+/-))/2i for (IA)/(AI)c
c------------------------------------------------------------------------c

         do 115 j = 1, nin
            sgp = 1 - (1+si(i))*(1-si(j))/2
            u1c(j,i) =  (uc(j,1,2)+uc(j,2,1))/2/ci
            u2c(j,i) = -(-uc(j,1,2)+uc(j,2,1))/2
            u3c(j,i) =  (uc(j,1,1)-uc(j,2,2))/2/ci
            u4c(j,i) = -(uc(j,1,1)+uc(j,2,2))/2*sgp
            ofac(j,i)=  (u1c(j,i)*dis(j,i,1) + u2c(j,i)*dis(j,i,2)
     1                 + u3c(j,i)*dis(j,i,3))
  115    continue

  150 continue

c------------------------------------------------------------------------c
c     for orientation update, correct reverse order                      c
c------------------------------------------------------------------------c

      do 117 i = il, iu
         do 116 j = 1, nin
            u1c(i,j) =-conjg(u1c(j,i))
            u2c(i,j) =-conjg(u2c(j,i))
            u3c(i,j) =-conjg(u3c(j,i))
            u4c(i,j) =-conjg(u4c(j,i))
            ofac(i,j)=+conjg(ofac(j,i))
  116    continue
  117 continue

c------------------------------------------------------------------------c
c     orientation update finished                                        c
c------------------------------------------------------------------------c

      end if

c------------------------------------------------------------------------c
c     big loop : calculate overlaps and gauge field action               c
c------------------------------------------------------------------------c

      do 12 i = il, iu
      do 80 j = 1, nin

c------------------------------------------------------------------------c
c     everything in units of sqrt(rho_I rho_J)                           c
c------------------------------------------------------------------------c

c        betas=bett**2/zr(i,5)/zr(j,5)
c        bets = betas
c        pibi = pi/sqrt(betas)
c        betn=sqrt(betas)

c------------------------------------------------------------------------c
c     factors gauge action and fermionic overlap                         c
c------------------------------------------------------------------------c
c     1-3: IA; 4-6: II; 7: f1; 8: kappa1; 9: f2; 10: kappa2              c
c------------------------------------------------------------------------c

c        tfac(1)=betas/(betas+5.21)
c        tfac(2)=betas/(betas+.75)
c        tfac(3)=1./(betas+1.73)
c        tfac(4)=betas/(betas+5.33)
c        tfac(5)=betas/(betas+1.17)
c        tfac(6)=1./(betas+2.08)
c        tfac(7)=bets/(.191*bets+1.)**2
c        tfac(8)=1.07/(bets/pi**2+0.533)
c        tfac(9)=1.+.069*bets/(0.011*bets+1.)**2
c        tfac(10)=1./(bets/pi**2+.69)

c------------------------------------------------------------------------c
c        minimal fix for low T behavior                                  c
c------------------------------------------------------------------------c

c        tfac(8) = 1.00/(bets/pi**2+0.533)
c        tfac(9) = 1.00

c------------------------------------------------------------------------c
c        variables for fermionic overlap, units of rho_I*rho_J           c
c------------------------------------------------------------------------c

         rhon=sqrt(zr(i,5)*zr(j,5))
         rat = zr(i,5)/zr(j,5)
         fact= 4/(rat+1/rat)**2
         r2  = r2s(j,i)+r2t(j,i)  
         rn  = r2/zr(i,5)/zr(j,5)
         r2si= sij(j,i)/(r2s(j,i)+eps)

         r1t = dis(j,i,4)/rhon
         r1s = sqrt(r2s(j,i))/rhon
         r1s = r1s + eps
c        sn = sin(pibi*r1t)
c        cs = cos(pibi*r1t)
c        sh = sinh(pibi*r1s)
c        ch = cosh(pibi*r1s)

c------------------------------------------------------------------------c
c        T_IJ = u_4*f1+...; pcl1=f1, pcl2=f2                             c
c------------------------------------------------------------------------c

c        fac1 = 1/((exp(-pibi*r1s/2) + pibi*r1s/2)/pibi**2
c    1                +2*(1-0.69*exp(-1.75*r1s/betn))/pi)
c        fac2 = 1+ 0.76*tfac(7)/(1+0.82*r1s*r1s)**2
c        fac2 = fac2*(1+(pibi*r1t)**2*0.178/(1+0.123*r1s*r1s))
c        fac3 = 1/((exp(-2.06*r1s/betn) + pibi*r1s/2)/pibi**2
c    1                +2*(1+0.42*exp(-.34*r1s/betn))/pi)
c        fac4 = tfac(9)

c        pcl1(j,i) =-pibi*sn*ch/(ch-cs+tfac(8))**2  *fac1 *fac2
c        pcl2(j,i) =-pibi*cs*sh/(ch-cs+tfac(10))**2 *fac3 *fac4

c        pcl1(j,i) = 1/rhon * pcl1(j,i)
c        pcl2(j,i) = 1/rhon * pcl2(j,i) * sqrt(r2si)

c------------------------------------------------------------------------c
c     check: zero temperature limit, mu=0                                c
c------------------------------------------------------------------------c

c        pcl1(j,i)= 4./(r1s**2+r1t**2+2)**2/rhon**2 * dis(j,i,4)
c        pcl2(j,i)= 4./(r1s**2+r1t**2+2)**2/rhon**2

c        pcl1p(j,i) = pcl1(j,i)
c        pcl2p(j,i) = pcl2(j,i)

c------------------------------------------------------------------------c
c     parametrizayion of T_IA (PAR5)                                     c
c------------------------------------------------------------------------c
  
         a1 =    1.33486101158204    
         a2 =    2.14268452272993     
         a3 =    2.87226236180785
         a4 =   0.756877352766242
         a5 =    2.37827701743211
         a6 =   0.893544815183667
         a7 =    1.06035128604853
         a8 =  -0.845176794766352     

         b1 =    2.36189735387699    
         b2 =    1.56110117795777     
         b3 =   7.230022501857282E-002
         b4 =    1.02709648541289     
         b5 =   0.818070505825908     
         b6 =   0.844197464583235    
         b7 =    3.27035365678753     
         b8 =    1.31391429361376    

c------------------------------------------------------------------------c
c     overlap matrix elements, mu!=0                                     c
c------------------------------------------------------------------------c

         xmun = xmu*rhon

         pcl1(j,i) = 2.d0/(2.0+rn)**2/rhon 
     &     * ( (2.d0*r1t + a1*xmun*r1t**2
     &     + a2*xmun*r1s**2 + a3*xmun**2*r1t )*cos(xmun*r1s) 
     &     - ( a4*xmun**2 + a5*r1s**2 - a6*r1t**2
     &     - a7*xmun*r1t*rn )*sin(xmun*r1s)/r1s 
     &     + a8*xmun**2*sin(xmun*r1s) )

         pcl2(j,i) = 2.0/(2.0+rn)**2/rhon**2/r1s 
     &     * ( (2.d0*r1s + b3*xmun**2*r1s**2
     &     - (b1*xmun*r1t**3 + b2*xmun*r1t*r1s**2)/r1s )*cos(xmun*r1s) 
     &     + (b1*r1t**3 + b4*3.d0*r1s**2*r1t  
     &     + b5*xmun*r1s**4 + b6*xmun*r1s**2*r1t**2 )
     &       *sin(xmun*r1s)/r1s**2 
     &     + (b7*xmun + b8*xmun**2*r1t)*sin(xmun*r1s) )

         pcl1p(j,i) = 2.d0/(2.0+rn)**2/rhon 
     &     * ( (2.d0*r1t + a1*(-xmun)*r1t**2
     &     + a2*(-xmun)*r1s**2 + a3*(-xmun)**2*r1t )*cos((-xmun)*r1s) 
     &     - ( a4*(-xmun)**2 + a5*r1s**2 - a6*r1t**2
     &     - a7*(-xmun)*r1t*rn )*sin((-xmun)*r1s)/r1s 
     &     + a8*(-xmun)**2*sin((-xmun)*r1s) )

         pcl2p(j,i) = 2.0/(2.0+rn)**2/rhon**2/r1s 
     &     * ( (2.d0*r1s + b3*(-xmun)**2*r1s**2
     &     - (b1*(-xmun)*r1t**3 + b2*(-xmun)*r1t*r1s**2)/r1s )
     &                                              *cos((-xmun)*r1s) 
     &     + (b1*r1t**3 + b4*3.d0*r1s**2*r1t  
     &     + b5*(-xmun)*r1s**4 + b6*(-xmun)*r1s**2*r1t**2 )
     &       *sin((-xmun)*r1s)/r1s**2 
     &     + (b7*(-xmun) + b8*(-xmun)**2*r1t)*sin((-xmun)*r1s) )

c-------------------------------------------------------------------------c
c     parametrization of M_IA (PAR5)                                      c
c-------------------------------------------------------------------------c

         a1 =    2.82924308859922     
         a2 =  -0.448009160464461    
         a3 =    1.56645720437148     
         a4 =   -2.07319129866595     
         a5 =   0.818929806872882    
         a6 =   -2.26578172373269     
         a7 =   0.296170345220125     
         a8 =   0.920621492331132    
       
         b1 =    3.76894634894014    
         b2 =   0.738194066223963    
         b3 =   -2.43307971010394    
         b4 =   0.738194066146792     
         b5 =  -0.799541520962931     
         b6 =   -5.65146975342066    
         b7 =    1.13260635638102     
         b8 =   0.647605956655825     

c------------------------------------------------------------------------c
c     gamma_4 matrix elements, mu=0                                      c
c------------------------------------------------------------------------c

c        rcl1(j,i) = 0.500/(1.0+0.1224*r1t**2+0.7175*r1s**4)
c        rcl2(j,i) = 0.169*r1t*r1s**2/(1.0+0.606*r1t**4+0.891*r1s**4) 
c        rcl2(j,i) = rcl2(j,i) * sqrt(r2si)
        
c        rcl1p(j,i) = rcl1(j,i)
c        rcl2p(j,i) = rcl2(j,i)

c------------------------------------------------------------------------c
c     gamma_4 matrix elements, mu!=0                                     c
c------------------------------------------------------------------------c
  
         rcl1(j,i) = 0.50/(1.d0 + 0.71*r1s**4 + 0.12*r1t**2)*
     &     ( (1.0 + a1*xmun*r1t + a2*xmun**2*r1s**2)*cos(xmun*r1s) 
     &     + (a3*xmun*r1s + a4*xmun*r1s**2 + a5*r1t)*sin(xmun*r1s) 
     &     + (a6 + a7*xmun*r1t + a8*xmun*r1t*r1s) 
     &            *r1t/r1s*sin(xmun*r1s)) 

         rcl2(j,i) = 0.170*r1s/(1.d0 + 0.89*r1s**4 + 0.61*r1t**4)*
     &     ( ( r1t*r1s**2 + b1*xmun*r1t + b2*xmun**2*r1s**2  
     &     + b3*xmun**2*r1t**2 + b4*xmun**2*r1s**2)*cos(xmun*r1s)
     &     + (b5 + b6*xmun*r1t + b7*xmun**2*r1s**2 + b8*xmun**2*r1s**3)
     &             *r1t/r1s*sin(xmun*r1s)) 
         rcl2(j,i) = rcl2(j,i) * sqrt(r2si)

         rcl1p(j,i)= 0.50/(1.d0 + 0.71*r1s**4 + 0.12*r1t**2)*
     &     ( (1.0 + a1*(-xmun)*r1t + a2*(-xmun)**2*r1s**2)
     &                                             *cos((-xmun)*r1s) 
     &     + (a3*(-xmun)*r1s + a4*(-xmun)*r1s**2 + a5*r1t)
     &                                             *sin((-xmun)*r1s) 
     &     + (a6 + a7*(-xmun)*r1t + a8*(-xmun)*r1t*r1s) 
     &            *r1t/r1s*sin((-xmun)*r1s)) 

         rcl2p(j,i)= 0.170*r1s/(1.d0 + 0.89*r1s**4 + 0.61*r1t**4)*
     &     ( ( r1t*r1s**2 + b1*(-xmun)*r1t + b2*(-xmun)**2*r1s**2  
     &     + b3*(-xmun)**2*r1t**2 + b4*(-xmun)**2*r1s**2)
     &                                              *cos((-xmun)*r1s)
     &     + (b5 + b6*(-xmun)*r1t + b7*(-xmun)**2*r1s**2 
     &     + b8*(-xmun)**2*r1s**3)
     &             *r1t/r1s*sin((-xmun)*r1s)) 
         rcl2p(j,i)= rcl2p(j,i) * sqrt(r2si)
           
c------------------------------------------------------------------------c
c     orientation invariants                                             c
c------------------------------------------------------------------------c

         r1  = sqrt(rn)
         r1i = sij(j,i)/sqrt(rn+eps)
         rr2i= sij(j,i)/(r2s(j,i)+r2t(j,i)+eps)
         tau = dis(j,i,4)
         sgp = 1 - (1+si(i))*(1-si(j))/2
         sge = si(i)*si(j)

         us1   = cabs(u1c(j,i))**2
         us23  = cabs(u2c(j,i))**2 + cabs(u3c(j,i))**2
         us4   = cabs(u4c(j,i))**2
         us    = us1 + us23 + us4
         us123 = us1 + us23
         us123s= us123*us123

         udotr = ofac(j,i)
         cusr  = (udotr + tau*u4c(j,i))*sqrt(rr2i)
         cusrs = cusr*conjg(cusr)
         cusrt = udotr*sqrt(rr2i)
         cusrst= cusrt*conjg(cusrt)
         utrans= us - cusrst - us4

c------------------------------------------------------------------------c
c     collect terms for II and IA interaction                            c
c------------------------------------------------------------------------c

c        uii1 = 1/(1+0.43*rn)**3*tfac(4)
c        uii2 = alog(rn+eps)/(1+1.17*rn)**4 * tfac(5)
c        uii3 = alog(1+betn*r1i)*tfac(6)
c        uia1 = 4.03/(rn+2.10)**2 *tfac(1)
c        uia2 = -(1.66/(1+1.68*rn)**3 +
c    1             0.72*alog(rn+eps)/(1+0.42*rn)**4)*tfac(2)
c        uia3 = (16.16/(rn+2.03)**2 - 2.73/(1+0.33*rn)**3)*
c    2             betas/(betas+0.24+11.5*rn/(1+1.14*rn))
c        uia4 = 0.36*alog(1+betn*r1i)/(1+0.013*rn)**4*tfac(3)

c------------------------------------------------------------------------c
c     fix small T behavior                                               c
c------------------------------------------------------------------------c

c        uia1 = 4.00/(rn+2.00)**2 *tfac(1)
c        uia3 = (16.00/(rn+2.00)**2 - 2.73/(1+0.33*rn)**3)*
c    2             betas/(betas+0.24+11.5*rn/(1+1.14*rn))

c------------------------------------------------------------------------c
c     T=0 limit                                                          c
c------------------------------------------------------------------------c

         uii10 = 1/(1+0.43*rn)**3
         uii20 = alog(rn+eps)/(1+1.17*rn)**4 
         uia10 = 4.00/(rn+2.00)**2 
         uia20 = -(1.66/(1+1.68*rn)**3 +
     1             0.72*alog(rn+eps)/(1+0.42*rn)**4)
         uia30 = 16.00/(rn+2.00)**2 - 2.73/(1+0.33*rn)**3

c------------------------------------------------------------------------c
c     anisotropic part of IA interaction (not used!)                     c
c------------------------------------------------------------------------c

c        uiap4 = -0.66*betn/(betas+0.58)*sin(pibi*r1t)

c------------------------------------------------------------------------c
c     du2 : IA interaction; du1 : II interaction                         c
c------------------------------------------------------------------------c
c     sij = 0 for I=J; (1-sge)=0 for (II), (AA); (1+sge)=0 for (IA),(AI) c
c------------------------------------------------------------------------c

c        fact   = 1.0
c        du2(j) = (uia1*us + uia2*us - uia3*cusrs + uia4*utrans
c    1                + uiap4*us4*0.0) * fact * sij(j,i)*(1-sge)/2
c        du1(j) = ((0.63*us123 + 0.071*us123s)*uii1
c    1              - (0.05*us123 + 0.47*us123s)*uii2
c    2              + (0.07*us123 + 0.05*us123s)*uii3)
c    3               * fact * sij(j,i)*(1+sge)/2

c------------------------------------------------------------------------c
c     T=0 limit                                                          c
c------------------------------------------------------------------------c

         du2(j) = (uia10*us + uia20*us - uia30*cusrs) 
     1               * sij(j,i)*(1-sge)/2
         du1(j) = ((0.63*us123 + 0.071*us123s)*uii10
     1              - (0.05*us123 + 0.47*us123s)*uii20)
     3               * sij(j,i)*(1+sge)/2

c------------------------------------------------------------------------c
c     hard core, controlled by parameter fcut                            c
c------------------------------------------------------------------------c

c        acf   = rn+(zr2(i)+zr2(j))/rhon**2
c        disc  = sqrt(abs(acf**2-4))
c        rlam  = (acf+disc)/2.0
c        core  = fcut*us/rlam**4
c        du1(j)= du1(j) + core*sij(j,i)
c        du2(j)= du2(j) + core*sij(j,i)

c------------------------------------------------------------------------c
c     include periodic image for anisotropic part of IA interaction      c
c------------------------------------------------------------------------c

c        taup = tau - sign(bett,tau)
c        r1tp = taup/sqrt(zr(i,5)*zr(j,5))
c        rnp  = sqrt(r1tp**2+r1s**2)
c        r1ip = sij(j,i)/sqrt(rnp+eps)
c        rr2ip= sij(j,i)/(r2s(j,i)+taup**2+eps)

c        cusr  = (udotr + taup*u4c(j,i))*sqrt(rr2ip)
c        cusrs = cusr*conjg(cusr)
c        cusrt = udotr*sqrt(rr2ip)
c        cusrst= cusrt*conjg(cusrt)
c        utrans= us - cusrst - us4

c------------------------------------------------------------------------c
c     collect terms for IA interaction                                   c
c------------------------------------------------------------------------c

c        uia1 = 4.03/(rnp+2.10)**2 *tfac(1)
c        uia2 = -(1.66/(1+1.68*rnp)**3 +
c    1             0.72*alog(rnp+eps)/(1+0.42*rnp)**4)*tfac(2)
c        uia3 = (16.16/(rnp+2.03)**2 - 2.73/(1+0.33*rnp)**3)*
c    2             betas/(betas+0.24+11.5*rnp/(1+1.14*rnp))
c        uia4 = 0.36*alog(1+betn*r1ip)/(1+0.013*rnp)**4*tfac(3)

c------------------------------------------------------------------------c
c     add interaction with mirror image                                  c
c------------------------------------------------------------------------c

c        du2(j) = du2(j) + (uia1*us + uia2*us - uia3*cusrs +
c    1               uia4*utrans ) * fact * sij(j,i)*(1-sge)/2

c------------------------------------------------------------------------c
c     end of loop over instantons j                                      c
c------------------------------------------------------------------------c

 80   continue

c------------------------------------------------------------------------c
c     collect gauge interaction; check: switch off II interaction        c
c------------------------------------------------------------------------c

      do 105 j = 1, nin
         uint(j,i) = du1(j) + du2(j)
c        uint(j,i) = du2(j)
  105 continue

c------------------------------------------------------------------------c
c     end of big loop over instantons I                                  c
c------------------------------------------------------------------------c

   12 continue

c------------------------------------------------------------------------c
c     reverse order (needed in update mode only)                         c
c------------------------------------------------------------------------c

      do 18 i = il, iu
      do 18 j = 1, nin
        uint(i,j) = uint(j,i)
        pcl1(i,j) =-pcl1(j,i)
        pcl2(i,j) = pcl2(j,i)
        pcl1p(i,j)=-pcl1p(j,i)
        pcl2p(i,j)= pcl2p(j,i)
        rcl1(i,j) = rcl1(j,i)
        rcl2(i,j) =-rcl2(j,i)
        rcl1p(i,j)= rcl1p(j,i)
        rcl2p(i,j)=-rcl2p(j,i)
   18 continue

c------------------------------------------------------------------------c
c     finally: nih*nih matrix of IA-overlap matrix elements              c
c------------------------------------------------------------------------c

      do 311 i=1, nin
      do 311 j=1, nin
         clm(i,j) = (0.0,0.0)
 311  continue

      do 312 i = 1, nih
      do 100 j = nih+1, nin
          ctij = (u4c(j,i)*pcl1(j,i) +ofac(j,i)*pcl2(j,i) )*ci
          ctji = (u4c(j,i)*pcl1p(j,i)+ofac(j,i)*pcl2p(j,i))*ci
          ctji = conjg(ctji)
          c4ij = (u4c(j,i)*rcl1(j,i) +ofac(j,i)*rcl2(j,i) )*ci
          c4ji = (u4c(j,i)*rcl1p(j,i)+ofac(j,i)*rcl2p(j,i))*ci
          c4ji = conjg(c4ji)
          clp(i,j-nih) = ctij 
          clm(i,j) = ctij 
          clm(j,i) = ctji 
          cl4(i,j) = c4ij
          cl4(j,i) = c4ji
c         clm(i,j) = ctij  + xmu*c4ij
c         clm(j,i) = ctji  - xmu*c4ji 
 100  continue
 312  continue

c------------------------------------------------------------------------c
c     renormalized size rfc(i) to be used for running coupling           c
c------------------------------------------------------------------------c

      do 410 i = 1, nin
      do 420 j = 1, nin
        ra(j) = (r2s(j,i)+r2t(j,i))/zr2(j)
        ffr(j) = sij(j,i)/(ra(j)+eps)
        rfc(i) = rfc(i) + ffr(j)
  420 continue
  410 continue

c------------------------------------------------------------------------c
c     log of density distribution, T!=0 switched off                     c
c------------------------------------------------------------------------c

      do 15 i = 1, nin
         ro = zr2(i)
         rop= zr2(i)/rfc(i)
         bet(i) = -b/2*alog(ro)
         betp(i)= -b/2*alog(rop)
c        z = pibi*zr(i,5) / rhon
c        t1     = 1.0/bett
c        fpis   = (tanh((t1-tc)/delt)+1.0)/2.0
c        tpisar = - tpar1*z**2
c    &          -tpar2*(-alog(1+z*z/3)
c    &          + 0.155*exp(-8*alog(1+0.159/z**1.5) ) )
c        drfc(i)= (p1+itwo*p1*p2/bet(i))*alog(bet(i))-bet(i)
c    &            -5*alog(zr(i,5))+ alc + tpisar*fpis*ipis
         xmupis = -nf*(xmu*zr(i,5))**2
         drfc(i)= (p1+itwo*p1*p2/bet(i))*alog(bet(i))-bet(i)
     &            -5*alog(zr(i,5))+ alc + xmupis*ipis
   15 continue

c------------------------------------------------------------------------c
c     collect density distribution                                       c
c------------------------------------------------------------------------c

      do 16 i = 1, nin
        rhofc = rhofc + drfc(i)
   16 continue

c------------------------------------------------------------------------c
c     collect gauge interaction, S_0=bet(i)                              c
c------------------------------------------------------------------------c

      do 17 i = 1, nin
      do 17 j = 1, nin
        uintt = uintt + bet(i)*uint(j,i)
   17 continue
      uintt = uintt/2.0

c------------------------------------------------------------------------c
c     size renormalization                                               c
c------------------------------------------------------------------------c

c     if(isize .eq. 1)then
c        do 19 i=1,nin
c           uintt = uintt + betp(i)-bet(i)
c  19    continue
c     endif

c------------------------------------------------------------------------c
c     log of bosonic weight function n(rho)*exp(-S)                      c
c------------------------------------------------------------------------c

      bfc = rhofc-uintt

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      function fs(r2)
      common /param/ al(4), alpha,rh0,sg,dz(4),drh, nc, nf, rmu, rms

      fs = 2/(r2+alpha)**2

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      function ff(x)

      data c /0.2/

      arg = c*x
      if (arg.gt.20) then
        ff = 0.0
      else
        ff = exp(-arg)
      end if

      return
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      subroutine setup(n,nin,nd,zr,e1,e2,iread,icold)
c------------------------------------------------------------------------c
c     initialize instanton distribution                                  c
c------------------------------------------------------------------------c
c     n,nin,nd   nc, number of instantons, array dimension               c
c     zr,e1,e2   instanton configuration                                 c
c     icold      icold=1: use cold start (polarized molecules)           c
c------------------------------------------------------------------------c
c     iread = 1  input from file                                         c
c     iread = 0  random configuration                                    c
c------------------------------------------------------------------------c

      parameter(ndd=256)

      complex u0, tau
      complex u(ndd,3,3), d(3,3)

      dimension zr(nd,5), er1(6), er2(6), fal(4)
      dimension x(6), y(6), z(6), e1(nd,6), e2(nd,6)
      dimension fdz(4), uv(4)

      common /param/ al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /taum/ tau(2,4,2,2)

      if (iread.eq.1) then

c------------------------------------------------------------------------c
c     input from file infile                                             c
c------------------------------------------------------------------------c

      open(unit=4,file='infile.dat',status='unknown')
      read(4,204) inc, inf, inin,(fal(k),k=1,4),frh0
      read(4,206) fsg, (fdz(k),k=1,4), fdrh, falpha
      read(4,205) frmu,frms
      write(2,204) inc, inf, inin, fal, frh0
      write(2,206) fsg, (fdz(k),k=1,4), fdrh, falpha
      write(2,205) frmu,frms
      read(4,*) nconfig
      do 25 i = 1, nin
        read(4,101) (zr(i,k),k=1,5)
        read(4,101) (e1(i,k), k=1,6)
        read(4,101) (e2(i,k), k=1,6)
   25 continue
      write (6,*) ' configuration has been read'


  101 format(1x,6f12.5)
  204 format (1x,3i5,5f10.4)
  205 format (1x,2f10.4)
  206 format (1x,7f10.4)

      else

c------------------------------------------------------------------------c
c     random configuration                                               c
c------------------------------------------------------------------------c

      do 100 ip = 1, nin
         do 10 i = 1, 4
           zr(ip,i) = al(i)*rang( )
c          zr(ip,i) = ip*al(i)/(nin+1)
  10     continue

c------------------------------------------------------------------------c
c     but fixed size                                                     c
c------------------------------------------------------------------------c

c-jac    zr(ip,5) = rh0+drh*(rang( )-0.5)
         zr(ip,5) = rh0

c------------------------------------------------------------------------c
c     random orientation                                                 c
c------------------------------------------------------------------------c

         do 20 i = 1, 6
            e1(ip,i) = 0.0
            e2(ip,i) = 0.0
            er1(i) = 0.0
            er2(i) = 0.0
  20     continue
         e1(ip,1) = 1.0
         e2(ip,3) = 1.0
         nu2 = 2*n
         call rsu(nu2,er1,er2)
         do 70 i = 1, nu2
            e1(ip,i) = er1(i)
            e2(ip,i) = er2(i)
  70     continue

  100 continue

      if (icold .eq. 1) then

c--------------------------------------------------------------------------c
c     replace instantons by molecules; reconstruct orientation             c
c--------------------------------------------------------------------------c

      call su3(3,nin,nd,e1,e2,u)

c--------------------------------------------------------------------------c
c     instanton coordinates unchanged, antiinstanton shifted + rotated     c
c--------------------------------------------------------------------------c

      drho = 2.00*rh0
      nin2 = nin/2.0
      do 200 ii=1,nin2
         ia = ii + nin2

c--------------------------------------------------------------------------c
c     polarized molecules                                                  c
c--------------------------------------------------------------------------c

         zr(ia,1) = zr(ii,1)
         zr(ia,2) = zr(ii,2)
         zr(ia,3) = zr(ii,3)
         zr(ia,4) = zr(ii,4)+drho*sign(1.0,rang()-0.5)
         uv(1) = 0.0
         uv(2) = 0.0
         uv(3) = 0.0
         uv(4) = (zr(ia,4)-zr(ii,4))/drho
         if(zr(ia,4) .gt. al(4)) zr(ia,4)=zr(ia,4)-al(4)
         if(zr(ia,4) .lt.  0.0 ) zr(ia,4)=zr(ia,4)+al(4)

c--------------------------------------------------------------------------c
c     orientation is taken to be U(ia) = U(i) * uv_\mu (\tau^(-))_\mu      c
c--------------------------------------------------------------------------c

         do 220 k1=1,3
         do 220 k2=1,3
            d(k1,k2) = (0.0,0.0)
 220     continue

         do 230 k1=1,2
         do 230 k2=1,2
            do 230 mu=1,4
               d(k1,k2) = d(k1,k2) + uv(mu)*tau(2,mu,k1,k2)
 230     continue
         d(3,3) = (1.0,0.0)

c--------------------------------------------------------------------------c
c        clear old, assign new orientation                                 c
c--------------------------------------------------------------------------c

         do 240 k1=1,3
         do 240 k2=1,3
            u(ia,k1,k2) = (0.0,0.0)
 240     continue

         do 250 k1=1,3
         do 250 k2=1,3
         do 250 m=1,3
            u(ia,k1,k2) = u(ia,k1,k2)+u(ii,k1,m)*d(m,k2)
 250     continue

c--------------------------------------------------------------------------c
c     change e1,e2 accordingly                                             c
c--------------------------------------------------------------------------c

         do 260 k=1,3
            e1(ia,2*k-1) = real(u(ia,k,1))
            e1(ia,2*k)   =aimag(u(ia,k,1))
            e2(ia,2*k-1) = real(u(ia,k,2))
            e2(ia,2*k)   =aimag(u(ia,k,2))
 260     continue

 200  continue

c------------------------------------------------------------------------c
c     end of molecule setup                                              c
c------------------------------------------------------------------------c

      end if

c------------------------------------------------------------------------c
c     end of iread                                                       c
c------------------------------------------------------------------------c

      end if

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine su3(n,nin,nd,x,y,u)
c------------------------------------------------------------------------c
c     reconstruct nin complex su(3) matrices u(i,3,3) from first two     c
c     rows  x(i,6),  y(i,6) (real, imaginary parts)                      c
c------------------------------------------------------------------------c

      parameter(n2=256)
      complex u
      dimension x(nd,6), y(nd,6), z(n2,6), u(nd,3,3)

      if (n.eq.2) goto 100

      do 10 i = 1, nin
         xl1 = x(i,1)
         xl2 = x(i,2)
         xl3 = x(i,3)
         xl4 = x(i,4)
         xl5 = x(i,5)
         xl6 = x(i,6)
         yl1 = y(i,1)
         yl2 = y(i,2)
         yl3 = y(i,3)
         yl4 = y(i,4)
         yl5 = y(i,5)
         yl6 = y(i,6)
         z(i,1) = xl3*yl5-xl4*yl6-xl5*yl3+xl6*yl4
         z(i,2) = -xl3*yl6-xl4*yl5+xl5*yl4+xl6*yl3
         z(i,3) = xl5*yl1-xl6*yl2-xl1*yl5+xl2*yl6
         z(i,4) = -xl5*yl2-xl6*yl1+xl1*yl6+xl2*yl5
         z(i,5) = xl1*yl3-xl2*yl4-xl3*yl1+xl4*yl2
         z(i,6) = -xl1*yl4-xl2*yl3+xl3*yl2+xl4*yl1
 10   continue

      do 20 i = 1, nin
         u(i,1,1) = cmplx(x(i,1), x(i,2))
         u(i,2,1) = cmplx(x(i,3), x(i,4))
         u(i,3,1) = cmplx(x(i,5), x(i,6))
         u(i,1,2) = cmplx(y(i,1), y(i,2))
         u(i,2,2) = cmplx(y(i,3), y(i,4))
         u(i,3,2) = cmplx(y(i,5), y(i,6))
         u(i,1,3) = cmplx(z(i,1), z(i,2))
         u(i,2,3) = cmplx(z(i,3), z(i,4))
         u(i,3,3) = cmplx(z(i,5), z(i,6))
  20  continue

      return

 100  continue
      do 50 i = 1, nin
         u(i,1,1) = cmplx(x(i,1), x(i,2))
         u(i,2,1) = cmplx(x(i,3), x(i,4))
         u(i,1,2) = cmplx(-x(i,3), x(i,4))
         u(i,2,2) = cmplx(x(i,1),-x(i,2))
  50  continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine store(nin,nd,uold,rolds,roldt,diso,
     &   pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &   pclo1p,pclo2p,rclo1p,rclo2p)
c------------------------------------------------------------------------c
c     store distances, orientations, overlaps and gauge interaction from c
c     commom /d/ in arrays uold,rolds, etc.                              c
c------------------------------------------------------------------------c

      parameter (n2=256)
      complex  ctr, utr, ctro, utro, uc1, uc2, uc3, uc4
      complex uc1o, uc2o, uc3o, uc4o, ofac, ofaco
      dimension uold(nd,nd), diso(nd,nd,4), pclo1(nd,nd)
      dimension rolds(nd,nd), roldt(nd,nd), pclo2(nd,nd)
      dimension uc1o(nd,nd), uc2o(nd,nd), uc3o(nd,nd)
      dimension uc4o(nd,nd), ofaco(nd,nd)
      dimension rclo1(nd,nd), rclo2(nd,nd)
      dimension pclo1p(nd,nd), pclo2p(nd,nd)
      dimension rclo1p(nd,nd), rclo2p(nd,nd)

      common /d/ uint(n2,n2),  dis(n2,n2,4), pcl1(n2,n2),
     & pcl2(n2,n2), r2s(n2,n2), r2t(n2,n2), ofac(n2,n2),
     & uc1(n2,n2), uc2(n2,n2), uc3(n2,n2), uc4(n2,n2),
     & rcl1(n2,n2), rcl2(n2,n2), pcl1p(n2,n2), pcl2p(n2,n2),
     & rcl1p(n2,n2), rcl2p(n2,n2)

      do 910 i = 1, nin
         do 911 ir = 1, 4
         do 911 j = 1, nin
            diso(j,i,ir) = dis(j,i,ir)
  911    continue
         do 910 j = 1, nin
            uold(j,i)  = uint(j,i)
            rolds(j,i) = r2s(j,i)
            roldt(j,i) = r2t(j,i)
            ofaco(j,i) = ofac(j,i)
            uc1o(j,i)  = uc1(j,i)
            uc2o(j,i)  = uc2(j,i)
            uc3o(j,i)  = uc3(j,i)
            uc4o(j,i)  = uc4(j,i)
            pclo1(j,i) = pcl1(j,i)
            pclo2(j,i) = pcl2(j,i)
            rclo1(j,i) = rcl1(j,i)
            rclo2(j,i) = rcl2(j,i)
            pclo1p(j,i) = pcl1p(j,i)
            pclo2p(j,i) = pcl2p(j,i)
            rclo1p(j,i) = rcl1p(j,i)
            rclo2p(j,i) = rcl2p(j,i)
  910 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--
                                                                       
      subroutine refresh(nin,nd,uold,rolds,roldt,diso,
     &   pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &   pclo1p,pclo2p,rclo1p,rclo2p)
c------------------------------------------------------------------------c
c     restore old distances, orientations, overlaps and gauge actions    c
c     from arrays uold, rolds, etc. to common /d/.                       c
c------------------------------------------------------------------------c

      parameter (n2=256)
      complex  ctr, utr, ctro, utro, uc1, uc2, uc3, uc4
      complex uc1o, uc2o, uc3o, uc4o, ofac, ofaco
      dimension uold(nd,nd), diso(nd,nd,4), pclo1(nd,nd)
      dimension rolds(nd,nd), roldt(nd,nd), pclo2(nd,nd)
      dimension uc1o(nd,nd), uc2o(nd,nd), uc3o(nd,nd)
      dimension uc4o(nd,nd), ofaco(nd,nd)
      dimension rclo1(nd,nd), rclo2(nd,nd)
      dimension pclo1p(nd,nd), pclo2p(nd,nd)
      dimension rclo1p(nd,nd), rclo2p(nd,nd)

      common /d/ uint(n2,n2),  dis(n2,n2,4), pcl1(n2,n2),
     & pcl2(n2,n2), r2s(n2,n2), r2t(n2,n2), ofac(n2,n2),
     & uc1(n2,n2), uc2(n2,n2), uc3(n2,n2), uc4(n2,n2),
     & rcl1(n2,n2), rcl2(n2,n2), pcl1p(n2,n2), pcl2p(n2,n2),
     & rcl1p(n2,n2), rcl2p(n2,n2)

      do 910 i = 1, nin
         do 911 ir = 1, 4
         do 911 j = 1, nin
            dis(j,i,ir) = diso(j,i,ir)
  911    continue
         do 910 j = 1, nin
            uint(j,i) = uold(j,i)
            r2s(j,i) = rolds(j,i)
            r2t(j,i) = roldt(j,i)
            uc1(j,i) = uc1o(j,i)
            uc2(j,i) = uc2o(j,i)
            uc3(j,i) = uc3o(j,i)
            uc4(j,i) = uc4o(j,i)
            ofac(j,i) = ofaco(j,i)
            pcl1(j,i) = pclo1(j,i)
            pcl2(j,i) = pclo2(j,i)
            rcl1(j,i) = rclo1(j,i)
            rcl2(j,i) = rclo2(j,i)
            pclo1p(j,i) = pcl1p(j,i)
            pclo2p(j,i) = pcl2p(j,i)
            rclo1p(j,i) = rcl1p(j,i)
            rclo2p(j,i) = rcl2p(j,i)
  910 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine iter(n,nin,nd,zr,e1,e2)
c------------------------------------------------------------------------c
c     perform one metropolis hit on every instanton. Need at least one   c
c     call to spectr() before iter() is called for the first time.       c
c     Change in action due to change in parameters of a single instanton c
c     is calculated in act(). Type of update controlled by parameter iop.c
c------------------------------------------------------------------------c
c     n,nin,nd   nc, number of instantons, array dimensions              c
c     zr,e1,e2   instanton configuration                                 c
c------------------------------------------------------------------------c
c     size in metropolis steps controlled by common/param/ ..sg,dz,drh,  c
c------------------------------------------------------------------------c
c     im controls updating strategy:                                     c
c     im = 1      metropolis algorithm                                   c
c     im = 2      cooling                                                c
c     im = 3      random model                                           c
c------------------------------------------------------------------------c

      parameter(n2=256)
      complex ctr, utr, ctro, utro, uc1o, uc2o, uc3o, uc4o
      complex ofaco
      dimension zr(nd,5), e1(nd,6), e2(nd,6), uold(n2,n2)
      dimension rolds(n2,n2), roldt(n2,n2)
      dimension x(6), y(6), zc(5), er1(6), er2(6) 
      dimension diso(n2,n2,4), ofaco(n2,n2)
      dimension uc1o(n2,n2), uc2o(n2,n2), uc3o(n2,n2), uc4o(n2,n2)
      dimension pclo1(n2,n2), pclo2(n2,n2)
      dimension rclo1(n2,n2), rclo2(n2,n2)
      dimension pclo1p(n2,n2), pclo2p(n2,n2)
      dimension rclo1p(n2,n2), rclo2p(n2,n2)

      common /param/al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /metr/ itot,ifail1,ifail2,ifail3,actt
      common /s1/   s1
      common /acti/ acold

      im    = 1
      itot  = 0
      ifail = 0

c------------------------------------------------------------------------c
c     loop over instantons                                               c
c------------------------------------------------------------------------c

      do 100 ip = 1, nin

c------------------------------------------------------------------------c
c     rho update, reject sizes larger than lambda^-1                     c
c------------------------------------------------------------------------c

         zrh = zr(ip,5)
 220     continue

         zr(ip,5) = zrh + drh*(rang( )-0.5)
c        zr(ip,5) = rang( )
         if(zr(ip,5) .le. 0.0 .or. zr(ip,5) .ge. 1.0) goto 220

c------------------------------------------------------------------------c
c     save old matrixelements                                            c
c------------------------------------------------------------------------c

         call store(nin,n2,uold,rolds,roldt,diso,
     &    pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &    pclo1p,pclo2p,rclo1p,rclo2p)

         acto = acold

c------------------------------------------------------------------------c
c     calculate change in action, iop=1 means only rho(ip) changed       c
c------------------------------------------------------------------------c

         iop=1
         call act(n,nin,nd,zr,e1,e2,acnew,ip,ip,iop)
         acold = acnew
         itot = itot + 1

c------------------------------------------------------------------------c
c     accept with propability exp(-(S_new-S_old))                        c
c------------------------------------------------------------------------c
c     xdec = 1 corresponds to cooling (->instanton crystal)              c
c------------------------------------------------------------------------c

         rat = exp((acnew-acto))
         if (im .eq. 1 .or. im .eq. 3) then
            xdec = rang()
         else if (im .eq. 2) then
            xdec = 1.0
         endif

c------------------------------------------------------------------------c
c     if metropolis hit is rejected, restore old matrixelements          c
c------------------------------------------------------------------------c

         if (rat .lt. xdec) then
            ifail1 = ifail1 + 1
            zr(ip,5) = zrh
            acold = acto
            call refresh(nin,n2,uold,rolds,roldt,
     &       diso,pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &       pclo1p,pclo2p,rclo1p,rclo2p)

         end if

c------------------------------------------------------------------------c
c     position update                                                    c
c------------------------------------------------------------------------c

         do 10 i = 1, 4
            zc(i) = zr(ip,i)
            rdm = rang( )
            zcc = zr(ip,i) + dz(i)*(rdm-.5)
c           zcc = al(i)*rdm
            if(zcc .gt. al(i)) zcc = zcc - al(i)
            if(zcc .lt.  0.0 ) zcc = zcc + al(i)
            zr(ip,i) = zcc
   10   continue

        call store(nin,n2,uold,rolds,roldt,diso,
     &   pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &   pclo1p,pclo2p,rclo1p,rclo2p)

        acto = acold

c------------------------------------------------------------------------c
c     calculate change in action, ip=2 means only zr(ip,mu) changed      c
c------------------------------------------------------------------------c

         iop=2
         call act(n, nin, nd, zr, e1, e2,acnew,ip,ip,iop)
         acold = acnew
         itot = itot + 1

c------------------------------------------------------------------------c
c        accept with propability exp(-(S_new-S_old))                     c
c------------------------------------------------------------------------c

         rat = exp((acnew-acto))
         if (im .eq. 1) then
            xdec = rang()
         else if (im .eq. 2) then
            xdec = 1.0
         else if (im .eq. 3) then
            xdec = 0.0
         endif

         if (rat .lt. xdec) then
            ifail2 = ifail2 + 1
            do 60 i = 1, 4
               zr(ip,i) = zc(i)
   60       continue
            acold = acto
            call refresh(nin,n2,uold,rolds,roldt,diso,
     &       pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &       pclo1p,pclo2p,rclo1p,rclo2p)

         end if

c------------------------------------------------------------------------c
c     orientation update                                                 c
c------------------------------------------------------------------------c

         do 5 i = 1, 6
            x(i) = e1(ip,i)
            y(i) = e2(ip,i)
            er1(i) = e1(ip,i)
            er2(i) = e2(ip,i)
   5     continue
         nu2 = n*2
         call rgu(nu2,er1,er2,sg)
c        call rsu(nu2,er1,er2)
         do 6 i = 1, 6
            e1(ip,i) = er1(i)
            e2(ip,i) = er2(i)
   6     continue
         call store(nin,n2,uold,rolds,roldt,
     &    diso,pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &    pclo1p,pclo2p,rclo1p,rclo2p)
         acto = acold

c------------------------------------------------------------------------c
c     calculate change in action, iop=3 means only e1,e2(6,ip) updated   c
c------------------------------------------------------------------------c

         iop=3
         call act(n, nin, nd, zr, e1, e2,acnew,ip,ip,iop)
         acold = acnew
         itot = itot + 1

c------------------------------------------------------------------------c
c     accept with propability exp(-(S_new-S_old))                        c
c------------------------------------------------------------------------c

         rat = exp((acnew-acto))
         if (im .eq. 1) then
            xdec = rang()
         else if (im .eq. 2) then
            xdec = 1.0
         else if (im .eq. 3) then
            xdec = 0.0
         endif

         if (rat.lt.xdec) then
            ifail3 = ifail3 + 1
            do 70 i = 1, 6
               e1(ip,i) = x(i)
               e2(ip,i) = y(i)
   70       continue
            acold = acto
            call refresh(nin,n2,uold,rolds,roldt,diso,
     &       pclo1,pclo2,rclo1,rclo2,uc1o,uc2o,uc3o,uc4o,ofaco,
     &       pclo1p,pclo2p,rclo1p,rclo2p)

         end if

c------------------------------------------------------------------------c
c      end of big loop over instantons                                   c
c------------------------------------------------------------------------c

  100 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine act(n,nin,nd,zr,e1,e2,sdet,il,iu,iop)
c------------------------------------------------------------------------c
c     calculate total action of instanton configuration. iop controls    c
c     mode of operation: full calculation or update only.                c
c------------------------------------------------------------------------c
c     input: n,nin,nd    nc, number of instantons, array dimension       c
c            zr,e1,e2    instanton configuration                         c
c            sdet        log of weight function, log(det)-S              c
c            il,iu       index of first, last instanton to be updated    c
c            iop         mode of operation (see below)                   c
c------------------------------------------------------------------------c
c     iop = 0 :  full calculation                                        c
c     iop = 1 :  only rho(il),..,rho(iu) is updated                      c
c     iop = 2 :  only positions updated                                  c
c     iop = 3 :  only orientations updated                               c
c------------------------------------------------------------------------c

      parameter(ni=256,ni2=128)
      complex clp, cdet, r, cdu, cds, cl2, clt, clf, cjc, fac
      complex clm, clr, ci, cl4
c     complex det1
      dimension clp(ni2,ni2), zr(nd,5), r(1000), cl2(ni2,ni2)
      dimension clf(ni2,ni2), cjc(ni2,ni2), fac(ni2,ni2),iptv(ni2)
      dimension ar(ni2,ni2), ai(ni2,ni2), wr(ni2), wi(ni2)
      dimension e1(nd,6), e2(nd,6)
      dimension clm(ni,ni), clr(ni,ni), cl4(ni,ni)

      common /param/ al(4), alpha,rh0,sg,dz(4),drh, nc, nf, rmu, rms
      common /pi/ pi, eps
      common /musim/ xmu, phase, imu

c------------------------------------------------------------------------c
c     calculate gauge field action and fermionic overlap matrix elements c
c------------------------------------------------------------------------c
c     if iop .ne. 0, only updates performed                              c
c------------------------------------------------------------------------c

      nih = nin/2
      nd2 = nd/2
      ci  = (0.0,1.0)
      call overl(n,nin,nd,nd2,zr,e1,e2,clp,clm,cl4,
     2               bofac,rhofc,uintt,il,iu,iop)

      rmul = 0.0
      do 15 i = 1, nih
         rhp  = zr(i,5)*zr(i+nih,5)
         rmul = rmul + alog(rhp)
   15 continue

c------------------------------------------------------------------------c
c     calculate effective action                                         c
c------------------------------------------------------------------------c

      if (nf .eq. 0) then

         sdet = bofac

      else if (imu .eq. 1) then

c------------------------------------------------------------------------c
c     fermionic determinant, (nf-1) light quarks, one heavy quark        c
c------------------------------------------------------------------------c
c     start with light quark                                             c
c------------------------------------------------------------------------c

         do 10 i = 1, nih
         do 10 j = 1, nih
            cjc(j,i) = conjg(clp(j,i))
   10    continue

c------------------------------------------------------------------------c
c     clt = T*T^(+), clf = T*T^(+)+m^2                                   c
c------------------------------------------------------------------------c

         do 12 i = 1, nih
         do 12 j = 1, nih
            clt = cmplx(0.0,0.0)
            do 14 k = 1, nih
               clt = clt + clp(k,i)*cjc(k,j)
   14       continue
            cl2(j,i) = clt
   12    continue

         do 18 i = 1, nih
            do 28 j = 1, nih
               clf(j,i) = cl2(j,i)
   28       continue
            clf(i,i) = cl2(i,i) + rmu*rmu
   18    continue

c------------------------------------------------------------------------c
c     determinant for u quark                                            c
c------------------------------------------------------------------------c

         call logdet(nih,clf,aludet)

c------------------------------------------------------------------------c
c     same for strange quark                                             c
c------------------------------------------------------------------------c

         do 19 i = 1, nih
            do 29 j = 1, nih
              clf(j,i) = cl2(j,i)
   29       continue
            clf(i,i) = cl2(i,i) +rms*rms
   19    continue

         call logdet(nih,clf,alsdet)

c------------------------------------------------------------------------c
c     collect light and heavy quark contributions to sdet = log(det)     c
c------------------------------------------------------------------------c

         sdet = (nf-1)*aludet + alsdet + nf*rmul
         rdet = sdet
         sdet = sdet + bofac

      else 

c------------------------------------------------------------------------c
c     same stuff for non-hermitean matrix with chemical potential        c
c------------------------------------------------------------------------c

         do 51 i=1,nih
         do 51 j=1,nih
            cl2(i,j) = (0.0,0.0)
            do 52 k=1,nih
               cl2(i,j) = cl2(i,j) + clm(i,k+nih)*clm(k+nih,j)
 52         continue
 51      continue

         do 53 i=1,nih
            do 55 j=1,nih
               clf(i,j) = cl2(i,j)
 55         continue
            clf(i,i) = cl2(i,i) + rmu*rmu 
 53      continue            

c------------------------------------------------------------------------c
c     determinant for u quark                                            c
c------------------------------------------------------------------------c

         call clogdet(nih,clf,aludet,phaseu)

c------------------------------------------------------------------------c
c     same for strange quark                                             c
c------------------------------------------------------------------------c

         do 63 i=1,nih
            do 65 j=1,nih
               clf(i,j) = cl2(i,j)
 65         continue
            clf(i,i) = cl2(i,i) + rms*rms 
 63      continue            

         call clogdet(nih,clf,alsdet,phases)

c------------------------------------------------------------------------c
c     collect light and heavy quark contributions to sdet = log(det)     c
c------------------------------------------------------------------------c

         sdet = (nf-1)*aludet + alsdet + nf*rmul
         rdet = sdet
         sdet = sdet + bofac
         phase= (nf-1)*phaseu + phases

      end if

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine logdet(nih,clf,aldet)
c------------------------------------------------------------------------c
c     calculate log(det(clf)) using imsl/lapack subroutines              c
c------------------------------------------------------------------------c
c     nih       dimension of clf                                         c
c     clf       complex hermitean matrix clf(nih,nih)                    c
c     aldet     log(det(clf))                                            c
c------------------------------------------------------------------------c
c     lapack workspace requirements as follows :                         c
c     lwork = 2*n-1, lrwork = 3*n-2                                      c
c     careful, seem to need:                                             c
c     lwork = 3*n,   lrwork = 4*n                                        c
c------------------------------------------------------------------------c

      parameter(ni=256, ni2=128, lwork=512, lrwork=512)
      complex clf, fac
      complex work
      dimension clf(ni2,ni2), fac(ni2,ni2),iptv(ni2)
      dimension work(lwork), rwork(lrwork), wr(ni2)
      complex*16 zwork(lwork), zlf(ni2,ni2)
      real*8     zrwork(lwork), zwr(ni2)

c------------------------------------------------------------------------c
c     imsl version                                                       c
c------------------------------------------------------------------------c

c     call lfthf(nih,clf,ni2,fac,ni2,iptv)
c     call lfdhf(nih,fac,ni2,iptv,det1,det2)

c     aldet = alog(det1) + det2*alog(10.0)

c------------------------------------------------------------------------c
c     lapack version                                                     c
c------------------------------------------------------------------------c

c     call cheev('n','u',nih,clf,ni2,wr,work,lwork,rwork,info)
c     if(info .ne. 0) stop

c     aldet = 0.0
c     do 5 i=1,nih
c        aldet = aldet + alog(abs(wr(i)))
c 5   continue

c------------------------------------------------------------------------c
c     diagonalize, lapack double precision                               c
c------------------------------------------------------------------------c

      do 3 i=1,nih
      do 3 j=1,nih
           zlf(i,j) = clf(i,j)
  3   continue

      call zheev('n','u',nih,zlf,ni2,zwr,zwork,lwork,zrwork,info)
      if(info .ne. 0) stop

      aldet = 0.0
      do 15 i=1,nih
         aldet = aldet + log(abs(zwr(i)))
  15  continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine logdetnh(nih,clf,aldet)
c------------------------------------------------------------------------c
c     calculate log(det(clf)) using imsl/lapack subroutines.             c
c------------------------------------------------------------------------c
c     uses the fact that log(det(D))=log(det(TT^+))                      c
c     Note: factor 2 supplied in calling routine                         c
c------------------------------------------------------------------------c
c     nih       dimension of clf                                         c
c     clf       complex matrix clf(nih,nih)                              c
c     aldet     log(det(clf))                                            c
c------------------------------------------------------------------------c
c     lapack workspace requirements as follows :                         c
c     lwork = 2*n-1, lrwork = 3*n-2                                      c
c     careful, seem to need:                                             c
c     lwork = 3*n,   lrwork = 4*n                                        c
c------------------------------------------------------------------------c

      parameter(ni=256, ni2=128, lwork=512, lrwork=512)
      complex clf, fac, det1, wr, vl, vr
      complex work
      dimension clf(ni2,ni2), fac(ni2,ni2),iptv(ni2)
      dimension work(lwork), rwork(lrwork), wr(ni2)
      dimension vl(ni2,ni2), vr(ni2,ni2)
      real*8     zrwork(lrwork)
      complex*16 zwork(lwork),zwr(ni2),zlf(ni2,ni2)
      complex*16 zvl(ni2,ni2),zvr(ni2,ni2)

c------------------------------------------------------------------------c
c     imsl version                                                       c
c------------------------------------------------------------------------c

c     call lftcg(nih,clf,ni2,fac,ni2,iptv)
c     call lfdcg(nih,fac,ni2,iptv,det1,det2)

c     aldet = alog(cabs(det1)) + det2*alog(10.0)

c------------------------------------------------------------------------c
c     lapack version                                                     c
c------------------------------------------------------------------------c

c     call cgeev('n','n',nih,clf,ni2,wr,vl,ni2,vr,ni2,
c    1           work,lwork,rwork,info)
c     if(info .ne. 0) stop

c     aldet = 0.0
c     do 5 i=1,nih
c        aldet = aldet + alog(cabs(wr(i)))
c 5   continue

c------------------------------------------------------------------------c
c     lapack double precision                                            c
c------------------------------------------------------------------------c

      do 3 i=1,nih
      do 3 j=1,nih
           zlf(i,j) = clf(i,j)
  3   continue
     
      call zgeev('n','n',nih,zlf,ni2,zwr,zvl,ni2,zvr,ni2,
     1           zwork,lwork,zrwork,info)
      if(info .ne. 0) stop

      aldet = 0.0
      do 15 i=1,nih
         aldet = aldet + log(abs(zwr(i)))
 15   continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--
                                                                       
      subroutine clogdet(nih,clf,aldet,phase)
c------------------------------------------------------------------------c
c     calculate log(det(clf)) using imsl/lapack subroutines.             c
c------------------------------------------------------------------------c
c     nin       dimension of clf                                         c
c     clf       complex matrix clf(ni,ni)                                c
c     aldet     log(|det(clf)|)                                          c
c     phase     arg(det(clf))                                            c
c------------------------------------------------------------------------c
c     lapack workspace requirements as follows :                         c
c     lwork = 2*n-1, lrwork = 3*n-2                                      c
c     careful, seem to need:                                             c
c     lwork = 3*n,   lrwork = 4*n                                        c
c------------------------------------------------------------------------c

      parameter(ni=256, ni2=128, lwork=512, lrwork=512)
      complex   clf
      dimension clf(ni2,ni2)
      real*8     zrwork(lrwork)
      complex*16 zwork(lwork), zwr(ni2), zlf(ni2,ni2)
      complex*16 zvl(ni2,ni2), zvr(ni2,ni2)

      common /pi/    pi, eps

c------------------------------------------------------------------------c
c     lapack double precision                                            c
c------------------------------------------------------------------------c

      do 3 i=1,nih
      do 3 j=1,nih
           zlf(i,j) = clf(i,j)
  3   continue
    
      call zgeev('n','n',nih,zlf,ni2,zwr,zvl,ni2,zvr,ni2,
     1           zwork,lwork,zrwork,info)
      if(info .ne. 0) stop

      aldet = 0.0
      phase = 0.0
      do 15 i=1,nih
         aldet = aldet + log(abs(zwr(i)))
         phase = phase + dimag(cdlog(zwr(i)))
 15   continue

c     phase = phase - nih*pi

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine spect(n,nin,nd,zr,e1,e2,rmdet,wr,zuwr,sdet,rv,uin,rlp,
     2            rmulp,rcos,rcii,rnii,rnia,r4ii,r4ia,rdii,rdia,mol,p)
c------------------------------------------------------------------------c
c     calculate spectrum, complete measure log(det)-S, and distribution  c
c     of relative orientations for instanton configuration. subroutine   c
c     assumes that configuration is calculated for the first time.       c
c------------------------------------------------------------------------c
c     n,nin,nd    number of colors, number of inst, array dimension      c
c     zr,e1,e2    instanton configuration                                c
c     rmdet       fermionic determinant                                  c
c     wr(nin/2)   fermionic eigenvalues                                  c
c     zuwr(nin)   complex eigenvalues                                    c   
c     sdet        log(measure)                                           c
c     rv          log of density distribution                            c
c     uin         gauge interaction                                      c
c     rlp(i,j)    absolute value of IA overlap matrix elements           c
c     rmulp(i,j)  same for non-hermitean part                            c
c     rcos(i,j)   cos^2 of relative IA angle                             c
c     rcii(i,j)   cos^2 of relative II angle                             c
c     rnii(i,j)   cos^2 of spatial IA angle                              c
c     rnia(i,j)   cos^2 of spatial II angle                              c
c     r4ii(i,j)   |u_4|^2/|u|^2 for IA                                   c
c     r4ia(i,j)   |u_4|^2/|u|^2 for II                                   c
c     mol(j)      I=mol(J) such that T_IJ maximal                        c
c------------------------------------------------------------------------c
c     lapack workspace requirements as follows :                         c
c     lwork = 2*n-1, lrwork = 3*n-2                                      c
c     careful, seem to need more than that:                              c
c     lwork = 3*n,   lrwork = 4*n                                        c
c------------------------------------------------------------------------c
c     In musim dirac operator is diagonalized both with (->zuwr(nin))    c
c     and without (-> wr(1,..,nih)) the chemical potential. For imu=1    c
c     determinant and effectice action are calculated from wr, else      c
c     from zwr.                                                          c
c------------------------------------------------------------------------c 

      parameter(ni=256, ni2=128, n2=256,lwork=512,lrwork=512)
      complex   clp, clt, cjc, cl2, ofac
      complex   ctr, utr, uc1, uc2, uc3, uc4, eval(ni2)
      complex   ctria, ctrii, ct4ia, ct4ii
      complex   work                          
      complex   clm, ci, ceig
      dimension clp(ni2,ni2), zr(nd,5)  
      dimension clm(ni,ni)
      dimension cl2(ni2,ni2), cjc(ni2,ni2)
      dimension rlp(ni2,ni2), rcos(ni2,ni2)
      dimension rnii(ni2,ni2),rnia(ni2,ni2)
      dimension r4ii(ni2,ni2),r4ia(ni2,ni2)
      dimension ruii(ni2,ni2),ruia(ni2,ni2)
      dimension rdii(ni2,ni2),rdia(ni2,ni2)
      dimension rcii(ni2,ni2),rmulp(ni2,ni2) 
      dimension mol(ni2), p(ni2)
      dimension ar(ni2,ni2), ai(ni2,ni2), wr(ni2), wi(ni2)
      dimension e1(nd,6), e2(nd,6)
      dimension work(lwork), rwork(lrwork)
      real*8     zrwork(lwork), zwr(ni2)
      complex*16 zwork(lwork), zl2(ni2,ni2)
      complex    clr(ni2,ni2), zuwr(ni2), zswr(ni2)
      complex*16 zlr(ni2,ni2), zdwr(ni2)        
      complex*16 zvl(ni2,ni2), zvr(ni2,ni2)
      complex    zzv(ni,ni)

c------------------------------------------------------------------------c
c     note : u1c,.. is renamed as uc1,...,uc4                            c
c------------------------------------------------------------------------c

      common /d/     uint(n2,n2),  dis(n2,n2,4), pcl1(n2,n2),
     2 pcl2(n2,n2), r2s(n2,n2), r2t(n2,n2), ofac(n2,n2),
     3 uc1(n2,n2), uc2(n2,n2), uc3(n2,n2), uc4(n2,n2),
     4 rcl1(n2,n2), rcl2(n2,n2)
      common /param/ al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /pi/    pi, eps
      common /sij/   sij(ni,ni)
      common /musim/ xmu, phase, imu

c------------------------------------------------------------------------c
c     do full calculation (iop=0) of overlaps and bosonic measure        c
c------------------------------------------------------------------------c

      nih = nin/2
      nd2 = nd/2
      nar = nin
      ci  = (0.0,1.0)
      call overl(n,nin,nd,nd2,zr,e1,e2,clp,clm,cl4,
     &             bofac,rv,uintt,1,nar,0)
      uin = uintt

c------------------------------------------------------------------------c
c     for statistics only: abs(clp) and IA orientation angles            c
c------------------------------------------------------------------------c

      do 10 i = 1, nih
         mol(i) =  1
         cmax   = 0.0
         do 10 j = 1, nih
            jp = j+nih

c------------------------------------------------------------------------c
c     orientation invariants for (IJ) and (IA)                           c
c------------------------------------------------------------------------c

             rr2ia = dis(jp,i,1)**2 + dis(jp,i,2)**2
     2              +dis(jp,i,3)**2 + dis(jp,i,4)**2
             rr2ii = dis(j,i,1)**2  + dis(j,i,2)**2
     2              +dis(j,i,3)**2  + dis(j,i,4)**2
             ctria = ofac(jp,i) + dis(jp,i,4)*uc4(jp,i)
             ctrii = ofac(j,i)  + dis(j,i,4) *uc4(j,i)
             ctuia =(cabs(uc1(jp,i))**2 + cabs(uc2(jp,i))**2
     2              +cabs(uc3(jp,i))**2 + cabs(uc4(jp,i))**2)
             ctuii =(cabs(uc1(j,i))**2  + cabs(uc2(j,i))**2
     2              +cabs(uc3(j,i))**2  + cabs(uc4(j,i))**2)
             ct4ia = uc4(jp,i)
             ct4ii = uc4(j,i)

c------------------------------------------------------------------------c
c     color and spatial orientation angles                               c
c------------------------------------------------------------------------c

             rdia(j,i) = sqrt(rr2ia)
             rdii(j,i) = sqrt(rr2ii)
             rcos(j,i) = cabs(ctria)**2/ctuia/rr2ia
             rcii(j,i) = cabs(ctrii)**2/ctuii/(rr2ii+eps)
             rnia(j,i) = dis(jp,i,4)**2/rr2ia
             rnii(j,i) = dis(j,i,4)**2/(rr2ii+eps)
             ruia(j,i) = ctuia
             ruii(j,i) = ctuii
             r4ia(j,i) = cabs(ct4ia)**2/ctuia
             r4ii(j,i) = cabs(ct4ii)**2/ctuii
             dipole    = ctuia-4.0*rcos(j,i)

c------------------------------------------------------------------------c
c     overlaps and pairs                                                 c
c------------------------------------------------------------------------c

             cjc(i,j)  = conjg(clp(i,j))
             rlp(j,i)  = cabs(clp(i,j))
             rmulp(j,i)= cabs(clm(i,j+nih)-clp(i,j))
             if (rlp(j,i) .gt. cmax) then
                cmax   = rlp(j,i)
                mol(i) = j
             endif

   10 continue

c------------------------------------------------------------------------c
c     calculate clt = T*T^(+)                                            c
c------------------------------------------------------------------------c

      do 12 i = 1, nih
      do 12 j = 1, nih
         clt = cmplx(0.0,0.0)
         do 14 k = 1, nih
            clt = clt + clp(k,i)*cjc(k,j)
   14    continue
         cl2(j,i) = clt
         zl2(j,i) = clt
   12 continue

c------------------------------------------------------------------------c
c     diagonalize T*T^(+), imsl                                          c
c------------------------------------------------------------------------c

c     call evlhf(nih,cl2,ni2,wr)

c------------------------------------------------------------------------c
c     diagonalize, lapack                                                c
c------------------------------------------------------------------------c

c     call cheev('n','l',nih,cl2,ni2,wr,work,lwork,rwork,info)
c     if(info .ne. 0) stop

c------------------------------------------------------------------------c
c     lapack, double precision                                           c
c------------------------------------------------------------------------c
 
      call zheev('n','l',nih,zl2,ni2,zwr,zwork,lwork,zrwork,info)
      if(info .ne. 0) stop

      do 7 i=1,nih
           wr(i) = zwr(i)
  7   continue

c------------------------------------------------------------------------c
c     calculate log(det) for (nf-1) light plus one heavy flavor          c
c------------------------------------------------------------------------c
      
      if (imu .eq. 1) then 
       
         rdet = 0.0
         sdet = 0.0
         srho = 0.0
         do 20 i = 1, nih
            rhp  = zr(i,5)*zr(i+nih,5)
            rms2 = rms*rms
            rmu2 = rmu*rmu
            eig  = abs(wr(i))
            if(wr(i) .lt. 0) write(6,*) 'warning: lam <0'
            srho = srho+alog(rhp)
            sdet = sdet+(nf-1)*alog(eig+rmu2) + alog(rms2+eig)
            rdet = rdet+(nf-1)*alog(eig) + alog(eig)
   20    continue
         sdet = sdet + nf*srho
         rdet = rdet + nf*srho

c------------------------------------------------------------------------c
c      calculate log(measure)                                            c
c------------------------------------------------------------------------c

         if (nf.eq.0) then
            rmdet = 0.0
         else
            rmdet = exp(sdet/(nf*nin))
         end if

         if (nf .eq. 0) then
           sdet = bofac
         else
           sdet = sdet + bofac
         end if              
      
      endif

      do 30 i = 1, nih
         wr(i) = sqrt(abs(wr(i)))
   30 continue

c------------------------------------------------------------------------c
c     same for non-hermitean matrix with chemical potential              c
c------------------------------------------------------------------------c

         do 51 i=1,nih
         do 51 j=1,nih
            clr(i,j) = (0.0,0.0)
            do 52 k=1,nih
               clr(i,j) = clr(i,j) + clm(i,k+nih)*clm(k+nih,j)
 52         continue
            zlr(i,j) = clr(i,j) 
 51      continue

      call zgeev('n','v',nih,zlr,ni2,zdwr,zvl,ni2,zvr,ni2,
     1           zwork,lwork,zrwork,info)
      if(info .ne. 0) stop

      do 57 i=1,nin
           zuwr(i) = sqrt(zdwr(i))
  57  continue

c------------------------------------------------------------------------c
c     calculate participation ratios                                     c
c------------------------------------------------------------------------c

      do 58 i=1,nih
         p(i) = 0.0
         test = 0.0 
         do 59 j=1,nih
            p(i) = p(i) + abs(zvr(j,i))**4
            test = test + abs(zvr(j,i))**2
  59     continue
         p(i) = 1.0/p(i)
  58  continue 

c------------------------------------------------------------------------c
c     reconstruct full eigenvector, calculate participation ratio        c
c------------------------------------------------------------------------c

      do 61 i=1,nih
         do 65 j=1,nih
            zzv(j,i) = zvr(j,i)
  65     continue
         do 66 j=1,nih
            jp = j+nih
            zzv(jp,i) = (0.d0,0.d0)
            do 67 k=1,nih
               zzv(jp,i) = zzv(jp,i)+clm(nih+j,k)*zvr(k,i)
  67        continue
            zzv(jp,i) = zzv(jp,i)/zuwr(i)
  66     continue
         xnorm2 = 0.0
         do 68 j=1,2*nih
            xnorm2 = xnorm2 + abs(zzv(j,i))**2
  68     continue
         xnorm = sqrt(xnorm2)
         do 69 j=1,2*nih
            zzv(j,i) = zzv(j,i)/xnorm
  69     continue
  61  continue

      do 78 i=1,nih
         pnew = 0.0
         test = 0.0 
         do 79 j=1,2*nih
            pnew = pnew + abs(zzv(j,i))**4
            test = test + abs(zzv(j,i))**2
  79     continue
         p(i) = 1.0/pnew
  78  continue 

c------------------------------------------------------------------------c
c     calculate log(|det|) for non-hermitean matrix                      c
c------------------------------------------------------------------------c
                        
      if (imu .ne. 1) then
                        
         rdet = 0.0
         sdet = 0.0
         srho = 0.0
         phase= 0.0
         do 70 i = 1, nih
            rhp  = zr(i,5)*zr(i+nih,5)
            ceig = zdwr(i)
            phase= phase + (nf-1)*aimag(log(ceig+rmu**2))
     1                   + aimag(log(ceig+rms**2))
            srho = srho  + alog(rhp)
            sdet = sdet  + (nf-1)*alog(abs(ceig)+rmu**2) 
     1                   + alog(abs(ceig)+rms**2)
            rdet = rdet  + nf*alog(abs(ceig)) 
   70    continue
         sdet = sdet + nf*srho
         rdet = rdet + nf*srho
c        phase= phase- nf*nin/2.0*pi

c------------------------------------------------------------------------c
c     calculate log(measure)                                             c
c------------------------------------------------------------------------c

         if (nf.eq.0) then
            rmdet = 0.0
         else
            rmdet = exp(sdet/(nf*nin))
         end if

         if (nf .eq. 0) then
            sdet = bofac
         else
            sdet = sdet + bofac
         endif

c------------------------------------------------------------------------c
c     end of imu block if                                                c
c------------------------------------------------------------------------c

      endif    
 
      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine dens(n,nin,nd,zr,e1,e2,crho)
c------------------------------------------------------------------------c
c     calculate baryon density. 
c------------------------------------------------------------------------c
c     n,nin,nd    number of colors, number of inst, array dimension      c
c     zr,e1,e2    instanton configuration                                c
c     rho         baryon density                                         c
c------------------------------------------------------------------------c

      parameter(ni=256, ni2=128, n2=256,lwork=2048)
      complex   clp, clm, cl4, cli, ci, ceig, crho
      dimension clp(ni2,ni2), zr(nd,5)  
      dimension clm(ni,ni), cl4(ni,ni), cli(ni,ni)
      dimension work(lwork), ipiv(ni)
  
      common /param/ al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /pi/    pi, eps
      common /musim/ xmu, phase, imu

c------------------------------------------------------------------------c
c     do full calculation (iop=0) of overlaps and bosonic measure        c
c------------------------------------------------------------------------c

      nih = nin/2
      nd2 = nd/2
      nar = nin
  
      call overl(n,nin,nd,nd2,zr,e1,e2,clp,clm,cl4,bfc,rv,ui,1,nar,0)
     
c-----------------------------------------------------------------------c
c     include mass term                                                 c
c-----------------------------------------------------------------------c

      do 10 i=1,nin
      do 10 j=1,nin
         cli(i,j) = (0.0,0.0)
  10  continue

      do 20 i=1, nin
         cli(i,i) = cmplx(0.0,-rmu)
         do 20 j=1, nin
            cli(i,j) = cli(i,j) + clm(i,j)
  20  continue

c-----------------------------------------------------------------------c
c     calculate inverse, use LAPACK                                     c
c-----------------------------------------------------------------------c

      call cgetrf(nin,nin,cli,ni,ipiv,info)
      if(info .ne. 0) stop
      call cgetri(nin,cli,ni,ipiv,work,lwork,info)
      if(info .ne. 0) stop

c-----------------------------------------------------------------------c
c     baryon density, rho= tr(S\gam_4)                                  c
c-----------------------------------------------------------------------c

      crho = (0.0,0.0)
      do 50 i=1, nin
      do 50 j=1, nin
         crho = crho + cli(i,j)*cl4(j,i)
  50  continue

      return
      end 

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      subroutine lev(xmin,st,n,nsm,ist)
c------------------------------------------------------------------------c
c     plot simple histogram in output file ftn02                         c
c------------------------------------------------------------------------c
c     input: xmin    smallest x value                                    c
c            st      bin width                                           c
c            n       number of bins                                      c
c            nsm     plot symbol                                         c
c            ist(n)  histogram array                                     c
c------------------------------------------------------------------------c

      dimension ist(n),a(6),k(6),t(50)
      data a    /1hx,1h*,1h+,1hc,1he,1h /
      data nhgt /50/

      j=0
      m=0
      s=0.
      d=xmin-st*0.5
      x=d
      do 1 i=1,n
      x=x+st
      s=s+ist(i)*x
      if(ist(i).lt.0)go to 9
      if(i.eq.1.or.i.eq.n)go to 1
      if(ist(i).gt.m)m=ist(i)
    1 j=j+ist(i)
      if(j.lt.1)go to 10
      s=s/j
      m=m/nhgt+1
      do 2 i=1,6
    2 k(i)=m*10*(i-1)
      write(2, 30)k
   30 format(5x,i20,5i10)
      write(2,40)
   40 format (2x,3h---,7(9(1h-),1h+))
      x=d
      d=0.
      do6 i=1,n
      do3 l=1,nhgt
    3 t(l)=a(6)
      nn=ist(i)/m
      if(nn.eq.0) go to 5
         if(nn.gt.nhgt)nn=nhgt
      do 4 l=1,nn
    4 t(l)=a(nsm)
    5 x=x+st
      write(2, 70)i,ist(i),x,t
   70 format(1x,1hi,i4,i7,1pe11.3,1x,50a1,1hi)
    6 d=d+ist(i)*(x-s)**2
      if(j.lt.2)goto 7
      d=sqrt(d/(j-1))
    7 write(2, 40)
      write(2, 30)k
   8  write(2, 90)j,s,d
   90 format(//2x,18hnumber of events =,i8,5x,
     *10haverage = ,e11.4,5x,8hsigma = ,e11.4//)
      return
    9 write(2, 80)i
   80 format(//2x,3hlev,i8,23h channel  less  than  0//)
      return
10    d=0.
      go to 8
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine lens(a,amin,st,m,ist)
c------------------------------------------------------------------------c
c     include value a in histogram array ist(n)                          c
c------------------------------------------------------------------------c
c     a      value to be added to histogram array                        c
c     amin   minimum value in histogram                                  c
c     st     bin width                                                   c
c     m      number of bins                                              c
c     ist(n) histogram array                                             c
c------------------------------------------------------------------------c

      dimension ist(150)
      j=(a-amin)/st+1.000001
      if(j.lt.1)j=1
      if(j.gt.m)j=m
      ist(j)=ist(j)+1
      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine tdlens(ax,ay,axmin,aymin,stx,sty,mx,my,ist,ndim)
c------------------------------------------------------------------------c
c     include value a in 2d histogram array ist(ndim,ndim)               c
c------------------------------------------------------------------------c
c     ax,ay       value to be added to histogram array                   c
c     axmin,aymin minimum value in histogram                             c
c     stx,sty     bin width                                              c
c     mx,my       number of bins                                         c
c     ist(nx,ny)  histogram array                                        c
c------------------------------------------------------------------------c

      dimension ist(ndim,ndim)

      jx=(ax-axmin)/stx+1.000001
      if(jx.lt.1)  jx=1
      if(jx.gt.mx) jx=mx
 
      jy=(ay-aymin)/sty+1.000001
      if(jy.lt.1)  jy=1
      if(jy.gt.my) jy=my

      ist(jx,jy)=ist(jx,jy)+1

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine ctdlens(ax,ay,c,axmin,aymin,stx,sty,mx,my,cist,ndim)
c------------------------------------------------------------------------c
c     same with complex phase c                                          c
c------------------------------------------------------------------------c
c     ax,ay       value to be added to histogram array                   c
c     axmin,aymin minimum value in histogram                             c
c     stx,sty     bin width                                              c
c     mx,my       number of bins                                         c
c     ist(nx,ny)  histogram array                                        c
c------------------------------------------------------------------------c

      complex c,cist(ndim,ndim)

      jx=(ax-axmin)/stx+1.000001
      if(jx.lt.1)  jx=1
      if(jx.gt.mx) jx=mx
 
      jy=(ay-aymin)/sty+1.000001
      if(jy.lt.1)  jy=1
      if(jy.gt.my) jy=my

      cist(jx,jy)=cist(jx,jy)+c

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine binadd(f,x,farray,f2array,narray,xmin,xbin,nbin)
c----------------------------------------------------------------------------c
c     add measurement f(x) to averages summed in array farray. variable x is c
c     binned in xmin,xmin+xbin,...,xmin+nbin*xbin.                           c
c----------------------------------------------------------------------------c
c     f     result of measurement                                            c
c     x     value of variable for this measurement                           c
c     farr  array of fuction values in bin k                                 c
c     f2arr array of function values squared                                 c
c     narr  array of number of measurements in bin k                         c
c     xmin  smallest x bin                                                   c
c     xbin  binsize                                                          c
c     nbin  number of bins                                                   c
c----------------------------------------------------------------------------c
      real farray(nbin),f2array(nbin)
      integer narray(nbin)

      kx = int((x-xmin)/xbin)+1
      kx = min(kx,nbin)
      kx = max(kx,1)
 
      farray(kx) = farray(kx) + f
      f2array(kx)= f2array(kx)+ f**2
      narray(kx) = narray(kx) + 1

      return
      end

c----------------------------------------------------------------------+-----
c----------------------------------------------------------------------+-----

      subroutine binav(farray,f2array,narray,nbin,fav,ferr)
c----------------------------------------------------------------------------c
c     determine average and uncertainity of binned function                  c
c----------------------------------------------------------------------------c
c     farr  array of fuction values in bin k                                 c
c     f2arr array of function values squared                                 c
c     narr  array of number of measurements in bin k                         c
c     nbin  number of bins                                                   c
c     fav   array of average function values                                 c
c     ferr  array of uncertainties                                           c
c----------------------------------------------------------------------------c
      real farray(nbin),f2array(nbin),fav(nbin),ferr(nbin)
      integer narray(nbin)

      do 10 i=1,nbin
         if(narray(i) .ne. 0) then
         fav(i) = farray(i)/float(narray(i))
         ferr(i)= disp(narray(i),farray(i),f2array(i))
         else 
         fav(i) = 0.0
         ferr(i)= 1.0
         endif
 10   continue

      return 
      end
 
c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine zero(n,ia)
c------------------------------------------------------------------------c
c     clear integer array ia(n)                                          c
c------------------------------------------------------------------------c

      dimension ia(n)

      do 10 i = 1, n
        ia(i) = 0.0
   10 continue

      return
      end      

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine zero2(nx,ny,ia)
c------------------------------------------------------------------------c
c     clear integer array ia(nx,ny)                                      c
c------------------------------------------------------------------------c

      dimension ia(nx,ny)

      do 10 i=1,nx
      do 10 j=1,ny
        ia(i,j) = 0.0
   10 continue

      return
      end      

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine zero2c(nx,ny,ca)
c------------------------------------------------------------------------c
c     clear complex array ca(nx,ny)                                      c
c------------------------------------------------------------------------c

      complex ca(nx,ny)

      do 10 i=1,nx
      do 10 j=1,ny
        ca(i,j) = (0.0,0.0)
   10 continue

      return
      end      
      
c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine zerof(n,a)                                          
c------------------------------------------------------------------------c
c     clear array a(n)                                                   c
c------------------------------------------------------------------------c

      dimension a(n)

      do 10 i = 1, n
        a(i) = 0.0
   10 continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      subroutine norm(nc,icos,nbin,xmin,st)
c------------------------------------------------------------------------c
c     normalize distribution of angles to uniform distribution.          c
c------------------------------------------------------------------------c
c     nc      dimension of group                                         c
c     icos(n) histogram array                                            c
c     nbin    number of bins                                             c
c     xmin    minimum value in histogram                                 c
c     st      bin width                                                  c
c------------------------------------------------------------------------c
      dimension icos(nbin)
      real c(20)

      npar = 8
      c(1) = 0.99994
      c(2) =-0.54220
      c(3) =-0.34209
      c(4) = 0.01935
      c(5) =-0.13441
      c(6) = 0.02259
      c(7) =-0.05701
      c(8) = 0.00732

      x0 = xmin+0.5*st
      f0 = sqrt((1.0-x0)/x0)
      nt = 0
      np = 0

      do 5  i=1,nbin
         nt = nt + icos(i)
   5  continue

      do 10 i=1,nbin
         x = xmin+(i-0.5)*st
         x = min(x,1.0-st)
         f = sqrt((1.0-x)/x)
         if(nc .eq. 3)then

         fit = 0.0
         do 20 j=1,npar
            l  = j-1
            y  = 2*x-1
            fit= fit + c(j)*p(l,y)
  20     continue
         f = f*fit

         endif
         icos(i) = icos(i)*f0/f
         np= np+icos(i)
  10  continue

      if(np .eq. 0) np = 1

      r = float(nt)/float(np)
      do 15 i=1,nbin
         icos(i) = icos(i)*r
  15  continue

      return
      end

c----------------------------------------------------------------------+--
c----------------------------------------------------------------------+--

      function p(l,x)
c------------------------------------------------------------------------c
c     legendre polynomials                                               c
c------------------------------------------------------------------------c

      real pp(0:100)

      pp(0) = 1.0
      pp(1) =  x

      do 10 i=2,l
         pp(i) = ((2*i-1)*x*pp(i-1)-(i-1)*pp(i-2))/float(I)
 10   continue

      p = pp(l)

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine taumat
c-------------------------------------------------------------------------c
c     initilize four-dim. tau matrices, the first number is index         c
c-------------------------------------------------------------------------c
      complex tau
      common /taum/tau(2,4,2,2)

      do 1 ind=1,2
         index=3-2*ind
         do 2 mu=1,4
         do 2 k1=1,2
         do 2 k2=1,2
            tau(ind,mu,k1,k2)=(0.,0.)
  2      continue
         tau(ind,4,1,1)=(0.,-1.)*index
         tau(ind,4,2,2)=(0.,-1.)*index
         tau(ind,3,1,1)=(1.,0.)
         tau(ind,3,2,2)=(-1.,0.)
         tau(ind,1,1,2)=(1.,0.)
         tau(ind,1,2,1)=(1.,0.)
         tau(ind,2,1,2)=(0.,-1.)
         tau(ind,2,2,1)=(0.,1.)
1     continue
      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------
