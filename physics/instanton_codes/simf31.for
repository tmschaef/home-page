
      program sim
c-------------------------------------------------------------------------c
c     determine partition function of an interacting instanton ensemble   c
c     with streamline interaction.                                        c
c-------------------------------------------------------------------------c
c     version:             3.1                                            c
c     creation date:     01-21-95                                         c
c     last modification: 02-04-96                                         c
c-------------------------------------------------------------------------c
c     based on version 3.1 of isim.for                                    c
c-------------------------------------------------------------------------c
c     version 3.0: new version based on main program from tsim. IMSL sub- c
c     routines EVLHF, LFTHF replaced by LAPACK counterpart CHEEV. In or-  c
c     der to convert back to IMSL uncomment RNSET in main, change to RNUNFc
c     in rang, change to EVLHF/LFTHF/LFDHF in logdet and spect.           c
c-------------------------------------------------------------------------c
c     version 3.1: subroutine overl includes correct ratio ansatz. Also:  c
c     included size renormalisation (optional), redefined uc4r (no effect)c
c     and changed definition of uin (includes usame).                     c
c-------------------------------------------------------------------------c
c     version 1.0: first version of simf. introduced parameter alp to     c
c     use integration in coupling constant.                               c
c-------------------------------------------------------------------------c
c     version 2.0: various minor bugs corrected, new version of random    c
c     number generator ran2, additional equilibration steps before integ  c
c     in coupling constant, use variational estimate of initial rho.      c
c-------------------------------------------------------------------------c
c     version 2.1: include error estimates for partition function, quark  c
c     condensate and average action.                                      c
c-------------------------------------------------------------------------c
c     version 2.2: improved error estimate, add statitical and discreti-  c
c     sation errors in free energy.                                       c
c-------------------------------------------------------------------------c
c     version 2.3: use simple interpolation to estimate various obser-    c
c     vables at optimum density.                                          c
c-------------------------------------------------------------------------c
c     version 2.4: hard core in ratio ansatz interaction, controlled by   c
c     parameter fcut.                                                     c
c-------------------------------------------------------------------------c
c     version 3.0: adapted from new version 2.0 of tsimf. uses hysteresis c
c     integration to minimize error assocaited with slow equilibration.   c
c-------------------------------------------------------------------------c
c     version 3.1: minor changes in order to study nonperturbative beta   c
c     functions, multi-instanton interactions etc.                        c
c-------------------------------------------------------------------------c
c     input file: insim.dat                                               c
c     nc,nf   number of colors,flavors                                    c
c     nin     number of instantons                                        c
c     rh0     average size                                                c
c     sg      metropolis step in group space                              c
c     dz(4)   metropolis step in position                                 c
c     drh     metropolis step in size                                     c
c     ninit   number of iterations at full interaction                    c
c     niter   number of iterations at each alpha                          c
c     nalp    number of couplings alpha                                   c
c     aa      parameter in Shuryaks overlap                               c
c     rmu,rms fermion masses                                              c
c     ieq,kp1 number of ieration to equilibrate, print                    c
c     iread   input from file (iread=1), random start (iread=0)           c
c     steig   bin size for eigenvalue plot                                c
c     strho   same for size distribution                                  c
c     acut    not used                                                    c
c     fcut    hard core in streamline interaction                         c
c     scut    cutoff in beta function                                     c
c     vcut    not used                                                    c
c     dens0   minimum density                                             c
c     dens1   maximum density                                             c
c     ndens   number of different densities                               c
c-------------------------------------------------------------------------c
c     input file: insimf.dat                                              c
c     determines start configuration if iread=1 is used. may use old      c
c     configuration as stored in outfile.                                 c
c-------------------------------------------------------------------------c
c     output file: outsimf.dat                                            c
c     contains full statistics as well as simple histograms.              c
c-------------------------------------------------------------------------c
c     Note: Program works correctly only for nc=2 and nc=3. Also: nf=0 is c
c     treated correctly, for nf=1 m=rms is used (!). In general, have     c
c     (nf-1) quarks of mass rmu and one quark with mass rms.              c
c-------------------------------------------------------------------------c
c     note: everything in units of lambda_QCD^(-1)                        c
c-------------------------------------------------------------------------c

      parameter(n=3, ni=64, ni2=ni/2, nbin=150, nstep=100)

      dimension rlp(ni2,ni2), rcos(ni2,ni2), rcii(ni2,ni2)
      dimension rnii(ni2, ni2), rnia(ni2, ni2), mol(ni2)
      dimension zr(ni,5),e1(ni,6),e2(ni,6), wr(ni2)
      dimension irho(nbin), ieig(nbin), ilap(nbin), icos(nbin)
      dimension irii(nbin), iria(nbin), icii(nbin), iact(nbin)
      dimension mcos(nbin), mcii(nbin), mlap(nbin)

c------------------------------------------------------------------------c
c     arrays for results at different densities                          c
c------------------------------------------------------------------------c

      dimension f0(nstep), f1(nstep), xqq(nstep), sav(nstep)
      dimension xrho(nstep), xd(nstep)
      dimension df1(nstep), df1p(nstep), dxqq(nstep), dsav(nstep)
      dimension dxrho(nstep)

c------------------------------------------------------------------------c
c     arrays for loop over couplings                                     c
c------------------------------------------------------------------------c

      dimension dfh(2,0:nstep), ddfh(2,0:nstep)
      dimension df(0:nstep), ddf(0:nstep)

      common /alpha/alp,rhob,xnu
      common /param/al(4),aa,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /cden/ b, alc, p1, p2
      common /metr/ itot, ifail1, ifail2, ifail3, actt
      common /acti/ acold
      common /pi/   pi, eps
      common /sij/  sij(ni,ni)
      common /seed/ iseed
      common /cut/  acut,fcut,scut,vcut
      common /c1c2/ c1,c2
      common /acc/  expmax

c------------------------------------------------------------------------c
c     read input file                                                    c
c------------------------------------------------------------------------c

      open(unit=1,file='insimf.dat',status='unknown')
      open(unit=2,file='outsimf.dat',status='unknown')
c     open(unit=14,file='outfile.dat',status='unknown')

      read(1,*) nc,nf,nin,rh0
      read(1,*) sg,(dz(k),k=1,4),drh
      read(1,*) ninit,niter,nalp
      read(1,*) aa,rmu,rms
      read(1,*) ieq,kp1
      read(1,*) iread, steig, strho
      read(1,*) acut,fcut,scut,vcut
      read(1,*) dens0,dens1,ndens
      close(unit=1)

c------------------------------------------------------------------------c
c     echo input                                                         c
c------------------------------------------------------------------------c

      nalp   = 2*((nalp+1)/2)
      nalp   = max(nalp,2)
      nlast  = (niter/kp1)*kp1
      nwrite = (nlast-ieq)/kp1+1

c     write(14,204) nc, nf, nin, al, rh0
c     write(14,206) sg, (dz(k),k=1,4), drh, aa
c     write(14,207) rmu,rms,nwrite

      write(2,801) nc,nf,nin
      write(2,802) ninit,niter,nalp
      write(2,803) iread,ieq,kp1
      write(2,804) rmu,rms
      write(2,805) rh0,aa
      write(2,806) steig,strho
      write(2,855) sg,drh
      write(2,856) (dz(k),k=1,4)
      write(2,807) acut,fcut,scut,vcut
      write(2,808) dens0,dens1,ndens

 801  format(1x,' N_c  = ',i5,' N_f   = ',i5,' N_in = ', i5)
 802  format(1x,' Ninit= ',i5,' Niter = ',i5,' Nalpha = ',i5)
 803  format(1x,' Iread= ',i5,' Nieq  = ',i5,' Nprint= ',i5)
 804  format(1x,' m_u  = ',f10.4,' m_s  = ',f10.4)
 805  format(1x,' Rh0  = ',f10.4,' aa   = ',f10.4)
 806  format(1x,' steig= ',f10.4,' strh = ',f10.4)
 855  format(1x,' sg   = ',f10.4,' drh  = ',f10.4)
 856  format(1x,' dz1  = ',f10.4,' dz2  = ',f10.4,
     2          ' dz3  = ',f10.4,' dz4  = ',f10.4)
 807  format(1x,' acut = ',f10.4,' fcut = ',f10.4,
     2          ' scut = ',f10.4,' vcut = ',f10.4)
 808  format(1x,' dens0= ',f10.4,' den1 = ',f10.4,
     2          ' ndens= ',i5)

c------------------------------------------------------------------------c
c     parameters for beta function, streamline, etc.                     c
c------------------------------------------------------------------------c

      nd   = ni
      pi   = 4*atan(1.0)
      eps  = exp(-40*alog(2.0))
      c1   = 3.0*pi/8.0
      c2   = (3.0*pi/32.0)**1.333333333
      b    = 11.0/3.0*nc -2.0/3.0*nf
      bp   = 34.0/3.0*nc*nc - 13.0/3.0*nc*nf +nf/float(nc)
      xnu  = (b-4.0)/2.0
      p1   = 2*nc-bp/2/b
      p2   = bp/2/b
      cnc  = 1.34**nf *4.66*exp(-1.68*nc)/pi/pi/(nc-1)*(b/2)**(bp/2/b)
      alc  = alog(cnc)
      gam2 = 27.0/4.0*nc/float(nc**2-1)*pi**2
      beta = 10.0
      expmax=50.0

c     iseed =-476
c     iseed =-9234
      iseed =-56789

      write(2,809)  b,bp,p1,p2,cnc,iseed
 809  format(1x,' b    = ',f10.4,' bp   = ',f10.4,/,
     2       1x,' p1   = ',f10.4,' p2   = ',f10.4,/,
     3       1x,' c_Nc = ',f10.4,/,
     4       1x,' iseed= ',i7,/)

c     call rnset(iseed)
      dum = rang( )

c------------------------------------------------------------------------c
c     loop over densities                                                c
c------------------------------------------------------------------------c

      ddens = (dens1-dens0)/float(ndens-1)
      xld0  = log(dens0)
      xld1  = log(dens1)
      dxld  = (xld1-xld0)/float(ndens-1)

      do 1000 id=1,ndens

c     dens  = dens0+(id-1)*ddens
      dens  = exp(xld0+(id-1)*dxld)
      xd(id)= dens
      vol   = nin/dens
      all   = vol**0.25
      do 1010 i=1,4
         al(i)= all
 1010 continue

      write(2,*)   ' Simulation at fixed density '
      write(2,999)
      write(2,*)
      write(2,810) dens
      write(2,811) al(1),al(4)
      write(2,*)
  810 format(1x,'density  nin/vol = ',f10.4)
  811 format(1x,'box      V_3*tau = ',f10.4,'^3 * ',f10.4)

c------------------------------------------------------------------------c
c     estimate initial size                                              c
c------------------------------------------------------------------------c

      rh0 = (xnu/beta/gam2/dens)**0.25
      rh0 = min(rh0,0.75)

c------------------------------------------------------------------------c
c     step1: simulation at full interaction strength                     c
c------------------------------------------------------------------------c

      alp = 1.0
      rhob= rh0

      write(2,*)   ' Simulation at full interaction strength'
      write(2,999)
      write(2,*)

c------------------------------------------------------------------------c
c     clear arrays etc.                                                  c
c------------------------------------------------------------------------c

      call zero(nbin,ieig)
      call zero(nbin,irho)
      call zero(nbin,irii)
      call zero(nbin,iria)
      call zero(nbin,icii)
      call zero(nbin,ilap)

      itott = 0
      stott = 0.0
      rmdett= 0.0
      srhot = 0.0
      uintt = 0.0
      sdett = 0.0
      qqt   = 0.0
      sintt = 0.0
      rho2t = 0.0

      ifailt1 = 0
      ifailt2 = 0
      ifailt3 = 0
      itel = 0
      nconf= 0

c------------------------------------------------------------------------c
c     initial configuration is made or read                              c
c------------------------------------------------------------------------c

      call setup(nc,nin,nd,zr,e1,e2, iread)

      do 5 i = 1, nin
      do 5 j = 1, nin
         sij(j,i) = abs(sign(1.0,j-i+0.5)+sign(1.0,j-i-0.5))/2
    5 continue

c------------------------------------------------------------------------c
c     calculate action and fermionic determinant for first configuration c
c------------------------------------------------------------------------c

      call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,act,rv,uint,rlp,
     2                 rcos,rcii,rnii,rnia,mol)
      acold= act
      stot =-act
      if(nf .ne. 0) then
         sdet =-nf*nin*log(rmdet)
      else
         sdet = 0.0
      endif
      srho =-rv
      sint = sdet+uint

      qq = 0.0
      do 7 k=1,nin/2
         qq = qq + 2.0*rmu/(rmu**2+wr(k)**2)
   7  continue
      qq = qq/vol

      rho2 = 0.0
      do 8 k=1,nin
         rho2 = rho2+zr(k,5)**2
  8   continue
      rho2 = rho2/nin

c------------------------------------------------------------------------c
c     write observables for initial configuration                        c
c------------------------------------------------------------------------c

      write(2,*)   ' Starting values '
      write(2,611) stot, sint
      write(2,612) srho, sdet, uint
      write(2,613) rmdet, qq, sqrt(rho2)
      write(2,*)
      write(2,104)

  611 format(1x,' S   = ',f10.4,' Sint = ',f10.4)
  612 format(1x,' rv  = ',f10.4,' Sdet = ',f10.4,' Sglu = ',f10.4)
  613 format(1x,' rmd = ',f10.4,' <qq> = ',f10.4,' rho  = ',f10.4)
  103 format(1x,i4,i5,'(',i3,',',i3,',',i3,')',8f8.3)
  104 format(1x,'   i itot(fl1,fl2,fl3)   S       w_I     S_int',
     2   '   S_qu    S_gl    rho     <qq>')

c------------------------------------------------------------------------c
c------------------------------------------------------------------------c
c     iterations start                                                   c
c------------------------------------------------------------------------c
c------------------------------------------------------------------------c

      do 10 i = 1, ninit

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
            sdett = 0.0
            srhot = 0.0
            uintt = 0.0
            stott = 0.0
            stott2= 0.0
            qqt   = 0.0
            qqt2  = 0.0
            sintt = 0.0
            sintt2= 0.0
            rho2t = 0.0
            rho2t2= 0.0
            itel  =  0
         end if
         itel = itel + 1

c------------------------------------------------------------------------c
c     calculate action and fermionic determinant                         c
c------------------------------------------------------------------------c

         call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,act,rv,uint,rlp,
     2                 rcos,rcii,rnii,rnia,mol)

c------------------------------------------------------------------------c
c     average n(rho), gauge interaction, ferm. det., and total action    c
c------------------------------------------------------------------------c

         srho  =-rv
         stot  =-act
         if(nf .ne. 0) then
            sdet =-nf*nin*log(rmdet)
         else
            sdet = 0.0
         endif
         sdett = sdett + sdet
         sint  = sdet  + uint
         sintt = sintt + sint
         sintt2= sintt2+ sint**2
         srhot = srhot + srho
         uintt = uintt + uint
         rmdett= rmdett+ rmdet
         stott = stott + stot
         stott2= stott2+ stot**2

c------------------------------------------------------------------------c
c     quark condensate, average size                                     c
c------------------------------------------------------------------------c

         qq = 0.0
         do 17 k=1,nin/2
            qq = qq + 2.0*rmu/(rmu**2+wr(k)**2)
  17     continue
         qq  = qq/vol
         qqt = qqt + qq
         qqt2= qqt2+ qq**2

         rho2 = 0.0
         do 18 k=1,nin
            rho2 = rho2 + zr(k,5)**2
  18     continue
         rho2  = rho2/nin
         rho2t = rho2t + rho2
         rho2t2= rho2t2+ rho2**2

c------------------------------------------------------------------------c
c     histogram action, size, eigenvalues, relative orientation          c
c------------------------------------------------------------------------c

         if (i .ge. ieq) then
            stlap = 0.030
            stcos = 0.025
            strii = 0.025
            sact = 0.10
            stc4 = 0.025
            call lens(stot/nin,-6.0,sact,40,iact)
            do 50 k = 1, nin
               call lens(zr(k,5),0.0,strho,40,irho)
   50       continue
            do 60 k = 1, nin/2
               call lens(wr(k), 0.0,steig, 60, ieig)
               call lens(rlp(mol(k),k),0.0,stlap,50,mlap)
               call lens(rcos(mol(k),k),0.0,stcos,42,mcos)
               call lens(rcii(mol(k),k),0.0,strii,42,mcii)
               do 60 l = 1, nin/2
                  call lens(rlp(k,l), 0.0, stlap,50, ilap)
                  call lens(rcos(k,l),-stcos, stcos,42,icos)
                  call lens(rcii(k,l),-strii, strii,42,icii)
c                 call lens(rnia(k,l),-strii, strii,42,iria)
c                 if (l .ne. k) then
c                    call lens(rnii(k,l),-strii, strii,42,irii)
c                 end if
   60       continue
         end if

c------------------------------------------------------------------------c
c     every kp1 sweeps do statistics, save configuration                 c
c------------------------------------------------------------------------c

         if ((i/kp1)*kp1 .eq. i) then

            write(6,700) id,0,i
            write(6,701) (ifail1+ifail2+ifail3)/(3.0*nin)
            write(6,702) stot/nin
            write(6,703) qq/vol
            write(6,704) sint/nin
            write(6,705) sqrt(rho2)
            write(6,*)

            write(2,103) i,itott,ifailt1,ifailt2,ifailt3,
     1           stott/(itel*nin),srhot/(nin*itel),sintt/(nin*itel),
     2           sdett/(itel*nin),uintt/(itel*nin),sqrt(rho2t/itel),
     3           qqt/itel
            ifailt1 = 0
            ifailt2 = 0
            ifailt3 = 0
            itott = 0

            if(i .ge. ieq) then

            nconf = nconf + 1
c           write(14,*) nconf
c           do 25 j = 1, nin
c              write(14,102) (zr(j,k),k=1,5)
c              write(14,102) (e1(j,k), k=1,6)
c              write(14,102) (e2(j,k), k=1,6)
c  25       continue

            endif

         end if

 700     format(1x,'iteration (id,ia,in)',3(1x,i3))
 701     format(1x,'failure rate     ',f10.4)
 702     format(1x,'average action   ',f10.4)
 703     format(1x,'quark condensate ',f10.4)
 704     format(1x,'interaction      ',f10.4)
 705     format(1x,'average size     ',f10.4)

c------------------------------------------------------------------------c
c     next iteration                                                     c
c------------------------------------------------------------------------c

   10 continue

c------------------------------------------------------------------------c
c     determine partition function for variational ansatz                c
c------------------------------------------------------------------------c

      sav0 = stott/(nin*itel)
      rhob = sqrt(rho2t/itel)
      call logz0(nin,vol,f0e,xmu0)
      f0(id) = f0e
      ftot   = f0e

      write(6,*)
      write(6,710) nin/vol
      write(6,711) xmu0
      write(6,712) f0e
      write(6,*)

      write(2,*)
      write(2,*)   'Variational estimate '
      write(2,999)
      write(2,713) nin/vol,rhob,f0(id),xmu0,sav0
      write(2,*)

 710  format(1x,'density          ',f10.4)
 711  format(1x,'mu0              ',f10.4)
 712  format(1x,'free energy, a=0 ',f10.4)
 713  format(1x,'nin/vol = ',f10.4,' rho = ',f10.4,' f0 = ',f10.4,/,
     1       1x,'mu0     = ',f10.4,' Sav = ',f10.4)

c------------------------------------------------------------------------c
c     averages at alpha = 1                                              c
c------------------------------------------------------------------------c

      xqq(id)  = qqt/itel
      dxqq(id) = disp(itel,qqt,qqt2)
      sav(id)  = stott/(itel*nin)
      dsav(id) = disp(itel*nin,stott,stott2)
      xrho(id) = sqrt(rho2t/itel)
      drho2    = disp(itel,rho2t,rho2t2)
      dxrho(id)= drho2/(2.0*xrho(id))

c------------------------------------------------------------------------c
c     also: calculate contribution to coupling constant integral         c
c------------------------------------------------------------------------c

      write(2,*)   'Integrate over coupling alpha'
      write(2,999)
      write(2,*)

      rho2 = rho2t/itel
      drho2= disp(itel,rho2t,rho2t2)
      sint = sintt/itel
      dsint= disp(itel,sintt,sintt2)

      dint = xnu*rho2/rhob**2-sint/nin
      ddint= xnu*drho2/rhob**2-dsint/nin

      dfh(1,nalp)  =-nin/vol * dint
      ddfh(1,nalp) =-nin/vol *ddint
      dfh(2,nalp)  = dfh(1,nalp)
      ddfh(2,nalp) = ddfh(1,nalp)

c------------------------------------------------------------------------c
c     alpha = 1 is a boundary point                                      c
c------------------------------------------------------------------------c

      dalp = 1.0/float(nalp)
      dint = dint/2.0
      fsum = f0(id)
      dfi  =-nin/vol*dint*dalp
      fsum = fsum + dfi

      write(2,*)
      write(2,715) alp
      write(2,716) sqrt(rho2)
      write(2,717) sint/nin
      write(2,718) dint
      write(2,719) f0(id),dfi,fsum

 715  format(1x,'coupling          alp = ',f10.4)
 716  format(1x,'average size      rho = ',f10.4)
 717  format(1x,'interaction / I   Sint= ',f10.4)
 718  format(1x,'integrand / I     Dlog= ',f10.4)
 719  format(1x,'free energy F0,DF,Ftot= ',3(f10.4,1x))

c------------------------------------------------------------------------c
c------------------------------------------------------------------------c
c     loop over coupling strengths                                       c
c------------------------------------------------------------------------c
c------------------------------------------------------------------------c

      do 2000 ia=nalp-1,0,-1
      alp =  ia*dalp

      write(2,*)
      write(2,714)  alp
      write(2,999)
      write(2,*)
 714  format(1x,'simulation at alpha = ',f10.4)

c------------------------------------------------------------------------c
c     clear arrays etc.                                                  c
c------------------------------------------------------------------------c

      call zero(nbin,ieig)
      call zero(nbin,irho)
      call zero(nbin,irii)
      call zero(nbin,iria)
      call zero(nbin,icii)
      call zero(nbin,ilap)

      itott = 0
      stott = 0.0
      rmdett= 0.0
      srhot = 0.0
      uintt = 0.0
      sdett = 0.0
      qqt   = 0.0
      sintt = 0.0
      rho2t = 0.0

      ifailt1 = 0
      ifailt2 = 0
      ifailt3 = 0
      itel = 0
      nconf= 0

c------------------------------------------------------------------------c
c     start from last configuration, calculate starting values           c
c------------------------------------------------------------------------c

      call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,act,rv,uint,rlp,
     2                 rcos,rcii,rnii,rnia,mol)
      acold= act
      stot =-act
      srho =-rv
      if(nf .ne. 0) then
         sdet =-nf*nin*log(rmdet)
      else
         sdet = 0.0
      endif
      sint = sdet+uint

c------------------------------------------------------------------------c
c     condensate, average size                                           c
c------------------------------------------------------------------------c

      qq = 0.0
      do 2007 k=1,nin/2
         qq = qq + 2.0*rmu/(rmu**2+wr(k)**2)
2007  continue
      qq = qq/vol

      rho2 = 0.0
      do 2008 k=1,nin
         rho2 = rho2+zr(k,5)**2
2008  continue
      rho2 = rho2/nin

c------------------------------------------------------------------------c
c     write observables for initial configuration                        c
c------------------------------------------------------------------------c

      write(2,*)   ' Starting values '
      write(2,611) stot, sint
      write(2,612) srho, sdet, uint
      write(2,613) rmdet, qq, sqrt(rho2)
      write(2,*)
      write(2,104)

c------------------------------------------------------------------------c
c     iterations start                                                   c
c------------------------------------------------------------------------c

      do 2010 i = 1, niter

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
            sdett = 0.0
            srhot = 0.0
            uintt = 0.0
            stott = 0.0
            qqt   = 0.0
            qqt2  = 0.0
            sintt = 0.0
            sintt2= 0.0
            rho2t = 0.0
            rho2t2= 0.0
            itel  =  0
         end if
         itel = itel + 1

c------------------------------------------------------------------------c
c     calculate action and fermionic determinant                         c
c------------------------------------------------------------------------c

         call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,act,rv,uint,rlp,
     2                 rcos,rcii,rnii,rnia,mol)

c------------------------------------------------------------------------c
c     average n(rho), gauge interaction, ferm. det., and total action    c
c------------------------------------------------------------------------c

         srho  =-rv
         stot  =-act
         if(nf .ne. 0) then
            sdet =-nf*nin*log(rmdet)
         else
            sdet = 0.0
         endif
         sdett = sdett + sdet
         sint  = sdet  + uint
         sintt = sintt + sint
         sintt2= sintt2+ sint**2
         srhot = srhot + srho
         uintt = uintt + uint
         rmdett= rmdett+ rmdet
         stott = stott + stot

c------------------------------------------------------------------------c
c     quark condensate, average size                                     c
c------------------------------------------------------------------------c

         qq = 0.0
         do 2017 k=1,nin/2
            qq = qq + 2.0*rmu/(rmu**2+wr(k)**2)
2017     continue
         qq  = qq/vol
         qqt = qqt + qq
         qqt2= qqt2+ qq**2

         rho2 = 0.0
         do 2018 k=1,nin
            rho2 = rho2 + zr(k,5)**2
2018     continue
         rho2  = rho2/nin
         rho2t = rho2t + rho2
         rho2t2= rho2t2+ rho2**2

c------------------------------------------------------------------------c
c     histogram action, size, eigenvalues, relative orientation          c
c------------------------------------------------------------------------c

         if (i .ge. ieq) then
            stlap = 0.030
            stcos = 0.025
            strii = 0.025
            sact = 0.10
            stc4 = 0.025
            call lens(stot/nin,-6.0,sact,40,iact)
            do 2050 k = 1, nin
               call lens(zr(k,5),0.0,strho,40,irho)
2050        continue
            do 2060 k = 1, nin/2
               call lens(wr(k), 0.0,steig, 60, ieig)
               call lens(rlp(mol(k),k),0.0,stlap,50,mlap)
               call lens(rcos(mol(k),k),0.0,stcos,42,mcos)
               call lens(rcii(mol(k),k),0.0,strii,42,mcii)
               do 2060 l = 1, nin/2
                  call lens(rlp(k,l), 0.0, stlap,50, ilap)
                  call lens(rcos(k,l),-stcos, stcos,42,icos)
                  call lens(rcii(k,l),-strii, strii,42,icii)
c                 call lens(rnia(k,l),-strii, strii,42,iria)
c                 if (l .ne. k) then
c                    call lens(rnii(k,l),-strii, strii,42,irii)
c                 end if
2060        continue
         end if

c------------------------------------------------------------------------c
c     every kp1 sweeps do statistics, save configuration                 c
c------------------------------------------------------------------------c

         if ((i/kp1)*kp1 .eq. i) then

            write(6,700) id,ia,i
            write(6,701) (ifail1+ifail2+ifail3)/(3.0*nin)
            write(6,702) stot/nin
            write(6,703) qq
            write(6,704) sint/nin
            write(6,705) sqrt(rho2)
            write(6,*)

            write(2,103) i,itott,ifailt1,ifailt2,ifailt3,
     1           stott/(itel*nin),srhot/(nin*itel),sintt/(nin*itel),
     2           sdett/(itel*nin),uintt/(itel*nin),sqrt(rho2t/itel),
     3           qqt/itel
            ifailt1 = 0
            ifailt2 = 0
            ifailt3 = 0
            itott = 0

            if(i .ge. ieq) then

            nconf = nconf + 1
c           write(14,*) nconf
c           do 25 j = 1, nin
c              write(14,102) (zr(j,k),k=1,5)
c              write(14,102) (e1(j,k), k=1,6)
c              write(14,102) (e2(j,k), k=1,6)
c  25       continue

            endif

         end if

c------------------------------------------------------------------------c
c     next iteration                                                     c
c------------------------------------------------------------------------c

2010  continue

c------------------------------------------------------------------------c
c     calculate contribution to partition function                       c
c------------------------------------------------------------------------c

      rho2 = rho2t/itel
      drho2= disp(itel,rho2t,rho2t2)
      sint = sintt/itel
      dsint= disp(itel,sintt,sintt2)

      dint = xnu*rho2/rhob**2-sint/nin
      ddint= xnu*drho2/rhob**2-dsint/nin

      dfh(1,ia)  =-nin/vol * dint
      ddfh(1,ia) =-nin/vol *ddint

c------------------------------------------------------------------------c
c     for ia.ne.0 average over up and down sweep                         c
c------------------------------------------------------------------------c

      dint = dint/2.0
      dfi  =-nin/vol*dint*dalp
      fsum = fsum + dfi

      write(2,*)
      write(2,715) alp
      write(2,716) sqrt(rho2)
      write(2,717) sint/nin
      write(2,718) dint
      write(2,719) f0(id),dfi,fsum

c------------------------------------------------------------------------c
c     next coupling strength alpha                                       c
c------------------------------------------------------------------------c

2000  continue

      dfh(2,0)  = dfh(1,0)
      ddfh(2,0) = ddfh(1,0)

c------------------------------------------------------------------------c
c------------------------------------------------------------------------c
c     loop over couplings in different direction                         c
c------------------------------------------------------------------------c
c------------------------------------------------------------------------c

      do 3000 ia=1,nalp
      alp =  ia*dalp

      write(2,*)
      write(2,714)  alp
      write(2,999)
      write(2,*)

c------------------------------------------------------------------------c
c     clear arrays etc.                                                  c
c------------------------------------------------------------------------c

      call zero(nbin,ieig)
      call zero(nbin,irho)
      call zero(nbin,irii)
      call zero(nbin,iria)
      call zero(nbin,icii)
      call zero(nbin,ilap)

      itott = 0
      stott = 0.0
      rmdett= 0.0
      srhot = 0.0
      uintt = 0.0
      sdett = 0.0
      qqt   = 0.0
      sintt = 0.0
      rho2t = 0.0

      ifailt1 = 0
      ifailt2 = 0
      ifailt3 = 0
      itel = 0
      nconf= 0

c------------------------------------------------------------------------c
c     start from last configuration, calculate starting values           c
c------------------------------------------------------------------------c

      call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,act,rv,uint,rlp,
     2                 rcos,rcii,rnii,rnia,mol)
      acold= act
      stot =-act
      srho =-rv
      if(nf .ne. 0) then
         sdet =-nf*nin*log(rmdet)
      else
         sdet = 0.0
      endif
      sint = sdet+uint

c------------------------------------------------------------------------c
c     condensate, average size                                           c
c------------------------------------------------------------------------c

      qq = 0.0
      do 3007 k=1,nin/2
         qq = qq + 2.0*rmu/(rmu**2+wr(k)**2)
3007  continue
      qq = qq/vol

      rho2 = 0.0
      do 3008 k=1,nin
         rho2 = rho2+zr(k,5)**2
3008  continue
      rho2 = rho2/nin

c------------------------------------------------------------------------c
c     write observables for initial configuration                        c
c------------------------------------------------------------------------c

      write(2,*)   ' Starting values '
      write(2,611) stot, sint
      write(2,612) srho, sdet, uint
      write(2,613) rmdet, qq, sqrt(rho2)
      write(2,*)
      write(2,104)

c------------------------------------------------------------------------c
c     iterations start                                                   c
c------------------------------------------------------------------------c

      do 3010 i = 1, niter

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
            sdett = 0.0
            srhot = 0.0
            uintt = 0.0
            stott = 0.0
            qqt   = 0.0
            qqt2  = 0.0
            sintt = 0.0
            sintt2= 0.0
            rho2t = 0.0
            rho2t2= 0.0
            itel  =  0
         end if
         itel = itel + 1

c------------------------------------------------------------------------c
c     calculate action and fermionic determinant                         c
c------------------------------------------------------------------------c

         call spect(nc,nin,nd,zr,e1,e2,rmdet,wr,act,rv,uint,rlp,
     2                 rcos,rcii,rnii,rnia,mol)

c------------------------------------------------------------------------c
c     average n(rho), gauge interaction, ferm. det., and total action    c
c------------------------------------------------------------------------c

         srho  =-rv
         stot  =-act
         if(nf .ne. 0) then
            sdet =-nf*nin*log(rmdet)
         else
            sdet = 0.0
         endif
         sdett = sdett + sdet
         sint  = sdet  + uint
         sintt = sintt + sint
         sintt2= sintt2+ sint**2
         srhot = srhot + srho
         uintt = uintt + uint
         rmdett= rmdett+ rmdet
         stott = stott + stot

c------------------------------------------------------------------------c
c     quark condensate, average size                                     c
c------------------------------------------------------------------------c

         qq = 0.0
         do 3017 k=1,nin/2
            qq = qq + 2.0*rmu/(rmu**2+wr(k)**2)
3017     continue
         qq  = qq/vol
         qqt = qqt + qq
         qqt2= qqt2+ qq**2

         rho2 = 0.0
         do 3018 k=1,nin
            rho2 = rho2 + zr(k,5)**2
3018     continue
         rho2  = rho2/nin
         rho2t = rho2t + rho2
         rho2t2= rho2t2+ rho2**2

c------------------------------------------------------------------------c
c     histogram action, size, eigenvalues, relative orientation          c
c------------------------------------------------------------------------c

         if (i .ge. ieq) then
            stlap = 0.030
            stcos = 0.025
            strii = 0.025
            sact = 0.10
            stc4 = 0.025
            call lens(stot/nin,-6.0,sact,40,iact)
            do 3050 k = 1, nin
               call lens(zr(k,5),0.0,strho,40,irho)
3050        continue
            do 3060 k = 1, nin/2
               call lens(wr(k), 0.0,steig, 60, ieig)
               call lens(rlp(mol(k),k),0.0,stlap,50,mlap)
               call lens(rcos(mol(k),k),0.0,stcos,42,mcos)
               call lens(rcii(mol(k),k),0.0,strii,42,mcii)
               do 3060 l = 1, nin/2
                  call lens(rlp(k,l), 0.0, stlap,50, ilap)
                  call lens(rcos(k,l),-stcos, stcos,42,icos)
                  call lens(rcii(k,l),-strii, strii,42,icii)
c                 call lens(rnia(k,l),-strii, strii,42,iria)
c                 if (l .ne. k) then
c                    call lens(rnii(k,l),-strii, strii,42,irii)
c                 end if
3060        continue
         end if

c------------------------------------------------------------------------c
c     every kp1 sweeps do statistics, save configuration                 c
c------------------------------------------------------------------------c

         if ((i/kp1)*kp1 .eq. i) then

            write(6,700) id,ia,i
            write(6,701) (ifail1+ifail2+ifail3)/(3.0*nin)
            write(6,702) stot/nin
            write(6,703) qq
            write(6,704) sint/nin
            write(6,705) sqrt(rho2)
            write(6,*)

            write(2,103) i,itott,ifailt1,ifailt2,ifailt3,
     1           stott/(itel*nin),srhot/(nin*itel),sintt/(nin*itel),
     2           sdett/(itel*nin),uintt/(itel*nin),sqrt(rho2t/itel),
     3           qqt/itel
            ifailt1 = 0
            ifailt2 = 0
            ifailt3 = 0
            itott = 0

            if(i .ge. ieq) then

            nconf = nconf + 1
c           write(14,*) nconf
c           do 25 j = 1, nin
c              write(14,102) (zr(j,k),k=1,5)
c              write(14,102) (e1(j,k), k=1,6)
c              write(14,102) (e2(j,k), k=1,6)
c  25       continue

            endif

         end if

c------------------------------------------------------------------------c
c     next iteration                                                     c
c------------------------------------------------------------------------c

3010  continue

c------------------------------------------------------------------------c
c     calculate contribution to partition function                       c
c------------------------------------------------------------------------c

      rho2 = rho2t/itel
      drho2= disp(itel,rho2t,rho2t2)
      sint = sintt/itel
      dsint= disp(itel,sintt,sintt2)

      dint = xnu*rho2/rhob**2-sint/nin
      ddint= xnu*drho2/rhob**2-dsint/nin

c------------------------------------------------------------------------c
c     for ia.ne.nalp average over up and down sweep                      c
c------------------------------------------------------------------------c

      if (ia .ne. nalp) then

      dfh(2,ia)  =-nin/vol * dint
      ddfh(2,ia) =-nin/vol *ddint

      dint = dint/2.0
      dfi  =-nin/vol*dint*dalp
      fsum = fsum + dfi

      write(2,*)
      write(2,715) alp
      write(2,716) sqrt(rho2)
      write(2,717) sint/nin
      write(2,718) dint
      write(2,719) f0(id),dfi,fsum

      endif

c------------------------------------------------------------------------c
c     next coupling strength alpha                                       c
c------------------------------------------------------------------------c

3000  continue

c------------------------------------------------------------------------c
c------------------------------------------------------------------------c
c     integration in coupling finished, average over up and down sweeps  c
c------------------------------------------------------------------------c
c------------------------------------------------------------------------c

      do 4005 ia=0,nalp
         df(ia) = ( dfh(1,ia) + dfh(2,ia) )/2.0
         ddf(ia)= ( ddfh(1,ia)+ ddfh(2,ia))/2.0
4005  continue

c------------------------------------------------------------------------c
c     integrate with boundaries; add error in squres                     c
c------------------------------------------------------------------------c

      fn     = ( df(0)+ df(nalp))/2.0*dalp
      dftot2 = (ddf(0)/2.0*dalp)**2 + (ddf(nalp)/2.0*dalp)**2

      do 4010 ia=1,nalp-1
         fn     = fn     +  df(ia)*dalp
         dftot2 = dftot2 +(ddf(ia)*dalp)**2
 4010 continue

c------------------------------------------------------------------------c
c     estimate discretization error from F(nalp/2)                       c
c------------------------------------------------------------------------c

      dalph = dalp*2.0
      fnh   =(df(0)+df(nalp))/2.0*dalph

      do 4020 ia=2,nalp-2,2
         fnh = fnh + df(ia)*dalph
 4020 continue

c------------------------------------------------------------------------c
c     estimate error from comparing up and down sweeps                   c
c------------------------------------------------------------------------c

      fnh1 =  ( df(0)+ df(nalp))/2.0*dalp

      do 4030 ia=1,nalp-1
         fnh1 = fnh1 +  dfh(1,ia)*dalp
 4030 continue

      fnh2 =  ( df(0)+ df(nalp))/2.0*dalp

      do 4040 ia=1,nalp-1
         fnh2 = fnh2 +  dfh(2,ia)*dalp
 4040 continue

c------------------------------------------------------------------------c
c     total free energy and error estimate                               c
c------------------------------------------------------------------------c

      f1(id)  = f0(id) + fn
      df1d    = abs(fn-fnh)/3.0
      df1s    = sqrt(dftot2)
      df1(id) = sqrt(df1d**2+df1s**2)
      df1p(id)= abs(fnh2-fnh1)/2.0

      write(2,*)
      write(2,*)   'integration in coupling finished'
      write(2,999)
      write(2,*)
      write(2,720) f0(id)
      write(2,721) f1(id)
      write(2,722) xqq(id)
      write(2,723) sav(id)
      write(2,724) xrho(id)

 720  format(1x,'free energy       F0  = ',f10.4)
 721  format(1x,'free energy       F1  = ',f10.4)
 722  format(1x,'quark condensate  <qq>= ',f10.4)
 723  format(1x,'average action    S_av= ',f10.4)
 724  format(1x,'average size      rho = ',f10.4)

c------------------------------------------------------------------------c
c     final output for full interaction                                  c
c------------------------------------------------------------------------c

c     open(unit=5,file='outplot.dat',status='unknown')
c     write(5,204)  nc, nf, nin, al, rh0
c     write(5,206)  sg, (dz(k),k=1,4), drh, aa
c     write(5,205)  rmu,rms

  204 format (1x,3i5,5f10.4)
  206 format (1x,7f10.4)
  205 format (1x,2f10.4)
  207 format (1x,2f10.4,1x,i5)
  102 format(6f12.5)

c     write (5,901) (irho(k), k = 1, 40)
c     write (5,901) (ieig(k), k = 1, 60)
c     write (5,901) (ilap(k), k = 1, 50)
c     write (5,901) (icos(k), k = 1, 42)
c     write (5,901) (icii(k), k = 1, 42)
c     write (5,901) (iria(k), k = 1, 42)
c     write (5,901) (irii(k), k = 1, 42)

  901 format(10(10i6,/))

c------------------------------------------------------------------------c
c     plot histograms for size, eigenvalues, overlaps, orientation, etc  c
c------------------------------------------------------------------------c

      write (2,113)
  113 format(1x,//,1x, '   Distribution of the sizes ')
      call lev(0.0,strho,40,2,irho)
      write (2,112)
  112 format(1x,//,1x, ' Distribution of the eigenvalues ')
      call lev(0.0,steig, 60,2, ieig)
c     write (2,117)
c 117 format(1x, ' Distribution of the overlap matrix elements ')
c     call lev(0.0, stlap,50,2,ilap)
c     write (2,117)
c 117 format(1x, ' Distribution of max overlap matrix elements ')
c     call lev(0.0, stlap,50,2,mlap)
c     write (2,118)
c 118 format(1x, ' Distribution of cos^2(th) ')
c     call lev(-stcos,stcos,42,2,icos)
c     write (2,119)
c 119 format(1x, ' Distribution of cos^2(alp) ')
c     call lev(-strii,strii,42,2,icii)

      call norm(icos,42,0.0,stcos)
      call norm(icii,42,0.0,strii)
c     write (2,124)
c 124 format(1x, ' normalized Distribution of cos^2(th) ')
c     call lev(-stcos,stcos,42,2,icos)
c     write (2,125)
c 125 format(1x, ' normalized Distribution of cos^2(alp) ')
c     call lev(-strii,strii,42,2,icii)

c     write(2,120)
c 120 format(1x, ' Distribution of cos^2(th) in pairs ')
c     call lev(-stcos,stcos,42,2,mcos)
c     write (2,121)
c 121 format(1x, ' Distribution of cos^2(alp) in pairs ')
c     call lev(-strii,strii,42,2,mcii)

      call norm(mcos,42,0.0,stcos)
      call norm(mcii,42,0.0,strii)
c     write(2,122)
c 122 format(1x, ' normalized distribution of cos^2(th) in molecules  ')
c     call lev(-stcos,stcos,42,2,mcos)
c     write (2,123)
c 123 format(1x, ' normalized distribution of cos^2(alp) in molecules ')
c     call lev(-strii,strii,42,2,mcii)

c      write (2,219)
c  219 format(1x, ' Distribution of the i-a norm  ')
c      call lev(-strii,strii,42,2,irii)
c      write (2,218)
c  218 format(1x, ' Distribution of the i-i norm  ')
c      call lev(-strii,strii,42,2,iria)
c      write (2,318)
c  318 format(1x, ' Distribution of the total action ')
c      call lev(-6.0, sact,40,2,iact)

c------------------------------------------------------------------------c
c     next density                                                       c
c------------------------------------------------------------------------c

1000  continue

c------------------------------------------------------------------------c
c------------------------------------------------------------------------c
c     loop over densities finished, final output                         c
c------------------------------------------------------------------------c
c------------------------------------------------------------------------c

      write(2,*)
      write(2,*)
      write(2,753)
      write(2,998)
      do 905 k=1,ndens
c        write(2,751) xd(k),f0(k),f1(k),xqq(k),sav(k)
         write(2,754) xd(k),f0(k),f1(k),df1(k),df1p(k),xqq(k),dxqq(k)
  905 continue

      write(2,*)
      write(2,755)
      write(2,999)
      do 910 k=1,ndens
         write(2,756) xd(k),f1(k),sav(k),dsav(k),xrho(k),dxrho(k)
  910 continue

c------------------------------------------------------------------------c
c     estimate for equilibrium density                                   c
c------------------------------------------------------------------------c

      if(ndens .ge. 3) then

      xmin = xd(1)
      fmin = f1(1)
      kmin = 1

      do 920 k=2,ndens
         if (f1(k) .lt. fmin) then
            kmin = k
            xmin = xd(k)
            fmin = f1(k)
         endif
  920 continue

      if(kmin .eq. 1) then
         k1 = 1
         k2 = 2
         k3 = 3
      else if(kmin .eq. ndens) then
         k1 = ndens-2
         k2 = ndens-1
         k3 = ndens
      else
         k1 = kmin-1
         k2 = kmin
         k3 = kmin+1
      endif

c------------------------------------------------------------------------c
c     simple parabolic interpolation, determine optimum density          c
c------------------------------------------------------------------------c

      x1 = xd(k1)
      x2 = xd(k2)
      x3 = xd(k3)
      y1 = f1(k1)
      y2 = f1(k2)
      y3 = f1(k3)

      den = (x2-x1)*(x3-x1)*(x3-x2)
      p0 = (-(x2**2*x3*y1) + x2*x3**2*y1 + x1**2*x3*y2
     1      - x1*x3**2*y2 - x1**2*x2*y3 + x1*x2**2*y3)/den
      p1 = (-(x2**2*y1) + x3**2*y1 + x1**2*y2 - x3**2*y2
     1      - x1**2*y3 + x2**2*y3)/(-den)
      p2 = (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)/den


      xmin =-(-(x2**2*y1) + x3**2*y1 + x1**2*y2 - x3**2*y2
     1        - x1**2*y3 + x2**2*y3)/
     2       (2*(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3))
      fmin= p0+p1*xmin+p2*xmin**2

c------------------------------------------------------------------------c
c     other observables: quark condensate                                c
c------------------------------------------------------------------------c

      y1 = xqq(k1)
      y2 = xqq(k2)
      y3 = xqq(k3)

      p0 = (-(x2**2*x3*y1) + x2*x3**2*y1 + x1**2*x3*y2
     1      - x1*x3**2*y2 - x1**2*x2*y3 + x1*x2**2*y3)/den
      p1 = (-(x2**2*y1) + x3**2*y1 + x1**2*y2 - x3**2*y2
     1      - x1**2*y3 + x2**2*y3)/(-den)
      p2 = (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)/den

      qqmin= p0+p1*xmin+p2*xmin**2

c------------------------------------------------------------------------c
c     other observables: average action                                  c
c------------------------------------------------------------------------c

      y1 = sav(k1)
      y2 = sav(k2)
      y3 = sav(k3)

      p0 = (-(x2**2*x3*y1) + x2*x3**2*y1 + x1**2*x3*y2
     1      - x1*x3**2*y2 - x1**2*x2*y3 + x1*x2**2*y3)/den
      p1 = (-(x2**2*y1) + x3**2*y1 + x1**2*y2 - x3**2*y2
     1      - x1**2*y3 + x2**2*y3)/(-den)
      p2 = (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)/den

      savmin= p0+p1*xmin+p2*xmin**2

c------------------------------------------------------------------------c
c     other observables: average size                                    c
c------------------------------------------------------------------------c

      y1 = xrho(k1)
      y2 = xrho(k2)
      y3 = xrho(k3)

      p0 = (-(x2**2*x3*y1) + x2*x3**2*y1 + x1**2*x3*y2
     1      - x1*x3**2*y2 - x1**2*x2*y3 + x1*x2**2*y3)/den
      p1 = (-(x2**2*y1) + x3**2*y1 + x1**2*y2 - x3**2*y2
     1      - x1**2*y3 + x2**2*y3)/(-den)
      p2 = (-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)/den

      rhomin= p0+p1*xmin+p2*xmin**2

c------------------------------------------------------------------------c
c     output                                                             c
c------------------------------------------------------------------------c

      write(2,*)
      write(2,*)   'equilibrium density'
      write(2,999)
      write(2,752) xmin,fmin,qqmin
      write(2,757) savmin,rhomin

      endif

  750 format(1x,'   N/V       F0        F1        <qq>      S_av')
  751 format(1x,5(f10.4))
  752 format(1x,' N/V = ',f10.4,' F = ',f10.4,' <qq> = ',f10.4)
  753 format(4x,'N/V',6x,'F0',7x,'F1',5x,'(stat)',3x,'(equ)',3x,'<qq>')
  754 format(1x,3(f9.4),2(' %',f7.4),(f9.4,' %',f7.4))
  755 format(1x,3x,'N/V',6x,'F1',7x,'S_av',7x,'rho')
  756 format(1x,2(f9.4),2(f9.4,' %',f7.4))
  757 format(1x,' Sav = ',f10.4,' rh= ',f10.4)
  998 format(1x,70('-'))
  999 format(1x,50('-'))


      stop
      end


c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine logz0(nin,vol,f0e,xmu0)
c------------------------------------------------------------------------c
c     calculate variational estimate for partition function              c
c------------------------------------------------------------------------c
c     nin    number of instantons                                        c
c     vol    four-volume                                                 c
c     f0e    variational free energy                                     c
c     xmu0   integral of d(rho)                                          c
c------------------------------------------------------------------------c

      common /alpha/ alp,rhob,xnu
      common /cden/  b, alc, p1, p2

c------------------------------------------------------------------------c
c     calculate mu0 = \int d\rho d(\rho)                                 c
c------------------------------------------------------------------------c

      rhom = 0.999
      nrho = 100
      drho = rhom/float(nrho)
      xmu0 = 0.0

      do 10 i=1, nrho-1, 2
         rho1 = i*drho
         rho2 = (i+1)*drho
         d1   = dvar(rho1)
         d2   = dvar(rho2)
         xmu0 = xmu0 + 4.0*d1 + 2.0*d2
  10  continue
      xmu0 = xmu0 - d2
      xmu0 = drho*xmu0/3.0

c------------------------------------------------------------------------c
c     one loop estimate                                                  c
c------------------------------------------------------------------------c

      nc   = 3
      bet  =-b*log(rhob)
      cnc  = exp(alc)
      xmu0a= cnc*bet**(2.0*nc)*exp(gammln(xnu))*(xnu/rhob**2)**(-xnu)

c------------------------------------------------------------------------c
c     free energy                                                        c
c------------------------------------------------------------------------c

      xlog2p = 1.837877067
      xlogn2 = (nin+1)*log(nin/2.0) - nin + xlog2p
      f0e    =-1.0/vol*( nin*log(vol*xmu0) - xlogn2 )
      f0a    =-1.0/vol*( nin*log(vol*xmu0a)- xlogn2 )

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function dvar(rho)
c------------------------------------------------------------------------c
c     variational ansatz for one instanton distribution                  c
c------------------------------------------------------------------------c
c     make sure that this agrees with whats being used in overl !        c
c------------------------------------------------------------------------c

      common /alpha/ alp,rhob,xnu
      common /cden/  b, alc, p1, p2
      common /cut/   acut,fcut,scut,vcut

      bet  =-b*log(rho)
      betnp= sqrt(bet**2+scut**2) 
      xlogd= alc - 5.0*log(rho) - betnp 
     2      + (p1+p1*p2/bet)*log(bet) - xnu*(rho/rhob)**2
      dvar = exp(xlogd)

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function gammln(xx)
c------------------------------------------------------------------------c
c     gamma function, numerical recipes                                  c
c------------------------------------------------------------------------c
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     2    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/

      x   = xx-one
      tmp = x+fpf
      tmp =(x+half)*log(tmp)-tmp
      ser = one

      do 11 j=1,6
         x   = x + one
         ser = ser + cof(j)/x
 11   continue

      gammln = tmp+log(stp*ser)

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function disp(n,xtot,x2tot)
c------------------------------------------------------------------------c
c     determine rms fluctuation, disp=del x                              c
c------------------------------------------------------------------------c
c     n     number of events                                             c
c     xtot  sum of x_i                                                   c
c     x2tot sum of X_i^2                                                 c
c------------------------------------------------------------------------c

      xav  = xtot/n
      del2 = x2tot/n**2 - xav**2/n
      del2 = max(del2,0.0)
      disp = sqrt(del2)

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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
      common /pi/ pi,eps0

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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
      common /pi/ pi,eps0

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function rang()
c------------------------------------------------------------------------c
c     various random number generators                                   c
c------------------------------------------------------------------------c
      common /seed/ iseed
c     rang = rnunf()
c     rang = ranf()
      rang = ran2(iseed)
c     rang = ran(iseed)

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function ran2old(idum)
c---------------------------------------------------------------------------c
c     numerical recipes random number generator, set iseed to a negative    c
c     value to initialize sequence. on first call, sequence is initialized  c
c     automatically.                                                        c
c---------------------------------------------------------------------------c
c     note: unix machines require save statement !                          c
c---------------------------------------------------------------------------c
      parameter (m=714025,ia=1366,ic=150889,rm=1.4005112E-6)
      dimension ir(97)
      save ir,iy,iff
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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine overl(n,nin,nd,nd2,zr,e1,e2,clp,bfc,rhofc,uintt,il,iu,
     1                 iop)
c------------------------------------------------------------------------c
c     Calculate overlap matrix elements, gauge field action and instan-  c
c     ton measure. Can be used in two different modes: evaluate all ma-  c
c     trix elements (iop=0), or update all matrix elements affected by   c
c     changing the collective coordinates of a single instanton.         c
c------------------------------------------------------------------------c
c     streamline version, includes ratio ansatz as comments.             c
c------------------------------------------------------------------------c
c     input:  n            number of colors                              c
c             nin          number of instantons                          c
c             nd           max number of instantons                      c
c             nd2          nd/2                                          c
c             zr(i,5)      position, size of instanton i                 c
c             e1(i,6)      first row of orientation matrix of inst i     c
c             e2(i,6)      second row                                    c
c     output: clp(i,j)     overlap matrix element                        c
c             bfc          instanton measure                             c
c             rhofc        t'Hooft factor                                c
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
c     common /d/      store positions for update (iop=1,2,3)             c
c     common /ucr/    store orientations for update                      c
c     common /cut/    parameters for hard core                           c
c     common /c1c2/   parameters for streamline m.e. parametrization     c
c------------------------------------------------------------------------c

      parameter(n2=64,ni=n2/2)
      complex u1, u2, tr, cl, clp, ctr, uc, utr, ut4, uts
      complex uc1r, uc2r, uc3r, uc4r, ci, ur, uu4, usus4, u4s
      dimension zr(nd,5), U1(n2,3,3),u2(n2,3,3),e1(nd,6), e2(nd,6)
      dimension tr(n2,3,3), clp(nd2,nd2), cl(ni,ni)
      dimension z1(5), z2(5), uc(n2,2,2), ur(n2,2,2)
      dimension ra(n2), rai(n2), si(n2)
      dimension ffr(n2), smni(ni), smna(ni), reff2(n2)
      dimension du1(n2), du2(n2), uinta(n2)
      dimension drfc(n2), rfc(n2), zr2(n2), bet(n2), betnp(n2)

      common /alpha/ alp,rhob,xnu
      common /param/ alb(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /cden/  b,alc,p1,p2
      common /d/     uint(n2,n2), r2(n2,n2), dis(n2,n2,4), pcl(n2,n2)
      common /pi/    pi, eps
      common /sij/   sij(n2,n2)
      common /cut/   acut, fcut, scut, vcut
      common /ucr/   uc1r(n2,n2), uc2r(n2,n2), uc3r(n2,n2), uc4r(n2,n2)
      common /c1c2/  c1, c2

      ci = cmplx(0.0,1.0)
      nih = nin/2

c------------------------------------------------------------------------c
c     si = -/+ for instantons/antiinstantons                             c
c------------------------------------------------------------------------c

      do 5 j = 1, nin
         si(j) = sign(1.0,j-nih-0.5)
         zr2(j)= zr(j,5)*zr(j,5)
    5 continue

c------------------------------------------------------------------------c
c     reconstruct orientation U1, U2=U1^\dagger                          c
c------------------------------------------------------------------------c

      call su3(n,nin,nd,e1,e2,u1)

      do 30 i1 = 1, n
      do 30 i2 = 1, n
         do 31 i = 1, nin
            u2(i,i1,i2) = conjg(u1(i,i2,i1))
   31    continue
   30 continue

      rhofc = 0.0
      uintt = 0.0
      do 216 i = 1, nin
        rfc(i) = 1.0
  216 continue

c------------------------------------------------------------------------c
c     for full calculation or position update recalculate distances      c
c------------------------------------------------------------------------c

      if(iop.eq.0 .or. iop.eq.2) then

      do 121 i = il, iu
         do 21 j = 1, nin
            r2(j,i) = 0.0
  21     continue
  121 continue

      do 212 i = il, iu
         do 20 ir = 1, 4
         do 22 j = 1, nin
            ds = zr(i,ir) - zr(j,ir)
            ads = abs(ds)
            asg = 1+sign(1.0,ads-alb(ir)/2)
            dis(j,i,ir) = ds - alb(ir)/2*asg*sign(1.0,ds)
            r2(j,i) = r2(j,i) + dis(j,i,ir)**2
   22    continue
   20    continue
  212 continue

c------------------------------------------------------------------------c
c     also correct reverse order                                         c
c------------------------------------------------------------------------c

      do 215 i = il, iu
         do 214 ir = 1, 4
         do 214 j = 1, nin
            dis(i,j,ir) = -dis(j,i,ir)
  214    continue
         do 215 j = 1, nin
            r2(i,j) = r2(j,i)
  215 continue

      end if

c------------------------------------------------------------------------c
c     calculate orientation vectors (except for rho update)              c
c------------------------------------------------------------------------c

      if(iop.eq.0 .or. iop.eq.2 .or. iop.eq.3) then

      do 150 i = il, iu

c------------------------------------------------------------------------c
c     tr = i(z_J-z_I).\tau^(+\-) for (IJ)=(IA,II,AA)/(AI)                c
c------------------------------------------------------------------------c

         do 90 j = 1, nin
            sgp = 1 - (1+si(i))*(1-si(j))/2
            tr(j,1,1) = ci*cmplx(dis(j,i,3),-sgp*dis(j,i,4))
            tr(j,2,2) = ci*cmplx(-dis(j,i,3),-sgp*dis(j,i,4))
            tr(j,1,2) = ci*cmplx(dis(j,i,1),-dis(j,i,2))
            tr(j,2,1) = ci*cmplx(dis(j,i,1), dis(j,i,2))
   90    continue

c------------------------------------------------------------------------c
c     relative orientation matrix UC=U_I^(+)*U_J                         c
c------------------------------------------------------------------------c

         do 35 k = 1, 2
         do 35 l = 1, 2
            do 95 j = 1, nin
               uc(j,k,l) = cmplx(0.0, 0.0)
   95       continue

            do 35 m = 1, n
            do 97 j = 1, nin
               uc(j,k,l) = uc(j,k,l) + u2(i,k,m)*u1(j,m,l)
   97       continue
   35    continue

c------------------------------------------------------------------------c
c     UR = R_IA.\tau^(+)*U_I^(+)*U_A  /  U_A^(+)*U_J*R_IA.\tau^(-)       c
c------------------------------------------------------------------------c

         do 36 k = 1, 2
         do 36 l = 1, 2
            do 105 j = 1, nin
               ur(j,k,l) = cmplx(0.0,0.0)
  105       continue

            do 112 m = 1, 2
            do 110 j = 1, nin
               sgp = 1 - (1+si(i))*(1-si(j))/2
               ur(j,k,l) = ur(j,k,l) + tr(j,k,m)*uc(j,m,l)*(1+sgp)/2
     2                               + uc(j,k,m)*tr(j,m,l)*(1-sgp)/2
  110       continue
  112       continue
   36    continue

c------------------------------------------------------------------------c
c     ucir_JI = 1/(2i) tr(ur*\tau_i^(-)), Note: (+/-) does not matter    c
c------------------------------------------------------------------------c

         do 115 j = 1, nin
            sgp = 1 - (1+si(i))*(1-si(j))/2
            uc1r(j,i) =  (ur(j,1,2)+ur(j,2,1))/2/ci
            uc2r(j,i) =  (ur(j,1,2)-ur(j,2,1))/2
            uc3r(j,i) =  (ur(j,1,1)-ur(j,2,2))/2/ci
            uc4r(j,i) =  (ur(j,1,1)+ur(j,2,2))/2*sgp
  115    continue
  150 continue

c------------------------------------------------------------------------c
c     orientation update, correct reverse order                          c
c------------------------------------------------------------------------c

      do 152 i = il, iu
      do 151 j = 1, nin
         uc1r(i,j) = -conjg(uc1r(j,i))
         uc2r(i,j) = -conjg(uc2r(j,i))
         uc3r(i,j) = -conjg(uc3r(j,i))
         uc4r(i,j) = -conjg(uc4r(j,i))
  151 continue
  152 continue

      end if

c------------------------------------------------------------------------c
c     big loop : calculate overlaps and gauge interaction of (IJ)        c
c------------------------------------------------------------------------c

      do 12 i = il, iu
      ihal = nint((-si(i)+1)/2.0) * nih
      do 80 jlp = 1, nih

c------------------------------------------------------------------------c
c     note : in streamline ansatz there is no (II) and (AA) interaction  c
c------------------------------------------------------------------------c

         j = jlp + ihal

c------------------------------------------------------------------------c
c     conformal parameter lambda                                         c
c------------------------------------------------------------------------c

         acf  = (r2(j,i)+zr2(j)+zr2(i))/zr(j,5)/zr(i,5)
         disc = sqrt(acf*acf-4)
         rlam = (acf+disc)/2
         rl2  = rlam*rlam
         den  = (rl2-1)**3
         rn   = r2(j,i)/(zr(j,5)*zr(i,5))

c------------------------------------------------------------------------c
c     parametrization of fermionic overlap matrix element                c
c------------------------------------------------------------------------c

         pcl(j,i) = c1*rlam*sqrt(rlam)/(1+1.25*(rl2-1)
     2            +c2*(rl2-1)**2)**0.75

c------------------------------------------------------------------------c
c     ratio ansatz result                                                c
c------------------------------------------------------------------------c

c        pcl(j,i) = 4.0*sqrt(rn)/(rn+2.0)**2

c------------------------------------------------------------------------c
c     orientation invariants                                             c
c------------------------------------------------------------------------c

         r2i  = sij(j,i)/(r2(j,i)+eps)
         uus  = (cabs(uc1r(j,i))**2 + cabs(uc2r(j,i))**2
     2         + cabs(uc3r(j,i))**2)*r2i
         u4s  = cabs(uc4r(j,i))**2*r2i
         uus4 = uus + u4s
         uu4  = ((uc1r(j,i))**2 + (uc2r(j,i))**2
     2         + (uc3r(j,i))**2 + (uc4r(j,i))**2)*r2i
         usus4= cmplx(uu4)
         d    = uus - 3*u4s
         troto= real(uus4*uus4+ 2*uu4*usus4)

c------------------------------------------------------------------------c
c     parametrization of IA gauge interaction (units of S_0)             c
c------------------------------------------------------------------------c

         du2(j) = -4*d*(1-rl2*rl2+4*rl2*alog(rlam))/den
     2           + 2*(d*d+troto)*(1-rl2+(1+rl2)*alog(rlam))/den

c------------------------------------------------------------------------c
c     include hard core, parameter fcut                                  c
c------------------------------------------------------------------------c

         du2(j) = du2(j) + fcut*uus4/rlam**4

c------------------------------------------------------------------------c
c     yet another core                                                   c
c------------------------------------------------------------------------c

c        if( du2(j) .lt. -vcut) du2(j) = 100.0

c------------------------------------------------------------------------c
c     alternatively : ratio ansatz interaction                           c
c------------------------------------------------------------------------c

c        duin1 = 4.00/(rn+2.00)**2
c        duin2 =-1.66/(1.0+1.68*rn)**3-0.72*log(rn)/(1.0+0.42*rn)**4
c        duin3 = 2.73/(1.0+0.33*rn)**3-16.00/(rn+2.00)**2
c        du2(j)= (duin1+duin2)*uus4+duin3*u4s
c        du2(j) = du2(j) + fcut*uus4/rlam**4

c------------------------------------------------------------------------c
c     end of loop over instanton J                                       c
c------------------------------------------------------------------------c

 80   continue

c------------------------------------------------------------------------c
c     store interaction for given instanton I                            c
c------------------------------------------------------------------------c

      ihaii = nint((si(i)+1)/2.0) * nih
      do 106 jlp = 1, nih
         j = jlp + ihal
         uint(j,i) = du2(j)
  106 continue

c------------------------------------------------------------------------c
c     end of big loop over instanton I                                   c
c------------------------------------------------------------------------c

   12 continue

c------------------------------------------------------------------------c
c     include hard core for (II) and (AA) interaction                    c
c------------------------------------------------------------------------c

      do 1012 i = il, iu
      ihaii = nint((si(i)+1)/2.0) * nih
      do 1080 jlp = 1, nih
         j = jlp + ihaii

c------------------------------------------------------------------------c
c     conformal parameter and orientation invariants as above            c
c------------------------------------------------------------------------c

         acf  = (r2(j,i)+zr2(j)+zr2(i))/zr(j,5)/zr(i,5)
         disc = sqrt(abs(acf*acf-4))
         rlam = (acf+disc)/2
         r2i  = sij(j,i)/(r2(j,i)+eps)
         uus  = (cabs(uc1r(j,i))**2 + cabs(uc2r(j,i))**2
     2          +cabs(uc3r(j,i))**2)*r2i
         uus4 = uus + cabs(uc4r(j,i))**2 *r2i

c------------------------------------------------------------------------c
c     only hard core                                                     c
c------------------------------------------------------------------------c

         uint(j,i) = fcut*uus4/rlam**4

c------------------------------------------------------------------------c
c     alternatively : ratio ansatz interaction                           c
c------------------------------------------------------------------------c

c        ut4 = (uc1r(j,i)*dis(j,i,1)+uc2r(j,i)*dis(j,i,2)
c    1         +uc3r(j,i)*dis(j,i,3)+uc4r(j,i)*dis(j,i,4))*r2i
c        uts = uus4 - cabs(ut4)**2
c        rn = r2(j,i)/zr(i,5)/zr(j,5)+eps
c        duin1 = 0.63/(1.00+0.43*rn)**3
c        duin2 =-0.05*log(rn)/(1.00+1.17*rn)**4
c        duin3 = 0.071/(1.00+0.43*rn)**3
c        duin4 =-0.47*log(rn)/(1.00+1.17*rn)**4
c        uint(j,i) = (duin1+duin2)*uts + (duin3+duin4)*uts**2
c        uint(j,i) = uint(j,i) + fcut*uus4/rlam**4

1080  continue
1012  continue

c------------------------------------------------------------------------c
c     for update: correct reverse order                                  c
c------------------------------------------------------------------------c

      do 18 i = il, iu
      do 18 j = 1, nin
         uint(i,j) = uint(j,i)
         pcl(i,j) = pcl(j,i)
   18 continue

c------------------------------------------------------------------------c
c     finally: nih*nih matrix of (IA)-overlap matrix elements            c
c------------------------------------------------------------------------c

      do 312 i = 1, nih
      do 100 j = nih+1, nin
         clp(i,j-nih) = -ci*uc4r(j,i)*pcl(j,i)
     2                  /sqrt(r2(j,i)*zr(i,5)*zr(j,5))
 100  continue
 312  continue

c------------------------------------------------------------------------c
c     not used                                                           c
c------------------------------------------------------------------------c

      do 177  i = 1, nih
         smni(i) = 0.0
         smna(i) = 0.0
  177 continue

      do 178 i = 1, nih
         reff2(i)     = zr2(i)
         reff2(i+nih) = zr2(i+nih)
  178 continue

      do 179 i = 1, nih
         smni(i) = exp(smni(i))
         smna(i) = exp(smna(i))
  179 continue

c------------------------------------------------------------------------c
c     Edward's size renormalisation                                      c
c------------------------------------------------------------------------c

      do 410 i=1, nin
      do 420 j=1, nin
         ra(j)  = r2(j,i)/zr2(j)
         ffr(j) = sij(j,i)/(ra(j)+eps)
         rfc(i) = rfc(i) + ffr(j)
 420  continue
 410  continue

c------------------------------------------------------------------------c
c     log of density distribution, bet(i)=g^2/(8\pi)                     c
c------------------------------------------------------------------------c
c     use one/two loop beta function, nonperturbative corrections        c
c------------------------------------------------------------------------c

      do 15 i = 1, nin
         ro = zr2(i)
c        ro = zr2(i)/rfc(i)
         bet(i)  = -b/2*alog(ro)
         betnp(i)= sqrt(bet(i)*bet(i) + scut*scut)
         drfc(i) =(p1+p1*p2/bet(i))*alog(bet(i))
     1            -betnp(i)
     2            -5*alog(zr(i,5))+alc
     3            -xnu*(1.0-alp)*zr2(i)/rhob**2

   15 continue

c------------------------------------------------------------------------c
c     collect density distribution                                       c
c------------------------------------------------------------------------c

      rhofc = 0.0

      do 16 i = 1, nin
         rhofc = rhofc + drfc(i)
   16 continue

c------------------------------------------------------------------------c
c     (IA) gauge interaction                                             c
c------------------------------------------------------------------------c

      do 17 i = 1, nih
      uinta(i) = 0.0
      usum = 0.0
      uinta(i+nih) = 0.0
      usumi = 0.0
      do 19 j = nih+1, nin
         uinta(i)   = uinta(i) + (betnp(i)+betnp(j))/2 * uint(j,i)/2
         uinta(i+nih) = uinta(i+nih) + (betnp(i+nih)+betnp(j-nih))/2
     &                     * uint(j-nih,i+nih)/2
  19  continue
  17  continue

c------------------------------------------------------------------------c
c     (II) and (AA) gauge interaction (hard core only)                   c
c------------------------------------------------------------------------c

      usame = 0.0

      do 1017 i = 1, nih
      do 1017 j = 1, nih
         in = i + nih
         jn = j + nih
         usame = usame + (betnp(i)+betnp(j))/2 * uint(j,i)/2
         usame = usame + (betnp(in)+betnp(jn))/2 * uint(jn,in)/2
 1017 continue

c------------------------------------------------------------------------c
c     log(measure), t'Hooft plus gauge interaction                       c
c------------------------------------------------------------------------c

      uintt = 0.0
      do 188 i = 1, nin
         if(uinta(i) .lt. -betnp(i)) write(77,*) betnp(i),uinta(i)
         uintai= max(-betnp(i),uinta(i))
         uintt = uintt + uintai
 188  continue

      uintt= uintt + usame
      bfc  = rhofc - alp*uintt

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function fs(r2)
      common /param/ al(4), alpha,rh0,sg,dz(4),drh, nc, nf, rmu, rms

      fs = 2/(r2+alpha)**2

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

c---------------------------------------------------------------------+----
c---------------------------------------------------------------------+----

      subroutine setup(n,nin,nd,zr,e1,e2,iread)
c------------------------------------------------------------------------c
c     initialize instanton distribution                                  c
c------------------------------------------------------------------------c
c     n,nin,nd   nc, number of instantons, array dimension               c
c     zr,e1,e2   instanton configuration                                 c
c------------------------------------------------------------------------c
c     iread = 1  input from file                                         c
c     iread = 0  random configuration                                    c
c------------------------------------------------------------------------c

      complex u0
      dimension zr(nd,5), er1(6), er2(6),fal(4)
      dimension x(6), y(6), z(6), e1(nd,6), e2(nd,6)
      dimension fdz(4)
      common /param/ al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms

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

      end if

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine su3(n,nin,nd,x,y,u)
c------------------------------------------------------------------------c
c     reconstruct nin complex su(3) matrices u(i,3,3) from first two     c
c     rows  x(i,6),  y(i,6) (real, imaginary parts)                      c
c------------------------------------------------------------------------c

      parameter(n2=64)
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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine store(nin,nd,uold,rold,diso,pclo,
     2                 uc1o,uc2o,uc3o,uc4o)
c------------------------------------------------------------------------c
c     store distances, orientations, overlaps and gauge interaction from c
c     commom blocks /d,ucr/ in arrays uold,rolds, etc.                   c
c------------------------------------------------------------------------c
      parameter (n2=64)

      complex uc1r, uc2r, uc3r, uc4r
      complex uc1o, uc2o, uc3o, uc4o
      dimension uold(nd,nd), rold(nd,nd), diso(nd,nd,4),pclo(nd,nd)
      dimension uc1o(nd,nd), uc2o(nd,nd), uc3o(nd,nd), uc4o(nd,nd)

      common /d/   uint(n2,n2), r2(n2,n2), dis(n2,n2,4), pcl(n2,n2)
      common /ucr/ uc1r(n2,n2), uc2r(n2,n2), uc3r(n2,n2), uc4r(n2,n2)

      do 910 i = 1, nin
        do 911 ir = 1, 4
        do 911 j = 1, nin
           diso(j,i,ir) = dis(j,i,ir)
  911   continue
        do 910 j = 1, nin
          uold(j,i) = uint(j,i)
          rold(j,i) = r2(j,i)
          pclo(j,i) = pcl(j,i)
          uc1o(j,i) = uc1r(j,i)
          uc2o(j,i) = uc2r(j,i)
          uc3o(j,i) = uc3r(j,i)
          uc4o(j,i) = uc4r(j,i)
  910   continue

       return
       end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine refresh(nin,nd,uold,rold,diso,pclo,
     2                   uc1o,uc2o,uc3o,uc4o)
c------------------------------------------------------------------------c
c     restore old distances, orientations, overlaps and gauge actions    c
c     from arrays uold, rolds, etc. to common /d,ucr/.                   c
c------------------------------------------------------------------------c
      parameter (n2=64)

      complex uc1r, uc2r, uc3r, uc4r
      complex uc1o, uc2o, uc3o, uc4o

      dimension uold(nd,nd), rold(nd,nd), diso(nd,nd,4),pclo(nd,nd)
      dimension uc1o(nd,nd), uc2o(nd,nd), uc3o(nd,nd), uc4o(nd,nd)

      common /d/   uint(n2,n2), r2(n2,n2), dis(n2,n2,4), pcl(n2,n2)
      common /ucr/ uc1r(n2,n2), uc2r(n2,n2), uc3r(n2,n2), uc4r(n2,n2)

      do 910 i = 1, nin
        do 911 ir = 1, 4
        do 911 j = 1, nin
          dis(j,i,ir) = diso(j,i,ir)
  911   continue
        do 910 j = 1, nin
          uint(j,i) = uold(j,i)
          r2(j,i) = rold(j,i)
          pcl(j,i) = pclo(j,i)
          uc1r(j,i) = uc1o(j,i)
          uc2r(j,i) = uc2o(j,i)
          uc3r(j,i) = uc3o(j,i)
          uc4r(j,i) = uc4o(j,i)
  910   continue

       return
       end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

      parameter(n2=64)

      complex  uc1o, uc2o, uc3o, uc4o

      dimension zr(nd,5), e1(nd,6), e2(nd,6)
      dimension x(6), y(6), zc(5), er1(6), er2(6)
      dimension uold(n2,n2), rold(n2,n2), diso(n2,n2,4), pclo(n2,n2)
      dimension uc1o(n2,n2), uc2o(n2,n2), uc3o(n2,n2), uc4o(n2,n2)

      common /param/al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /metr/ itot,ifail1,ifail2,ifail3,actt
      common /s1/   s1
      common /acti/ acold
      common /acc/  expmax

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

         call store(nin,n2,uold,rold,diso,pclo,uc1o,uc2o,uc3o,uc4o)

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

         dels= acnew-acto
         dels= min(dels,expmax)
         dels= max(dels,-expmax)
         rat = exp(dels)
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
            call refresh(nin,n2,uold,rold,diso,pclo,
     2           uc1o,uc2o,uc3o,uc4o)
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

        call store(nin,n2,uold,rold,diso,pclo,uc1o,uc2o,uc3o,uc4o)
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

         dels= acnew-acto
         dels= min(dels,expmax)
         dels= max(dels,-expmax)
         rat = exp(dels)
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
            call refresh(nin,n2,uold,rold,diso,pclo,
     2           uc1o,uc2o,uc3o,uc4o)
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

         call store(nin,n2,uold,rold,diso,pclo,uc1o,uc2o,uc3o,uc4o)
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

         dels= acnew-acto
         dels= min(dels,expmax)
         dels= max(dels,-expmax)
         rat = exp(dels)
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
            call refresh(nin,n2,uold,rold,diso,pclo,
     2           uc1o,uc2o,uc3o,uc4o)
         end if

c------------------------------------------------------------------------c
c      end of big loop over instantons                                   c
c------------------------------------------------------------------------c

  100 continue

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

      parameter(ni=64, ni2=ni/2)
      complex clp, cdet, r, cdu, cds, cl2, clt, clf, cjc,fac
c     complex det1
      dimension clp(ni2,ni2), zr(nd,5), r(1000), cl2(ni2,ni2)
      dimension clf(ni2,ni2), cjc(ni2,ni2), fac(ni2,ni2),iptv(ni2)
      dimension ar(ni2,ni2), ai(ni2,ni2), wr(ni2), wi(ni2)
      dimension e1(nd,6), e2(nd,6)

      common /alpha/ alp,rhob,xnu
      common /param/ al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /pi/ pi, eps

c------------------------------------------------------------------------c
c     calculate gauge field action and fermionic overlap matrix elements c
c------------------------------------------------------------------------c
c     if iop .ne. 0, only updates performed                              c
c------------------------------------------------------------------------c

      nih = nin/2
      nd2 = nd/2
      call overl(n,nin,nd,nd2,zr,e1,e2,clp,bofac,rhofc,uintt,il,iu,iop)

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

      else if (rmu.eq.0.0 .and. rms.eq. 0.0) then

c------------------------------------------------------------------------c
c     fermionic determinant, massless quarks                             c
c------------------------------------------------------------------------c

         do 5  i = 1, nih
         do 5  j = 1, nih
            clf(j,i) = clp(j,i)
   5     continue

c------------------------------------------------------------------------c
c     determinant from LU factorization                                  c
c------------------------------------------------------------------------c

         call logdetnh(nih,clf,aldet)
         sdet = 2.0*nf*aldet + nf*rmul
         rdet = sdet
         sdet = alp*sdet + bofac

      else

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
         sdet = alp*sdet + bofac

      end if

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

      parameter(ni=64, ni2=ni/2, lwork=96, lrwork=128)
      complex clf, fac
      complex work
      dimension clf(ni2,ni2), fac(ni2,ni2),iptv(ni2)
      dimension work(lwork), rwork(lrwork), wr(ni2)

c------------------------------------------------------------------------c
c     imsl version                                                       c
c------------------------------------------------------------------------c

c     call lfthf(nih,clf,ni2,fac,ni2,iptv)
c     call lfdhf(nih,fac,ni2,iptv,det1,det2)

c     aldet = alog(det1) + det2*alog(10.0)

c------------------------------------------------------------------------c
c     lapack version                                                     c
c------------------------------------------------------------------------c

      call cheev('n','u',nih,clf,ni2,wr,work,lwork,rwork,info)
      if(info .ne. 0) stop

      aldet = 0.0
      do 5 i=1,nih
         aldet = aldet + alog(abs(wr(i)))
  5   continue

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine logdetnh(nih,clf,aldet)
c------------------------------------------------------------------------c
c     calculate log(det(clf)) using imsl/lapack subroutines.             c
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

      parameter(ni=64, ni2=ni/2, lwork=96, lrwork=128)
      complex clf, fac, det1, wr, vl, vr
      complex work
      dimension clf(ni2,ni2), fac(ni2,ni2),iptv(ni2)
      dimension work(lwork), rwork(lrwork), wr(ni2)
      dimension vl(ni2,ni2), vr(ni2,ni2)

c------------------------------------------------------------------------c
c     imsl version                                                       c
c------------------------------------------------------------------------c

c     call lftcg(nih,clf,ni2,fac,ni2,iptv)
c     call lfdcg(nih,fac,ni2,iptv,det1,det2)

c     aldet = alog(cabs(det1)) + det2*alog(10.0)

c------------------------------------------------------------------------c
c     lapack version                                                     c
c------------------------------------------------------------------------c

      call cgeev('n','n',nih,clf,ni2,wr,vl,ni2,vr,ni2,
     1           work,lwork,rwork,info)
      if(info .ne. 0) stop

      aldet = 0.0
      do 5 i=1,nih
         aldet = aldet + alog(cabs(wr(i)))
  5   continue

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine spect(n,nin,nd,zr,e1,e2,rmdet,wr,sdet,rv,uin,rlp,
     2                 rcos,rcii,rnii,rnia,mol)
c------------------------------------------------------------------------c
c     calculate spectrum, complete measure log(det)-S, and distribution  c
c     of relative orientations for instanton configuration. subroutine   c
c     assumes that configuration is calculated for the first time.       c
c------------------------------------------------------------------------c
c     note that spectrum is calculated even for nf=0.                    c
c------------------------------------------------------------------------c
c     n,nin,nd    number of colors, number of inst, array dimension      c
c     zr,e1,e2    instanton configuration                                c
c     rmdet       fermionic determinant                                  c
c     wr(nin/2)   fermionic eigenvalues                                  c
c     sdet        log(measure)                                           c
c     rv          log of density distribution                            c
c     uin         gauge interaction                                      c
c     rlp(i,j)    absolute value of IA overlap matrix elements           c
c     rcos(i,j)   cos of relative IA angle                               c
c     rcii(i,j)   |u_4|^2 for IA orientation                             c
c     rnii,rnia   not used                                               c
c------------------------------------------------------------------------c
c     lapack workspace requirements as follows :                         c
c     lwork = 2*n-1, lrwork = 3*n-2                                      c
c     careful, seem to need more than that:                              c
c     lwork = 3*n,   lrwork = 4*n                                        c
c------------------------------------------------------------------------c

      parameter(ni=64, ni2=ni/2, n2=ni, lwork=96,lrwork=128)
      complex   clp, clt, cjc, cl2, ctu, ctr, ctu4
      complex   uc1r, uc2r, uc3r, uc4r, eval(ni2)
      complex   work
      dimension clp(ni2,ni2), zr(nd,5), rlp(ni2,ni2), rcos(ni2,ni2)
      dimension cl2(ni2,ni2), cjc(ni2,ni2), rnii(ni2,ni2),rnia(ni2,ni2)
      dimension rcii(ni2, ni2), mol(ni2)
      dimension ar(ni2,ni2), ai(ni2,ni2), wr(ni2), wi(ni2)
      dimension e1(nd,6), e2(nd,6)
      dimension work(lwork), rwork(lrwork)

      common /alpha/ alp,rhob,xnu
      common /param/ al(4),alpha,rh0,sg,dz(4),drh,nc,nf,rmu,rms
      common /d/     uint(n2,n2),r2(n2,n2),dis(n2,n2,4),pcl(n2,n2)
      common /ucr/   uc1r(n2,n2),uc2r(n2,n2),uc3r(n2,n2),uc4r(n2,n2)
      common /pi/    pi, eps
      common /sij/   sij(ni,ni)

c------------------------------------------------------------------------c
c     do full calculation (iop=0) of overlaps and bosonic measure        c
c------------------------------------------------------------------------c

      nih = nin/2
      nd2 = nd/2
      nar = nin
      call overl(n,nin,nd,nd2,zr,e1,e2,clp,bofac,rv,uintt,1,nar,0)
      uin = uintt

c------------------------------------------------------------------------c
c     for statistics only: abs(clp) and IA orientation angles            c
c------------------------------------------------------------------------c

      do 10 i = 1, nih
         mol(i) =  1
         cmax   = 0.0
         do 10 j = 1, nih
             jp  = j+nih
             rr2 = dis(jp,i,1)**2 + dis(jp,i,2)**2
     2            +dis(jp,i,3)**2 + dis(jp,i,4)**2
             ctu =(cabs(uc1r(jp,i))**2 + cabs(uc2r(jp,i))**2
     2            +cabs(uc3r(jp,i))**2 + cabs(uc4r(jp,i))**2)/rr2
             ctr = uc4r(jp,i)
             ctu4=(uc1r(jp,i)*dis(jp,i,1)+uc1r(jp,i)*dis(jp,i,2)
     2            +uc3r(jp,i)*dis(jp,i,3)+uc4r(jp,i)*dis(jp,i,4))/
     3             rr2
             rcii(j,i) = cabs(ctu4)**2/ctu
             rcos(j,i) = cabs(ctr)**2/ctu/rr2
             dipole    = ctu-4.0*rcos(j,i)
             cjc(j,i)  = conjg(clp(j,i))
             rlp(j,i)  = cabs(clp(j,i))
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
   12 continue

c------------------------------------------------------------------------c
c     diagonalize T*T^(+), imsl                                          c
c------------------------------------------------------------------------c

c     call evlhf(nih,cl2,ni2,wr)

c------------------------------------------------------------------------c
c     diagonalize, lapack                                                c
c------------------------------------------------------------------------c

      call cheev('n','l',nih,cl2,ni2,wr,work,lwork,rwork,info)
      if(info .ne. 0) stop

c------------------------------------------------------------------------c
c     calculate log(det) for (nf-1) light plus one heavy flavor          c
c------------------------------------------------------------------------c

      rdet = 0.0
      sdet = 0.0
      srho = 0.0
      do 20 i = 1, nih
         rhp  = zr(i,5)*zr(i+nih,5)
         rms2 = rms*rms
         rmu2 = rmu*rmu
         eig  = wr(i)
         if (eig .le. 0.0) then
            write(6,*) 'danger: negative eigenvalue found'
            eig = abs(wr(i))
         endif
         srho = srho + alog(rhp)
         sdet = sdet+(nf-1)*alog(eig+rmu2) + alog(rms2+eig)
         rdet = rdet+(nf-1)*alog(eig) + alog(eig)
   20 continue
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
        sdet = alp*sdet + bofac
      end if

      do 30 i = 1, nih
         wr(i) = sqrt(abs(wr(i)))
   30 continue

      return
      end

c---------------------------------------------------------------------+----
c---------------------------------------------------------------------+----

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine norm(icos,nbin,xmin,st)
c------------------------------------------------------------------------c
c     normalize distribution of angles to sin^2(x).                      c
c------------------------------------------------------------------------c
c     icos(n) histogram array                                            c
c     nbin    number of bins                                             c
c     xmin    minimum value in histogram                                 c
c     st      bin width                                                  c
c------------------------------------------------------------------------c
      dimension icos(nbin)
      x0 = xmin+0.5*st
      f0 = sqrt((1.0-x0)/x0)
      do 10 i=1,nbin
         x = xmin+(i-0.5)*st
         x = min(x,1.0-st)
         f = sqrt((1.0-x)/x)
         icos(i) = icos(i)*f0/f
  10  continue
      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

