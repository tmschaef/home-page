
      program random
c---------------------------------------------------------------------------c
c     evaluate finite T correlation functions.                              c
c---------------------------------------------------------------------------c
c     version           :        2.7                                        c
c     last modification :      12-12-95                                     c
c---------------------------------------------------------------------------c
c     this version of tcor is based on version 3.1 of newcor. version 3.0   c
c     is completely revised and simplified, in version 3.1 all correlation  c
c     functions are calcultated in 1-direction in order to make the 4-dir   c
c     periodic.                                                             c
c---------------------------------------------------------------------------c
c     version 1.0 has new finite temperature propagator ptfull, but the     c
c     overlap matrix elements are unchanged. version 2.0 has the full       c
c     temperature dependence in the fermionic overlap matrix elements.      c
c---------------------------------------------------------------------------c
c     version 2.1 has revised parametrizations of the fermionic overlap     c
c     matrix element at finite temperature.                                 c
c---------------------------------------------------------------------------c
c     version 2.2 has new finite temperature setup. new input parameter     c
c     pfrac controls the number of molecules that are polarized in the      c
c     time direction.                                                       c
c---------------------------------------------------------------------------c
c     version 2.3: added omega and f1 correlation functions.                c
c---------------------------------------------------------------------------c
c     version 2.4: adapt program for use with both imsl (cray) and cern     c
c     (alpha) libraries. add independent random number generator. Also      c
c     modified setup to read configurations from infile.dat.                c
c---------------------------------------------------------------------------c
c     version 2.5: calculate different components of vector correlators.    c
c---------------------------------------------------------------------------c
c     version 2.6: adapt program for both temporal and spatila correlators. c
c---------------------------------------------------------------------------c
c     This version uses LAPACK subroutines for matrix diagonalization and   c
c     inversion. In order to go back to IMSL uncomment RNSET in main,       c
c     switch to RNUNF() in rang, change to EVLHF/LINCG in spect and rminv,  c
c     change to BSK1/BSKS in scalcor, masscor and stfree.                   c
c---------------------------------------------------------------------------c
c     input file: tincor.dat                                                c
c     nc,nf    number of colors, flavors                                    c
c     nin      number of instantons                                         c
c     nmol     number of molecules (cocktail model, iread=0)                c
c     rh0      instanton size                                               c
c     rmu,rms  quark masses                                                 c
c     idu,idu  not used                                                     c
c     iread    mode (0=random model, 1=input from infile.dat)               c
c     nconfig  number of configurations                                     c
c     nitc     number of points/configuration                               c
c     dummy    not used (->wave fct)                                        c
c     delt     step size for correlators                                    c
c     ndelt    number of steps                                              c
c     idu      not used                                                     c
c     ipexp    pathexponent (0/1=off/on)                                    c
c     iseed    random number seed                                           c
c     idum     not used                                                     c
c     idir     direction in which correlators are calculated                c
c     alb(4)   box size                                                     c
c---------------------------------------------------------------------------c
c     note: alb(4) is interpreted as inverse temperature, except if not=1   c
c     in main (then beta=1000).                                             c
c---------------------------------------------------------------------------c
c     parameters : n    number of colors                                    c
c                  ni   number of instantons                                c
c                  ncf  number of configurations                            c
c                  ndl  number of deltas for correlator                     c
c---------------------------------------------------------------------------c

      parameter(n=3, ni=256, ni2=ni/2, nbin=150, ncf=100, ndl = 20)

      dimension zr(ni,5),e1(ni,6),e2(ni,6), wr(ni2), xxsp(4,2)
      dimension xxsm(4,2)
      complex g, gm5
      complex lambda, lf
      complex clp, clps, clpu(ni,ni)
      complex cotr, cotrs
      complex u, u1
      complex cprod, ceig, u0, v, csom
      complex ca0, cb0, cat, cbt, catt, cbtt, clap
      complex cj, cj0, cjt, cjtt, uus, snew, sold
      complex sfree, szm, uuu, ca, cb, s, sm, condt
      complex szxx, szyy, conddd
      complex as(16,3,3), sas(16), sast(16,ndl,ncf), aas(ndl,16)
      complex sab(ndl,ncf)
      complex cl(6), cdl(4), cdiq(5)
      complex clt(6), cdlt(10), cdiqt(5), cdiqtt(5,ndl,ncf)

c-------------------------------------------------------------------------c
c     propagators                                                         c
c-------------------------------------------------------------------------c

      complex sfr(3,4,3,4)
      complex sfullm(3,4,3,4),sfullp(3,4,3,4)
      complex sxy(3,4,3,4),syx(3,4,3,4)
      complex sax(3,4,3,4),say(3,4,3,4),saz(3,4,3,4)
      complex ssxy(3,4,3,4),ssyx(3,4,3,4)
      complex sxx(3,4,3,4), syy(3,4,3,4)
      complex pax(3,3),pyx(3,3),pexp(3,3),pexp2(3,3)
      dimension point(4), axx(4,2), ayy(4,2)
      dimension ndis(ndl), ist(100), shu(ndl)

c-------------------------------------------------------------------------c
c     arrays needed for meson correlators                                 c
c-------------------------------------------------------------------------c

      real kpmy,kmmy,kp2my,km2my,kpmya,kmmya
      real kpmye,kmmye,kpmyep,kmmyep
      dimension scmy(ndl,ncf),psmy(ndl,ncf)
      dimension axmy(ndl,ncf),vemy(ndl,ncf)
      dimension sc2my(ndl,ncf),ps2my(ndl,ncf)
      dimension ax2my(ndl,ncf),ve2my(ndl,ncf)
      dimension scmya(ndl),psmya(ndl),axmya(ndl),vemya(ndl)
      dimension scmye(ndl),psmye(ndl),axmye(ndl),vemye(ndl)
      dimension scmyep(ndl),psmyep(ndl),axmyep(ndl),vemyep(ndl)
      dimension pemy(ndl),pe2my(ndl),pemya(ndl),pemye(ndl)
      dimension kpmy(ndl,ncf),kmmy(ndl,ncf)
      dimension kp2my(ndl,ncf),km2my(ndl,ncf)
      dimension kpmya(ndl),kmmya(ndl)
      dimension kpmye(ndl),kmmye(ndl)
      dimension kpmyep(ndl),kmmyep(ndl)

c------------------------------------------------------------------------c
c     spatial, temporal components of vector correlators                 c
c------------------------------------------------------------------------c

      dimension avmy(ndl,ncf),vvmy(ndl,ncf)
      dimension av2my(ndl,ncf),vv2my(ndl,ncf)
      dimension avmya(ndl),vvmya(ndl)
      dimension avmye(ndl),vvmye(ndl)
      dimension avmyep(ndl),vvmyep(ndl)

      dimension a4my(ndl,ncf),v4my(ndl,ncf)
      dimension a42my(ndl,ncf),v42my(ndl,ncf)
      dimension a4mya(ndl),v4mya(ndl)
      dimension a4mye(ndl),v4mye(ndl)
      dimension a4myep(ndl),v4myep(ndl)

      dimension a1my(ndl,ncf),v1my(ndl,ncf)
      dimension a12my(ndl,ncf),v12my(ndl,ncf)
      dimension a1mya(ndl),v1mya(ndl)
      dimension a1mye(ndl),v1mye(ndl)
      dimension a1myep(ndl),v1myep(ndl)

c------------------------------------------------------------------------c
c     unflavored mesons                                                  c
c------------------------------------------------------------------------c

      dimension sigmy(ndl,ncf),sig2my(ndl,ncf)
      dimension sigmya(ndl),sigmye(ndl),sigmyep(ndl)
      dimension etamy(ndl,ncf),eta2my(ndl,ncf)
      dimension etamya(ndl),etamye(ndl),etamyep(ndl)
      dimension f1my(ndl,ncf),f12my(ndl,ncf)
      dimension f1mya(ndl),f1mye(ndl),f1myep(ndl)
      dimension ommy(ndl,ncf),om2my(ndl,ncf)
      dimension ommya(ndl),ommye(ndl),ommyep(ndl)

c------------------------------------------------------------------------c
c     off diagonal correlators                                           c
c------------------------------------------------------------------------c

      dimension vsmy(ndl,ncf),pamy(ndl,ncf)
      dimension vtmy(ndl,ncf),atmy(ndl,ncf)
      dimension vs2my(ndl,ncf),pa2my(ndl,ncf)
      dimension vt2my(ndl,ncf),at2my(ndl,ncf)
      dimension vsmya(ndl),pamya(ndl)
      dimension vtmya(ndl),atmya(ndl)
      dimension vsmye(ndl),pamye(ndl)
      dimension vtmye(ndl),atmye(ndl)
      dimension vsmyep(ndl),pamyep(ndl)
      dimension vtmyep(ndl),atmyep(ndl)

c------------------------------------------------------------------------c
c     diquarks, baryons                                                  c
c------------------------------------------------------------------------c

      dimension xc(ndl,ncf),xc2(ndl,ncf)
      dimension xca(ndl),xce(ndl)
      dimension xcep(ndl)
      dimension clmy(6,ndl,ncf),cdlmy(4,ndl,ncf)
      dimension cdiqmy(5,ndl,ncf)
      dimension cl2my(6,ndl,ncf),cdl2my(4,ndl,ncf)
      dimension cdiq2my(5,ndl,ncf)
      dimension clmya(6,ndl),cdlmya(4,ndl)
      dimension cdiqmya(5,ndl)
      dimension clmye(6,ndl),cdlmye(4,ndl)
      dimension cdiqmye(5,ndl)
      dimension clmyep(6,ndl),cdlmyep(4,ndl)
      dimension cdiqmyep(5,ndl)
      dimension cdfree(5),clfree(6),cdlfree(4)

c------------------------------------------------------------------------c
c     multi strange mesons                                               c
c------------------------------------------------------------------------c

      dimension sscmy(ndl,ncf),spsmy(ndl,ncf)
      dimension saxmy(ndl,ncf),svemy(ndl,ncf)
      dimension ssc2my(ndl,ncf),sps2my(ndl,ncf)
      dimension sax2my(ndl,ncf),sve2my(ndl,ncf)
      dimension sscmya(ndl),spsmya(ndl),saxmya(ndl),svemya(ndl)
      dimension sscmye(ndl),spsmye(ndl),saxmye(ndl),svemye(ndl)
      dimension sscmyep(ndl),spsmyep(ndl),saxmyep(ndl),svemyep(ndl)

c------------------------------------------------------------------------c
c     multi strange baryons                                              c
c------------------------------------------------------------------------c

      dimension sclmy(6,ndl,ncf),scdlmy(4,ndl,ncf)
      dimension sscdlmy(4,ndl,ncf),scdiqmy(5,ndl,ncf)
      dimension scl2my(6,ndl,ncf),scdl2my(4,ndl,ncf)
      dimension sscdl2my(4,ndl,ncf),scdiq2my(5,ndl,ncf)
      dimension sclmya(6,ndl),scdlmya(4,ndl)
      dimension sscdlmya(4,ndl),scdiqmya(5,ndl)
      dimension sclmye(6,ndl),scdlmye(4,ndl)
      dimension sscdlmye(4,ndl),scdiqmye(5,ndl)
      dimension sclmyep(6,ndl),scdlmyep(4,ndl)
      dimension sscdlmyep(4,ndl),scdiqmyep(5,ndl)

c------------------------------------------------------------------------c
c     lambdas                                                            c
c------------------------------------------------------------------------c

      dimension cclmy(6,ndl,ncf),dclmy(6,ndl,ncf)
      dimension ccl2my(6,ndl,ncf),dcl2my(6,ndl,ncf)
      dimension cclmya(6,ndl),dclmya(6,ndl)
      dimension cclmye(6,ndl),dclmye(6,ndl)
      dimension cclmyep(6,ndl),dclmyep(6,ndl)

c----------------------------------------------------------------------c
c     common blocks                                                    c
c----------------------------------------------------------------------c

      common /gamf5/ gm5(4,4), nm5(4,4)
      common /gamf/  g(5,4), nu(5,4), nv(5,4), c(4), nuc(4), nic(4)
      common /lam/   lambda(8,3,3)
      common /lamf/  lf(8,3), lu(8,3), li(8,3)
      common /param/ a, alpha,rh0,sg,dz,drh, nc, nf, rmu, rms
      common /cden/  b,alc, p1, p2
      common /metr/  itot, ifail, actt
      common /acti/  acold
      common /pi/    pi, eps
      common /sij/   sij(ni,ni)
      common /clp/   clp(ni,ni), clps(ni,ni)
      common /wol/   wol, itow, isuc
      common /xo/    xo(4)
      common /grint/ wdt, ddx
      common /temp/  keyt,al4,temp
      common /u1/    u1(ni,3,3)
      common /sweight/ sweight
      common /const/ const
      common /nconf/ nconfigB
      common /box/   alb(4)
      common /clap/  clap(ni,ni)
      common /counter/ icount
      common /seed/  iseed
      common /c1c2/  c1, c2
      common /scalmass/ scalmass
c     common /worksp/rwksp
      common /cold/  not
      common /direc/ idir

c----------------------------------------------------------------------------c
c     allocate workspace for diagonalization (66068 on the cray)             c
c----------------------------------------------------------------------------c

c     real rwksp(263188)
c     call iwkin(263188)

c----------------------------------------------------------------------------c
c     input parameters                                                       c
c----------------------------------------------------------------------------c

      not  = 0
      nbins= 40
      xmin = 0.0
      xmax = 2.0

      open (unit=1, file='tincor.dat', status='old')
      read (1,*) nc, nf, nin, nmol, pfrac, rh0, scalmass
      read (1,*) rmu,rms, ieq, kp1
      read (1,*) iread, nconfig
      read (1,*) nitc, delfix, delt , ndelt
      read (1,*) ntrans, ipexp, iseed
      read (1,*) igeom, idir
      read (1,*) (alb(k), k = 1, 4)
      close (unit=1)

c----------------------------------------------------------------------------c
c     echo input parameters                                                  c
c----------------------------------------------------------------------------c

      open (unit=2, file='toutcor.dat', status='unknown')
      write(2,*) ' tcor, version 2.7'
      write(2,*) ' -----------------'
      write(2,501) nc, nf, nin
      write(2,502) nmol, pfrac, rh0
      write(2,503) rmu, rms, scalmass
      write(2,504) ieq, kp1, iread, nconfig
      write(2,505) nitc, delt, ndelt
      write(2,506) ipexp, idir, iseed
      write(2,507) (alb(k),k=1,4)
      write(2,*)

  501 format(1x,' N_c = ',i5,5x,' N_f = ',i5,5x,' N_in = ',i5)
  502 format(1x,' N_mo= ',i5,5x,' pfr = ',f10.4,' rho  = ',f10.4)
  503 format(1x,' m_u = ',f10.4,' m_s = ',f10.4,' M_sc = ',f10.4)
  504 format(1x,' ieq = ',i5,5x,' ikp = ',i5,5x,
     1          ' irea= ',i5,5x,' ncon= ',i5)
  505 format(1x,' nitc= ',i5,5x,' delt= ',f10.4,' N_del= ',i5)
  506 format(1x,' ipex= ',i5,5x,' idir= ',i5,5x,' iseed= ',i10)
  507 format(1x,' a_1 = ',f10.4,' a_2 = ',f10.4,
     2          ' a_3 = ',f10.4,' a_4 = ',f10.4)

c-------------------------------------------------------------------------c
c     flush random number generator                                       c
c-------------------------------------------------------------------------c

      jseed = iseed
c     call ranget(iseed)
c     call rnset(iseed)
      do 45 i = 1, jseed
        dum = rang( )
  45  continue
      print *, 'jseed =',jseed
      keyt = 0
      al4 = a
      volume = alb(1)*alb(2)*alb(3)*alb(4)
      ndd = ni

c---------------------------------------------------------------------------c
c     initialize constants                                                  c
c---------------------------------------------------------------------------c

      pi = 3.1415926
      c1 = 3*pi/8
      c2 = (3*pi/32)**1.333333333
      const=1/(2*pi**2)
      eps = exp(-40*alog(2.0))
      b = 11.0/3.0*nc -2.0/3.0*nf
      bp = 34.0/3.0*nc*nc - 13.0/3.0*nc*nf +nf/float(nc)
      pp1 = 2*nc-bp/2/b
      pp2 = bp/2/b
      cnc = 1.34**nf *4.66*exp(-1.68*nc)/pi/pi/(nc-1)*(b/2)**(bp/2/b)
      alc = alog(cnc)
      sweight = 1.0
      write(2,121)  jseed
  121 format(1x, ' iseed = ', i12,/)

      nih = nin/2
c---------------------------------------------------------------------------c
c     initialize gamma and tau matrices                                     c
c---------------------------------------------------------------------------c

      call gammat
      call gamfast
      call gamlr
      call lammat
      call taumat

      icon = 0
      ncon = 0

c---------------------------------------------------------------------------c
c     clear summation arrays used for correlators                           c
c---------------------------------------------------------------------------c

      call zero2(ndl,ncf,scmy)
      call zero2(ndl,ncf,psmy)
      call zero2(ndl,ncf,axmy)
      call zero2(ndl,ncf,vemy)

      call zero2(ndl,ncf,sc2my)
      call zero2(ndl,ncf,ps2my)
      call zero2(ndl,ncf,ax2my)
      call zero2(ndl,ncf,ve2my)

      call zero2(ndl,ncf,vvmy)
      call zero2(ndl,ncf,avmy)
      call zero2(ndl,ncf,v4my)
      call zero2(ndl,ncf,a4my)
      call zero2(ndl,ncf,v1my)
      call zero2(ndl,ncf,a1my)

      call zero2(ndl,ncf,vv2my)
      call zero2(ndl,ncf,av2my)
      call zero2(ndl,ncf,v42my)
      call zero2(ndl,ncf,a42my)
      call zero2(ndl,ncf,v12my)
      call zero2(ndl,ncf,a12my)

      call zero2(ndl,ncf,sigmy)
      call zero2(ndl,ncf,etamy)
      call zero2(ndl,ncf,f1my)
      call zero2(ndl,ncf,ommy)

      call zero2(ndl,ncf,sig2my)
      call zero2(ndl,ncf,eta2my)
      call zero2(ndl,ncf,f12my)
      call zero2(ndl,ncf,om2my)

      call zero2(ndl,ncf,sscmy)
      call zero2(ndl,ncf,spsmy)
      call zero2(ndl,ncf,saxmy)
      call zero2(ndl,ncf,svemy)

      call zero2(ndl,ncf,ssc2my)
      call zero2(ndl,ncf,sps2my)
      call zero2(ndl,ncf,sax2my)
      call zero2(ndl,ncf,sve2my)

      call zero2(ndl,ncf,vsmy)
      call zero2(ndl,ncf,pamy)
      call zero2(ndl,ncf,vtmy)
      call zero2(ndl,ncf,atmy)

      call zero2(ndl,ncf,vs2my)
      call zero2(ndl,ncf,pa2my)
      call zero2(ndl,ncf,vt2my)
      call zero2(ndl,ncf,at2my)

      call zero2(ndl,ncf,kpmy)
      call zero2(ndl,ncf,kmmy)
      call zero2(ndl,ncf,kp2my)
      call zero2(ndl,ncf,km2my)

      call zero3(6,ndl,ncf,clmy)
      call zero3(4,ndl,ncf,cdlmy)
      call zero3(5,ndl,ncf,cdiqmy)

      call zero3(6,ndl,ncf,cl2my)
      call zero3(4,ndl,ncf,cdl2my)
      call zero3(5,ndl,ncf,cdiq2my)

      call zero3(6,ndl,ncf,sclmy)
      call zero3(6,ndl,ncf,cclmy)
      call zero3(6,ndl,ncf,dclmy)
      call zero3(4,ndl,ncf,scdlmy)
      call zero3(4,ndl,ncf,sscdlmy)
      call zero3(5,ndl,ncf,scdiqmy)

      call zero3(6,ndl,ncf,scl2my)
      call zero3(4,ndl,ncf,scdl2my)
      call zero3(4,ndl,ncf,sscdl2my)
      call zero3(5,ndl,ncf,scdiq2my)

      call zero1(ndl,pemy)
      call zero1(ndl,pe2my)

      call izero(nbins,ist)

      xtij  = 0.0
      xtij2 = 0.0
      xdet  = 0.0
      xdet2 = 0.0

      xqq   = 0.0
      xqq2  = 0.0
      xuu   = 0.0
      xuu2  = 0.0
      xuudd = 0.0
      xuudd2= 0.0
      xuuuu = 0.0
      xuuuu2= 0.0

      xua1 = 0.0
      xua12= 0.0
      xqg  = 0.0
      xqg2 = 0.0

      xo1  = 0.0
      xo12 = 0.0
      xo2  = 0.0
      xo22 = 0.0
      xo3  = 0.0
      xo32 = 0.0
      xo4  = 0.0
      xo42 = 0.0

      xss  = 0.0
      xss2 = 0.0
      xpp  = 0.0
      xpp2 = 0.0
      xvv  = 0.0
      xvv2 = 0.0
      xaa  = 0.0
      xaa2 = 0.0
      xtt  = 0.0
      xtt2 = 0.0

c---------------------------------------------------------------------------c
c     clear counter for correlator statistics                               c
c---------------------------------------------------------------------------c

      do 344 k = 1, ndelt
         ndis(k)=0
344   continue

c---------------------------------------------------------------------------c
c     loop counters for config,                                             c
c---------------------------------------------------------------------------c

      icount = 0
      ic     = 0

c---------------------------------------------------------------------------c
c     initialize pexps to 1 (used in heavy-light,diqu.)                     c
c---------------------------------------------------------------------------c

      do 342 i=1,3
      do 342 j=1,3
         pax(i,j)  = (0.0,0.0)
         pyx(i,j)  = (0.0,0.0)
         pexp(i,j) = (0.0,0.0)
         pexp2(i,j)= (0.0,0.0)
 342  continue
      do 343 i=1,3
         pax(i,i)  = (1.0,0.0)
         pyx(i,i)  = (1.0,0.0)
         pexp(i,i) = (1.0,0.0)
         pexp2(i,i)= (1.0,0.0)
 343  continue

c---------------------------------------------------------------------------c
c     loop over different configurations                                    c
c---------------------------------------------------------------------------c

      do 180 ic=1,nconfig

c---------------------------------------------------------------------------c
c     generate a new configuration                                          c
c---------------------------------------------------------------------------c

      call setup(nc,nin,nmol,pfrac,zr,e1,e2, iread, ic)

c---------------------------------------------------------------------------c
c     calculate overlap matrix elements, save clpu                          c
c---------------------------------------------------------------------------c

      call rmtinv(nc,nin,ndd,zr,e1,e2,rmu,rms)

      do 818 i = 1, nin
      do 818 j = 1, nin
         clpu(j,i) = clp(j,i)
 818  continue

      cotr = cmplx(0.0, 0.0)
      cotrs= cmplx(0.0, 0.0)
      do 330 i = 1, nin
         cotr = cotr + clp(i,i)
         cotrs= cotrs+ clps(i,i)
 330  continue

c---------------------------------------------------------------------------c
c     calculate light, strange quark condensates from tr((T+im)^-1)         c
c---------------------------------------------------------------------------c

      qbarq = aimag(cotr)/volume
      call myaddto(qbarq,xqq,xqq2)
      write (2,901) cotr/volume
      print 901, cotr/volume
      write(2,902) cotrs/volume
      print 902, cotrs/volume
 901  format(1x, ' Tr(1/(T+im))/vol:   ', 2f12.5)
 902  format(1x, ' Tr(1/(T+im_s))/vol: ', 2f12.5)

c--------------------------------------------------------------------------c
c     calculate spectrum, average fermion determinant                      c
c--------------------------------------------------------------------------c

      call tspect(nc,nin,ndd,zr,e1,e2,rmdet,wr,sdet,t2av,uin)
      write(2,903) sdet
 903  format(1x, ' log det(iD+im):    ', f12.5)
      call myaddto(sdet,xdet,xdet2)
      write(2,904) t2av
 904  format(1x, ' <\sum_j |T_ij|^2>: ', f12.5)
      call myaddto(t2av,xtij,xtij2)

c--------------------------------------------------------------------------c
c     plot eigenvalue distribution                                         c
c--------------------------------------------------------------------------c
c
c      nbins = 20
c      call bin(nih,nbins,wr,ist,xmin,xmax,1)
c      st   = (xmax-xmin)/float(nbins)
c      call lev(xmin,st,nbins,1,ist)
c
c---------------------------------------------------------------------------c
c     or: include eigenvalues in histogram, plot later                      c
c---------------------------------------------------------------------------c

      call addbin(nih,nbins,wr,xmin,xmax,ist)

c---------------------------------------------------------------------------c
c     loop over different points in a given configuration                   c
c---------------------------------------------------------------------------c

      do 100 i = 1, nitc
         itel = itel + 1
         delmax = ndelt*delt

         call getx(point,delmax)

c---------------------------------------------------------------------------c
c     local zero mode propagator at x (origin remains fixed in loop over k) c
c---------------------------------------------------------------------------c

         do 605 m=1,4
            axx(m,1) = point(m)
            axx(m,2) = point(m)
 605     continue
         call ztmodes(nc,nin,ndd,zr,e1,e2,axx,rmu,sfr,sxx,snorm)

c---------------------------------------------------------------------------c
c     loop over different separations at which correlator is measured       c
c---------------------------------------------------------------------------c

         do 200 k = 1, ndelt
         delfix = k*delt

c---------------------------------------------------------------------------c
c     select endpoints for propagators                                      c
c---------------------------------------------------------------------------c

         call getsplit(point,xxsp,delfix)

c---------------------------------------------------------------------------c
c     calculate propagators (non-strange and strange) and pyx               c
c---------------------------------------------------------------------------c

         do 817 il = 1, nin
         do 817 jl = 1, nin
             clp(jl,il) = clpu(jl,il)
  817    continue

         call ptfull(nin,ndd,zr,e1,e2,xxsp,rmu,sfr,sfullp)

c         do 816 il = 1, nin
c         do 816 jl = 1, nin
c             clp(jl,il) = clps(jl,il)
c  816    continue
c
c         call ptfull(nin,ndd,zr,e1,e2,xxsp,rms,sfr,sfullm)

         if ( ipexp .eq. 0 ) goto 209
         call pathexp(nin,ndd,zr,e1,e2,xxsp,xxsp,pyx)
 209     continue

c---------------------------------------------------------------------------c
c     local zero mode propagator at y for disconnected correlators          c
c---------------------------------------------------------------------------c

         do 610 m = 1, 4
            ayy(m,1) = xxsp(m,2)
            ayy(m,2) = xxsp(m,2)
  610    continue

         call ztmodes(nc,nin,ndd,zr,e1,e2,ayy,rmu,sfr,syy,snorm)
c        do 611 k1=1,3
c        do 611 k2=1,3
c        do 611 i1=1,4
c        do 611 i2=1,4
c           syy(k1,i1,k2,i2) = syy(k1,i1,k2,i2)/10000.0
c 611    continue
c        do 612 k1=1,3
c        do 612 i1=1,4
c           syy(k1,i1,k1,i1) = (0.0,1.0)/12.0
c 612    continue

c---------------------------------------------------------------------------c
c     check : use free propagator                                           c
c---------------------------------------------------------------------------c
c
c        call pfree(xxsp,0.0,sfullp)
c        call pfree(xxsp,rms,sfullm)
c
c---------------------------------------------------------------------------c
c     propagators for meson/baryon correlators                              c
c---------------------------------------------------------------------------c

         do 920 k1 = 1, 3
         do 920 l1 = 1, 4
         do 920 k2 = 1, 3
         do 920 l2 = 1, 4
            sxy(k1,l1,k2,l2) = sfullp(k1,l1,k2,l2)
c            ssxy(k1,l1,k2,l2)= sfullm(k1,l1,k2,l2)
  920    continue

         do 930 k1 = 1, 3
         do 930 l1 = 1, 4
         do 930 k2 = 1, 3
         do 930 l2 = 1, 4
            syx(k1,l1,k2,l2) =-g(5,l1)
     1         *conjg(sxy(k2,nv(5,l2),k1,nu(5,l1)))*g(5,nv(5,l2))
c            ssyx(k1,l1,k2,l2)=-g(5,l1)
c     1         *conjg(ssxy(k2,nv(5,l2),k1,nu(5,l1)))*g(5,nv(5,l2))
  930    continue

         do 932 k1 = 1, 3
         do 932 k2 = 1, 3
            pax(k1,k2) = conjg(pyx(k2,k1))
  932    continue

c---------------------------------------------------------------------------c
c     calculate higher condensates                                          c
c---------------------------------------------------------------------------c

         call cond(syy,uu,uudd,uuuu,ua1)
         call myaddto(uu,xuu,xuu2)
         call myaddto(uudd,xuudd,xuudd2)
         call myaddto(uuuu,xuuuu,xuuuu2)
         call myaddto(ua1,xua1,xua12)

c---------------------------------------------------------------------------c
c     weinberg operators                                                    c
c---------------------------------------------------------------------------c

         call moreops(syy,o1,o2,o3,o4)
         call myaddto(o1,xo1,xo12)
         call myaddto(o2,xo2,xo22)
         call myaddto(o3,xo3,xo32)
         call myaddto(o4,xo4,xo42)

c---------------------------------------------------------------------------c
c     colored four quark condensates                                        c
c---------------------------------------------------------------------------c

         call fourq(syy,ss,pp,vv,aa,tt)
         call myaddto(ss,xss,xss2)
         call myaddto(pp,xpp,xpp2)
         call myaddto(vv,xvv,xvv2)
         call myaddto(aa,xaa,xaa2)

c---------------------------------------------------------------------------c
c     calculate correlators, include into statistical sample                c
c---------------------------------------------------------------------------c

         ndis(k) = ndis(k)+1

         call prmeson(sxy,syx,pexp,prsc,prps,praxt,prvect,
     1                prav,prvv,pra4,prv4,pra1,prv1)

         call myaddto(prsc,scmy(k,ic),sc2my(k,ic))
         call myaddto(prps,psmy(k,ic),ps2my(k,ic))
         call myaddto(praxt,axmy(k,ic),ax2my(k,ic))
         call myaddto(prvect,vemy(k,ic),ve2my(k,ic))

         call myaddto(prav,avmy(k,ic),av2my(k,ic))
         call myaddto(prvv,vvmy(k,ic),vv2my(k,ic))
         call myaddto(pra4,a4my(k,ic),a42my(k,ic))
         call myaddto(prv4,v4my(k,ic),v42my(k,ic))
         call myaddto(pra1,a1my(k,ic),a12my(k,ic))
         call myaddto(prv1,v1my(k,ic),v12my(k,ic))

c--------------------------------------------------------------------------c
c     disconnected part of scalar correlators                              c
c--------------------------------------------------------------------------c

         call twoloop(sxx,syy,prt,prt5,prv,pra)

         prsig = prsc - 2.0*prt
         preta = prps - 2.0*prt5
         prf1  = praxt - 2.0*pra
         prom  = prvect- 2.0*prv

         call myaddto(prsig,sigmy(k,ic),sig2my(k,ic))
         call myaddto(preta,etamy(k,ic),eta2my(k,ic))
         call myaddto(prf1,f1my(k,ic),f12my(k,ic))
         call myaddto(prom,ommy(k,ic),om2my(k,ic))

c--------------------------------------------------------------------------c
c     trace of schwinger factor                                            c
c--------------------------------------------------------------------------c

         prtr = real(pax(1,1)+pax(2,2)+pax(3,3))/3.0
         call myaddto(prtr,pemy(k),pe2my(k))

c--------------------------------------------------------------------------c
c     doubly strange mesons                                                c
c--------------------------------------------------------------------------c
c
c         call prmeson(ssxy,ssyx,pexp,prsc,prps,praxt,prvect)
c
c         call myaddto(prsc,sscmy(k,ic),ssc2my(k,ic))
c         call myaddto(prps,spsmy(k,ic),sps2my(k,ic))
c         call myaddto(praxt,saxmy(k,ic),sax2my(k,ic))
c         call myaddto(prvect,svemy(k,ic),sve2my(k,ic))
c
c--------------------------------------------------------------------------c
c     use off-diagonal correlators                                         c
c--------------------------------------------------------------------------c

         call prmix(sxy,syx,pexp,pvs,ppa,pvt,pat)

         call myaddto(pvs,vsmy(k,ic),vs2my(k,ic))
         call myaddto(ppa,pamy(k,ic),pa2my(k,ic))
         call myaddto(pvt,vtmy(k,ic),vt2my(k,ic))
         call myaddto(pat,atmy(k,ic),at2my(k,ic))

c--------------------------------------------------------------------------c
c     heavy-light mesons                                                   c
c--------------------------------------------------------------------------c

         call prhlm(sxy,pyx,xkp,xkm)

         call myaddto(xkp,kpmy(k,ic),kp2my(k,ic))
         call myaddto(xkm,kmmy(k,ic),km2my(k,ic))

c--------------------------------------------------------------------------c
c     diquark correlator                                                   c
c--------------------------------------------------------------------------c

         call prdiq(sxy,sxy,pax,pexp,cdiq)

         do 935 idiq=1,5
         xcdiq = real(cdiq(idiq))
         call myaddto(xcdiq,cdiqmy(idiq,k,ic),cdiq2my(idiq,k,ic))
  935    continue

c--------------------------------------------------------------------------c
c     doubly strange diquark                                               c
c--------------------------------------------------------------------------c
c
c         call prdiq(ssxy,ssxy,pax,pexp,cdiq)
c
c         do 936 idiq=1,5
c         xcdiq = real(cdiq(idiq))
c         call myaddto(xcdiq,scdiqmy(idiq,k,ic),scdiq2my(idiq,k,ic))
c  936    continue
c
c--------------------------------------------------------------------------c
c     nucleon correlator                                                   c
c--------------------------------------------------------------------------c

         call nprop(sxy,sxy,sxy,pexp,pexp,cl)

         do 940 in=1,6
         xcl = aimag(cl(in))
         call myaddto(xcl,clmy(in,k,ic),cl2my(in,k,ic))
  940    continue

c--------------------------------------------------------------------------c
c     lambda correlator                                                     c
c--------------------------------------------------------------------------c
c
c         call lprop(sxy,sxy,ssxy,pexp,pexp,cl)
c
c         do 941 in=1,6
c         xcl = aimag(cl(in))
c         call myaddto(xcl,cclmy(in,k,ic),ccl2my(in,k,ic))
c  941    continue
c
c--------------------------------------------------------------------------c
c     more lambda correlators                                              c
c--------------------------------------------------------------------------c
c
c         call lprop2(sxy,sxy,ssxy,pexp,pexp,cl)
c
c         do 942 in=1,2
c         xcl = aimag(cl(in))
c         call myaddto(xcl,dclmy(in,k,ic),dcl2my(in,k,ic))
c  942    continue
c
c--------------------------------------------------------------------------c
c     xi correlator                                                        c
c--------------------------------------------------------------------------c
c
c         call nprop(ssxy,sxy,ssxy,pexp,pexp,cl)
c
c         do 943 in=1,6
c         xcl = aimag(cl(in))
c         call myaddto(xcl,sclmy(in,k,ic),scl2my(in,k,ic))
c  943    continue
c
c--------------------------------------------------------------------------c
c     delta correlator                                                     c
c--------------------------------------------------------------------------c

         call dprop(sxy,sxy,sxy,pexp,pexp,cdl)

         do 950 id=1,4
         xcdl = aimag(cdl(id))
         call myaddto(xcdl,cdlmy(id,k,ic),cdl2my(id,k,ic))
  950    continue

c--------------------------------------------------------------------------c
c     xi^* correlator                                                      c
c--------------------------------------------------------------------------c
c
c         call delbs(ssxy,sxy,pexp,cdl)
c
c         do 951 id=1,4
c         xcdl = aimag(cdl(id))
c         call myaddto(xcdl,scdlmy(id,k,ic),scdl2my(id,k,ic))
c  951    continue
c
c--------------------------------------------------------------------------c
c     omega correlator                                                     c
c--------------------------------------------------------------------------c
c
c         call delbs(ssxy,ssxy,pexp,cdl)
c
c         do 952 id=1,4
c         xcdl = aimag(cdl(id))
c         call myaddto(xcdl,sscdlmy(id,k,ic),sscdl2my(id,k,ic))
c  952    continue
c
c----------------------------------------------------------------------------c
c     end of loop over splittings                                            c
c----------------------------------------------------------------------------c

  200    continue

c----------------------------------------------------------------------------c
c     end of loop over points in fixed configuration                         c
c----------------------------------------------------------------------------c

  100 continue

c----------------------------------------------------------------------------c
c     generate next configuration                                            c
c----------------------------------------------------------------------------c

      ncon = ncon + 1
 180  continue

c----------------------------------------------------------------------------c
c     monte carlo finished, do statistics; start with condensates            c
c----------------------------------------------------------------------------c

      call mydisp(nconfig,xdet,xdet2,xdeta,xdete)
      call mydisp(nconfig,xtij,xtij2,xtija,xtije)

      call mydisp(nconfig,xqq,xqq2,qq,qqe)
      qbarq = qq**(1./3.)*197.3
      qbarqe= qbarq/3.*qqe/qq

      ntot = nconfig*nitc*ndelt
      call mydisp(ntot,xuu,xuu2,uua,uue)
      call mydisp(ntot,xuudd,xuudd2,uudda,uudde)
      call mydisp(ntot,xuuuu,xuuuu2,uuuua,uuuue)
      call mydisp(ntot,xua1,xua12,ua1a,ua1e)

      uu   = sign(1.0,uua)*abs(uua)**(1.0/3.0)*197.3
      uue  = uu/3.0*uue/uua
      uudd = sign(1.0,uudda)*abs(uudda)**(1.0/6.0)*197.3
      uudde= uudd/6.0*uudde/uudda
      uuuu = sign(1.0,uuuua)*abs(uuuua)**(1.0/6.0)*197.3
      uuuue= uuuu/6.0*uuuue/uuuua
      ua1  = sign(1.0,ua1a)*abs(ua1a)**(1.0/6.0)*197.3
      ua1e = ua1/6.0*ua1e/ua1a

c----------------------------------------------------------------------------c
c     weinberg operators                                                     c
c----------------------------------------------------------------------------c

      call mydisp(ntot,xo1,xo12,o1a,o1e)
      call mydisp(ntot,xo2,xo22,o2a,o2e)
      call mydisp(ntot,xo3,xo32,o3a,o3e)
      call mydisp(ntot,xo4,xo42,o4a,o4e)

      o1 = sign(1.0,o1a)*abs(o1a)**(1./6.)*197.3
      o1e= o1/6.0*o1e/o1a
      o2 = sign(1.0,o2a)*abs(o2a)**(1./6.)*197.3
      o2e= o2/6.0*o2e/o2a
      o3 = sign(1.0,o3a)*abs(o3a)**(1./6.)*197.3
      o3e= o3/6.0*o3e/o3a
      o4 = sign(1.0,o4a)*abs(o4a)**(1./6.)*197.3
      o4e= o4/6.0*o4e/o4a

c----------------------------------------------------------------------------c
c     four quark condensates                                                 c
c----------------------------------------------------------------------------c

      call mydisp(ntot,xss,xss2,ssa,sse)
      call mydisp(ntot,xpp,xpp2,ppa,ppe)
      call mydisp(ntot,xvv,xvv2,vva,vve)
      call mydisp(ntot,xaa,xaa2,aaa,aae)

      ss = sign(1.0,ssa)*abs(ssa)**(1./6.)*197.3
      sse= ss/6.0*sse/ssa
      pp = sign(1.0,ppa)*abs(ppa)**(1./6.)*197.3
      ppe= pp/6.0*ppe/ppa
      vv = sign(1.0,vva)*abs(vva)**(1./6.)*197.3
      vve= vv/6.0*vve/vva
      aa = sign(1.0,aaa)*abs(aaa)**(1./6.)*197.3
      aae= aa/6.0*aae/aaa

c----------------------------------------------------------------------------c
c     meson correlators                                                      c
c----------------------------------------------------------------------------c

      call disp(ndelt,nconfig,nitc,scmy,sc2my,scmya,scmye,scmyep)
      call disp(ndelt,nconfig,nitc,psmy,ps2my,psmya,psmye,psmyep)
      call disp(ndelt,nconfig,nitc,axmy,ax2my,axmya,axmye,axmyep)
      call disp(ndelt,nconfig,nitc,vemy,ve2my,vemya,vemye,vemyep)

c----------------------------------------------------------------------------c
c     spatial, temporal components                                           c
c----------------------------------------------------------------------------c

      call disp(ndelt,nconfig,nitc,avmy,av2my,avmya,avmye,avmyep)
      call disp(ndelt,nconfig,nitc,vvmy,vv2my,vvmya,vvmye,vvmyep)
      call disp(ndelt,nconfig,nitc,a4my,a42my,a4mya,a4mye,a4myep)
      call disp(ndelt,nconfig,nitc,v4my,v42my,v4mya,v4mye,v4myep)
      call disp(ndelt,nconfig,nitc,a1my,a12my,a1mya,a1mye,a1myep)
      call disp(ndelt,nconfig,nitc,v1my,v12my,v1mya,v1mye,v1myep)

c----------------------------------------------------------------------------c
c     unflavored mesons                                                      c
c----------------------------------------------------------------------------c

      call disp(ndelt,nconfig,nitc,sigmy,sig2my,sigmya,sigmye,sigmyep)
      call disp(ndelt,nconfig,nitc,etamy,eta2my,etamya,etamye,etamyep)
      call disp(ndelt,nconfig,nitc,f1my,f12my,f1mya,f1mye,f1myep)
      call disp(ndelt,nconfig,nitc,ommy,om2my,ommya,ommye,ommyep)

c----------------------------------------------------------------------------c
c     multi strange mesons                                                   c
c----------------------------------------------------------------------------c
c
c      call disp(ndelt,nconfig,nitc,sscmy,ssc2my,sscmya,sscmye,sscmyep)
c      call disp(ndelt,nconfig,nitc,spsmy,sps2my,spsmya,spsmye,spsmyep)
c      call disp(ndelt,nconfig,nitc,saxmy,sax2my,saxmya,saxmye,saxmyep)
c      call disp(ndelt,nconfig,nitc,svemy,sve2my,svemya,svemye,svemyep)
c
c---------------------------------------------------------------------------c
c     off diagonal correlators                                              c
c---------------------------------------------------------------------------c

      call disp(ndelt,nconfig,nitc,vsmy,vs2my,vsmya,vsmye,vsmyep)
      call disp(ndelt,nconfig,nitc,pamy,pa2my,pamya,pamye,pamyep)
      call disp(ndelt,nconfig,nitc,vtmy,vt2my,vtmya,vtmye,vtmyep)
      call disp(ndelt,nconfig,nitc,atmy,at2my,atmya,atmye,atmyep)

c---------------------------------------------------------------------------c
c     heavy-light mesons                                                    c
c---------------------------------------------------------------------------c

      call disp(ndelt,nconfig,nitc,kpmy,kp2my,kpmya,kpmye,kpmyep)
      call disp(ndelt,nconfig,nitc,kmmy,km2my,kmmya,kmmye,kmmyep)

c---------------------------------------------------------------------------c
c     diquarks, need additional loop over dirac structures                  c
c---------------------------------------------------------------------------c

      do 960 idiq=1,5
         do 961 i=1,ndelt
         do 961 j=1,nconfig
            xc(i,j) = cdiqmy(idiq,i,j)
            xc2(i,j)= cdiq2my(idiq,i,j)
 961     continue
         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
         do 962 i=1,ndelt
            cdiqmya(idiq,i) = xca(i)
            cdiqmye(idiq,i) = xce(i)
            cdiqmyep(idiq,i)= xcep(i)
 962     continue
 960  continue

c---------------------------------------------------------------------------c
c     multi strange diquarks                                                c
c---------------------------------------------------------------------------c
c
c      do 963 idiq=1,5
c         do 964 i=1,ndelt
c         do 964 j=1,nconfig
c            xc(i,j) = scdiqmy(idiq,i,j)
c            xc2(i,j)= scdiq2my(idiq,i,j)
c 964     continue
c         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
c         do 965 i=1,ndelt
c            scdiqmya(idiq,i) = xca(i)
c            scdiqmye(idiq,i) = xce(i)
c            scdiqmyep(idiq,i)= xcep(i)
c 965     continue
c 963  continue
c
c--------------------------------------------------------------------------c
c     nucleons                                                             c
c--------------------------------------------------------------------------c

      do 970 in=1,6
         do 971 i=1,ndelt
         do 971 j=1,nconfig
            xc(i,j) = clmy(in,i,j)
            xc2(i,j)= cl2my(in,i,j)
 971     continue
         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
         do 972 i=1,ndelt
            clmya(in,i) = xca(i)
            clmye(in,i) = xce(i)
            clmyep(in,i)= xcep(i)
 972     continue
 970  continue

c--------------------------------------------------------------------------c
c     lambdas                                                              c
c--------------------------------------------------------------------------c
c
c      do 1970 in=1,6
c         do 1971 i=1,ndelt
c         do 1971 j=1,nconfig
c            xc(i,j) = cclmy(in,i,j)
c            xc2(i,j)= ccl2my(in,i,j)
c1971     continue
c         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
c         do 1972 i=1,ndelt
c            cclmya(in,i) = xca(i)
c            cclmye(in,i) = xce(i)
c            cclmyep(in,i)= xcep(i)
c1972     continue
c1970  continue
c
c--------------------------------------------------------------------------c
c     more lambdas                                                         c
c--------------------------------------------------------------------------cc
c
c      do 2970 in=1,2
c         do 2971 i=1,ndelt
c         do 2971 j=1,nconfig
c            xc(i,j) = dclmy(in,i,j)
c            xc2(i,j)= dcl2my(in,i,j)
c2971     continue
c         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
c         do 2972 i=1,ndelt
c            dclmya(in,i) = xca(i)
c            dclmye(in,i) = xce(i)
c            dclmyep(in,i)= xcep(i)
c2972     continue
c2970  continue
c
c--------------------------------------------------------------------------c
c     xi                                                                   c
c--------------------------------------------------------------------------c
c
c      do 973 in=1,6
c         do 974 i=1,ndelt
c         do 974 j=1,nconfig
c            xc(i,j) = sclmy(in,i,j)
c            xc2(i,j)= scl2my(in,i,j)
c 974     continue
c         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
c         do 975 i=1,ndelt
c            sclmya(in,i) = xca(i)
c            sclmye(in,i) = xce(i)
c            sclmyep(in,i)= xcep(i)
c 975     continue
c 973  continue
c
c--------------------------------------------------------------------------c
c     delta                                                                c
c--------------------------------------------------------------------------c

      do 980 id=1,4
         do 981 i=1,ndelt
         do 981 j=1,nconfig
            xc(i,j) = cdlmy(id,i,j)
            xc2(i,j)= cdl2my(id,i,j)
 981     continue
         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
         do 982 i=1,ndelt
            cdlmya(id,i) = xca(i)
            cdlmye(id,i) = xce(i)
            cdlmyep(id,i)= xcep(i)
 982     continue
 980  continue

c--------------------------------------------------------------------------c
c     xi^*                                                                 c
c--------------------------------------------------------------------------c
c
c      do 983 id=1,4
c         do 984 i=1,ndelt
c         do 984 j=1,nconfig
c            xc(i,j) = scdlmy(id,i,j)
c            xc2(i,j)= scdl2my(id,i,j)
c 984     continue
c         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
c         do 985 i=1,ndelt
c            scdlmya(id,i) = xca(i)
c            scdlmye(id,i) = xce(i)
c            scdlmyep(id,i)= xcep(i)
c 985     continue
c 983  continue
c
c--------------------------------------------------------------------------c
c     omega                                                                c
c--------------------------------------------------------------------------c
c
c      do 986 id=1,4
c         do 987 i=1,ndelt
c         do 987 j=1,nconfig
c            xc(i,j) = sscdlmy(id,i,j)
c            xc2(i,j)= sscdl2my(id,i,j)
c 987     continue
c         call disp(ndelt,nconfig,nitc,xc,xc2,xca,xce,xcep)
c         do 988 i=1,ndelt
c            sscdlmya(id,i) = xca(i)
c            sscdlmye(id,i) = xce(i)
c            sscdlmyep(id,i)= xcep(i)
c 988     continue
c 986  continue
c
c--------------------------------------------------------------------------c
c     output, start with condensates                                       c
c--------------------------------------------------------------------------c

      write(2,*)
      write(2,*)
      write(2,1109) xtija,xtije
1109  format(1x, ' <\sum_j |T_ij|^2>      = ',f10.6,'+/-',f10.5)
      write(2,*)
      write(2,1110) xdeta,xdete
1110  format(1x, ' log(det_f(iD+im))      = ',f10.3,'+/-',f10.5)
      write(2,*)
      write(2,1111) qbarq , qbarqe
1111  format(1x, ' quark condensate  <qq> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1112) uu,uue
1112  format(1x, ' average local     <uu> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1113) uudd,uudde
1113  format(1x, ' 4 quark cond.   <uudd> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1114) uuuu,uuuue
1114  format(1x, ' 4 quark cond.   <uuuu> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1108) ua1,ua1e
1108  format(1x, ' t Hooft op.   <det qq> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)

c--------------------------------------------------------------------------c
c     colored 4-quark condensates                                          c
c--------------------------------------------------------------------------c

      write(2,*)
      write(2,1115) ss,sse
1115  format(1x, ' 4 quark cond. <uSddSu> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1116) pp,uue
1116  format(1x, ' 4 quark cond. <uPddPu> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1117) vv,vve
1117  format(1x, ' 4 quark cond. <uVddVu> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1118) aa,aae
1118  format(1x, ' 4 quark cond. <uAddAu> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)

c--------------------------------------------------------------------------c
c     weinberg operators                                                   c
c--------------------------------------------------------------------------c

      write(2,*)
      write(2,1119) o1,o1e
1119  format(1x, ' 4 quark cond.    <O_1> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1120) o2,o2e
1120  format(1x, ' 4 quark cond.    <O_2> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1121) o3,o3e
1121  format(1x, ' 4 quark cond.    <O_3> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)
      write(2,1122) o4,o4e
1122  format(1x, ' 4 quark cond.    <O_4> = ',f10.5,'+/-',f10.5,' MeV')
      write(2,*)

c--------------------------------------------------------------------------c
c     compare with free correlators                                        c
c--------------------------------------------------------------------------c

      cdfree(1) =-6.0/pi**4
      cdfree(2) = 6.0/pi**4
      cdfree(3) = 12./pi**4
      cdfree(4) = 12./pi**4
      cdfree(5) = cdfree(3)

      clfree(1) =-96.0/pi**6
      clfree(2) = clfree(1)
      clfree(3) =-576.0/pi**6
      clfree(4) = clfree(3)
      clfree(5) =-clfree(1)
      clfree(6) =-clfree(1)

      cdlfree(1) =-72.0/pi**6
      cdlfree(2) =-cdlfree(1)
      cdlfree(3) = 36.0/pi**6
      cdlfree(4) =-cdlfree(3)

c---------------------------------------------------------------------------c
c     equal r shuryak factor                                                c
c---------------------------------------------------------------------------c

      if(idir .eq. 4) then
       do 549 k=1,ndelt
         delfix = delt*k
         y = pi*delfix/alb(4)
         shu(k) = y**3/2.0*(1.0+cos(y)**2)/sin(y)**3
549    continue
      else
       do 548 k=1,ndelt
         shu(k) = 1.0
548    continue
      endif

c---------------------------------------------------------------------------c
c     mesonic correlators                                                   c
c---------------------------------------------------------------------------c

      write(2,550)
550   format(1x, '  scalar, Pi/Pi0')
      do 551 k=1,ndelt
         delfix = delt*k
         scfree = 3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,scmya(k)/scfree,scmye(k)/abs(scfree),
     1   scmyep(k)/abs(scfree)
551   continue

      write(2,552)
552   format(1x, '  pseudoscalar, Pi/Pi0 ')
      do 553 k=1,ndelt
         delfix = delt*k
         psfree =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,psmya(k)/psfree,psmye(k)/abs(psfree),
     1   psmyep(k)/abs(psfree)
553   continue

      write(2,554)
554   format(1x, '  axial vector, Pi/Pi0 ')
      do 555 k=1,ndelt
         delfix = delt*k
         axfree =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,axmya(k)/axfree,axmye(k)/abs(axfree),
     1   axmyep(k)/abs(axfree)
555   continue

      write(2,556)
556   format(1x, '  vector,  Pi/Pi0 ')
      do 557 k=1,ndelt
         delfix = delt*k
         vefree =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,vemya(k)/vefree,vemye(k)/abs(vefree),
     1   vemyep(k)/abs(vefree)
557   continue

c---------------------------------------------------------------------------c
c     spatial, temporal components                                          c
c---------------------------------------------------------------------------c

      write(2,9550)
9550  format(1x, ' axial vv, Pi/Pi0')
      do 9551 k=1,ndelt
         delfix = delt*k
         avfree =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,avmya(k)/avfree,avmye(k)/abs(avfree),
     1   avmyep(k)/abs(avfree)
9551  continue

      write(2,9552)
9552  format(1x, '  vector vv, Pi/Pi0 ')
      do 9553 k=1,ndelt
         delfix = delt*k
         vvfree =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,vvmya(k)/vvfree,vvmye(k)/abs(vvfree),
     1   vvmyep(k)/abs(vvfree)
9553  continue

      write(2,9554)
9554  format(1x, '  axial 44, Pi/Pi0 ')
      do 9555 k=1,ndelt
         delfix = delt*k
         a4free =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,a4mya(k)/a4free,a4mye(k)/abs(a4free),
     1   a4myep(k)/abs(a4free)
9555  continue

      write(2,9556)
9556  format(1x, '  vector 44,  Pi/Pi0 ')
      do 9557 k=1,ndelt
         delfix = delt*k
         v4free =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,v4mya(k)/v4free,v4mye(k)/abs(v4free),
     1   v4myep(k)/abs(v4free)
9557  continue

      write(2,9558)
9558  format(1x, '  axial 11, Pi/Pi0 ')
      do 9559 k=1,ndelt
         delfix = delt*k
         a1free =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,a1mya(k)/a1free,a1mye(k)/abs(a1free),
     1   a1myep(k)/abs(a1free)
9559  continue

      write(2,9560)
9560  format(1x, '  vector 11,  Pi/Pi0 ')
      do 9561 k=1,ndelt
         delfix = delt*k
         v1free =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,v1mya(k)/v1free,v1mye(k)/abs(v1free),
     1   v1myep(k)/abs(v1free)
9561  continue

c---------------------------------------------------------------------------c
c     unflavored mesons                                                     c
c---------------------------------------------------------------------------c

      write(2,1550)
1550  format(1x, '  sigma, Pi/Pi0')
      do 1551 k=1,ndelt
         delfix = delt*k
         sigfree = 3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,sigmya(k)/sigfree,sigmye(k)/abs(sigfree),
     1   sigmyep(k)/abs(sigfree)
1551  continue

      write(2,1552)
1552  format(1x, '  eta prime, Pi/Pi0 ')
      do 1553 k=1,ndelt
         delfix = delt*k
         etafree =-3.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,etamya(k)/etafree,etamye(k)/abs(etafree),
     1   etamyep(k)/abs(etafree)
1553  continue

      write(2,1554)
1554  format(1x, '  f1, Pi/Pi0')
      do 1555 k=1,ndelt
         delfix = delt*k
         f1free =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,f1mya(k)/f1free,f1mye(k)/abs(f1free),
     1   f1myep(k)/abs(f1free)
1555  continue

      write(2,1556)
1556  format(1x, '  omega, Pi/Pi0 ')
      do 1557 k=1,ndelt
         delfix = delt*k
         omfree =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,ommya(k)/omfree,ommye(k)/abs(omfree),
     1   ommyep(k)/abs(omfree)
1557  continue

c--------------------------------------------------------------------------c
c     multi strange mesons                                                 c
c--------------------------------------------------------------------------c
c
c      write(2,540)
c540   format(1x, ' s=2 scalar, Pi/Pi0')
c      do 541 k=1,ndelt
c         delfix = delt*k
c         scfree = 3.0/pi**4/delfix**6
c         write(2,444) delfix,sscmya(k)/scfree,sscmye(k)/abs(scfree),
c     1   sscmyep(k)/abs(scfree)
c541   continue
c
c      write(2,542)
c542   format(1x, ' s=2 pseudoscalar, Pi/Pi0 ')
c      do 543 k=1,ndelt
c         delfix = delt*k
c         psfree =-3.0/pi**4/delfix**6
c         write(2,444) delfix,spsmya(k)/psfree,spsmye(k)/abs(psfree),
c     1   spsmyep(k)/abs(psfree)
c543   continue
c
c      write(2,544)
c544   format(1x, ' s=2 axial vector, Pi/Pi0 ')
c      do 545 k=1,ndelt
c         delfix = delt*k
c         axfree =-6.0/pi**4/delfix**6
c         write(2,444) delfix,saxmya(k)/axfree,saxmye(k)/abs(axfree),
c     1   saxmyep(k)/abs(axfree)
c545   continue
c
c      write(2,546)
c546   format(1x, ' s=2 vector,  Pi/Pi0 ')
c      do 547 k=1,ndelt
c         delfix = delt*k
c         vefree =-6.0/pi**4/delfix**6
c         write(2,444) delfix,svemya(k)/vefree,svemye(k)/abs(vefree),
c     1   svemyep(k)/abs(vefree)
c547   continue
c
c-------------------------------------------------------------------------c
c     off diagonal correlators                                            c
c-------------------------------------------------------------------------c

      write(2,561)
561   format(1x, ' scalar-vector,  Pi/Piv ')
      do 562 k=1,ndelt
         delfix = delt*k
         vefree =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,vsmya(k)/vefree,vsmye(k)/abs(vefree),
     1   vsmyep(k)/abs(vefree)
562   continue

      write(2,563)
563   format(1x, ' axial-pseudoscalar,  Pi/Pia ')
      do 564 k=1,ndelt
         delfix = delt*k
         axfree =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,pamya(k)/axfree,pamye(k)/abs(axfree),
     1   pamyep(k)/abs(axfree)
564   continue

      write(2,565)
565   format(1x, ' tensor-vector,  Pi/Piv ')
      do 566 k=1,ndelt
         delfix = delt*k
         vefree =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,vtmya(k)/vefree,vtmye(k)/abs(vefree),
     1   vtmyep(k)/abs(vefree)
566   continue

      write(2,567)
567   format(1x, ' tensor-axialvector,  Pi/Pia ')
      do 568 k=1,ndelt
         delfix = delt*k
         axfree =-6.0/pi**4/delfix**6*shu(k)**2
         write(2,444) delfix,atmya(k)/axfree,atmye(k)/abs(axfree),
     1   atmyep(k)/abs(axfree)
568   continue

c----------------------------------------------------------------------------c
c     heavy-light mesons                                                     c
c----------------------------------------------------------------------------c

      write(2,530)
530   format(1x, ' k plus,  Pi/Pi0 ')
      do 531 k=1,ndelt
         delfix = delt*k
         xkfree =-3.0/pi**2/delfix**3*shu(k)
         write(2,444) delfix,kpmya(k)/xkfree,kpmye(k)/abs(xkfree),
     1   kpmyep(k)/abs(xkfree)
531   continue

      write(2,532)
532   format(1x, ' k minus,  Pi/Pi0 ')
      do 533 k=1,ndelt
         delfix = delt*k
         xkfree = 3.0/pi**2/delfix**3*shu(k)
         write(2,444) delfix,kmmya(k)/xkfree,kmmye(k)/abs(xkfree),
     1   kmmyep(k)/abs(xkfree)
533   continue

c----------------------------------------------------------------------------c
c     diquarks                                                               c
c----------------------------------------------------------------------------c

      do 570 idiq=1,5

      write(2,571) idiq
 571  format(1x, ' diquark cdiq(',i1,')/cfree')
      do 572 k=1,ndelt
         delfix = delt*k
         cd0 = cdfree(idiq)/delfix**6*shu(k)**2
         write(2,444) delfix,cdiqmya(idiq,k)/cd0,
     1   cdiqmye(idiq,k)/abs(cd0),cdiqmyep(idiq,k)/abs(cd0)
 572  continue

 570  continue

c----------------------------------------------------------------------------c
c     multi strange diquarks (only idiq=3,5 exist)                           c
c----------------------------------------------------------------------------c
c
c      do 575 idiq=3,5,2
c
c      write(2,576) idiq
c 576  format(1x, ' s=2 strange diquark cdiq(',i1,')/cfree')
c      do 577 k=1,ndelt
c         delfix = delt*k
c         cd0 = cdfree(idiq)/delfix**6
c         write(2,444) delfix,scdiqmya(idiq,k)/cd0,
c     1   scdiqmye(idiq,k)/abs(cd0),scdiqmyep(idiq,k)/abs(cd0)
c 577  continue
c
c 575  continue
c
c----------------------------------------------------------------------------c
c     nucleon                                                                c
c----------------------------------------------------------------------------c

      do 580 in=1,6

      write(2,581) in
 581  format(1x, ' nucleon cl(',i1,')/clfree')
      do 582 k=1,ndelt
         delfix = delt*k
         cl0 = clfree(in)/delfix**9*shu(k)**3
         write(2,444) delfix,clmya(in,k)/cl0,clmye(in,k)/abs(cl0),
     1   clmyep(in,k)/abs(cl0)
 582  continue

 580  continue

c----------------------------------------------------------------------------c
c     lambda                                                                 c
c----------------------------------------------------------------------------c
c
c      do 1580 in=1,6
c
c      write(2,1581) in
c1581  format(1x, ' lambda p/s cl(',i1,')/clfree')
c      do 1582 k=1,ndelt
c         delfix = delt*k
c         cl0 = clfree(in)/delfix**9
c         write(2,444) delfix,cclmya(in,k)/cl0,cclmye(in,k)/abs(cl0),
c     1   cclmyep(in,k)/abs(cl0)
c1582  continue
c
c1580  continue
c
c----------------------------------------------------------------------------c
c     more lambdas                                                           c
c----------------------------------------------------------------------------c
c
c      do 2580 in=1,2
c
c      write(2,2581) in
c2581  format(1x, ' lambda v cl(',i1,')/clfree')
c      do 2582 k=1,ndelt
c         delfix = delt*k
c         cl0 = clfree(in)/delfix**9
c         write(2,444) delfix,dclmya(in,k)/cl0,dclmye(in,k)/abs(cl0),
c     1   dclmyep(in,k)/abs(cl0)
c2582  continue
c
c2580  continue
c
c----------------------------------------------------------------------------c
c     xi                                                                     c
c----------------------------------------------------------------------------c
c
c      do 583 in=1,6
c
c      write(2,584) in
c 584  format(1x, ' xi cl(',i1,')/clfree')
c      do 585 k=1,ndelt
c         delfix = delt*k
c         cl0 = clfree(in)/delfix**9
c         write(2,444) delfix,sclmya(in,k)/cl0,sclmye(in,k)/abs(cl0),
c     1   sclmyep(in,k)/abs(cl0)
c 585  continue
c
c 583  continue
c
c----------------------------------------------------------------------------c
c     delta                                                                  c
c----------------------------------------------------------------------------c

      do 590 id=1,4

      write(2,591) id
 591  format(1x, ' delta cdl(',i1,')/cdlfree')
      do 592 k=1,ndelt
         delfix = delt*k
         cdl0 = cdlfree(id)/delfix**9*shu(k)**3
         write(2,444) delfix,cdlmya(id,k)/cdl0,cdlmye(id,k)/abs(cdl0),
     1   cdlmyep(id,k)/abs(cdl0)
 592  continue

 590  continue

c----------------------------------------------------------------------------c
c     xi^*                                                                   c
c----------------------------------------------------------------------------c
c
c      do 593 id=1,4
c
c      write(2,594) id
c 594  format(1x, ' xi^* cdl(',i1,')/cdlfree')
c      do 595 k=1,ndelt
c         delfix = delt*k
c         cdl0 = cdlfree(id)/delfix**9
c         write(2,444) delfix,scdlmya(id,k)/cdl0,
c     1   scdlmye(id,k)/abs(cdl0),scdlmyep(id,k)/abs(cdl0)
c 595  continue
c
c 593  continue
c
c----------------------------------------------------------------------------c
c     omega                                                                  c
c----------------------------------------------------------------------------c
c
c      do 596 id=1,4
c
c      write(2,597) id
c 597  format(1x, ' omega  cdl(',i1,')/cdlfree')
c      do 598 k=1,ndelt
c         delfix = delt*k
c         cdl0 = cdlfree(id)/delfix**9
c         write(2,444) delfix,sscdlmya(id,k)/cdl0,
c     1   sscdlmye(id,k)/abs(cdl0),sscdlmyep(id,k)/abs(cdl0)
c 598  continue
c
c 596  continue
c
c----------------------------------------------------------------------------c
c     trace of schwinger factor                                              c
c----------------------------------------------------------------------------c

      write(2,*) 'trace of pathexp'
      do 593 k=1,ndelt
         delfix = delt*k
         call mydisp(ndis(k),pemy(k),pe2my(k),pea,pee)
         write(2,333) delfix,pea,pee
 593  continue

c----------------------------------------------------------------------------c
c     plot eigenvalue distribution                                           c
c----------------------------------------------------------------------------c

      st   = (xmax-xmin)/float(nbins)
      call lev(xmin,st,nbins,1,ist)

      stop

  222 format(2(1x,f12.5))
  333 format(3(1x,f12.5))
  444 format(4(1x,f12.5))
 1021 format(i5,f8.3,f9.4,f8.4,f9.4,f8.4,f9.4,f8.4)
 1023 format(i5,f6.3,3(f12.4,f11.3))
 1022 format(1x, i5, a70)

      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine gammat
c----------------------------------------------------------------------------c
c     initialize euclidean gamma matrixes, spatial=-i*minkowski              c
c----------------------------------------------------------------------------c
      complex gamma, gamu5
      common /gam/   gamma(5,4,4), gr(4,4), gl(4,4), c(4,4)
      common /gamu5/ gamu5(4,4,4)
      dimension un(4,4)

      do 1  mu= 1,5
      do 1  m1= 1,4
      do 1  m2= 1,4
        gamma(mu,m1,m2)=cmplx(0.,0.)
        c(m1,m2) = 0.0
        gr(m1,m2) = 0.0
        gl(m1,m2) = 0.0
        un(m1,m2) = 0.0
1     continue

       gamma(1,1,4)=(0.,-1.)
       gamma(1,2,3)=(0.,-1.)
       gamma(1,3,2)=(0.,1.)
       gamma(1,4,1)=(0.,1.)
       gamma(2,1,4)=(-1.,0.)
       gamma(2,2,3)=(1.,0.)
       gamma(2,3,2)=(1.,0.)
       gamma(2,4,1)=(-1.,0.)
       gamma(3,1,3)=(0.,-1.)
       gamma(3,2,4)=(0.,1.)
       gamma(3,3,1)=(0.,1.)
       gamma(3,4,2)=(0.,-1.)
       gamma(4,1,1)=(1.,0.)
       gamma(4,2,2)=(1.,0.)
       gamma(4,3,3)=(-1.,0.)
       gamma(4,4,4)=(-1.,0.)

       gamma(5,1,3)=(-1.,0.)
       gamma(5,2,4)=(-1.,0.)
       gamma(5,3,1)=(-1.,0.)
       gamma(5,4,2)=(-1.,0.)

       c(1,4) = 1.0
       c(2,3) = -1.0
       c(3,2) = 1.0
       c(4,1) = -1.0

       un(1,1) = 1.0
       un(2,2) = 1.0
       un(3,3) = 1.0
       un(4,4) = 1.0

       do 10 i = 1, 4
       do 10 j = 1, 4
         gr(i,j) = (un(i,j)+real(gamma(5,i,j)))/2
         gl(i,j) = (un(i,j)-real(gamma(5,i,j)))/2
         do 20 ip = 1, 4
           gamu5(ip,i,j) = cmplx(0.0,0.0)
         do 20 k = 1, 4
           gamu5(ip,i,j) = gamu5(ip,i,j) + gamma(ip,i,k)*gamma(5,k,j)
   20    continue
   10  continue

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine gamfast
c----------------------------------------------------------------------------c
c     initialize arrays needed for fast gamma multiplication                 c
c     of euclidean gamma matrixes, spatial=-i*minkowski                      c
c----------------------------------------------------------------------------c
c     fast gamma trick:                                                      c
c     -----------------                                                      c
c     \sum_jk (\gam_a)_ij (A)_jk (\gam_b)_kl                                 c
c       = \sum_jk gamma(a,i,j)*a(j,k)*gamma(b,k,l)                           c
c       = g(a,i)*a(nu(a,i),ni(b,l))*g(b,ni(b,l))                             c
c----------------------------------------------------------------------------c

      complex g, gm5, guv
      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)
      common /gamf5/ gm5(4,4), nm5(4,4)
      common /guv/ guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)

c       gamma(1,1,4)=(0.,-1.)
c       gamma(1,2,3)=(0.,-1.)
c       gamma(1,3,2)=(0.,1.)
c       gamma(1,4,1)=(0.,1.)
       g(1,1) = (0.0,-1.0)
       nu(1,1) = 4
       ni(1,4) = 1
       g(1,2) = (0.0,-1.0)
       nu(1,2) = 3
       ni(1,3) = 2
       g(1,3) = (0.0,1.0)
       nu(1,3) = 2
       ni(1,2) = 3
       g(1,4) = (0.0,1.0)
       nu(1,4) = 1
       ni(1,1) = 4

c       gamma(2,1,4)=(-1.,0.)
c       gamma(2,2,3)=(1.,0.)
c       gamma(2,3,2)=(1.,0.)
c       gamma(2,4,1)=(-1.,0.)
       g(2,1) = (-1.0,0.0)
       nu(2,1) = 4
       ni(2,4) = 1
       g(2,2) = (1.0,0.0)
       nu(2,2) = 3
       ni(2,3) = 2
       g(2,3) = (1.0,0.0)
       nu(2,3) = 2
       ni(2,2) = 3
       g(2,4) = (-1.0,0.0)
       nu(2,4) = 1
       ni(2,1) = 4

c       gamma(3,1,3)=(0.,-1.)
c       gamma(3,2,4)=(0.,1.)
c       gamma(3,3,1)=(0.,1.)
c       gamma(3,4,2)=(0.,-1.)
       g(3,1) = (0.0,-1.0)
       nu(3,1) = 3
       ni(3,3) = 1
       g(3,2) = (0.0,1.0)
       nu(3,2) = 4
       ni(3,4) = 2
       g(3,3) = (0.0,1.0)
       nu(3,3) = 1
       ni(3,1) = 3
       g(3,4) = (0.0,-1.0)
       nu(3,4) = 2
       ni(3,2) = 4

c       gamma(4,1,1)=(1.,0.)
c       gamma(4,2,2)=(1.,0.)
c       gamma(4,3,3)=(-1.,0.)
c       gamma(4,4,4)=(-1.,0.)
       g(4,1) = (1.0,0.0)
       nu(4,1) = 1
       ni(4,1) = 1
       g(4,2) = (1.0,0.0)
       nu(4,2) = 2
       ni(4,2) = 2
       g(4,3) = (-1.0,0.0)
       nu(4,3) = 3
       ni(4,3) = 3
       g(4,4) = (-1.0,0.0)
       nu(4,4) = 4
       ni(4,4) = 4

c       gamma(5,1,3)=(-1.,0.)
c       gamma(5,2,4)=(-1.,0.)
c       gamma(5,3,1)=(-1.,0.)
c       gamma(5,4,2)=(-1.,0.)
       g(5,1) = (-1.0,0.0)
       nu(5,1) = 3
       ni(5,3) = 1
       g(5,2) = (-1.0,0.0)
       nu(5,2) = 4
       ni(5,4) = 2
       g(5,3) = (-1.0,0.0)
       nu(5,3) = 1
       ni(5,1) = 3
       g(5,4) = (-1.0,0.0)
       nu(5,4) = 2
       ni(5,2) = 4

c       c(1,4) = 1.0
c       c(2,3) = -1.0
c       c(3,2) = 1.0
c       c(4,1) = -1.0
       c(1) = 1.0
       nuc(1) = 4
       nic(4) = 1
       c(2) = -1.0
       nuc(2) = 3
       nic(3) = 2
       c(3) = 1.0
       nuc(3) = 2
       nic(2) = 3
       c(4) = -1.0
       nuc(4) = 1
       nic(1) = 4

      do 10 mu = 1, 4
      do 10 m = 1, 4
        gm5(mu,m) = g(mu,m)*g(5,nu(mu,m))
        nm5(mu,m) = nu(5,nu(mu,m))
   10 continue

      i = 0
      do 20 mu = 1, 4
      do 20 mv = mu+1, 4
        i = i+1
        mui(i) = mu
        mvi(i) = mv
        do 30 m = 1, 4
          guv(i,m) = g(mu,m)*g(mv,nu(mu,m))*cmplx(0.0,1.0)
          nuv(i,m) = nu(mv,nu(mu,m))
          nvu(i,m) = ni(mu,ni(mv,m))
   30   continue
   20 continue

      return
      end

c----------------------------------------------------------------------+-------
c----------------------------------------------------------------------+-------

      subroutine gamlr
c--------------------------------------------------------------------------c
c     initialize l/r projection of gamma matrices                          c
c--------------------------------------------------------------------------c
      complex gamplu,gammin,gamma
      common/lr/ gamplu(4,4,4),gammin(4,4,4)
      common /gam/ gamma(5,4,4), gr(4,4), gl(4,4), cccc(4,4)
          do 1 mu=1,4
          do 1 m1=1,4
          do 1 m2=1,4
           gamplu(mu,m1,m2)=0.5*gamma(mu,m1,m2)
           gammin(mu,m1,m2)=0.5*gamma(mu,m1,m2)
          do 1 m3=1,4
            gamplu(mu,m1,m2)= gamplu(mu,m1,m2)+
     *    gamma(mu,m1,m3)*gamma(5,m3,m2)*0.5
            gammin(mu,m1,m2)= gammin(mu,m1,m2)-
     *    gamma(mu,m1,m3)*gamma(5,m3,m2)*0.5
1     continue
      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine lammat
c-------------------------------------------------------------------------c
c     initialize Gell-Mann matrices and arrays for fast lambda algebra    c
c-------------------------------------------------------------------------c
c     basic fast lambda formula:                                          c
c     --------------------------                                          c
c     \sum_jk (\lam^a)_ij (A)_jk (\lam^b)_kl                              c
c           = \sum_jk lam(a,i,j)*a(j,k)*lam(b,k,l)                        c
c           = lf(a,i)*a(lu(a,i),li(b,l))*lf(b,li(b,l))                    c
c-------------------------------------------------------------------------c

      complex lambda,lf
      common /lam/  lambda(8,3,3)
      common /lamf/ lf(8,3), lu(8,3), li(8,3)

      do 10 k=1,8
      do 10 i=1,3
      do 10 j=1,3
         lambda(k,i,j) = (0.0,0.0)
 10   continue

      lambda(1,1,2) = (1.0,0.0)
      lambda(1,2,1) = (1.0,0.0)

      lambda(2,1,2) = (0.0,-1.)
      lambda(2,2,1) = (0.0,1.0)

      lambda(3,1,1) = (1.0,0.0)
      lambda(3,2,2) = (-1.,0.0)

      lambda(4,1,3) = (1.0,0.0)
      lambda(4,3,1) = (1.0,0.0)

      lambda(5,1,3) = (0.0,-1.)
      lambda(5,3,1) = (0.0,1.0)

      lambda(6,2,3) = (1.0,0.0)
      lambda(6,3,2) = (1.0,0.0)

      lambda(7,2,3) = (0.0,-1.)
      lambda(7,3,2) = (0.0,1.0)

      lambda(8,1,1) = (1.0,0.0)/sqrt(3.0)
      lambda(8,2,2) = (1.0,0.0)/sqrt(3.0)
      lambda(8,3,3) = (-2.,0.0)/sqrt(3.0)

c-------------------------------------------------------------------------c
c     lf(a,i) = non zero entry in i-th row of lambda^a                    c
c-------------------------------------------------------------------------c

      lf(1,1) = (1.0,0.0)
      lu(1,1) = 2
      li(1,2) = 1
      lf(1,2) = (1.0,0.0)
      lu(1,2) = 1
      li(1,1) = 2
      lf(1,3) = (0.0,0.0)
      lu(1,3) = 3
      li(1,3) = 3

      lf(2,1) = (0.0,-1.)
      lu(2,1) = 2
      li(2,2) = 1
      lf(2,2) = (0.0,1.0)
      lu(2,2) = 1
      li(2,1) = 2
      lf(2,3) = (0.0,0.0)
      lu(2,3) = 3
      li(2,3) = 3

      lf(3,1) = (1.0,0.0)
      lu(3,1) = 1
      li(3,1) = 1
      lf(3,2) = (-1.,0.0)
      lu(3,2) = 2
      li(3,2) = 2
      lf(3,3) = (0.0,0.0)
      lu(3,3) = 3
      li(3,3) = 3

      lf(4,1) = (1.0,0.0)
      lu(4,1) = 3
      li(4,3) = 1
      lf(4,2) = (0.0,0.0)
      lu(4,2) = 2
      li(4,2) = 2
      lf(4,3) = (1.0,0.0)
      lu(4,3) = 1
      li(4,1) = 3

      lf(5,1) = (0.0,-1.)
      lu(5,1) = 3
      li(5,3) = 1
      lf(5,2) = (0.0,0.0)
      lu(5,2) = 2
      li(5,2) = 2
      lf(5,3) = (0.0,1.0)
      lu(5,3) = 1
      li(5,1) = 3

      lf(6,1) = (0.0,0.0)
      lu(6,1) = 1
      li(6,1) = 1
      lf(6,2) = (1.0,0.0)
      lu(6,2) = 3
      li(6,3) = 2
      lf(6,3) = (1.0,0.0)
      lu(6,3) = 2
      li(6,2) = 3

      lf(7,1) = (0.0,0.0)
      lu(7,1) = 1
      li(7,1) = 1
      lf(7,2) = (0.0,-1.)
      lu(7,2) = 3
      li(7,3) = 2
      lf(7,3) = (0.0,1.0)
      lu(7,3) = 2
      li(7,2) = 3

      lf(8,1) = (1.0,0.0)/sqrt(3.0)
      lu(8,1) = 1
      li(8,1) = 1
      lf(8,2) = (1.0,0.0)/sqrt(3.0)
      lu(8,2) = 2
      li(8,2) = 2
      lf(8,3) = (-2.,0.0)/sqrt(3.0)
      lu(8,3) = 3
      li(8,3) = 3

      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine getx( point, delmax)
c---------------------------------------------------------------------------c
c     select first point for correlator measurement. make sure that         c
c     everything fits in the box.                                           c
c---------------------------------------------------------------------------c
c     input : delmax    maximum distance in time direction                  c
c     output: point(k)  coordinates of point                                c
c---------------------------------------------------------------------------c
      dimension point(4)
      common /box/   alb(4)
      common /direc/ idir

      do 10 k = 2, 4
         point(k) = alb(k)*rang()
   10 continue

      point(idir) = (alb(idir)-delmax)*rang()

      return
      end

c---------------------------------------------------------------------+-----
c---------------------------------------------------------------------+-----

      subroutine getsplit(point, xxsp, delfix)
c---------------------------------------------------------------------------c
c     select end points of propagators for meson correlators                c
c---------------------------------------------------------------------------c
c     point(4)   coordinates of starting point                              c
c     xxsp(4,2)  end points of propagator                                   c
c     delfix     distance in tau direction                                  c
c---------------------------------------------------------------------------c
      dimension xxsm(4,2), xxsp(4,2),point(4)
      common /box/   alb(4)
      common /direc/ idir
      do 10 k = 1, 4
         r = point(k)
         xxsp(k,1) = r
         xxsp(k,2) = r
   10 continue

      xxsp(idir,2) = point(idir) + delfix

      return
      end

c----------------------------------------------------------------------+-----
c----------------------------------------------------------------------+-----

      subroutine su3(n,nin,nd,x,y,u)
c----------------------------------------------------------------------------c
c     determine third row of su(3) matrix if the first two are given.        c
c----------------------------------------------------------------------------c
c     input :  n      number of colors                                       c
c              nin    number of instantons                                   c
c              nd     maximum number of instantons                           c
c              x(i,j) real,imaginary parts in first row of matrix u(i)       c
c              y(i,j) same for second row                                    c
c     output:  z(i,j) real,imaginary parts of third row                      c
c----------------------------------------------------------------------------c
      parameter(n2=256)
      complex u
      dimension x(nd,6), y(nd,6), z(n2,6), u(nd,n,n)

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
20    continue
      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine setup(n,nin,nmol,pfrac,zr,e1,e2,iread,ic)
c--------------------------------------------------------------------------c
c     generate a new configuration                                         c
c--------------------------------------------------------------------------c
c     this is the most general version of the cocktail setup. this         c
c     version initializes nin instantons, nin-2*nmol are random, the re-   c
c     maining 2*nmol are molecules. a fraction pfrac of the molecules is   c
c     oriented in the time direction.                                      c
c--------------------------------------------------------------------------c
c     input : n       number of colors                                     c
c             nin     number of instantons                                 c
c             nmol    number of molecules                                  c
c             iread   input from file/random configuration (=1/0)          c
c     output: zr(i,5) position,size of instanton i                         c
c             e1(i,6) first column of rotation matrix for instanton i      c
c             e2(i,6) second column   --- " ---                            c
c             icon    number of configuration                              c
c--------------------------------------------------------------------------c

      parameter( nd=256 )

      complex u0, tau
      complex u(nd,3,3), d(3,3)
      character*9 a8
      dimension zr(nd,5), er1(6), er2(6)
      dimension x(6), y(6), z(6), e1(nd,6), e2(nd,6)
      dimension falb(4), uv(4)
      common /param/ a,alpha,rh0, sg, dz, drh, nc, nf, rmu, rms
      common /counter/ icount
      common /box/ alb(4)
      common /nconf/ nconfig
      common /taum/ tau(2,4,2,2)

      if (iread.eq.1) then

c------------------------------------------------------------------------c
c     input from file                                                    c
c------------------------------------------------------------------------c

      if (ic .eq. 1) then

c------------------------------------------------------------------------c
c     initialize                                                         c
c------------------------------------------------------------------------c

         open(unit=4,file='infile.dat',status='old')
         read(4,*) inc, inf, inin, (falb(k),k=1,4), frho
         read(4,*) fsg, fdx1, fdx2, fdx3, fdx4, fdrho, falpha
         read(4,*) frmu, frms, fnconf

         if(inin .ne. nin) stop
         if(inc  .ne. nc)  stop
         if(fnconf .lt. nconfig) stop
         do 5 j=1,4
         if(falb(j) .ne. alb(j)) print 997
 5       continue
         if(inf .ne. nf) print 998
         if(frmu .ne. rmu .or. frms .ne. rms) print 999

      end if

 997  format(' Warning: box size has changed')
 998  format(' Warning: nf has changed')
 999  format(' Warning: masses have changed')

c------------------------------------------------------------------------c
c     read new configuration                                             c
c------------------------------------------------------------------------c

      read (4,*) icon
      write(6,208) icon
      write(2,208) icon
      do 25 i = 1, nin
         read(4,101) (zr(i,k),k=1,5)
         read(4,102) (e1(i,k),k=1,6)
         read(4,102) (e2(i,k),k=1,6)
   25 continue

  208 format(1x, ' configuration ', i5, '   has been read ',/)
  101 format(1x,5f12.5)
  102 format(1x,6f12.5)

c--------------------------------------------------------------------------c
c     iread = 0 : generate random configuration (but fixed size)           c
c--------------------------------------------------------------------------c

      else

c--------------------------------------------------------------------------c
c     optional: write header to outfile                                    c
c--------------------------------------------------------------------------c

      if ( ic .eq. 1 ) then
         open(unit=14,file='outfile.dat',status='unknown')
         write(14,204) nc,nf,nin,alb,rh0
         write(14,206) 1.0,1.0,1.0,1.0,1.0,1.0,1.0
         write(14,207) rmu,rms,nconfig
      endif

      icount = icount + 1
      write(6,701) icount
      write(2,701) icount
  701 format(1x,/,1x,'  Configuration: ', i5)
  204 format(1x,3i5,5f10.4)
  206 format(1x,7f10.4)
  207 format(1x,2f10.4,1x,i5)

c--------------------------------------------------------------------------c
c     loop over instantons                                                 c
c--------------------------------------------------------------------------c

      nin2 = nin/2
      do 100 ip = 1, nin
         do 10 i = 1, 4
          zr(ip,i) = alb(i)*rang( )
 10      continue
         zr(ip,5) = rh0
         do 20 i = 1, 6
            e1(ip,i) = 0.0
            e2(ip,i) = 0.0
 20      continue
         e1(ip,1) = 1.0
         e2(ip,3) = 1.0
         call rsu(6,er1,er2)
         do 70 i = 1, 6
            e1(ip,i) = er1(i)
            e2(ip,i) = er2(i)
 70      continue
 100  continue

c--------------------------------------------------------------------------c
c     replace first nmol instantons by molecules; reconstruct orientation  c
c--------------------------------------------------------------------------c

      call su3(3,nin,nd,e1,e2,u)

c--------------------------------------------------------------------------c
c     instanton coordinates unchanged, antiinstanton shifted + rotated     c
c--------------------------------------------------------------------------c

      drho = 2.00*rh0

      do 200 ii=1,nmol
         ia = ii + nin2

c--------------------------------------------------------------------------c
c     polarized molecules                                                  c
c--------------------------------------------------------------------------c

         if(ii .lt. pfrac*nmol) then
            zr(ia,1) = zr(ii,1)
            zr(ia,2) = zr(ii,2)
            zr(ia,3) = zr(ii,3)
            zr(ia,4) = zr(ii,4)+drho*sign(1.0,rang()-0.5)
            uv(1) = 0.0
            uv(2) = 0.0
            uv(3) = 0.0
            uv(4) = (zr(ia,4)-zr(ii,4))/drho
            if(zr(ia,4) .gt. alb(4)) zr(ia,4)=zr(ia,4)-alb(4)
            if(zr(ia,4) .lt. alb(4)) zr(ia,4)=zr(ia,4)+alb(4)
         else

c--------------------------------------------------------------------------c
c     random relative orientation                                          c
c--------------------------------------------------------------------------c

 215        continue

            call randuv(uv,4)
            do 210 m=1,4
               zr(ia,m) = zr(ii,m) + drho*uv(m)
 210        continue
            if (zr(ia,1) .gt. alb(1) .or. zr(ia,1) .lt. 0) goto 215
            if (zr(ia,2) .gt. alb(2) .or. zr(ia,2) .lt. 0) goto 215
            if (zr(ia,3) .gt. alb(3) .or. zr(ia,3) .lt. 0) goto 215
            if (zr(ia,4) .gt. alb(4) .or. zr(ia,4) .lt. 0) goto 215

         endif

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

c--------------------------------------------------------------------------c
c     optional: save configuration to outfile                              c
c--------------------------------------------------------------------------c

      write(14,*) ic
      do 300 j=1,nin
         write(14,101) (zr(j,k),k=1,5)
         write(14,102) (e1(j,k),k=1,6)
         write(14,102) (e2(j,k),k=1,6)
 300  continue

c--------------------------------------------------------------------------c
c     end of random + molecules part                                       c
c--------------------------------------------------------------------------c

      end if

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine rsu(n, x, y)
c----------------------------------------------------------------------------c
c     generate random su(3) matrix                                           c
c----------------------------------------------------------------------------c
c     input : n    twice the number of colors                                c
c     output: x(i) real,imaginary parts of first row (i=1,...,2n)            c
c             y(i) same for second row                                       c
c----------------------------------------------------------------------------c

      dimension x(n), y(n), xs(6)

      call randuv(x, n)

      xs(1) = x(2)
      xs(2) = - x(1)
      xs(3) = x(4)
      xs(4) = -x(3)
      xs(5) = x(6)
      xs(6) = -x(5)

      call randuv(y,n)

      call rdot(n, x, y, xdy)
      call rdot(n, xs, y, xsdy)

      r2 = 0.0
      do 10 i = 1, 6
        y(i) = y(i) - xdy*x(i) -xsdy*xs(i)
        r2 = r2 + y(i)*y(i)
   10 continue
      r = sqrt(r2)
      do 20 i = 1, 6
        y(i) = y(i)/r
   20 continue

      return
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      subroutine rdot(n, x, y, xdy)

      dimension x(n), y(n)

      xdy = 0.0
      do 10 i = 1, n
        xdy = xdy + x(i)*y(i)
   10 continue

      return
      end

c----------------------------------------------------------------------+--------
c----------------------------------------------------------------------+--------

      subroutine randuv(x,n)
c-----------------------------------------------------------------------------c
c     generate random n-dimensional unit vector                               c
c-----------------------------------------------------------------------------c

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

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      function rang()
c------------------------------------------------------------------------c
c     various random number generators                                   c
c------------------------------------------------------------------------c
c     rnunf (imsl), ranf (cray), ran2 (numerical recipes), ran (vax)     c
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
      ran2 = iy*rm
      idum = mod(ia*idum+ic,m)
      ir(j)= idum

      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine myaddto(x,xtot,x2tot)
c---------------------------------------------------------------------------c
c     include result x in counters xtot and x2tot                           c
c---------------------------------------------------------------------------c
      xtot=xtot+x
      x2tot=x2tot+x*x
      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine mydisp(n,xtot,x2tot,xav,xerr)
c---------------------------------------------------------------------------c
c     estimate average and error from xtot and x2tot                        c
c---------------------------------------------------------------------------c
c     input : n     number of measurements                                  c
c             xtot  sum of x_i                                              c
c             x2tot sum of x**2                                             c
c     output: xav   average                                                 c
c             xerr  error estimate                                          c
c---------------------------------------------------------------------------c
      if(n.lt.1)goto 10
      xav=xtot/(n*1.)
      del2=x2tot/(1.*n*n)-xav*xav/(1.*n)
      if(del2.lt.0.)del2=0.
      xerr=sqrt(del2)
10    continue
      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine disp(ndelt,nconfig,nitc,x,x2,xa,xe,xep)
c--------------------------------------------------------------------------c
c     calculate ensemble averages and error estimates.                     c
c--------------------------------------------------------------------------c
c     ndelt        number of separations                                   c
c     nconfig      number of configurations                                c
c     nitc         number of points in each configuration                  c
c     x(k,ic)      sum of measurements in config ic                        c
c     x2(k,ic)     sum of squares                                          c
c     xa(k)        ensemble average                                        c
c     xe(k)        error estimate                                          c
c     xep(k)       improved error estimate                                 c
c--------------------------------------------------------------------------c
      parameter( ncf=100 )
      parameter( ndl=20  )
      real x(ndl,ncf),x2(ndl,ncf)
      real xa(ndl),xe(ndl),xep(ndl)
c--------------------------------------------------------------------------c
c     first method : all measurements are uncorrelated                     c
c--------------------------------------------------------------------------c

      do 20 k=1,ndelt
         xtot = 0.0
         x2tot= 0.0
         do 30 ic=1, nconfig
               xtot = xtot + x(k,ic)
               x2tot= x2tot+ x2(k,ic)
 30      continue
         n     = nconfig*nitc
         xav   = xtot/float(n)
         xvar  = abs( x2tot-n*xav**2 )/float(n)
         xa(k) = xav
         xe(k) = sqrt(xvar/float(n))
 20   continue

c----------------------------------------------------------------------------c
c     second method : independent error estimates for each point and config. c
c----------------------------------------------------------------------------c

      do 60 k=1,ndelt
         xtot = 0.0
         x2tot= 0.0
         dxtot= 0.0
         do 70 ic=1,nconfig

c----------------------------------------------------------------------------c
c     average and error for given config                                     c
c----------------------------------------------------------------------------c

            xav   = x(k,ic)/float(nitc)
            xvar  = abs( x2(k,ic)-nitc*xav**2 )/float(nitc)
            xerr  = sqrt(xvar/float(nitc))
            xtot  = xtot + xav
            x2tot = x2tot+ xav**2
            dxtot = dxtot+ xerr
 70      continue

c----------------------------------------------------------------------------c
c     ensemble average and error                                             c
c----------------------------------------------------------------------------c

         xav    = xtot/float(nconfig)
         xvar   = abs(x2tot - nconfig*xav**2)/float(nconfig)
         xep(k) = sqrt(xvar/float(nconfig)) + dxtot/float(nconfig)
 60   continue
      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine zero(n,ria)
c-------------------------------------------------------------------------c
c     clear array ria(n)                                                  c
c-------------------------------------------------------------------------c
      dimension ria(n)
      do 10 i = 1, n
        ria(i) = 0.0
   10 continue
      return
      end

      subroutine izero(n,narr)
c-------------------------------------------------------------------------c
c     clear integer array narr(n)                                         c
c-------------------------------------------------------------------------c
      dimension narr(n)
      do 1 k=1,n
         narr(k)=0.
1     continue
      return
      end

      subroutine zero1(n,arr)
c-------------------------------------------------------------------------c
c     clear array arr(n)                                                  c
c-------------------------------------------------------------------------c
      dimension arr(n)
      do 1 k=1,n
         arr(k)=0.
1     continue
      return
      end

      subroutine zero2(n1,n2,arr)
c-------------------------------------------------------------------------c
c     clear array arr(n1,n2)                                              c
c-------------------------------------------------------------------------c
      dimension arr(n1,n2)
      do 1 j=1,n1
      do 1 k=1,n2
         arr(j,k)=0.
1     continue
      return
      end

      subroutine zero3(n1,n2,n3,arr)
c-------------------------------------------------------------------------c
c     clear arr(n1,n2,n3)                                                 c
c-------------------------------------------------------------------------c
      dimension arr(n1,n2,n3)
      do 10 i=1,n1
      do 10 j=1,n2
      do 10 k=1,n3
         arr(i,j,k) = 0.0
 10   continue
      return
      end

      subroutine zero4(n1,n2,n3,n4,arr)
c-------------------------------------------------------------------------c
c     clear arr(n1,n2,n3)                                                 c
c-------------------------------------------------------------------------c
      dimension arr(n1,n2,n3,n4)
      do 10 i=1,n1
      do 10 j=1,n2
      do 10 k=1,n3
      do 10 l=1,n4
         arr(i,j,k,l) = 0.0
 10   continue
      return
      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine pfree(xx,xmass,sfree)
c----------------------------------------------------------------------c
c     calculate free or vacuum dominance propagator from point xx(.,1) c
c     to xx(.,2). adjust coef,coef2 as needed.                         c
c----------------------------------------------------------------------c
c     input : xx(k,index)    endpoints; k=1,..,4 coordinates,          c
c                            index=1,2 for beginning, end point        c
c     output: sfree(a,i,b,j) free (vac dom) propagator; a,b=1,..3      c
c                            color, i,j=1,..,4 dirac        .          c
c----------------------------------------------------------------------c

      implicit real (a-z)
      real xx(4,2), dx(4)
      real gr(4,4), gl(4,4), c(4,4), u(4,4)
      complex sfree(3,4,3,4), gamma(5,4,4), coef, coef2
      integer a,b,i,j,k

      common /gam/ gamma,gr,gl,c

      do 10 a=1,3
      do 10 b=1,3
         do 10 i=1,4
         do 10 j=1,4
            sfree(a,i,b,j) = (0.0,0.0)
  10  continue
      do 15 i=1,4
      do 15 j=1,4
         U(i,j) = 0.0
  15  continue
      do 16 i=1,4
         U(i,i) = 1.0
  16  continue

c---------------------------------------------------------------------c
c     distance between end points                                     c
c---------------------------------------------------------------------c

      tau2 = 0.0
      do 20 k=1,4
         dx(k)= xx(k,2)-xx(k,1)
         tau2 = tau2 + dx(k)**2
  20  continue
      tau  = sqrt(tau2)
      Pi   = 3.1415926
      coef = (0.0,1.0)/(2.0*Pi**2*tau**4)
      coef2= (0.0,1.0)/12.0*(250.0/197.3)**3 * 0.0
     1      +(0.0,1.0)*xmass/(4.0*Pi**2*tau**2) * 0.0

c----------------------------------------------------------------------c
c     only color diagonal part is non zero                             c
c----------------------------------------------------------------------c

      do 30 i=1,4
      do 30 j=1,4
         do 30 a=1,3
            sfree(a,i,a,j) = coef * ( gamma(1,i,j)*dx(1)
     2        + gamma(2,i,j)*dx(2) + gamma(3,i,j)*dx(3)
     3        + gamma(4,i,j)*dx(4) ) + coef2 * U(i,j)
  30  continue
      return
      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine bin(n,nbin,val,ist,xmin,xmax,ix)
c-----------------------------------------------------------------------c
c     bin data set val(n) into nbin segments                            c
c-----------------------------------------------------------------------c
c     input : n         number of events                                c
c             nbin      number of bins                                  c
c             val(n)    events val(1)...val(n)                          c
c             ist(nbin) histo array                                     c
c             ix.eq.1   xmin,xmax determined from data (output)         c
c             ix.ne.1   xmin,xmax given on input                        c
c-----------------------------------------------------------------------c
      real val(n)
      integer ist(nbin)
      if (ix .eq. 1) then
         xmin = val(1)
         xmax = xmin
         do 10 i=1,n
            xmin = min(xmin,val(i))
            xmax = max(xmax,val(i))
 10      continue
      endif
      st = (xmax-xmin)/float(nbin)
      do 15 j=1,nbin
         ist(j) = 0
 15   continue
      do 20 i=1,n
         j = int( (val(i)-xmin)/st ) + 1
         if (j .lt. 1)    j=1
         if (j .gt. nbin) j=nbin
         ist(j) = ist(j) + 1
 20   continue
      end


c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine addbin(n,nbin,val,xmin,xmax,ist)
c-----------------------------------------------------------------------c
c     include data set val(n) in histogram array ist(nbin)              c
c-----------------------------------------------------------------------c
c     input : n         number of events                                c
c             nbin      number of bins                                  c
c             val(n)    events val(1)...val(n)                          c
c             ist(nbin) histo array                                     c
c             xmin      lower limit of histogram                        c
c             xmax      upper bound                                     c
c-----------------------------------------------------------------------c
      real val(n)
      integer ist(nbin)
      st = (xmax-xmin)/float(nbin)
      do 20 i=1,n
         j = int( (val(i)-xmin)/st ) + 1
         if (j .lt. 1)    j=1
         if (j .gt. nbin) j=nbin
         ist(j) = ist(j) + 1
 20   continue
      end

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine lev(xmin,st,n,nsm,ist)
c-----------------------------------------------------------------------c
c     plot simple histogram in output file ftn02                        c
c-----------------------------------------------------------------------c
c     input:  xmin    smallest x value                                  c
c             st      bin width                                         c
c             n       number of bins                                    c
c             nsm     plot symbol                                       c
c             ist(n)  histogram array                                   c
c-----------------------------------------------------------------------c
      dimension ist(n),a(6),k(6),t(50)
      data a/1hx,1h*,1h+,1hc,1he,1h /
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

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine scalcor(nin,deltas,fak)
c----------------------------------------------------------------------------c
c     calculate screening factor fak(I)=D_M(x,z_I)/D_0(x,z_I)                c
c----------------------------------------------------------------------------c
c     nin       number of instantons                                         c
c     deltas(i) distance (x-z_i)^2                                           c
c     fak(i)    screening factor                                             c
c----------------------------------------------------------------------------c
c     screening mass stored in common scalmass, not flavor dependent !       c
c----------------------------------------------------------------------------c

      common /scalmass/ amass
      dimension deltas(nin), fak(nin)

      do 10 is = 1, nin
         t = sqrt(deltas(is))
         e = amass
         x = e*t
c------------------------------------------------------------------------c
c     imsl                                                               c
c------------------------------------------------------------------------c
c        dprop=e/t*bsk1(x)
c------------------------------------------------------------------------c
c     cern                                                               c
c------------------------------------------------------------------------c
         dprop=e/t*besk1(x)
         fak(is)=dprop*deltas(is)
   10 continue

      return
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      subroutine masscorr(amass,x,corg,corm)
c-------------------------------------------------------------------------c
c     mass corrections for free propagators                               c
c-------------------------------------------------------------------------c
c     amass  current quark mass                                           c
c     corg   S_m(x)/S_0(x), gamma matrix part                             c
c     corm   D_m(x)/D_0(x)                                                c
c-------------------------------------------------------------------------c
      real k(0:2)

      arg = amass*x
      if ( arg .gt. 86.0) then
         k(0) = 0.0
         k(1) = 0.0
         k(2) = 0.0
      else
c------------------------------------------------------------------------c
c     imsl                                                               c
c------------------------------------------------------------------------c
c        call bsks(0.0,arg,3,k)
c------------------------------------------------------------------------c
c     cern                                                               c
c------------------------------------------------------------------------c
         k(0) = besk0(arg)
         k(1) = besk1(arg)
         k(2) = k(0)+2.0/arg*k(1)
      endif
      corm = arg*k(1)
      corg = arg**2*k(2)/2.0

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine mult2(nin,a,b,c)
c---------------------------------------------------------------------------c
c     multiply nin 2*2 matrices a(i,2,2)*b(i,2,2)-->c(i,2,2)                c
c---------------------------------------------------------------------------c

      parameter(ni=256)
      complex a(ni,2,2), b(ni,2,2),c(ni,2,2)

      do 10 is = 1, nin
         c(is,1,1) = a(is,1,2)*b(is,2,1)+a(is,1,1)*b(is,1,1)
         c(is,1,2) = a(is,1,2)*b(is,2,2)+a(is,1,1)*b(is,1,2)
         c(is,2,1) = a(is,2,2)*b(is,2,1)+a(is,2,1)*b(is,1,1)
         c(is,2,2) = a(is,2,2)*b(is,2,2)+a(is,2,1)*b(is,1,2)
   10 continue

      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine mult3(nin,a,b,c)
c---------------------------------------------------------------------------c
c     multiply nin 3*3 matrices a(i,3,3)*b(i,3,3)-->c(i,3,3)                c
c---------------------------------------------------------------------------c

      parameter (ni=256)
      complex a(ni,3,3),b(ni,3,3),c(ni,3,3)

      do 10 is = 1, nin
           c(is,1,1) = a(is,1,1)*b(is,1,1)+a(is,1,2)*b(is,2,1)+
     2                 a(is,1,3)*b(is,3,1)
           c(is,1,2) = a(is,1,1)*b(is,1,2)+a(is,1,2)*b(is,2,2)+
     2                 a(is,1,3)*b(is,3,2)
           c(is,1,3) = a(is,1,1)*b(is,1,3)+a(is,1,2)*b(is,2,3)+
     2                 a(is,1,3)*b(is,3,3)
           c(is,2,1) = a(is,2,1)*b(is,1,1)+a(is,2,2)*b(is,2,1)+
     2                 a(is,2,3)*b(is,3,1)
           c(is,2,2) = a(is,2,1)*b(is,1,2)+a(is,2,2)*b(is,2,2)+
     2                 a(is,2,3)*b(is,3,2)
           c(is,2,3) = a(is,2,1)*b(is,1,3)+a(is,2,2)*b(is,2,3)+
     2                 a(is,2,3)*b(is,3,3)
           c(is,3,1) = a(is,3,1)*b(is,1,1)+a(is,3,2)*b(is,2,1)+
     2                 a(is,3,3)*b(is,3,1)
           c(is,3,2) = a(is,3,1)*b(is,1,2)+a(is,3,2)*b(is,2,2)+
     2                 a(is,3,3)*b(is,3,2)
           c(is,3,3) = a(is,3,1)*b(is,1,3)+a(is,3,2)*b(is,2,3)+
     2                 a(is,3,3)*b(is,3,3)
   10 continue

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine taumat
c----------------------------------------------------------------------------c
c     initilize four-dim. tau matrices, the first number is index            c
c----------------------------------------------------------------------------c
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

      subroutine sigvec(nin,x,ind,s)
c----------------------------------------------------------------------------c
c     putting 4-vector into 2*2 matrix, sig(index,mu)*x(mu), index is + or -.c
c     sweep over all instantons, index is reversed for second half.          c
c----------------------------------------------------------------------------c
      parameter (ni=256)

      complex s(ni,2,2)
      dimension x(ni,4)

      nih = nin/2
      do 10 is = 1, nin
         indmax = (is-1)/nih + 1
         index = (3-2*indmax)*ind
c----------------------------------------------------------------------------c
c     \sig_\mu^(+) = (\vec\sig,-i) !!                                        c
c----------------------------------------------------------------------------c
         s(is,1,1)= index*x(is,4)*(0.,-1.)+x(is,3)
         s(is,2,2)= index*x(is,4)*(0.,-1.)-x(is,3)
         s(is,1,2)= x(is,1) - x(is,2)*(0.,1.)
         s(is,2,1)= x(is,1) + x(is,2)*(0.,1.)
   10 continue

      return
      end

c----------------------------------------------------------------------+-----
c----------------------------------------------------------------------+-----

      subroutine tordis(nin,x,y,z,rr)
c----------------------------------------------------------------------------c
c     determine shortest path from point x(k) to all instnatons in           c
c     box with periodic boundary conditions.                                 c
c----------------------------------------------------------------------------c
c     input : nin     number of instantons                                   c
c             x(k)    coordinates of point x                                 c
c             y(i,k)  coordinates of i-th instanton                          c
c                     note that y(i,5) is instanton size                     c
c     output: z(i,k)  vector of shortest path to i                           c
c             rr(i)   distance squared to instanton i                        c
c----------------------------------------------------------------------------c
      parameter (n2=256)
      common /box/ alb(4)
      dimension x(4),y(n2,5),z(n2,4), rr(n2)

      do 5 is = 1, nin
        rr(is) = 0.
   5  continue
       do 1 m=1,4
         do 10 is = 1, nin
           dis=x(m)-y(is,m)
           adis=abs(dis)
           asg = 1+sign(1.0,adis-alb(m)/2)
           dis = dis - alb(m)/2*asg*sign(1.0,dis)
           rr(is) = rr(is) + dis**2
           z(is,m)=dis
  10    continue
   1  continue

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine prmeson(sxy,syx,pexp,prsc,prps,praxt,prvect,
     1                   prav,prvv,pra4,prv4,pra1,prv1)
c---------------------------------------------------------------------------c
c     calculate meson bs amplitudes (including pexp)                        c
c---------------------------------------------------------------------------c
c     input : sxy(a,i,b,j)  first propagator (a,b=1,..,3;i,j=1,..,4)        c
c             syx(a,i,b,j)  second propagator                               c
c             pexp(a,b)     pathexp (a,b=1,..,3)                            c
c             prsc          scalar (isovector) bs amplitude                 c
c             prps          peudoscalar bs amplitude                        c
c             praxt         axialvector bs amplitude (traced)               c
c             prvect        vector bs amplitude (traced)                    c
c             prav,pra4     axial vector, spatial and temporal part         c
c             prvv,prv4     vector, spatial and temporal part               c
c---------------------------------------------------------------------------c

      parameter(nc=3)
      complex pexp(3,3)
      complex syxexp(3,4,3,4)
      dimension prax(4,4), prvec(4,4)
      complex sgg(nc,4,nc,4,4), sgg5(nc,4,nc,4,4)
      complex sxy(3,4,3,4), syx(3,4,3,4), s5(3,4,3,4)
      complex gamma, gamu5
      complex cprax(4,4), cprvec(4,4), cprsc, cprps
      complex g
      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)
      common /direc/ idir

c---------------------------------------------------------------------------c
c     multiply pexp by second propagator -->  syxexp(a,i;b,j)               c
c---------------------------------------------------------------------------c

      do 220 m1=1,4
      do 220 m2=1,4
      do 220 k1=1,3
      do 220 k2=1,3
         syxexp(k1,m1,k2,m2)=(0.,0.)
         do 221 kk=1,3
            syxexp(k1,m1,k2,m2)=syxexp(k1,m1,k2,m2)
     2                  + pexp(kk,k2)*syx(k1,m1,kk,m2)
221      continue
220   continue

c---------------------------------------------------------------------------c
c     calculate (gamma*sxy*gamma) for gamma=(s,p,v,a)                       c
c---------------------------------------------------------------------------c

      do 10 i = 1, 3
      do 10 j = 1, 3
        do 20 m = 1, 4
        do 20 n = 1, 4

c--------------------------------------------------------------------------c
c     for vectors do fixed index ip (no sum yet)                           c
c--------------------------------------------------------------------------c

          do 130 ip = 1, 4
          sgg(i,m,j,n,ip) = g(ip,m)*sxy(i,nu(ip,m),j,ni(ip,n))
     1                     *g(ip,ni(ip,n))
          sgg5(i,m,j,n,ip) = g(ip,m)*g(5,nu(ip,m))
     1              *sxy(i,nu(5,nu(ip,m)),j,ni(ip,ni(5,n)))
     2                       *g(ip,ni(ip,ni(5,n)))*g(5,ni(5,n))
 130      continue

c--------------------------------------------------------------------------c
c     pseudoscalar, scalar is just sxy                                     c
c--------------------------------------------------------------------------c

          s5(i,m,j,n) = g(5,m)*sxy(i,nu(5,m),j,ni(5,n))*g(5,ni(5,n))
  20    continue
  10  continue

c--------------------------------------------------------------------------c
c     contract vectors with second propagator                              c
c--------------------------------------------------------------------------c

      do 70 ip = 1, 4
         iq = ip

c--------------------------------------------------------------------------c
c     might look at non-diagonal correlator (not implemented)              c
c--------------------------------------------------------------------------c

         cprax(ip,iq) = 0.0
         cprvec(ip,iq) = 0.0
         do 80 i = 1, 4
         do 80 j = 1, 4
            do 90 k = 1, 3
            do 90 l = 1, 3
            cprax(ip,iq) = cprax(ip,iq)
     2                  + sgg5(k,i,l,j,ip)*(syxexp(l,j,k,i))
            cprvec(ip,iq) = cprvec(ip,iq)
     2                   + sgg(k,i,l,j,ip)*(syxexp(l,j,k,i))
   90       continue
   80    continue
   70 continue

c--------------------------------------------------------------------------c
c     contract scalars with second propagator                              c
c--------------------------------------------------------------------------c

      cprps = cmplx(0.0,0.0)
      cprsc = cmplx(0.0,0.0)
      do 85 i = 1, 4
      do 85 j = 1, 4
          do 95 k = 1, 3
          do 95 l = 1, 3
            cprps = cprps + s5(k,i,l,j)*(syxexp(l,j,k,i))
            cprsc = cprsc + sxy(k,i,l,j)*(syxexp(l,j,k,i))
   95     continue
   85 continue

      prps = real(cprps)
      prsc = real(cprsc)

c---------------------------------------------------------------------------c
c     trace vector correlators                                              c
c---------------------------------------------------------------------------c

      prav = 0.0
      prvv = 0.0
      do 100 ip = 1, 3
c        prax(ip,ip) = real(cprax(ip,ip))
c        prvec(ip,ip) = real(cprvec(ip,ip))
        prav = prav + real(cprax(ip,ip))
        prvv = prvv + real(cprvec(ip,ip))
  100 continue

      pra4 = real(cprax(4,4))
      prv4 = real(cprvec(4,4))

      pra1 = real(cprax(idir,idir))
      prv1 = real(cprvec(idir,idir))

      praxt = prav + pra4
      prvect= prvv + prv4

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine prmix(sxy,syx,pexp,pvs,ppa,pvt,pat)
c--------------------------------------------------------------------------c
c     calcualte meson bs amplitudes using non-diagonal correlators         c
c--------------------------------------------------------------------------c
c     input : sxy(a,i,b,j)  first propagator (a,b=1,..,3;i,j=1,..,4)       c
c             syx(a,i,b,j)  second propagator                              c
c             pexp(a,b)     pathexp (a,b=1,..,3)                           c
c     output: pvs           vector bs amplitude with scalar source         c
c             ppa           pseudoscalar with axial vector source          c
c             pvt           vector with tensor source                      c
c             pat           axial vector with tensor source                c
c--------------------------------------------------------------------------c
      parameter(nc=3)

      complex sgg(nc,4,nc,4)
      complex pexp(3,3)
      complex svt(3,4,3,4), svt5(3,4,3,4), s4(3,4,3,4)
      complex sxy(3,4,3,4), syx(3,4,3,4), s5(3,4,3,4)
      complex syxexp(3,4,3,4)
      complex gamma, gamu5
      complex cvs, cpa, cvt, cat
      complex g, gm5, guv

      common /gamf/  g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)
      common /gamf5/ gm5(4,4), nm5(4,4)
      common /guv/   guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)
      common /direc/ idir

c-------------------------------------------------------------------------c
c     multiply pexp by second propagator -->  syxexp(a,i;b,j)             c
c-------------------------------------------------------------------------c

      do 220 m1=1,4
      do 220 m2=1,4
      do 220 k1=1,3
      do 220 k2=1,3
         syxexp(k1,m1,k2,m2)=(0.,0.)
         do 221 kk=1,3
            syxexp(k1,m1,k2,m2)=syxexp(k1,m1,k2,m2)
     2                  + pexp(kk,k2)*syx(k1,m1,kk,m2)
221      continue
220   continue

c-------------------------------------------------------------------------c
c     sgg = \vec{gamma} S(x-y) \vec{gamma}                                c
c-------------------------------------------------------------------------c

      do 1 i = 1, 3
      do 1 j = 1, 3
        do 2 m = 1, 4
        do 2 n = 1, 4
          sgg(i,m,j,n) = cmplx(0.0,0.0)
  2     continue
  1   continue

      do 10 i = 1, 3
      do 10 j = 1, 3
        do 20 m = 1, 4
        do 20 n = 1, 4
          do 130 ip = 1, 3
          sgg(i,m,j,n) = sgg(i,m,j,n) +
     1                    g(ip,m)*sxy(i,nu(ip,m),j,ni(ip,n))
     1                     *g(ip,ni(ip,n))
 130      continue
  20    continue
  10  continue

c-------------------------------------------------------------------------c
c     contract \Gamma_1 S(x-y) \Gamma_2 (1=antenna, 2=source)             c
c-------------------------------------------------------------------------c

      do 510 i = 1, 3
      do 510 j = 1, 3
        do 520 m = 1, 4
        do 520 n = 1, 4
          s5(i,m,j,n)  = g(5,m)*sxy(i,nu(5,m),j,ni(5,ni(idir,n)))
     2                    *g(5,ni(5,ni(idir,n)))*g(idir,ni(idir,n))
          s4(i,m,j,n)  = sxy(i,m,j,ni(idir,n))
     2                    *g(idir,ni(idir,n))
          svt(i,m,j,n) = sgg(i,m,j,ni(idir,n))
     2                    *g(idir,ni(idir,n))
 520    continue
 510  continue

c-------------------------------------------------------------------------c
c     \gamma_5\gamma_0\vec\gamma = \eps_{ijk}\sigma_{jk} for av           c
c-------------------------------------------------------------------------c

      do 610 i = 1, 3
      do 610 j = 1, 3
        do 620 m = 1, 4
        do 620 n = 1, 4
          svt5(i,m,j,n) = g(5,m)*svt(i,nu(5,m),j,ni(5,n))
     1        * g(5,ni(5,n))
 620    continue
 610  continue

c-------------------------------------------------------------------------c
c     contract with second propagator                                     c
c-------------------------------------------------------------------------c

      cvs = cmplx(0.0,0.0)
      cpa = cmplx(0.0,0.0)
      cvt = cmplx(0.0,0.0)
      cat = cmplx(0.0,0.0)
      do 85 i = 1, 4
      do 85 j = 1, 4
         do 95 k = 1, 3
         do 95 l = 1, 3
            cvs = cvs + s4(k,i,l,j) * syxexp(l,j,k,i)
            cpa = cpa + s5(k,i,l,j) * syxexp(l,j,k,i)
            cvt = cvt + svt(k,i,l,j)* syxexp(l,j,k,i)
            cat = cat +svt5(k,i,l,j)* syxexp(l,j,k,i)
   95    continue
   85 continue

c-------------------------------------------------------------------------c
c     note that off-diagonal correlators might be imaginary               c
c-------------------------------------------------------------------------c

      pvs =aimag(cvs)
      ppa = real(cpa)
      pvt = real(cvt)
      pat =aimag(cat)

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine prhlm(sxy,pexp,xkp,xkm)
c--------------------------------------------------------------------------c
c     calcualate heavy-light meson correlators (incl. pexp)                c
c--------------------------------------------------------------------------c
c     input : sxy(a,i,b,j)  propagator (a,b=1,..,3;i,j=1,..,4)             c
c             pexp(a,b)     pathexp (a,b=1,..,3)                           c
c     output: xkp           1+\gamma_0 structure                           c
c             xkm           1-\gamma_0 structure                           c
c--------------------------------------------------------------------------c

      complex pexp(3,3)
      complex sxy(3,4,3,4), sxyexp(3,4,3,4), sxy0(3,4,3,4)
      complex ckp,ckm
      complex gamma, gamu5
      complex g, gm5, guv

      common /gamf/  g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)
      common /gamf5/ gm5(4,4), nm5(4,4)
      common /guv/   guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)
      common /direc/ idir

c-------------------------------------------------------------------------c
c     multiply propagator by pexp  -->  sxyexp(a,i;b,j)                   c
c-------------------------------------------------------------------------c

      do 10 m1=1,4
      do 10 m2=1,4
      do 10 k1=1,3
      do 10 k2=1,3
         sxyexp(k1,m1,k2,m2)=(0.,0.)
         do 11 kk=1,3
            sxyexp(k1,m1,k2,m2) = sxyexp(k1,m1,k2,m2)
     2                  + sxy(k1,m1,kk,m2)*pexp(kk,k2)
 11      continue
 10   continue

c--------------------------------------------------------------------------c
c     contract with \gamma_1                                               c
c--------------------------------------------------------------------------c

      do 20 i=1,3
      do 20 j=1,3
         do 20 m=1,4
         do 20 n=1,4
            sxy0(i,m,j,n) = g(idir,m)*sxyexp(i,nu(idir,m),j,n)
   20 continue

c---------------------------------------------------------------------------c
c     calculate traces                                                      c
c---------------------------------------------------------------------------c

      ckp = cmplx(0.0,0.0)
      ckm = cmplx(0.0,0.0)

      do 30 j = 1, 4
         do 40 k = 1, 3
            ckp = ckp + ( sxyexp(k,j,k,j) + sxy0(k,j,k,j) )/2.0
            ckm = ckm + ( sxyexp(k,j,k,j) - sxy0(k,j,k,j) )/2.0
  40     continue
  30  continue

      xkp = aimag(ckp)
      xkm = aimag(ckm)

      return
      end

c----------------------------------------------------------------------+-----
C----------------------------------------------------------------------+-----

      subroutine prdiq(qax,qay,pax,pexp,cdiq)
c----------------------------------------------------------------------------c
c     diquark correlators and wavefunctions                                  c
c----------------------------------------------------------------------------c
c     input :  qax(a,i;b,j)  quark propagator from origin to x               c
c              qay(a,i;b,j)  quark propagator from origin to y               c
c              pax(a,b)      pathexp connecting origin and x                 c
c              pexp(a,b)     pathexp connecting y and x                      c
c     output:  cdiq(m)       diquak correlators (m=1,..,5)                   c
c----------------------------------------------------------------------------c
c     cdiq(1)  scalar (Gamma) diquark                                        c
c     cdiq(2)  pseudoscalar diquark                                          c
c     cdiq(3)  vector diquark                                                c
c     cdiq(4)  axial vector diquark                                          c
c     cdiq(5)  tensor diquark                                                c
c----------------------------------------------------------------------------c

      complex cdiq(5)
      complex qax(3,4,3,4),qay(3,4,3,4)
      complex pax(3,3),pexp(3,3)
      complex qayexp(3,4,3,4)

      complex sc(3,4,3,4), s55(3,4,3,4)
      complex smumu(3,4,3,4), sm5m5(3,4,3,4)
      complex suvuv(3,4,3,4)

      complex guv,g

      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4),nic(4)
      common /guv/ guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)

c----------------------------------------------------------------------------c
c     \sum_{abc} \eps_{abc} F(a,b,c) = \sum_i F(eps(1,i),eps(2,i),eps(3,i))  c
c----------------------------------------------------------------------------c

      integer eps(3,6)
      data eps/1,2,3, 2,3,1, 3,1,2, 3,2,1, 2,1,3, 1,3,2/

c----------------------------------------------------------------------------c
c     multiply qay by pathexponent                                           c
c----------------------------------------------------------------------------c

      do 5 m1=1,4
      do 5 m2=1,4
      do 5 k1=1,3
      do 5 k2=1,3
         qayexp(k1,m1,k2,m2) = (0.,0.)
         do 7 kk=1,3
            qayexp(k1,m1,k2,m2) = qayexp(k1,m1,k2,m2)
     2                       + qay(k1,m1,kk,m2)*pexp(kk,k2)
  7      continue
  5   continue

c----------------------------------------------------------------------------c
c     calculate (\Gamma S(x-a) \Gamma)                                       c
c----------------------------------------------------------------------------c

      do 10 i = 1, 3
      do 10 j = 1, 3
      do 20 m = 1, 4
      do 20 n = 1, 4

      s55(i,m,j,n) = g(5,m)*qax(i,nu(5,m),j,ni(5,n))*g(5,ni(5,n))
      smumu(i,m,j,n) = g(1,m)*qax(i,nu(1,m),j,ni(1,n))*g(1,ni(1,n))
     2              + g(2,m)*qax(i,nu(2,m),j,ni(2,n))*g(2,ni(2,n))
     3              + g(3,m)*qax(i,nu(3,m),j,ni(3,n))*g(3,ni(3,n))
     4              + g(4,m)*qax(i,nu(4,m),j,ni(4,n))*g(4,ni(4,n))
      suvuv(i,m,j,n)=guv(1,m)*qax(i,nuv(1,m),j,nvu(1,n))*guv(1,nvu(1,n))
     1              +guv(2,m)*qax(i,nuv(2,m),j,nvu(2,n))*guv(2,nvu(2,n))
     2              +guv(3,m)*qax(i,nuv(3,m),j,nvu(3,n))*guv(3,nvu(3,n))
     3              +guv(4,m)*qax(i,nuv(4,m),j,nvu(4,n))*guv(4,nvu(4,n))
     4              +guv(5,m)*qax(i,nuv(5,m),j,nvu(5,n))*guv(5,nvu(5,n))
     5              +guv(6,m)*qax(i,nuv(6,m),j,nvu(6,n))*guv(6,nvu(6,n))

c----------------------------------------------------------------------------c
c     calculate (C S(y-a)^T C)                                               c
c----------------------------------------------------------------------------c

      sc(i,m,j,n) = c(m)*qayexp(i,nic(n),j,nuc(m))*c(nic(n))

  20  continue
  10  continue

c----------------------------------------------------------------------------c
c     (\gamma_5\gamma_mu S(y-a) \gamma_mu\gamma_5)                           c
c----------------------------------------------------------------------------c

      do 12 i = 1, 3
      do 12 j = 1, 3
        do 22 m = 1, 4
        do 22 n = 1, 4
          sm5m5(i,m,j,n) = g(5,m)*smumu(i,nu(5,m),j,ni(5,n))
     2                     *g(5,ni(5,n))*(-1.0)
  22    continue
  12  continue

c----------------------------------------------------------------------------c
c     clear arrays                                                           c
c----------------------------------------------------------------------------c

      do 13 i = 1, 5
        cdiq(i) = cmplx(0.0,0.0)
   13 continue

c----------------------------------------------------------------------------c
c     diquark correlators, sum over color (\sum_{ij})                        c
c----------------------------------------------------------------------------c

      do 15 i = 1, 6
      do 15 j = 1, 6

c----------------------------------------------------------------------------c
c     i,j>3 correspond to odd permutations of (abc)=(123)                    c
c----------------------------------------------------------------------------c

      sign = 1
      if (i .gt. 3) sign=-sign
      if (j .gt. 3) sign=-sign

c----------------------------------------------------------------------------c
c     contract with second propagator and pexp                               c
c----------------------------------------------------------------------------c

      do 14 m = 1, 4
      do 14 n = 1, 4
        cdiq(1) = cdiq(1) + sign * pax(eps(1,i),eps(1,j)) *
     2      qax(eps(2,i),m,eps(2,j),n) * sc(eps(3,i),n,eps(3,j),m)
        cdiq(2) = cdiq(2) + sign * pax(eps(1,i),eps(1,j)) *
     2      s55(eps(2,i),m,eps(2,j),n) * sc(eps(3,i),n,eps(3,j),m)
        cdiq(3) = cdiq(3) + sign * pax(eps(1,i),eps(1,j)) *
     2      smumu(eps(2,i),m,eps(2,j),n) * sc(eps(3,i),n,eps(3,j),m)
        cdiq(4) = cdiq(4) + sign * pax(eps(1,i),eps(1,j)) *
     2      sm5m5(eps(2,i),m,eps(2,j),n) * sc(eps(3,i),n,eps(3,j),m)
        cdiq(5) = cdiq(5) + sign * pax(eps(1,i),eps(1,j)) *
     2      suvuv(eps(2,i),m,eps(2,j),n) * sc(eps(3,i),n,eps(3,j),m)
   14 continue

   15 continue

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine nprop(qax,qay,qaz,pexp,pexp2,cl)
c----------------------------------------------------------------------------c
c     nucleon correlators and wavefunctions. all propagators are assumed to  c
c     be different. note that qay is the d-quark propagator.                 c
c----------------------------------------------------------------------------c
c     input :  qax(a,i;b,j)  quark propagator from origin to x               c
c              qay(a,i;b,j)  quark propagator from origin to y               c
c              qaz(a,i;b,j)  quark propagator from origin to z               c
c              pexp(a,b)     pathexp connecting y and x                      c
c              pexp2(a,b)    pathexp connecting z and x                      c
c     output:  cl(m)         nucleon correlators (m=1,..,6)                  c
c----------------------------------------------------------------------------c
c     cl(1)    <\eta_1\bar\eta_1>                                            c
c     cl(2)    <\eta_1\bar\eta_1\gamma_4>                                    c
c     cl(3)    <\eta_2\bar\eta_2>                                            c
c     cl(4)    <\eta_2\bar\eta_2\gamma_4>                                    c
c     cl(5)    <\eta_1\bar\eta_2>                                            c
c     cl(6)    <\eta_1\bar\eta_2\gamma_4>                                    c
c----------------------------------------------------------------------------c

      complex cl(6)
      complex qax(3,4,3,4),qay(3,4,3,4),qaz(3,4,3,4)
      complex pexp(3,3),pexp2(3,3)
      complex qayexp(3,4,3,4),qazexp(3,4,3,4)

      complex sxr5(3,4,3,4), sxl5(3,4,3,4), sx55(3,4,3,4)
      complex szr5(3,4,3,4), szl5(3,4,3,4), sz55(3,4,3,4)
      complex sc(3,4,3,4)

      complex sxr(3,4,3,4), sxl(3,4,3,4)
      complex sxrl(3,4,3,4), sxlr(3,4,3,4)
      complex szr(3,4,3,4), szl(3,4,3,4)
      complex szrl(3,4,3,4), szlr(3,4,3,4)

      complex prllr(4,4), plrrl(4,4), prrrr(4,4), pllll(4,4)
      complex prlll(4,4), plrrr(4,4)
      complex grrll(4,4), gllrr(4,4), grlrl(4,4), glrlr(4,4)
      complex grrlr(4,4), gllrl(4,4)

      complex srllr, slrrl, srrrr, sllll
      complex srlll, slrrr
      complex srrll, sllrr, srlrl, slrlr
      complex srrlr, sllrl

      complex ssrr, ssll, ssrl, sslr
      complex ssrl4, sslr4, sll1, srr1

      complex sl4(4,4), sr4(4,4), slr4(4,4), srl4(4,4)
      complex tlr4(4,4), trl4(4,4)

      complex g, guv

      common /gamf/  g(5,4), nu(5,4), ni(5,4), c(4), nuc(4),nic(4)
      common /guv/   guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)
      common /direc/ idir

c----------------------------------------------------------------------------c
c     \sum_{abc} \eps_{abc} F(a,b,c) = \sum_i F(eps(1,i),eps(2,i),eps(3,i))  c
c----------------------------------------------------------------------------c

      integer eps(3,6)
      data eps/1,2,3, 2,3,1, 3,1,2, 3,2,1, 2,1,3, 1,3,2/

c----------------------------------------------------------------------------c
c     multiply qay and qaz by pathexponents                                  c
c----------------------------------------------------------------------------c

      do 5 m1=1,4
      do 5 m2=1,4
      do 5 k1=1,3
      do 5 k2=1,3
         qayexp(k1,m1,k2,m2) = (0.,0.)
         qazexp(k1,m1,k2,m2) = (0.,0.)
         do 7 kk=1,3
            qayexp(k1,m1,k2,m2) = qayexp(k1,m1,k2,m2)
     2                       + qay(k1,m1,kk,m2)*pexp(kk,k2)
            qazexp(k1,m1,k2,m2) = qazexp(k1,m1,k2,m2)
     2                       + qaz(k1,m1,kk,m2)*pexp2(kk,k2)
  7      continue
  5   continue

c----------------------------------------------------------------------------c
c     multiply qax and qaz with \gamma_5 frm l,r                             c
c----------------------------------------------------------------------------c

      do 10 i = 1, 3
      do 10 j = 1, 3
        do 20 m = 1, 4
        do 20 n = 1, 4

          sxl5(i,m,j,n) = g(5,m)*qax(i,nu(5,m),j,n)
          sxr5(i,m,j,n) = qax(i,m,j,ni(5,n))*g(5,ni(5,n))
          sx55(i,m,j,n) = g(5,m)*qax(i,nu(5,m),j,ni(5,n))
     2                                            *g(5,ni(5,n))

          szl5(i,m,j,n) = g(5,m)*qazexp(i,nu(5,m),j,n)
          szr5(i,m,j,n) = qazexp(i,m,j,ni(5,n))*g(5,ni(5,n))
          sz55(i,m,j,n) = g(5,m)*qazexp(i,nu(5,m),j,ni(5,n))
     2                                            *g(5,ni(5,n))

c----------------------------------------------------------------------------c
c     sc = C S(a,y)^T C                                                      c
c----------------------------------------------------------------------------c

          sc(i,m,j,n) = c(m)*qayexp(i,nic(n),j,nuc(m))*c(nic(n))

  20    continue

c----------------------------------------------------------------------------c
c     calculate l/r projections of qax and qaz                               c
c----------------------------------------------------------------------------c

        do 30 m = 1, 4
        do 30 n = 1, 4

          sxr(i,m,j,n) = (qax(i,m,j,n)-sxl5(i,m,j,n)-sxr5(i,m,j,n)
     2                  +sx55(i,m,j,n))/4.
          sxl(i,m,j,n) = (qax(i,m,j,n)+sxl5(i,m,j,n)+sxr5(i,m,j,n)
     2                  +sx55(i,m,j,n))/4.
          sxrl(i,m,j,n)= (qax(i,m,j,n)-sxl5(i,m,j,n)+sxr5(i,m,j,n)
     2                  -sx55(i,m,j,n))/4.
          sxlr(i,m,j,n)= (qax(i,m,j,n)+sxl5(i,m,j,n)-sxr5(i,m,j,n)
     2                  -sx55(i,m,j,n))/4.

          szr(i,m,j,n) = (qazexp(i,m,j,n)-szl5(i,m,j,n)-szr5(i,m,j,n)
     2                  +sz55(i,m,j,n))/4.
          szl(i,m,j,n) = (qazexp(i,m,j,n)+szl5(i,m,j,n)+szr5(i,m,j,n)
     2                  +sz55(i,m,j,n))/4.
          szrl(i,m,j,n)= (qazexp(i,m,j,n)-szl5(i,m,j,n)+szr5(i,m,j,n)
     2                  -sz55(i,m,j,n))/4.
          szlr(i,m,j,n)= (qazexp(i,m,j,n)+szl5(i,m,j,n)-szr5(i,m,j,n)
     2                  -sz55(i,m,j,n))/4.

   30   continue

  10  continue

c----------------------------------------------------------------------------c
c     clear arrays for nucleon correlators                                   c
c----------------------------------------------------------------------------c

      do 40 i = 1, 6
        cl(i) = cmplx(0.0,0.0)
 40   continue

c----------------------------------------------------------------------------c
c     sum over color \sum_{abc}\sum_{a'b'c'}->\sum_{ij}                      c
c----------------------------------------------------------------------------c

      do 100 i=1,6
      do 100 j=1,6

c----------------------------------------------------------------------------c
c     i,j>3 correspond to odd permutations of (abc)=(123)                    c
c----------------------------------------------------------------------------c

      sign = 1
      if(i.gt.3)  sign=-sign
      if(j.gt.3)  sign=-sign

c----------------------------------------------------------------------------c
c     multilpy l/r projections with \gamma_1                                 c
c----------------------------------------------------------------------------c

      do 35 m = 1, 4
      do 35 n = 1, 4

        sl4(m,n)  = szl(eps(1,i),m,eps(3,j),ni(idir,n))
     2                *g(idir,ni(idir,n))
        sr4(m,n)  = szr(eps(1,i),m,eps(3,j),ni(idir,n))
     2                *g(idir,ni(idir,n))
        slr4(m,n) = szlr(eps(1,i),m,eps(3,j),ni(idir,n))
     2                *g(idir,ni(idir,n))
        srl4(m,n) = szrl(eps(1,i),m,eps(3,j),ni(idir,n))
     2                *g(idir,ni(idir,n))

c----------------------------------------------------------------------------c
c     tlr4 has different color structure                                     c
c----------------------------------------------------------------------------c

        tlr4(m,n) = szlr(eps(3,i),m,eps(3,j),ni(idir,n))
     2                *g(idir,ni(idir,n))
        trl4(m,n) = szrl(eps(3,i),m,eps(3,j),ni(idir,n))
     2                *g(idir,ni(idir,n))

 35   continue

c----------------------------------------------------------------------------c
c     multiply by second propagator                                          c
c----------------------------------------------------------------------------c

         do 60 m = 1, 4
         do 60 n = 1, 4

           prllr(m,n) = cmplx(0.0,0.0)
           plrrl(m,n) = cmplx(0.0,0.0)
           prrrr(m,n) = cmplx(0.0,0.0)
           pllll(m,n) = cmplx(0.0,0.0)
           prlll(m,n) = cmplx(0.0,0.0)
           plrrr(m,n) = cmplx(0.0,0.0)

           grrll(m,n) = cmplx(0.0,0.0)
           gllrr(m,n) = cmplx(0.0,0.0)
           grlrl(m,n) = cmplx(0.0,0.0)
           glrlr(m,n) = cmplx(0.0,0.0)
           grrlr(m,n) = cmplx(0.0,0.0)
           gllrl(m,n) = cmplx(0.0,0.0)

           do 50 k = 1, 4

c----------------------------------------------------------------------------c
c     prllr = r*S(z-a)*l l*S(x-a)*r   etc.                                   c
c----------------------------------------------------------------------------c

           prllr(m,n) = prllr(m,n) + szrl(eps(1,i),m,eps(3,j),k)
     2                 *sxlr(eps(3,i),k,eps(1,j),n)
           plrrl(m,n) = plrrl(m,n) + szlr(eps(1,i),m,eps(3,j),k)
     2                 *sxrl(eps(3,i),k,eps(1,j),n)
           prrrr(m,n) = prrrr(m,n) + szr(eps(1,i),m,eps(3,j),k)
     2                 *sxr(eps(3,i),k,eps(1,j),n)
           pllll(m,n) = pllll(m,n) + szl(eps(1,i),m,eps(3,j),k)
     2                 *sxl(eps(3,i),k,eps(1,j),n)
           prlll(m,n) = prlll(m,n) + szrl(eps(1,i),m,eps(3,j),k)
     2                 *sxl(eps(3,i),k,eps(1,j),n)
           plrrr(m,n) = plrrr(m,n) + szlr(eps(1,i),m,eps(3,j),k)
     2                 *sxr(eps(3,i),k,eps(1,j),n)

c----------------------------------------------------------------------------c
c     grrll = r*S(z-a)*r \gamma_0 l*S(x-a)*l                                 c
c----------------------------------------------------------------------------c

           grrll(m,n) = grrll(m,n) + sr4(m,k)
     2                 *sxl(eps(3,i),k,eps(1,j),n)
           gllrr(m,n) = gllrr(m,n) + sl4(m,k)
     2                 *sxr(eps(3,i),k,eps(1,j),n)
           grlrl(m,n) = grlrl(m,n) + srl4(m,k)
     2                 *sxrl(eps(3,i),k,eps(1,j),n)
           glrlr(m,n) = glrlr(m,n) + slr4(m,k)
     2                 *sxlr(eps(3,i),k,eps(1,j),n)
           grrlr(m,n) = grrlr(m,n) + sr4(m,k)
     2                 *sxlr(eps(3,i),k,eps(1,j),n)
           gllrl(m,n) = gllrl(m,n) + sl4(m,k)
     2                 *sxrl(eps(3,i),k,eps(1,j),n)

   50      continue
   60    continue

c----------------------------------------------------------------------------c
c    clear traces                                                            c
c----------------------------------------------------------------------------c

         sslr4 = cmplx(0.0, 0.0)
         ssrl4 = cmplx(0.0, 0.0)
         srr1  = cmplx(0.0, 0.0)
         sll1  = cmplx(0.0, 0.0)

         srllr = cmplx(0.0, 0.0)
         slrrl = cmplx(0.0, 0.0)
         srrrr = cmplx(0.0, 0.0)
         sllll = cmplx(0.0, 0.0)
         srlll = cmplx(0.0, 0.0)
         slrrr = cmplx(0.0, 0.0)
         srrll = cmplx(0.0, 0.0)
         sllrr = cmplx(0.0, 0.0)
         srlrl = cmplx(0.0, 0.0)
         slrlr = cmplx(0.0, 0.0)
         srrlr = cmplx(0.0, 0.0)
         sllrl = cmplx(0.0, 0.0)

         ssrr = cmplx(0.0, 0.0)
         ssll = cmplx(0.0, 0.0)
         ssrl = cmplx(0.0, 0.0)
         sslr = cmplx(0.0, 0.0)

         do 70 m = 1, 4

c----------------------------------------------------------------------------c
c     Tr(L*S(z-a)*R*gamma_4), .. Tr(R*S(z-a)*R), ..                          c
c----------------------------------------------------------------------------c

           sslr4 = sslr4 + tlr4(m,m)
           ssrl4 = ssrl4 + trl4(m,m)
           srr1 = srr1 + szr(eps(3,i),m,eps(3,j),m)
           sll1 = sll1 + szl(eps(3,i),m,eps(3,j),m)

           do 80 k = 1, 4

c----------------------------------------------------------------------------c
c     multply by third propagator and trace                                  c
c----------------------------------------------------------------------------c

             srllr = srllr + prllr(m,k)*sc(eps(2,i),k,eps(2,j),m)
             slrrl = slrrl + plrrl(m,k)*sc(eps(2,i),k,eps(2,j),m)
             srrrr = srrrr + prrrr(m,k)*sc(eps(2,i),k,eps(2,j),m)
             sllll = sllll + pllll(m,k)*sc(eps(2,i),k,eps(2,j),m)
             srlll = srlll + prlll(m,k)*sc(eps(2,i),k,eps(2,j),m)
             slrrr = slrrr + plrrr(m,k)*sc(eps(2,i),k,eps(2,j),m)

             srrll = srrll + grrll(m,k)*sc(eps(2,i),k,eps(2,j),m)
             sllrr = sllrr + gllrr(m,k)*sc(eps(2,i),k,eps(2,j),m)
             srlrl = srlrl + grlrl(m,k)*sc(eps(2,i),k,eps(2,j),m)
             slrlr = slrlr + glrlr(m,k)*sc(eps(2,i),k,eps(2,j),m)
             srrlr = srrlr + grrlr(m,k)*sc(eps(2,i),k,eps(2,j),m)
             sllrl = sllrl + gllrl(m,k)*sc(eps(2,i),k,eps(2,j),m)

c----------------------------------------------------------------------------c
c     Tr(R*S(x-a)*R*C*S(y-a)^T*C) ...                                        c
c----------------------------------------------------------------------------c

             ssrr = ssrr +
     2        sxr(eps(1,i),m,eps(1,j),k)*sc(eps(2,i),k,eps(2,j),m)
             ssll = ssll +
     2        sxl(eps(1,i),m,eps(1,j),k)*sc(eps(2,i),k,eps(2,j),m)
             ssrl = ssrl+
     2        sxrl(eps(1,i),m,eps(1,j),k)*sc(eps(2,i),k,eps(2,j),m)
             sslr = sslr +
     2        sxlr(eps(1,i),m,eps(1,j),k)*sc(eps(2,i),k,eps(2,j),m)


   80      continue
   70    continue

c----------------------------------------------------------------------------c
c     nucleon correlators                                                    c
c----------------------------------------------------------------------------c

         cl(1) = cl(1) + (srllr - ssrr*sll1 + slrrl
     2               - ssll*srr1)*sign
         cl(2) = cl(2) + (srrll - ssrl*sslr4
     2               + sllrr - sslr*ssrl4)*sign
         cl(3) = cl(3) + (srrrr - ssrr*srr1
     2               + sllll - ssll*sll1)*sign
         cl(4) = cl(4) + (srlrl - ssrl*ssrl4
     2               + slrlr - sslr*sslr4)*sign
         cl(5) = cl(5) + (srlll - ssrl*sll1
     2               + slrrr - sslr*srr1)*sign
         cl(6) = cl(6) + (srrlr - ssrr*sslr4
     2               + sllrl - ssll*ssrl4)*sign

c----------------------------------------------------------------------------c
c     end of sum over color indeces                                          c
c----------------------------------------------------------------------------c

  100 continue

      cl(1) = 16.0*cl(1)
      cl(2) = 16.0*cl(2)
      cl(3) = 64.0*cl(3)
      cl(4) = 64.0*cl(4)
      cl(5) = 32.0*cl(5)
      cl(6) = 32.0*cl(6)

      return
      end

c----------------------------------------------------------------------+-----
c----------------------------------------------------------------------+-----

      subroutine dprop(qax,qay,qaz,pexp,pexp2,cdl)
c----------------------------------------------------------------------------c
c     delta correlators and wavefunctions                                    c
c----------------------------------------------------------------------------c
c     input :  qax(a,i;b,j)  quark propagator from origin to x               c
c              qay(a,i;b,j)  quark propagator from origin to y               c
c              qaz(a,i;b,j)  quark propagator from origin to z               c
c              pexp(a,b)     pathexp connecting x and y                      c
c              pexp2(a,b)    pathexp connecting z and y                      c
c     output:  cdl(m)        delta correlators (m=1,..,4)                    c
c----------------------------------------------------------------------------c
c     cl(1)    Tr(\Pi_{mu mu})                                               c
c     cl(2)    Tr(\Pi_{mu mu}\gamma_4)                                       c
c     cl(3)    Tr(\Pi_{44})                                                  c
c     cl(4)    Tr(\Pi_{44}\gamma_4)                                          c
c----------------------------------------------------------------------------c

      complex cdl(4)
      complex qax(3,4,3,4),qay(3,4,3,4),qaz(3,4,3,4)
      complex pexp(3,3),pexp2(3,3)
      complex qayexp(3,4,3,4),qazexp(3,4,3,4)

      complex sr0, sr1
      complex sxmnc(4,4,4,4),symnc(4,4,4,4)
      complex txmnc,tymnc

      complex sx0(4,4), sy0(4,4), t0(4,4)
      complex s001(4,4), s000(4,4), sdl1(4,4), sdl0(4,4)
      complex t001, t000, tdl1, tdl0
      complex t200, t2sc

      complex g
      complex gamma, guv

      common /gamf/  g(5,4), nu(5,4), ni(5,4), c(4), nuc(4),nic(4)
      common /guv/   guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)
      common /direc/ idir


c----------------------------------------------------------------------------c
c     \sum_{abc} \eps_{abc} F(a,b,c) = \sum_i F(eps(1,i),eps(2,i),eps(3,i))  c
c----------------------------------------------------------------------------c

      integer eps(3,6)
      data eps/1,2,3, 2,3,1, 3,1,2, 3,2,1, 2,1,3, 1,3,2/

c----------------------------------------------------------------------------c
c     multiply qay and qaz by pathexponents                                  c
c----------------------------------------------------------------------------c

      do 5 m1=1,4
      do 5 m2=1,4
      do 5 k1=1,3
      do 5 k2=1,3
         qayexp(k1,m1,k2,m2) = (0.,0.)
         qazexp(k1,m1,k2,m2) = (0.,0.)
         do 7 kk=1,3
            qayexp(k1,m1,k2,m2) = qayexp(k1,m1,k2,m2)
     2                       + qay(k1,m1,kk,m2)*pexp(kk,k2)
            qazexp(k1,m1,k2,m2) = qazexp(k1,m1,k2,m2)
     2                       + qaz(k1,m1,kk,m2)*pexp2(kk,k2)
  7      continue
  5   continue

c----------------------------------------------------------------------------c
c     clear arrays for delta correlators                                     c
c----------------------------------------------------------------------------c

      do 6 i = 1, 4
        cdl(i) = cmplx(0.0,0.0)
    6 continue

c----------------------------------------------------------------------------c
c     sum over color \sum_{abc}\sum_{a'b'c'}->\sum_{ij}                      c
c----------------------------------------------------------------------------c

      do 100 i=1,6
      do 100 j=1,6

c----------------------------------------------------------------------------c
c     i,j>3 correspond to odd permutations of (abc)=(123)                    c
c----------------------------------------------------------------------------c

      sign = 1
      if(i.gt.3)  sign=-sign
      if(j.gt.3)  sign=-sign

      do 35 m = 1, 4
      do 35 n = 1, 4

        sx0(m,n) = g(idir,m)*qax(eps(3,i),nu(idir,m),eps(1,j),n)
        sy0(m,n) = g(idir,m)*qayexp(eps(3,i),nu(idir,m),eps(1,j),n)
        t0(m,n)  = g(idir,m)*qazexp(eps(3,i),nu(idir,m),eps(3,j),n)

   35 continue

c----------------------------------------------------------------------------c
c     sxmnc(n,m,k,l) = (\gamma_{m} C S(x-a)^T C \gamma_{n})_{k,l}            c
c----------------------------------------------------------------------------c

      do 213 k = 1, 4
      do 213 l = 1, 4
      do 214 m = 1, 4
      do 214 n = 1, 4
         sxmnc(n,m,k,l) = g(m,k)*c(nu(m,k))
     2     * qax(eps(2,i),nic(ni(n,l)),eps(2,j),nuc(nu(m,k)))
     3     * c(nic(ni(n,l)))*g(n,ni(n,l))
         symnc(n,m,k,l) = g(m,k)*c(nu(m,k))
     2     * qayexp(eps(2,i),nic(ni(n,l)),eps(2,j),nuc(nu(m,k)))
     3     * c(nic(ni(n,l)))*g(n,ni(n,l))
 214  continue
 213  continue

c----------------------------------------------------------------------------c
c     multiply smnc by second propagator                                     c
c----------------------------------------------------------------------------c

      do 65 m = 1, 4
      do 65 n = 1, 4
         s000(m,n) = cmplx(0.0,0.0)
         s001(m,n) = cmplx(0.0,0.0)
         sdl0(m,n) = cmplx(0.0,0.0)
         sdl1(m,n) = cmplx(0.0,0.0)


         do 55 k = 1, 4

c----------------------------------------------------------------------------c
c     s000(m,n) = (\gamma_4 S(x-a)\gamma_4 C S(y-a)^T C \gamma_4)_{mn}       c
c       +  ( x <-> y )                                                       c
c----------------------------------------------------------------------------c

            s001(m,n) = s001(m,n) +
     1       (  qax(eps(3,i),m,eps(1,j),k) * symnc(idir,idir,k,n)
     1       +qayexp(eps(3,i),m,eps(1,j),k)* sxmnc(idir,idir,k,n) )/2.0
            s000(m,n) = s000(m,n) +
     1           ( sx0(m,k)*symnc(idir,idir,k,n)
     1           + sy0(m,k)*sxmnc(idir,idir,k,n) )/2.0

c----------------------------------------------------------------------------c
c     sdl0(m,n) = (\gamma_4 S(x-a)\gamma_mu C S(y-a)^T gamma_mu)_{mn}        c
c       +  ( x <-> y )                                                       c
c----------------------------------------------------------------------------c

            txmnc =  sxmnc(1,1,k,n)+
     1           sxmnc(2,2,k,n)+sxmnc(3,3,k,n)+sxmnc(4,4,k,n)
            tymnc =  symnc(1,1,k,n)+
     1           symnc(2,2,k,n)+symnc(3,3,k,n)+symnc(4,4,k,n)

            sdl1(m,n) = sdl1(m,n) +
     1           (   qax(eps(3,i),m,eps(1,j),k) * tymnc
     2           + qayexp(eps(3,i),m,eps(1,j),k)* txmnc )/2.0
            sdl0(m,n) = sdl0(m,n) +
     1           (   sx0(m,k) * tymnc
     2           +   sy0(m,k) * txmnc )/2.0

  55     continue
  65  continue

c----------------------------------------------------------------------------c
c     clear traces                                                           c
c----------------------------------------------------------------------------c

      sr1 = cmplx(0.0,0.0)
      sr0 = cmplx(0.0,0.0)

      t000 = cmplx(0.0,0.0)
      t001 = cmplx(0.0,0.0)
      tdl0 = cmplx(0.0,0.0)
      tdl1 = cmplx(0.0,0.0)
      t200 = cmplx(0.0,0.0)
      t2sc = cmplx(0.0,0.0)

      do 70 m = 1, 4

c----------------------------------------------------------------------------c
c     sr1 = Tr(S(z-a)\gamma_4)                                               c
c----------------------------------------------------------------------------c

         sr1 = sr1 + qazexp(eps(3,i),m,eps(3,j),m)
         sr0 = sr0 + t0(m,m)

         do 80 k = 1, 4

c----------------------------------------------------------------------------c
c      txxx = Tr(S(z-a) sxxx)                                                c
c----------------------------------------------------------------------------c

            tdl1 = tdl1 + qazexp(eps(1,i),m,eps(3,j),k)*sdl1(k,m)
            tdl0 = tdl0 + qazexp(eps(1,i),m,eps(3,j),k)*sdl0(k,m)

            t000 = t000 + qazexp(eps(1,i),m,eps(3,j),k)*s000(k,m)
            t001 = t001 + qazexp(eps(1,i),m,eps(3,j),k)*s001(k,m)

c-    ---------------------------------------------------------------------------c
c     t200 = Tr(S(x-a)\gamma_4 C S(y-a)^T C \gamma_4)                        c
c----------------------------------------------------------------------------c

            t200 = t200 + qax(eps(1,i),m,eps(1,j),k)*
     2           symnc(idir,idir,k,m)
            t2sc = t2sc + qax(eps(1,i),m,eps(1,j),k)*
     2          (symnc(1,1,k,m)+symnc(2,2,k,m)+symnc(3,3,k,m)
     2              +symnc(4,4,k,m))


  80     continue
  70  continue

c----------------------------------------------------------------------------c
c     delta correlators                                                      c
c----------------------------------------------------------------------------c

         cdl(1) = cdl(1) + (4.0*tdl1 - 2.0*sr1*t2sc)*sign
         cdl(2) = cdl(2) + (4.0*tdl0 - 2.0*sr0*t2sc)*sign
         cdl(3) = cdl(3) + (4.0*t001 - 2.0*sr1*t200)*sign
         cdl(4) = cdl(4) + (4.0*t000 - 2.0*sr0*t200)*sign

c----------------------------------------------------------------------------c
c     end of sum over color indices                                          c
c----------------------------------------------------------------------------c

  100 continue

      return
      end

c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      subroutine delbs(sax,say,pexp,cdl)
c----------------------------------------------------------------------------c
c     symmetrized delta wave function. say is the d quark propagator in del+ c
c----------------------------------------------------------------------------c
c     input :  sax(a,i;b,j)  u-quark propagator from origin to x             c
c              say(a,i;b,j)  d-quark propagator from origin to y             c
c              pexp(a,b)     pathexp connecting y and x                      c
c     output:  cdl(m)        delta correlators (see dprop)                   c
c----------------------------------------------------------------------------c
      complex  sax(3,4,3,4),say(3,4,3,4)
      complex  pexp(3,3),pexp2(3,3)
      complex  cdl(4),cdl1(4),cdl2(4)

c----------------------------------------------------------------------------c
c     second pathexp is unit matrix                                          c
c----------------------------------------------------------------------------c

      do 10 i=1,3
      do 10 j=1,3
         pexp2(i,j) = (0.0,0.0)
  10  continue
      do 20 i=1,3
         pexp2(i,i) = (1.0,0.0)
  20  continue

c----------------------------------------------------------------------------c
c     calculate possible permutations of correlator                          c
c----------------------------------------------------------------------------c

      call dprop(sax,say,sax,pexp,pexp2,cdl1)
      call dprop(sax,sax,say,pexp2,pexp,cdl2)

      do 30 i=1,4
         cdl(i) = (2.0*cdl1(i)+cdl2(i))/3.0
  30  continue

      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine lprop(su,sd0,ss0,pexp,pexp2,cl)
c----------------------------------------------------------------------------c
c     lambda correlators and wavefunctions. this subroutine considers        c
c     scalar and pseudoscalar diquark lamda currents.                        c
c----------------------------------------------------------------------------c
c     input :  su(a,i;b,j)   u quark propagator from origin to x             c
c              sd0(a,i;b,j)  d quark propagator from origin to y             c
c              ss0(a,i;b,j)  s quark propagator from origin to z             c
c              pexp(a,b)     pathexp connecting y and x                      c
c              pexp2(a,b)    pathexp connecting z and x                      c
c     output:  cl(m)         lambda correlators (m=1,..,6)                   c
c----------------------------------------------------------------------------c
c     cl(1)    <\eta_1\bar\eta_1>                                            c
c     cl(2)    <\eta_1\bar\eta_1\gamma_4>                                    c
c     cl(3)    <\eta_2\bar\eta_2>                                            c
c     cl(4)    <\eta_2\bar\eta_2\gamma_4>                                    c
c     cl(5)    <\eta_1\bar\eta_2>                                            c
c     cl(6)    <\eta_1\bar\eta_2\gamma_4>                                    c
c----------------------------------------------------------------------------c

      complex cl(6)
      complex su(3,4,3,4), sd0(3,4,3,4), ss0(3,4,3,4)
      complex pexp(3,3), pexp2(3,3)
      complex sd(3,4,3,4), ss(3,4,3,4)

      complex sc(3,4,3,4)
      complex psu5(3,4,3,4), p5su5(3,4,3,4)
      complex p5ss(3,4,3,4), p5ss5(3,4,3,4)
      complex ss4(3,4,3,4), p5ss4(3,4,3,4), p5ss54(3,4,3,4)

      complex t, t4, t5, t54, t55, t554
      complex tss, tpp, tsp

      complex g, guv

      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4),nic(4)
      common /guv/  guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)
      common /direc/ idir

c----------------------------------------------------------------------------c
c     \sum_{abc} \eps_{abc} F(a,b,c) = \sum_i F(eps(1,i),eps(2,i),eps(3,i))  c
c----------------------------------------------------------------------------c

      integer eps(3,6)
      data eps/1,2,3, 2,3,1, 3,1,2, 3,2,1, 2,1,3, 1,3,2/

c----------------------------------------------------------------------------c
c     multiply sd0 and ss0 by pathexponents                                  c
c----------------------------------------------------------------------------c

      do 5 m1=1,4
      do 5 m2=1,4
      do 5 k1=1,3
      do 5 k2=1,3
         sd(k1,m1,k2,m2) = (0.,0.)
         ss(k1,m1,k2,m2) = (0.,0.)
         do 7 kk=1,3
            sd(k1,m1,k2,m2) = sd(k1,m1,k2,m2)
     2                       + sd0(k1,m1,kk,m2)*pexp(kk,k2)
            ss(k1,m1,k2,m2) = ss(k1,m1,k2,m2)
     2                       + ss0(k1,m1,kk,m2)*pexp2(kk,k2)
  7      continue
  5   continue

c----------------------------------------------------------------------------c
c     multiply ss and su with \gamma_5 frm l,r                               c
c----------------------------------------------------------------------------c

      do 10 i = 1, 3
      do 10 j = 1, 3
        do 20 m = 1, 4
        do 20 n = 1, 4

          p5ss(i,m,j,n) = g(5,m)*ss(i,nu(5,m),j,n)
          psu5(i,m,j,n) = su(i,m,j,ni(5,n))*g(5,ni(5,n))
          p5ss5(i,m,j,n)= g(5,m)*ss(i,nu(5,m),j,ni(5,n))
     2                                            *g(5,ni(5,n))
          p5su5(i,m,j,n)= g(5,m)*su(i,nu(5,m),j,ni(5,n))
     2                                            *g(5,ni(5,n))

c----------------------------------------------------------------------------c
c     sc = C Sd(a,y)^T C                                                     c
c----------------------------------------------------------------------------c

          sc(i,m,j,n) = c(m)*sd(i,nic(n),j,nuc(m))*c(nic(n))

  20    continue

  10  continue

c----------------------------------------------------------------------------c
c     multiply by gamma_4                                                    c
c----------------------------------------------------------------------------c

      do 30 i = 1, 3
      do 30 j = 1, 3
        do 30 m = 1, 4
        do 30 n = 1, 4

          ss4(i,m,j,n)   =    ss(i,m,j,ni(idir,n))*g(idir,ni(idir,n))
          p5ss4(i,m,j,n) =  p5ss(i,m,j,ni(idir,n))*g(idir,ni(idir,n))
          p5ss54(i,m,j,n)= p5ss5(i,m,j,ni(idir,n))*g(idir,ni(idir,n))

30    continue

c----------------------------------------------------------------------------c
c     clear arrays for lambda correlators                                    c
c----------------------------------------------------------------------------c

      do 40 i = 1, 6
        cl(i) = cmplx(0.0,0.0)
 40   continue

c----------------------------------------------------------------------------c
c     sum over color \sum_{abc}\sum_{a'b'c'}->\sum_{ij}                      c
c----------------------------------------------------------------------------c

      do 100 i=1,6
      do 100 j=1,6

c----------------------------------------------------------------------------c
c     i,j>3 correspond to odd permutations of (abc)=(123)                    c
c----------------------------------------------------------------------------c

      sign = 1
      if(i.gt.3)  sign=-sign
      if(j.gt.3)  sign=-sign

c----------------------------------------------------------------------------c
c     clear traces                                                           c
c----------------------------------------------------------------------------c

      t   = (0.0,0.0)
      t4  = (0.0,0.0)
      t5  = (0.0,0.0)
      t54 = (0.0,0.0)
      t55 = (0.0,0.0)
      t554= (0.0,0.0)
      tss = (0.0,0.0)
      tpp = (0.0,0.0)
      tsp = (0.0,0.0)

c----------------------------------------------------------------------------c
c     take traces                                                            c
c----------------------------------------------------------------------------c

      do 50 m=1,4

         t   = t   +    ss(eps(1,i),m,eps(1,j),m)
         t4  = t4  +   ss4(eps(1,i),m,eps(1,j),m)
         t5  = t5  +  p5ss(eps(1,i),m,eps(1,j),m)
         t54 = t54 + p5ss4(eps(1,i),m,eps(1,j),m)
         t55 = t55 + p5ss5(eps(1,i),m,eps(1,j),m)
         t554= t554+p5ss54(eps(1,i),m,eps(1,j),m)

c----------------------------------------------------------------------------c
c     multiply by second propagator and trace                                c
c----------------------------------------------------------------------------c

         do 60 n=1,4

         tss = tss + su(eps(2,i),m,eps(2,j),n)
     1               * sc(eps(3,i),n,eps(3,j),m)
         tpp = tpp + p5su5(eps(2,i),m,eps(2,j),n)
     1               * sc(eps(3,i),n,eps(3,j),m)
         tsp = tsp + psu5(eps(2,i),m,eps(2,j),n)
     1               * sc(eps(3,i),n,eps(3,j),m)

 60      continue
 50      continue

c----------------------------------------------------------------------------c
c     lambda correlators                                                     c
c----------------------------------------------------------------------------c

         cl(1) = cl(1) + t   * tpp * sign
         cl(2) = cl(2) + t4  * tpp * sign
         cl(3) = cl(3) + t55 * tss * sign
         cl(4) = cl(4) + t554* tss * sign
         cl(5) = cl(5) + t5  * tsp * sign
         cl(6) = cl(6) + t54 * tsp * sign

c----------------------------------------------------------------------------c
c     end of sum over color indices                                          c
c----------------------------------------------------------------------------c

  100 continue

c----------------------------------------------------------------------------c
c     normalize like nucleon                                                 c
c----------------------------------------------------------------------------c

      cl(1) =- 8.0*cl(1)
      cl(2) =  8.0*cl(2)
      cl(3) =-48.0*cl(3)
      cl(4) = 48.0*cl(4)
      cl(5) = 8.0*sqrt(6.0)*cl(5)
      cl(6) = 8.0*sqrt(6.0)*cl(6)

      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine lprop2(su,sd0,ss0,pexp,pexp2,cl)
c----------------------------------------------------------------------------c
c     lambda correlators and wavefunctions. this subroutine uses Ioffe       c
c     currents generalized to the lambda case.                               c
c----------------------------------------------------------------------------c
c     input :  su(a,i;b,j)   u quark propagator from origin to x             c
c              sd0(a,i;b,j)  d quark propagator from origin to y             c
c              ss0(a,i;b,j)  s quark propagator from origin to z             c
c              pexp(a,b)     pathexp connecting y and x                      c
c              pexp2(a,b)    pathexp connecting z and x                      c
c     output:  cl(m)         lambda correlators (m=1,..,6)                   c
c----------------------------------------------------------------------------c
c     cl(1)    <\eta_1\bar\eta_1>                                            c
c     cl(2)    <\eta_1\bar\eta_1\gamma_4>                                    c
c     cl(3)    <\eta_2\bar\eta_2>                                            c
c     cl(4)    <\eta_2\bar\eta_2\gamma_4>                                    c
c     cl(5)    <\eta_1\bar\eta_2>                                            c
c     cl(6)    <\eta_1\bar\eta_2\gamma_4>                                    c
c----------------------------------------------------------------------------c

      complex cl(6)
      complex su(3,4,3,4), sd0(3,4,3,4), ss0(3,4,3,4)
      complex pexp(3,3), pexp2(3,3)
      complex sd(3,4,3,4), ss(3,4,3,4)

      complex sc(3,4,3,4)
      complex psu5(3,4,3,4), p5sd(3,4,3,4), p5sd5(3,4,3,4)

      complex pm5sn5(4,4,3,4,3,4), pmsn(4,4,3,4,3,4)
      complex pm5sn(4,4,3,4,3,4), pmsn5(4,4,3,4,3,4)
      complex pmsn54(4,4,3,4,3,4), pm5sn4(4,4,3,4,3,4)
      complex pmsns(4,4,4,4)

      complex tt1, tt2, tt3, tt4
      complex tmsn(4,4), tmsn4(4,4), tmsns(4,4)

      complex g, guv

      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4),nic(4)
      common /guv/  guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)
      common /direc/idir

c----------------------------------------------------------------------------c
c     \sum_{abc} \eps_{abc} F(a,b,c) = \sum_i F(eps(1,i),eps(2,i),eps(3,i))  c
c----------------------------------------------------------------------------c

      integer eps(3,6)
      data eps/1,2,3, 2,3,1, 3,1,2, 3,2,1, 2,1,3, 1,3,2/

c----------------------------------------------------------------------------c
c     multiply sd0 and ss0 by pathexponents                                  c
c----------------------------------------------------------------------------c

      do 5 m1=1,4
      do 5 m2=1,4
      do 5 k1=1,3
      do 5 k2=1,3
         sd(k1,m1,k2,m2) = (0.,0.)
         ss(k1,m1,k2,m2) = (0.,0.)
         do 7 kk=1,3
            sd(k1,m1,k2,m2) = sd(k1,m1,k2,m2)
     2                       + sd0(k1,m1,kk,m2)*pexp(kk,k2)
            ss(k1,m1,k2,m2) = ss(k1,m1,k2,m2)
     2                       + ss0(k1,m1,kk,m2)*pexp2(kk,k2)
  7      continue
  5   continue

c----------------------------------------------------------------------------c
c     multiply ss and su with \gamma_5 frm l,r                               c
c----------------------------------------------------------------------------c

      do 10 i = 1, 3
      do 10 j = 1, 3
        do 20 m = 1, 4
        do 20 n = 1, 4

          p5sd(i,m,j,n) = g(5,m)*sd(i,nu(5,m),j,n)
          psu5(i,m,j,n) = su(i,m,j,ni(5,n))*g(5,ni(5,n))
          p5sd5(i,m,j,n)= g(5,m)*sd(i,nu(5,m),j,ni(5,n))
     2                                            *g(5,ni(5,n))

c----------------------------------------------------------------------------c
c     sc = C Ss(a,y)^T C                                                     c
c----------------------------------------------------------------------------c

          sc(i,m,j,n) = c(m)*ss(i,nic(n),j,nuc(m))*c(nic(n))

  20    continue

  10  continue

c----------------------------------------------------------------------------c
c     construct (\gamma_\mu S(x) \gamma_\nu)                                 c
c----------------------------------------------------------------------------c

      do 30 na=1,4
      do 30 nb=1,4
         do 30 i=1,3
         do 30 j=1,3
         do 30 m=1,4
         do 30 n=1,4

         pmsn(na,nb,i,m,j,n)  = g(na,m)*su(i,nu(na,m),j,ni(nb,n))
     1                                    *g(nb,ni(nb,n))
         pmsn5(na,nb,i,m,j,n) =-g(na,m)*psu5(i,nu(na,m),j,ni(nb,n))
     1                                    *g(nb,ni(nb,n))
         pm5sn(na,nb,i,m,j,n) = g(na,m)*p5sd(i,nu(na,m),j,ni(nb,n))
     1                                    *g(nb,ni(nb,n))
         pm5sn5(na,nb,i,m,j,n)=-g(na,m)*p5sd5(i,nu(na,m),j,ni(nb,n))
     1                                    *g(nb,ni(nb,n))

 30   continue

c----------------------------------------------------------------------------c
c     multiply by gamma_4                                                    c
c----------------------------------------------------------------------------c

      do 35 na=1,4
      do 35 nb=1,4
         do 35 i = 1, 3
         do 35 j = 1, 3
         do 35 m = 1, 4
         do 35 n = 1, 4

         pm5sn4(na,nb,i,m,j,n)= pm5sn5(na,nb,i,m,j,ni(idir,n))
     1                                    *g(idir,ni(idir,n))
         pmsn54(na,nb,i,m,j,n)=  pmsn5(na,nb,i,m,j,ni(idir,n))
     1                                    *g(idir,ni(idir,n))

 35   continue

c----------------------------------------------------------------------------c
c     clear arrays for lambda correlators                                    c
c----------------------------------------------------------------------------c

      do 40 i = 1, 6
        cl(i) = cmplx(0.0,0.0)
 40   continue

c----------------------------------------------------------------------------c
c     sum over color \sum_{abc}\sum_{a'b'c'}->\sum_{ij}                      c
c----------------------------------------------------------------------------c

      do 100 i=1,6
      do 100 j=1,6

c----------------------------------------------------------------------------c
c     i,j>3 correspond to odd permutations of (abc)=(123)                    c
c----------------------------------------------------------------------------c

      sign = 1
      if(i.gt.3)  sign=-sign
      if(j.gt.3)  sign=-sign

c----------------------------------------------------------------------------c
c     clear traces                                                           c
c----------------------------------------------------------------------------c

      tt1 = (0.0,0.0)
      tt2 = (0.0,0.0)
      tt3 = (0.0,0.0)
      tt4 = (0.0,0.0)

      do 50 na=1,4
      do 50 nb=1,4
         tmsn(na,nb)  = (0.0,0.0)
         tmsn4(na,nb) = (0.0,0.0)
         tmsns(na,nb) = (0.0,0.0)
 50   continue

c----------------------------------------------------------------------------c
c     take traces                                                            c
c----------------------------------------------------------------------------c

      do 60 na=1,4
      do 60 nb=1,4

      do 70 k=1,4

         tmsn(na,nb)   = tmsn(na,nb)
     1               + pm5sn5(na,nb,eps(1,i),k,eps(1,j),k)
         tmsn4(na,nb)  = tmsn4(na,nb)
     1               + pm5sn4(na,nb,eps(1,i),k,eps(1,j),k)

c----------------------------------------------------------------------------c
c     multiply by C S^T(x) C and trace                                       c
c----------------------------------------------------------------------------c

         do 70 l=1,4

         tmsns(na,nb) = tmsns(na,nb) +
     1     pmsn(na,nb,eps(2,i),k,eps(2,j),l)*sc(eps(3,i),l,eps(3,j),k)

c----------------------------------------------------------------------------c
c     multiply by C S^T(x) C, no trace                                       c
c----------------------------------------------------------------------------c

         pmsns(na,nb,k,l) = (0.0,0.0)

         do 70 m=1,4

         pmsns(na,nb,k,l) = pmsns(na,nb,k,l) +
     1     pm5sn(na,nb,eps(1,i),k,eps(1,j),m)*sc(eps(2,i),m,eps(2,j),l)

 70      continue
 60      continue

c----------------------------------------------------------------------------c
c     contract vector indices                                                c
c----------------------------------------------------------------------------c

         do 110 na=1,4
         do 110 nb=1,4

         tt1 = tt1 + tmsn(na,nb) * tmsns(na,nb)
         tt2 = tt2 + tmsn4(na,nb)* tmsns(na,nb)

c----------------------------------------------------------------------------c
c     .. and contract with third propagator                                  c
c----------------------------------------------------------------------------c

         do 120 k=1,4
         do 120 l=1,4

         tt3 = tt3 + pmsns(na,nb,k,l)
     1                 * pmsn5(na,nb,eps(3,i),l,eps(3,j),k)
         tt4 = tt4 + pmsns(na,nb,k,l)
     1                 * pmsn54(na,nb,eps(3,i),l,eps(3,j),k)

 120     continue
 110     continue

c----------------------------------------------------------------------------c
c     lambda correlators, normalize to nucleon                               c
c----------------------------------------------------------------------------c

         cl(1) = cl(1) + 4.0/3.0*sign * ( tt1 - tt3 )
         cl(2) = cl(2) - 4.0/3.0*sign * ( tt2 - tt4 )

c----------------------------------------------------------------------------c
c     end of sum over color indices                                          c
c----------------------------------------------------------------------------c

  100 continue

      return
      end

c----------------------------------------------------------------------+-----
c----------------------------------------------------------------------+-----

      subroutine twoloop(sxx,syy,prt,prt5,prv,pra)
c-------------------------------------------------------------------------c
c     calculate twoloop part of scalar and pseudoscalar correlators       c
c-------------------------------------------------------------------------c
c     sxx(a,i,b,j)  propagator s(x,x)^{ab}_{ij}                           c
c     syy(a,i,b,j)  propagator s(y,y)^{ab}_{ij}                           c
c     prt           Tr(S(x))Tr(S(y))                                      c
c     prt5          Tr(S(x)\gam_5)Tr(S(y)\gam_5)                          c
c     prv           Tr(S(x)\gam_m)Tr(S(y)\gam_m)                          c
c     pra           Tr(S(x)\gam_5\gam_m)Tr(S(y)\gam_m\gam_5)              c
c-------------------------------------------------------------------------c

      complex sxx(3,4,3,4), syy(3,4,3,4)
      complex ctrx, ctry, ctrx5, ctry5, ctrv, ctra
      complex cvx(4), cvy(4), cax(4), cay(4)
      complex g

      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)

      ctrx = cmplx(0.0,0.0)
      ctry = cmplx(0.0,0.0)
      ctrx5= cmplx(0.0,0.0)
      ctry5= cmplx(0.0,0.0)
      ctrv = cmplx(0.0,0.0)
      ctra = cmplx(0.0,0.0)

      do 5 mu=1,4
         cvx(mu) = cmplx(0.0,0.0)
         cvy(mu) = cmplx(0.0,0.0)
         cax(mu) = cmplx(0.0,0.0)
         cay(mu) = cmplx(0.0,0.0)
   5  continue

c-------------------------------------------------------------------------c
c     calculate traces                                                    c
c-------------------------------------------------------------------------c

      do 15 k = 1, 3
         ctrx = ctrx
     2    + sxx(k,1,k,1)+sxx(k,2,k,2)+sxx(k,3,k,3)+sxx(k,4,k,4)
         ctry = ctry
     2    + syy(k,1,k,1)+syy(k,2,k,2)+syy(k,3,k,3)+syy(k,4,k,4)
         ctrx5= ctrx5
     2    - sxx(k,1,k,3)-sxx(k,2,k,4)-sxx(k,3,k,1)-sxx(k,4,k,2)
         ctry5= ctry5
     2    - syy(k,1,k,3)-syy(k,2,k,4)-syy(k,3,k,1)-syy(k,4,k,2)
         do 20 mu= 1, 4
         do 20 i = 1, 4
            cvx(mu) = cvx(mu) + g(mu,i)*sxx(k,nu(mu,i),k,i)
            cvy(mu) = cvy(mu) + g(mu,i)*syy(k,nu(mu,i),k,i)
            cax(mu) = cax(mu) + g(5,i)*g(mu,nu(5,i))
     1                  * sxx(k,nu(mu,nu(5,i)),k,i)
            cay(mu) = cay(mu) + g(5,i)*g(mu,nu(5,i))
     1                  * syy(k,nu(mu,nu(5,i)),k,i)
 20      continue
 15   continue

c-------------------------------------------------------------------------c
c     disconnected correlators                                            c
c-------------------------------------------------------------------------c

      prt  = real(ctrx*ctry)
      prt5 = real(ctrx5*ctry5)

      do 30 mu=1,4
         ctrv = ctrv + cvx(mu)*cvy(mu)
         ctra = ctra + cax(mu)*cay(mu)
 30   continue

      prv = real(ctrv)
      pra = real(ctra)

      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine cond(sxx,uu,uudd,uuuu,ua1)
c-------------------------------------------------------------------------c
c     calculate (scalar) four quark condensates from zero mode propagator c
c-------------------------------------------------------------------------c
c     sxx(a,i,b,j)  propagator s(x,x)^{ab}_{ij}                           c
c     uu,uudd       Tr(S(x,x)), Tr(S(x,x))^2                              c
c     uuuu          Tr(S(x,x))^2-Tr(S(x,x)S(x,x))                         c
c     ua1           det_f(\bar q_Lq_R+(L<->R))                            c
c-------------------------------------------------------------------------c

      complex sxx(3,4,3,4), mp(3,4,3,4)
      complex ctrs, ctrs2, cu5, cu5u5, cua1
      complex g

      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)

      ctrs = (0.0,0.0)
      ctrs2= (0.0,0.0)
      cu5  = (0.0,0.0)
      cu5u5= (0.0,0.0)

c-------------------------------------------------------------------------c
c     local value of quark condensate                                     c
c-------------------------------------------------------------------------c

      do 10 k=1,3
      do 10 i=1,4
         ctrs = ctrs + sxx(k,i,k,i)
 10   continue

c-------------------------------------------------------------------------c
c     scalar four quark condensates                                       c
c-------------------------------------------------------------------------c

      do 20 k1=1,3
      do 20 k2=1,3
      do 20 i1=1,4
      do 20 i2=1,4
         mp(k1,i1,k2,i2) =
     2          g(5,i1)*sxx(k1,nu(5,i1),k2,ni(5,i2))*g(5,ni(5,i2))
         ctrs2= ctrs2+ sxx(k1,i1,k2,i2)*sxx(k2,i2,k1,i1)
 20   continue

      do 30 k=1,3
      do 30 i=1,4
         cu5 = cu5 + mp(k,i,k,i)
 30   continue

      do 40 k1=1,3
      do 40 k2=1,3
      do 40 i1=1,4
      do 40 i2=1,4
         cu5u5 = cu5u5 + mp(k1,i1,k2,i2)*sxx(k2,i2,k1,i1)
 40   continue

c-------------------------------------------------------------------------c
c     normalize with respect to factorization                             c
c-------------------------------------------------------------------------c

      uu  = aimag(ctrs)
      uudd= real(ctrs**2)
      uuuu= real(ctrs**2-ctrs2)*12.0/11.0
      ua1 = real(ctrs**2 + ctrs2 + cu5**2 + cu5u5)*6.0/7.0
      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine mixed(sxx,gmunu,qsgq)
c-------------------------------------------------------------------------c
c     calculate mixed condensate <\bar q\sig_\mu\nu G_\mu\nu q>.          c
c-------------------------------------------------------------------------c
c     sxx(3,4,3,4)   local zero mode propagator                           c
c     gmunu(3,3,4,4) field strength tensor                                c
c     qsgq           mixed condensate                                     c
c-------------------------------------------------------------------------c

      complex sxx(3,4,3,4), gmunu(3,3,4,4), cqsgq
      complex guv

      common /guv/ guv(6,4), nuv(6,4), mui(6), mvi(6), nvu(6,4)

      cqsgq = (0.0,0.0)

      do 10 k=1,3
      do 10 l=1,3
      do 10 i=1,6
      do 10 j=1,4
         cqsgq = cqsgq + sxx(k,j,l,nvu(i,j))*guv(i,nvu(i,j))
     2                    *gmunu(l,k,mui(i),mvi(i))
 10   continue

      qsgq = 2.0*real(cqsgq)

      return
      end

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine moreops(sxx,o1,o2,o3,o4)
c-------------------------------------------------------------------------c
c     calculate chiral symmetry breaking four quark operators defined     c
c     kapusta and shuryak.                                                c
c-------------------------------------------------------------------------c
c     input :  sxx(a,i,b,j)  local zero mode propagator                   c
c     output:  o1,...,o4     four quark operators                         c
c-------------------------------------------------------------------------c

      complex sxx(3,4,3,4), saa(3,4,3,4)
      complex co1v, co2v, co3v, co4v
      complex co1a, co2a, co3a, co4a
      complex mp(3,4,3,4), mpc(3,4,3,4)
      complex mv(3,4,3,4), ma(3,4,3,4)
      complex mvc(3,4,3,4), mac(3,4,3,4)
      complex mv0(3,4,3,4), ma0(3,4,3,4)
      complex mv0c(3,4,3,4), ma0c(3,4,3,4)

      complex lambda, lf, g

      common /lam/  lambda(8,3,3)
      common /lamf/ lf(8,3), lu(8,3), li(8,3)
      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)

c-------------------------------------------------------------------------c
c     contract (\lambda^a S(x,x) \lambda^a)                               c
c-------------------------------------------------------------------------c

      do 10 k=1,3
      do 10 l=1,3
      do 10 i=1,4
      do 10 j=1,4
         saa(k,i,l,j) = (0.0,0.0)
         do 10 n=1,8
            saa(k,i,l,j) = saa(k,i,l,j) +
     1         lf(n,k)*sxx(lu(n,k),i,li(n,l),j)*lf(n,li(n,l))
 10   continue

c------------------------------------------------------------------------c
c     contract  (\Gamma Saa(x,x) \Gamma)                                 c
c------------------------------------------------------------------------c

      do 20 k=1,3
      do 20 l=1,3
      do 20 i=1,4
      do 20 j=1,4
         mp(k,i,l,j) = g(5,i)*sxx(k,nu(5,i),l,ni(5,j))*g(5,ni(5,j))
         mpc(k,i,l,j)= g(5,i)*saa(k,nu(5,i),l,ni(5,j))*g(5,ni(5,j))
 20   continue

      do 30 k=1,3
      do 30 l=1,3
      do 30 i=1,4
      do 30 j=1,4

         mu = 4

         mv0(k,i,l,j) =
     1        g(mu,i)*sxx(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
         ma0(k,i,l,j) =
     1       -g(mu,i)* mp(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))

         mv0c(k,i,l,j) =
     1        g(mu,i)*saa(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
         ma0c(k,i,l,j) =
     1       -g(mu,i)*mpc(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))

         mv(k,i,l,j) = (0.0,0.0)
         ma(k,i,l,j) = (0.0,0.0)
         mvc(k,i,l,j)= (0.0,0.0)
         mac(k,i,l,j)= (0.0,0.0)

         do 30 mu=1,3
         mv(k,i,l,j) = mv(k,i,l,j) +
     1        g(mu,i)*sxx(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
         ma(k,i,l,j) = ma(k,i,l,j) -
     1        g(mu,i)* mp(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
         mvc(k,i,l,j)= mvc(k,i,l,j) +
     1        g(mu,i)*saa(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
         mac(k,i,l,j)= mac(k,i,l,j) -
     1        g(mu,i)*mpc(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
 30   continue

c------------------------------------------------------------------------c
c     contract with second propagator                                    c
c------------------------------------------------------------------------c

      co1v = (0.0,0.0)
      co2v = (0.0,0.0)
      co3v = (0.0,0.0)
      co4v = (0.0,0.0)
      co1a = (0.0,0.0)
      co2a = (0.0,0.0)
      co3a = (0.0,0.0)
      co4a = (0.0,0.0)

      do 40 k=1,3
      do 40 l=1,3
      do 40 i=1,4
      do 40 j=1,4
         co1v = co1v + mv0(k,i,l,j)*sxx(l,j,k,i)
         co2v = co2v +  mv(k,i,l,j)*sxx(l,j,k,i)
         co3v = co3v +mv0c(k,i,l,j)*sxx(l,j,k,i)
         co4v = co4v + mvc(k,i,l,j)*sxx(l,j,k,i)
         co1a = co1a + ma0(k,i,l,j)*sxx(l,j,k,i)
         co2a = co2a +  ma(k,i,l,j)*sxx(l,j,k,i)
         co3a = co3a +ma0c(k,i,l,j)*sxx(l,j,k,i)
         co4a = co4a + mac(k,i,l,j)*sxx(l,j,k,i)
 40   continue

c------------------------------------------------------------------------c
c     normalize to naive factorization                                   c
c------------------------------------------------------------------------c

      o1 = real(co1v-co1a)/2.0*12.0
      o2 = real(co2v-co2a)/2.0*4.0
      o3 = real(co3v-co3a)/2.0*9.0/4.0
      o4 = real(co4v-co4a)/2.0*3.0/4.0

      return
      end

c---------------------------------------------------------------------+-----
c---------------------------------------------------------------------+-----

      subroutine fourq(sxx,ss,pp,vv,aa,tt)
c-------------------------------------------------------------------------c
c     calculate colored four quark condensates, normalize with respect to c
c     result from naive factorization ss=<s>^2, etc.                      c
c-------------------------------------------------------------------------c
c     input :  sxx(a,i,b,j)  local zero mode propagator                   c
c     output:  gg={ss,pp,..} Tr(S(x,x)\lam\gS(x,x)\lam\g) with            c
c                            g=1,\gamma_5,\gamma_\mu,\gamma_5\gamma_\mu   c
c-------------------------------------------------------------------------c

      complex sxx(3,4,3,4), saa(3,4,3,4)
      complex css, cpp, cvv, caa, ctt
      complex mp(3,4,3,4), mv(3,4,3,4), ma(3,4,3,4), mt(3,4,3,4)
      complex lambda, lf, g

      common /lam/  lambda(8,3,3)
      common /lamf/ lf(8,3), lu(8,3), li(8,3)
      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)

c-------------------------------------------------------------------------c
c     contract (\lambda^a S(x,x) \lambda^a)                               c
c-------------------------------------------------------------------------c

      do 10 k=1,3
      do 10 l=1,3
      do 10 i=1,4
      do 10 j=1,4
         saa(k,i,l,j) = (0.0,0.0)
         do 10 n=1,8
            saa(k,i,l,j) = saa(k,i,l,j) +
     1         lf(n,k)*sxx(lu(n,k),i,li(n,l),j)*lf(n,li(n,l))
 10   continue

c------------------------------------------------------------------------c
c     contract  (\Gamma Saa(x,x) \Gamma)                                 c
c------------------------------------------------------------------------c

      do 20 k=1,3
      do 20 l=1,3
      do 20 i=1,4
      do 20 j=1,4
         mp(k,i,l,j) = g(5,i)*saa(k,nu(5,i),l,ni(5,j))*g(5,ni(5,j))
 20   continue

      do 30 k=1,3
      do 30 l=1,3
      do 30 i=1,4
      do 30 j=1,4
         mv(k,i,l,j) = (0.0,0.0)
         ma(k,i,l,j) = (0.0,0.0)
         do 30 mu=1,4
         mv(k,i,l,j) = mv(k,i,l,j) +
     1        g(mu,i)*saa(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
         ma(k,i,l,j) = ma(k,i,l,j) +
     1        g(mu,i)* mp(k,nu(mu,i),l,ni(mu,j))*g(mu,ni(mu,j))
 30   continue

c------------------------------------------------------------------------c
c     contract with second propagator, sign in caa from permuting \g_5   c
c------------------------------------------------------------------------c

      css = (0.0,0.0)
      cpp = (0.0,0.0)
      cvv = (0.0,0.0)
      caa = (0.0,0.0)

      do 40 k=1,3
      do 40 l=1,3
      do 40 i=1,4
      do 40 j=1,4
         css = css + saa(k,i,l,j)*sxx(l,j,k,i)
         cpp = cpp +  mp(k,i,l,j)*sxx(l,j,k,i)
         cvv = cvv +  mv(k,i,l,j)*sxx(l,j,k,i)
         caa = caa -  ma(k,i,l,j)*sxx(l,j,k,i)
 40   continue

c------------------------------------------------------------------------c
c     normalize to naive factorization                                   c
c------------------------------------------------------------------------c

      ss = real(css)*9.0/4.0
      pp = real(cpp)*9.0/4.0
      vv = real(cvv)*9.0/16.0
      aa =-real(caa)*9.0/16.0
      tt = 0.0

      return
      end

c---------------------------------------------------------------------+-----
c---------------------------------------------------------------------+-----

      subroutine potsu3(ni,ndd,zr,e1,e2,xa,a)
c-------------------------------------------------------------------------c
c     evaluate vector potential a(i,j,mu) at point xa(k). collective      c
c     variables in zr,e1,e2. first ni/2 pseudoparticles are instantons    c
c     second half antiinstantons. vector potential is calculated in the   c
c     ratio ansatz.                                                       c
c-------------------------------------------------------------------------c
c     this version uses antihermitean generators.                         c
c-------------------------------------------------------------------------c
c     ni        number of instantons                                      c
c     ndd       maximum number of instantons (for array dimensions)       c
c     zr(i,5)   position, size of instanton i                             c
c     e1(i,6)   real, complex parts of first row of su(3) matrix          c
c     e2(i,6)   same for second row                                       c
c     xa(k)     coordinates of point                                      c
c     a(i,j,mu) vector potential at xa(k) (output)                        c
c-------------------------------------------------------------------------c

      dimension  zr(ndd,5),e1(ndd,6),e2(ndd,6),zz(256,4),zs(256)
      dimension y(4),z(4),a(3,3,4),xa(4),ai3(3,3,4),ai(3,4)
      complex a,ai3,uuu(256,3,3),u1(3,3),u2(3,3)

c-------------------------------------------------------------------------c
c     modify ratio ansatz by exponential tail                             c
c-------------------------------------------------------------------------c

        fdecay=0.2

c-------------------------------------------------------------------------c
c     initialize a(i,j,mu)=0                                              c
c-------------------------------------------------------------------------c

      do 1 k1=1,3
      do 1 k2=1,3
      do 1 m=1,4
         a(k1,k2,m)=0.
1     continue

c-------------------------------------------------------------------------c
c     shortest path to all instantons ( --> zz(i,4) )                     c
c-------------------------------------------------------------------------c

      call tordis(ni,xa,zr,zz,zs)

c-------------------------------------------------------------------------c
c     construct su(3) rotation matrix uuu(i,3,3) from e1,e2               c
c-------------------------------------------------------------------------c

       call su3(3,ni,ndd,e1,e2,uuu)

c-------------------------------------------------------------------------c
c      loop over all instantons                                           c
c-------------------------------------------------------------------------c

       f=1.
       do 2 in=1,ni

       rr = zs(in)
       do 3   m= 1,4
          z(m) = zz(in,m)
3      continue

       do 5 k1=1,3
          do 5 m=1,4
             ai(k1,m)=0.
             do 5 k2=1,3
                ai3(k1,k2,m)=0.
5      continue

c-------------------------------------------------------------------------c
c      factors for ratio ansatz                                           c
c-------------------------------------------------------------------------c

       ros = zr(in,5)**2
       ff  = exp(-fdecay*rr/ros)
       f   = f + ros*ff/rr
       z4  = 2.0*ff*ros/rr**2

c-------------------------------------------------------------------------c
c     signum=+- for instanton/antiinstanton                               c
c-------------------------------------------------------------------------c

      signum = 1 - 2*int( (2.0*in-1.0)/float(ni) )

c-------------------------------------------------------------------------c
c     ai(a,mu) = eta/etabar(a,mu,nu)*x(nu)                                c
c-------------------------------------------------------------------------c

      do 4 k=1,3
         ai(k,4) = signum*z(k)*z4
         ai(k,k) =-z(4)*z4*signum
4     continue

      ai(1,2) = z(3)*z4
      ai(2,3) = z(1)*z4
      ai(3,1) = z(2)*z4
      ai(1,3) =-z(2)*z4
      ai(3,2) =-z(1)*z4
      ai(2,1) =-z(3)*z4

c------------------------------------------------------------------------c
c     ai3(i,j,mu) = ai(a,mu) * 1/(2i) * tau(a)_(i,j)                     c
c------------------------------------------------------------------------c

      do 41 mm=1,4
         ai3(1,1,mm)= ai(3,mm)*0.5 *(0.,-1.)
         ai3(2,2,mm)=-ai(3,mm)*0.5 *(0.,-1.)
         ai3(1,2,mm)=(ai(1,mm) - ai(2,mm)*(0.,1.))*0.5*(0.,-1.)
         ai3(2,1,mm)=(ai(1,mm) + ai(2,mm)*(0.,1.))*0.5*(0.,-1.)
41    continue

c------------------------------------------------------------------------c
c     a(i,j,mu) -> a(i,j,mu) + u+(i,k) ai(k,l,mu) u(l,j)                 c
c------------------------------------------------------------------------c

      do 60  m= 1,4
         do 61 k1= 1,3
            do 61 k2= 1,3
               u1(k1,k2)=(0.,0.)
61    continue

      do 62 k1= 1,3
         do 62 k2= 1,3
            do 62 kk= 1,3
            u1(k1,k2)= u1(k1,k2)+ai3(k1,kk,m)*uuu(in,kk,k2)
62    continue

      do 63 k1= 1,3
         do 63 k2= 1,3
            do 63 kk= 1,3
            a(k1,k2,m)=a(k1,k2,m)+u1(kk,k2)*conjg(uuu(in,kk,k1))
63    continue

60    continue

c-------------------------------------------------------------------------c
c     end of loop over all instantons                                     c
c-------------------------------------------------------------------------c

2     continue

c-------------------------------------------------------------------------c
c     ratio ansatz                                                        c
c-------------------------------------------------------------------------c

      fact = 1/f
      do 7  m= 1,4
      do 7  k1= 1,3
      do 7  k2= 1,3
         a(k1,k2,m)= a(k1,k2,m)*fact
7     continue
      return
      end

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine sumsu3(ni,ndd,zr,e1,e2,xa,a)
c-------------------------------------------------------------------------c
c     evaluate vector potential a(i,j,mu) at point xa(k). collective      c
c     variables in zr,e1,e2. first ni/2 pseudoparticles are instantons    c
c     second half antiinstantons. vector potential is calculated in the   c
c     sum ansatz.                                                         c
c-------------------------------------------------------------------------c
c     this version assumes antihermitean generators.                      c
c-------------------------------------------------------------------------c
c     ni        number of instantons                                      c
c     ndd       maximum number of instantons (for array dimensions)       c
c     zr(i,5)   position, size of instanton i                             c
c     e1(i,6)   real, complex parts of first row of su(3) matrix          c
c     e2(i,6)   same for second row                                       c
c     xa(k)     coordinates of point                                      c
c     a(i,j,mu) vector potential at xa(k) (output)                        c
c-------------------------------------------------------------------------c

      dimension  zr(ndd,5),e1(ndd,6),e2(ndd,6),zz(256,4),zs(256)
      dimension y(4),z(4),a(3,3,4),xa(4),ai3(3,3,4),ai(3,4)
      complex a,ai3,uuu(256,3,3),u1(3,3),u2(3,3)

c-------------------------------------------------------------------------c
c     initialize a(i,j,mu)=0                                              c
c-------------------------------------------------------------------------c

      do 1 k1=1,3
      do 1 k2=1,3
      do 1 m=1,4
         a(k1,k2,m)=0.
1     continue

c-------------------------------------------------------------------------c
c     shortest path to all instantons ( --> zz(i,4) )                     c
c-------------------------------------------------------------------------c

      call tordis(ni,xa,zr,zz,zs)

c-------------------------------------------------------------------------c
c     construct su(3) rotation matrix uuu(i,3,3) from e1,e2               c
c-------------------------------------------------------------------------c

       call su3(3,ni,ndd,e1,e2,uuu)

c-------------------------------------------------------------------------c
c      loop over all instantons                                           c
c-------------------------------------------------------------------------c

       f=1.
       do 2 in=1,ni

       rr = zs(in)
       do 3   m= 1,4
          z(m) = zz(in,m)
3      continue

       do 5 k1=1,3
          do 5 m=1,4
             ai(k1,m)=0.
             do 5 k2=1,3
                ai3(k1,k2,m)=0.
5      continue

c-------------------------------------------------------------------------c
c      instanton profile                                                  c
c-------------------------------------------------------------------------c

       ros = zr(in,5)**2
       z4  = 2.0*ros/rr/(ros+rr)

c-------------------------------------------------------------------------c
c     signum=+- for instanton/antiinstanton                               c
c-------------------------------------------------------------------------c

      signum = 1 - 2*int( (2.0*in-1.0)/float(ni) )
      if (ni .eq. 1) signum=1

c-------------------------------------------------------------------------c
c     ai(a,mu) = etabar/eta(a,mu,nu)*x(nu) for I,A                        c
c-------------------------------------------------------------------------c

      do 4 k=1,3
         ai(k,4) = signum*z(k)*z4
         ai(k,k) =-z(4)*z4*signum
4     continue

      ai(1,2) = z(3)*z4
      ai(2,3) = z(1)*z4
      ai(3,1) = z(2)*z4
      ai(1,3) =-z(2)*z4
      ai(3,2) =-z(1)*z4
      ai(2,1) =-z(3)*z4

c------------------------------------------------------------------------c
c     ai3(i,j,mu) = ai(a,mu) * 1/(2i) * tau(a)_(i,j)                     c
c------------------------------------------------------------------------c

      do 41 mm=1,4
         ai3(1,1,mm)= ai(3,mm)*0.5*(0.,-1.)
         ai3(2,2,mm)=-ai(3,mm)*0.5*(0.,-1.)
         ai3(1,2,mm)=(ai(1,mm) - ai(2,mm)*(0.,1.))*0.5*(0.,-1.)
         ai3(2,1,mm)=(ai(1,mm) + ai(2,mm)*(0.,1.))*0.5*(0.,-1.)
41    continue

c------------------------------------------------------------------------c
c     a(i,j,mu) -> a(i,j,mu) + u+(i,k) ai(k,l,mu) u(l,j)                 c
c------------------------------------------------------------------------c

      do 60  m= 1,4
         do 61 k1= 1,3
            do 61 k2= 1,3
               u1(k1,k2)=(0.,0.)
61    continue

      do 62 k1= 1,3
         do 62 k2= 1,3
            do 62 kk= 1,3
            u1(k1,k2)= u1(k1,k2)+ai3(k1,kk,m)*uuu(in,kk,k2)
62    continue

      do 63 k1= 1,3
         do 63 k2= 1,3
            do 63 kk= 1,3
            a(k1,k2,m)=a(k1,k2,m)+u1(kk,k2)*conjg(uuu(in,kk,k1))
63    continue

60    continue

c-------------------------------------------------------------------------c
c     end of loop over all instantons                                     c
c-------------------------------------------------------------------------c

2     continue

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine pathexp(ni,ndd,zr,e1,e2,xxsm,xxsp,u)
c-------------------------------------------------------------------------c
c     calculate pathexp connecting endpoint of xxsm and starting point of c
c     xxsp. uses subroutine potsu3 to calculate vectorpotential.          c
c-------------------------------------------------------------------------c
c     this version assumes antihermitean generators.                      c
c-------------------------------------------------------------------------c
c     ni        number of instantons                                      c
c     ndd       maximum number of instantons (for array dimensions)       c
c     zr(i,5)   position, size of instanton i                             c
c     e1(i,6)   real, complex parts of first row of su(3) matrix          c
c     e2(i,6)   same for second row                                       c
c     xxsp(k,i) endpoints of first propagator (k=1,..,4 i=1,2)            c
c     xxsm(k,i) same for second propagator                                c
c     u(a,b)    pathexp (a,b=1,..,3)                                      c
c-------------------------------------------------------------------------c

      complex u(3,3),uu(3,3), a(3,3,4),
     1        adx(3,3),adx2(3,3),adx3,dua(3,3),dua1(3,3)
      dimension xa(4),xb(4),xx(4),dx(4),xi(4),xf(4),z(4)
      dimension zr(ndd,5),e1(ndd,6),e2(ndd,6),xxsm(4,2),xxsp(4,2)

c-------------------------------------------------------------------------c
c     initialize u=1                                                      c
c-------------------------------------------------------------------------c

      nc=3
      do 1 k1=1,nc
      do 1 k2=1,nc
         u(k2,k1)=(0.,0.)
         if(k1.eq.k2) u(k2,k1)=(1.,0.)
1     continue

c-------------------------------------------------------------------------c
c     determine vector z(k) connecting endpoints                          c
c-------------------------------------------------------------------------c

      rr = 0.0
      do 50 m=1,4
         xi(m)= xxsm(m,2)
         z(m) = xxsp(m,1) - xxsm(m,2)
         rr   = rr + z(m)**2
50    continue

c-------------------------------------------------------------------------c
c     determine stepsize dx(k)  (step used to be 0.2 fm)                  c
c-------------------------------------------------------------------------c

      step = 0.1
      n    = sqrt(rr)/step
      if (n.lt.2) n=2
      do 40 mu=1,4
         dx(mu)=z(mu)/float(n-1)
40    continue

c------------------------------------------------------------------------c
c     loop over intermediate steps                                       c
c------------------------------------------------------------------------c

      do 10 i=1,n

      do 2 mu=1,4
         xa(mu) = xi(mu) + (i-1)*dx(mu)
         xb(mu) = xa(mu) + dx(mu)
         xx(mu) = 0.5*( xb(mu) + xa(mu) )
2     continue

c------------------------------------------------------------------------c
c     calculate vector potential a(i,j,mu) at point xa(mu)               c
c------------------------------------------------------------------------c

      call potsu3(ni,ndd,zr,e1,e2,xa,a)

c------------------------------------------------------------------------c
c     calculate adx=-a(mu)*dx(mu)                                        c
c------------------------------------------------------------------------c
c     for hermitean generators change to adx=i*a(mu)*dx(mu)              c
c------------------------------------------------------------------------c

      znorm=0.
      do 4 k1=1,nc
      do 4 k2=1,nc
         adx3=(0.,0.)
         do 3 mu=1,4
            adx3=adx3+a(k2,k1,mu)*dx(mu)
3        continue
         adx(k2,k1)= (-1.0,0.0)*adx3
         znorm=znorm+adx3*conjg(adx3)
4     continue

c------------------------------------------------------------------------c
c     if norm(adx) is large, split dx in 2**knorm intervalls             c
c------------------------------------------------------------------------c

      knorm=0
      zmax=0.5
      if(znorm.lt.zmax)go to 42
      knorm=0.5*alog(znorm/zmax)/0.68+1
      fnorm=0.5**knorm

c------------------------------------------------------------------------c
c     splitting only used to evaluate exp, a(mu) is assumed to be const. c
c------------------------------------------------------------------------c

      do 41 k1=1,nc
         do 41 k2=1,nc
            adx(k2,k1)=adx(k2,k1)*fnorm
41       continue
42    continue

c------------------------------------------------------------------------c
c     calculation of  (a*dx)**2                                          c
c------------------------------------------------------------------------c

      do 5 k1=1,nc
         do 5 k2=1,nc
            adx2(k2,k1)=(0.,0.)
            do 5 k3=1,nc
              adx2(k2,k1)=adx2(k2,k1)+adx(k2,k3)*adx(k3,k1)
5     continue

c------------------------------------------------------------------------c
c    calculation of  du=exp(i*a*dx)   up to third order                  c
c------------------------------------------------------------------------c

      do 7 k1=1,nc
         do 7 k2=1,nc
            adx3=(0.,0.)
            do 6 k3=1,nc
               adx3=adx3+adx2(k2,k3)*adx(k3,k1)
6           continue
            dua(k2,k1)=adx(k2,k1)+0.5*adx2(k2,k1)+0.1666667*adx3
            if(k2.eq.k1) dua(k2,k1)=dua(k2,k1)+(1.,0.)
7     continue

c------------------------------------------------------------------------c
c     if interval was split, evaluate dua**(2**knorm)                    c
c------------------------------------------------------------------------c

      if(knorm.eq.0) go to 74
      do 73 inorm=1,knorm
         do 71 k1=1,nc
            do 71 k2=1,nc
               dua1(k2,k1)=(0.,0.)
               do 71 k3=1,nc
                  dua1(k2,k1)=dua1(k2,k1)+dua(k2,k3)*dua(k3,k1)
71       continue

         do 72 k1=1,nc
            do 72 k2=1,nc
               dua(k2,k1)=dua1(k2,k1)
72       continue
73    continue

74    continue

c------------------------------------------------------------------------c
c     collect contribution from dx, u -> u*dua                           c
c------------------------------------------------------------------------c

      do 91 k1=1,nc
         do 91 k2=1,nc
            uu(k2,k1)=(0.,0.)
            do 91 k3=1,nc
               uu(k2,k1)=uu(k2,k1)+dua(k2,k3)*u(k3,k1)
91    continue

      do 92 k1=1,nc
         do 92 k2=1,nc
            u(k2,k1)=uu(k2,k1)
92    continue

c---------------------------------------------------------------------c
c     end of loop over intervals dx                                   c
c---------------------------------------------------------------------c

10    continue

      return
      end



c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine fmunu(ni,zr,e1,e2,x,f,f2,e,b,ee,bb)
c-------------------------------------------------------------------------c
c     evaluate field strength tensor for arbitrary gauge potential given  c
c     subroutine potsu3(). gauge field is assumed to be an antihermitean  c
c     vector field A_\mu=A_\mu^a(\tau^a/2i) with generators normalized    c
c     as tr(\tau^a\tau^b)=2\delta^{ab}.                                   c
c-------------------------------------------------------------------------c
c     ni        number of instantons                                      c
c     zr(i,5)   position, size of instanton i                             c
c     e1(i,6)   real, complex parts of first row of orientation matrix    c
c     e2(i,6)   same for second row                                       c
c     f(a,b,m,n)field strength tensor (a,b=1,3 color; m,n=1,4 vector)     c
c     f2        F^a_{\mu\nu}F^a_{\mu\nu})                                 c
c     e(a,b,i)  electric field (a,b=1,3 color; i=1,3 vector)              c
c     b(a,b,i)  magnetic field                                            c
c     ee,bb     (E^a_i)^2, (B^a_i)^2                                      c
c-------------------------------------------------------------------------c

      parameter(nd=256)
      real zr(nd,5), e1(nd,6), e2(nd,6), x(4), xx(4)
      complex f(3,3,4,4), ff(3,3,4,4)
      complex e(3,3,3), b(3,3,3)
      complex a(3,3,4), aa(3,3,-4:4,4), da(3,3,4)
      complex cf2, cee, cbb
      integer eps(3,2)

      data eps/ 2,3,1, 3,1,2 /

c-------------------------------------------------------------------------c
c     vector potential at x                                               c
c-------------------------------------------------------------------------c

      call potsu3(ni,nd,zr,e1,e2,x,a)

c-------------------------------------------------------------------------c
c     calculate vector potential on star around central point             c
c-------------------------------------------------------------------------c
c     mu labels direction of dx                                           c
c-------------------------------------------------------------------------c

      del = 1.0e-5
      do 10 mu=-4,4
         do 20 nu=1,4
            xx(nu) = x(nu)
 20      continue
         mup = abs(mu)
         xx(mup) = xx(mup) + del*sign(1,mu)
         call potsu3(ni,nd,zr,e1,e2,xx,da)
         do 30 nu=1,4
         do 30 k=1,3
         do 30 l=1,3
            aa(k,l,mu,nu) = da(k,l,nu)
 30      continue
 10   continue

c-------------------------------------------------------------------------c
c     calculate f_{\mu\nu} (not antisymmetrized)                          c
c-------------------------------------------------------------------------c

      do 40 mu=1,4
      do 40 nu=1,4
         do 40 k=1,3
         do 40 l=1,3
            mup = -mu
            ff(k,l,mu,nu) = 1./(2.*del)*(aa(k,l,mu,nu)-aa(k,l,mup,nu))
            do 50 m=1,3
            ff(k,l,mu,nu) = ff(k,l,mu,nu) + a(k,m,mu)*a(m,l,nu)
 50         continue
 40   continue

c-------------------------------------------------------------------------c
c     antisymmetrize                                                      c
c-------------------------------------------------------------------------c

      do 60 mu=1,4
      do 60 nu=1,4
         do 60 k=1,3
         do 60 l=1,3
            f(k,l,mu,nu) = ff(k,l,mu,nu) - ff(k,l,nu,mu)
 60   continue

c-------------------------------------------------------------------------c
c     electric and magnetic fields                                        c
c-------------------------------------------------------------------------c

      do 70 i=1,3
         do 70 k=1,3
         do 70 l=1,3
            e(k,l,i) = f(k,l,4,i)
            b(k,l,i) = f(k,l,eps(i,1),eps(i,2))
 70   continue

c-------------------------------------------------------------------------c
c     f2, e2, b2                                                          c
c-------------------------------------------------------------------------c

      cf2 = 0.0
      cee = 0.0
      cbb = 0.0

      do 80 k=1,3
      do 80 l=1,3
         do 90 mu=1,4
         do 90 nu=mu,4
            cf2 = cf2 + 2.0*f(k,l,mu,nu)*f(l,k,mu,nu)
 90      continue
         do 100 i=1,3
            cee = cee + e(k,l,i)*e(l,k,i)
            cbb = cbb + b(k,l,i)*b(l,k,i)
 100     continue
 80   continue

      f2 =-2.0*real(cf2)
      ee =-2.0*real(cee)
      bb =-2.0*real(cbb)

      return
      end
c----------------------------------------------------------------------+---
c----------------------------------------------------------------------+---

      subroutine ptfull(nin,ndd,zr,e1,e2,xx,rmass,ptfr,sxy)
c---------------------------------------------------------------------------c
c     calculate finite temperature quark propagator in instanton vacuum.    c
c---------------------------------------------------------------------------c
c     in this version the finite temperature non zero mode propagator is    c
c     approximated by nmode terms in a mode sum of the zero t propagator.   c
c     there is no correction for the modified instanton profile in the      c
c     non zero mode part.                                                   c
c---------------------------------------------------------------------------c
c     temperature is controlled by box dimension beta=alb(4) except when    c
c     parameter /cold/ not is set to 1 (then beta=1000) .                   c
c---------------------------------------------------------------------------c
c     parameter iopt controls how the propagator is calculated:             c
c     iopt = 1  full propagator as described above                          c
c     iopt = 2  only zero modes plus finite T free propagator               c
c     iopt = 3  zero modes plus finite T free plus zero T non zero modes    c
c     iopt = 4  free finite T propagator                                    c
c---------------------------------------------------------------------------c
c     nin       number of instantons (nin/2 inst. plus nin/2 antiinst.)     c
c     ndd       max number of instantons                                    c
c     zr(i,5)   position, size of instanton i                               c
c     e1(i,6)   first row of orientation matrix                             c
c     e2(i,6)   second row of orientation matrix                            c
c     xx(4,2)   endpoints of propagator                                     c
c     rmass     current mass                                                c
c     sfr(3,4,.)free propagator                                             c
c     sxy(3,4,.)full propagator                                             c
c---------------------------------------------------------------------------c
c     before calling pfull, the following commons have to be initialized:   c
c     /param/     current masses                                            c
c     /scalmass/  screening mass                                            c
c     /const/     normalization of free propagator                          c
c     /pi/        guess what                                                c
c     /c1c2/      parametrization of stream line overlap m.e.               c
c     /gamma/     gamma matrices, fast gamma multiplication                 c
c     /box/       box dimensions                                            c
c---------------------------------------------------------------------------c
c     for each configuartion need one call to rminv or spect before pfull   c
c     is called for the first time.                                         c
c---------------------------------------------------------------------------c
c     for correct flavor dependence copy clpu (light quark) or clps         c
c     (strange quark) in /clp/ clp before calling pfull with corres-        c
c     ponding rmass.                                                        c
c---------------------------------------------------------------------------c

      parameter(n=3, ni=256, ni2=ni/2, nbin=150, ncf=100, ndl = 10)

      dimension zr(ni,5), e1(ni,6), e2(ni,6), xx(4,2)
      dimension zold(4), znew(4)

      complex clp, clps
      complex sfr(3,4,3,4), sxy(3,4,3,4)
      complex szm(3,4,3,4), snzm(3,4,3,4)
      complex s(3,4,3,4), sfrxy(3,4,3,4)
      complex sm(3,4,3,4), smfr(3,4,3,4)
      complex smunc(3,4,3,4)
      complex stfr(3,4,3,4), dtfr(3,4,3,4), ptfr(3,4,3,4)

c---------------------------------------------------------------------------c
c     common blocks: scalmass gives screening mass, clp overlap m.e.        c
c---------------------------------------------------------------------------c

      common /scalmass/ scalmass
      common /param/ a, alpha,rh0,sg,dz,drh, nc, nf, rmu, rms
      common /clp/   clp(ni,ni), clps(ni,ni)
      common /box/   alb(4)
      common /cold/  not

      iopt = 1
      beta = alb(4)
      if (not .eq. 1) beta = 1000.0

c---------------------------------------------------------------------------c
c     zero mode part of propagator at finite temperature                    c
c---------------------------------------------------------------------------c

      if (iopt .ne. 4) then
      call ztmodes(nc,nin,ndd,zr,e1,e2,xx,rmass,sfr,szm,snorm)
      endif

c---------------------------------------------------------------------------c
c     free finite temperature propagators                                   c
c---------------------------------------------------------------------------c

      if(iopt .ne. 1) then
c     call ptfree(xx,rmass,ptfr)
c     call dtfree(xx,rmass,dtfr)
      call stfree(xx,rmass,stfr)
      endif

c---------------------------------------------------------------------------c
c     non zero mode part at zero temperature                                c
c---------------------------------------------------------------------------c

      if (iopt .eq. 1 .or. iopt .eq. 3) then
      rs = 0.0
      do 110 m = 1, 3
         zold(m) = xx(m,1)
         znew(m) = xx(m,2)
         rs = rs + (zold(m)-znew(m))**2
  110 continue
      zold(4) = xx(4,1)
      znew(4) = xx(4,2)
      zs = rs + (znew(4)-zold(4))**2
      zs = sqrt(zs)

      call prtnz(nin,zold,znew,zr,s,sfrxy,sm,smfr,rmass,smunc)

      endif

c---------------------------------------------------------------------------c
c     mass corrections for free scalar and vector part of spin 1/2 prop.    c
c---------------------------------------------------------------------------c

      if (iopt .eq. 1) then

      call masscorr(rmass,zs,corg,corm)

c---------------------------------------------------------------------------c
c     replace free prop in S and SM=im*D by massive counterparts            c
c---------------------------------------------------------------------------c

      do 1121 k1 = 1, 3
      do 1121 l1 = 1, 4
      do 1121 k2 = 1, 3
      do 1121 l2 = 1, 4
         s(k1,l1,k2,l2)  = s(k1,l1,k2,l2)
     1                  - (1-corg)*sfrxy(k1,l1,k2,l2)
         sm(k1,l1,k2,l2) = sm(k1,l1,k2,l2)-(1-corm)*smfr(k1,l1,k2,l2)
         snzm(k1,l1,k2,l2) = s(k1,l1,k2,l2) + sm(k1,l1,k2,l2)
 1121 continue

c---------------------------------------------------------------------------c
c     mode sum over zero T background field propagator                      c
c---------------------------------------------------------------------------c
c     note : nmodes should be even, nm=-1,+1,-2,+2,.. counts modes          c
c---------------------------------------------------------------------------c

      nmodes = 6

      do 300 i=1,nmodes
         i2 = int(i/2.0+0.75)
         nm = i2*(-1)**i
         znew(4) = xx(4,2)+nm*beta
         zs = rs + (znew(4)-zold(4))**2
         zs = sqrt(zs)

         call prtnz(nin,zold,znew,zr,s,sfrxy,sm,smfr,rmass,smunc)
         call masscorr(rmass,zs,corg,corm)

         do 300 k1=1,3
         do 300 l1=1,4
         do 300 k2=1,3
         do 300 l2=1,4

         s(k1,l1,k2,l2)  = s(k1,l1,k2,l2)
     1                  - (1-corg)*sfrxy(k1,l1,k2,l2)
         sm(k1,l1,k2,l2) = sm(k1,l1,k2,l2)-(1-corm)*smfr(k1,l1,k2,l2)

c---------------------------------------------------------------------------c
c     add contribution from mode nm                                         c
c---------------------------------------------------------------------------c

         snzm(k1,l1,k2,l2) = snzm(k1,l1,k2,l2)
     1       + (s(k1,l1,k2,l2) + sm(k1,l1,k2,l2))*(-1)**abs(nm)

 300  continue
      endif

c---------------------------------------------------------------------------c
c     construct full propagator; different models:                          c
c---------------------------------------------------------------------------c
c     full prop = zero modes + non zero modes                               c
c---------------------------------------------------------------------------c

      if (iopt .eq. 1) then
      do 120 k1 = 1, 3
      do 120 l1 = 1, 4
      do 120 k2 = 1, 3
      do 120 l2 = 1, 4
          sxy(k1,l1,k2,l2) = szm(k1,l1,k2,l2)+ snzm(k1,l1,k2,l2)
  120 continue

c---------------------------------------------------------------------------c
c     simpler model: zero modes plus free part                              c
c---------------------------------------------------------------------------c

      else if(iopt .eq. 2) then
      do 130 k1 = 1, 3
      do 130 l1 = 1, 4
      do 130 k2 = 1, 3
      do 130 l2 = 1, 4
          sxy(k1,l1,k2,l2) = szm(k1,l1,k2,l2)+stfr(k1,l1,k2,l2)
  130 continue

c---------------------------------------------------------------------------c
c     somewhat better: zero modes plus free plus zero T non zero modes      c
c---------------------------------------------------------------------------c

      else if (iopt .eq. 3) then
      do 140 k1 = 1, 3
      do 140 l1 = 1, 4
      do 140 k2 = 1, 3
      do 140 l2 = 1, 4
          sxy(k1,l1,k2,l2) = szm(k1,l1,k2,l2) + stfr(k1,l1,k2,l2)
     1           + s(k1,l1,k2,l2)-sfrxy(k1,l1,k2,l2)
     2           +sm(k1,l1,k2,l2)- smfr(k1,l1,k2,l2)
  140 continue

c---------------------------------------------------------------------------c
c     even simpler: free finite T propagator                                c
c---------------------------------------------------------------------------c

      else if (iopt .eq. 4) then
      do 150 k1 = 1, 3
      do 150 l1 = 1, 4
      do 150 k2 = 1, 3
      do 150 l2 = 1, 4
          sxy(k1,l1,k2,l2) = stfr(k1,l1,k2,l2)
  150 continue
      endif

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine ztmodes(n,nin,nd,zr,e1,e2,xx,rms,sfr,szm,snorm)
c---------------------------------------------------------------------------c
c     zero mode part of propagator at finite temperature.                   c
c---------------------------------------------------------------------------c
c     temperature is given by inverse box size in 4-direction (as specified c
c     in common /box/ alb(4)). propagator is not periodic in the other      c
c     directions.                                                           c
c---------------------------------------------------------------------------c
c     n         number of colors                                            c
c     nin       number of instantons (nin/2 inst. plus nin/2 antiinst.)     c
c     nd        max number of instantons                                    c
c     zr(i,5)   position, size of instanton i                               c
c     e1(i,6)   first row of orientation matrix                             c
c     e2(i,6)   second row of orientation matrix                            c
c     xx(4,2)   endpoints of propagator                                     c
c     rmass     current mass                                                c
c     sfr(3,4,.)free propagator                                             c
c     szm(3,4,.)zero mode propagator                                        c
c     snorm     not used ?                                                  c
c---------------------------------------------------------------------------c
c     im = 1    calculate profiles analytically                             c
c     im = 2    calculate derivative of scalar profile numerically          c
c     im = 3    zero temperature profile                                    c
c---------------------------------------------------------------------------c

      parameter(n2=256)
      parameter(im=1)
      complex u1, clm, gamma, sfr, clap, clps
      complex phi, szm, epu, gmx, clp, cphi, ccphi
      complex g, gm5
      dimension phi(n2,3,4,2), zr(nd,5),e1(nd,6),e2(nd,6)
      dimension si(n2), dis(n2,4), xx(4,2), cphi(n2,3,4), ccphi(n2,3,4)
      dimension sfr(3,4,3,4), szm(3,4,3,4), epu(n2,3,4), gmx(n2,4,4)
      dimension r2(n2), scfact(n2), scfacr(n2), rnor(n2)

c---------------------------------------------------------------------------c
c     common /clp/ contains overlap matrix elements calculated in rminv     c
c---------------------------------------------------------------------------c
c     common /u1/  contains orientation matrix (also from rminv)            c
c---------------------------------------------------------------------------c

      common /gamf/ g(5,4), nu(5,4), ni(5,4), c(4), nuc(4), nic(4)
      common /gamf5/ gm5(4,4), nm5(4,4)
      common /clap/ clap(n2,n2)
      common /param/ a, alpha
      common /pi/ pi
      common /gam/ gamma(5,4,4), gr(4,4), gl(4,4), cccc(4,4)
      common /u1/ u1(n2,3,3)
      common /clp/ clp(n2,n2), clps(n2,n2)
      common /box/ alb(4)
      common /cold/  not

c---------------------------------------------------------------------------c
c     si=-/+ for instanton/antiinstanton                                    c
c---------------------------------------------------------------------------c

      nih = nin/2
      ahl = a/2
      beta = alb(4)
      if (not .eq. 1) beta=1000.0
      t = 1.0/beta

      coef = 1.0/(2.0*sqrt(2.0)*pi)
      dtau = 1.0e-5
      dr   = 1.0e-5

      do 5 j = 1, nin
        si(j) = sign(1.0,j-nih-0.5)
  5   continue
      if (nin.eq.1) then
        si(1) = -1.0
      end if

c---------------------------------------------------------------------------c
c     ix=1,2 for |x> and <y|, construct wavefuntions \phi(x),\phi(y)        c
c---------------------------------------------------------------------------c

      do 100 ix = 1, 2

c---------------------------------------------------------------------------c
c     dis(j,mu)=(x-z_j)_\mu along shortest path                             c
c---------------------------------------------------------------------------c

        call zero(nin,r2)

        do 20 ir = 1, 3
           do 910 j = 1, nin
              ds = xx(ir,ix)-zr(j,ir)
              ads = abs(ds)
              asg = 1+sign(1.0,ads-alb(ir)/2)
              dis(j,ir) = ds - alb(ir)/2*asg*sign(1.0,ds)
              r2(j) = r2(j) + dis(j,ir)**2
  910      continue
   20   continue

c---------------------------------------------------------------------------c
c    along 4 direction, wavefct. takes care of antiperiodic boundary cond.  c
c---------------------------------------------------------------------------c

        if (not .eq. 1) then
        do 915 j = 1, nin
           ds = xx(4,ix)-zr(j,4)
           ads = abs(ds)
           asg = 1+sign(1.0,ads-alb(4)/2)
           dis(j,4) = ds - alb(4)/2*asg*sign(1.0,ds)
  915   continue
        else
        do 916 j = 1, nin
           ds = xx(4,ix)-zr(j,4)
           dis(j,4) = ds
  916   continue
        endif

c---------------------------------------------------------------------------c
c     zero mode profile                                                     c
c---------------------------------------------------------------------------c

        rnorm = 0.0
        do 920 j = 1, nin
           rho = zr(j,5)
           r   = sqrt(r2(j))
           tau = dis(j,4)
           fact= coef/rho*sqrt( 1 + Pi*rho**2*t*Sinh(2*Pi*r*t)/
     1       (r*(-Cos(2*Pi*t*tau) + Cosh(2*Pi*r*t))) )

           if (im .eq. 1) then

           scfact(j) =-fact * (
     1     -2*Pi**2*rho**2*t**2*Sin(Pi*t*tau)*Sinh(Pi*r*t)*
     2     (2*r + r*Cos(2*Pi*t*tau) + r*Cosh(2*Pi*r*t) +
     3      Pi*rho**2*t*Sinh(2*Pi*r*t))
     4      /(r*Cos(2*Pi*t*tau) - r*Cosh(2*Pi*r*t) -
     5      Pi*rho**2*t*Sinh(2*Pi*r*t))**2 )

           scfacr(j) =-fact * (
     1      Pi*rho**2*t*Cos(Pi*t*tau)*(3*Pi*r*t*Cosh(Pi*r*t) -
     2      2*Pi*r*t*Cos(2*Pi*t*tau)*Cosh(Pi*r*t) -
     3      Pi*r*t*Cosh(3*Pi*r*t) +
     4      Sinh(Pi*r*t) + 3*Pi**2*rho**2*t**2*Sinh(Pi*r*t) +
     5      2*Cos(2*Pi*t*tau)*Sinh(Pi*r*t) - Sinh(3*Pi*r*t) -
     6      Pi**2*rho**2*t**2*Sinh(3*Pi*r*t))/
     7     (r*Cos(2*Pi*t*tau) - r*Cosh(2*Pi*r*t) -
     8      Pi*rho**2*t*Sinh(2*Pi*r*t))**2 )

c---------------------------------------------------------------------------c
c     or : calculate derivatives numerically                                c
c---------------------------------------------------------------------------c

            else if (im .eq. 2) then

            tau1 = tau + dtau
            tau2 = tau - dtau
            rr1  = r + dr
            rr2  = r - dr

            f1 = 2*Pi*rho**2*t*Cos(Pi*t*tau)*Sinh(Pi*rr1*t)/
     1        (-(rr1*Cos(2*Pi*t*tau)) + rr1*Cosh(2*Pi*rr1*t)
     2         + Pi*rho**2*t*Sinh(2*Pi*rr1*t))
            f2 = 2*Pi*rho**2*t*Cos(Pi*t*tau)*Sinh(Pi*rr2*t)/
     1        (-(rr2*Cos(2*Pi*t*tau)) + rr2*Cosh(2*Pi*rr2*t)
     2         + Pi*rho**2*t*Sinh(2*Pi*rr2*t))
            f3 = 2*Pi*rho**2*t*Cos(Pi*t*tau1)*Sinh(Pi*r*t)/
     1        (-(r*Cos(2*Pi*t*tau1)) + r*Cosh(2*Pi*r*t)
     2         + Pi*rho**2*t*Sinh(2*Pi*r*t))
            f4 = 2*Pi*rho**2*t*Cos(Pi*t*tau2)*Sinh(Pi*r*t)/
     1        (-(r*Cos(2*Pi*t*tau2)) + r*Cosh(2*Pi*r*t)
     2         + Pi*rho**2*t*Sinh(2*Pi*r*t))
            scfacr(j) =-fact * (f1-f2)/(2.0*dr)
            scfact(j) =-fact * (f3-f4)/(2.0*dtau)

c---------------------------------------------------------------------------c
c     zero temperature profile                                              c
c---------------------------------------------------------------------------c

            else if (im .eq. 3) then

            x2 = r**2+tau**2
            ff = coef*2.0*rho/sqrt(x2)/(x2+rho**2)**1.5
            scfacr(j) = ff * r
            scfact(j) = ff *tau

            endif

            rnor(j) = scfacr(j)*scfacr(j)*r2(j)
  920   continue

        if (ix.eq.1) then
           rnorm = 0.0
           do 925 j = 1, nin
              rnorm = rnorm + rnor(j)
  925      continue
           snorm = rnorm*4
        end if

c---------------------------------------------------------------------------c
c       gmx(j,k,l) = (\gamma_4)_{kl}\del_4 F(x,z_j) + (\gamma_i)_{kl}*..    c
c---------------------------------------------------------------------------c

        do 36 k = 1, 4
           do 931 l = 1, 4
           do 930 j = 1, nin
              gmx(j,k,l) = cmplx(0.0,0.0)
  930      continue
  931      continue
           do 37 m = 1, 3
           do 940 j = 1, nin
              gmx(j,k,nu(m,k)) = gmx(j,k,nu(m,k))
     1          + dis(j,m)*g(m,k)/sqrt(r2(j))*scfacr(j)
  940      continue
   37      continue
           do 945 j = 1, nin
              gmx(j,k,nu(4,k)) = gmx(j,k,nu(4,k))
     1          + g(4,k)*scfact(j)
  945      continue
   36   continue

c---------------------------------------------------------------------------c
c       spin-color: epu(j,k,m)=(U^j)_{ka}\Omega_{am}                        c
c---------------------------------------------------------------------------c

        do 50 k = 1, 3
           do 950 j = 1, nin
              epu(j,k,1) = u1(j,k,2)
              epu(j,k,2) = -u1(j,k,1)
              epu(j,k,3) = u1(j,k,2)*si(j)
              epu(j,k,4) = -u1(j,k,1)*si(j)
  950      continue
   50   continue

c---------------------------------------------------------------------------c
c       wavefunction phi(j,k,l,ix) = (\phi_J)^k_l  k:color l:spin           c
c---------------------------------------------------------------------------c

        do 60 k = 1, 3
        do 60 l = 1, 4
            do 960 j = 1, nin
              phi(j,k,l,ix) =
     1            gmx(j,l,1)*epu(j,k,1) + gmx(j,l,2)*epu(j,k,2)
     2          + gmx(j,l,3)*epu(j,k,3) + gmx(j,l,4)*epu(j,k,4)
  960       continue
   60   continue

c---------------------------------------------------------------------------c
c     end of loop over ix (|x>,<y|)                                         c
c---------------------------------------------------------------------------c

  100 continue

c---------------------------------------------------------------------------c
c     adjoint of phi(y)                                                     c
c---------------------------------------------------------------------------c

      do 460 k2 = 1, 3
      do 460 l2 = 1, 4
        do 460 i = 1, nin
          ccphi(i,k2,l2) = conjg(phi(i,k2,l2,2))
 460  continue

c---------------------------------------------------------------------------c
c     rotate \phi_I --> (T+im)^(-1)_{IJ} \phi_J                             c
c---------------------------------------------------------------------------c

      do 300 k2 = 1, 3
      do 300 l2 = 1, 4
        do 410 i = 1, nin
          cphi(i,k2,l2) = cmplx(0.0,0.0)
          do 410 m = 1, nin
          cphi(i,k2,l2) = cphi(i,k2,l2) + clp(i,m)*ccphi(m,k2,l2)
 410  continue

c---------------------------------------------------------------------------c
c     propagator  S(x,y)=\phi_I(x) S_{IJ} \phi_J(y)^\dagger                 c
c---------------------------------------------------------------------------c

      do 302 k1 = 1, 3
      do 302 l1 = 1, 4
        szm(k1,l1,k2,l2) = cmplx(0.0,0.0)
  302 continue
      do 400 j = 1, nin
        szm(1,1,k2,l2) = szm(1,1,k2,l2)+phi(j,1,1,1)*cphi(j,k2,l2)
        szm(1,2,k2,l2) = szm(1,2,k2,l2)+phi(j,1,2,1)*cphi(j,k2,l2)
        szm(1,3,k2,l2) = szm(1,3,k2,l2)+phi(j,1,3,1)*cphi(j,k2,l2)
        szm(1,4,k2,l2) = szm(1,4,k2,l2)+phi(j,1,4,1)*cphi(j,k2,l2)
        szm(2,1,k2,l2) = szm(2,1,k2,l2)+phi(j,2,1,1)*cphi(j,k2,l2)
        szm(2,2,k2,l2) = szm(2,2,k2,l2)+phi(j,2,2,1)*cphi(j,k2,l2)
        szm(2,3,k2,l2) = szm(2,3,k2,l2)+phi(j,2,3,1)*cphi(j,k2,l2)
        szm(2,4,k2,l2) = szm(2,4,k2,l2)+phi(j,2,4,1)*cphi(j,k2,l2)
        szm(3,1,k2,l2) = szm(3,1,k2,l2)+phi(j,3,1,1)*cphi(j,k2,l2)
        szm(3,2,k2,l2) = szm(3,2,k2,l2)+phi(j,3,2,1)*cphi(j,k2,l2)
        szm(3,3,k2,l2) = szm(3,3,k2,l2)+phi(j,3,3,1)*cphi(j,k2,l2)
        szm(3,4,k2,l2) = szm(3,4,k2,l2)+phi(j,3,4,1)*cphi(j,k2,l2)
  400 continue
  300 continue



      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine prtnz(nin,zold,znew,zr,s,sfree,sm,smfr,amass,smunc)
c----------------------------------------------------------------------------c
c     calculation of the propagator due to nonzero modes                     c
c----------------------------------------------------------------------------c
c     nin         number of instantons                                       c
c     zold(4)     x(4)                                                       c
c     znew(4)     y(4)                                                       c
c     zr(ni,5)    collective coordinates                                     c
c     s(3,4,3,4)  non zero mode propagator                                   c
c     sfree(3,.)  free propagator                                            c
c     sm(3,4,3,)  linear mass correction m*D(x,y)                            c
c     smfr(3,4,)  same for free propagaror m*D_0(x,y)                        c
c     amass       current quark mass                                         c
c     smunc(3,.)                                                             c
c----------------------------------------------------------------------------c

      parameter(ni=256)

      complex s(3,4,3,4), sfree(3,4,3,4),uuu,cccl,cccr,corr
      complex sm(3,4,3,4),sxy(ni,2,2)
      complex smfr(3,4,3,4), dm(ni,3,3),ddm(ni,3,3),dddm(ni,3,3)
      complex dmunc(ni,3,3), smunc(3,4,3,4)
      complex gamma,tau
      complex ucross(ni,3,3), cfactor(ni)
      complex ul(ni,3,3),ur(ni,3,3)
      complex uul(ni,3,3),uur(ni,3,3)
      complex uuul(ni,3,3),uuur(ni,3,3)
      complex al(ni,2,2),ar(ni,2,2)
      complex bl(ni,2,2),br(ni,2,2)
      complex bbl(ni,2,2), bbr(ni,2,2)
      complex cl(ni,2,2),cr(ni,2,2)
      complex dl(ni,2,2),dr(ni,2,2)
      complex ddl(ni,2,2),ddr(ni,2,2)
      complex uup(ni,3,3), uum(ni,3,3)
      complex sx(ni,2,2),sy(ni,2,2), sdell(ni,2,2), sdelr(ni,2,2)
      complex gamplu,gammin,delhat(4,4), uu

c---------------------------------------------------------------------------c
c     we assume that gamma was already called by the master !!!             c
c---------------------------------------------------------------------------c

      common /gam/ gamma(5,4,4), gr(4,4), gl(4,4), cccc(4,4)
      complex g, gm5
      common /gamf5/ gm5(4,4), nm5(4,4)
      common /gamf/ g(5,4), nu(5,4), nv(5,4), c(4), nuc(4), nic(4)
      common /taum/tau(2,4,2,2)
      common /lr/ gamplu(4,4,4),gammin(4,4,4)
      common /param/ a, alpha, rh0
      common /box/ alb(4)
      common /const/ const
      common /u1/ uu(ni,3,3)
      dimension znew(4),zold(4),delta(4),z(5), uuu(3,3)
      dimension xx(ni,4), yy(ni,4), xs(ni), ys(ni)
      dimension zr(ni,5),e1(ni,6),e2(ni,6)
      dimension xr(ni), yr(ni), hx(ni), hy(ni), factor(ni)
      dimension factm(ni), delt(ni,4), ros(ni), fact(ni)
      dimension corr(ni), indmx(ni), indsg(ni), fakx(ni), faky(ni)
      complex ss(3,4,3,4)
      dimension sccor(ni)

      ndd = ni
      nih = nin/2

c---------------------------------------------------------------------------c
c     indmx = 1,2  indsg=+1,-1  for inst/antiinstanton                      c
c---------------------------------------------------------------------------c

      do 820 is = 1, nin
        indmx(is) = (is-1)/nih + 1
        indsg(is) = 3-2*indmx(is)
  820 continue

c---------------------------------------------------------------------------c
c     calculate shortest distance around torus                              c
c---------------------------------------------------------------------------c
c     asg = 0,2 if shortest path is direct/backwards                        c
c---------------------------------------------------------------------------c

      dis = zold(4) -znew(4)
      delta(4) = dis
c     adis=abs(dis)
c     asg = 1+sign(1.0,adis-alb(4)/2)
c     delta(4) = dis - alb(4)/2*asg*sign(1.0,dis)

      dis = zold(3) -znew(3)
      adis=abs(dis)
      asg = 1+sign(1.0,adis-alb(3)/2)
      delta(3) = dis - alb(3)/2*asg*sign(1.0,dis)

      dis = zold(2) -znew(2)
      adis=abs(dis)
      asg = 1+sign(1.0,adis-alb(2)/2)
      delta(2) = dis - alb(2)/2*asg*sign(1.0,dis)

      dis = zold(1) -znew(1)
      adis=abs(dis)
      asg = 1+sign(1.0,adis-alb(1)/2)
      delta(1) = dis - alb(1)/2*asg*sign(1.0,dis)

      deltas = delta(1)*delta(1) + delta(2)*delta(2)
     1       + delta(3)*delta(3) + delta(4)*delta(4)

c---------------------------------------------------------------------------c
c     delhat = gamma_\mu (x-y)_\mu                                          c
c---------------------------------------------------------------------------c

      do 21  m1= 1,4
      do 21  m2= 1,4
        delhat(m1,m2)=gamma(1,m1,m2)*delta(1)
     *              + gamma(2,m1,m2)*delta(2)
     *              + gamma(3,m1,m2)*delta(3)
     *              + gamma(4,m1,m2)*delta(4)
21    continue

c---------------------------------------------------------------------------c
c     clear propagators                                                     c
c---------------------------------------------------------------------------c

      do 1  m1= 1,4
      do 1  m2= 1,4
      do 1  k1= 1,3
      do 1  k2= 1,3
        sfree(k1,m1,k2,m2)=(0.,0.)
        s(k1,m1,k2,m2)=(0.,0.)
        sm(k1,m1,k2,m2)=(0.,0.)
        smunc(k1,m1,k2,m2)=(0.,0.)
        smfr(k1,m1,k2,m2) = cmplx(0.0,0.0)
1     continue

c---------------------------------------------------------------------------c
c     free propagators                                                      c
c---------------------------------------------------------------------------c

      do 2 k1 = 1, 3
      do 2 m1 = 1, 4
        smfr(k1,m1,k1,m1) = amass*const/deltas*cmplx(0.0,1.0)/2
      do 2 m2 = 1, 4
        sfree(k1,m1,k1,m2)=-const*delhat(m1,m2)/deltas**2
     2             *cmplx(0.0,-1.0)
  2   continue

c---------------------------------------------------------------------------c
c     sweep over all instantons                                             c
c---------------------------------------------------------------------------c
c     xx(i,4) = (x-z_i)_\mu ; xs(i) = xx(i)**2                              c
c---------------------------------------------------------------------------c

      call tordis3(nin,zold,zr,xx,xs)
      call tordis3(nin,znew,zr,yy,ys)

c---------------------------------------------------------------------------c
c     sx(i,2,2) = sig(-)_{ij}\cdot (x-z_i) ; sy = sig(+)\cdot (y-z_i)       c
c---------------------------------------------------------------------------c
c     note : sigvec reverses index for antiinstanton                        c
c---------------------------------------------------------------------------c

      call sigvec(nin,xx,-1,sx)
      call sigvec(nin,yy,1,sy)

c---------------------------------------------------------------------------c
c     screening : fakx(i)= D_M(x,z_i)/D_0(x,z_i)  faky(i) = ...             c
c---------------------------------------------------------------------------c

      call scalcor(nin,xs,fakx)
      call scalcor(nin,ys,faky)

c---------------------------------------------------------------------------c
c     \sum_i (fakx(i----)*faky(i)) , needed to uncorrect free propagator        c
c---------------------------------------------------------------------------c

      corfree = 0.0
      do 509 k = 1, nin
        corfree = corfree + fakx(k)*faky(k)
  509 continue

c---------------------------------------------------------------------------c
c     various factors appearing in scal prop D_I                            c
c-----------------------------------------------------------------------c

      do 510 is = 1, nin
         ros(is) = zr(is,5)*zr(is,5)
         xr(is)  =    xs(is)+ros(is)
         yr(is)  =    ys(is)+ros(is)
         hx(is)  = 1./sqrt(1.+ros(is)/xs(is))
         hy(is)  = 1./sqrt(1.+ros(is)/ys(is))
         factor(is) =-hx(is)*hy(is)*const/deltas**2
         cfactor(is)= factor(is)*cmplx(0.0,-1.0)*fakx(is)*faky(is)
         factm(is)  = const*0.5/deltas*amass*fakx(is)*faky(is)
  510 continue

c---------------------------------------------------------------------------c
c     orientation uu determined in main, ucross = uu^\dagger                c
c---------------------------------------------------------------------------c

      do 10 k1=1,3
      do 10 k2=1,3
         do 520 is = 1, nin
            ucross(is,k2,k1)=conjg(uu(is,k1,k2))
  520    continue
   10 continue

c----------------------------------------------------------------------------c
c     sxy(i,2,2) = sig(-)*(x-z_i) sig(+)*(y-z_i)                             c
c----------------------------------------------------------------------------c

      call mult2(nin,sx,sy,sxy)

c----------------------------------------------------------------------------c
c     embed scalar propagator into su(3); dm(i,3,3) = D_I(3,3)               c
c----------------------------------------------------------------------------c
c     note: dm does not include 1/(x-y)^2                                    c
c----------------------------------------------------------------------------c

      do 530 is = 1, nin
         fact(is) = ros(is)/(xs(is)*ys(is))*hx(is)*hy(is)
         sccor(is) = fact(is)*deltas/2*(1-fakx(is)*faky(is))
         sccor(is) = 0.0
         dm(is,1,1)=sxy(is,1,1)*fact(is) + hx(is)*hy(is)
     2                     + sccor(is)
         dm(is,1,2)=sxy(is,1,2)*fact(is)
         dm(is,2,1)=sxy(is,2,1)*fact(is)
         dm(is,2,2)=sxy(is,2,2)*fact(is) + hx(is)*hy(is)
     2                     + sccor(is)
         dm(is,1,3)=(0.,0.)
         dm(is,2,3)=(0.,0.)
         dm(is,3,3)=(1.,0.)
         dm(is,3,2)=(0.,0.)
         dm(is,3,1)=(0.,0.)
  530 continue

c----------------------------------------------------------------------------c
c     orientation : D_I = U*D_I*U^+                                          c
c----------------------------------------------------------------------------c

      call mult3(nin,uu,dm,ddm)
      call mult3(nin,ddm,ucross,dddm)

c----------------------------------------------------------------------------c
c     not used ?                                                             c
c----------------------------------------------------------------------------c

      do 531 k1 = 1, 3
      do 531 k2 = 1, 3
      do 531 is = 1, nin
         dmunc(is,k1,k2) = dddm(is,k1,k2)
  531 continue

      do 532 is = 1, nin
         dmunc(is,1,1) = dddm(is,1,1)-sccor(is)
         dmunc(is,2,2) = dddm(is,2,2)-sccor(is)
  532 continue

c----------------------------------------------------------------------------c
c     sm = i*m*\sum_I D_I ; including screening corrections                  c
c----------------------------------------------------------------------------c

      do 31 k1=1,3
      do 31 k2=1,3
      do 31 m1=1,4
         do 540 is = 1, nin
            sm(k1,m1,k2,m1) = sm(k1,m1,k2,m1) +
     1             dddm(is,k1,k2) * factm(is)*cmplx(0.0,1.0)
            smunc(k1,m1,k2,m1) = smunc(k1,m1,k2,m1) +
     1             dmunc(is,k1,k2)* factm(is)*cmplx(0.0,1.0)
 540     continue
 31   continue

c----------------------------------------------------------------------------c
c     sm = sm0 + i*m*\sum_I (D_I-D_0)                                        c
c----------------------------------------------------------------------------c

      do 650 k1 = 1, 3
      do 650 m1 = 1, 4
        sm(k1,m1,k1,m1) = sm(k1,m1,k1,m1)+smfr(k1,m1,k1,m1)
     2                   - corfree*smfr(k1,m1,k1,m1)
        smunc(k1,m1,k1,m1) = smunc(k1,m1,k1,m1)
     2                       -(nin-1)*smfr(k1,m1,k1,m1)
 650  continue

c----------------------------------------------------------------------------c
c     scalar propagator finished                                             c
c----------------------------------------------------------------------------c

c----------------------------------------------------------------------------c
c     decompose  S = (sl_\mu \gam_5^+ sr_\mu \gam_5^-) \gamma_mu             c
c----------------------------------------------------------------------------c
c     loop over mu, calculate sl(mu) and sr(mu)                              c
c----------------------------------------------------------------------------c

      do 5 mu=1,4

c----------------------------------------------------------------------------c
c     initialize al/ar= sig^(+/-)_\mu                                        c
c----------------------------------------------------------------------------c

         do 550 is = 1, nin
            indmax = indmx(is)
            al(is,1,1)=tau(indmax,mu,1,1)
            al(is,1,2)=tau(indmax,mu,1,2)
            al(is,2,1)=tau(indmax,mu,2,1)
            al(is,2,2)=tau(indmax,mu,2,2)
            ar(is,1,1)=tau(3-indmax,mu,1,1)
            ar(is,1,2)=tau(3-indmax,mu,1,2)
            ar(is,2,1)=tau(3-indmax,mu,2,1)
            ar(is,2,2)=tau(3-indmax,mu,2,2)
  550    continue

         do 560 is = 1, nin
            delt(is,1) = delta(1)
            delt(is,2) = delta(2)
            delt(is,3) = delta(3)
            delt(is,4) = delta(4)
  560    continue

c----------------------------------------------------------------------------c
c     sdel = (x-y)\cdot\sig^(+/-)                                            c
c----------------------------------------------------------------------------c

         call sigvec(nin, delt,-1,sdell)
         call sigvec(nin, delt, 1,sdelr)

c----------------------------------------------------------------------------c
c     bl = \sig_\mu * sdel , br = ...                                        c
c----------------------------------------------------------------------------c

         call mult2(nin,al,sdell,bl)
         call mult2(nin,sdelr,ar,br)

c----------------------------------------------------------------------------c
c     structure inside \sig\dot x ( ... ) \sig\dot y                         c
c----------------------------------------------------------------------------c

         do 580 is = 1, nin
            bl(is,1,1) = bl(is,1,1)*deltas*0.5/xr(is)
     2                  + cmplx(delta(mu),0.0)
            bl(is,1,2) = bl(is,1,2)*deltas*0.5/xr(is)
            bl(is,2,1) = bl(is,2,1)*deltas*0.5/xr(is)
            bl(is,2,2) = bl(is,2,2)*deltas*0.5/xr(is)
     2                  + cmplx(delta(mu),0.0)
            br(is,1,1) = br(is,1,1)*deltas*0.5/yr(is)
     2                  + cmplx(delta(mu),0.0)
            br(is,1,2) = br(is,1,2)*deltas*0.5/yr(is)
            br(is,2,1) = br(is,2,1)*deltas*0.5/yr(is)
            br(is,2,2) = br(is,2,2)*deltas*0.5/yr(is)
     2                  + cmplx(delta(mu),0.0)
  580    continue

c----------------------------------------------------------------------------c
c     multiply by \sig\dot x,\sig\dot y                                      c
c----------------------------------------------------------------------------c

         call mult2(nin,sx,bl,cl)
         call mult2(nin,sx,br,cr)
         call mult2(nin,cl,sy,dl)
         call mult2(nin,cr,sy,dr)

c----------------------------------------------------------------------------c
c     add free propagator                                                    c
c----------------------------------------------------------------------------c

         do 590 is = 1, nin
            roxy=ros(is)/(xs(is)*ys(is))
            dl(is,1,1) = dl(is,1,1)*roxy+cmplx(delta(mu),0.0)
            dl(is,1,2) = dl(is,1,2)*roxy
            dl(is,2,1) = dl(is,2,1)*roxy
            dl(is,2,2) = dl(is,2,2)*roxy+cmplx(delta(mu),0.0)
            dr(is,1,1) = dr(is,1,1)*roxy+cmplx(delta(mu),0.0)
            dr(is,1,2) = dr(is,1,2)*roxy
            dr(is,2,1) = dr(is,2,1)*roxy
            dr(is,2,2) = dr(is,2,2)*roxy+cmplx(delta(mu),0.0)
  590    continue

c----------------------------------------------------------------------------c
c     embed in su(3)                                                         c
c----------------------------------------------------------------------------c

         do 600 is = 1, nin
            ul(is,1,1) = dl(is,1,1)
            ul(is,1,2) = dl(is,1,2)
            ul(is,2,1) = dl(is,2,1)
            ul(is,2,2) = dl(is,2,2)
            ur(is,1,1) = dr(is,1,1)
            ur(is,1,2) = dr(is,1,2)
            ur(is,2,1) = dr(is,2,1)
            ur(is,2,2) = dr(is,2,2)
  600    continue

c----------------------------------------------------------------------------c
c     free part should note be multiplied by hx,hy; correct in advance       c
c----------------------------------------------------------------------------c

         do 610 is = 1, nin
            ur(is,3,3)=cmplx(delta(mu)/(hx(is)*hy(is)),0.0)
            ul(is,3,3)=cmplx(delta(mu)/(hx(is)*hy(is)),0.0)
            ul(is,1,3)=(0.,0.)
            ul(is,2,3)=(0.,0.)
            ul(is,3,2)=(0.,0.)
            ul(is,3,1)=(0.,0.)
            ur(is,1,3)=(0.,0.)
            ur(is,2,3)=(0.,0.)
            ur(is,3,2)=(0.,0.)
            ur(is,3,1)=(0.,0.)
  610    continue

c----------------------------------------------------------------------------c
c     orientation                                                            c
c----------------------------------------------------------------------------c

         call mult3(nin,uu,ul,uul)
         call mult3(nin,uu,ur,uur)
         call mult3(nin,uul,ucross,uuul)
         call mult3(nin,uur,ucross,uuur)

c----------------------------------------------------------------------------c
c     go from \gamma_5^(+/-) to (1,\gamma_5) decomposition                   c
c----------------------------------------------------------------------------c
c     note that antiinstanton has opposite chirality                         c
c----------------------------------------------------------------------------c

         do 750 k1 = 1, 3
         do 750 k2 = 1, 3
            do 760 is = 1, nin
                  isg = indsg(is)
                  uup(is,k1,k2) = (uuul(is,k1,k2)+uuur(is,k1,k2))/2
                  uum(is,k1,k2) = (uuul(is,k1,k2)-uuur(is,k1,k2))/2
     2                            *isg
  760       continue
  750    continue

c---------------------------------------------------------------------------c
c     contract with \gamma_\mu                                              c
c---------------------------------------------------------------------------c

         do 9 m1=1,4
         do 9 k1=1,3
         do 9 k2=1,3
            do 770 is = 1, nin
               s(k1,m1,k2,nu(mu,m1)) = s(k1,m1,k2,nu(mu,m1)) +
     2             uup(is,k1,k2)*g(mu,m1)*cfactor(is)
               s(k1,m1,k2,nm5(mu,m1)) = s(k1,m1,k2,nm5(mu,m1))
     2          +  uum(is,k1,k2)*gm5(mu,m1)*cfactor(is)
  770    continue
    9 continue

c----------------------------------------------------------------------------c
c     this is the  end of the long loop over mu                              c
c----------------------------------------------------------------------------c

5     continue

c----------------------------------------------------------------------------c
c     uncorrect free part for screening                                      c
c----------------------------------------------------------------------------c

      do 640 k1 = 1,3
      do 640 k2 = 1,3
      do 640 m1 = 1,4
      do 640 m2 = 1,4
         s(k1,m1,k2,m2) = s(k1,m1,k2,m2) + sfree(k1,m1,k2,m2)
     2                   -corfree*sfree(k1,m1,k2,m2)
  640 continue

      return
      end

c----------------------------------------------------------------------+-----
c----------------------------------------------------------------------+-----

      subroutine tordis3(nin,x,y,z,rr)
c----------------------------------------------------------------------------c
c     determine shortest path from point x(k) to all instantons in box with  c
c     periodic boundary conditions. x(4) is reflected back in first period.  c
c----------------------------------------------------------------------------c
c     input : nin     number of instantons                                   c
c             x(k)    coordinates of point x                                 c
c             y(i,k)  coordinates of i-th instanton                          c
c                     note that y(i,5) is instanton size                     c
c     output: z(i,k)  vector of shortest path to i                           c
c             rr(i)   distance squared to instanton i                        c
c----------------------------------------------------------------------------c
      parameter (n2=256)
      common /box/ alb(4)
      dimension x(4),y(n2,5),z(n2,4), rr(n2)

      x(4) = x(4) - alb(4)*int(x(4)/alb(4))

      do 5 is = 1, nin
        rr(is) = 0.
   5  continue

       do 1 m=1,3
          do 10 is = 1, nin
             dis=x(m)-y(is,m)
             adis=abs(dis)
             asg = 1+sign(1.0,adis-alb(m)/2)
             dis = dis - alb(m)/2*asg*sign(1.0,dis)
             rr(is) = rr(is) + dis**2
             z(is,m)=dis
  10     continue
   1  continue

      m=4
      do 20 is = 1, nin
         dis=x(m)-y(is,m)
         adis=abs(dis)
         asg = 1+sign(1.0,adis-alb(m)/2)
         dis = dis - alb(m)/2*asg*sign(1.0,dis)
         z(is,m)=dis
         rr(is) = rr(is) + dis**2
  20  continue

      return
      end

c---------------------------------------------------------------------+-----
c---------------------------------------------------------------------+-----

      subroutine rmtinv(n,nin,nd,zr,e1,e2,rmu,rms)
c----------------------------------------------------------------------------c
c     calculate fermionic overlap matrix elements at finite temperature      c
c     and invert (T+im). the temperature is controlled by /box/alb(4).       c
c     if /cold/not is set to 1, beta=1000.                                   c
c----------------------------------------------------------------------------c
c     this version uses jac's parametrzation of the finite T overlap matrix  c
c     elements. for comparison, the zero T sum ansatz and stream line matrix c
c     elements are also included. revised version has modified f2 to give    c
c     smoother zero temperature limit.                                       c
c----------------------------------------------------------------------------c
c     n         number of colors                                             c
c     nin       number of instantons                                         c
c     nd        max number of instantons (also in param.!)                   c
c     zr(ni,5)  collective coordinates                                       c
c     e1(i,6)   first row of orientation matrix                              c
c     e2(i,6)   second row of orientation matrix                             c
c     rmu       u quark mass                                                 c
c     rms       s quark mass                                                 c
c----------------------------------------------------------------------------c

      parameter(n2=256,ni=n2/2,lwork=2048)
      complex work
      complex cl, clps, clpu
      dimension rbuf(500)
      complex clp(n2,n2)
      complex u1, u2, tr, ctr, uc, utr, u, utest
      complex ci, ur, uu4, usus4
      complex u1c(n2,n2), u2c(n2,n2), u3c(n2,n2), u4c(n2,n2)
      complex ofac(n2,n2)
      complex pcl1(n2,n2), pcl2(n2,n2)
      dimension dis(n2,n2,4)
      dimension zr(nd,5), u(n2,3,3),u2(n2,3,3),e1(nd,6), e2(nd,6)
      dimension tr(n2,2,2)
      dimension z1(5), z2(5), uc(n2,2,2), ur(n2,2,2)
      dimension ra(n2), rai(n2), si(n2)
      dimension du1(n2), du2(n2)
      dimension drfc(n2), rfc(n2), zr2(n2), bet(n2)
      dimension r2s(n2,n2), r2t(n2,n2)
      dimension tfac(10)
      dimension work(lwork), ipiv(n2)

c---------------------------------------------------------------------------c
c     overlaps clpu and orientation u1 are passed to pfull                  c
c---------------------------------------------------------------------------c
c     common /clspect/ cl used in spect                                     c
c---------------------------------------------------------------------------c

      common /u1/ u1(n2,3,3)
      common /clspect/ cl(n2,n2)
      common /clp/ clpu(n2,n2), clps(n2,n2)
      common /param/ a, alpha
      common /cden/ b,alc,p1,p2
      common /pi/ pi, eps
      common /sij/ sij(n2,n2)
      common /box/ alb(4)
      common /c1c2/ c1, c2
      common /cold/ not

c---------------------------------------------------------------------------c
c     check zero tempearture limit                                          c
c---------------------------------------------------------------------------c

      bett = alb(4)
      if (not .eq. 1) bett=1000.0
      tinv = bett

c---------------------------------------------------------------------------c
c     sij=0,1 for i=j, i<j; si=-/+ for instanton/antiinstanton              c
c---------------------------------------------------------------------------c

      ci = cmplx(0.0,1.0)
      nih = nin/2
      ahl = a/2

      do 2 i = 1, nin
      do 2 j = 1, nin
           sij(j,i) = abs(sign(1.0,j-i+0.5)+sign(1.0,j-i-0.5))/2.0
  2   continue

      do 5 j = 1, nin
        si(j) = sign(1.0,j-nih-0.5)
        zr2(j) = zr(j,5)*zr(j,5)
  5   continue

c---------------------------------------------------------------------------c
c     reconstruct orientation u1, u2=u1^\dagger                             c
c---------------------------------------------------------------------------c

      call su3(n,nin,nd,e1, e2,u)

      do 30 i1 = 1, n
      do 30 i2 = 1, n
         do 31 i = 1, nin
            u1(i,i2,i1) = u(i,i2,i1)
            u2(i,i1,i2) = conjg(u1(i,i2,i1))
   31    continue
   30 continue

c---------------------------------------------------------------------------c
c     calculate distances: dis(j,i,mu) = (z_I-z_J)_\mu                      c
c---------------------------------------------------------------------------c
c     note order : antiinstanton j, instanton i                             c
c---------------------------------------------------------------------------c

      iop = 0
      if(iop.eq.0 .or. iop.eq.2) then

      do 121 i = 1, nin
         do 21 j = 1, nin
            r2s(j,i) = 0.0
            r2t(j,i) = 0.0
  21    continue
  121 continue

      do 212 i = 1, nih
         do 20 ir = 1, 3
         do 22 j = nih+1, nin
            ds = zr(i,ir) - zr(j,ir)
            ads = abs(ds)
            asg = 1+sign(1.0,ads-alb(ir)/2.0)
            dis(j,i,ir) = ds - alb(ir)/2.0*asg*sign(1.0,ds)
            r2s(j,i) = r2s(j,i) + dis(j,i,ir)**2
   22    continue
   20    continue

c------------------------------------------------------------------------c
c     symmetrize 4-direction ?                                           c
c------------------------------------------------------------------------c

         do 922 j = nih+1, nin
            ds = zr(i,4) - zr(j,4)
c           dis(j,i,4) = ds
c           r2t(j,i) = dis(j,i,4)**2
            ads = abs(ds)
            asg = 1+sign(1.0,ads-alb(4)/2.0)
            dis(j,i,4) = ds - alb(4)/2.0*asg*sign(1.0,ds)
            r2t(j,i) =  dis(j,i,4)**2
  922    continue
  212 continue
      endif

c---------------------------------------------------------------------------c
c     loop over instantons I, determine orientation m.e. with antiinst. J   c
c---------------------------------------------------------------------------c

      if(iop.eq.0 .or. iop.eq.2 .or. iop.eq.3) then
      do 150 i = 1, nih

c---------------------------------------------------------------------------c
c     tr(j,a,b) = i\tau^(+)\cdot R_{IJ}                                     c
c---------------------------------------------------------------------------c

        do 90 j = nih+1, nin
          sg = 1 - (1+si(i))*(1-si(j))/2
          tr(j,1,1) = ci*cmplx(dis(j,i,3),-sg*dis(j,i,4))
          tr(j,2,2) = ci*cmplx(-dis(j,i,3),-sg*dis(j,i,4))
          tr(j,1,2) = ci*cmplx(dis(j,i,1),-dis(j,i,2))
          tr(j,2,1) = ci*cmplx(dis(j,i,1), dis(j,i,2))
   90   continue

c---------------------------------------------------------------------------c
c     relative orientation matrix UC = U_I^(+) U_J                          c
c---------------------------------------------------------------------------c

        do 35 k = 1, 2
        do 35 l = 1, 2
        do 95 j = nih+1, nin
          uc(j,k,l) = cmplx(0.0, 0.0)
   95   continue

        do 35 m = 1, n
          do 97 j = nih+1, nin
            uc(j,k,l) = uc(j,k,l) + u2(i,k,m)*u1(j,m,l)
   97     continue
   35   continue

c------------------------------------------------------------------------c
c     orientation vector u1c_JI = tr(UC_IJ \tau_1^-)/2i                  c
c------------------------------------------------------------------------c
c     note: see definition of cl, effectively define u with \tau^+       c
c------------------------------------------------------------------------c

         do 115 j = nih+1, nin
            u1c(j,i) =  (uc(j,1,2)+uc(j,2,1))/2/ci
            u2c(j,i) = -(-uc(j,1,2)+uc(j,2,1))/2
            u3c(j,i) =  (uc(j,1,1)-uc(j,2,2))/2/ci
            u4c(j,i) =  (uc(j,1,1)+uc(j,2,2))/2
            ofac(j,i)=  (u1c(j,i)*dis(j,i,1) + u2c(j,i)*dis(j,i,2)
     1                   + u3c(j,i)*dis(j,i,3))
            utest    =  (ofac(j,i)-u4c(j,i)*dis(j,i,4))/
     1                     sqrt(r2s(j,i)+r2t(j,i))
c           write(56,*) j,i,ofac(j,i),u4c(j,i),utest
  115    continue

c---------------------------------------------------------------------------c
c     end of loop over instantons I                                         c
c---------------------------------------------------------------------------c

  150 continue

      end if

c---------------------------------------------------------------------------c
c     calculate overlap matrix elements cl(i,j)                             c
c---------------------------------------------------------------------------c

      do 12 i = 1, nih
      ihal = nint((-si(i)+1)/2.0) * nih
      do 80 jlp = 1, nih
         j = jlp + ihal

c------------------------------------------------------------------------c
c     everything in units of sqrt(rho_I rho_J)                           c
c------------------------------------------------------------------------c

         betas=bett**2/zr(i,5)/zr(j,5)
         bets = betas
         pibi = pi/sqrt(betas)
         betn=sqrt(betas)
         rhon=sqrt(zr(i,5)*zr(j,5))

c------------------------------------------------------------------------c
c     coefficients:  7: f1; 8: kappa1; 9: f2; 10: kappa2                 c
c------------------------------------------------------------------------c

         tfac(7)=bets/(.191*bets+1.)**2
c        tfac(8)=1.07/(bets/pi**2+0.533)
c        tfac(9)=1.+.069*bets/(0.011*bets+1.)**2
         tfac(10)=1./(bets/pi**2+.69)

c------------------------------------------------------------------------c
c     minimial fix for improved small T behavior                         c
c------------------------------------------------------------------------c

         tfac(8)=1.00/(bets/pi**2+0.533)
         tfac(9)=1.

c------------------------------------------------------------------------c
c        variables for fermionic overlap, units of rho_I*rho_J           c
c------------------------------------------------------------------------c

         rat = zr(i,5)/zr(j,5)
         fact= 4/(rat+1/rat)**2
         rn  = (r2s(j,i)+r2t(j,i))/zr(i,5)/zr(j,5)
         r2i = sij(j,i)/(r2s(j,i)+eps)

         r1t = dis(j,i,4)/rhon
         r1s = sqrt(r2s(j,i)/zr(i,5)/zr(j,5))
         sn = sin(pibi*r1t)
         cs = cos(pibi*r1t)
         sh = sinh(pibi*r1s)
         ch = cosh(pibi*r1s)

c------------------------------------------------------------------------c
c     jac's parametrization                                              c
c------------------------------------------------------------------------c

         fac1 = 1/((exp(-pibi*r1s/2) + pibi*r1s/2)/pibi**2
     1                +2*(1-0.69*exp(-1.75*r1s/betn))/pi)
         fac2 = 1+ 0.76*tfac(7)/(1+0.82*r1s*r1s)**2
         fac2 = fac2*(1+(pibi*r1t)**2*0.178/(1+0.123*r1s*r1s))
         fac3 = 1/((exp(-2.06*r1s/betn) + pibi*r1s/2)/pibi**2
     1                +2*(1+0.42*exp(-.34*r1s/betn))/pi)
         fac4 = tfac(9)

c------------------------------------------------------------------------c
c     my simplified version                                              c
c------------------------------------------------------------------------c

c        aa = 0.5*pi**2*(2.0+2.0/3.0)
c        xx = exp(-aa/bets)
c        yy = aa/(bets+aa)
c        fac1 = pibi**2*xx + pi/2.0*yy*(1.0-xx)
c        fac2 = 1.0
c        fac3 = fac1
c        fac4 = 1.0

c------------------------------------------------------------------------c
c        T_IJ = u_4*f1+...; pcl1= f1, pcl2= f2                           c
c------------------------------------------------------------------------c

         pcl1(j,i) = pibi*sn*ch/(ch-cs+tfac(8))**2 *fac1 *fac2
         pcl2(j,i) = pibi*cs*sh/(ch-cs+tfac(10))**2 *fac3 *fac4

         pcl1(j,i) = 1/rhon * pcl1(j,i)
         pcl2(j,i) = 1/rhon * pcl2(j,i) * sqrt(r2i)

c------------------------------------------------------------------------c
c     check: zero temperature limit                                      c
c------------------------------------------------------------------------c

c        pcl1(j,i)= 4./(r1s**2+r1t**2+2.)**2/rhon**2 * dis(j,i,4)
c        pcl2(j,i)= 4./(r1s**2+r1t**2+2.)**2/rhon**2

c------------------------------------------------------------------------c
c     zero temperature streamline                                        c
c------------------------------------------------------------------------c

         x2  = r2s(j,i)+r2t(j,i)
         acf = (x2+zr2(j)+zr2(i))/rhon**2
         disc= sqrt(acf**2-4.0)
         rlam= (acf+disc)/2.0
         rl2 = rlam**2
         fstr= c1*rlam*sqrt(rlam)/
     2           (1.0+1.25*(rl2-1.0)+c2*(rl2-1.0)**2)**0.75

c        pcl1(j,i) = fstr/sqrt(x2)/rhon * dis(j,i,4)
c        pcl2(j,i) = fstr/sqrt(x2)/rhon

c------------------------------------------------------------------------c
c     end of loop over instantons j                                      c
c------------------------------------------------------------------------c

 80   continue
 12   continue

c------------------------------------------------------------------------c
c     construct overlap matrix cl                                        c
c------------------------------------------------------------------------c
c     note: cl has indices in strange order                              c
c------------------------------------------------------------------------c

      do 275 i = 1, nin
      do 275 j = 1, nin
        cl(j,i) = cmplx(0.0,0.0)
  275 continue

      do 280 i = 1, nih
      do 280 j = nih+1, nin

c---------------------------------------------------------------------------c
c     T_IA as defined in DP Streamline ansatz                               c
c---------------------------------------------------------------------------c

         cl(j,i) = (-u4c(j,i)*pcl1(j,i)+ofac(j,i)*pcl2(j,i))*ci
         cl(i,j) = conjg(cl(j,i))
c        write(55,*) i,j,cl(i,j)

 280  continue

c---------------------------------------------------------------------------c
c     add current mass i(T+im)=(iT-m), copy result in clp                   c
c---------------------------------------------------------------------------c
c     note : clp_IJ has indices in correct order                            c
c---------------------------------------------------------------------------c

      do 290 i = 1, nin
         cl(i,i) = cmplx(0.0,-rmu)
         do 295 j = 1, nin
            clp(i,j) = cl(j,i)
            clpu(i,j)= cl(j,i)
  295    continue
  290 continue

c---------------------------------------------------------------------------c
c     invert light quark (iT-m)^(-1) --> clpu, imsl                         c
c---------------------------------------------------------------------------c

c     call lincg(nin,clp,n2,clpu,n2)

c---------------------------------------------------------------------------c
c     invert light quark (iT-m), lapack (clpu -> clpu)                      c
c---------------------------------------------------------------------------c

      call cgetrf(nin,nin,clpu,n2,ipiv,info)
      if(info .ne. 0) stop
      call cgetri(nin,clpu,n2,ipiv,work,lwork,info)
      if(info .ne. 0) stop

c---------------------------------------------------------------------------c
c     same for ms (iT-ms)^(-1) --> clps                                     c
c---------------------------------------------------------------------------c

      do 490 i = 1, nin
         cl(i,i) = cmplx(0.0,-rms)
         do 495 j = 1, nin
            clp(j,i) = cl(i,j)
            clps(j,i)= cl(i,j)
  495    continue
  490 continue

c---------------------------------------------------------------------------c
c     imsl                                                                  c
c---------------------------------------------------------------------------c

c     call lincg(nin,clp,n2,clps,n2)

c---------------------------------------------------------------------------c
c     lapack                                                                c
c---------------------------------------------------------------------------c

      call cgetrf(nin,nin,clps,n2,ipiv,info)
      if(info .ne. 0) stop
      call cgetri(nin,clps,n2,ipiv,work,lwork,info)
      if(info .ne. 0) stop

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine tspect(n,nin,nd,zr,e1,e2,rmdet,wr,sdet,t2av,uin)
c----------------------------------------------------------------------------c
c     determine spectrum of (T+im)^(-1), calculate det and trace             c
c----------------------------------------------------------------------------c
c     n         number of colors                                             c
c     nin       number of instantons                                         c
c     nd        max number of instantons (also in param.!)                   c
c     zr(ni,5)  collective coordinates                                       c
c     e1(i,6)   first row of orientation matrix                              c
c     e2(i,6)   second row of orientation matrix                             c
c     rmdet     exp(sdet/(nc*nf))                                            c
c     wr(i)     eigenvalues of iT                                            c
c     sdet      log det(T^+T+m^2)                                            c
c     t2av      average overlap |T_IJ|^2                                     c
c     uin       not used                                                     c
c----------------------------------------------------------------------------c

      parameter(ni=256, ni2=ni/2, n2=ni, lwork=768, lrwork=1024)
      complex work
      complex clp, clt, cjc, cl2, clap, cl
      complex ctr, utr, uc1r, uc2r, uc3r, uc4r
      dimension clp(ni2,ni2), zr(nd,5), rlp(ni2,ni2), rcos(ni2,ni2)
      dimension cl2(ni2,ni2), cjc(ni2,ni2), rnii(ni2,ni2),rnia(ni2,ni2)
      dimension rcii(ni2, ni2)
      dimension ar(ni2,ni2), ai(ni2,ni2), wr(ni2), wi(ni2)
      dimension e1(nd,6), e2(nd,6)
      dimension work(lwork), rwork(lrwork)

c---------------------------------------------------------------------------c
c     overlap matrix cl calculated in rminv                                 c
c---------------------------------------------------------------------------c
c     masses in common /param/                                              c
c---------------------------------------------------------------------------c

      common /d/ uint(n2,n2), r2(n2,n2), dis(n2,n2,4), pcl(n2,n2)
      common /ucr/ uc1r(n2,n2), uc2r(n2,n2), uc3r(n2,n2), uc4r(n2,n2)
      common /param/ a, alpha,rh0,sg,dz,drh, nc, nf, rmu, rms
      common /pi/ pi, eps
      common /sij/ sij(ni,ni)
      common /clspect/ cl(n2,n2)

c---------------------------------------------------------------------------c
c     call rminv (in case it has not been done yet)                         c
c---------------------------------------------------------------------------c

      nih = nin/2
      nd2 = nd/2
      nar = nin

      call rmtinv(n,nin,nd,zr,e1,e2,rmu,rms)

c---------------------------------------------------------------------------c
c     clp (in this subr. only!) is upper right block of cl.                 c
c---------------------------------------------------------------------------c
c     note : remember that indices in cl are interchanged                   c
c---------------------------------------------------------------------------c

      do 5 i = 1, nih
      do 5 j = 1, nih
        clp(i,j) = cl(j+nih,i)
   5  continue

c---------------------------------------------------------------------------c
c     calculate average overlap matrix element                              c
c---------------------------------------------------------------------------c

      t2av = 0.0
      do 7 i=1, nih
      do 7 j=1, nih
         t2av = t2av + cabs(clp(i,j))**2
   7  continue
      t2av = t2av/nih

c---------------------------------------------------------------------------c
c     calculate clt = T*T^\dagger                                           c
c---------------------------------------------------------------------------c

      do 10 i = 1, nih
      do 10 j = 1, nih
        cjc(j,i) = conjg(clp(j,i))
   10 continue
      do 12 i = 1, nih
      do 12 j = 1, nih
        clt = cmplx(0.0,0.0)
        do 214 k = 1, nih
          clt = clt + clp(k,i)*cjc(k,j)
  214   continue
        cl2(j,i) = clt
   12 continue
      do 16 i = 1, nih
      do 16 j = 1, nih
        ar(i,j) = real(cl2(i,j))
        ai(i,j) = aimag(cl2(i,j))
   16 continue

c---------------------------------------------------------------------------c
c     diagonalize T*T^\dagger, imsl                                         c
c---------------------------------------------------------------------------c

c     call evlhf(nih,cl2,ni2,wr)

c---------------------------------------------------------------------------c
c     diagonalize, lapack                                                   c
c---------------------------------------------------------------------------c

      call cheev('n','l',nih,cl2,ni2,wr,work,lwork,rwork,info)
      if(info .ne. 0) stop

c---------------------------------------------------------------------------c
c     determinants                                                          c
c---------------------------------------------------------------------------c

      rdet = 0.0
      sdet = 0.0
      do 20 i = 1, nih
         rhp = zr(i,5)*zr(i+nih,5)
c        rhp = 1
         rms2 = rms*rms*rhp
         rmu2 = rmu*rmu*rhp
         eig = wr(i)*rhp
         sdet = sdet+(nf-1)*alog(eig+rmu2) + alog(rms2+eig)
         rdet = rdet+(nf-1)*alog(eig) + alog(eig)
   20 continue
      if (nf.eq.0) then
        rmdet = 0.0
      else
        rmdet = exp(sdet/(nf*nin))
      end if
      bofac = 0.0
      if (nf .eq. 0) then
        sdet = bofac
      else
        sdet = sdet + bofac
      end if

c---------------------------------------------------------------------------c
c     wr(i) : positive eigenvalues of T                                     c
c---------------------------------------------------------------------------c
c     calculate <qq> (in fm^(-3) for unit rho = 1 fm^(-4))                  c
c---------------------------------------------------------------------------c

      som = 0.0
      do 30 i = 1, nih
        if(wr(i) .le. 0) write(6,*) i,wr(i)
        wr(i) = sqrt(abs(wr(i)))
        som = som + 2*rmu/(rmu*rmu+wr(i)*wr(i))
   30 continue

      print 101, som/nin
  101 format(1x, 'sum_i=1,nin/2  2*rmu/(rmu*rmu+wr(i)*wr(i)) ', f10.4)

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine ptfree(xx,xmass,sfree)
c--------------------------------------------------------------------------c
c     calculate free finite temperature quark propagator from point        c
c     xx(.,1) to xx(.,2).                                                  c
c--------------------------------------------------------------------------c
c     temperature given by inverse length of box in 4-direction (from      c
c     common /box/), propagator satifies no b.c. in other directions.      c
c--------------------------------------------------------------------------c
c     input : xx(k,index)    endpoints; k=1,..,4 coordinates,              c
c                            index=1,2 labels beginning, endpoint          c
c     output: sfree(a,i,b,j) free massless propagator; a,b=1,..3           c
c                            color, i,j=1,..,4 dirac indices.              c
c--------------------------------------------------------------------------c

      real xx(4,2),dx(4)
      real alb(4)
      real gr(4,4),gl(4,4),c(4,4),u(4,4)
      complex sfree(3,4,3,4),gamma(5,4,4),coef,coef2
      integer a,b

      common /gam/ gamma,gr,gl,c
      common /box/ alb
      common /cold/not

c--------------------------------------------------------------------------c
c     inverse temperature from box dimension                               c
c--------------------------------------------------------------------------c

      beta = alb(4)
      if(not .eq. 1) beta=1000.0
      t = 1.0/beta

      do 10 a=1,3
      do 10 b=1,3
         do 10 i=1,4
         do 10 j=1,4
            sfree(a,i,b,j) = (0.0,0.0)
  10  continue

c------------------------------------------------------------------------c
c     distance between end points, dx defined as delta in prnz           c
c------------------------------------------------------------------------c

      tau  = xx(4,1)-xx(4,2)
      dx(4)= tau
      r2   = 0.0
      do 20 k=1,3
         dx(k)= xx(k,1)-xx(k,2)
         r2   = r2 + dx(k)**2
  20  continue
      r  = sqrt(r2)
      pi = 3.1415926
      coef = (0.0,1.0)/(4.0*pi**2)

c------------------------------------------------------------------------c
c     sfree is color diagonal, two dirac structures                      c
c------------------------------------------------------------------------c

      st =  2*Pi**2*T**2*(2 + Cos(2*Pi*T*tau)
     1  + Cosh(2*Pi*r*T))*Sin(Pi*T*tau)*
     2    Sinh(Pi*r*T)/(r*(Cos(2*Pi*T*tau) - Cosh(2*Pi*r*T))**2)

      sr = -Pi*T*Cos(Pi*T*tau)*(3*Pi*r*T*Cosh(Pi*r*T) -
     1 2*Pi*r*T*Cos(2*Pi*T*tau)*Cosh(Pi*r*T) - Pi*r*T*Cosh(3*Pi*r*T) +
     2 Sinh(Pi*r*T) + 2*Cos(2*Pi*T*tau)*Sinh(Pi*r*T) - Sinh(3*Pi*r*T))/
     3 (r**2*(Cos(2*Pi*T*tau) - Cosh(2*Pi*r*T))**2)

c------------------------------------------------------------------------c
c     T->0 limit                                                         c
c------------------------------------------------------------------------c

c     st = 2.0*tau/(r**2+tau**2)**2
c     sr = 2.0* r /(r**2+tau**2)**2

c------------------------------------------------------------------------c
c     massless propagator                                                c
c------------------------------------------------------------------------c

      do 30 i=1,4
      do 30 j=1,4
         do 30 a=1,3
            sfree(a,i,a,j) =  coef*( gamma(4,i,j)*st
     2         + (gamma(1,i,j)*dx(1)+gamma(2,i,j)*dx(2)+
     3            gamma(3,i,j)*dx(3))*sr/r )
  30  continue

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine dtfree(xx,xmass,dfree)
c--------------------------------------------------------------------------c
c     calculate i*m*D(x) where D(x) is the free finite temperature scalar  c
c     propagator satisfying antiperiodic (fermionic) boundary conditions.  c
c--------------------------------------------------------------------------c
c     temperature given by inverse length of box in 4-direction (from      c
c     common /box/), propagator satifies no b.c. in other directions.      c
c--------------------------------------------------------------------------c
c     input : xx(k,index)    endpoints; k=1,..,4 coordinates,              c
c                            index=1,2 labels beginning, endpoint          c
c             xmass          quark current mass                            c
c     output: dfree(a,i,b,j) i*m*D_sc(x); a,b=1,..3 color                  c
c                            i,j=1,..,4 dirac indices.                     c
c--------------------------------------------------------------------------c
c     note that S(x) = sfree + dfree has the full temperature dependence   c
c     but only the leading (short distance) mass term.                     c
c--------------------------------------------------------------------------c

      real xx(4,2),dx(4)
      real alb(4)
      complex dfree(3,4,3,4),coef,coef2
      integer a,b

      common /box/ alb
      common /cold/not

c--------------------------------------------------------------------------c
c     inverse temperature from box dimension                               c
c--------------------------------------------------------------------------c

      beta = alb(4)
      if(not .eq. 1) beta=1000.0
      t = 1.0/beta

      do 10 a=1,3
      do 10 b=1,3
         do 10 i=1,4
         do 10 j=1,4
            dfree(a,i,b,j) = (0.0,0.0)
  10  continue

c------------------------------------------------------------------------c
c     distance between end points                                        c
c------------------------------------------------------------------------c

      tau = xx(4,1)-xx(4,2)
      r2 = 0.0
      do 20 k=1,3
         dx(k)= xx(k,2)-xx(k,1)
         r2   = r2 + dx(k)**2
  20  continue
      r  = sqrt(r2)
      pi = 3.1415926
      coef = (0.0,1.0)/(4.0*pi**2)*xmass

c------------------------------------------------------------------------c
c     scalar propagator, last factor makes is antiperiodic               c
c------------------------------------------------------------------------c

      dt = Pi*T/r*Sinh(2*Pi*r*T)/(Cosh(2*Pi*r*T)-Cos(2*Pi*T*tau))
     1    * Cos(Pi*T*tau)/Cosh(Pi*r*T)

c------------------------------------------------------------------------c
c     T->0 limit                                                         c
c------------------------------------------------------------------------c

c     dt = 1.0/(r**2+tau**2)

c------------------------------------------------------------------------c
c     dfree is color, dirac diagonal                                     c
c------------------------------------------------------------------------c

      do 30 i=1,4
         do 30 a=1,3
            dfree(a,i,a,i) = coef * dt
  30  continue

      return
      end

c----------------------------------------------------------------------+------
c----------------------------------------------------------------------+------

      subroutine stfree(xx,xmass,sfree)
c--------------------------------------------------------------------------c
c     calculate free finite temperature massive quark propagator from      c
c     xx(.,1) to xx(.,2). calculation uses exact zero T propagator, finite c
c     temperature is introduced by brute force, doing a truncated mode sum.c
c--------------------------------------------------------------------------c
c     temperature given by inverse length of box in 4-direction (from      c
c     common /box/), propagator satifies no b.c. in other directions.      c
c--------------------------------------------------------------------------c
c     input : xx(k,index)    endpoints; k=1,..,4 coordinates,              c
c                            index=1,2 labels beginning, endpoint          c
c     output: sfree(a,i,b,j) free massless propagator; a,b=1,..3           c
c                            color, i,j=1,..,4 dirac indices.              c
c--------------------------------------------------------------------------c

      real xx(4,2),dx(4)
      real alb(4)
      real gr(4,4),gl(4,4),c(4,4),u(4,4)
      real k(0:2)
      complex sfree(3,4,3,4),gamma(5,4,4),coef,coef2,coef3
      integer a,b

      common /gam/ gamma,gr,gl,c
      common /box/ alb
      common /cold/not

c--------------------------------------------------------------------------c
c     inverse temperature from box dimension                               c
c--------------------------------------------------------------------------c

      beta = alb(4)
      if( not .eq. 1) beta=1000.0
      t = 1.0/beta

      do 10 a=1,3
      do 10 b=1,3
         do 10 i=1,4
         do 10 j=1,4
            sfree(a,i,b,j) = (0.0,0.0)
  10  continue

      do 15 i=1,4
      do 15 j=1,4
         u(i,j) = 0.0
  15  continue
      do 16 i=1,4
         u(i,i) = 1.0
  16  continue

c------------------------------------------------------------------------c
c     distance between end points, dx defined as delta in prnz           c
c------------------------------------------------------------------------c

      r2 = 0.0
      do 20 j=1,3
         dx(j)= xx(j,1)-xx(j,2)
         r2   = r2 + dx(j)**2
  20  continue
      tau  = xx(4,1)-xx(4,2)
      dx(4)= tau
      x2 = r2+tau**2
      x  = sqrt(x2)
      r  = sqrt(r2)
      pi = 3.1415926

      arg  = x*xmass
      if (arg .ge. 86.6) then
         k(0)=0.0
         k(1)=0.0
         k(2)=0.0
      else
c------------------------------------------------------------------------c
c     imsl                                                               c
c------------------------------------------------------------------------c
c        call bsks(0.0,arg,3,k)
c------------------------------------------------------------------------c
c     cern                                                               c
c------------------------------------------------------------------------c
         k(0) = besk0(arg)
         k(1) = besk1(arg)
         k(2) = k(0)+2.0/arg*k(1)
      endif
      coef = (0.0,1.0)/(4.0*pi**2)
      coef2= coef*xmass**2/x2*k(2)
      coef3= coef*xmass**2/x*k(1)

c------------------------------------------------------------------------c
c     zero temperature                                                   c
c------------------------------------------------------------------------c

      do 30 i=1,4
      do 30 j=1,4
         do 30 a=1,3
            sfree(a,i,a,j) =
     1       coef2 * ( gamma(4,i,j)*dx(4)
     2        + gamma(1,i,j)*dx(1) + gamma(2,i,j)*dx(2)
     3        + gamma(3,i,j)*dx(3) ) + coef3*u(i,j)
  30  continue

c------------------------------------------------------------------------c
c     mode sum, nmodes= number of modes, nm=-1,+1,-2,.. index of mode    c
c------------------------------------------------------------------------c

      nmodes = 6

      do 100 i1=1,nmodes
         i2 = int(i1/2.0+0.75)
         nm = i2*(-1)**i1

         tnew = xx(4,2)+nm*beta
         tau  = xx(4,1)-tnew
         dx(4)= tau
         x2 = r2+tau**2
         x  = sqrt(x2)

c------------------------------------------------------------------------c
c     (-1)**nm takes care of antiperiodic boundary conditions            c
c------------------------------------------------------------------------c

         arg  = x*xmass
         if (arg .ge. 86.6) then
            k(0)=0.0
            k(1)=0.0
            k(2)=0.0
         else
c------------------------------------------------------------------------c
c     imsl                                                               c
c------------------------------------------------------------------------c
c           call bsks(0.0,arg,3,k)
c------------------------------------------------------------------------c
c     cern                                                               c
c------------------------------------------------------------------------c
            k(0) = besk0(arg)
            k(1) = besk1(arg)
            k(2) = k(0)+2.0/arg*k(1)
         endif
         coef2= coef*xmass**2/x2*k(2) * (-1)**abs(nm)
         coef3= coef*xmass**2/x *k(1) * (-1)**abs(nm)

c------------------------------------------------------------------------c
c     add mode to zero T piece                                           c
c------------------------------------------------------------------------c

         do 130 i=1,4
         do 130 j=1,4
            do 130 a=1,3
               sfree(a,i,a,j) = sfree(a,i,a,j) +
     1            coef2 * ( gamma(4,i,j)*dx(4)
     2          + gamma(1,i,j)*dx(1) + gamma(2,i,j)*dx(2)
     3          + gamma(3,i,j)*dx(3) ) + coef3*u(i,j)
 130     continue

 100  continue

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

