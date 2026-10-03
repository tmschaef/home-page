      program corglue
c--------------------------------------------------------------------------c
c     glueball correlators and wave functions. This program uses the ran-  c
c     dom instanton model. Other configurations can be read in from infile c
c--------------------------------------------------------------------------c
c     version           :      1.5                                         c
c     creation date     :    06-06-94                                      c
c     last modification :    10-27-94                                      c
c--------------------------------------------------------------------------c
c     version 1.2 : modified setup can read configurations from infile.dat.c
c     real part of gluonic operators is passed to main.                    c
c--------------------------------------------------------------------------c
c     version 1.3 : introduced running coupling constant g(x). Fixed some  c
c     of the problem with infinite check loops. Corrected error estimates. c
c--------------------------------------------------------------------------c
c     version 1.4 : included double precision version of fmunu, potsu3,    c
c     sumsu3 and tordis.                                                   c
c--------------------------------------------------------------------------c
c     version 1.5 : corrected subtraction procedure and error estimate for c
c     scalar correlator.                                                   c
c--------------------------------------------------------------------------c
c     input : inglue.dat                                                   c
c             infile.dat                                                   c
c     output: outglue.dat                                                  c
c--------------------------------------------------------------------------c
      parameter(n=3,ni=16,nbin=150,ncf=100,ndl = 20)

      common /seed/  iseed
      common /param/ a, alpha,rh0,sg,dz,drh, nc, nf, rmu, rms
      common /nconf/ nconfig
      common /box/   alb(4)
      common /counter/ icount

      dimension zr(ni,5),e1(ni,6),e2(ni,6)
      complex pexp(3,3)
      dimension point(4),point1(4),point2(4)
      dimension ndis(ndl),scmy(ndl),psmy(ndl),t1my(ndl),t2my(ndl)
      dimension sc2my(ndl),ps2my(ndl),t12my(ndl),t22my(ndl)
      dimension scmya(ndl),psmya(ndl),t1mya(ndl),t2mya(ndl)
      dimension scmye(ndl),psmye(ndl),t1mye(ndl),t2mye(ndl)

      nran = 0

c--------------------------------------------------------------------------c
c     input parameters                                                     c
c--------------------------------------------------------------------------c

      open (unit=1, file='inglue.dat', status='old')
      read (1,*) nc, nf, nin, rh0,iread, nconfig,index
      read (1,*) nitc, delt, splitin, ndelt
      read (1,*) (alb(k), k = 1, 4)
      close (unit=1)

      open (unit=2, file='outglue.dat', status='new')
      write (2,221) nc,nf,nin
      write (2,222) rh0,iread
      write (2,223) nconfig,nitc
      write (2,224) delt,splitin,ndelt
      write (2,225) (alb(k), k = 1, 4)

 221  format('nc = ',i3,' nf = ',i3,' nin = ',i3)
 222  format('rho = ',f8.3,' iread = ',i2)
 223  format('nconf = ',i6,' nitc = ',i6)
 224  format('del = ',f8.3,' split = ',f8.3,' nd = ',i3)
 225  format(' box : ',4(f8.3))

      iseed=9234
      ni2=ni/2

      keyt = 0
      volume = alb(1)*alb(2)*alb(3)*alb(4)
      ndd = ni
      pi  = 3.1415926
      nih = nin/2

      if (index .eq. 1) then
         taumax = splitin*ndelt
      else if (index .eq. 2) then
         taumax = delt
         splitmax = splitin*ndelt
      endif

      icon = 0
      ncon = 0
      nch  = 0
      ncht = 0
      nrej = 0
      nmax = 10
      eps  = 1.e-5

c--------------------------------------------------------------------------c
c     clear summation arrays                                               c
c--------------------------------------------------------------------------c

      call myzero(ndl,scmy)
      call myzero(ndl,psmy)
      call myzero(ndl,t1my)
      call myzero(ndl,t2my)
      call myzero(ndl,sc2my)
      call myzero(ndl,ps2my)
      call myzero(ndl,t12my)
      call myzero(ndl,t22my)
      call myzero(ndl,scmya)
      call myzero(ndl,psmya)
      call myzero(ndl,t1mya)
      call myzero(ndl,t2mya)
      call myzero(ndl,scmye)
      call myzero(ndl,psmye)
      call myzero(ndl,t1mye)
      call myzero(ndl,t2mye)

c--------------------------------------------------------------------------c
c     local quantities                                                     c
c--------------------------------------------------------------------------c

      nloc=0
      scl  = 0.0
      psl  = 0.0
      t1l  = 0.0
      t2l  = 0.0
      sc2l = 0.0
      ps2l = 0.0
      t12l = 0.0
      t22l = 0.0

      scl2  = 0.0
      psl2  = 0.0
      t1l2  = 0.0
      t2l2  = 0.0
      sc2l2 = 0.0
      ps2l2 = 0.0
      t12l2 = 0.0
      t22l2 = 0.0

      do 544 k = 1, ndelt
         ndis(k)=0
 544  continue

      icount = 0

c--------------------------------------------------------------------------c
c     new configuration                                                    c
c--------------------------------------------------------------------------c

      do 180 ic=1,nconfig

      call setup(nc,nin,ndd,zr,e1,e2,iread,icon)

      itel = 0

      do 61 m1=1,3
      do 61 m2=1,3
         pexp(m1,m2)=(0.,0.)
         if(m1.eq.m2) pexp(m1,m2)=(1.,0.)
61    continue

c--------------------------------------------------------------------------c
c     average over initial point                                           c
c--------------------------------------------------------------------------c

      do 100 i = 1, nitc
         itel = itel + 1
         nch  = 0

777      continue

c--------------------------------------------------------------------------c
c     index = 1,2 corresponds to correlator/wavefunction                   c
c--------------------------------------------------------------------------c
c     correlator measured in 4-direction, wavefct split in 1-direction     c
c--------------------------------------------------------------------------c

         if(index. eq.1) then
           point(1)=alb(1)*rang()
         else
           point(1)=alb(1)*0.5
c          point(1)=(alb(1)-splitmax)*rang()+splitmax/2.0
         endif
         point(2)= alb(2)*rang()
         point(3)= alb(3)*rang()
c        point(4)= 0.35
         point(4)=(alb(4)-taumax)*rang()

         do 15 in=1,nin
         do 15 ip=1,4
            dis = abs(zr(in,ip)-point(ip))
            a1  = alb(ip)/2.0-eps
            a2  = alb(ip)/2.0+eps
            if( dis .gt. a1 .and. dis .lt. a2)then
                nrej = nrej+1
                goto 777
            endif
 15      continue

c--------------------------------------------------------------------------c
c     local gluonic observables                                            c
c--------------------------------------------------------------------------c

         call glue(nin,ndd,zr,e1,e2,point,sc,ps,t1,t2)

c--------------------------------------------------------------------------c
c     outrageous stuff is eliminated                                       c
c--------------------------------------------------------------------------c

         check=sc*rh0**4
         if(check.gt.500.0) then
           nch = nch + 1
           ncht= ncht+ 1
           if(nch .gt. nmax) goto 180
           go to 777
         endif


         nloc=nloc+1

         sc2 = sc**2
         ps2 = ps**2
         t12 = t1**2
         t22 = t2**2

         call myaddto(sc,scl,scl2)
         call myaddto(ps,psl,psl2)
         call myaddto(t1,t1l,t1l2)
         call myaddto(t2,t2l,t2l2)

         call myaddto(sc2,sc2l,sc2l2)
         call myaddto(ps2,ps2l,ps2l2)
         call myaddto(t12,t12l,t12l2)
         call myaddto(t22,t22l,t22l2)

c--------------------------------------------------------------------------c
c     loop over endpoints                                                  c
c--------------------------------------------------------------------------c

         do 200 k = 1, ndelt

            do 710 kk=1,4
               point1(kk)=point(kk)
               point2(kk)=point(kk)
  710       continue

            if(index. eq.1)  then

c--------------------------------------------------------------------------c
c     endpoint for correlator                                              c
c--------------------------------------------------------------------------c

            delt=(k-1)*splitin
            point1(4) = point(4) + delt
            point2(4) = point(4) + delt

            call glue(nin,ndd,zr,e1,e2,point1,sc2,ps2, t12,t22)

            else

c-------------------------------------------------------------------------c
c     two endpoints for wavefunction                                      c
c-------------------------------------------------------------------------c

            split = splitin*(k-1)
            point1(1)= point(1)+0.5*split
            point2(1)= point(1)-0.5*split
            point1(4) = point(4) + delt
            point2(4) = point(4) + delt

            call sglue(nin,ndd,zr,e1,e2,point1,point2,sc2,ps2, t12,t22)

            endif

c-------------------------------------------------------------------------c
c     again, eliminate results that are too big                           c
c-------------------------------------------------------------------------c

            check=sc2*rh0**4
            if(check.gt.500.0) then
              nch = nch + 1
              ncht= ncht+ 1
              if(nch .gt. nmax) goto 180
              go to 777
            endif

            prsc=sc*sc2
            prps=ps*ps2
            prt1=t1*t12
            prt2=t2*t22

            ndis(k)=ndis(k)+1
            call myaddto(prsc,scmy(k),sc2my(k))
            call myaddto(prps,psmy(k),ps2my(k))
            call myaddto(prt1,t1my(k),t12my(k))
            call myaddto(prt2,t2my(k),t22my(k))

c-------------------------------------------------------------------------c
c     end of loop over endpoint                                           c
c-------------------------------------------------------------------------c

  200    continue

c-------------------------------------------------------------------------c
c     end of loop over initial point                                      c
c-------------------------------------------------------------------------c

  100 continue

c-------------------------------------------------------------------------c
c     new configuration                                                   c
c-------------------------------------------------------------------------c

      ncon = ncon + 1

 180  continue

c-------------------------------------------------------------------------c
c     calculate averages                                                  c
c-------------------------------------------------------------------------c

       call mydisp(nloc,scl,scl2,sca,sce)
       call mydisp(nloc,psl,psl2,psa,pse)
       call mydisp(nloc,t1l,t1l2,t1a,t1e)
       call mydisp(nloc,t2l,t2l2,t2a,t2e)

       call mydisp(nloc,sc2l,sc2l2,sc2a,sc2e)
       call mydisp(nloc,ps2l,ps2l2,ps2a,ps2e)
       call mydisp(nloc,t12l,t12l2,t12a,t12e)
       call mydisp(nloc,t22l,t22l2,t22a,t22e)

      do 545 k=1,ndelt
         call mydisp(ndis(k),scmy(k),sc2my(k),scmya(k),scmye(k))
         call mydisp(ndis(k),psmy(k),ps2my(k),psmya(k),psmye(k))
         call mydisp(ndis(k),t1my(k),t12my(k),t1mya(k),t1mye(k))
         call mydisp(ndis(k),t2my(k),t22my(k),t2mya(k),t2mye(k))
 545  continue

c-----------------------------------------------------------------------c
c     output, local average                                             c
c-----------------------------------------------------------------------c

      write(6,*) 'number of rejected points',ncht
      write(6,*) 'rejceted in advance      ',nrej
      write(2,*)
      write(2,*) 'scalar       :',sca,'+/-',sce
      write(2,*) 'pseudoscalar :',psa,'+/-',pse
      write(2,*) 'tensor 1     :',t1a,'+/-',t1e
      write(2,*) 'tensor 2     :',t2a,'+/-',t2e

      write(2,*)
      write(2,*) 'sc^2         :',sc2a,'+/-',sc2e
      write(2,*) 'ps^2         :',ps2a,'+/-',ps2e
      write(2,*) 't1^2         :',t12a,'+/-',t12e
      write(2,*) 't2^2         :',t22a,'+/-',t22e

c-----------------------------------------------------------------------c
c     subtraction constant for scalar                                   c
c-----------------------------------------------------------------------c
c     various possibilities: global <G^2>, correlator on last point,etc.c
c-----------------------------------------------------------------------c

      scconst = sca**2/sc2a
      scdisc  = sca**2
      scdisce = 2.0*sca*sce
      scdisc2 = scmya(ndelt)
      scdisc2e= scmye(ndelt)
      smin    = scmya(1)
      do ix=2,ndelt
         smin = min(smin,scmya(ix))
      enddo
      scdisc3 = smin
      scdisc3e= scdisce
      write(2,*)
      write(2,521) 1./scconst
 521  format(1x,' G^4/(G^2)^2 = ',f10.3)

c-----------------------------------------------------------------------c
c     scalar correlator                                                 c
c-----------------------------------------------------------------------c

      write(2,*)
      if (index.eq.1) then
         write(2,550)
      else
         write(2,1550) scmya(1)
      endif
 550  format(1x, ' scalar ')
1550  format(1x, ' scalar        Pi(x) = ',g12.5)
      delstep=splitin

c-----------------------------------------------------------------------c
c     note that only correlator is subtracted                           c
c-----------------------------------------------------------------------c

      do 551 k=1,ndelt
         delt=(k-1)*delstep
         if (index.eq.1) then
            scal = (scmya(k)-scdisc)/(scmya(1)-scdisc)
            write(2,*) delt,scal,scmye(k)/scmya(1)
         else
            write(2,*) delt,scmya(k)/scmya(1),scmye(k)/abs(scmya(1))
         endif
 551  continue

c-----------------------------------------------------------------------c
c     peudoscalar correlator                                            c
c-----------------------------------------------------------------------c

      write(2,*)
      if (index.eq.1) then
         write(2,552)
      else
         write(2,1552) psmya(1)
      endif
 552  format(1x, ' pseudoscalar ')
1552  format(1x, ' pseudoscalar  Pi(x) = ',g12.5)

      do 553 k=1,ndelt
         delt=(k-1)*delstep
         write(2,*)delt,psmya(k)/psmya(1),psmye(k)/abs(psmya(1))
553   continue

c-----------------------------------------------------------------------c
c     tensor correlator                                                 c
c-----------------------------------------------------------------------c

      write(2,*)
      if (index.eq.1) then
         write(2,554)
      else
         write(2,1554) t1mya(1)
      endif
 554  format(1x, ' tensor 1 ')
1554  format(1x, ' tensor 1      Pi(x) = ',g12.5)

      do 555 k=1,ndelt
         delt=(k-1)*delstep
         write(2,*)delt,t1mya(k)/t1mya(1),t1mye(k)/abs(t1mya(1))
555   continue

      write(2,*)
      write(2,556)
556   format(1x, ' tensor 2 ')
      do 557 k=1,ndelt
         delt=(k-1)*delstep
         write(2,*) delt,t2mya(k)/t2mya(1),t2mye(k)/abs(t2mya(1))
557   continue
558   continue

      if(index.eq.2)stop

c------------------------------------------------------------------------c
c     comparison with perturbative part                                  c
c------------------------------------------------------------------------c

      const=3.1415926**4/3./2.**7

      write(2,*)
      write(2,*) 'scalar correlator'
      write(2,*) ' x           s(x)        dels(x)     log(s(x))'
      do k=2,ndelt
         x=(k-1)*delstep
         cc=const*x**8/g(x)**4
         sub =scmya(k)-scdisc
         sube=scmye(k)+scdisce
         scal=1.+ sub*cc
         sp=1+sube/(abs(sub)+1.0/cc)
         sm=1-sube/(abs(sub)+1.0/cc)
         sclog=reglog(scal)
         scl1=reglog(sp)
         scl2=-reglog(sm)
         write(2,572) x,scal,scmye(k)*cc,sclog,scl1,scl2
      enddo

c------------------------------------------------------------------------c
c     pseudoscalar, psfree=-scfree                                       c
c------------------------------------------------------------------------c

      write(2,*)
      write(2,*) 'pseudoscalar correlator'
      write(2,*) ' x           p(x)        delp(x)     log(p(x))'

      do  k=2,ndelt
          x=(k-1)*delstep
          cc=const*x**8/g(x)**4
          pscal=1.- psmya(k)*cc
          psp=1+psmye(k)/(psmya(k)+1.0/cc)
          psm=1-psmye(k)/(psmya(k)+1.0/cc)
          psclog=reglog(pscal)
          pscl1=reglog(psp)
          pscl2=-reglog(psm)
          write(2,572) x,pscal,psmye(k)*cc,psclog,pscl1,pscl2
       enddo

c------------------------------------------------------------------------c
c     tensor, tfree=1/4*scfree                                           c
c------------------------------------------------------------------------c

      write(2,*)
      write(2,*) 'tensor correlator'
      write(2,*) ' x           t(x)        delt(x)     log(t(x))'

      do k=2,ndelt
         x=(k-1)*delstep
         cc=const*x**8/g(x)**4
         cc=cc/4.
         t1cal=1.+ t1mya(k)*cc
         tsp=1+t1mye(k)/(t1mya(k)+1.0/cc)
         tsm=1-t1mye(k)/(t1mya(k)+1.0/cc)
         t1clog=reglog(t1cal)
         t1cl1=reglog(tsp)
         t1cl2=-reglog(tsm)
         write(2,572) x,t1cal,t1mye(k)*cc,t1clog,t1cl1,t1cl2
      enddo

      write(2,*)
      write(2,*) 'tensor 2'
      write(2,*) ' x           t2(x)       delt2(x)    log(t2(x))'

      do k=2,ndelt
         x=(k-1)*delstep
         cc=const*x**8/g(x)**4
         cc=cc/4.
         t2cal= t2mya(k)*cc
         tsp=1+t2mye(k)/t2mya(k)
         tsm=1-t2mye(k)/t2mya(k)
         t2clog=reglog(t2cal)
         t2cl1=reglog(tsp)
         t2cl2=-reglog(tsm)
         write(2,572) x,t2cal,t2mye(k)*cc,t2clog,t2cl1,t2cl2
      enddo

      write(2,*)
      write(2,*) 'scalar correlator (modified subtraction)'
      write(2,*) ' x           s(x)        dels(x)     log(s(x))'
      do k=2,ndelt
         x=(k-1)*delstep
         cc=const*x**8/g(x)**4
         sub =scmya(k)-scdisc2
         sube=scmye(k)+scdisc2e
         scal=1.+ sub*cc
         sp=1+sube/(abs(sub)+1.0/cc)
         sm=1-sube/(abs(sub)+1.0/cc)
         sclog=reglog(scal)
         scl1=reglog(sp)
         scl2=-reglog(sm)
         write(2,572) x,scal,scmye(k)*cc,sclog,scl1,scl2
      enddo

      write(2,*)
      write(2,*) 'scalar correlator (modified subtraction)'
      write(2,*) ' x           s(x)        dels(x)     log(s(x))'
      do k=2,ndelt
         x=(k-1)*delstep
         cc=const*x**8/g(x)**4
         sub =scmya(k)-scdisc3
         sube=scmye(k)+scdisc3e
         scal=1.+ sub*cc
         sp=1+sube/(abs(sub)+1.0/cc)
         sm=1-sube/(abs(sub)+1.0/cc)
         sclog=reglog(scal)
         scl1=reglog(sp)
         scl2=-reglog(sm)
         write(2,572) x,scal,scmye(k)*cc,sclog,scl1,scl2
      enddo

      write(2,*)
      write(2,*) 'perturbative'
      write(2,*) ' x           pi(x)       pi(x)/g^4   g(x)'

      do k=2,ndelt
         x=(k-1)*delstep
         cc0= 1.0/(const*x**8)
         gg = g(x)
         cc = cc0*gg**4
         write(2,573) x,cc,cc0,gg
      enddo

572   format(1x,7(1x,g11.5))
573   format(1x,4(1x,g11.5))
571   continue

      stop

 1021 format(i5,f8.3,f9.4,f8.4,f9.4,f8.4,f9.4,f8.4)
 1023 format(i5,f6.3,3(f12.4,f11.3))
 1022 format(1x, i5, a70)

      end

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      function g(x)
c------------------------------------------------------------------------c
c     modified running coupling constant                                 c
c------------------------------------------------------------------------c
      common /param/ a, alpha,rh0,sg,dz,drh, nc, nf, rmu, rms
      xlam = 1.0
      pi = 3.1415926
      b  = 11.0/2.0-nf/3.0
      alp= pi/b/log(6.7+1.0/(xlam*x))
      g  = sqrt(4.0*pi*alp)
      return
      end

c------------------------------------------------------------------------
c------------------------------------------------------------------------


      function reglog(x)
      if (x.gt.0.) then
         reglog=alog10(x)
      else
         reglog=-20.
      endif
      return
      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine mult2(nin,a,b,c)
      parameter(ni=16)

      complex a(ni,2,2), b(ni,2,2),c(ni,2,2)

      do 10 is = 1, nin
         c(is,1,1) = a(is,1,2)*b(is,2,1)+a(is,1,1)*b(is,1,1)
         c(is,1,2) = a(is,1,2)*b(is,2,2)+a(is,1,1)*b(is,1,2)
         c(is,2,1) = a(is,2,2)*b(is,2,1)+a(is,2,1)*b(is,1,1)
         c(is,2,2) = a(is,2,2)*b(is,2,2)+a(is,2,1)*b(is,1,2)
   10 continue

      return
      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine mult3(nin,a,b,c)

      parameter (ni=16)
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

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine taumat
c------------------------------------------------------------------------c
c     four-dim. tau matrices, the first number is index                  c
c------------------------------------------------------------------------c
      complex tau
      common /taum/tau(2,4,2,2)
       do 1 ind=1,2
       index=3-2*ind
         do 2 mu=1,4
         do 2 k1=1,2
         do 2 k2=1,2
         tau(ind,mu,k1,k2)=(0.,0.)
2        continue
       tau(ind,4,1,1)=(0.,-1.)*index
       tau(ind,4,2,2)=(0.,-1.)*index
       tau(ind,3,1,1)=(1.,0.)
       tau(ind,3,2,2)=(-1.,0.)
       tau(ind,1,1,2)=(1.,0.)
       tau(ind,1,2,1)=(1.,0.)
       tau(ind,2,1,2)=(0.,-1.)
       tau(ind,2,2,1)=(0.,1.)
1      continue
      return
      end

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine tordis(nin,x,y,z,rr)

      parameter (n2=16)
      common /box/ alb(4)
      dimension x(4),y(n2,5),z(n2,4),rr(n2)

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
   10   continue
 1    continue

      return
      end

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine dtordis(nin,x,y,z,rr)

      parameter (n2=16)
      common /box/ alb(4)
      real*4 y(n2,5)
      real*8 x(4),z(n2,4),rr(n2)
      real*8 dis,adis,asg,al2

      do 5 is = 1, nin
        rr(is) = 0.
   5  continue
       do 1 m=1,4
         do 10 is = 1, nin
           dis=x(m)-y(is,m)
           adis=abs(dis)
           al2 = dble(alb(m))/2.d0
           asg = 1+sign(1.d0,adis-al2)
           dis = dis - alb(m)/2*asg*sign(1.d0,dis)
           rr(is) = rr(is) + dis**2
           z(is,m)=dis
   10   continue
 1    continue

      return
      end

c---------------------------------------------------------------------------
c---------------------------------------------------------------------------

      subroutine su3(n,nin,nd,x,y,u)

      parameter(n2=16)
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

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine setup(n,nin,nd,zr,e1,e2,iread,icon)
c------------------------------------------------------------------------c
c     generate new instanton configuration                               c
c------------------------------------------------------------------------c
c     n          number of colors                                        c
c     nin        number of instantons                                    c
c     nd         array dimension (nin<nd)                                c
c     zr,e1,e2   instanton configuration                                 c
c     iread      random config or input from file (0,1)                  c
c     icon       switch                                                  c
c------------------------------------------------------------------------c

      complex u0
      character*9 a8
      dimension falb(4)
      dimension zr(nd,5), er1(6), er2(6)
      dimension x(6), y(6), z(6), e1(nd,6), e2(nd,6)

      common /param/   a,alpha,rh0, sg, dz, drh, nc, nf, rmu, rms
      common /counter/ icount
      common /box/     alb(4)
      common /nconf/   nconfig

      if (iread.eq.1) then

c------------------------------------------------------------------------c
c     input from file                                                    c
c------------------------------------------------------------------------c

      if (icon.eq. 0) then

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

 997  format('Warning: box size has changed')
 998  format('Warning: nf has changed')
 999  format('Warning: masses have changed')

c------------------------------------------------------------------------c
c     read new configuration                                             c
c------------------------------------------------------------------------c

      read (4,*) icon
      print 207, icon
      do 25 i = 1, nin
         read(4,101) (zr(i,k),k=1,5)
         read(4,101) (e1(i,k), k=1,6)
         read(4,101) (e2(i,k), k=1,6)
   25 continue

  207 format(1x, ' configuration ', i5, '   has been read ',/)
  101 format(1x,6f12.5)

      else

c------------------------------------------------------------------------c
c     random configuration                                               c
c------------------------------------------------------------------------c

      icount = icount + 1
      print *, ' Configuration: ', icount
  701 format(1x,/,1x,'  Configuration: ', i5)

c------------------------------------------------------------------------c
c     loop over instantons                                               c
c------------------------------------------------------------------------c

      do 100 ip = 1, nin

         do 10 i = 1, 4
            zr(ip,i) = alb(i)*rang( )
   10    continue

         zr(ip,5) = rh0
         do 20 i = 1, 6
            e1(ip,i) = 0.0
            e2(ip,i) = 0.0
   20    continue
         e1(ip,1) = 1.0
         e2(ip,3) = 1.0
         call rsu(6,er1,er2)
         do 70 i = 1, 6
            e1(ip,i) = er1(i)
            e2(ip,i) = er2(i)
   70    continue

  100 continue

      end if

      return
      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine rsu(n, x, y)

      dimension x(n), y(n), xs(6)

      call r6(x, n)

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

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine rdot(n, x, y, xdy)

      dimension x(n), y(n)

      xdy = 0.0
      do 10 i = 1, n
        xdy = xdy + x(i)*y(i)
   10 continue

      return
      end

c-------------------------------------------------------------------------
c-------------------------------------------------------------------------

      subroutine r6(x,n)

      dimension x(n)

 100  continue
      r2 = 0.0
      do 10 i = 1, 6
        x(i) = 2*rang( ) - 1
        r2 = r2+x(i)*x(i)
  10  continue
      if (r2 .gt. 1) goto 100
      r = sqrt(r2)
      do 20 i = 1, 6
        x(i) = x(i)/r
   20 continue

      return
      end

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine myaddto(x,xtot,x2tot)
      xtot=xtot+x
      x2tot=x2tot+x*x
      return
      end

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine mydisp(n,xtot,x2tot,xav,xerr)
      if(n.lt.1)goto 10
      xav=xtot/(n*1.)
      del2=x2tot/(1.*n*n)-xav*xav/(1.*n)
      if(del2.lt.0.)del2=0.
      xerr=sqrt(del2)
10    continue
      return
      end

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine myzero(n,arr)
      dimension arr(n)
      do 1 k=1,n
         arr(n)=0.
1     continue
      return
      end

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine glue(ni,ndd,zr,e1,e2,x,sc,ps,t1,t2)
c-------------------------------------------------------------------------c
c     calculate gluonic observables                                       c
c-------------------------------------------------------------------------c
c     ni,ndd    number of instantons, array dimension                     c
c     zr,e1,e2  configuration                                             c
c     x(4)      point where operators are evaluated                       c
c     sc        G^2                                                       c
c     ps        G\tilde G                                                 c
c     t1        (E^2-B^2)                                                 c
c     t2        2B_3^2-B_1^2-B_2^2                                        c
c-------------------------------------------------------------------------c

      complex ct2,es,bs,eb
      complex g(3,3,4,4)
      dimension  zr(ndd,5),e1(ndd,6),e2(ndd,6),x(4)
      complex e(3,3,3),b(3,3,3)

c-------------------------------------------------------------------------c
c     calculate field strength tensor                                     c
c-------------------------------------------------------------------------c

      call dfmunu(ni,zr,e1,e2,x,g,f2,e,b,ee,bb)

      es = ( 0.,0.)
      bs = ( 0.,0.)
      eb = ( 0.,0.)
      ct2= ( 0.,0.)

      do 1 k1=1,3
      do 1 k2=1,3
         ct2=ct2+2*g(k1,k2,1,2)*g(k2,k1,1,2)
     2            -g(k1,k2,3,2)*g(k2,k1,3,2)
     3            -g(k1,k2,3,1)*g(k2,k1,3,1)
         do 2 m=1,3
            eb=eb+e(k1,k2,m)*b(k2,k1,m)
2        continue
1     continue

      ps = 4.0*real(eb)
      sc = real(f2)
      t1 = real(ee-bb)
      t2 = real(ct2)

c thus scalar is 2(e^2+b^2) while tensor  t1 is without 2, and with minus

      return
      end

c---------------------------------------------------------------------------
c---------------------------------------------------------------------------

      subroutine sglue(ni,ndd,zr,e1,e2,x1,x2,sc,ps,t1,t2)
c-------------------------------------------------------------------------c
c     calculate point split gluonic operators for glueball wavefunctions  c
c-------------------------------------------------------------------------c
c     ni,ndd    number of instantons, array dimension                     c
c     zr,e1,e2  configuration                                             c
c     x1,x2     points where operators are evaluated                      c
c     sc        G^2                                                       c
c     ps        G\tilde G                                                 c
c     t1        (E^2-B^2)                                                 c
c     t2        2B_3^2-B_1^2-B_2^2                                        c
c-------------------------------------------------------------------------c

      complex ct2,es,bs,eb
      complex g(3,3,4,4)
      dimension  zr(ndd,5),e1(ndd,6),e2(ndd,6),x1(4),x2(4)
      complex ef1(3,3,3),bf1(3,3,3),ef2(3,3,3),bf2(3,3,3)

c-------------------------------------------------------------------------c
c     calculate field strength at x1,x2                                   c
c-------------------------------------------------------------------------c

      call dfmunu(ni,zr,e1,e2,x1,g,f2,ef1,bf1,ee1,bb1)
      call dfmunu(ni,zr,e1,e2,x2,g,f2,ef2,bf2,ee2,bb2)

      es = ( 0.,0.)
      bs = ( 0.,0.)
      eb = ( 0.,0.)
      ct2= ( 0.,0.)

      do 1 k1=1,3
      do 1 k2=1,3
         do 2 m=1,3
            eb=eb+ef1(k1,k2,m)*bf2(k2,k1,m)+ef2(k1,k2,m)*bf1(k2,k1,m)
            es=es+ef1(k1,k2,m)*ef2(k2,k1,m)
            bs=bs+bf1(k1,k2,m)*bf2(k2,k1,m)
2        continue
         ct2=ct2+2*bf1(k1,k2,1)*bf2(k2,k1,1)
     2         - bf1(k1,k2,2)*bf2(k2,k1,2)-bf1(k1,k2,3)*bf2(k2,k1,3)
1     continue

      sc = real(es+bs)
      ps = 4.0*real(eb)
      t1 = real(es-bs)
      t2 = real(ct2)

      return
      end

c---------------------------------------------------------------------------
c---------------------------------------------------------------------------

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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      function ran2(idum)
c---------------------------------------------------------------------------c
c     numerical recipes random number generator, set iseed to a negative    c
c     value to initialize sequence. on first call, sequence is initialized  c
c     automatically.                                                        c
c---------------------------------------------------------------------------c
c     cray requires save statement !                                        c
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
      parameter(n2=16)
c no difference between ndd and n2 should be there!
      dimension  zr(ndd,5),e1(ndd,6),e2(ndd,6),zz(n2,4),zs(n2)
      dimension y(4),z(4),a(3,3,4),xa(4),ai3(3,3,4),ai(3,4)
      complex a,ai3,uuu(n2,3,3),u1(3,3),u2(3,3)

c-------------------------------------------------------------------------c
c     modify ratio ansatz by exponential tail                             c
c-------------------------------------------------------------------------c

c        fdecay=0.2
c temporally absent
         fdecay=0.

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

c---------------------------------------------------------------------+-----
c---------------------------------------------------------------------+-----

      subroutine dpotsu3(ni,ndd,zr,e1,e2,xa,a)
c-------------------------------------------------------------------------c
c     evaluate vector potential a(i,j,mu) at point xa(k). collective      c
c     variables in zr,e1,e2. first ni/2 pseudoparticles are instantons    c
c     second half antiinstantons. vector potential is calculated in the   c
c     ratio ansatz.                                                       c
c-------------------------------------------------------------------------c
c     this version uses antihermitean generators.                         c
c-------------------------------------------------------------------------c
c              double precision !                                         c
c-------------------------------------------------------------------------c
c     note: on in/output only xa,a are double precision !                 c
c-------------------------------------------------------------------------c
c     ni        number of instantons                                      c
c     ndd       maximum number of instantons (for array dimensions)       c
c     zr(i,5)   position, size of instanton i                             c
c     e1(i,6)   real, complex parts of first row of su(3) matrix          c
c     e2(i,6)   same for second row                                       c
c     xa(k)     coordinates of point                                      c
c     a(i,j,mu) vector potential at xa(k) (output)                        c
c-------------------------------------------------------------------------c
c     note: make sure that n2=ndd !                                       c
c-------------------------------------------------------------------------c
      parameter(n2=16)

      real*4 zr(ndd,5),e1(ndd,6),e2(ndd,6)
      real*8 zz(n2,4),zs(n2)
      real*8 y(4),z(4),xa(4),ai(3,4)
      real*8 ros, ff, f, z4, fdecay, fact, rr
      complex*8  uuu(n2,3,3)
      complex*16 a(3,3,4),ai3(3,3,4),u1(3,3),u2(3,3)

c-------------------------------------------------------------------------c
c     modify ratio ansatz by exponential tail                             c
c-------------------------------------------------------------------------c

c     fdecay=0.2
      fdecay=0.d0

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

      call dtordis(ni,xa,zr,zz,zs)

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

      dimension  zr(ndd,5),e1(ndd,6),e2(ndd,6),zz(512,4),zs(512)
      dimension y(4),z(4),a(3,3,4),xa(4),ai3(3,3,4),ai(3,4)
      complex a,ai3,uuu(512,3,3),u1(3,3),u2(3,3)

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

c--------------------------------------------------------------------------
c--------------------------------------------------------------------------

      subroutine dsumsu3(ni,ndd,zr,e1,e2,xa,a)
c-------------------------------------------------------------------------c
c     evaluate vector potential a(i,j,mu) at point xa(k). collective      c
c     variables in zr,e1,e2. first ni/2 pseudoparticles are instantons    c
c     second half antiinstantons. vector potential is calculated in the   c
c     sum ansatz.                                                         c
c-------------------------------------------------------------------------c
c     this version assumes antihermitean generators.                      c
c-------------------------------------------------------------------------c
c              double precision !                                         c
c-------------------------------------------------------------------------c
c     note: on in/output only xa,a are double precision !                 c
c-------------------------------------------------------------------------c
c     ni        number of instantons                                      c
c     ndd       maximum number of instantons (for array dimensions)       c
c     zr(i,5)   position, size of instanton i                             c
c     e1(i,6)   real, complex parts of first row of su(3) matrix          c
c     e2(i,6)   same for second row                                       c
c     xa(k)     coordinates of point                                      c
c     a(i,j,mu) vector potential at xa(k) (output)                        c
c-------------------------------------------------------------------------c
c     note: make sure that n2=ndd !                                       c
c-------------------------------------------------------------------------c
      parameter(n2=16)

      real*4 zr(ndd,5),e1(ndd,6),e2(ndd,6)
      real*8 zz(n2,4),zs(n2)
      real*8 y(4),z(4),xa(4),ai(3,4)
      real*8 ros, ff, f, z4, fdecay, fact, rr
      complex*8  uuu(n2,3,3)
      complex*16 a(3,3,4),ai3(3,3,4),u1(3,3),u2(3,3)

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

      call dtordis(ni,xa,zr,zz,zs)

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
c     call sumsu3(ni,ndd,zr,e1,e2,xa,a)
c switch to sum ansatz
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

      parameter(nd=16)
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

c     call sumsu3(ni,nd,zr,e1,e2,x,a)
      call potsu3(ni,nd,zr,e1,e2,x,a)

c-------------------------------------------------------------------------c
c     calculate vector potential on star around central point             c
c-------------------------------------------------------------------------c
c     mu labels direction of dx                                           c
c-------------------------------------------------------------------------c

      del = 1.0e-4
      do 10 mu=-4,4
      if(mu.eq.0)go to 10

         do 20 nu=1,4
            xx(nu) = x(nu)
 20      continue
         mup = abs(mu)
         xx(mup) = xx(mup) + del*sign(1,mu)
c        call sumsu3(ni,nd,zr,e1,e2,xx,da)
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

c---------------------------------------------------------------------+------
c---------------------------------------------------------------------+------

      subroutine dfmunu(ni,zr,e1,e2,x,f,f2,e,b,ee,bb)
c-------------------------------------------------------------------------c
c     evaluate field strength tensor for arbitrary gauge potential given  c
c     subroutine potsu3(). gauge field is assumed to be an antihermitean  c
c     vector field A_\mu=A_\mu^a(\tau^a/2i) with generators normalized    c
c     as tr(\tau^a\tau^b)=2\delta^{ab}.                                   c
c-------------------------------------------------------------------------c
c             double precision !                                          c
c-------------------------------------------------------------------------c
c     Note: double precision only internally, in/out in single precision  c
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

      parameter(nd=16)

      integer eps(3,2)
      real*4  zr(nd,5), e1(nd,6), e2(nd,6), x(4)
      real*4  f2, ee, bb
      real*8  xx(4), del
      complex*8  f(3,3,4,4), e(3,3,3), b(3,3,3)
      complex*8  cf2, cee, cbb
      complex*16 a(3,3,4), aa(3,3,-4:4,4), da(3,3,4)
      complex*16 ff(3,3,4,4)

      data eps/ 2,3,1, 3,1,2 /

c-------------------------------------------------------------------------c
c     vector potential at x                                               c
c-------------------------------------------------------------------------c

      do 5 nu=1,4
         xx(nu) = x(nu)
  5   continue

c     call sumsu3(ni,nd,zr,e1,e2,xx,a)
      call dpotsu3(ni,nd,zr,e1,e2,xx,a)

c-------------------------------------------------------------------------c
c     calculate vector potential on star around central point             c
c-------------------------------------------------------------------------c
c     mu labels direction of dx                                           c
c-------------------------------------------------------------------------c

      del = 1.0d-5
      do 10 mu=-4,4
      if(mu.eq.0)go to 10

         do 20 nu=1,4
            xx(nu) = x(nu)
 20      continue
         mup = abs(mu)
         xx(mup) = xx(mup) + del*sign(1,mu)
c        call sumsu3(ni,nd,zr,e1,e2,xx,da)
         call dpotsu3(ni,nd,zr,e1,e2,xx,da)
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


