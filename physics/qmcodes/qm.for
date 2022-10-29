      program main
c----------------------------------------------------------------------------c
c     lattice calculation in quantum mechanics.                              c
c----------------------------------------------------------------------------c
c     version :            1.0                                               c
c     creation date :   03-15-98                                             c
c     last modified :   03-25-98                                             c
c----------------------------------------------------------------------------c
c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
c----------------------------------------------------------------------------c
c     lattice x(i=0,n), peridic b.c. x(0)=x(n)                               c
c----------------------------------------------------------------------------c
c     nmc    number of mc sweeps                                             c
c     nc     number of correlator measurements in a single configuration     c
c     neq    number of equlibration sweeps                                   c
c     kp     number of sweeps between writeout of complete configuration     c
c----------------------------------------------------------------------------c

      parameter(npar=10002,nmax=100)
      real x(-1:npar)
      integer ix(100)
      real xcor_sum(nmax), xcor2_sum(nmax)
      real xcor_av(nmax),xcor_er(nmax)
      real x2cor_sum(nmax),x2cor2_sum(nmax)
      real x2cor_av(nmax), x2cor_er(nmax)
      real x3cor_sum(nmax),x3cor2_sum(nmax)
      real x3cor_av(nmax), x3cor_er(nmax)
      real x2sub_av(nmax), x2sub_er(nmax)
      
      common /par/ f,a
      common /seed/iseed

      open(unit=16,file='qm.dat',status='unknown')
      open(unit=17,file='config.dat',status='unknown')
      open(unit=18,file='trajectory.dat',status='unknown')
      open(unit=19,file='qmdist.dat',status='unknown')
      open(unit=20,file='qmcor.dat',status='unknown')
      open(unit=21,file='qmcor2.dat',status='unknown')
      open(unit=22,file='qmcor3.dat',status='unknown') 

      write(6,*) 'separation of wells f (f=1.4)'
      read(5,*)  f
      write(6,*) 'grid size n<1000 (n=100)'
      read(5,*)  n
      write(6,*) 'grid spacing a (dtau=0.05)'
      read(5,*)  a
      write(6,*) 'cold/hot start (0,1)'
      read(5,*)  icold
      write(6,*) 'equilibration sweeps'
      read(5,*)  neq
      write(6,*) 'monte carlo sweeps'
      read(5,*)  nmc   
      write(6,*) 'update x (delx)'
      read(5,*)  delx     
      write(6,*) 'number of points in correlator'
      read(5,*)  np 
      write(6,*) 'number of measurements per config'
      read(5,*)  nc 
      write(6,*) 'write every kth config'
      read(5,*)  kp
      tmax = n*a

c----------------------------------------------------------------------------c
c     echo input parameters                                                  c
c----------------------------------------------------------------------------c

      write(16,*)   'lattice qm 1.0'
      write(16,*)   '--------------'
      write(16,101) f,n,a
      write(16,102) nmc,neq
      write(16,103) np,nc
      write(16,104) delx,icold
c     write(17,*)   n,nmc/kp,n*a,f
      
 101  format(1x,' f    = ',f8.2,' n   = ',i5,' a   = ',f5.4)
 102  format(1x,' nmc  = ',i8  ,' neq = ',i5)
 103  format(1x,' np    = ',i8  ,' nc  = ',i5)
 104  format(1x,' delx = ',f8.2,' icol= ',i5)
      
c----------------------------------------------------------------------------c
c     parameters for histograms                                              c
c----------------------------------------------------------------------------c
      
      nxhist = 50
      xhist_min = -1.5*f
      stxhist= 3.0*f/float(nxhist)

c----------------------------------------------------------------------------c
c     clear summation arrays                                                 c
c----------------------------------------------------------------------------c

      stot_sum = 0.0
      stot2_sum= 0.0
      vtot_sum = 0.0
      vtot2_sum= 0.0
      ttot_sum = 0.0
      ttot2_sum= 0.0
      tvir_sum = 0.0
      tvir2_sum= 0.0
      x_sum = 0.0
      x2_sum= 0.0
      x4_sum= 0.0
      x8_sum= 0.0
    
      call zero(xcor_sum,nmax)
      call zero(xcor2_sum,nmax)
      call zero(x2cor_sum,nmax)
      call zero(x2cor2_sum,nmax)
      call zero(x3cor_sum,nmax)
      call zero(x3cor2_sum,nmax)
      call izero(ix,nxhist)

c----------------------------------------------------------------------------c
c     initialize                                                             c
c----------------------------------------------------------------------------c

      iseed =-1234

      do 10 i=1,n
         if (icold .eq.0) then
            x(i) = -f
         else 
            x(i) = 2.0*ran2(iseed)*f-f
         endif
  10  continue       
  
c----------------------------------------------------------------------------c
c     periodic boundary conditions                                           c
c----------------------------------------------------------------------------c

      x(n)   = x(1)
      x(n+1) = x(2)
      x(0)   = x(n-1)
      x(-1)  = x(n-2)
      
      stot = 0.0

c----------------------------------------------------------------------------c
c     initial action                                                         c
c----------------------------------------------------------------------------c
           
      do 50 i=1,n 
         xp = (x(i+1)-x(i))/a
         t  = 1.0/4.0*xp**2
         v  = (x(i)**2-f**2)**2
         s  = a*(t+v)
         stot = stot + s
  50  continue   
  
      nacc = 0
      nhit = 0    
      nconf= 0
      ncor = 0

c----------------------------------------------------------------------------c
c     monte carlo sweeps                                                     c
c----------------------------------------------------------------------------c

      do 100 i=1,nmc
         
         nconf = nconf+1
              
         if(i .eq. neq)then
            ncor = 0
            nconf = 0
            stot_sum = 0.0
            stot2_sum= 0.0
            vtot_sum = 0.0
            vtot2_sum= 0.0
            ttot_sum = 0.0
            ttot2_sum= 0.0
            tvir_sum = 0.0
            tvir2_sum= 0.0
            x_sum = 0.0
            x2_sum= 0.0
            x4_sum= 0.0
            x8_sum= 0.0
            call zero(xcor_sum,nmax)
            call zero(xcor2_sum,nmax)
            call zero(x2cor_sum,nmax)
            call zero(x2cor2_sum,nmax)
            call zero(x3cor_sum,nmax)
            call zero(x3cor2_sum,nmax)
            call izero(ix,nxhist)
            call izero(iz,nzhist)
         endif
         
c----------------------------------------------------------------------------c
c     one sweep thorough configuration                                       c
c----------------------------------------------------------------------------c
 
         do 200 j=1,n-1
            
            nhit = nhit+1  
            
            xpm = (x(j)-x(j-1))/a
            xpp = (x(j+1)-x(j))/a
            t = 1.0/4.0*(xpm**2+xpp**2)
            v = (x(j)**2-f**2)**2
            sold = a*(t+v)

            xnew = x(j) + delx*(2.0*ran2(iseed)-1.0)
             
            xpm = (xnew-x(j-1))/a
            xpp = (x(j+1)-xnew)/a
            t = 1.0/4.0*(xpm**2+xpp**2)
            v = (xnew**2-f**2)**2
            snew = a*(t+v)
            dels = snew-sold           
                       
            p  = ran2(iseed)
            dels = min(dels,70.0)
            dels = max(dels,-70.0)
            
            if (exp(-dels) .gt. p) then
               x(j) = xnew
               nacc = nacc + 1
            endif
            
 200     continue        
 
         x(n)   = x(1)
         x(n+1) = x(2)
         x(0)   = x(n-1)
         x(-1)  = x(n-2)
         
c----------------------------------------------------------------------------c
c     calculate action etc.                                                  c
c----------------------------------------------------------------------------c

         stot = 0.0
         ttot = 0.0
         tvtot= 0.0
         vtot = 0.0  
         do 150 j=1,n 
            xp = (x(j+1)-x(j))/a
            t  = 1.0/4.0*xp**2
            v  = (x(j)**2-f**2)**2
            tv = 2.0*x(j)**2*(x(j)**2-f**2)
            s  = a*(t+v)
            ttot = ttot +a*t
            vtot = vtot +a*v
            tvtot= tvtot+a*tv
            stot = stot + s
 150     continue   
         write(18,*) i,stot,ttot,vtot
                                           
         if (mod(i,kp) .eq. 0) then
            write(6,*)
            write(6,*) 'configuration   ',i
            write(6,*) 'acceptance rate ',float(nacc)/float(nhit)
            write(6,*) 'action (T,V)    ',stot,ttot,vtot         
c           write(17,*) 'configuration',i
            write(17,*)
            do 155 k=1,n
               write(17,*) k*a,x(k)
 155        continue
         endif

c----------------------------------------------------------------------------c
c     include in sample                                                      c
c----------------------------------------------------------------------------c
         
         stot_sum = stot_sum + stot
         stot2_sum= stot2_sum+ stot**2
         vtot_sum = vtot_sum + vtot
         vtot2_sum= vtot2_sum+ vtot**2
         ttot_sum = ttot_sum + ttot
         ttot2_sum= ttot2_sum+ ttot**2
         tvir_sum = tvir_sum + tvtot
         tvir2_sum= tvir2_sum+ tvtot**2
                                                   
         do 160 k=1,n
            call lens(x(k),xhist_min,stxhist,nxhist,ix)
            x_sum = x_sum + x(k)
            x2_sum= x2_sum+ x(k)**2
            x4_sum= x4_sum+ x(k)**4
            x8_sum= x8_sum+ x(k)**8
 160     continue
         
c----------------------------------------------------------------------------c
c     correlation function                                                   c
c----------------------------------------------------------------------------c

         do 250 ic=1,nc
            
            ncor = ncor + 1 
            ip0 = int( (n-np)*ran2(iseed) ) 
            x0  = x(ip0) 
            
            do 280 ip=1,np
            
               x1 = x(ip0+ip)
               xcor = x0*x1
               x2cor= xcor**2
               x3cor= xcor**3
               
               xcor_sum(ip)  = xcor_sum(ip)  + xcor
               xcor2_sum(ip) = xcor2_sum(ip) + xcor**2
               x2cor_sum(ip) = x2cor_sum(ip) + x2cor
               x2cor2_sum(ip)= x2cor2_sum(ip)+ x2cor**2
               x3cor_sum(ip) = x3cor_sum(ip) + x3cor
               x3cor2_sum(ip)= x3cor2_sum(ip)+ x3cor**2

 280        continue
         
 250     continue
         
c----------------------------------------------------------------------------c
c     next configuration                                                     c
c----------------------------------------------------------------------------c

 100  continue
                   
c----------------------------------------------------------------------------c
c     averages                                                               c
c----------------------------------------------------------------------------c

      call disp(nconf,stot_sum,stot2_sum,stot_av,stot_err)
      call disp(nconf,vtot_sum,vtot2_sum,vtot_av,vtot_err)
      call disp(nconf,ttot_sum,ttot2_sum,ttot_av,ttot_err)
      call disp(nconf,tvir_sum,tvir2_sum,tvir_av,tvir_err)
      call disp(nconf*n,x_sum,x2_sum,x_av,x_err)
      call disp(nconf*n,x2_sum,x4_sum,x2_av,x2_err)
      call disp(nconf*n,x4_sum,x8_sum,x4_av,x4_err)
      
      do 300 ip=1,np
      call disp(ncor,xcor_sum(ip),xcor2_sum(ip),
     &          xcor_av(ip),xcor_er(ip))
      call disp(ncor,x2cor_sum(ip),x2cor2_sum(ip),
     &          x2cor_av(ip),x2cor_er(ip))
      call disp(ncor,x3cor_sum(ip),x3cor2_sum(ip),
     &          x3cor_av(ip),x3cor_er(ip))
 300  continue

      v_av  = vtot_av/tmax
      v_err = vtot_err/tmax
      t_av  = ttot_av/tmax
      t_err = ttot_err/tmax
      tv_av = tvir_av/tmax
      tv_err= tvir_err/tmax
      e_av  = v_av+tv_av
      e_err = sqrt(v_err**2+tv_err**2)

c----------------------------------------------------------------------------c
c     output                                                                 c
c----------------------------------------------------------------------------c

      write(16,*)
      write(16,901) stot_av,stot_err
      write(16,902) v_av,v_err
      write(16,903) t_av,t_err
      write(16,904) tv_av,tv_err
      write(16,905) e_av,e_err
      write(16,906) x_av,x_err
      write(16,907) x2_av,x2_err
      write(16,908) x4_av,x4_err
      write(16,*) 
      
c----------------------------------------------------------------------------c
c     correlation function, log derivative                                   c
c----------------------------------------------------------------------------c

      write(16,*) '# x correlation function'
      write(20,*) '# tau       x(tau)       dx(tau)     dlog'
      do 500 ip=1,np-1
         dx  = (xcor_av(ip)-xcor_av(ip+1))/xcor_av(ip)/a
         dxe2= (xcor_er(ip+1)/xcor_av(ip))**2
     &        +(xcor_er(ip)*xcor_av(ip+1)/xcor_av(ip)**2)**2
         dxe = sqrt(dxe2)/a
         write(16,555) ip*a,xcor_av(ip),xcor_er(ip),dx,dxe
         write(20,555) ip*a,xcor_av(ip),xcor_er(ip),dx,dxe
 500  continue        

c----------------------------------------------------------------------------c
c     subtracted x^2 correlation function, log derivative                    c
c----------------------------------------------------------------------------c

      xx_sub = x2cor_av(np)
      xx_er  = x2cor_er(np)
      do 505 ip=1,np
         x2sub_av(ip) = x2cor_av(ip)-xx_sub
         x2sub_er(ip) = sqrt(x2cor_er(ip)**2+xx_er**2)
 505  continue

      write(16,*) '# x2 correlation function'
      write(21,*) '# tau       x2(tau)      dx2(tau)     dlog'
      do 510 ip=1,np-1
         dx  = (x2sub_av(ip)-x2sub_av(ip+1))/x2sub_av(ip)/a
         dxe2= (x2sub_er(ip+1)/x2sub_av(ip))**2
     &        +(x2sub_er(ip)*x2sub_av(ip+1)/x2sub_av(ip)**2)**2
         dxe = sqrt(dxe2)/a
         write(16,555) ip*a,x2cor_av(ip),x2cor_er(ip),dx,dxe
         write(21,555) ip*a,x2cor_av(ip),x2cor_er(ip),dx,dxe
 510  continue  

c----------------------------------------------------------------------------c
c     x^3 correlation function, log derivative                               c
c----------------------------------------------------------------------------c

      write(16,*) '# x^3 correlation function'
      write(22,*) '# tau       x(tau)       dx(tau)     dlog'
      do 520 ip=1,np-1
         dx  = (x3cor_av(ip)-x3cor_av(ip+1))/x3cor_av(ip)/a
         dxe2= (x3cor_er(ip+1)/x3cor_av(ip))**2
     &        +(x3cor_er(ip)*x3cor_av(ip+1)/x3cor_av(ip)**2)**2
         dxe = sqrt(dxe2)/a
         write(16,555) ip*a,x3cor_av(ip),x3cor_er(ip),dx,dxe
         write(22,555) ip*a,x3cor_av(ip),x3cor_er(ip),dx,dxe
 520  continue   

c----------------------------------------------------------------------------c
c     wave function                                                          c
c----------------------------------------------------------------------------c

      write(16,*)
      write(16,*) ' x distribution '
      call lev(xhist_min,stxhist,nxhist,2,ix,16)
      write(19,*) ' x distribution '
      xnorm = 0
      do 600 i=1,nxhist 
         xnorm = xnorm + ix(i)*stxhist
 600  continue
      do 610 i=1,nxhist 
         xx = xhist_min + (i-1)*stxhist
         write(19,*) xx,ix(i)/xnorm
 610  continue
      
 901  format(1x,'stot  = ',f12.5,' +/- ',f12.5)
 902  format(1x,'v_av  = ',f12.5,' +/- ',f12.5)
 903  format(1x,'t_av  = ',f12.5,' +/- ',f12.5)
 904  format(1x,'t(vir)= ',f12.5,' +/- ',f12.5)
 905  format(1x,'e_av  = ',f12.5,' +/- ',f12.5)
 906  format(1x,'x     = ',f12.5,' +/- ',f12.5)
 907  format(1x,'x^2   = ',f12.5,' +/- ',f12.5)
 908  format(1x,'x^4   = ',f12.5,' +/- ',f12.5)
 222  format(1x,2(f12.5))
 333  format(1x,3(f12.5))
 444  format(1x,4(f12.5))
 555  format(1x,5(f12.5))
 666  format(1x,6(f12.5))
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
     &ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     &ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
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
   
c----------------------------------------------------------------------+-----
c----------------------------------------------------------------------+-----
 
      subroutine disp(n,xtot,x2tot,xav,xerr)
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

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine smooth(x,xs,npar,n,ns)
c----------------------------------------------------------------------------c
c     naive data smoothing with sliding window                               c
c----------------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      real*8 x(-npar:npar),xs(-npar:npar)
      
      do 10 i=-n,n
         xsum = 0.0
         do 20 j=i-ns,i+ns
            xsum = xsum + x(j)
  20     continue
         xav = xsum/float(2*ns+1)
         xs(i)= xav
  10  continue

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine cool(x,xs,npar,n,ncool)
c----------------------------------------------------------------------------c
c     local cooling algorithm                                                c
c----------------------------------------------------------------------------c
      implicit real (a-h,o-z)
      real x(-npar:npar),xs(-npar:npar)
      common /par/  f,dtau
      common /seed/ iseed
      
      nhit = 10
      delx = 0.001

      do 2  i=-n,n
         xs(i) = x(i)
  2   continue   

      do 5  k=1,ncool

      do 10 i=-n,n
         do 20 j=1,nhit
            xdot = (xs(i+1)-xs(i))/dtau
            v    = (xs(i)**2-f**2)**2
            sold = xdot**2/4. + v
            xnew = xs(i)+delx*(ran2(iseed)-0.5)
            xdot = (xs(i+1)-xnew)/dtau
            v    = (xnew**2-f**2)**2
            snew = xdot**2/4. + v
            if (snew .lt. sold) xs(i)=xnew
  20     continue
  10  continue

   5  continue

      return
      end
   
c---------------------------------------------------------------------+----
c---------------------------------------------------------------------+----

      subroutine lev(xmin,st,n,nsm,ist,idev)
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
      write(idev,30)k
   30 format(5x,i20,5i10)
      write(idev,40)
   40 format(2x,3h---,7(9(1h-),1h+))
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
      write(idev,70)i,ist(i),x,t
   70 format(1x,1x,i4,i7,1pe11.3,1x,50a1,1hi)
    6 d=d+ist(i)*(x-s)**2
      if(j.lt.2)goto 7
      d=sqrt(d/(j-1))
    7 write(idev,40)
      write(idev,30)k
   8  write(idev,90)j,s,d
   90 format(//2x,18hnumber of events =,i8,5x,
     *10haverage = ,e11.4,5x,8hsigma = ,e11.4//)
      return
    9 write(idev,80)i
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

      subroutine zero(a,n)
c------------------------------------------------------------------------c
c     clear array a(n)                                                   c
c------------------------------------------------------------------------c
      real a(n)

      do 10 i=1,n
         a(i) = 0.0
 10   continue

      return
      end
c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine izero(ia,n)
c------------------------------------------------------------------------c
c     clear integer array ia(n)                                          c
c------------------------------------------------------------------------c
      integer ia(n)

      do 10 i=1,n
         ia(i) = 0
 10   continue

      return
      end

