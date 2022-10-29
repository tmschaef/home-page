      program main
c----------------------------------------------------------------------------c
c     streamline configurations in quantum mechanical double well potential. c
c----------------------------------------------------------------------------c
c     version :            1.0                                               c
c     creation date :   07-23-96                                             c
c     last modified :   07-23-96                                             c
c----------------------------------------------------------------------------c
c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
c----------------------------------------------------------------------------c

      parameter(npar=1002)
      implicit real*8 (a-h,o-z) 
      real*8 x(-npar:npar),xn(-npar:npar),sd(-npar:npar)
      real*8 xs(-npar:npar),for(-npar:npar),fors(-npar:npar)
      common /par/ f,dtau
      common /seed/ iseed

      open(unit=16,file='stream.dat',status='unknown')
      open(unit=17,file='sint_stream.dat',status='unknown')
      write(6,*) 'separation of wells f (f=1.4)'
      read(5,*)  f
      write(6,*) 'initial pair separation r0 (r0=1.8)'
      read(5,*)  r0
      write(6,*) 'streamline time step dxi (dxi=0.001)'
      read(5,*)  dxi
      write(6,*) 'grid size n<1000 (n=100)'
      read(5,*)  n
      write(6,*) 'grid spacing dtau (dtau=0.05)'
      read(5,*)  dtau
      iseed = -1234

c----------------------------------------------------------------------------c
c     ns: smoothing window, ncool: cooling sweeps                            c
c----------------------------------------------------------------------------c

      n1 = n+2
      ns = 2
      ncool = 2
      
c----------------------------------------------------------------------------c
c     initialize single instanton, sum, product                              c
c----------------------------------------------------------------------------c
  
      taui =-r0/2.0
      taua = r0/2.0 

      do 10 i=-n1,n1
         tau  = i*dtau
         x(i) = f*(tanh(2*f*(tau-taui)) - tanh(2*f*(tau-taua)) - 1.0 )
c        x(i) = f*tanh(2*f*(tau-taui))*tanh(2*f*(tau-taua))
c        x(i) = f*tanh(2*f*(tau-taui))
c        x(i) =-f+0.1*exp(-tau**2)
         xn(i)= x(i)      
         nin = 2
  10  continue

c----------------------------------------------------------------------------c
c     calculate action                                                       c
c----------------------------------------------------------------------------c

      s0   = 4.0/3.0*f**3
      stot = 0.0
      nconf= 1

      do 20 i=-n,n
         tau  = i*dtau
         xdot = (x(i+1)-x(i-1))/(2.0*dtau)
         v    = (x(i)**2-f**2)**2
         sd(i)= xdot**2/4. + v
         stot = stot + sd(i)*dtau
c        xddot = (x(i+1)-2*x(i)+x(i-1))/dtau**2
         xddot = (-x(i+2)+16*x(i+1)-30*x(i)+16*x(i-1)-x(i-2))
     1          /(12*dtau**2)     
         for(i)= -xddot/2. + 4*x(i)*(x(i)**2-f**2)
  20  continue

      call smooth(for,fors,npar,n,ns)

      sint = stot - s0*nin

      write(16,999) j,nconf,stot/s0,sint/s0  
      write(16,998)
      do 30 i=-n,n 
         tau = i*dtau
         write(16,666) tau,x(i),x(i),sd(i),for(i),fors(i)
  30  continue

      write(6,*) 's = ',stot/s0

c----------------------------------------------------------------------------c
c     streamline evolution                                                   c
c----------------------------------------------------------------------------c

      sinit = stot
      slast = stot
      sstep = stot
      nmax  = 1000000

      do 100 j=1,nmax

c----------------------------------------------------------------------------c
c     streamline step                                                        c
c----------------------------------------------------------------------------c

         do 50 i=-n,n
c           xddot = (x(i+1)-2*x(i)+x(i-1))/dtau**2
            xddot = (-x(i+2)+16*x(i+1)-30*x(i)+16*x(i-1)-x(i-2))
     1             /(12*dtau**2)     
            for(i)= -xddot/2. + 4*x(i)*(x(i)**2-f**2)
  50     continue  

         call smooth(for,fors,npar,n,ns)

         do 55 i=-n,n
            xn(i) = x(i) - for(i)*dxi
  55     continue
         xn(-n-2) = xn(-n)
         xn(-n-1) = xn(-n)
         xn(n+1)  = xn(n)
         xn(n+2)  = xn(n)

         call smooth(xn,xs,npar,n,ns)

         stot = 0.0
         do 60 i=-n,n
            tau  = i*dtau
            x(i) = xn(i)
            xdot = (xn(i+1)-xn(i-1))/(2.0*dtau)
            v    = (xn(i)**2-f**2)**2
            sd(i)= xdot**2/4. + v
            stot = stot + sd(i)*dtau
  60     continue
         x(-n-2) = xn(-n-2)
         x(-n-1) = xn(-n-1)
         x(n+1)  = xn(n+1)
         x(n+2)  = xn(n+2)

         ds  = sold - stot
         sold= stot
         sint= stot - s0*nin
         write(6,*) 's,ds,sint',stot/s0,ds/s0,sint/s0

c----------------------------------------------------------------------------c
c     write configuration                                                    c
c----------------------------------------------------------------------------c

         if (abs(slast-stot)/sinit .gt. 0.1) then
            slast = stot
            nconf = nconf + 1
            write(16,999) j,nconf,stot/s0,sint/s0
            write(16,998)  
            do 70 i=-n,n 
               tau = i*dtau
               write(16,666) tau,xn(i),xs(i),sd(i),for(i),fors(i)         
  70        continue
            call dist(n,dtau,xs,sd,dis,dis2)
            write(17,333) dis,dis2,sint/s0 
         endif
         if (stot/sinit .lt. 0.05) then
            nconf = nconf + 1
            write(16,999) j,nconf,stot/s0 
            do 80 i=-n,n 
               tau = i*dtau
               write(16,666) tau,xn(i),xs(i),sd(i),for(i),fors(i)
  80        continue
            call dist(n,dtau,xs,sd,dis,dis2)
            write(17,333) dis,dis2,sint/s0 
            stop
          endif

 100  continue

 222  format(1x,2(f12.5))
 333  format(1x,3(f12.5))
 444  format(1x,4(f12.5))
 555  format(1x,5(f12.5))
 666  format(1x,6(f12.5))
 998  format(5x,'# tau',9x,'xnew',8x,'xcool',7x,'s_dens',6x,'force',
     1          7x,'f_smooth')
 999  format(1x,'# j =',i5,' nconf =',i3,' stot =',f12.5,
     1          ' sint =',f12.5)
      end

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

      subroutine cool(x,xs,npar,n,ncool)
c----------------------------------------------------------------------------c
c     local cooling algorithm                                                c
c----------------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      real*8 x(-npar:npar),xs(-npar:npar)
      common /par/ f,dtau
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

      subroutine dist(n,dtau,xs,sd,dis,dis2)
c------------------------------------------------------------------------c
c     estimate distance                                                  c
c------------------------------------------------------------------------c
      parameter(npar=1002)
      implicit real*8 (a-h,o-z)
      real*8 xs(-npar:npar),sd(-npar:npar)
      imax = n 
      smax = sd(n)
      icr  = 0 
      do 10 i=n,0,-1
         ss = sd(i) 
         if(ss .gt. smax) then
            smax = ss
            imax = i
         endif
         if(xs(i) .lt. 0 .and. xs(i-1) .gt. 0) icr=i
 10   continue
      dis = 2*imax*dtau
      dis2= 2*icr*dtau
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
   
   
