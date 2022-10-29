      program main
c----------------------------------------------------------------------------c
c     lattice calculation in quantum mechanics.                              c
c----------------------------------------------------------------------------c
c     this version applies cooling to MC configurations.                     c
c----------------------------------------------------------------------------c
c     version :            2.0                                               c
c     creation date :   07-23-96                                             c
c     last modified :   07-23-96                                             c
c----------------------------------------------------------------------------c
c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
c----------------------------------------------------------------------------c
c     lattice x(i=0,n), peridic b.c. x(0)=x(n)                               c
c----------------------------------------------------------------------------c
c     nmc    number of mc sweeps                                             c
c     nc     number of correlator measurements in a single configuration     c
c     neq    number of equlibration sweeps                                   c
c     kp     number of sweeps between cooling                                c
c     ncool  number of cooling sweeps in a single configuration              c
c     kp2    number of sweeps between writeout of complete configuration     c
c----------------------------------------------------------------------------c

      parameter(npar=10002,nmax=5000)
      real x(-1:npar),xs(-1:npar),xi(npar),xa(npar),z(npar)
      integer ix(100),iz(100)
      integer ipa(nmax)

c----------------------------------------------------------------------------c
c     correlators <x(#0)x(t)>                                                 c
c----------------------------------------------------------------------------c

      real xcor_sum(nmax), xcor2_sum(nmax)
      real xcor_av(nmax),  xcor_er(nmax)
      real xcool_sum(nmax),xcool2_sum(nmax)
      real xcool_av(nmax), xcool_er(nmax)

c----------------------------------------------------------------------------c
c     correlators <x^2(0)x^2(t)>                                             c
c----------------------------------------------------------------------------c

      real x2cor_sum(nmax),x2cor2_sum(nmax)
      real x2cor_av(nmax), x2cor_er(nmax)
      real x2sub_av(nmax), x2sub_er(nmax)
      real x2cool_sum(nmax),x2cool2_sum(nmax)
      real x2cool_av(nmax), x2cool_er(nmax)
      real x2cool_sub_av(nmax), x2cool_sub_er(nmax)

c----------------------------------------------------------------------------c
c     correlators <x^3(0)x^3(t)>                                             c
c----------------------------------------------------------------------------c

      real x3cor_sum(nmax),x3cor2_sum(nmax)
      real x3cor_av(nmax), x3cor_er(nmax)
      real x3cool_sum(nmax),x3cool2_sum(nmax)
      real x3cool_av(nmax), x3cool_er(nmax)

c----------------------------------------------------------------------------c
c     number of instantons, cooled action                                    c
c----------------------------------------------------------------------------c

      real nin_sum(0:nmax),nin2_sum(0:nmax)
      real nin_av(0:nmax), nin_er(0:nmax)    
      real scool_sum(0:nmax),scool2_sum(0:nmax)
      real scool_av(0:nmax), scool_er(0:nmax)  

      common /par/ f,a,delx
      common /seed/iseed

      open(unit=16,file='qm.dat',status='unknown')
      open(unit=17,file='config.dat',status='unknown')
      open(unit=18,file='trajectory.dat',status='unknown')
      open(unit=19,file='qmdist.dat',status='unknown')
      open(unit=20,file='coolconfig.dat',status='unknown')
      open(unit=21,file='cor.dat',status='unknown')
      open(unit=22,file='coolcor.dat',status='unknown')
      open(unit=23,file='nin.dat',status='unknown')
      open(unit=24,file='scool.dat',status='unknown')
      open(unit=25,file='sinst.dat',status='unknown')
      open(unit=26,file='cor2.dat',status='unknown')
      open(unit=27,file='coolcor2.dat',status='unknown')
      open(unit=28,file='cor3.dat',status='unknown')
      open(unit=29,file='coolcor3.dat',status='unknown')
      open(unit=30,file='zdist.dat',status='unknown')

c----------------------------------------------------------------------------c
c     input                                                                  c
c----------------------------------------------------------------------------c

      write(6,*) 'separation of wells f (f=1.4)'
      read(5,*)  f
      write(6,*) 'grid size n<10000 (n=100)'
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
      read(5,*)  kp2
      write(6,*) 'number of sweeps between cooling'
      read(5,*)  kp
      write(6,*) 'number of cooling sweeps (ncool<5000)'
      read(5,*)  ncool            
      tmax = n*a
      
c----------------------------------------------------------------------------c
c     echo input parameters                                                  c
c----------------------------------------------------------------------------c

      pi = 3.1415926
      s0 = 4.0/3.0*f**3
      de = 8*sqrt(2.0/pi)*f**2.5*exp(-s0)
      de2= de*(1.0-71.0/72.0/s0)
      write(16,*)   'lattice qm 1.1'
      write(16,*)   '--------------'
      write(16,101) f,n,a
      write(16,102) nmc,neq
      write(16,103) np,nc
      write(16,104) delx,icold,ncool
      write(16,105) s0,de,de*n*a
      write(16,106) s0,de2,de2*n*a
c     write(17,*)   '#',n,nmc/kp,n*a,f
c     write(20,*)   '#',n,nmc/kp,n*a,f
      
 101  format(1x,' f    = ',f8.2,' n   = ',i8,  ' a   = ',f8.4)
 102  format(1x,' nmc  = ',i8  ,' neq = ',i8)
 103  format(1x,' np   = ',i8  ,' nc  = ',i8)
 104  format(1x,' delx = ',f8.2,' icol= ',i8,  ' ncoo= ',i8)
 105  format(1x,' S_0  = ',f8.2,' dE  = ',f8.2,' dE*L= ',f8.2)
 106  format(1x,' S_0  = ',f8.2,' dE_2= ',f8.2,' dE2*L=',f8.2)
      
c----------------------------------------------------------------------------c
c     parameters for histograms                                              c
c----------------------------------------------------------------------------c
      
      nxhist = 50
      xhist_min = -1.5*f
      stxhist= 3.0*f/float(nxhist)
      nzhist = 40
      stzhist= 4.01/float(nzhist)

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
      call zero(xcool_sum,nmax)
      call zero(xcool2_sum,nmax)

      call zero(x2cor_sum,nmax)
      call zero(x2cor2_sum,nmax)
      call zero(x2cool_sum,nmax)
      call zero(x2cool2_sum,nmax)

      call zero(x3cor_sum,nmax)
      call zero(x3cor2_sum,nmax)
      call zero(x3cool_sum,nmax)
      call zero(x3cool2_sum,nmax)

      call zero(nin_sum,nmax+1)
      call zero(nin2_sum,nmax+1)
      call zero(scool_sum,nmax+1)
      call zero(scool2_sum,nmax+1)
      call izero(ix,nxhist)
      call izero(iz,nzhist)

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
      ncoolconf = 0
      ncoolcor  = 0

c----------------------------------------------------------------------------c
c     monte carlo sweeps                                                     c
c----------------------------------------------------------------------------c

      do 100 i=1,nmc
         
         nconf = nconf+1
              
         if(i .eq. neq)then
            nconf = 0
            ncor  = 0
            ncoolconf = 0
            ncoolcor  = 0
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
            call zero(xcool_sum,nmax)
            call zero(xcool2_sum,nmax)
            call zero(x2cor_sum,nmax)
            call zero(x2cor2_sum,nmax)
            call zero(x2cool_sum,nmax)
            call zero(x2cool2_sum,nmax)
            call zero(x3cor_sum,nmax)
            call zero(x3cor2_sum,nmax)
            call zero(x3cool_sum,nmax)
            call zero(x3cool2_sum,nmax)
            call zero(nin_sum,nmax)
            call zero(nin2_sum,nmax)
            call zero(scool_sum,nmax+1)
            call zero(scool2_sum,nmax+1)
            call izero(ix,nxhist)
            call izero(iz,nzhist)
         endif
         
c----------------------------------------------------------------------------c
c     one sweep through configuration                                        c
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
            xs(j)= x(j)
            ttot = ttot +a*t
            vtot = vtot +a*v
            tvtot= tvtot+a*tv
            stot = stot + s
 150     continue 
         if(i .le. 10000) write(18,*) i,stot,ttot,vtot
                  
         xs(n)   = x(n)
         xs(n+1) = x(n+1)
         xs(0)   = x(0)
         xs(-1)  = x(-1)
                         
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
            ipa(ic)= ip0
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
c     cooling and topological charge                                         c
c----------------------------------------------------------------------------c

         if (mod(i,kp) .eq. 0)then
         ncoolconf = ncoolconf + 1
         call inst(xs,n,ni,na,xi,xa,z)
         call act(n,xs,ss,ts,vs)
         nin = ni+na
         nin_sum(0)   = nin_sum(0)   + nin
         nin2_sum(0)  = nin2_sum(0)  + nin**2
         scool_sum(0) = scool_sum(0) + ss
         scool2_sum(0)= scool2_sum(0)+ ss**2

         do 300 icool=1,ncool
  
            call cool(xs,npar,n,1)
            call inst(xs,n,ni,na,xi,xa,z)
            call act(n,xs,ss,ts,vs)
            nin = ni+na
            nin_sum(icool)   = nin_sum(icool)   + nin
            nin2_sum(icool)  = nin2_sum(icool)  + nin**2
            scool_sum(icool) = scool_sum(icool) + ss
            scool2_sum(icool)= scool2_sum(icool)+ ss**2

 300     continue      

c----------------------------------------------------------------------------c
c     cooled configuration: instanton distribution                           c
c----------------------------------------------------------------------------c

         do 310 ii=1,nin-1,2
            if(ii .eq. 1) then
               zm = z(nin)-tmax
            else
               zm = z(ii-1)
            endif
            z0 = z(ii)
            zp = z(ii+1)
            zia= min(zp-z0,z0-zm)
            call lens(zia,0.0,stzhist,nzhist,iz)
 310     continue

c----------------------------------------------------------------------------c
c     cooled correlator                                                      c
c----------------------------------------------------------------------------c

         do 350 ic=1,nc
            ncoolcor = ncoolcor + 1 
            ip0 = ipa(ic)
            x0  = xs(ip0) 
            do 380 ip=1,np

               x1 = xs(ip0+ip)
               xcor = x0*x1
               x2cor= xcor**2
               x3cor= xcor**3

               xcool_sum(ip)  = xcool_sum(ip)  + xcor
               xcool2_sum(ip) = xcool2_sum(ip) + xcor**2
               xcool_sum(ip)  = xcool_sum(ip)  + xcor
               xcool2_sum(ip) = xcool2_sum(ip) + xcor**2
               x2cool_sum(ip) = x2cool_sum(ip) + x2cor
               x2cool2_sum(ip)= x2cool2_sum(ip)+ x2cor**2
               x3cool_sum(ip) = x3cool_sum(ip) + x3cor
               x3cool2_sum(ip)= x3cool2_sum(ip)+ x3cor**2

 380        continue
 350     continue
         endif

c----------------------------------------------------------------------------c
c     write configuration                                                    c
c----------------------------------------------------------------------------c
   
         if (mod(i,kp2) .eq. 0) then
            write(6,*)
            write(6,*) 'configuration   ',i
            write(6,*) 'acceptance rate ',float(nacc)/float(nhit)
            write(6,*) 'action (T,V)    ',stot,ttot,vtot         
c           write(17,*) '# configuration',i
c           write(20,*) '# configuration',i
            do 155 k=1,n
               write(17,*) k*a,x(k)
               write(20,*) k*a,xs(k)
 155        continue
         endif

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
 
c----------------------------------------------------------------------------c
c     correlators                                                            c
c----------------------------------------------------------------------------c

      do 400 ip=1,np
      call disp(ncor,xcor_sum(ip),xcor2_sum(ip),
     &          xcor_av(ip),xcor_er(ip))
      call disp(ncor,x2cor_sum(ip),x2cor2_sum(ip),
     &          x2cor_av(ip),x2cor_er(ip))
      call disp(ncor,x3cor_sum(ip),x3cor2_sum(ip),
     &          x3cor_av(ip),x3cor_er(ip))
      call disp(ncoolcor,xcool_sum(ip),xcool2_sum(ip),
     &          xcool_av(ip),xcool_er(ip))
      call disp(ncoolcor,x2cool_sum(ip),x2cool2_sum(ip),
     &          x2cool_av(ip),x2cool_er(ip))
      call disp(ncoolcor,x3cool_sum(ip),x3cool2_sum(ip),
     &          x3cool_av(ip),x3cool_er(ip))
 400  continue
 
c----------------------------------------------------------------------------c
c     instanton density, cooled action                                       c
c----------------------------------------------------------------------------c

      do 410 ic=0,ncool
      call disp(ncoolconf,nin_sum(ic),nin2_sum(ic),
     &          nin_av(ic),nin_er(ic))
      call disp(ncoolconf,scool_sum(ic),scool2_sum(ic),
     &          scool_av(ic),scool_er(ic))
 410  continue

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
      write(16,*) 'nconf = ',nconf
      write(16,*) 'ncoolc= ',ncoolconf
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
c     correlators etc                                                        c
c----------------------------------------------------------------------------c
      
      write(16,*) 
      write(16,*) ' <x(0)x(t)> correlation function'
      write(21,*) '#<x(0)x(t)> correlation function'
      do 500 ip=1,np
         dx = dl(xcor_av(ip),xcor_av(ip+1),a)
         dxe=dle(xcor_av(ip),xcor_av(ip+1),
     &             xcor_er(ip),xcor_er(ip+1),a)
         write(16,555) ip*a,xcor_av(ip),xcor_er(ip),dx,dxe
         write(21,555) ip*a,xcor_av(ip),xcor_er(ip),dx,dxe
 500  continue        

      write(16,*) 
      write(16,*) ' <x(0)x(t)> cooled correlation function'
      write(22,*) '#<x(0)x(t)> cooled correlation function'
      do 501 ip=1,np
         dx = dl(xcool_av(ip),xcool_av(ip+1),a)
         dxe=dle(xcool_av(ip),xcool_av(ip+1),
     &             xcool_er(ip),xcool_er(ip+1),a)
         write(16,555) ip*a,xcool_av(ip),xcool_er(ip),dx,dxe
         write(22,555) ip*a,xcool_av(ip),xcool_er(ip),dx,dxe
 501  continue        

c----------------------------------------------------------------------------c
c     <x^2(0)x^2(t) correlator requires subtraction                          c
c----------------------------------------------------------------------------c
 
      xx_sub = x2cor_av(np)
      xx_er  = x2cor_er(np)
      do 502 ip=1,np
         x2sub_av(ip) = x2cor_av(ip)-xx_sub
         x2sub_er(ip) = sqrt(x2cor_er(ip)**2+xx_er**2)
 502  continue

      write(16,*) 
      write(16,*) ' <x^2(0)x^2(t)> correlation function'
      write(26,*) '#<x^2(0)x^2(t)> correlation function'
      do 503 ip=1,np
         dx = dl(x2sub_av(ip),x2sub_av(ip+1),a)
         dxe=dle(x2sub_av(ip),x2sub_av(ip+1),
     &             x2sub_er(ip),x2sub_er(ip+1),a)
         write(16,555) ip*a,x2cor_av(ip),x2cor_er(ip),dx,dxe
         write(26,555) ip*a,x2cor_av(ip),x2cor_er(ip),dx,dxe
 503  continue        

      xx_sub = x2cool_av(np)
      xx_er  = x2cool_er(np)
      do 504 ip=1,np
         x2cool_sub_av(ip) = x2cool_av(ip)-xx_sub
         x2cool_sub_er(ip) = sqrt(x2cool_er(ip)**2+xx_er**2)
 504  continue

      write(16,*) 
      write(16,*) ' <x^2(0)x^2(t)> cooled correlation function'
      write(27,*) '#<x^2(0)x^2(t)> cooled correlation function'
      do 505 ip=1,np
         dx = dl(x2cool_sub_av(ip),x2cool_sub_av(ip+1),a)
         dxe=dle(x2cool_sub_av(ip),x2cool_sub_av(ip+1),
     &             x2cool_sub_er(ip),x2cool_sub_er(ip+1),a)
         write(16,555) ip*a,x2cool_av(ip),x2cool_er(ip),dx,dxe
         write(27,555) ip*a,x2cool_av(ip),x2cool_er(ip),dx,dxe
 505  continue        

c----------------------------------------------------------------------------c
c     x^3(0)x^3(t) correlator                                                c
c----------------------------------------------------------------------------c

      write(16,*) 
      write(16,*) ' <x^3(0)x^3(t)> correlation function'
      write(28,*) '#<x^3(0)x^3(t)> correlation function'
      do 506 ip=1,np
         dx = dl(x3cor_av(ip),x3cor_av(ip+1),a)
         dxe=dle(x3cor_av(ip),x3cor_av(ip+1),
     &             x3cor_er(ip),x3cor_er(ip+1),a)
         write(16,555) ip*a,x3cor_av(ip),x3cor_er(ip),dx,dxe
         write(28,555) ip*a,x3cor_av(ip),x3cor_er(ip),dx,dxe
 506  continue        

      write(16,*) 
      write(16,*) ' <x^3(0)x^3(t)> cooled correlation function'
      write(29,*) '#<x^3(0)x^3(t)> cooled correlation function'
      do 507 ip=1,np
         dx = dl(x3cool_av(ip),x3cool_av(ip+1),a)
         dxe=dle(x3cool_av(ip),x3cool_av(ip+1),
     &             x3cool_er(ip),x3cool_er(ip+1),a)
         write(16,555) ip*a,x3cool_av(ip),x3cool_er(ip),dx,dxe
         write(29,555) ip*a,x3cool_av(ip),x3cool_er(ip),dx,dxe
 507  continue        

c----------------------------------------------------------------------------c
c     instanton density                                                      c
c----------------------------------------------------------------------------c

      write(16,*) 
      write(16,*) ' number of instantons'
      write(23,*) '#number of instantons'
      do 510 ic=0,ncool
         write(16,556) ic,nin_av(ic),nin_er(ic),de*tmax,de2*tmax
         write(23,556) ic,nin_av(ic),nin_er(ic),de*tmax,de2*tmax
 510  continue        

      write(16,*) 
      write(16,*) ' action vs cooling sweeps'
      write(24,*) '#action vs cooling sweeps'
      do 515 ic=0,ncool
         sin = nin_av(ic)*s0
         write(16,443) ic,scool_av(ic),scool_er(ic),sin
         write(24,443) ic,scool_av(ic),scool_er(ic),sin
 515  continue                            
 
      write(16,*) 
      write(16,*) ' action per instanton, S_0 = ',4.0/3.0*f**3
      write(25,*) '#action per instanton, S_0 = ',4.0/3.0*f**3
      do 520 ic=0,ncool          
         si_av= scool_av(ic)/nin_av(ic)                    
         del2 =(scool_er(ic)/scool_av(ic))**2+(nin_er(ic)/nin_av(ic))**2
         si_er= si_av*sqrt(del2)
         write(16,443) ic,si_av,si_er,s0
         write(25,443) ic,si_av,si_er,s0
 520  continue        


c----------------------------------------------------------------------------c
c     histograms                                                             c
c----------------------------------------------------------------------------c
         
      write(16,*)
      write(16,*) ' x distribution '
      call lev(xhist_min,stxhist,nxhist,2,ix,16)

      write(16,*)
      write(16,*) ' z distribution '
      call lev(0.0,stzhist,nzhist,2,iz,16)
      do 610 i=1,nzhist 
         xx = (i+0.5)*stzhist
         write(30,*) xx,iz(i)
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
 332  format(1x,i4,2(f12.5))
 333  format(1x,3(f12.5))
 443  format(1x,i4,3(f12.5))
 444  format(1x,4(f12.5))
 555  format(1x,5(f12.5))
 556  format(1x,i4,4(f12.5))
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

      subroutine act(n,x,stot,ttot,vtot)
c-------------------------------------------------------------------------c
c     discretized action for configuration x(n)                           c
c-------------------------------------------------------------------------c
      real x(-1:n)
      
      common /par/ f,a,delx

      stot = 0.0
      ttot = 0.0
      vtot = 0.0  
      do 50 j=1,n 
         xp = (x(j+1)-x(j))/a
         t  = 1.0/4.0*xp**2
         v  = (x(j)**2-f**2)**2
         s  = a*(t+v)
         ttot = ttot +a*t
         vtot = vtot +a*v
         stot = stot + s
  50  continue   

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---
 
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

      function dl(xcor1,xcor2,a)
c----------------------------------------------------------------------------c
c     log derivative                                                         c
c----------------------------------------------------------------------------c
      dl = (xcor1-xcor2)/(xcor1*a)
      end

      function dle(xcor1,xcor2,xcor1e,xcor2e,a)
c----------------------------------------------------------------------------c
c     log derivative, error                                                  c
c----------------------------------------------------------------------------c
      dle2 = (xcor2e/xcor1)**2+(xcor1e*xcor2/xcor1**2)**2
      dle  = sqrt(dle2)
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine smooth(x,xs,npar,n,ns)
c----------------------------------------------------------------------------c
c     naive data smoothing with sliding window                               c
c----------------------------------------------------------------------------c
     
      real x(-1:npar),xs(-1:npar)
      
      do 10 i=1,n
         xsum = 0.0
         do 20 j=i-ns,i+ns
            jp = j
            if(j .lt. 1) jp=n+j
            if(j .gt. n) jp=j-n
            xsum = xsum + x(jp)
  20     continue
         xav = xsum/float(2*ns+1)
         xs(i)= xav
  10  continue

      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

      subroutine cool(xs,npar,n,ncool)
c----------------------------------------------------------------------------c
c     local cooling algorithm                                                c
c----------------------------------------------------------------------------c
    
      real xs(-1:npar)
      common /par/  f,a,delx
      common /seed/ iseed
      
      nhit = 10
      delxp= 0.1*delx

      do 5  k=1,ncool

      do 10 i=1,n

         xpm = (xs(i)-xs(i-1))/a
         xpp = (xs(i+1)-xs(i))/a
         t = 1.0/4.0*(xpm**2+xpp**2)
         v = (xs(i)**2-f**2)**2
         sold = a*(t+v)

         do 20 j=1,nhit
            
            xnew = xs(i) + delxp*(2.0*ran2(iseed)-1.0)
             
            xpm = (xnew-xs(i-1))/a
            xpp = (xs(i+1)-xnew)/a
            t = 1.0/4.0*(xpm**2+xpp**2)
            v = (xnew**2-f**2)**2
            snew = a*(t+v)

            if (snew .lt. sold) xs(i)=xnew

  20     continue

  10  continue

   5  continue

      return
      end

c---------------------------------------------------------------------+----
c---------------------------------------------------------------------+----

      subroutine inst(x,n,ni,na,xi,xa,z)
c------------------------------------------------------------------------c
c     return number and location of (anti) instantons                    c
c------------------------------------------------------------------------c
      real x(n),xi(n),xa(n),z(n)

      common /par/ f,a,delx

      ni = 0
      na = 0
      nin= 0
 
      ix = int(sign(1.0,x(1)))
 
      do 10 i=2,n
         tau = a*i
         ixp = int(sign(1.0,x(i)))
         if(ixp .gt. ix) then
            ni = ni+1
            nin= nin+1
            xi(ni) = tau
            z(nin) = tau
         else if(ixp .lt. ix) then
            na = na+1
            nin= nin+1
            xa(na) = tau
            z(nin) = tau
         endif
         ix = ixp
  10  continue
      
      return
      end

c---------------------------------------------------------------------+---
c---------------------------------------------------------------------+---

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

