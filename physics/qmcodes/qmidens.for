      program main
c----------------------------------------------------------------------------c
c     lattice calculation in quantum mechanics.                              c
c----------------------------------------------------------------------------c
c     calculate instanton density by adiabatically switching from gaussian   c
c     approximation to full potential. Have                                  c
c      d(inst)=d(gaussian)*exp(-S(non-gaussian))                             c
c     where S(non-gaussian)=\int d\alpha <S-S(gaussian)>_\alpha. perform     c
c     reference calculation of fluctuations around trivial vacuum.           c
c----------------------------------------------------------------------------c
c     instanton is placed at beta/2. anti-symmetric boundary conditions are  c
c     used. position is fixed during update by requiring x(beta/2)=0.        c
c----------------------------------------------------------------------------c
c     version :            1.0                                               c
c     creation date :   07-23-96                                             c
c     last modified :   07-23-96                                             c
c----------------------------------------------------------------------------c
c     action m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                       c
c----------------------------------------------------------------------------c
      parameter(npar=10002,nmax=100)
      real x(-1:npar),x0(-1:npar),w(-1:npar),vi(-1:npar)
      real va_av(0:nmax),va_err(0:nmax)

      common /par/ f,a
      common /seed/iseed

      open(unit=16,file='qmidens.dat',status='unknown')
      open(unit=17,file='idens_conf.dat',status='unknown')
      open(unit=18,file='vac_conf.dat',status='unknown')

      write(6,*) 'separation of wells f (f=1.4)'
      read(5,*)  f
      write(6,*) 'grid size n<1000, n even (n=100)'
      read(5,*)  n
      write(6,*) 'grid spacing a (dtau=0.05)'
      read(5,*)  a
      write(6,*) 'equilibration sweeps'
      read(5,*)  neq
      write(6,*) 'monte carlo sweeps'
      read(5,*)  nmc 
      write(6,*) 'write every kp configuration'
      read(5,*)  kp   
      write(6,*) 'update x (delx)'
      read(5,*)  delx 
      write(6,*) 'number of steps in switching process'
      read(5,*)  nalpha   
      
      dalpha = 1.0/float(nalpha)
      beta   = n*a
      tau0   = beta/2.0
      n0     = n/2
      n0p    = n0+1
      n0m    = n0-1
      sng    = 0.0
      svacng = 0.0
      pi = 3.1415926
      
c----------------------------------------------------------------------------c
c     echo input parameters                                                  c
c----------------------------------------------------------------------------c

      write(16,*)   'lattice qm idens 1.0'
      write(16,*)   '--------------------'
      write(16,101) f,n,a
      write(16,102) nmc,neq
      write(16,104) delx
      write(16,105) nalpha
      write(16,*) 
      write(16,*)   'fluctuations around instanton path'
      write(16,*)   '----------------------------------'
      
 101  format(1x,' f    = ',f8.2,' n   = ',i5,' a   = ',f5.4)
 102  format(1x,' nmc  = ',i8  ,' neq = ',i5)
 104  format(1x,' delx = ',f8.2) 
 105  format(1x,' nalp = ',i5)

c----------------------------------------------------------------------------c
c     clear summation arrays                                                 c
c----------------------------------------------------------------------------c

      stot_sum   = 0.0
      stot2_sum  = 0.0
      vav_sum    = 0.0
      vav2_sum   = 0.0 
      valpha_sum = 0.0
      valpha2_sum= 0.0

c----------------------------------------------------------------------------c
c     initialize instanton and gaussian potential                            c
c----------------------------------------------------------------------------c

      iseed =-1234

      do 10 i=1,n
         tau  = i*a
         x(i) = f*tanh(2.0*f*(tau-tau0))
         x0(i)= x(i)
         w(i) =-4.0*(f**2-3.0*x0(i)**2)
         vi(i)= (x0(i)**2-f**2)**2
         write(17,222) tau,x(i)
  10  continue       
  
c----------------------------------------------------------------------------c
c     anti-periodic boundary conditions                                      c
c----------------------------------------------------------------------------c

      x(n)   =-x(1)
      x(n+1) =-x(2)
      x(0)   =-x(n-1)
      x(-1)  =-x(n-2)
      
c----------------------------------------------------------------------------c
c     initial action                                                         c
c----------------------------------------------------------------------------c

      stot = 0.0
           
      do 50 i=1,n 
         xp = (x(i+1)-x(i))/a
         t  = 1.0/4.0*xp**2    
         v  = (x(i)**2-f**2)**2
         s  = a*(t+v)
         stot = stot + s
  50  continue
   
      s0   = 4.0/3.0*f**3
      dens = 8*sqrt(2.0/pi)*f**2.5*exp(-s0)
      f0   = dens

      write(16,*) 
      write(16,801) f0,stot,s0
 801  format(' f = ',f8.5,' s0(num) = ',f8.5,' s(ana) = ',f8.5)
 
c----------------------------------------------------------------------------c
c     loop over coupling constant (up/down)                                  c
c----------------------------------------------------------------------------c

      do 500 ialpha=0,2*nalpha      
      
      if(ialpha .le. nalpha)then
         alpha = ialpha*dalpha
      else
         alpha = 2.0-ialpha*dalpha
      endif          
                
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
            nconf = 0
            stot_sum   = 0.0
            stot2_sum  = 0.0
            vav_sum    = 0.0
            vav2_sum   = 0.0 
            valpha_sum = 0.0
            valpha2_sum= 0.0
         endif
         
c----------------------------------------------------------------------------c
c     one sweep thorough configuration                                       c
c----------------------------------------------------------------------------c
 
         do 200 j=1,n-1
            
            nhit = nhit+1  

c----------------------------------------------------------------------------c
c     fix instanton center                                                   c
c----------------------------------------------------------------------------c

            if(j .eq. n0)goto 200

c----------------------------------------------------------------------------c
c     old action                                                             c
c----------------------------------------------------------------------------c
 
            xpm = (x(j)-x(j-1))/a
            xpp = (x(j+1)-x(j))/a
            t  = 1.0/4.0*(xpm**2+xpp**2)
            v0 = 1.0/2.0*w(j)*(x(j)-x0(j))**2+vi(j)
            v1 = (x(j)**2-f**2)**2
            v  = alpha*(v1-v0) + v0
            sold = a*(t+v)

c----------------------------------------------------------------------------c
c     jakobian of constraint                                                 c
c----------------------------------------------------------------------------c

            if(j .eq. n0m .or. j .eq. n0p)then
               vel  = (x(n0p)-x(n0m))/(2.0*a)
               sjak =-log(abs(vel))
               sold = sold + sjak
            endif

c----------------------------------------------------------------------------c
c     MC hit                                                                 c
c----------------------------------------------------------------------------c
 
            xnew = x(j) + delx*(2.0*ran2(iseed)-1.0)

c----------------------------------------------------------------------------c
c     new action                                                             c
c----------------------------------------------------------------------------c
             
            xpm = (xnew-x(j-1))/a
            xpp = (x(j+1)-xnew)/a
            t  = 1.0/4.0*(xpm**2+xpp**2)
            v0 = 1.0/2.0*w(j)*(xnew-x0(j))**2+vi(j)
            v1 = (xnew**2-f**2)**2
            v  = alpha*(v1-v0) + v0
            snew = a*(t+v)

c----------------------------------------------------------------------------c
c     jakobian                                                               c
c----------------------------------------------------------------------------c

            if(j .eq. n0m)then
               vel  = (x(n0p)-xnew)/(2.0*a)
               sjak =-log(abs(vel))
               snew = snew + sjak
            else if(j .eq. n0p)then
               vel  = (xnew-x(n0m))/(2.0*a)
               sjak =-log(abs(vel))
               snew = snew + sjak
            endif
                       
c----------------------------------------------------------------------------c
c     accept/reject                                                          c
c----------------------------------------------------------------------------c

            dels = snew-sold      
            p    = ran2(iseed)
            dels = min(dels,70.0)
            dels = max(dels,-70.0)
            
            if (exp(-dels) .gt. p) then
               x(j) = xnew
               nacc = nacc + 1
            endif
            
 200     continue        
 
         x(n)   =-x(1)
         x(n+1) =-x(2)
         x(0)   =-x(n-1)
         x(-1)  =-x(n-2)
         
c----------------------------------------------------------------------------c
c     calculate action etc.                                                  c
c----------------------------------------------------------------------------c

         stot = 0.0
         ttot = 0.0
         vtot = 0.0 
         ptot = 0.0
          
         do 150 j=1,n 
            xp = (x(j+1)-x(j))/a
            t  = 1.0/4.0*xp**2       
            v0 = 1.0/2.0*w(j)*(x(j)-x0(j))**2+vi(j)
            v1 = (x(j)**2-f**2)**2
            v  = alpha*(v1-v0) + v0
            s  = a*(t+v)
            ttot = ttot +a*t
            vtot = vtot +a*v
            stot = stot + s    
            ptot = ptot +a*(v1-v0)
 150     continue
                                           
         if (mod(i,kp) .eq. 0) then
            write(6,*)
            write(6,*) 'configuration   ',i  
            write(6,*) 'coupling        ',alpha
            write(6,*) 'acceptance rate ',float(nacc)/float(nhit)
            write(6,*) 'action (T,V)    ',stot,ttot,vtot
         endif

c----------------------------------------------------------------------------c
c     include in sample                                                      c
c----------------------------------------------------------------------------c
         
         stot_sum = stot_sum + stot
         stot2_sum= stot2_sum+ stot**2
         vav_sum  = vav_sum  + vtot/beta
         vav2_sum = vav2_sum + vtot**2/beta
         valpha_sum = valpha_sum + ptot/beta
         valpha2_sum= valpha2_sum+ ptot**2/beta

c----------------------------------------------------------------------------c
c     next configuration                                                     c
c----------------------------------------------------------------------------c

 100  continue
                   
c----------------------------------------------------------------------------c
c     averages                                                               c
c----------------------------------------------------------------------------c

      call disp(nconf,stot_sum,stot2_sum,stot_av,stot_err)
      call disp(nconf,vav_sum,vav2_sum,v_av,v_err)
      call disp(nconf,valpha_sum,valpha2_sum,valpha_av,valpha_err)      
                             
      va_av(ialpha) = valpha_av
      va_err(ialpha)= valpha_err
      
      if(mod(ialpha,2*nalpha) .eq. 0) then
         da = dalpha/4.0
      else 
         da = dalpha/2.0
      endif
      dsng = da*valpha_av
      sng  = sng + dsng
      
c----------------------------------------------------------------------------c
c     output                                                                 c
c----------------------------------------------------------------------------c

      write(16,*)  
      write(16,900) alpha
      write(16,901) stot_av,stot_err
      write(16,905) v_av,v_err
      write(16,906) valpha_av,valpha_err
      write(16,907) sng,dsng
      write(16,*) 
      do 170 j=1,n
         write(17,222) j*a,x(j)
 170  continue

c----------------------------------------------------------------------------c
c     end of loop over coupling constants                                    c
c----------------------------------------------------------------------------c 
      
 500  continue

c----------------------------------------------------------------------------c
c     final estimate of integral over coupling                               c
c----------------------------------------------------------------------------c
                                  
      sup_sum = 0.0
      sup_err = 0.0
      sup_hal = 0.0
      sdw_sum = 0.0
      sdw_err = 0.0
      sdw_hal = 0.0
      
c----------------------------------------------------------------------------c
c     have sum=1/2(up+down) and up = 1/2*f0+f1+...+1/2*fn, down=...          c
c----------------------------------------------------------------------------c

      do 600 ia=0,nalpha   
         if(mod(ia,nalpha) .eq. 0)then
            da = dalpha/4.0
         else
            da = dalpha/2.0
         endif
         iap = ia+nalpha 
         sup_sum = sup_sum + da*va_av(ia)
         sup_err = sup_err + da*va_err(ia)**2
         sdw_sum = sdw_sum + da*va_av(iap)
         sdw_err = sdw_err + da*va_err(iap)**2
 600  continue

      do 620 ia=0,nalpha,2   
         if(mod(ia,nalpha) .eq. 0)then
            da = dalpha/2.0
         else
            da = dalpha
         endif
         iap = ia+nalpha 
         sup_hal = sup_hal + da*va_av(ia)
         sdw_hal = sdw_hal + da*va_av(iap)
 620  continue          
      
c----------------------------------------------------------------------------c
c     uncertainties                                                          c
c----------------------------------------------------------------------------c
      
      ds     = sup_sum + sdw_sum
      sng    = ds
      dens_ng= dens*exp(-sng)
      ds_err = sqrt(sup_err + sdw_err)
      ds_hal = sup_hal + sdw_hal
      ds_dif = abs(sup_sum - sdw_sum)
      ds_dis = abs(ds - ds_hal)/2.0
      ds_tot = sqrt(ds_err**2 + ds_dif**2 + ds_dis**2)
      dens_er= dens_ng*ds_tot 
                                                                               
c----------------------------------------------------------------------------c
c     output                                                                 c
c----------------------------------------------------------------------------c
          
      do 700 i=6,16,10
      write(i,*)
      write(i,911) sng,ds_tot
      write(i,912) ds_err,ds_dif,ds_dis
      write(i,*) 
      write(i,*) 'reference calculation: trivial vacuum'
      write(i,*) '-------------------------------------'
      write(i,*)
 700  continue 

c----------------------------------------------------------------------------c
c     repeat calculation for fluctuations around trivial vacuum              c
c----------------------------------------------------------------------------c
c----------------------------------------------------------------------------c
c     clear summation arrays                                                 c
c----------------------------------------------------------------------------c

      stot_sum   = 0.0
      stot2_sum  = 0.0
      vav_sum    = 0.0
      vav2_sum   = 0.0 
      valpha_sum = 0.0
      valpha2_sum= 0.0

c----------------------------------------------------------------------------c
c     initialize                                                             c
c----------------------------------------------------------------------------c

      do 1010 i=1,n
         tau  = i*a
         x(i) = f
         x0(i)= x(i)
         w(i) = 8.0*f**2
         write(18,222) tau,x(i)
 1010 continue       
  
c----------------------------------------------------------------------------c
c     periodic boundary conditions                                           c
c----------------------------------------------------------------------------c

      x(n)   = x(1)
      x(n+1) = x(2)
      x(0)   = x(n-1)
      x(-1)  = x(n-2)
 
c----------------------------------------------------------------------------c
c     loop over coupling constant (up/down)                                  c
c----------------------------------------------------------------------------c

      do 1500 ialpha=0,2*nalpha      
      
      if(ialpha .le. nalpha)then
         alpha = ialpha*dalpha
      else
         alpha = 2.0-ialpha*dalpha
      endif          
                
      nacc = 0
      nhit = 0    
      nconf= 0
      ncor = 0

c----------------------------------------------------------------------------c
c     monte carlo sweeps                                                     c
c----------------------------------------------------------------------------c

      do 1100 i=1,nmc
         
         nconf = nconf+1
              
         if(i .eq. neq)then
            nconf = 0
            stot_sum   = 0.0
            stot2_sum  = 0.0
            vav_sum    = 0.0
            vav2_sum   = 0.0 
            valpha_sum = 0.0
            valpha2_sum= 0.0
         endif
         
c----------------------------------------------------------------------------c
c     one sweep thorough configuration                                       c
c----------------------------------------------------------------------------c
 
         do 1200 j=1,n-1
            
            nhit = nhit+1  

c----------------------------------------------------------------------------c
c     old action                                                             c
c----------------------------------------------------------------------------c
 
            xpm = (x(j)-x(j-1))/a
            xpp = (x(j+1)-x(j))/a
            t  = 1.0/4.0*(xpm**2+xpp**2)
            v0 = 1.0/2.0*w(j)*(x(j)-x0(j))**2
            v1 = (x(j)**2-f**2)**2
            v  = alpha*(v1-v0) + v0
            sold = a*(t+v)

c----------------------------------------------------------------------------c
c     MC hit                                                                 c
c----------------------------------------------------------------------------c
 
            xnew = x(j) + delx*(2.0*ran2(iseed)-1.0)

c----------------------------------------------------------------------------c
c     new action                                                             c
c----------------------------------------------------------------------------c
             
            xpm = (xnew-x(j-1))/a
            xpp = (x(j+1)-xnew)/a
            t  = 1.0/4.0*(xpm**2+xpp**2)
            v0 = 1.0/2.0*w(j)*(xnew-x0(j))**2
            v1 = (xnew**2-f**2)**2
            v  = alpha*(v1-v0) + v0
            snew = a*(t+v)
                       
c----------------------------------------------------------------------------c
c     accept/reject                                                          c
c----------------------------------------------------------------------------c

            dels = snew-sold      
            p    = ran2(iseed)
            dels = min(dels,70.0)
            dels = max(dels,-70.0)
            
            if (exp(-dels) .gt. p) then
               x(j) = xnew
               nacc = nacc + 1
            endif
            
 1200    continue        
 
         x(n)   = x(1)
         x(n+1) = x(2)
         x(0)   = x(n-1)
         x(-1)  = x(n-2)
         
c----------------------------------------------------------------------------c
c     calculate action etc.                                                  c
c----------------------------------------------------------------------------c

         stot = 0.0
         ttot = 0.0
         vtot = 0.0 
         ptot = 0.0
          
         do 1150 j=1,n 
            xp = (x(j+1)-x(j))/a
            t  = 1.0/4.0*xp**2       
            v0 = 1.0/2.0*w(j)*(x(j)-x0(j))**2
            v1 = (x(j)**2-f**2)**2
            v  = alpha*(v1-v0) + v0
            s  = a*(t+v)
            ttot = ttot +a*t
            vtot = vtot +a*v
            stot = stot + s    
            ptot = ptot +a*(v1-v0)
1150     continue
                                           
         if (mod(i,kp) .eq. 0) then
            write(6,*)
            write(6,*) 'configuration   ',i  
            write(6,*) 'coupling        ',alpha
            write(6,*) 'acceptance rate ',float(nacc)/float(nhit)
            write(6,*) 'action (T,V)    ',stot,ttot,vtot
         endif

c----------------------------------------------------------------------------c
c     include in sample                                                      c
c----------------------------------------------------------------------------c
         
         stot_sum = stot_sum + stot
         stot2_sum= stot2_sum+ stot**2
         vav_sum  = vav_sum  + vtot/beta
         vav2_sum = vav2_sum + vtot**2/beta
         valpha_sum = valpha_sum + ptot/beta
         valpha2_sum= valpha2_sum+ ptot**2/beta

c----------------------------------------------------------------------------c
c     next configuration                                                     c
c----------------------------------------------------------------------------c

 1100 continue
                   
c----------------------------------------------------------------------------c
c     averages                                                               c
c----------------------------------------------------------------------------c

      call disp(nconf,stot_sum,stot2_sum,stot_av,stot_err)
      call disp(nconf,vav_sum,vav2_sum,v_av,v_err)
      call disp(nconf,valpha_sum,valpha2_sum,valpha_av,valpha_err)      
                             
      va_av(ialpha) = valpha_av
      va_err(ialpha)= valpha_err
      
      if(mod(ialpha,2*nalpha) .eq. 0) then
         da = dalpha/4.0
      else 
         da = dalpha/2.0
      endif
      dsvacng = da*valpha_av
      svacng  = svacng + dsvacng
      
c----------------------------------------------------------------------------c
c     output                                                                 c
c----------------------------------------------------------------------------c

      write(16,*)  
      write(16,900) alpha
      write(16,901) stot_av,stot_err
      write(16,905) v_av,v_err
      write(16,906) valpha_av,valpha_err
      write(16,907) svacng,dsvacng
      write(16,*) 
      do 1170 j=1,n
         write(18,222) j*a,x(j)
1170  continue

c----------------------------------------------------------------------------c
c     end of loop over coupling constants                                    c
c----------------------------------------------------------------------------c 
      
 1500 continue

c----------------------------------------------------------------------------c
c     final estimate of integral over coupling                               c
c----------------------------------------------------------------------------c
                                  
      svacup_sum = 0.0
      svacup_err = 0.0
      svacup_hal = 0.0
      svacdw_sum = 0.0
      svacdw_err = 0.0
      svacdw_hal = 0.0
      
c----------------------------------------------------------------------------c
c     have sum=1/2(up+down) and up = 1/2*f0+f1+...+1/2*fn, down=...          c
c----------------------------------------------------------------------------c

      do 1600 ia=0,nalpha   
         if(mod(ia,nalpha) .eq. 0)then
            da = dalpha/4.0
         else
            da = dalpha/2.0
         endif
         iap = ia+nalpha 
         svacup_sum = svacup_sum + da*va_av(ia)
         svacup_err = svacup_err + da*va_err(ia)**2
         svacdw_sum = svacdw_sum + da*va_av(iap)
         svacdw_err = svacdw_err + da*va_err(iap)**2
1600  continue

      do 1620 ia=0,nalpha,2   
         if(mod(ia,nalpha) .eq. 0)then
            da = dalpha/2.0
         else
            da = dalpha
         endif
         iap = ia+nalpha 
         svacup_hal = svacup_hal + da*va_av(ia)
         svacdw_hal = svacdw_hal + da*va_av(iap)
1620  continue          
      
c----------------------------------------------------------------------------c
c     uncertainties                                                          c
c----------------------------------------------------------------------------c
      
      dsvac     = svacup_sum + svacdw_sum
      svacng    = dsvac
      fvac      = exp(-svacng)
      dsvac_err = sqrt(svacup_err + svacdw_err)
      dsvac_hal = svacup_hal + svacdw_hal
      dsvac_dif = abs(svacup_sum - svacdw_sum)
      dsvac_dis = abs(dsvac - dsvac_hal)/2.0
      dsvac_tot = sqrt(dsvac_err**2 + dsvac_dif**2 + dsvac_dis**2)
      fvac_er   = fvac*dsvac_tot 
                                                                               
c----------------------------------------------------------------------------c
c     output                                                                 c
c----------------------------------------------------------------------------c
          
      do 1700 i=6,16,10
      write(i,*)
      write(i,1911) svacng,dsvac_tot
      write(i,1912) dsvac_err,dsvac_dif,dsvac_dis
      write(i,*) 
      write(i,1914) fvac,fvac_er
1700  continue 

c----------------------------------------------------------------------------c
c     final answer                                                           c
c----------------------------------------------------------------------------c

      seff    = sng - svacng
      seff_er = sqrt(ds_tot**2+dsvac_tot**2)
      dens_ng = dens*exp(-seff)
      dens_er = dens_ng*seff_er

      do 1800 i=6,16,10
      write(i,*)
      write(i,915) seff,sng,svacng
      write(i,916) seff,seff_er
      write(i,*) 
      write(i,917) dens_ng,dens
      write(i,918) dens_ng,dens_er
1800  continue 

 900  format(1x,'alpha = ',f12.5)      
 901  format(1x,'stot  = ',f12.5,' +/- ',f12.5)
 905  format(1x,'v     = ',f12.5,' +/- ',f12.5)
 906  format(1x,'v_alp = ',f12.5,' +/- ',f12.5)
 907  format(1x,'S_ng  = ',f12.5,' dS= ',f12.5) 
 911  format(1x,'S_ng  = ',f12.5,' +/- ',f12.5)
 912  format(1x,'stat  = ',f12.5,' up/d',f12.5,' disc= ',f12.5)
 913  format(1x,'d_gs  = ',f12.5)
 914  format(1x,'d_ng  = ',f12.5,' +/- ',f12.5)
1911  format(1x,'Svac_n= ',f12.5,' +/- ',f12.5)
1912  format(1x,'stat  = ',f12.5,' up/d',f12.5,' disc= ',f12.5)
1914  format(1x,'Fvac_n= ',f12.5,' +/- ',f12.5)
 915  format(1x,'seff  = ',f12.5,' ins ',f12.5,' vac = ',f12.5)
 916  format(1x,'seff  = ',f12.5,' +/- ',f12.5)
 917  format(1x,'dens  = ',f12.5,' gaus',f12.5)
 918  format(1x,'dens  = ',f12.5,' +/- ',f12.5)


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

