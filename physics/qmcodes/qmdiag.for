      program main
c----------------------------------------------------------------------------c
c     direct diagonalization of quantum mechanical anharmonic oscillator.    c
c----------------------------------------------------------------------------c
c     version :            1.0                                               c
c     creation date :   07-23-96                                             c
c     last modified :   07-23-96                                             c
c----------------------------------------------------------------------------c
c     hamiltonian m/2(\dot x)^2+k(x^2-f^2)^2, units 2m=k=1.                  c
c----------------------------------------------------------------------------c
c     harmonic oscillator: H_0=m/2(\dot x)^2+m/2*w^2x^2                      c
c     perturbation:        H_1=A*x^4+B*x^2+C                                 c
c----------------------------------------------------------------------------c
      parameter(nmax=1000)
      real h(0:nmax-1,0:nmax-1),e(0:nmax-1)
      real v(0:nmax-1,0:nmax-1)
      real psi(0:nmax-1),rho(0:nmax-1)
      real rho2(0:nmax-1),rho3(0:nmax-1)
      real xcorp(0:nmax,0:nmax)
      real x3corp(0:nmax,0:nmax),x2corp(0:nmax,0:nmax)
      real m
       
      common /par/ m,w 
       
      write(6,*) 'parameter f'
      read(5,*)   f
      write(6,*) 'dimension of matrix (n>4)' 
      read(5,*)   ndim   
      write(6,*)  'unperturbed oscillator frequency (4*f)'
      read(5,*)   w0

c----------------------------------------------------------------------------c
c     output files, parameters                                               c
c----------------------------------------------------------------------------c

      open(unit=16,file='qmdiag.dat',status='unknown')
      open(unit=17,file='cor.dat',status='unknown')
      open(unit=18,file='cor2.dat',status='unknown')
      open(unit=19,file='z.dat',status='unknown')
      open(unit=20,file='psi.dat',status='unknown')
      open(unit=21,file='dcor.dat',status='unknown')

      eps = 1.e-30
      pi  = 3.1415926
 
      taumax = 2.5
      ntau = 100
      dtau = taumax/float(ntau)
            
      xmax = 2.0*f
      nx = 100
      dx = 2.0*xmax/float(nx)  
      
      write(16,*) 'qmdiag 1.0'
      write(16,*) '----------'
      write(16,601) f,ndim
      write(16,602) taumax,ntau
      write(16,603) xmax,nx

c----------------------------------------------------------------------------c
c     initialize parameters, clear arrays                                    c
c----------------------------------------------------------------------------c

      do 10 i=1,nmax
      do 10 j=1,nmax
         h(i,j) = 0.0
  10  continue
  
      m = 0.5
      w = w0
      
      a = 1.0
      b =-2.0*f**2-m*w**2/2.0
      c = f**4

      cw= 1.0/sqrt(m*w)
      
      c22 = cw**2/2.0 
      c44 = cw**4/4.0
            
c----------------------------------------------------------------------------c
c     build up h                                                             c
c----------------------------------------------------------------------------c

      do 100 n = 0,ndim-1
         
c----------------------------------------------------------------------------c
c     <n|h|n>                                                                c
c----------------------------------------------------------------------------c

         x4 = c44*3.0*((n+1)**2+n**2)
         x2 = c22*(2*n+1)
         e0 = w*(n+0.5) + c
         
         h(n,n) = a*x4 + b*x2 + e0

c----------------------------------------------------------------------------c
c     <n|h|n+2>                                                              c
c----------------------------------------------------------------------------c
         
         x4 = c44*sqrt((n+1.0)*(n+2))*(4*n+6)
         x2 = c22*sqrt((n+1.0)*(n+2))
         
         hh = a*x4 + b*x2
         h(n,n+2) = hh
         h(n+2,n) = hh
         
c----------------------------------------------------------------------------c
c     <n|h|n+4>                                                              c
c----------------------------------------------------------------------------c
              
         x4 = c44*sqrt((n+1.0)*(n+2)*(n+3)*(n+4))
         
         hh = a*x4     
         h(n,n+4) = hh     
         h(n+4,n) = hh
         
 100  continue
 
c----------------------------------------------------------------------------c
c     diagonalize                                                            c
c----------------------------------------------------------------------------c

      call jacobi(h,ndim,nmax,e,v,nrot)
      call eigsrt(e,v,ndim,nmax)

c----------------------------------------------------------------------------c
c     energy eigenvalues and matrix elements <0|x|n>                         c
c----------------------------------------------------------------------------c
              
      write(6,*) 
      write(16,*) 
      write(6,901) 
      write(16,901) 
      do 200 n=0,ndim-1
         cn = 0.0
         dn = 0.0
         en = 0.0
         do 205 k=0,ndim-1
            km3 = max(k-3,0)
            km2 = max(k-2,0)
            km1 = max(k-1,0)
            kp1 = min(k+1,ndim-1)
            kp2 = min(k+2,ndim-1)
            kp3 = min(k+3,ndim-1)
            cn = cn + ( sqrt(float(k))*v(km1,0)
     &          +sqrt(float(k+1))*v(kp1,0) )*v(k,n)
            dn = dn + ( sqrt(float(k*(k-1)))*v(km2,0)
     &          +(2*k+1)*v(k,0)
     &          +sqrt(float((k+1)*(k+2)))*v(kp2,0) )*v(k,n)  
            en = en + ( sqrt(float(k*(k-1)*(k-2)))*v(km3,0)
     &          +3*k*sqrt(float(k))*v(km1,0)
     &          +3*(k+1)*sqrt(float(k+1))*v(kp1,0)
     &          + sqrt(float((k+1)*(k+2)*(k+3)))*v(kp3,0) )*v(k,n)
 205     continue
         rho(n) = cw**2/2.0*cn**2
         rho2(n)= cw**4/4.0*dn**2
         rho3(n)= cw**6/8.0*en**2
         write(6,551)  n,e(n),rho(n),rho2(n),rho3(n)
         write(16,551) n,e(n),rho(n),rho2(n),rho3(n)
 200  continue

c----------------------------------------------------------------------------c
c     groundstate wave function                                              c
c----------------------------------------------------------------------------c

      write(16,*)
      write(16,902)    
      xnorm = 0.0
      xnorm2= 0.0
      do 210 k=0,nx
         x =-xmax+k*dx  
         psix = 0.0   
         call psiosc(ndim-1,x,psi)
         do 220 j=0,ndim-1
            psix = psix + v(j,0)*psi(j) 
 220     continue

c----------------------------------------------------------------------------c
c     compare to simple model                                                c
c----------------------------------------------------------------------------c

         psip = (2.0*f/pi)**0.25 * exp(-f*(x-f)**2)
         psim = (2.0*f/pi)**0.25 * exp(-f*(x+f)**2)
         psi0 = 1.0/sqrt(2.0)*(psip+psim)

c----------------------------------------------------------------------------c
c     check normalization                                                    c
c----------------------------------------------------------------------------c

         xnorm = xnorm + dx*psix**2
         xnorm2= xnorm2+ dx*psi(0)**2
         write(16,556) x,psix,psix**2,psi0,psi0**2
         write(20,556) x,psix,psix**2,psi0,psi0**2
 210  continue

      xnorm3 = 0.0
      do 230 j=0,ndim-1
         xnorm3 = xnorm3 + v(j,0)**2 
 230  continue

c     write(6,*)
c     write(6,*) 'norm ',xnorm,xnorm2,xnorm3
 
c----------------------------------------------------------------------------c
c     coordinate space correlator                                            c
c----------------------------------------------------------------------------c
      
      e0 = e(0)
      write(16,*)
      write(16,903)
      write(17,906)
      do 250 k=0,ntau   
         tau = k*dtau
         xcor= 0.0 
         do 260 j=1,n-1                   
            ej = e(j)
            xcor = xcor + rho(j)*exp(-(ej-e0)*tau)
            xcorp(k,j) = xcor
 260     continue
         write(16,555) tau,xcor,xcorp(k,1),xcorp(k,3),xcorp(k,5) 
         write(17,333) tau,xcor,0.01
 250  continue

c----------------------------------------------------------------------------c
c     log derivative                                                         c
c----------------------------------------------------------------------------c

      write(16,*)
      write(16,904)
      do 280 k=0,ntau
         tau = k*dtau
         xcor= xcorp(k,n-1)
         write(16,555) tau,log(xcor+eps),log(xcorp(k,1)+eps),
     &   log(xcorp(k,3)+eps),log(xcorp(k,5)+eps)
 280  continue

      write(16,*)
      write(16,905)
      do 290 k=0,ntau-1
         tau = k*dtau
         dlog =-(log(xcorp(k+1,n-1)+eps)-log(xcorp(k,n-1)+eps))/dtau
         dlog1=-(log(xcorp(k+1,1)+eps)-log(xcorp(k,1)+eps))/dtau
         dlog3=-(log(xcorp(k+1,3)+eps)-log(xcorp(k,3)+eps))/dtau
         dlog5=-(log(xcorp(k+1,5)+eps)-log(xcorp(k,5)+eps))/dtau
         write(16,555) tau,dlog,dlog1,dlog3,dlog5
 290  continue

c----------------------------------------------------------------------------c
c     x^2,x^3 correlator                                                     c
c----------------------------------------------------------------------------c

      write(16,*) 
      write(16,907)
      write(18,908)
      do 300 k=0,ntau   
         tau = k*dtau
         x1cor = 0.0
         x2cor = 0.0 
         x3cor = 0.0
         dx1cor= 0.0
         dx2cor= 0.0
         dx3cor= 0.0
         do 310 j=0,n-1                   
            ej = e(j)
            x1cor = x1cor + rho(j) *exp(-(ej-e0)*tau)
            x2cor = x2cor + rho2(j)*exp(-(ej-e0)*tau)
            x3cor = x3cor + rho3(j)*exp(-(ej-e0)*tau)
            dx1cor= dx1cor+ rho(j) *(ej-e0)*exp(-(ej-e0)*tau)
            dx2cor= dx2cor+ rho2(j)*(ej-e0)*exp(-(ej-e0)*tau)
            dx3cor= dx3cor+ rho3(j)*(ej-e0)*exp(-(ej-e0)*tau)
            x2corp(k,j) = x2cor
            x3corp(k,j) = x3cor
 310     continue
         dx1 = dx1cor/x1cor
         dx2 = dx2cor/(x2cor-rho2(0))
         dx3 = dx3cor/x3cor
         write(16,555) tau,x2cor,x2corp(k,0),x2corp(k,2),x2corp(k,4) 
         write(18,333) tau,x2cor,0.01
         write(21,777) tau,x1cor,x2cor,x3cor,dx1,dx2,dx3
 300  continue

      write(16,*)
      write(18,*) 
      write(16,909)
      write(18,910)
      do 320 k=0,ntau   
         tau = k*dtau
         x3cor = x3corp(k,n-1)
         write(16,555) tau,x3cor,x3corp(k,1),x3corp(k,3),x3corp(k,5) 
         write(18,333) tau,x3cor,0.01
 320  continue

c----------------------------------------------------------------------------c
c     partition function                                                     c
c----------------------------------------------------------------------------c

      xlmax = 100.0
      xlmin = 0.1
      xlogmax = log(xlmax)
      xlogmin = log(xlmin)
      nl = 50
      dlog = (xlogmax-xlogmin)/float(nl)
      
      do 400 il=0,nl
         xlog = xlogmin+il*dlog
         xl = exp(xlog)
         t  = 1.0/xl
         z  = 1.0
         do 410 i=1,ndim-1
            z = z + exp(-(e(i)-e(0))*xl)
 410     continue
         p = t*log(z)-e(0)
         write(19,333) t,xl,p
 400  continue
  

 222  format(2(1x,f12.5))
 331  format(1x,i4,2(f12.5))
 333  format(3(1x,f12.5))
 444  format(4(1x,f12.5))
 551  format(1x,i4,4(f12.5))
 555  format(5(1x,g12.5))
 556  format(5(1x,f12.5))
 777  format(7(1x,f12.5))
 601  format(' f   = ',f12.4,' ndim = ',i5)
 602  format(' t_m = ',f12.4,' ntau = ',i5)
 603  format(' x_m = ',f12.4,' nx   = ',i5)
 901  format(2x,'n',5x,'E_n',10x,'|c_n|^2')
 902  format(1x,'x',8x,'psi(x)',8x,'psi(x)^2')
 903  format(3x,'t',12x,'x(0)x(t)',5x,'1 state',5x,'3 states')
 906  format('#',2x,'t',12x,'x(0)x(t)')
 904  format(3x,'t',12x,'log x(0)x(t) ',1x,'1 state',5x,'3 states')
 905  format(3x,'t',12x,'d log(pi)/dt ',1x,'1 state',5x,'3 states')
 907  format(3x,'t',12x,'x^2(0)x^2(t)',1x,'1 state',5x,'3 states')
 908  format(3x,'t',12x,'x^2(0)x^2(t)')
 909  format(3x,'t',12x,'x^3(0)x^3(t)',1x,'1 state',5x,'3 states')
 910  format(3x,'t',12x,'x^3(0)x^3(t)')

      end

c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+-----

      subroutine hermite(n,x,p)
c----------------------------------------------------------------------------c
c     hermite polynomials p(i=0,..,n)=H_i(x)                                 c
c----------------------------------------------------------------------------c
      real p(0:n)
      
      p(0) = 1.0
      p(1) = 2.0*x
      
      do 10 i=2,n
         p(i) = 2.0*x*p(i-1) - 2*(i-1)*p(i-2)
  10  continue
      
      return
      end

c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+-----

      subroutine hermite2(n,x,p)
c----------------------------------------------------------------------------c
c     rescaled hermite polynomials H_i(x)/2^i                                c
c----------------------------------------------------------------------------c
      real p(0:n)
      
      p(0) = 1.0
      p(1) = x
      
      do 10 i=2,n
         p(i) = x*p(i-1) - (i-1)*p(i-2)/2.0
  10  continue
      
      return
      end

c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+-----

      subroutine hermite3(n,x,p)
c----------------------------------------------------------------------------c
c     rescaled hermite polynomials H_i(x)/2^i/sqrt(i!)                       c
c----------------------------------------------------------------------------c
      real p(0:n)
      
      p(0) = 1.0
      p(1) = x
      
      do 10 i=2,n
         p(i) = ( x*p(i-1) - sqrt(i-1.0)*p(i-2)/2.0 )/sqrt(float(i))
  10  continue
      
      return
      end
      
c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+-----

      subroutine psiosc(n,x,psi)
c----------------------------------------------------------------------------c
c     harmonic oscillator wave functions psi(i=0,..,n)=psi_i(x)              c
c----------------------------------------------------------------------------c
      parameter(nmax=1000)
      real psi(0:n),h(0:nmax-1)
      real m 
       
      common /par/ m,w
                    
      pi= 3.1415926    
      y = sqrt(m*w)*x
c     call hermite(n,y,h) 
c     call hermite2(n,y,h) 
      call hermite3(n,y,h) 
       
      do 10 i=0,n
c        xnorm  = (m*w/pi)**0.25/(2.0**(i/2.0)*snfac(i))  
c        xnorm  = (m*w/pi)**0.25*2.0**(i/2.0) /snfac(i)  
         xnorm  = (m*w/pi)**0.25*2.0**(i/2.0)  
         psi(i) = xnorm * h(i) * exp(-m*w/2.0*x**2)
  10  continue
  
      return
      end

c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+----- 
      
      function nfac(n)
c----------------------------------------------------------------------------c
c     n factorial                                                            c
c----------------------------------------------------------------------------c

      nfac = 1
      
      do 10 i=2,n
         nfac = nfac*i
  10  continue   
  
      return                                                           
      end                                                                     

c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+-----

      function snfac(n)
c----------------------------------------------------------------------------c
c     sqrt(n!)                                                               c
c----------------------------------------------------------------------------c

      snfac = 1
      
      do 10 i=2,n
         snfac = snfac*sqrt(float(i))
  10  continue   
  
      return                                                           
      end                                                                     

c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+----- 

      subroutine jacobi(a,n,np,d,v,nrot)
c----------------------------------------------------------------------------c
c     diagonalize real symmetric n*n matrix a using jacobi method.           c
c----------------------------------------------------------------------------c
c     a(n,n)    real, symmetric n*n matrix                                   c
c     n         dimension of matrix                                          c
c     np        dimension of arrays                                          c
c     d(n)      eigenvalues                                                  c
c     v(n,n)    eigenvectors v(i,1),...,v(n,i)                               c
c     nrot      number of jacobi rotations                                   c
c----------------------------------------------------------------------------c
c     note that the matrix a(n,n) is destoyed !                              c
c----------------------------------------------------------------------------c
      INTEGER n,np,nrot,NMAX
      REAL    a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      end

c-----------------------------------------------------------------------+-----
c-----------------------------------------------------------------------+-----
      
      subroutine eigsrt(d,v,n,np)
c----------------------------------------------------------------------------c
c     order eigenvectors and eigenvalues found by jacobi.                    c
c----------------------------------------------------------------------------c
      INTEGER n,np
      REAL d(np),v(np,np)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      end











