      program main
c--------------------------------------------------------------------------c
c     simple implementation of historic maximum entropy spectral function  c
c     reconstruction.                                                      c
c--------------------------------------------------------------------------c
c     version            :    1.0                                          c
c     creation date      :  02-05-97                                       c
c     last modification  :  02-09-97                                       c
c--------------------------------------------------------------------------c
c     reconstruct spectrum of anharmonic oscillator.                       c
c     input: file cor.dat                                                  c
c     title                                                                c
c     tau(i=1,..,ndat)  cor(i=1,..,ndat)  dcor(i=1,..,ndat)                c
c     output: file spec_rec.dat                                            c
c     title                                                                c
c     e(i=1,..,nspec)   rho(e)                                             c
c--------------------------------------------------------------------------c
      parameter(nsmax=100,ncmax=1000,nlmax=100)
      implicit real*8 (a-h,o-z)
      real*8 rho(nsmax),rho0(nsmax),rhog(nsmax)
      real*8 achi2(nlmax), aentr(nlmax), alam(nlmax)
      real*8 rholam(nlmax,nsmax)

      common /data/ taudat(ncmax),cordat(ncmax),cordate(ncmax)
      common /data/ ndat,nspec,e0
      common /ref/  rho0

      open(unit=17,file='cor.dat',status='unknown')
      open(unit=18,file='spec_rec.dat',status='unknown')
      open(unit=19,file='chi2.dat',status='unknown')
      open(unit=20,file='spec_lam.dat',status='unknown')

      pi= 3.1415926
   
c--------------------------------------------------------------------------c
c     input                                                                c
c--------------------------------------------------------------------------c
    
      write(6,*) 'number of data points (ndat)'
      read(5,*) ndat   
      write(6,*) 'maximum energy'
      read(5,*) e0
      write(6,*) 'number of points in spectral fct (nspec)'
      read(5,*) nspec 
      write(6,*) 'rescale errors'
      read(5,*) escale         
      write(6,*) 'default model: w0'
      read(5,*) w0
      write(6,*) 'default model: x^2'
      read(5,*) x2
      write(6,*) 'maximum value of lambda'
      read(5,*) xlam0
      write(6,*) 'minimum value of lambda'
      read(5,*) xlam1
      write(6,*) 'number of lambda values'
      read(5,*) nlam

c--------------------------------------------------------------------------c
c     read correlation function                                            c
c--------------------------------------------------------------------------c
   
      read(17,*)
      do 10 i=1,ndat 
         read(17,*) taudat(i),cordat(i),cordate(i)
         cordate(i) = escale*cordate(i)
  10  continue

c--------------------------------------------------------------------------c
c     grid for spectral function                                           c
c--------------------------------------------------------------------------c

      de= e0/float(nspec)

c--------------------------------------------------------------------------c
c     initial guess and default model                                      c
c--------------------------------------------------------------------------c 

      do 20 i=1,nspec
         e = i*de
         rho0(i) = x2*4.0/(pi*w0**3)*e**2/(1.0+(e/w0)**2)**2
         rhog(i) = rho0(i)
  20  continue

c--------------------------------------------------------------------------c
c     perform inversion for several values of lambda                       c
c--------------------------------------------------------------------------c

      dlam = (xlam0-xlam1)/float(nlam-1)

      ncross = 0

      write(6,*) ' xlam, chi2, entr'
      do 50 i=1,nlam
         xlam = xlam0-(i-1)*dlam
         call invert(rhog,rho,chi2,xlam,entr)
         write(6,333) xlam,chi2,entr 
         alam(i) = xlam
         achi2(i)= chi2
         aentr(i)= entr
         do 55 j=1,nspec
            rholam(i,j) = rho(j)
            rhog(j) = rho(j)
  55     continue
         if (i .gt. 1 .and. achi2(i-1) .gt. 1 .and.
     &       achi2(i) .lt. 1 .and. ncross .eq. 0) then
            xlambd1 = alam(i-1)
            xlambd2 = alam(i)
            chi2bd1 = achi2(i-1)
            chi2bd2 = achi2(i)
            ncross = 1
         endif
  50  continue

c--------------------------------------------------------------------------c
c     determine lambda from historic max ent principle chi2/ndf=1          c
c--------------------------------------------------------------------------c

      if (ncross .eq. 0 .and. achi2(1) .lt. 1) then
            xlambd1 = alam(1)
            xlambd2 = alam(2)
            chi2bd1 = achi2(1)
            chi2bd2 = achi2(2)
      else if (ncross .eq. 0) then
            xlambd1 = alam(nlam-1)
            xlambd2 = alam(nlam)
            chi2bd1 = achi2(nlam-1)
            chi2bd2 = achi2(nlam)
      endif 

      xlammag = xlambd1 + 
     &     (1.0-chi2bd1)/(chi2bd2-chi2bd1)*(xlambd2-xlambd1)
      xlammag = max(xlammag,0.0d0)
      xlammag = min(xlammag,2.0*xlam0)

c--------------------------------------------------------------------------c
c     final inversion for this value of lambda                             c
c--------------------------------------------------------------------------c

      xlam = xlammag
      call invert(rhog,rho,chi2,xlam,entr) 
      write(6,333) xlam,chi2,entr 

c--------------------------------------------------------------------------c
c     output, reconstructed spectrum                                       c
c--------------------------------------------------------------------------c
      
      write(18,*) nspec
      write(18,*) '  e           rho(e)'
      do 100 i=1,nspec
         e = i*de
         write(18,222) e,rho(i)
 100  continue

c--------------------------------------------------------------------------c
c     output, original data and fit from reconstructed spectrum            c
c--------------------------------------------------------------------------c

      write(18,*)      
      write(18,*) ndat  
      write(18,*) ' tau          cordat(tau)    cordate(i)',
     &            '    correc(tau)' 
      do 110 i=1,ndat
         tau  = taudat(i) 
         correc = cor(tau,rho,nspec,e0)
         write(18,444) tau,cordat(i),cordate(i),correc
 110  continue
      write(18,*) 
      write(18,*) ' xlam = ',xlam
      write(18,*) ' chi2 = ',chi2
      write(18,*) ' entr = ',entr

c--------------------------------------------------------------------------c
c     output, reconstructed spectrum for different lamdas                  c
c--------------------------------------------------------------------------c

      write(19,*) ' lam          chi2            entr' 
      write(20,*) ' nlam     nspec'
      write(20,*) nlam,nspec
      do 120 i=1,nlam
         write(19,333) alam(i),achi2(i),aentr(i)
         write(20,*) ' xlam = ',alam(i)
         write(20,*) ' e         rho(e)'  
         do 125 j=1,nspec
            e = j*de
            write(20,222) e,rholam(i,j)
 125     continue      
 120  continue
      write(20,*) ' xlam0 = ',xlam
      write(20,*) ' e         rho(e)'  
      do 128 j=1,nspec
         e = j*de
         write(20,222) e,rho(j)
 128  continue      

 111  format(1x,f12.5)
 222  format(2(1x,f12.5))
 333  format(3(1x,f12.5))
 444  format(4(1x,f12.5))
 555  format(5(1x,f12.5))
 666  format(6(1x,f12.5))
 999  format(70(1h-))
      end              
      
c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine invert(rhog,rho,chi20,xlam0,entr0)
c--------------------------------------------------------------------------c
c     perform regularized inversion for fixed value of lambda.             c
c--------------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      parameter(nsmax=100,ncmax=1000)

      real*8 xi(nsmax,nsmax)
      real*8 rhog(1),rho(1)

      common /data/ taudat(ncmax),cordat(ncmax),cordate(ncmax)
      common /data/ ndat,nspec,e0
      common /mem/  chi2,entr,xlam

      xlam = xlam0
      ndim = nsmax

c--------------------------------------------------------------------------c
c     initialize                                                           c
c--------------------------------------------------------------------------c

      do 10 i=1,nspec
      do 20 j=1,nspec
         xi(i,j) = 0.0 
  20  continue
      xi(i,i) = 1.0
  10  continue    

      do 30 i=1,nspec
         rho(i) = rhog(i)
  30  continue

      ftol = 1.d-5

c--------------------------------------------------------------------------c
c     inversion for fixed xlam, return chi2 and entropy                    c
c--------------------------------------------------------------------------c

      call powell(rho,xi,nspec,ndim,ftol,iter,fret)

      chi20 = chi2
      entr0 = entr

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      function func(rho)
c--------------------------------------------------------------------------c
c     function to be minimized for spectral fct reconstruction.            c
c--------------------------------------------------------------------------c
      parameter(nsmax=100,ncmax=1000)
      implicit real*8 (a-h,o-z)
      real*8 rho(1)
      real*8 taudat(ncmax),cordat(ncmax),cordate(ncmax)     
      real*8 rho0(nsmax)

      common /data/ taudat,cordat,cordate
      common /data/ ndat,nspec,e0
      common /ref/  rho0
      common /mem/  chi2,entr,xlam
      
      pi = 3.1415926d0
      de = e0/float(nspec)

c--------------------------------------------------------------------------c
c     chi2 for fit to correlator data                                      c
c--------------------------------------------------------------------------c 

      chi2 = 0.0
      do 10 i=1,ndat
         tau  = taudat(i) 
         corm = cor(tau,rho,nspec,e0)
         chi2 = chi2 + (cordat(i)-corm)**2/cordate(i)**2
  10  continue
      
      chi2 = chi2/float(nspec)

c--------------------------------------------------------------------------c
c     enforce positivity                                                   c
c--------------------------------------------------------------------------c

      xneg = 0.0
      do 20 i=1,ndat
         if (rho(i) .lt. 0.0) xneg = xneg + rho(i)**2
  20  continue
      xneg = xneg*10000.0

c--------------------------------------------------------------------------c
c     calculate entropy with respect to default model                      c
c--------------------------------------------------------------------------c

      entr = 0.0
      do 40 i=1,nspec
         entr = entr + rho(i)*log(abs(rho(i)/rho0(i)))
     &         - (rho(i)-rho0(i))
  40  continue
      entr = entr*de

c--------------------------------------------------------------------------c
c     function to be minimized                                             c
c--------------------------------------------------------------------------c 
     
      func = chi2 + xlam*entr + xneg
c     write(6,444) ' chi2 ',func,chi2,xneg,entr

 444  format(a5,4(1x,g12.5))
      return 
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      function cor(tau,rho,n,e0)
c--------------------------------------------------------------------------c
c     calculates coordinate space correlators for spectral function data   c
c     given on a grid rho(n).                                              c
c--------------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      real*8 rho(1)

      pi = 3.1415926
      de = e0/float(n)
      cor = 0.d0

      do 10 i=1,n
         e  = de*i
         cor = cor + de*rho(i)*exp(-e*tau)
 10   continue

      return
      end

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      subroutine powell(p,xi,n,np,ftol,iter,fret)
c--------------------------------------------------------------------------c
c     multidimensional minimization using powell's method, from num rec.   c
c--------------------------------------------------------------------------c
c     minimizes function func(x), where x is an n-dimesnsional vector.     c
c--------------------------------------------------------------------------c
c     p(n)      initial guess                                              c
c     xi(n,n)   initial direction (may take unit vectors)                  c
c     np        array size of xi                                           c
c     ftol      fractional tolerance                                       c
c     iter      number of iterations (output)                              c
c     fret      function value (output)                                    c
c--------------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      INTEGER iter,n,np,NMAX,ITMAX
      REAL*8 fret,ftol,p(np),xi(np,np),func
      EXTERNAL func
      PARAMETER (NMAX=100,ITMAX=500)
CU    USES func,linmin
      INTEGER i,ibig,j
      REAL*8 del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
      fret=func(p)
      do 11 j=1,n
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        call linmin(p,xit,n,fret)
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
      if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      SUBROUTINE linmin(p,xi,n,fret)
      implicit real*8 (a-h,o-z)
      INTEGER n,NMAX
      REAL*8 fret,p(n),xi(n),TOL
      PARAMETER (NMAX=100,TOL=1.d-6)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      REAL*8 ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      FUNCTION f1dim(x)
      implicit real*8 (a-h,o-z)
      INTEGER NMAX
      REAL*8 f1dim,func,x
      PARAMETER (NMAX=100)
CU    USES func
      INTEGER j,ncom
      REAL*8 pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      implicit real*8 (a-h,o-z)
      REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.d-30)
      REAL*8 dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END

c----------------------------------------------------------------------+----
c----------------------------------------------------------------------+----

      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      implicit real*8 (a-h,o-z)
      INTEGER ITMAX
      REAL*8 brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100,CGOLD=0.3819660d0,ZEPS=1.0d-10)
      INTEGER iter
      REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END
