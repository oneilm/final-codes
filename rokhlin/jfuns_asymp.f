        implicit real *8 (a-h,o-z)
        complex *16 fjs(100 000),zz,x,z,f66,ima,f1
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER rnu'
        READ *,rnu
        CALL PRIN2('rnu=*',rnu,1 )
c 
c        This code documents the asymptotic expansion for Bessel
c        functions (J_{nu}) for large nu, valid to 1/(nu)^5. The
c        said expansion is obtained by elementary analysis from
c        the formula on page 227 of Watson's "Bessel Functions"
c 
  
        nu=rnu+0.001
        x=3*ima+2
cccc        x=3*ima
c 
        done=1
        pi=atan(done)*4
        e=exp(done)
c 
        z=x/rnu
        t=1/rnu
c 
c       evaluate the Bessel function
c 
        zz=x
        eps=1.0d-3000
        call JBFUN(zZ,fjs(2),N,EPS)
  
  
        f1=fjs(nu+1)
  
        call prin2('f1 (Bessel function) =*',f1,1)
c 
c        evaluate the J-function via a 4-term expansion
c        and compare the resultwith J_{nu} evaluated directly
c 
        f66=1 - (t*(1 + 3*x**2))/12 + (t**2*(1 + 78*x**2 + 9*x**4))/288
     1        - (t**3*(-139 + 14085*x**2 + 4995*x**4 + 135*x**6))/51840
     2       + (t**4*(-571 + 674412*x **2 + 564030*x**4 + 39420*x**6 +
     3       405*x**8))/2488320
  
 1500 continue
  
        f66=f66*(x*e/2/rnu)**rnu/sqrt(2*pi*rnu)
  
  
  
  
  
        call prin2('and f66 =*',f66,2)
        call prin2('and f66-f1)/f1 =*',(f66-f1)/f1,2)
        call prin2('and abs((f66-f1)/f1) =*',abs((f66-f1)/f1),1)
  
        stop
        end
  
c 
c 
c 
c 
c 
        function facts(nu)
        implicit real *8 (a-h,o-z)
c 
        save
        d=1
        do 1200 i=1,nu
c 
        d=d*i
 1200 continue
c 
        facts=d
  
        return
        end
  
  
  
  
  
