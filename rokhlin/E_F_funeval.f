        implicit real *8 (a-h,o-z)
        dimension xs(12),ws(12),xs8(8),ws8(8)

c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER c'
        READ *,c
        CALL PRIN2('c=*',c,1 )
C 


        x=0.99999999
        a=0.999999
cc        a=0.97
cc        a=0.2
cccc        x=1
c
cccc        call Efun_eval(x,a,eps,f,dfdx)
        call Ffun_eval(x,a,eps,f,dfdx)



        call prin2('from efun_eval, f=*',f,1)

cccc        call Efun_eval2(x,a,f2,dfdx2)


cccc        call E_F_funeval(x,a,f2,dfdx2,fF,dfFdx)

        call E_F_funeval(x,a,f27,dfdx27,f2,dfdx2)

        call prin2('and f2=*',f2,1)
        call prin2('and f2-f=*',f2-f,1)
        call prin2('and f2/f=*',f2/f,1)



        call prin2('and dfdx=*',dfdx,1)
        call prin2('and dfdx2=*',dfdx2,1)
        call prin2('and dfdx2-dfdx=*',dfdx2-dfdx,1)



        stop
        end
c 
c 
c 
c 
c 
        subroutine E_F_funeval(x,a,fE,dfEdx,fF,dfFdx)
        implicit real *8 (a-h,o-z)
        save
        dimension xs9(9),ws9(9)
c
        data xs9/
     1  0.5638078146207924D-01,0.2160364885037746D+00,
     2  0.3933358751675277D+00,0.5573613589026810D+00,
     3  0.6986310146537872D+00,0.8132945846251771D+00,
     4  0.9002236800650369D+00,0.9596749768459339D+00,
     5  0.9923789681449862D+00/
c        
        data ws9/
     1  0.1281278118908609D+00,0.1774773860933696D+00,
     2  0.1730176845957494D+00,0.1536053653007898D+00,
     3  0.1283276930455867D+00,0.1008326420252267D+00,
     4  0.7308354814119854D-01,0.4594968670470709D-01,
     5  0.1957818220251128D-01/
c
c        This subroutine calculates the elliptic integrals
c
c        E(a,x)=\int_0^x sqrt(1-a^2*sin(t)^2) dt,               (1)
c
c
c        F(a,x)=\int_0^x 1/sqrt(1-a^2*sin(t)^2) dt,             (2)
c
c
c        together with the derivatives of the said elliptic 
c        integrals. The evaluation is accurate to the full
c        double precision.
c
        b=x
        fE=0
        fF=0
        do 2200 i=1,9
c
        d=sqrt(1-a**2*(sin(b*xs9(i)))**2)
        fE=fE+d*ws9(i)
        fF=fF+ws9(i)/d
 2200 continue
c
        fE=fE*x
        fF=fF*x
c
        d=1-a**2*sin(x)**2
        dfEdx=sqrt(d)
        dfFdx=1/dfEdx
c
        return
        end
c 
c 
c 
c 
c 
        subroutine Efun_eval2(x,a,f,dfdx)
        implicit real *8 (a-h,o-z)
        external fun17
        save
        dimension xs8(8),ws8(8)
c
        data xs8/
     1  0.1032871719713107D+00,0.3041758538824985D+00,
     2  0.4889882735353880D+00,0.6498439018087656D+00,
     3  0.7819071172835492D+00,0.8829752849658865D+00,
     4  0.9525722562021406D+00,0.9910228739984174D+00/
c        
        data ws8/
     1  0.2056092105909210D+00,0.1944094174036297D+00,
     2  0.1738933171454199D+00,0.1470201294455017D+00,
     3  0.1167555042514758D+00,0.8531763778642262D-01,
     4  0.5394108908269785D-01,0.2305369429393149D-01/
c
        b=x
        f=0
        do 2200 i=1,8
c
        f=f+sqrt(1-a**2*(sin(b*xs8(i)))**2)*ws8(i)
 2200 continue
c
        f=f*x
c
        d=1-a**2*sin(x)**2
        dfdx=sqrt(d)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine Ffun_eval(x,a,eps,f,dfdx)
        implicit real *8 (a-h,o-z)
        external fun19
        save
c
        eps=1.0d-15
cccc        eps=1.0d-22
        m=24
        done=1
        pi=atan(done)*4
c
        aa=0
        bb=x
c 
        call adapgaus(ier,aa,bb,fun19,a,par2,m,eps,
     1      rint,maxrec,numint)
        f=rint
c
        d=1-a**2*sin(x)**2
        dfdx=1/sqrt(d)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine Efun_eval(x,a,eps,f,dfdx)
        implicit real *8 (a-h,o-z)
        external fun17
        save
c
        eps=1.0d-15
cccc        eps=1.0d-22
        m=24
        done=1
        pi=atan(done)*4
c
        aa=0
        bb=x
c 
        call adapgaus(ier,aa,bb,fun17,a,par2,m,eps,
     1      rint,maxrec,numint)
        f=rint
c
        d=1-a**2*sin(x)**2
        dfdx=sqrt(d)
c
        return
        end
c 
c 
c 
c 
c 
        function fun17(x,a,dummy)
        implicit real *8 (a-h,o-z)
c 
        save
c
        d=1-a**2*sin(x)**2
        fun17=sqrt(d)
c
        return
        end
c 
c 
c 
c 
c 
        function fun19(x,a,dummy)
        implicit real *8 (a-h,o-z)
c 
        save
c
        d=1-a**2*sin(x)**2
        fun19=1/sqrt(d)
c
        return
        end

