        implicit double precision (a-h,o-z)
        double precision xs(40 000 000),ys(40 000 000),zs(40 000 000)
        DIMENSION vals(1000 000),bins(1000 000),fs(1000 000),
     1      ixs(1000 000)
        DIMENSION vals2(1000 000),bins2(1000 000),fs2(1000 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1)
c
        call corrand3(N,xs)
cccc        call corrand(N,xs)

        k=30
cccc        k=200
cccc        k=100

        iw=31
        call hystog0(N,xs,k,iw,bins,vals)

        iw=32
        call zquaplot(iw,xs,n/2,1,'area*')
c
c        construct the normal distributions
c
        call corrand_norm(N,xs,ys)
        call corrand_norm(N,ys,zs)


        iw=33
        call hystog0(N,xs,k,iw,bins,vals)
c
        iw=34
        call hystog0(N,ys,k,iw,bins2,vals2)
c
c        construct the distributions analytically and plot
c        them things together with the experimental ones
c
        done=1
        pi=atan(done)*4
c
        h=bins(2)-bins(1)
        rint=0
        do 2200 i=1,k
c
        fs(i)=exp(-bins(i)**2/2) /sqrt(pi*2)
        rint=rint+fs(i)
 2200 continue
c
        rint=rint*h
        call prin2('and rint=*',rint,1)
c
c
        iw2=53
        call quagraph2(iw2,bins,vals,k-1,3,bins,fs,k-1,3,
     1      'hystogram for vals vs. analytical result *')
c
        h=bins2(2)-bins2(1)
        rint2=0
        do 2400 i=1,k
c
        fs2(i)=exp(-bins2(i)**2/2) /sqrt(pi*2)
        rint2=rint2+fs2(i)
 2400 continue
c
        rint2=rint2*h
        call prin2('and rint2=*',rint2,1)
c
        iw2=54
        call quagraph2(iw2,bins2,vals2,k-1,3,bins2,fs2,k-1,3,
     1      'hystogram for vals2 vs. analytical result *')
c
        call arbscap(xs,xs,n,d)

        call prin2('(xs,xs)=*',d,1)

        call arbscap(xs,ys,n,d2)

        call prin2('(xs,ys)=*',d2,1)
c
c        test the integer number generator
c 
        call corrand_integer_knuth(n,ixs)
        call prinf('after _integr, ixs=*',ixs,n)

        call corrand_integer_knuth(n,ixs)
        call prinf('after _integr, ixs=*',ixs,n)
c
ccc        call bubba(N,ixs)

ccc        call prinf('after bubba, ixs=*',ixs,n)



        do 3200 i=1,n
c
        fs2(i)=ixs(i)
 3200 continue


        iw=36
        call zquaplot(iw,fs2,n/2,2,'area*')
c
c       test the gold-plated random number generator
c
        call corrand_norm(N,ys,zs)

        call prin2('ys as created*',ys,20)

        done=1
        pi=atan(done)*4

        dd=sqrt(2*done)

        do 4200 i=1,n
c
cc        zs(i)=exp(-ys(i)**2/2)/sqrt(2*pi)

cc        zs(i)=exp(-ys(i)**2/2)/sqrt(2*pi) *ys(i)

        zs(i)=-exp(-ys(i)**2/2)/sqrt(2*pi) *ys(i) 

        zs(i)=ys(i)/abs(zs(i))


       
        call qerrfun(ys(i)/dd,zs(i))

        zs(i)=zs(i)/2

cccc        zs(i)=1/zs(i)
 4200 continue


        call prin2('zs as created*',zs,20)

        iw=37
        k=300
        k=10
        call hystog0(N,zs,k,iw,bins2,vals2)





        stop
        end
C 
C 
C 
C 
C 
        SUBROUTINE bubba(N,ixs)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION ixs(1)
c
        do 1400 i=1,n
        do 1200 j=1,n-1
        if(ixs(j) .gt. ixs(j+1) ) then
            jj=ixs(j)
            ixs(j)=ixs(j+1)
            ixs(j+1)=jj
        endif
c
 1200 continue
 1400 continue


        return
        end
C 
C 
C 
C 
C 
        SUBROUTINE hystog0(N,xs,k,iw,bins,vals)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION xs(1),vals(1),bins(1)
c
c        determine the maximum value of the random variable
c
        xmax=-1.0d30
        xmin=1.0d30
        do 1200 i=1,n
c
        if(xmax .lt. xs(i)) xmax=xs(i)
        if(xmin .gt. xs(i)) xmin=xs(i)
 1200 continue
c
        call prin2('xmin=*',xmin,1)
        call prin2('xmax=*',xmax,1)
c
c       . . . bin them things
c
        h=(xmax-xmin)/k

        call prin2('h=*',h,1)
c
        done=1
        pi=atan(done)*4
c
        rint=0
        do 1400 i=1,k
c
        vals(i)=0
        bins(i)=xmin+h/2+(i-1)*h
c
 1400 continue


        call prin2('bins=*',bins,k)
c
        do 1600 i=1,n
c
        d=(xs(i)-xmin)/h
        j=d+1
        vals(j)=vals(j)+1
 1600 continue
c
        do 1800 i=1,k
c
        vals(i)=vals(i)/n*k/(xmax-xmin)
 1800 continue
c
        call prin2('vals=*',vals,k)
c
c       . . . plot
c
        call quagraph(iw,bins,vals,k-1,3,'hystogram*')
c
        return
        end
c
c 
c 
c 
c 
        subroutine arbscap(x,y,n,d)
        implicit double precision (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        d=0
        do 2000 i=1,n
        d=d+x(i)*y(i)
 2000 continue
        return
        end
  
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code, and the beginning 
c        of the random number code proper
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
C        This file contains four user-callable subroutines: corrand,
c        corrand2, corrand3, corrand_norm. Following is a brief
c        description of the said four subroutines.
c
c   corrand, corrand2, corrand3 - the functions of the 
c        subroutines corrand, corrand2, corrand3 are identical;
c        each of them returns to the user a pseudo-random vector 
c        distributed uniformly on the interval [0,1]. The algorithms
c        used by them are also identical, and are as follows: each
c        constructs three pseudo-random sequences using three 
c        congruential generators, averages them, and uses the 
c        obvious transformation to bring it back to the uniform
c        distribution. Each of the three uses its own three congru-
c        ential generators. Those in corrand2, corrand3 are copied
c        from the Numerical Recipies, and the ones in corrand are 
c        from nowhere, and are probably too short, anyway.
c
c   corrand_norm - returns to the user two random vectors y1, y2, 
c        each of which is distributed normally with the 
c        distribution density
c
c        1/sqrt(2 * \pi) * e^(y^2/2)                         (1)
c   
c   corrand_integer_knuth - for a user-specified n, returns to 
c        the user a random permutation of n integers 1,2,...,n
c
C 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c 
c 
        subroutine corrand_integer_knuth(n,ixs)
        implicit double precision (a-h,o-z)
        save
        dimension ixs(n),tt(12)
c
c        This subroutine returns to the user a random permutation 
c        of n integer numbers 1,2,...,n.
c
c              Input parameters:
c
c  n - random numbers in the array ixs will be the distributed
c        uniformly on the interval [1,n]
c
c              Output parameters:
c
c  ixs - a pseudo-random permutation of length n
c
        call corrand3(11,tt)
        call corrand3(11,tt)
        call corrand3(11,tt)
        do 1200 i=1,n
c
        ixs(i)=i
 1200 continue
c
        done=1
        do 1400 i=1,n-1
c
        call corrand3(1,tt)
c
        k=n-i+1

cccc        call prinf('k=*',k,1)

        h=done/k
        j=tt(1)/h+1

cccc        call prinf('and j=*',j,1)
c
        jj=ixs(k)
        ixs(k)=ixs(j)
        ixs(j)=jj
 1400 continue
c
        return
        end
c
c
c
c
c
        SUBROUTINE corrand_norm(N,Y1,y2)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION Y1(1),y2(1)
c
c        This subroutine returns to the user two random
c        vectors y1, y2, each of which is distiibuted
c        normally with the distribution density
c
c        1/sqrt(2 * \pi) * e^(y^2/2)                         (1)
c   
c
c              Input parameters:
c
c  n - the number of elements to be returned in each of the 
c        arrays y1, y2
c
c              Output parameters:
c  y1, y2 - two pseudo-random arrays distributed normally 
c        with the distribution density (1)
c
c
c        . . . construct vectors of variables distributed
c              uniformly on the interval [0,1]
c
        call corrand3(N,Y1)
        call corrand3(N,Y2)
c
c       combine the variables y1, y2 converting them
c       into variables distributed normally (Box-Muller 
c       algorithm)
c
        done=1
        pi=atan(done)*4
        do 1400 i=1,n
c
        z1=sqrt(-2*log(y1(i)))*cos(2*pi*y2(i))
        z2=sqrt(-2*log(y1(i)))*sin(2*pi*y2(i))
c
        y1(i)=z1
        y2(i)=z2
 1400 continue
c
        return
        end
c
c
c
c
c
        SUBROUTINE corrand3(N,Y)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
c
c        This subroutine returns to the user a collection of
c        reasonably random numbers uniformly distributed on
c        the interval [0,1]
c
C       . . . CONDUCT PRELIMINARY RANDOMIZATION
C 
        IF (IFCALL.EQ.1 ) GOTO  1300
        IFCALL =1
c
c        construct parameters for the first process
c
        LAMBDA=3661
        MU=30809
        IP=145800
        M=41
c
c        construct parameters for the second process
c
        LAMBDAq=8121
        MUq=28411
        IPq=134456
        Mq=43
c
c        construct parameters for the third process
c
        LAMBDAqq=3613
        MUqq=45289
        IPqq=214326
        Mqq=53
c
        DO 1200 I=1,100
        M1=M*LAMBDA  +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq  +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq  +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
 1200 CONTINUE
c
 1300 CONTINUE
c
        D=1
        D=D/IP
        Dq=1
        Dq=Dq/IPq
        Dqq=1
        Dqq=Dqq/IPqq
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        DO 1400 I=1,N
c
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
        Y(I)=M*D+Mq*Dq+mqq*dqq
c
        call corrand_comp(y(i),rint)
c
        y(i)=rint
 1400 CONTINUE
c
        RETURN
c
c
c
c
        entry corrand3_init(dummy)
c
        ifcall=0
        return
        END
c
c
c
c
c
        SUBROUTINE corrand2(N,Y)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
c
c        This subroutine returns to the user a collection of
c        reasonably random numbers uniformly distributed on
c        the interval [0,1]
c
C       . . . CONDUCT PRELIMINARY RANDOMIZATION
C 
        IF (IFCALL.EQ.1 ) GOTO  1300
        IFCALL =1
c
c        construct parameters for the first process
c
        LAMBDA=1255
        MU=6173
        IP=29282
        M=17
c
c        construct parameters for the second process
c
        LAMBDAq=421
        MUq=54773
        IPq=259200
        Mq=19
c
c        construct parameters for the third process
c
        LAMBDAqq=1366
        MUqq=150889
        IPqq=714025
        Mqq=31
c
        DO 1200 I=1,100
        M1=M*LAMBDA  +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq  +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq  +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
 1200 CONTINUE
c
 1300 CONTINUE
c
        D=1
        D=D/IP
        Dq=1
        Dq=Dq/IPq
        Dqq=1
        Dqq=Dqq/IPqq
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        DO 1400 I=1,N
c
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
        Y(I)=M*D+Mq*Dq+mqq*dqq
c
        call corrand_comp(y(i),rint)
c
        y(i)=rint
 1400 CONTINUE
c
        RETURN
        END
c
c
c
c
c
        SUBROUTINE corrand(N,Y)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
c
c        This subroutine returns to the user a collection of
c        reasonably random numbers uniformly distributed on
c        the interval [0,1]
c
C       . . . CONDUCT PRELIMINARY RANDOMIZATION
C 
        IF (IFCALL.EQ.1 ) GOTO  1300
        IFCALL =1
c
c        construct parameters for the first process
c
        LAMBDA=2**7 +9
        MU=1
        IP=2**20
        M=17
c
c        construct parameters for the second process
c
        LAMBDAq=2**5 +29
        MUq=3
        IPq=2**22
        Mq=19
c
c        construct parameters for the third process
c
        LAMBDAqq=2**3 +31
        MUqq=13
        IPqq=2**18
        Mqq=31
c
        DO 1200 I=1,100
        M1=M*LAMBDA  +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq  +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq  +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
 1200 CONTINUE
c
 1300 CONTINUE
c
        D=1
        D=D/IP
        Dq=1
        Dq=Dq/IPq
        Dqq=1
        Dqq=Dqq/IPqq
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        DO 1400 I=1,N
c
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
        Y(I)=M*D+Mq*Dq+mqq*dqq
c
        call corrand_comp(y(i),rint)
c
        y(i)=rint
 1400 CONTINUE
c
        RETURN
        END
c
c
c
c
c
        SUBROUTINE corrand_comp(x,rint)
        IMPLICIT double precision (A-H,O-Z)
        save
c
        if (x .lt. 1) then
            rint=x**3/6
            return
        endif
c      
        if ( (x .ge. 1) .and. (x .le. 2) ) then
            rint= 0.75d0*x-(x-1.5d0)**3/3 - 0.625d0
            return
        endif
c
        rint= (x-3)**3/6 +1
c
        return
        end
