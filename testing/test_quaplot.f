        implicit real *8 (a-h,o-z)
        dimension x1(1000),y1(1000),x2(1000),y2(1000),
     1      x3(1000),y3(1000),ts(1000)
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
        PRINT *, 'ENTER ntest'
        READ *,ntest
        CALL PRINf('ntest=*',ntest,1 )

c
c        construct the nodes to be plotted
c

        done=1
        pi=atan(done)*4
        h=2*pi/ntest
        do 1200 i=1,ntest
c
        ts(i)=(i-1)*h
c
        x1(i)=cos(ts(i))
        y1(i)=sin(ts(i))
c
        x2(i)=x1(i)*2
        y2(i)=y1(i)*2
c
        x3(i)=x2(i)*2
        y3(i)=y2(i)*2
 1200 continue
c
c        plot the first circle
c
        itype1=3
        iw=21
        call quaplot(iw,x1,y1,ntest,itype1,'first circle*')

        itype1=1
        itype2=2
        itype3=3
        iw=22

        n2=ntest/2

        call quaplot2(iw,x1,y1,ntest,itype1,x2,y2,n2,itype2,
     1      'first two curves*')

        iw=23
        n3=ntest/3

        call quaplot3(iw,x1,y1,ntest,itype1,x2,y2,n2,itype2,
     1      x3,y3,n3,itype3,'all three curves*')


c
c
c
        iw=31

        itype1=3


        call quagraph(iw,ts,x1,ntest,itype1,
     1      'first graph*')


        iw=32

        itype1=3
        itype2=2

        call quagraph2(iw,ts,x1,ntest,itype1,ts,y1,n2,itype2,
     1      'first two curves*')


        iw=33

        itype1=1
        itype2=2
        itype3=3

        call quagraph3(iw,ts,x1,ntest,itype1,ts,y1,n2,itype2,
     1      ts,x2,n3,itype3,'first two curves*')



        stop
        end
