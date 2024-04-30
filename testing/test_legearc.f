        implicit real *8 (a-h,o-z)
        real *8 xs(10 000),ys(10 000),ts(10 000),whtsout(10 000),
     1      tsout(10 000),cxs(0:10 000),cys(0:10 000),xsout(10 000),
     2      ysout(10 000)
        complex *16 czs(0:10000),zsout(10000),ima,zwhtsout(10000)
c
c
        done=1
        pi=4*atan(done)
        ima=(0,1)
c
        call prini(6,13)
        print *, 'enter n:'
        read *, n
        call prinf('n=*',n,1)
c
        do 1600 i=0,n
        cxs(i)=.4d0**i
        cys(i)=(-.3d0)**i
 1600 continue
c
        call prin2('x coefs=*',cxs,n+1)
        call prin2('y coefs=*',cys,n+1)
c
        m=1000
        a=-1
        b=1
        h=(b-a)/(m-1)
c
        do 2000 i=1,m
        t=a+h*(i-1)
        ts(i)=t
        call legeexev(t,xs(i),cxs,n)
        call legeexev(t,ys(i),cys,n)
 2000 continue
c
        iw=11
        itype=3
cccc        call quagraph(iw,xs,ys,m,itype,'parameterized curve*')
c
        iw=12
        itype=3
cccc        call quagraph(iw,ts,xs,m,itype,'x value of curve*')
c
        iw=13
        itype=3
cccc        call quagraph(iw,ts,ys,m,itype,'y value of curve*')

c
c       now calculate the legendre nodes in arc length along the curve
c
        k=16
        call legearcnodes(ier,n,cxs,cys,k,xsout,ysout,whtsout,
     1      tsout,rlout)
c
        print *
        call prin2('rlout=*',rlout,1)
        call prin2('tsout=*',tsout,k)
        call prin2('whtsout=*',whtsout,k)
        call prin2('xsout=*',xsout,k)
        call prin2('ysout=*',ysout,k)


        iw=21
        itype2=2
        itype3=3
        call quagraph2(iw,xs,ys,m,itype3,xsout,ysout,k,itype2,
     1      'parameterized curve and legendre nodes*')

c
c       make complex and compare
c
        do 2200 i=0,n
        czs(i)=cxs(i)+ima*cys(i)
 2200 continue
c
        call zlegearcnodes(ier,n,czs,k,zsout,zwhtsout,
     1      tsout,rlout)
c
        print *
        call prin2('from zlegearcnodes, rlout=*',rlout,1)
        call prin2('from zlegearcnodes, tsout=*',tsout,k)
        call prin2('from zlegearcnodes, zwhtsout=*',zwhtsout,2*k)
        call prin2('from zlegearcnodes, zsout=*',zsout,2*k)

        stop
        end
c
c
