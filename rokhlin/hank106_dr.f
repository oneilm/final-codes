        implicit real *8 (a-h,o-z)
        integer *4 istart(1000)
        complex *16 h0derss(100 000),h0s(10000),ima,z,
     1      h0,h1,h02,h12,h0stest(10000),h0stest2(10000),
     2      h1stest(10000),h1stest2(10000),diffs0(10000),
     3      diffs1(10000),u0,h1derss(100 000),
     4      h1s(10000),rk,rksav
c 
        dimension ab(2,10000),centers(2,10000),xstest(10000),
     1      rstest(10000)
  
        dimension w(1000 000)
c 
  
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
        ifexpon=1
C 
C       SET ALL PARAMETERS
C 
cc         PRINT *, 'ENTER n'
cc         READ *,n
cc         CALL PRINF('n=*',n,1 )
  
c 
c       construct all parameters
c 
        rk=ima/1000001
        rksav=rk/cdabs(rk)
c 
c 
c       initialize array
c 
c 
        nab = 10000
        lw = 1000000
ccc        call hank106datagen(rk,rmin,rmax,ab,nab,ninterval,w,lw,
ccc     1          istart,ier)
        call hank106datagen(rk,ier)
ccc        call prin2(' ab = *',ab,2*ninterval)
        call prinf(' ier = *',ier,1)
        r = 0.005
        z=r*rksav
ccc        z=ima*1.0d-6
c 
c        construct the scales that will convert the arrgument of
c        the taylor series to be evaluated into a real number
c 
  
  
ccc        call hank106(z,h0,h1,ab,ninterval,w,istart)
        call hank106(Z,H0,H1,ifexpon)
  
        call prin2('after hank106eval, h0=*',h0,2)
  
        call hank103(z,h02,h12,ifexpon)
  
        call prin2('and h02=*',h02,2)
        call prin2('and h02-h0=*',h02-h0,2)
  
  
        call prin2('after hank106eva, h1=*',h1,2)
  
        call prin2('and h12=*',h12,2)
        call prin2('and h12-h1=*',h12-h1,2)
  
ccc        stop
c 
c        conduct mass testing
c 
        rmin = 1.0d-6
        rmax = 2.0d-4
        ntest=1000
  
        call prin2('rmax =*',rmax,1)
        call prin2('rmin =*',rmin,1)
        htest=(rmax-rmin)/(ntest-1)
        call prin2('htest =*',htest,1)
c 
        do 2200 i=1,ntest
c 
cc        rstest(i)=(i-1)*htest+rmin+1.0d-12
        rstest(i)=(i-1)*htest+rmin
        call prin2('rstest =*',rstest(i),1)
 2200 continue
  
ccc        rstest(1)=rstest(1)+1.0d-12
ccc        rstest(ntest)=rstest(ntest)-1.0d-12
  
        call prin2('and rstest=*',rstest,ntest)
  
        time1=clotatim()
  
        do 2500 ijk=1,1000
c 
        do 2400 i=1,ntest
c 
        z=rstest(i)*rksav
  
  
ccc        call hank106(rstest(i),h0stest(i),h1stest(i),w)
ccc        call hank106(z,h0stest(i),h1stest(i),ab,ninterval,w,istart)
        call hank106(z,h0stest(i),h1stest(i),ifexpon)
ccc        write(17,*) rstest(i), dimag(h0stest(i))
c 
 2400 continue
 2500 continue
  
        time2=clotatim()
  
        call prin2('and time for hank106 is*',time2-time1,1)
  
        call prin2('h0stest=*',h0stest,ntest*2)
c 
  
        time3=clotatim()
  
        do 2700 ijk=1,1000
c 
        do 2600 i=1,ntest
c 
        z=rstest(i)*rksav
        call hank103(z,h0stest2(i),h1stest2(i),ifexpon)
  
 2600 continue
 2700 continue
  
        call prin2('h0stest2=*',h0stest2,ntest*2)
  
        time4=clotatim()
  
        call prin2('and time for hank103 is*',time4-time3,1)
  
  
        err0=0
        d0=0
        err1=0
        d1=0
  
        do 2800 i=1,ntest
  
        diffs0(i)=h0stest(i)-h0stest2(i)
        diffs1(i)=h1stest(i)-h1stest2(i)
c 
        d0=d0+abs(h0stest(i))**2
        call prin2(' d0 = *',d0,1)
        d1=d1+abs(h1stest(i))**2
c 
        err0=err0+abs(diffs0(i))**2
        call prin2(' err0 = *',err0,1)
        err1=err1+abs(diffs1(i))**2
        write(11,*) rstest(i), sqrt(abs(diffs0(i)))
 2800 continue
  
  
cccc        call prin2('diffs0=*',diffs0,ntest*2)
cccc        call prin2('diffs1=*',diffs1,ntest*2)
  
        errrel0=sqrt(err0/d0)
        errrel1=sqrt(err1/d1)
  
        call prin2('and errrel0=*',errrel0,1)
        call prin2('and errrel1=*',errrel1,1)
  
        call prinf('and, by the way, keep=*',keep,1)
  
        stop
        end
c 
c 
  
