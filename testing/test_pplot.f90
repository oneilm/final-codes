program test_pplot
  implicit real *8 (a-h,o-z)
  real *8 ts(10000),xs(10000),ys(10000),zs(10000)
  real *8, allocatable :: vexps(:,:)

  call prini(6,13)

  done=1
  pi=4*atan(done)
  ima=(0,1)

  print *, 'enter n:'
  read *, n
  call prinf('n=*',n,1)

  !c
  !c       test the plotting routine
  !c
  iw=51
  m1=500
  n1=200
  allocate( vexps(m1,n1) )

  x1=-5
  x2=5
  y1=-2
  y2=8

  xcen=-2
  ycen=7
  hx=(x2-x1)/n1
  hy=(y2-y1)/m1

  do  j=1,n1
     do  i=1,m1
        x=x1+hx*(j-1)
        y=y2-hy*(i-1)
        vexps(i,j)=exp(-((x-xcen)**2 + (y-ycen)**2)/2)
     end do
  end do

  call pyimage(iw,m1,n1,vexps,'incoming field*')

  !c
  !c       now plot the field as well as the scatterer
  !c
  x0=0
  y0=4
  do  i=1,n
     t=2*pi/n*(i-1)
     xs(i)=x0+cos(t)
     ys(i)=y0+sin(t)
  end do



  iw=52
  call pyimage2(iw,m1,n1,vexps,xs,ys,n,x1,x2,y1,y2, &
           'field and scatterer*')


  !c
  !c       now plot a curve
  !c
  a=0
  b=2*pi
  h=(b-a)/(n-1)

  do  i=1,n
     ts(i)=a+h*(i-1)
     xs(i)=cos(ts(i))+4*sin(ts(i))-3*cos(7*ts(i))
     ys(i)=exp(-ts(i)**2)
     zs(i)=exp(-ts(i))
  end do


  n2=n+23
  h=(b-a)/(n2-1)

  do  i=1,n2
     ts(i)=a+h*(i-1)
     ys(i)=exp(-ts(i)**2)
  end do

  
  iw=11
  itype=3
  call pyplot(iw,ts,xs,n,itype,'one test curve*')



  iw=12
  itype=1
  itype2=3
  call pyplot2(iw,ts,xs,n,itype,ts,ys,n2,itype2, &
       'two test curves*')


  stop
end program test_pplot
