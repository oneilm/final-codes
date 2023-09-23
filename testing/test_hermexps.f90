program test_hermexps
  implicit real *8 (a-h,o-z)
  dimension funs(1000),ders(1000),diffs(1000), &
       funsp(1000),funsm(1000),dersp(1000),dersm(1000)

  dimension w(10000),rints(1000),rintsp(1000), &
       rintsm(1000),funs2(1000)

  double precision :: ts(10000), vals(10000)

  dimension xs(1000),uu(100000),vv(100000),whts(1000), &
           ww(100000),coefs(1000)

  call prini(6,13)
  !
  !C       SET ALL PARAMETERS
  !C
  PRINT *, 'ENTER n'
  READ *,n
  CALL PRINf('n=*',n,1 )

  nmax=n

  h = 0.00001d0
  x=2
  !cccc        x=0
  ifinit=1

  call herminte(x, funs, ders, rints,n,w,ifinit)

  call prin2('funs=*',funs,nmax)
  call prin2('ders=*',ders,nmax)
  call prin2('rints=*',rints,2)

  call herminte(x+h,funsp,dersp,rintsp,n,w,ifinit)
  call herminte(x-h,funsm,dersm,rintsm,n,w,ifinit)

  funs2(1) = (rintsp(1)-rintsm(1))/h/2

  !cccc        call prin2('and rintsp=*',rintsp,1)
  !cccc        call prin2('and rintsm=*',rintsm,1)

  funs2(2)=(rintsp(2)-rintsm(2))/h/2

  call prin2('and funs2=*',funs2,2)
  call prin2('and funs2(1)/funs(1)=*',funs2(1)/funs(1),1)
  call prin2('and funs2(2)/funs(2)=*',funs2(2)/funs(2),1)

  do  i=1,n
     funs2(i)=(rintsp(i)-rintsm(i))/h/2
     diffs(i)=funs2(i)-funs(i)
  end do


  call prin2('and funs2=*',funs2,n)
  call prin2('while funs=*',funs,n)
  call prin2('and diffs=*',diffs,n)

  !c
  !c       check that the matrices u,v are inverses of each other
  !c
  itype=2
  call hermexps(itype,n,xs,uu,vv,whts,w)


  call prin2('uu=*',uu,n*n)
  !cccc        call matmat(uu,vv,n,ww)

  call matmat(vv,uu,n,ww)


  call prin2('uu vv =*',ww,n*n)

  !c
  !c       construct a test function and decompose it
  !c
  do  i=1,n
     funs(i)=exp(-xs(i)**2/2) * xs(i)
     !cccc        funs(i)=1
  end do


  call prin2('funs as created*',funs,n)

  call matvec(uu,funs,coefs,n)
  call prin2('and coefs=*',coefs,n)


  !
  ! plot a specific hermite function
  !
  m = 9
  a = -5
  b = 5
  nplot = 200

  ifinit = 1
  do i = 1,nplot
     ts(i) = a + (b-a)/(nplot-1)*(i-1)
     call hermfuns(ts(i) , m+10, funs, ders, ifinit, w)
     vals(i) = funs(m+1)
  end do

  iw = 21
  itype = 2
  call pyplot(iw, ts, vals, nplot, itype, 'hermite polynomial*')


  stop
end program test_hermexps





subroutine matvec(a,x,y,n)
  implicit real *8 (a-h,o-z)
  dimension a(n,n)
  real *8 x(n),y(n),d

  do  i=1,n
     d=0
     do  j=1,n
        d=d+a(i,j)*x(j)
     end do

     y(i)=d
  end do

  return
end subroutine matvec





subroutine matmat(a,b,n,c)
  implicit real *8 (a-h,o-z)
  save
  dimension a(n,n),b(n,n),c(n,n)

  do  i=1,n
     do  j=1,n
        d=0
        do  k=1,n
           d=d+a(i,k)*b(k,j)
        end do

        c(i,j)=d
     end do
  end do

  return
end subroutine matmat
