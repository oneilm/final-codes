c***********************************************************************
      subroutine alpert_sqrt(n,norder,xnodes,whts)
c***********************************************************************
      implicit none
      integer n,i,k,ia,norder,inc,m,ii,iasqrt
      real *8 xnodes(n),whts(n)
      real *8 eps,done,pi,h,t,xk,delta,deltold,der,pol,sum
      real *8 xsqrt(16),wsqrt(16)
      real *8 x(16),w(16)
c
c     This subroutine constructs the nodes and the
c     weights of the n-point gauss-trapezoidal quadrature of 
c     order norder on [0,1] with square root singularity at 0.
c
c     input parameters:
c
c     n - the number of nodes in the quadrature n > 2*n order+1
c     norder - Alpert correction order
c
c     output parameters:
c
c     xnodes - the nodes of the n-point quadrature
c     whts   - the weights of the n-point quadrature
c
c
      if (norder.eq.4) then
         ia = 2
         x(1) = 2.000000000000000D-1 
         x(2) = 1.000000000000000D0 
         w(1) = 5.208333333333333D-1 
         w(2) = 9.791666666666667D-1
      else if (norder.eq.8) then
         ia = 4
         x(1) = 2.087647422032129D-1 
         x(2) = 9.786087373714483D-1 
         x(3) = 1.989541386579751D0 
         x(4) = 3.000000000000000D0
         w(1) = 5.207988277246498D-1 
         w(2) = 9.535038018555888D-1 
         w(3) = 1.024871626402471D0 
         w(4) = 1.000825744017291D0
      else if (norder.eq.16) then
         ia = 7
         x(1) = 9.919337841451028D-2 
         x(2) = 5.076592669645529D-1 
         x(3) = 1.184972925827278D0 
         x(4) = 2.047493467134072D0 
         x(5) = 3.007168911869310D0 
         x(6) = 4.000474996776184D0 
         x(7) = 5.000007879022339D0 
         x(8) = 6.000000000000000D0
         w(1) = 2.528198928766921D-1 
         w(2) = 5.550158230159486D-1 
         w(3) = 7.852321453615224D-1 
         w(4) = 9.245915673876714D-1 
         w(5) = 9.839350200445296D-1 
         w(6) = 9.984463448413151D-1 
         w(7) = 9.999592378464547D-1 
         w(8) = 9.999999686258662D-1
      else
         return
      endif
c
      if (norder.eq.4) then
         iasqrt = 3
         xsqrt(1) = 1.189242434021285D-02
         xsqrt(2) = 2.578220434738662D-01
         xsqrt(3) = 1.007750064585281D0
         xsqrt(4) = 2.000000000000000D0
         wsqrt(1) = 5.927215035616424D-02
         wsqrt(2) = 4.955981740306228D-01
         wsqrt(3) = 9.427131290628058D-01
         wsqrt(4) = 1.002416546550407D0
      else if (norder.eq.8) then
         iasqrt = 5
         xsqrt(1) = 1.214130606523435D-03
         xsqrt(2) = 3.223952700027058D-02
         xsqrt(3) = 1.790935383649920D-01
         xsqrt(4) = 5.437663805244631D-01
         xsqrt(5) = 1.176116628396759D0
         xsqrt(6) = 2.031848210716014D0
         xsqrt(7) = 3.001961225690812D0
         xsqrt(8) = 4.000000000000000D0
         wsqrt(1) = 6.199844884297793D-03
         wsqrt(2) = 7.106286791720044D-02
         wsqrt(3) = 2.408930104410471D-01
         wsqrt(4) = 4.975929263668960D-01
         wsqrt(5) = 7.592446540441226D-01
         wsqrt(6) = 9.322446399614420D-01
         wsqrt(7) = 9.928171438160095D-01
         wsqrt(8) = 9.999449125689846D-01
      else if (norder.eq.16) then
         iasqrt = 10
         xsqrt(1) = 2.158438988280793D-04
         xsqrt(2) = 5.898432743709196D-03
         xsqrt(3) = 3.462795956896131D-02
         xsqrt(4) = 1.145586495070213D-01
         xsqrt(5) = 2.790344218856415D-01
         xsqrt(6) = 5.600113798653321D-01
         xsqrt(7) = 9.814091242883119D-01
         xsqrt(8) = 1.553594853974655D0
         xsqrt(9) = 2.270179114036658D0
         xsqrt(10) = 3.108234601715371D0
         xsqrt(11) = 4.032930893996553D0
         xsqrt(12) = 5.006803270228157D0
         xsqrt(13) = 6.000815466735179D0
         xsqrt(14) = 7.000045035079542D0
         xsqrt(15) = 8.000000738923901D0
         xsqrt(16) = 9.000000000000000D0
         wsqrt(1) = 1.105804873501181D-03
         wsqrt(2) = 1.324499944707956D-02
         wsqrt(3) = 4.899842307592144D-02
         wsqrt(4) = 1.165326192868815D-01
         wsqrt(5) = 2.178586693194957D-01
         wsqrt(6) = 3.481766016945031D-01
         wsqrt(7) = 4.964027915911545D-01
         wsqrt(8) = 6.469026189623831D-01
         wsqrt(9) = 7.823688971783889D-01
         wsqrt(10) = 8.877772445893361D-01
         wsqrt(11) = 9.551665077035583D-01
         wsqrt(12) = 9.876285579741800D-01
         wsqrt(13) = 9.979929183863017D-01
         wsqrt(14) = 9.998470620634641D-01
         wsqrt(15) = 9.999962891645340D-01
         wsqrt(16) = 9.999999946893169D-01
      else
         return
      endif
c
      m = n-norder-norder/2
ccc      write(6,*) 'm=',m
      h = 1.0d0/(m+ia+iasqrt-1.0d0)
      do i = 1,norder
         whts(i) = wsqrt(i)*h
         xnodes(i) = xsqrt(i)*h
      enddo
      do i = 1,m
         whts(norder+i) = h
         xnodes(norder+i) = iasqrt*h+(i-1)*h
      enddo
c
      ii = norder+m
      inc = 1
      do i = norder/2,1,-1
         xnodes(ii+inc) = 1.0d0-x(i)*h
         whts(ii+inc) = w(i)*h
         inc = inc+1
      enddo
      return
      end
c
c
c
ccc***********************************************************************
      subroutine alpert_sqrta(a,b,n,norder,xnodes,whts)
c***********************************************************************
      implicit none
      integer n,i,k,ia,norder,inc,m,ii,iasqrt
      real *8 xnodes(n),whts(n),a,b
      real *8 eps,done,pi,h,t,xk,delta,deltold,der,pol,sum
      real *8 xsqrt(16),wsqrt(16)
      real *8 x(16),w(16)
c
c     This subroutine constructs the nodes and the
c     weights of the n-point gauss-trapezoidal quadrature of 
c     order norder on [a,b] with sqrt singularity at a.
c
c     input parameters:
c
c     a,b  the interval
c     n - the number of nodes in the quadrature n > 2*n order+1
c     norder - Alpert correction order
c
c     output parameters:
c
c     xnodes - the nodes of the n-point quadrature
c     whts   - the weights of the n-point quadrature
c
c
      if (norder.eq.4) then
         ia = 2
         x(1) = 2.000000000000000D-1 
         x(2) = 1.000000000000000D0 
         w(1) = 5.208333333333333D-1 
         w(2) = 9.791666666666667D-1
      else if (norder.eq.8) then
         ia = 4
         x(1) = 2.087647422032129D-1 
         x(2) = 9.786087373714483D-1 
         x(3) = 1.989541386579751D0 
         x(4) = 3.000000000000000D0
         w(1) = 5.207988277246498D-1 
         w(2) = 9.535038018555888D-1 
         w(3) = 1.024871626402471D0 
         w(4) = 1.000825744017291D0
      else if (norder.eq.16) then
         ia = 7
         x(1) = 9.919337841451028D-2 
         x(2) = 5.076592669645529D-1 
         x(3) = 1.184972925827278D0 
         x(4) = 2.047493467134072D0 
         x(5) = 3.007168911869310D0 
         x(6) = 4.000474996776184D0 
         x(7) = 5.000007879022339D0 
         x(8) = 6.000000000000000D0
         w(1) = 2.528198928766921D-1 
         w(2) = 5.550158230159486D-1 
         w(3) = 7.852321453615224D-1 
         w(4) = 9.245915673876714D-1 
         w(5) = 9.839350200445296D-1 
         w(6) = 9.984463448413151D-1 
         w(7) = 9.999592378464547D-1 
         w(8) = 9.999999686258662D-1
      else
         return
      endif
c
      if (norder.eq.4) then
         iasqrt = 3
         xsqrt(1) = 1.189242434021285D-02
         xsqrt(2) = 2.578220434738662D-01
         xsqrt(3) = 1.007750064585281D0
         xsqrt(4) = 2.000000000000000D0
         wsqrt(1) = 5.927215035616424D-02
         wsqrt(2) = 4.955981740306228D-01
         wsqrt(3) = 9.427131290628058D-01
         wsqrt(4) = 1.002416546550407D0
      else if (norder.eq.8) then
         iasqrt = 5
         xsqrt(1) = 1.214130606523435D-03
         xsqrt(2) = 3.223952700027058D-02
         xsqrt(3) = 1.790935383649920D-01
         xsqrt(4) = 5.437663805244631D-01
         xsqrt(5) = 1.176116628396759D0
         xsqrt(6) = 2.031848210716014D0
         xsqrt(7) = 3.001961225690812D0
         xsqrt(8) = 4.000000000000000D0
         wsqrt(1) = 6.199844884297793D-03
         wsqrt(2) = 7.106286791720044D-02
         wsqrt(3) = 2.408930104410471D-01
         wsqrt(4) = 4.975929263668960D-01
         wsqrt(5) = 7.592446540441226D-01
         wsqrt(6) = 9.322446399614420D-01
         wsqrt(7) = 9.928171438160095D-01
         wsqrt(8) = 9.999449125689846D-01
      else if (norder.eq.16) then
         iasqrt = 10
         xsqrt(1) = 2.158438988280793D-04
         xsqrt(2) = 5.898432743709196D-03
         xsqrt(3) = 3.462795956896131D-02
         xsqrt(4) = 1.145586495070213D-01
         xsqrt(5) = 2.790344218856415D-01
         xsqrt(6) = 5.600113798653321D-01
         xsqrt(7) = 9.814091242883119D-01
         xsqrt(8) = 1.553594853974655D0
         xsqrt(9) = 2.270179114036658D0
         xsqrt(10) = 3.108234601715371D0
         xsqrt(11) = 4.032930893996553D0
         xsqrt(12) = 5.006803270228157D0
         xsqrt(13) = 6.000815466735179D0
         xsqrt(14) = 7.000045035079542D0
         xsqrt(15) = 8.000000738923901D0
         xsqrt(16) = 9.000000000000000D0
         wsqrt(1) = 1.105804873501181D-03
         wsqrt(2) = 1.324499944707956D-02
         wsqrt(3) = 4.899842307592144D-02
         wsqrt(4) = 1.165326192868815D-01
         wsqrt(5) = 2.178586693194957D-01
         wsqrt(6) = 3.481766016945031D-01
         wsqrt(7) = 4.964027915911545D-01
         wsqrt(8) = 6.469026189623831D-01
         wsqrt(9) = 7.823688971783889D-01
         wsqrt(10) = 8.877772445893361D-01
         wsqrt(11) = 9.551665077035583D-01
         wsqrt(12) = 9.876285579741800D-01
         wsqrt(13) = 9.979929183863017D-01
         wsqrt(14) = 9.998470620634641D-01
         wsqrt(15) = 9.999962891645340D-01
         wsqrt(16) = 9.999999946893169D-01
      else
         return
      endif
c
      m = n-norder-norder/2
ccc      write(6,*) 'm=',m
cccc      h = 1.0d0/(m+ia+iasqrt-1.0d0)
      h = (b-a)/(m+ia+iasqrt-1.0d0)
      do i = 1,norder
         whts(i) = wsqrt(i)*h
         xnodes(i) = a + xsqrt(i)*h
      enddo
      do i = 1,m
         whts(norder+i) = h
         xnodes(norder+i) = a + iasqrt*h+(i-1)*h
      enddo
c
      ii = norder+m
      inc = 1
      do i = norder/2,1,-1
         xnodes(ii+inc) = b-x(i)*h
         whts(ii+inc) = w(i)*h
         inc = inc+1
      enddo
      return
      end
c
c
c
c
c
ccc***********************************************************************
      subroutine alpert_sqrtb(a,b,n,norder,xnodes,whts)
c***********************************************************************
      implicit none
      integer n,i,k,ia,norder,inc,m,ii,iasqrt
      real *8 xnodes(n),whts(n),a,b
      real *8 eps,done,pi,h,t,xk,delta,deltold,der,pol,sum
      real *8 xsqrt(16),wsqrt(16)
      real *8 x(16),w(16)
c
c     This subroutine constructs the nodes and the
c     weights of the n-point gauss-trapezoidal quadrature of 
c     order norder on [a,b] with sqrt singularity at b.
c
c     input parameters:
c
c     a,b  the interval
c     n - the number of nodes in the quadrature n > 2*n order+1
c     norder - Alpert correction order
c
c     output parameters:
c
c     xnodes - the nodes of the n-point quadrature
c     whts   - the weights of the n-point quadrature
c
c
      if (norder.eq.4) then
         ia = 2
         x(1) = 2.000000000000000D-1 
         x(2) = 1.000000000000000D0 
         w(1) = 5.208333333333333D-1 
         w(2) = 9.791666666666667D-1
      else if (norder.eq.8) then
         ia = 4
         x(1) = 2.087647422032129D-1 
         x(2) = 9.786087373714483D-1 
         x(3) = 1.989541386579751D0 
         x(4) = 3.000000000000000D0
         w(1) = 5.207988277246498D-1 
         w(2) = 9.535038018555888D-1 
         w(3) = 1.024871626402471D0 
         w(4) = 1.000825744017291D0
      else if (norder.eq.16) then
         ia = 7
         x(1) = 9.919337841451028D-2 
         x(2) = 5.076592669645529D-1 
         x(3) = 1.184972925827278D0 
         x(4) = 2.047493467134072D0 
         x(5) = 3.007168911869310D0 
         x(6) = 4.000474996776184D0 
         x(7) = 5.000007879022339D0 
         x(8) = 6.000000000000000D0
         w(1) = 2.528198928766921D-1 
         w(2) = 5.550158230159486D-1 
         w(3) = 7.852321453615224D-1 
         w(4) = 9.245915673876714D-1 
         w(5) = 9.839350200445296D-1 
         w(6) = 9.984463448413151D-1 
         w(7) = 9.999592378464547D-1 
         w(8) = 9.999999686258662D-1
      else
         return
      endif
c
      if (norder.eq.4) then
         iasqrt = 3
         xsqrt(1) = 1.189242434021285D-02
         xsqrt(2) = 2.578220434738662D-01
         xsqrt(3) = 1.007750064585281D0
         xsqrt(4) = 2.000000000000000D0
         wsqrt(1) = 5.927215035616424D-02
         wsqrt(2) = 4.955981740306228D-01
         wsqrt(3) = 9.427131290628058D-01
         wsqrt(4) = 1.002416546550407D0
      else if (norder.eq.8) then
         iasqrt = 5
         xsqrt(1) = 1.214130606523435D-03
         xsqrt(2) = 3.223952700027058D-02
         xsqrt(3) = 1.790935383649920D-01
         xsqrt(4) = 5.437663805244631D-01
         xsqrt(5) = 1.176116628396759D0
         xsqrt(6) = 2.031848210716014D0
         xsqrt(7) = 3.001961225690812D0
         xsqrt(8) = 4.000000000000000D0
         wsqrt(1) = 6.199844884297793D-03
         wsqrt(2) = 7.106286791720044D-02
         wsqrt(3) = 2.408930104410471D-01
         wsqrt(4) = 4.975929263668960D-01
         wsqrt(5) = 7.592446540441226D-01
         wsqrt(6) = 9.322446399614420D-01
         wsqrt(7) = 9.928171438160095D-01
         wsqrt(8) = 9.999449125689846D-01
      else if (norder.eq.16) then
         iasqrt = 10
         xsqrt(1) = 2.158438988280793D-04
         xsqrt(2) = 5.898432743709196D-03
         xsqrt(3) = 3.462795956896131D-02
         xsqrt(4) = 1.145586495070213D-01
         xsqrt(5) = 2.790344218856415D-01
         xsqrt(6) = 5.600113798653321D-01
         xsqrt(7) = 9.814091242883119D-01
         xsqrt(8) = 1.553594853974655D0
         xsqrt(9) = 2.270179114036658D0
         xsqrt(10) = 3.108234601715371D0
         xsqrt(11) = 4.032930893996553D0
         xsqrt(12) = 5.006803270228157D0
         xsqrt(13) = 6.000815466735179D0
         xsqrt(14) = 7.000045035079542D0
         xsqrt(15) = 8.000000738923901D0
         xsqrt(16) = 9.000000000000000D0
         wsqrt(1) = 1.105804873501181D-03
         wsqrt(2) = 1.324499944707956D-02
         wsqrt(3) = 4.899842307592144D-02
         wsqrt(4) = 1.165326192868815D-01
         wsqrt(5) = 2.178586693194957D-01
         wsqrt(6) = 3.481766016945031D-01
         wsqrt(7) = 4.964027915911545D-01
         wsqrt(8) = 6.469026189623831D-01
         wsqrt(9) = 7.823688971783889D-01
         wsqrt(10) = 8.877772445893361D-01
         wsqrt(11) = 9.551665077035583D-01
         wsqrt(12) = 9.876285579741800D-01
         wsqrt(13) = 9.979929183863017D-01
         wsqrt(14) = 9.998470620634641D-01
         wsqrt(15) = 9.999962891645340D-01
         wsqrt(16) = 9.999999946893169D-01
      else
         return
      endif
c
      m = n-norder-norder/2
ccc      write(6,*) 'm=',m
      h = 1.0d0/(m+ia+iasqrt-1.0d0)
      h = (b-a)/(m+ia+iasqrt-1.0d0)
      do i = 1,norder/2
         whts(i) = w(i)*h
         xnodes(i) = a+ x(i)*h
      enddo
      do i = 1,m
         whts(norder/2+i) = h
         xnodes(norder/2+i) = a + ia*h+(i-1)*h
      enddo
c
      ii = norder/2+m
      inc = 1
      do i = norder,1,-1
         xnodes(ii+inc) = b-xsqrt(i)*h
         whts(ii+inc) = wsqrt(i)*h
         inc = inc+1
      enddo
      return
      end
c
c
c
c
c
ccc***********************************************************************
      subroutine alpert_sqrtab(a,b,n,norder,xnodes,whts)
c***********************************************************************
      implicit none
      integer n,i,k,ia,norder,inc,m,ii,iasqrt
      real *8 xnodes(n),whts(n),a,b
      real *8 eps,done,pi,h,t,xk,delta,deltold,der,pol,sum
      real *8 xsqrt(16),wsqrt(16)
c
c     This subroutine constructs the nodes and the
c     weights of the n-point gauss-trapezoidal quadrature of 
c     order norder on [a,b] with sqrt singularity at a and b.
c
c     input parameters:
c
c     a,b  the interval
c     n - the number of nodes in the quadrature n > 2*n order+1
c     norder - Alpert correction order
c
c     output parameters:
c
c     xnodes - the nodes of the n-point quadrature
c     whts   - the weights of the n-point quadrature
c
c
      if (norder.eq.4) then
         iasqrt = 3
         xsqrt(1) = 1.189242434021285D-02
         xsqrt(2) = 2.578220434738662D-01
         xsqrt(3) = 1.007750064585281D0
         xsqrt(4) = 2.000000000000000D0
         wsqrt(1) = 5.927215035616424D-02
         wsqrt(2) = 4.955981740306228D-01
         wsqrt(3) = 9.427131290628058D-01
         wsqrt(4) = 1.002416546550407D0
      else if (norder.eq.8) then
         iasqrt = 5
         xsqrt(1) = 1.214130606523435D-03
         xsqrt(2) = 3.223952700027058D-02
         xsqrt(3) = 1.790935383649920D-01
         xsqrt(4) = 5.437663805244631D-01
         xsqrt(5) = 1.176116628396759D0
         xsqrt(6) = 2.031848210716014D0
         xsqrt(7) = 3.001961225690812D0
         xsqrt(8) = 4.000000000000000D0
         wsqrt(1) = 6.199844884297793D-03
         wsqrt(2) = 7.106286791720044D-02
         wsqrt(3) = 2.408930104410471D-01
         wsqrt(4) = 4.975929263668960D-01
         wsqrt(5) = 7.592446540441226D-01
         wsqrt(6) = 9.322446399614420D-01
         wsqrt(7) = 9.928171438160095D-01
         wsqrt(8) = 9.999449125689846D-01
      else if (norder.eq.16) then
         iasqrt = 10
         xsqrt(1) = 2.158438988280793D-04
         xsqrt(2) = 5.898432743709196D-03
         xsqrt(3) = 3.462795956896131D-02
         xsqrt(4) = 1.145586495070213D-01
         xsqrt(5) = 2.790344218856415D-01
         xsqrt(6) = 5.600113798653321D-01
         xsqrt(7) = 9.814091242883119D-01
         xsqrt(8) = 1.553594853974655D0
         xsqrt(9) = 2.270179114036658D0
         xsqrt(10) = 3.108234601715371D0
         xsqrt(11) = 4.032930893996553D0
         xsqrt(12) = 5.006803270228157D0
         xsqrt(13) = 6.000815466735179D0
         xsqrt(14) = 7.000045035079542D0
         xsqrt(15) = 8.000000738923901D0
         xsqrt(16) = 9.000000000000000D0
         wsqrt(1) = 1.105804873501181D-03
         wsqrt(2) = 1.324499944707956D-02
         wsqrt(3) = 4.899842307592144D-02
         wsqrt(4) = 1.165326192868815D-01
         wsqrt(5) = 2.178586693194957D-01
         wsqrt(6) = 3.481766016945031D-01
         wsqrt(7) = 4.964027915911545D-01
         wsqrt(8) = 6.469026189623831D-01
         wsqrt(9) = 7.823688971783889D-01
         wsqrt(10) = 8.877772445893361D-01
         wsqrt(11) = 9.551665077035583D-01
         wsqrt(12) = 9.876285579741800D-01
         wsqrt(13) = 9.979929183863017D-01
         wsqrt(14) = 9.998470620634641D-01
         wsqrt(15) = 9.999962891645340D-01
         wsqrt(16) = 9.999999946893169D-01
      else
         return
      endif
c
      m = n-norder-norder
ccc      write(6,*) 'm=',m
      h = (b-a)/(m+2*iasqrt-1.0d0)
      do i = 1,norder
         whts(i) = wsqrt(i)*h
         xnodes(i) = a + xsqrt(i)*h
      enddo
      do i = 1,m
         whts(norder+i) = h
         xnodes(norder+i) = a + iasqrt*h+(i-1)*h
      enddo
c
      ii = norder+m
      inc = 1
      do i = norder,1,-1
         xnodes(ii+inc) = b-xsqrt(i)*h
         whts(ii+inc) = wsqrt(i)*h
         inc = inc+1
      enddo
      return
      end
c
c
c
c
c
c
ccc***********************************************************************
      subroutine alpert_smooth(a,b,n,norder,xnodes,whts)
c***********************************************************************
      implicit none
      integer n,i,k,ia,norder,inc,m,ii,iasqrt
      real *8 xnodes(n),whts(n),a,b
      real *8 eps,done,pi,h,t,xk,delta,deltold,der,pol,sum
      real *8 x(16),w(16)
c
c     This subroutine constructs the nodes and the
c     weights of the n-point gauss-trapezoidal quadrature of 
c     order norder on [a,b] for smooth functions
c
c     input parameters:
c
c     a,b  the interval
c     n - the number of nodes in the quadrature n > 2*norder+1
c     norder - Alpert correction order
c
c     output parameters:
c
c     xnodes - the nodes of the n-point quadrature
c     whts   - the weights of the n-point quadrature
c
c
      if (norder.eq.4) then
         ia = 2
         x(1) = 2.000000000000000D-1 
         x(2) = 1.000000000000000D0 
         w(1) = 5.208333333333333D-1 
         w(2) = 9.791666666666667D-1
      else if (norder.eq.8) then
         ia = 4
         x(1) = 2.087647422032129D-1 
         x(2) = 9.786087373714483D-1 
         x(3) = 1.989541386579751D0 
         x(4) = 3.000000000000000D0
         w(1) = 5.207988277246498D-1 
         w(2) = 9.535038018555888D-1 
         w(3) = 1.024871626402471D0 
         w(4) = 1.000825744017291D0
      else if (norder.eq.16) then
         ia = 7
         x(1) = 9.919337841451028D-2 
         x(2) = 5.076592669645529D-1 
         x(3) = 1.184972925827278D0 
         x(4) = 2.047493467134072D0 
         x(5) = 3.007168911869310D0 
         x(6) = 4.000474996776184D0 
         x(7) = 5.000007879022339D0 
         x(8) = 6.000000000000000D0
         w(1) = 2.528198928766921D-1 
         w(2) = 5.550158230159486D-1 
         w(3) = 7.852321453615224D-1 
         w(4) = 9.245915673876714D-1 
         w(5) = 9.839350200445296D-1 
         w(6) = 9.984463448413151D-1 
         w(7) = 9.999592378464547D-1 
         w(8) = 9.999999686258662D-1
      else
         return
      endif
c
      m = n-norder/2-norder/2
ccc      write(6,*) 'm=',m
      h = (b-a)/(m+ia+ia-1.0d0)
      do i = 1,norder/2
         whts(i) = w(i)*h
         xnodes(i) = a+ x(i)*h
      enddo
      do i = 1,m
         whts(norder/2+i) = h
         xnodes(norder/2+i) = a + ia*h+(i-1)*h
      enddo
c
      ii = norder/2+m
      inc = 1
      do i = norder/2,1,-1
         xnodes(ii+inc) = b-x(i)*h
         whts(ii+inc) = w(i)*h
         inc = inc+1
      enddo
      return
      end
c
c

