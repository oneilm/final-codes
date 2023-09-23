c 
c 
c 
c 
c 
        subroutine o3rrdall(ir,xs,whts,
     1      points1,weights1,numspts,coefslrg1,
     2      points2,weights2,coefslrg2,npts2,nregions2,
     3      points3,weights3,coefslrg3,npts3,nregions3,
     4      points4,weights4,coefslrg4,npts4,nregions4,
     5      points5,weights5,coefslrg5,npts5,nregions5)
c 
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(1),whts(1)
        dimension numspts(20),
     1      points1(24,20),weights1(24,20),coefslrg1(20,24,20),
     2      points2(27,21),weights2(27,21),coefslrg2(20,27,21),
     3      points3(31,21),weights3(31,21),coefslrg3(20,31,21),
     4      points4(38,21),weights4(38,21),coefslrg4(20,38,21),
     5      points5(47,42),weights5(47,42),coefslrg5(20,47,42)
c 
c        This subroutine reads from the FORTRAN unit ir the data
c        to be used in the discretization of polygonal surfaces in
c        two dimensions (hopefully, these have been previously stored
c        there by a prior call to the subroutine o3rwrall), and returns
c        them to the user.
c 
c                     Input parameters:
c 
c  ir - the FORTRAN unit number from which the subroutine is to
c        read the data.
c 
c                     Output parameters:
c 
c  xs - the 20 support nodes on the interval [-1,1]
c  ws - the 20 quadrature weights corresponding to the support nodes xs
c  points1,weights1,coefslrg1 - nodes, quadrature weights, and
c        interpolation matrices corresponding to the interactions
c        of the chunk with itself
c  points2,weights2,coefslrg2 - nodes, quadrature weights, and
c        interpolation matrices corresponding to the interactions
c        of the chunk with nodes in the zone 1 with respect to the
c        chunk (theta > pi/6)
c  points3,weights3,coefslrg3 - nodes, quadrature weights, and
c        interpolation matrices corresponding to the interactions
c        of the chunk with nodes in the zone 2 with respect to the
c        chunk (theta > pi/12)
c  points4,weights4,coefslrg4 - nodes, quadrature weights, and
c        interpolation matrices corresponding to the interactions
c        of the chunk with nodes in the zone 3 with respect to the
c        chunk (theta > pi/24)
c  points5,weights5,coefslrg5 - nodes, quadrature weights, and
c        interpolation matrices corresponding to the interactions
c        of the chunk with nodes in the zone 4 with respect to the
c        chunk (theta > pi/36)
c 
        rewind(ir)
c 
c       read from disk the data to be used in the construction of
c       interactions of a chunk with itself
c 
        call o1rread(ir,xs,whts,points1,weights1,
     1      numspts,coefslrg1)
c 
c       read from disk the data to be used in the
c       construction of interactions between a chunk and its neighbors
c 
c       . . . 30-degree regime
c 
        call o3rread0(ir,npts2,nregions2,xs,whts,
     1      points2,weights2,coefslrg2)
c 
c       . . . 15-degree regime
c 
        call o3rread0(ir,npts3,nregions3,xs,whts,
     1      points3,weights3,coefslrg3)
c 
c       . . . 7.5-degree regime
c 
        call o3rread0(ir,npts4,nregions4,xs,whts,
     1      points4,weights4,coefslrg4)
c 
c       . . . 5-degree regime
c 
        call o3rread0(ir,npts5,nregions5,xs,whts,
     1      points5,weights5,coefslrg5)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o3rread0(ir,npts,nregions,xs,whts,
     1      points,weights,coefslrg)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension weights(1),points(1),coefslrg(1)
        real *8 xs(20),whts(20)
c 
        read(ir) npts
        read(ir) nregions
        read(ir) (xs(i),i=1,20)
        read(ir) (whts(i),i=1,20)
c 
        read(ir) (points(i),i=1,nregions*npts)
        read(ir) (weights(i),i=1,nregions*npts)
        read(ir) (coefslrg(i),i=1,nregions*npts*20)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o3rwrall(iw,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(21),whts(21),w(1)
c 
c       This subroutine generates and stores on disk as unformatted
c       FORTRAN records the precomputed data to be used in the
c       construction of discretizations of polygonal surfaces in
c       two dimensions. The data can be subsegently retrieved via the
c       subroutines o3rrdall, o1rread (see).
c 
c                    Input parameters:
c 
c  iw - the FORTRAN unit number on which the data are to be written
c 
c                    Output parameters: none
c 
c                    work arrays
c 
c  w - must be at least 44000 real *8 locations long
c 
c 
c       . . . construct the data to be used in the construction of
c             interactions of a chunk with itself, and store it on disk
c 
        call o1rwrite(iw)
c 
c       construct and store on disk the data to be used in the
c       construction of interactions between a chunk and its neighbors
c 
        ipoints=1
        lpoints=2000
c 
        iweights=ipoints+lpoints
        lweights=2000
c 
        icoefs=iweights+lweights
c 
c       . . . 30-degree regime
c 
        npts=27
        nregions=21
c 
        call o2rnear7(w(ipoints),w(iweights) )
c 
        call o3rmatbs(w(icoefs),xs,whts,npts,nregions,
     1      w(ipoints),w(iweights))
c 
        call o3rwrt0(iw,npts,nregions,xs,whts,
     1      w(ipoints),w(iweights),w(icoefs))
c 
c       . . . 15-degree regime
c 
        npts=31
        nregions=21
c 
        call o3rnear7(w(ipoints),w(iweights))
c 
        call o3rmatbs(w(icoefs),xs,whts,npts,nregions,
     1      w(ipoints),w(iweights))
c 
        call o3rwrt0(iw,npts,nregions,xs,whts,
     1      w(ipoints),w(iweights),w(icoefs))
c 
c       . . . 7.5-degree regime
c 
        npts=38
        nregions=21
c 
        call o4rnear7(w(ipoints),w(iweights))
c 
        call o3rmatbs(w(icoefs),xs,whts,npts,nregions,
     1      w(ipoints),w(iweights))
c 
        call o3rwrt0(iw,npts,nregions,xs,whts,
     1      w(ipoints),w(iweights),w(icoefs))
c 
c       . . . 5-degree regime
c 
        npts=47
        nregions=42
c 
        call o5rnear7(w(ipoints),w(iweights))
c 
        call o3rmatbs(w(icoefs),xs,whts,npts,nregions,
     1      w(ipoints),w(iweights))
c 
        call o3rwrt0(iw,npts,nregions,xs,whts,
     1      w(ipoints),w(iweights),w(icoefs))
c 
c       construct and store on the same disk unit
c       the data for the "universal" diecretizations of
c       intervals
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine o3rwrt0(iw,npts,nregions,xs,whts,
     1      points,weights,coefslrg)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension weights(npts*nregions),points(npts*nregions),
     1      coefslrg(20*npts*nregions)
c 
        real *8 xs(20),whts(20)
c 
  
        call prinf('in o3rwrt0, npts=*',npts,1)
        call prinf('in o3rwrt0, nregions=*',nregions,1)
  
        write(iw) npts
        write(iw) nregions
        write(iw) xs
        write(iw) whts
        write(iw) (points(i),i=1,npts*nregions)
        write(iw) weights
        write(iw) coefslrg
        return
        end
c 
c 
c 
c 
c 
        subroutine o3rmatbs(coefslrg,xs,whts,npts,nregions,
     1      points,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension points(npts,nregions),weights(npts,nregions),
     1      xs(20),whts(20),dumpts(100),dumwhts(100),
     2      coefslrg(20,npts,nregions),forminte(100)
c 
c        . . . retrieve from the subroutine o1rnodes the support
c              nodes and corresponding weights
c 
        call o1rnodes(xs,whts,nn)
c 
c       for each of the 20 sets of quadrature nodes, construct the matrix
c       of interpolation coefficients connecting the values of a function
c       at the support nodes to its values at the quadrature nodes
c       PLEASE NOTE THAT DIFFERENT quadrature formulae have different
c       sets of nodes
c 
        nnn=nregions
        do 2000 k=1,nnn
c 
        numpts=npts
c 
        do 1800 i=1,numpts
c 
        x=points(i,k)
c 
        call o1rinter(x,forminte,dumpts,dumwhts,nn)
c 
        do 1600 j=1,nn
        coefslrg(j,i,k)=forminte(j)
 1600 continue
c 
 1800 continue
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ifinside(theta,x1,y1,x2,y2,x,y,ifin)
        implicit real *8 (a-h,o-z)
c 
c       Given an interval in R^2 and a point outside that interval,
c       this subroutine determines wheter the point is located
c       within a parallelogram of which the interval is the diagonal,
c       and the angles adjacent to that diagonal are 2*theta:
c 
c 
c 
c                                 *
c                                * *
c                               *   *
c                              *     *
c                             *       *
c                            *         *
c                           *           *
c                          *             *
c               (x1,y1)   *  theta        *  (x2,y2)
c                        *-----------------*
c                         *               *
c                          *             *
c                           *           *
c                            *         *
c                             *       *     * (x,y)
c                              *     *
c                               *   *
c                                * *
c                                 *
c 
c                      Input parameters:
c 
c  theta - half of the angle in the parallelogram adjacent to the
c          user-specified interval
c  (x1,y1) - the coordinates of one end of the interval
c  (x2,y2) - the coordinates of second end of the interval
c  (x,y) - the coordinates of the point whose location with respect
c          to the interval ((x1,y1),(x2,y2)) we are checking
c 
c                      Output parametes:
c 
c  ifin - the onteger parameter telling the used whether the point is
c          inside or outside the parallelogram:
c     ifin=0 means that the point is outside the parallelogram
c     ifin=1 means that the point is inside the parallelogram
c 
c       . . . scale
c 
        save
        d=sqrt( (x2-x1)**2+(y2-y1)**2)
        scale=2/d
c 
        xx1=x1*scale
        xx2=x2*scale
c 
        yy1=y1*scale
        yy2=y2*scale
c 
        xx=x*scale
        yy=y*scale
c 
c        . . . shift
c 
        xmiddle=(xx1+xx2)/2
        ymiddle=(yy1+yy2)/2
c 
        xx=xx-xmiddle
        yy=yy-ymiddle
c 
        xx1=xx1-xmiddle
        yy1=yy1-ymiddle
c 
        xx2=xx2-xmiddle
        yy2=yy2-ymiddle
c 
c       find the rotation that will convert the point(xx1,yy1)
c       into the point (1,0)
c 
        a11=xx1
        a12=yy1
        a21=-yy1
        a22=xx1
c 
        xx3=a11*xx1+a12*yy1
        yy3=a21*xx1+a22*yy1
c 
        xxx=a11*xx+a12*yy
        yyy=a21*xx+a22*yy
c 
        call ifinsid0(xxx,yyy,theta,ifin)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ifinsid0(x,y,theta,ifin)
        implicit real *8 (a-h,o-z)
c 
c       construct the equations of the four lines
c 
        save
        done=1
        pi=atan(done)*4
c 
        alpha1=tan(theta)
        beta1=alpha1
c 
        alpha2=-alpha1
        beta2=alpha2
c 
        alpha3=tan(pi-theta)
        beta3=-alpha3
c 
        alpha4=-alpha3
        beta4=-alpha4
c 
        d1=alpha1*x+beta1-y
        d2=alpha2*x+beta2-y
        d3=alpha3*x+beta3-y
        d4=alpha4*x+beta4-y
c 
        ifin=1
        if(d1 .lt.0) ifin=0
        if(d2 .gt.0) ifin=0
        if(d3 .lt.0) ifin=0
        if(d4 .gt.0) ifin=0
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine o1rauwrt(iw,n,w)
        implicit real *8 (a-h,o-z)
        save
        dimension points(21),weights(21),w(1)
c 
c       this subroutine createds and stores on disk on the FORTRAN
c       unit iw the nodes and weights of a "universal" quadrature,
c       and the matrix interpolating "appropriately behaved" functions
c       from the 20 support nodes to the n "universal" nodes. The data
c       stored by this subroutine can be retrieved from disk via the
c       subroutine o1raurd (see)
c 
c                        Input parameters:
c 
c  iw - the FORTRAN unit number on which the subroutine will store
c       the data it generates
c  n - the number of nodes in the universal quadrature
c 
c                        Output parameters: none
c 
c 
c                        Work arrays:
c 
c  w - must be at least n*22+100 real *8 locations long
c 
c 
c       . . . construct all the appropriate parameters
c 
        ndigits=15
c 
        ixs=1
        lxs=n+10
c 
        iwhts=ixs+lxs
        lwhts=n+10
c 
        iamatr=iwhts+lwhts
  
        call o1raupnt(ndigits,n,w(ixs),w(iwhts),w(iamatr),
     1      points,weights)
c 
c       store all of these parameters on disk
c 
        write(iw) n
        write(iw) (w(ixs+i-1),i=1,n)
        write(iw) (w(iwhts+i-1),i=1,n)
        write(iw) (w(iamatr+i-1),i=1,n*20)
c 
        write(iw) (points(i),i=1,20)
        write(iw) (weights(i),i=1,20)
        return
        end
c 
c 
c 
c 
c 
        subroutine o1raupnt(ndigits,npts,xs,whts,amatr,
     1      points,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),points(21),
     1      weights(21),amatr(20,npts)
c 
c       this subroutine returns to the user the nodes and weights
c       of a "universal" quadrature, and the matrix interpolating
c       "appropriately behaved" functions from the 20 support nodes
c       to the npts "universal" nodes. For convenience, it also
c       returns the standard support nodes and weights
c 
c                        Input parameters:
c 
c  ndigits - the accuracy (number of digits) with which the "universal"
c       quadrature is constructed
c  npts - the number of nodes in the universal quadrature
c 
c                        Output parameters:
c 
c  xs - the nodes of the universal quadrature
c  whts - the weights of the universal quadrature
c  amatr - the matrix whose adjoint converts the values of an
c       "appropriately behaved" function at the support nodes into its
c       values at the universal nodes
c  points - the support nodes (20 of them things)
c  weights - the weights corresponding to the 20 support nodes
c 
c 
c 
c        . . . construct the "universal" quadrature on the
c              interval [-1,1]
c 
        itype=1
  
        call probexps(ndigits,itype,npts,xs,u,v,whts)
  
  
        call prin2('after probexps, xs=*',xs,npts)
        call prin2('after probexps, whts=*',whts,npts)
  
c 
c       construct the interpolation matrix, converting the values of a
c       function at the support nodes into its values at the "universal:
c       quadrature nodes
c 
c 
c       one point after another, construct the matrix interpolating
c       functions from the 20 support nodes to the n "universal" nodes
c 
        do 1200 i=1,npts
c 
        call o1rinter(xs(i),amatr(1,i),points,weights,nn)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1raurd(ir,n,xs,whts,amatr,points,weights)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(1),whts(1),amatr(1),points(1),weights(1)
  
c 
c       this subroutine raed from the FORTRAN unit ir and returns
c       to the user the nodes and weights of a "universal" quadrature,
c       and the matrix interpolating "appropriately behaved" functions
c       from the 20 support nodes to the npts "universal" nodes.
c       For convenience, it also returns the standard support nodes and
c       weights. Hopefully, these data have been stored on disk via a
c       preceding call to the subroutine o1rauwrt (see)
c 
c 
c                        Input parameters:
c 
c  ir - the FORTRAN unit number from which the data are to be read
c 
c                        Output parameters:
c 
c  npts - the number of nodes in the universal quadrature
c  xs - the nodes of the universal quadrature
c  whts - the weights of the universal quadrature
c  amatr - the matrix whose adjoint converts the values of an
c       "appropriately behaved" function at the support nodes into its
c       values at the universal nodes
c  points - the support nodes (20 of them things)
c  weights - the weights corresponding to the 20 support nodes
c 
c 
c       . . . read all parameters from disk
c 
        read(ir) n
        read(ir) (xs(i),i=1,n)
        read(ir) (whts(i),i=1,n)
        read(ir) (amatr(i),i=1,n*20)
        read(ir) (points(i),i=1,20)
        read(ir) (weights(i),i=1,20)
c 
        return
        end
        subroutine o2rnear7(xsout,wsout)
        implicit real *8 (a-h,o-z)
        save
        dimension xsout(27,21),wsout(27,21)
c 
        dimension xs1(27),xs2(27),xs3(27),xs4(27),xs5(27),
     1      xs6(27),xs7(27),xs8(27),xs9(27),xs10(27),xs11(27),
     2      xs12(27),xs13(27),xs14(27),xs15(27),xs16(27),xs17(27),
     3      xs18(27),xs19(27),xs20(27),xs21(27),xz(27,21)
c 
        equivalence (xz(1,1),xs1(1)),(xz(1,2),xs2(1)),
     1      (xz(1,3),xs3(1)),(xz(1,4),xs4(1)),(xz(1,5),xs5(1)),
     2      (xz(1,6),xs6(1)),(xz(1,7),xs7(1)),(xz(1,8),xs8(1)),
     3      (xz(1,9),xs9(1)),(xz(1,10),xs10(1)),(xz(1,11),xs11(1))
  
        equivalence (xz(1,12),xs12(1)),
     1      (xz(1,13),xs13(1)),(xz(1,14),xs14(1)),(xz(1,15),xs15(1)),
     2      (xz(1,16),xs16(1)),(xz(1,17),xs17(1)),(xz(1,18),xs18(1)),
     3       (xz(1,19),xs19(1)),(xz(1,20),xs20(1)),(xz(1,21),xs21(1))
        dimension ws1(27),ws2(27),ws3(27),ws4(27),ws5(27),
     1      ws6(27),ws7(27),ws8(27),ws9(27),ws10(27),ws11(27),
     2      ws12(27),ws13(27),ws14(27),ws15(27),ws16(27),ws17(27),
     3      ws18(27),ws19(27),ws20(27),ws21(27),wz(27,21)
c 
        equivalence (wz(1,1),ws1(1)),(wz(1,2),ws2(1)),
     1      (wz(1,3),ws3(1)),(wz(1,4),ws4(1)),(wz(1,5),ws5(1)),
     2      (wz(1,6),ws6(1)),(wz(1,7),ws7(1)),(wz(1,8),ws8(1)),
     3      (wz(1,9),ws9(1)),(wz(1,10),ws10(1)),(wz(1,11),ws11(1))
  
        equivalence (wz(1,12),ws12(1)),
     1      (wz(1,13),ws13(1)),(wz(1,14),ws14(1)),(wz(1,15),ws15(1)),
     2      (wz(1,16),ws16(1)),(wz(1,17),ws17(1)),(wz(1,18),ws18(1)),
     3       (wz(1,19),ws19(1)),(wz(1,20),ws20(1)),(wz(1,21),ws21(1))
c 
        data xs1/
     1  -.9936675873447934E+00,-.9579391696926212E+00,
     2  -.8787060482145128E+00,-.7523150143005424E+00,
     3  -.5822976707784945E+00,-.3772696348217959E+00,
     4  -.1493615070749633E+00,0.8720532135734249E-01,
     5  0.3175242173527475E+00,0.5274040737016059E+00,
     6  0.7047411107676309E+00,0.8408743401809580E+00,
     7  0.9319809378780863E+00,0.9806564757234544E+00,
     8  0.9975503197911038E+00,0.9999012294131495E+00,
     9  0.9999875241314127E+00,0.9999950962438275E+00,
     *  0.9999970389760744E+00,0.9999978519236261E+00,
     *  0.9999983252846016E+00,0.9999986674018415E+00,
     *  0.9999989474392523E+00,0.9999992041688313E+00,
     *  0.9999994679538533E+00,0.9999997362498980E+00,
     *  0.9999999429563582E+00/
c 
        data xs2/
     1  -.9937480863543350E+00,-.9583368582955133E+00,
     2  -.8796105060565540E+00,-.7538378000099827E+00,
     3  -.5844726817812615E+00,-.3800541580130333E+00,
     4  -.1526413372285159E+00,0.8360351021750301E-01,
     5  0.3138148715941980E+00,0.5238206395888565E+00,
     6  0.7015138354278181E+00,0.8382095809757499E+00,
     7  0.9300394056112622E+00,0.9795166986090926E+00,
     8  0.9971248197819675E+00,0.9998353923182610E+00,
     9  0.9999758798550020E+00,0.9999903116137667E+00,
     *  0.9999941137475270E+00,0.9999957199473368E+00,
     *  0.9999966599743464E+00,0.9999973410814114E+00,
     *  0.9999978993863798E+00,0.9999984118282603E+00,
     *  0.9999989386179682E+00,0.9999994741174507E+00,
     *  0.9999998863109196E+00/
c 
        data xs3/
     1  -.9938417201259711E+00,-.9588061306444388E+00,
     2  -.8806895700971652E+00,-.7556696255401570E+00,
     3  -.5871057930598829E+00,-.3834423364804888E+00,
     4  -.1566496763217627E+00,0.7918343784603277E-01,
     5  0.3092430843606839E+00,0.5193823313928765E+00,
     6  0.6974932417780841E+00,0.8348655556413591E+00,
     7  0.9275787106383162E+00,0.9780473083445171E+00,
     8  0.9965492069166513E+00,0.9997278905701150E+00,
     9  0.9999537395512257E+00,0.9999809250474558E+00,
     *  0.9999883200251802E+00,0.9999914819010307E+00,
     *  0.9999933445927980E+00,0.9999946986406973E+00,
     *  0.9999958106181210E+00,0.9999968328360285E+00,
     *  0.9999978843611714E+00,0.9999989524828787E+00,
     *  0.9999997736665676E+00/
c 
        data xs4/
     1  -.9939506513934541E+00,-.9593609182505075E+00,
     2  -.8819802849259005E+00,-.7578791904545550E+00,
     3  -.5903015893070052E+00,-.3875742957309607E+00,
     4  -.1615578006841628E+00,0.7375028459063282E-01,
     5  0.3036004370857310E+00,0.5138787190550227E+00,
     6  0.6924795334189457E+00,0.8306660219211957E+00,
     7  0.9244578187471630E+00,0.9761498927214323E+00,
     8  0.9957665530975865E+00,0.9995535326920675E+00,
     9  0.9999120740431311E+00,0.9999625984515823E+00,
     *  0.9999768741624719E+00,0.9999830710964877E+00,
     *  0.9999867524179493E+00,0.9999894396398015E+00,
     *  0.9999916517272482E+00,0.9999936893274062E+00,
     *  0.9999957870797118E+00,0.9999979158962820E+00,
     *  0.9999995500183780E+00/
c 
        data xs5/
     1  -.9940778197142884E+00,-.9600200082694887E+00,
     2  -.8835321929208818E+00,-.7605576401040869E+00,
     3  -.5941978347576755E+00,-.3926335622715584E+00,
     4  -.1675891547716756E+00,0.6705042706171138E-01,
     5  0.2966157481194387E+00,0.5070358435350852E+00,
     6  0.6862122222607134E+00,0.8253805189917576E+00,
     7  0.9204908555714628E+00,0.9736917914919161E+00,
     8  0.9946957195920907E+00,0.9992724278534606E+00,
     9  0.9998345276943972E+00,0.9999270137534617E+00,
     *  0.9999543296019400E+00,0.9999664115460756E+00,
     *  0.9999736644657378E+00,0.9999789865001600E+00,
     *  0.9999833808059710E+00,0.9999874387810299E+00,
     *  0.9999916208279636E+00,0.9999958595473486E+00,
     *  0.9999991068476946E+00/
c 
        data xs6/
     1  -.9942271297726186E+00,-.9608081863664124E+00,
     2  -.8854100931892220E+00,-.7638231862509637E+00,
     3  -.5989719258232771E+00,-.3988552798600757E+00,
     4  -.1750293243083475E+00,0.5875972316626088E-01,
     5  0.2879418079445786E+00,0.4985019077585857E+00,
     6  0.6783555537884438E+00,0.8187103305799457E+00,
     7  0.9154340947343600E+00,0.9704941878615636E+00,
     8  0.9932208530514335E+00,0.9988217306608733E+00,
     9  0.9996918872714204E+00,0.9998583584374410E+00,
     *  0.9999100830166257E+00,0.9999334921748304E+00,
     *  0.9999477283716058E+00,0.9999582422883488E+00,
     *  0.9999669564815923E+00,0.9999750294398806E+00,
     *  0.9999833595853695E+00,0.9999917891316414E+00,
     *  0.9999982308592028E+00/
c 
        data xs7/
     1  -.9944036093369105E+00,-.9617568371436086E+00,
     2  -.8876949351260283E+00,-.7678219344462600E+00,
     3  -.6048415029147307E+00,-.4065269599571687E+00,
     4  -.1842277042767396E+00,0.4848007776643626E-01,
     5  0.2771496770988373E+00,0.4878393532750176E+00,
     6  0.6684888074725268E+00,0.8102774230268988E+00,
     7  0.9089735646125543E+00,0.9663183141652114E+00,
     8  0.9911772171863142E+00,0.9981031569437915E+00,
     9  0.9994327515619237E+00,0.9997268800457857E+00,
     *  0.9998236229520054E+00,0.9998686392253178E+00,
     *  0.9998964546577733E+00,0.9999171609971325E+00,
     *  0.9999344047125986E+00,0.9999504436632953E+00,
     *  0.9999670182534171E+00,0.9999837556566302E+00,
     *  0.9999965051324096E+00/
c 
        data xs8/
     1  -.9946103397027501E+00,-.9628872991292584E+00,
     2  -.8904438273199854E+00,-.7726589139470103E+00,
     3  -.6119659513778544E+00,-.4158647866011870E+00,
     4  -.1954565764397797E+00,0.3588764699668912E-01,
     5  0.2638737493907353E+00,0.4746568041245836E+00,
     6  0.6562160228278235E+00,0.7997042614127736E+00,
     7  0.9007702542579741E+00,0.9608776622088337E+00,
     8  0.9883414123089074E+00,0.9969662137338737E+00,
     9  0.9989684337167957E+00,0.9994773567252184E+00,
     *  0.9996556237088114E+00,0.9997413888328603E+00,
     *  0.9997954219771994E+00,0.9998360385207212E+00,
     *  0.9998700658031983E+00,0.9999018762317160E+00,
     *  0.9999348068134677E+00,0.9999679678269446E+00,
     *  0.9999931218595370E+00/
c 
        data xs9/
     1  -.9948344827996033E+00,-.9641315621151052E+00,
     2  -.8934957518020904E+00,-.7780596612247745E+00,
     3  -.6199582448878692E+00,-.4263912338482565E+00,
     4  -.2081864065737689E+00,0.2151640088790496E-01,
     5  0.2486033657277205E+00,0.4593542445397187E+00,
     6  0.6418110516371669E+00,0.7871129212100563E+00,
     7  0.8907865399289889E+00,0.9539983830670073E+00,
     8  0.9844734184985509E+00,0.9951976925191006E+00,
     9  0.9981515211578229E+00,0.9990095237689651E+00,
     *  0.9993317495959251E+00,0.9994931597452967E+00,
     *  0.9995973070204403E+00,0.9996765424579261E+00,
     *  0.9997434339589781E+00,0.9998063760623923E+00,
     *  0.9998716664794489E+00,0.9999371529741305E+00,
     *  0.9999865411940365E+00/
c 
        data xs10/
     1  -.9950387498651165E+00,-.9652753355984658E+00,
     2  -.8963176437405627E+00,-.7830800622883151E+00,
     3  -.6274324058215504E+00,-.4363073388775807E+00,
     4  -.2202852266103240E+00,0.7710793920593151E-02,
     5  0.2337450629948294E+00,0.4442302880433200E+00,
     6  0.6272874207915775E+00,0.7740667186404322E+00,
     7  0.8800154184225299E+00,0.9460822090727370E+00,
     8  0.9795242287462884E+00,0.9925648908755760E+00,
     9  0.9967584535284484E+00,0.9981492085994901E+00,
     *  0.9987152158508481E+00,0.9990135552424073E+00,
     *  0.9992119858581683E+00,0.9993653007659720E+00,
     *  0.9994960681612590E+00,0.9996201939923774E+00,
     *  0.9997492267888892E+00,0.9998778235606447E+00,
     *  0.9999739448394174E+00/
c 
        data xs11/
     1  -.9952445696673776E+00,-.9664236247032639E+00,
     2  -.8991237429070493E+00,-.7880116676629694E+00,
     3  -.6346769086245221E+00,-.4457874922124787E+00,
     4  -.2316934164307702E+00,-.5129375024510718E-02,
     5  0.2201089237501024E+00,0.4305176382229683E+00,
     6  0.6142335488405948E+00,0.7623431520057908E+00,
     7  0.8701471316060649E+00,0.9383847708911990E+00,
     8  0.9740672428101651E+00,0.9890649584388427E+00,
     9  0.9945560720448980E+00,0.9966398569653154E+00,
     *  0.9975784095287468E+00,0.9981100896607791E+00,
     *  0.9984793758027959E+00,0.9987712060772300E+00,
     *  0.9990241866970562E+00,0.9992676256401261E+00,
     *  0.9995209222649994E+00,0.9997694506064134E+00,
     *  0.9999513321602580E+00/
c 
        data xs12/
     1  -.9955743674046716E+00,-.9682926129732662E+00,
     2  -.9037027666609707E+00,-.7960269710998789E+00,
     3  -.6463649453743607E+00,-.4609370446563985E+00,
     4  -.2497118742925286E+00,-.2511926107877160E-01,
     5  0.1992590634813447E+00,0.4100229077630790E+00,
     6  0.5952787006097564E+00,0.7459280558719728E+00,
     7  0.8569255143101806E+00,0.9285236420602720E+00,
     8  0.9671881027392499E+00,0.9843857418382909E+00,
     9  0.9912562115501382E+00,0.9941523823227671E+00,
     *  0.9956046650782726E+00,0.9965072951054776E+00,
     *  0.9971722135916283E+00,0.9977164139024515E+00,
     *  0.9982042227314976E+00,0.9986883420838673E+00,
     *  0.9991887441888501E+00,0.9996449118631456E+00,
     *  0.9999333078702468E+00/
c 
        data xs13/
     1  -.9960692012515575E+00,-.9711645749892504E+00,
     2  -.9108511011120860E+00,-.8087057366788244E+00,
     3  -.6650930377153186E+00,-.4855363059698625E+00,
     4  -.2793775530052763E+00,-.5851908414295226E-01,
     5  0.1638481855640006E+00,0.3745413631514012E+00,
     6  0.5616843129860956E+00,0.7159875149539427E+00,
     7  0.8320070017052779E+00,0.9093686023955363E+00,
     8  0.9536223324320205E+00,0.9751887931897412E+00,
     9  0.9848055232798660E+00,0.9892908528623079E+00,
     *  0.9917406093226608E+00,0.9933612797293128E+00,
     *  0.9945983438076864E+00,0.9956339193311061E+00,
     *  0.9965829454612292E+00,0.9975405416761869E+00,
     *  0.9985260919250667E+00,0.9993970102302254E+00,
     *  0.9999019436259080E+00/
c 
        data xs14/
     1  -.9964491273766515E+00,-.9734726819292928E+00,
     2  -.9168250922221183E+00,-.8196715766546383E+00,
     3  -.6818061998941704E+00,-.5081490507406851E+00,
     4  -.3074560166884501E+00,-.9110439832577133E-01,
     5  0.1281387294112113E+00,0.3373844178478254E+00,
     6  0.5249067328895469E+00,0.6814181187034918E+00,
     7  0.8013133068248519E+00,0.8838247668879993E+00,
     8  0.9337205667225457E+00,0.9602791703546876E+00,
     9  0.9734976287457924E+00,0.9803702049644476E+00,
     *  0.9844874731224713E+00,0.9873870993457363E+00,
     *  0.9896782874996395E+00,0.9916413285176082E+00,
     *  0.9934780705923274E+00,0.9953498120938793E+00,
     *  0.9972550408619134E+00,0.9988977406656736E+00,
     *  0.9998250630804353E+00/
c 
        data xs15/
     1  -.9967322371174239E+00,-.9752683938876805E+00,
     2  -.9216474206973673E+00,-.8288057037837683E+00,
     3  -.6961279367782469E+00,-.5280640935877379E+00,
     4  -.3328900953862611E+00,-.1215305888768491E+00,
     5  0.9364405937216325E-01,0.3000669741011910E+00,
     6  0.4862513289946129E+00,0.6430636319813419E+00,
     7  0.7649563951272356E+00,0.8510535786697617E+00,
     8  0.9056500431387455E+00,0.9370775198056764E+00,
     9  0.9544650750111894E+00,0.9646247962961050E+00,
     *  0.9713674628025086E+00,0.9764436894483952E+00,
     *  0.9806124804437018E+00,0.9842920254061287E+00,
     *  0.9878257616451347E+00,0.9914515727317958E+00,
     *  0.9950606756228779E+00,0.9980573497044272E+00,
     *  0.9996956868034630E+00/
c 
        data xs16/
     1  -.9969587424228807E+00,-.9767485780983165E+00,
     2  -.9257269066614786E+00,-.8367086761221417E+00,
     3  -.7087812503789401E+00,-.5460300337981974E+00,
     4  -.3563458134915081E+00,-.1502790599954564E+00,
     5  0.6013859321144421E-01,0.2626256895604657E+00,
     6  0.4459192457903167E+00,0.6010558204517462E+00,
     7  0.7226336315855624E+00,0.8099155519180969E+00,
     8  0.8672067471199461E+00,0.9024676573623520E+00,
     9  0.9241700990245851E+00,0.9385997981200804E+00,
     *  0.9492897124905786E+00,0.9579209905436394E+00,
     *  0.9653563114547241E+00,0.9722307649118482E+00,
     *  0.9790792415448408E+00,0.9860448460018345E+00,
     *  0.9925126321090809E+00,0.9973012775795365E+00,
     *  0.9996081175627924E+00/
c 
        data xs17/
     1  -.9973224274243953E+00,-.9792059998487070E+00,
     2  -.9326841656675401E+00,-.8504577334814779E+00,
     3  -.7311232777462858E+00,-.5781070273427076E+00,
     4  -.3985672444118445E+00,-.2023101683556153E+00,
     5  -.6629588294904121E-03,0.1947101594080434E+00,
     6  0.3730341293178838E+00,0.5256741551655633E+00,
     7  0.6473848860230114E+00,0.7373906742862015E+00,
     8  0.7996749171027893E+00,0.8415225161949693E+00,
     9  0.8705317904254974E+00,0.8922121694464125E+00,
     *  0.9096627688893979E+00,0.9245293175476259E+00,
     *  0.9379219650429604E+00,0.9508195172578251E+00,
     *  0.9638663737498636E+00,0.9767748024678702E+00,
     *  0.9880938385941140E+00,0.9959299577312029E+00,
     *  0.9994414672913101E+00/
c 
        data xs18/
     1  -.9976955425554720E+00,-.9818145097568824E+00,
     2  -.9402924827897701E+00,-.8658765565669158E+00,
     3  -.7567480596214627E+00,-.6156912463551070E+00,
     4  -.4491079192272968E+00,-.2660088506889063E+00,
     5  -.7695239836011980E-01,0.1070873525169247E+00,
     6  0.2759337186743109E+00,0.4215228381523881E+00,
     7  0.5392910633609729E+00,0.6292737001578037E+00,
     8  0.6959284706845703E+00,0.7459262138004989E+00,
     9  0.7852219612175662E+00,0.8177323333307700E+00,
     *  0.8458053929800970E+00,0.8711272827776244E+00,
     *  0.8952188982796310E+00,0.9192135861986384E+00,
     *  0.9430795116528019E+00,0.9651995618577827E+00,
     *  0.9830131675806641E+00,0.9944194313144825E+00,
     *  0.9992538122636878E+00/
c 
        data xs19/
     1  -.9981423850478836E+00,-.9850400056700046E+00,
     2  -.9499662434084777E+00,-.8859335686365309E+00,
     3  -.7907428599444446E+00,-.6664733938898874E+00,
     4  -.5186645899929088E+00,-.3554089165278825E+00,
     5  -.1863425348079152E+00,-.2147687219218508E-01,
     6  0.1301989497793472E+00,0.2623285186694337E+00,
     7  0.3725870494378156E+00,0.4628200448936106E+00,
     8  0.5372807466196510E+00,0.6002674124228062E+00,
     9  0.6550137171745098E+00,0.7038994410418313E+00,
     *  0.7490175216299972E+00,0.7923844978224759E+00,
     *  0.8354293195432267E+00,0.8780458107365321E+00,
     *  0.9181963816439309E+00,0.9525979373498630E+00,
     *  0.9780994498585247E+00,0.9931817777049470E+00,
     *  0.9991299919297958E+00/
c 
        data xs20/
     1  -.9986540744095632E+00,-.9889061910394767E+00,
     2  -.9620765468436818E+00,-.9120112787934916E+00,
     3  -.8364420012459097E+00,-.7368346580568033E+00,
     4  -.6177324351676943E+00,-.4858378442462243E+00,
     5  -.3488231293339688E+00,-.2138860479584513E+00,
     6  -.8634248677143075E-01,0.3111193741219043E-01,
     7  0.1380473152739090E+00,0.2353032063944305E+00,
     8  0.3242079299763090E+00,0.4063348524379073E+00,
     9  0.4835512326919293E+00,0.5579461539920267E+00,
     *  0.6312922817366760E+00,0.7041273576982980E+00,
     *  0.7751216933289464E+00,0.8412709298383745E+00,
     *  0.8987736513350573E+00,0.9441567063605778E+00,
     *  0.9753590485445780E+00,0.9926402412373609E+00,
     *  0.9990923671468079E+00/
c 
        data xs21/
     1  -.9992363060487008E+00,-.9935060264362001E+00,
     2  -.9770653832623323E+00,-.9451976406284681E+00,
     3  -.8954769136260526E+00,-.8278977566937035E+00,
     4  -.7444283225284132E+00,-.6483123600257108E+00,
     5  -.5433134619204210E+00,-.4330467185925888E+00,
     6  -.3204919981332580E+00,-.2077137860852003E+00,
     7  -.9576670552774280E-01,0.1522277154094187E-01,
     8  0.1257636604763872E+00,0.2364752245389554E+00,
     9  0.3474446913746911E+00,0.4577285771656441E+00,
     *  0.5652103892966860E+00,0.6668265006708167E+00,
     *  0.7590329989935687E+00,0.8383800653777862E+00,
     *  0.9021100371002746E+00,0.9487326090509984E+00,
     *  0.9785325267843799E+00,0.9939094810734341E+00,
     *  0.9992819305311855E+00/
c 
c 
c 
c 
c 
        data ws1/
     1  0.1771180746264163E-01,0.5612651452264730E-01,
     2  0.1028872659089096E+00,0.1492734700434235E+00,
     3  0.1892999449341969E+00,0.2187106006840512E+00,
     4  0.2347117545548852E+00,0.2359136394781541E+00,
     5  0.2223360113569183E+00,0.1953840391207088E+00,
     6  0.1578237538963901E+00,0.1137797139445527E+00,
     7  0.6883480522334508E-01,0.3028595192195138E-01,
     8  0.6593993679462321E-02,0.3031996669242932E-03,
     9  0.1633965639685043E-04,0.3231601005165888E-05,
     *  0.1145792758093868E-05,0.5857421921088699E-06,
     *  0.3905469229907773E-06,0.3037483876299762E-06,
     *  0.2624587299024022E-06,0.2565514895945987E-06,
     *  0.2713171176653337E-06,0.2537653424637160E-06,
     *  0.1416292008526390E-06/
c 
        data ws2/
     1  0.1751013972637862E-01,0.5570230828922125E-01,
     2  0.1023115886072573E+00,0.1486253170285182E+00,
     3  0.1886564849314522E+00,0.2181470941921142E+00,
     4  0.2342945971862009E+00,0.2356938337080623E+00,
     5  0.2223439119377328E+00,0.1956273817124462E+00,
     6  0.1582887295829968E+00,0.1144323476247188E+00,
     7  0.6961569112181959E-01,0.3108249796099407E-01,
     8  0.7175935967504636E-02,0.4477990582212761E-03,
     9  0.3053893766042349E-04,0.6287041731518681E-05,
     *  0.2256748689696217E-05,0.1161092956926643E-05,
     *  0.7768327274224728E-06,0.6052028584908543E-06,
     *  0.5235923338946442E-06,0.5123101728107093E-06,
     *  0.5417572990292029E-06,0.5062476034035852E-06,
     *  0.2823026611111327E-06/
c 
        data ws3/
     1  0.1727428201591147E-01,0.5519570119703394E-01,
     2  0.1016145972257547E+00,0.1478324356324334E+00,
     3  0.1878626974867507E+00,0.2174457077242826E+00,
     4  0.2337680342623946E+00,0.2354056855582140E+00,
     5  0.2223328187002911E+00,0.1959052185287648E+00,
     6  0.1588420874252487E+00,0.1152237035750732E+00,
     7  0.7057550544265401E-01,0.3207292179271303E-01,
     8  0.7910615376579389E-02,0.6551100287833019E-03,
     9  0.5621541341685911E-04,0.1213715557581217E-04,
     *  0.4424635655929440E-05,0.2295341749677093E-05,
     *  0.1542564065857984E-05,0.1204404495344713E-05,
     *  0.1043715014864274E-05,0.1022538598634989E-05,
     *  0.1081202824658253E-05,0.1009120831595893E-05,
     *  0.5620972309135221E-06/
c 
        data ws4/
     1  0.1699817605272736E-01,0.5458899529966729E-01,
     2  0.1007682680814107E+00,0.1468602805098873E+00,
     3  0.1868821349577904E+00,0.2165724863763555E+00,
     4  0.2331042563303415E+00,0.2350302990356705E+00,
     5  0.2222947019903805E+00,0.1962217355505399E+00,
     6  0.1595009676980836E+00,0.1161835633706959E+00,
     7  0.7175420854628201E-01,0.3330102548207065E-01,
     8  0.8835606927527804E-02,0.9500294379634443E-03,
     9  0.1017257332404970E-03,0.2321484193338365E-04,
     *  0.8627464373743499E-05,0.4522832761031223E-05,
     *  0.3056818445543117E-05,0.2393440017178414E-05,
     *  0.2078541005166629E-05,0.2039726459874243E-05,
     *  0.2156432036722319E-05,0.2009544154644968E-05,
     *  0.1117746006196378E-05/
c 
        data ws5/
     1  0.1667360641486185E-01,0.5385848220249739E-01,
     2  0.9973552889098655E-01,0.1456637203110238E+00,
     3  0.1856676859669193E+00,0.2154838494388243E+00,
     4  0.2322676561635685E+00,0.2345433923953367E+00,
     5  0.2222190866390348E+00,0.1965821055218620E+00,
     6  0.1602872111047320E+00,0.1173498376390852E+00,
     7  0.7320209844674562E-01,0.3482110248565914E-01,
     8  0.9998232202192226E-02,0.1366424419321974E-02,
     9  0.1806690528113772E-03,0.4392146172304472E-04,
     *  0.1671270698694803E-04,0.8877484037519057E-05,
     *  0.6042748813156342E-05,0.4748210213036726E-05,
     *  0.4134707950706187E-05,0.4065926686196402E-05,
     *  0.4297642249773373E-05,0.3996966363492554E-05,
     *  0.2219131017874026E-05/
c 
        data ws6/
     1  0.1628967778952740E-01,0.5297308941067170E-01,
     2  0.9846845514062319E-01,0.1441851655341786E+00,
     3  0.1841597873076417E+00,0.2141248394772792E+00,
     4  0.2312130663542856E+00,0.2339136067554571E+00,
     5  0.2220916993253428E+00,0.1969918012808026E+00,
     6  0.1612275204088079E+00,0.1187694237055402E+00,
     7  0.7498092204564473E-01,0.3669906236477438E-01,
     8  0.1145737055501413E-01,0.1949730301020510E-02,
     9  0.3144999568857782E-03,0.8204334273291398E-04,
     *  0.3212335260839250E-04,0.1734409586516374E-04,
     *  0.1191006175252047E-04,0.9400122410437431E-05,
     *  0.8213661350836036E-05,0.8097950410502388E-05,
     *  0.8556721571411240E-05,0.7937995879652195E-05,
     *  0.4397016268017631E-05/
c 
        data ws7/
     1  0.1583244590930609E-01,0.5189394788276375E-01,
     2  0.9690815394917880E-01,0.1423547799571895E+00,
     3  0.1822865516924061E+00,0.2124286402488524E+00,
     4  0.2298843930713464E+00,0.2331004034543616E+00,
     5  0.2218921328195989E+00,0.1974545608183715E+00,
     6  0.1623519364498449E+00,0.1204969302373588E+00,
     7  0.7716245073775910E-01,0.3901124928341566E-01,
     8  0.1328409818170567E-01,0.2759555395647409E-02,
     9  0.5359176204672411E-03,0.1509842411329033E-03,
     *  0.6116790314732482E-04,0.3369484558210246E-04,
     *  0.2338875231151630E-04,0.1856147709376579E-04,
     *  0.1628890522403894E-04,0.1611106799648538E-04,
     *  0.1701548609088130E-04,0.1573402258472944E-04,
     *  0.8689665621290915E-05/
c 
        data ws8/
     1  0.1529290475846048E-01,0.5059338807232442E-01,
     2  0.9501173252704520E-01,0.1401206056032045E+00,
     3  0.1799918674719573E+00,0.2103382580804357E+00,
     4  0.2282270632071664E+00,0.2320561151421790E+00,
     5  0.2215855959526729E+00,0.1979545451999335E+00,
     6  0.1636674604876485E+00,0.1225617002432157E+00,
     7  0.7979327119936655E-01,0.4181537895327204E-01,
     8  0.1554759047591009E-01,0.3868446303802960E-02,
     9  0.8923280754513432E-03,0.2729667485982417E-03,
     *  0.1151364899441443E-03,0.6499854110743930E-04,
     *  0.4571571502692199E-04,0.3652870939288916E-04,
     *  0.3223268978186874E-04,0.3200839992646297E-04,
     *  0.3377882203043600E-04,0.3110208486271519E-04,
     *  0.1711097014574022E-04/
c 
        data ws9/
     1  0.1470412851215310E-01,0.4914779679398834E-01,
     2  0.9288576744019576E-01,0.1375993466940508E+00,
     3  0.1773805365036133E+00,0.2079271071324602E+00,
     4  0.2262701624906356E+00,0.2307603205383358E+00,
     5  0.2211040151731564E+00,0.1983824548546291E+00,
     6  0.1650329250964015E+00,0.1248041686695069E+00,
     7  0.8271867900225360E-01,0.4500199204078366E-01,
     8  0.1823055993562084E-01,0.5336306175983505E-02,
     9  0.1444497523590992E-02,0.4824108731636764E-03,
     *  0.2134675537882129E-03,0.1242084218419045E-03,
     *  0.8878718391827364E-04,0.7155720879247928E-04,
     *  0.6359271655663638E-04,0.6346886463525481E-04,
     *  0.6688972385766177E-04,0.6122802986247385E-04,
     *  0.3350650045804029E-04/
c 
        data ws10/
     1  0.1416560346226777E-01,0.4781114331637116E-01,
     2  0.9090501901311611E-01,0.1352278868340615E+00,
     3  0.1748893393822680E+00,0.2055758362331482E+00,
     4  0.2242922872054199E+00,0.2293544412607500E+00,
     5  0.2204268301160226E+00,0.1985417062126308E+00,
     6  0.1660765206192163E+00,0.1267028553320216E+00,
     7  0.8533567872372616E-01,0.4801978271666135E-01,
     8  0.2100848239299097E-01,0.7104691793172063E-02,
     9  0.2239689217513054E-02,0.8240848303356855E-03,
     *  0.3870228207649433E-03,0.2340714510632819E-03,
     *  0.1708026608624092E-03,0.1392210935753142E-03,
     *  0.1249294421504537E-03,0.1255000145187719E-03,
     *  0.1319133945405553E-03,0.1196516968625745E-03,
     *  0.6494052785705139E-04/
c 
        data ws11/
     1  0.1362300275213502E-01,0.4647753948863603E-01,
     2  0.8896694522517087E-01,0.1329623402009892E+00,
     3  0.1725734804401811E+00,0.2034589097071296E+00,
     4  0.2225846825699069E+00,0.2282215506132418E+00,
     5  0.2199818967545744E+00,0.1988364764636636E+00,
     6  0.1670887667185361E+00,0.1283264678392249E+00,
     7  0.8738694632949308E-01,0.5027039701435784E-01,
     8  0.2318164238130350E-01,0.8768561260364844E-02,
     9  0.3188106614283242E-02,0.1314838522476689E-02,
     *  0.6706932794833888E-03,0.4288130301220894E-03,
     *  0.3222895656696625E-03,0.2672607432687730E-03,
     *  0.2436019845672480E-03,0.2470039509666009E-03,
     *  0.2574928329348151E-03,0.2286187354973051E-03,
     *  0.1216450981857394E-03/
c 
        data ws12/
     1  0.1274651432989933E-01,0.4429216563096405E-01,
     2  0.8581122767093378E-01,0.1293283357296758E+00,
     3  0.1689369279335857E+00,0.2002361432495789E+00,
     4  0.2201207084060132E+00,0.2267786148608004E+00,
     5  0.2197225668817001E+00,0.1998013808133825E+00,
     6  0.1691729559982800E+00,0.1312638157886672E+00,
     7  0.9075427178264187E-01,0.5352972840881522E-01,
     8  0.2581408431181369E-01,0.1053108125533126E-01,
     9  0.4232456952074845E-02,0.1941840371732077E-02,
     *  0.1096397003376001E-02,0.7552557477012680E-03,
     *  0.5920002224485220E-03,0.5064928629984955E-03,
     *  0.4783067990912441E-03,0.4943091711040441E-03,
     *  0.4963755517620936E-03,0.3925951015097143E-03,
     *  0.1734284560406113E-03/
c 
        data ws13/
     1  0.1141799575381782E-01,0.4088034362203177E-01,
     2  0.8078541374604266E-01,0.1234174659799555E+00,
     3  0.1628656566473524E+00,0.1946699969256986E+00,
     4  0.2156426059751066E+00,0.2238541641361183E+00,
     5  0.2186615555617774E+00,0.2007241066873647E+00,
     6  0.1719925955371848E+00,0.1356856847381585E+00,
     7  0.9627336862179725E-01,0.5938460705896637E-01,
     8  0.3094662192475180E-01,0.1406289397153139E-01,
     9  0.6289812998526040E-02,0.3161479937848882E-02,
     *  0.1921544500384297E-02,0.1386503433819097E-02,
     *  0.1115061197069200E-02,0.9744942301673663E-03,
     *  0.9404796757965440E-03,0.9798323250564535E-03,
     *  0.9671088352710734E-03,0.7267608280126125E-03,
     *  0.2718305725660992E-03/
c 
        data ws14/
     1  0.1038064927342182E-01,0.3804092365220404E-01,
     2  0.7636789213614047E-01,0.1179480413840217E+00,
     3  0.1569416637739032E+00,0.1888932774662162E+00,
     4  0.2105760744160665E+00,0.2199781227478337E+00,
     5  0.2163385034638949E+00,0.2001772500206598E+00,
     6  0.1732986624324631E+00,0.1387700820005962E+00,
     7  0.1008777543385223E+00,0.6495172472789384E-01,
     8  0.3645774779169984E-01,0.1839978960335826E-01,
     9  0.9201492762153454E-02,0.5114758473934472E-02,
     *  0.3357307422966970E-02,0.2534017593492176E-02,
     *  0.2091258912485322E-02,0.1867832342894825E-02,
     *  0.1834924836599211E-02,0.1911256709457297E-02,
     *  0.1848192862152322E-02,0.1351036052120744E-02,
     *  0.4897533524853332E-03/
c 
        data ws15/
     1  0.9595608978091610E-02,0.3575956486661454E-01,
     2  0.7263929781826729E-01,0.1131229137537504E+00,
     3  0.1514706728397128E+00,0.1832582789458035E+00,
     4  0.2052489240626110E+00,0.2153874835800367E+00,
     5  0.2128459862601318E+00,0.1980640781062497E+00,
     6  0.1727616779102738E+00,0.1399194336141686E+00,
     7  0.1037089780673814E+00,0.6921728759544303E-01,
     8  0.4142498430922161E-01,0.2301122176986170E-01,
     9  0.1290650433069892E-01,0.8028441157304323E-02,
     *  0.5729794185240898E-02,0.4538650522694080E-02,
     *  0.3864422242059812E-02,0.3552436085455829E-02,
     *  0.3559564951399731E-02,0.3677302957025791E-02,
     *  0.3436397000750934E-02,0.2415690636748959E-02,
     *  0.8543974221896624E-03/
c 
        data ws16/
     1  0.8960921210334843E-02,0.3383739562630139E-01,
     2  0.6938747775214098E-01,0.1087788372153409E+00,
     3  0.1463755231591148E+00,0.1777928831713405E+00,
     4  0.1997975339821350E+00,0.2103101291789431E+00,
     5  0.2084570701561894E+00,0.1946185593122652E+00,
     6  0.1704534139046636E+00,0.1388998849260039E+00,
     7  0.1041234613863590E+00,0.7117992163769307E-01,
     8  0.4478867233408942E-01,0.2719500913855069E-01,
     9  0.1726422975395961E-01,0.1216515418871108E-01,
     *  0.9476515191080042E-02,0.7920221982853269E-02,
     *  0.7054029856387065E-02,0.6789654084752078E-02,
     *  0.6939166956750525E-02,0.6887063749203472E-02,
     *  0.5833928468428052E-02,0.3589899864310428E-02,
     *  0.1123438221795848E-02/
c 
        data ws17/
     1  0.7929784268438364E-02,0.3057067037165196E-01,
     2  0.6368106646788778E-01,0.1009864843062061E+00,
     3  0.1370932312709824E+00,0.1677283162142682E+00,
     4  0.1896981108689510E+00,0.2009084054058046E+00,
     5  0.2004214940620323E+00,0.1885080235763060E+00,
     6  0.1666797206158594E+00,0.1377067209159094E+00,
     7  0.1055776198813336E+00,0.7512193377523672E-01,
     8  0.5071708936548888E-01,0.3428906431657564E-01,
     9  0.2464991621107681E-01,0.1920766495009610E-01,
     *  0.1595186647217316E-01,0.1395930029074346E-01,
     *  0.1299270243853581E-01,0.1292207011844955E-01,
     *  0.1312919093769194E-01,0.1243327204850994E-01,
     *  0.9862156037676716E-02,0.5644161052964884E-02,
     *  0.1629962157706705E-02/
c 
        data ws18/
     1  0.6859665877955894E-02,0.2701826934542556E-01,
     2  0.5723398067502169E-01,0.9188527725941879E-01,
     3  0.1258912056740960E+00,0.1551387481103230E+00,
     4  0.1765084531654525E+00,0.1879111078574665E+00,
     5  0.1883526253933123E+00,0.1780088729191478E+00,
     6  0.1583161210196309E+00,0.1320806138554837E+00,
     7  0.1034616562934428E+00,0.7728687918927112E-01,
     8  0.5719096953491466E-01,0.4381879338585507E-01,
     9  0.3541325081567663E-01,0.2998202375037412E-01,
     *  0.2643617248163136E-01,0.2446479365019776E-01,
     *  0.2392306605948040E-01,0.2407013092278978E-01,
     *  0.2338035029732937E-01,0.2040768334358973E-01,
     *  0.1484364410402514E-01,0.7924336225116750E-02,
     *  0.2191308062536796E-02/
c 
        data ws19/
     1  0.5564442076781132E-02,0.2252669420193703E-01,
     2  0.4879413307696426E-01,0.7962705857108962E-01,
     3  0.1103939104931460E+00,0.1372067430647425E+00,
     4  0.1570454322729574E+00,0.1678406046901328E+00,
     5  0.1686003889429334E+00,0.1596051037327543E+00,
     6  0.1426609525841632E+00,0.1212057910936778E+00,
     7  0.9965207181479874E-01,0.8156810308952017E-01,
     8  0.6808754800029210E-01,0.5841831139561033E-01,
     9  0.5145995116292349E-01,0.4665577200279776E-01,
     *  0.4392545416166342E-01,0.4306555594718443E-01,
     *  0.4301446780971456E-01,0.4186605838550606E-01,
     *  0.3786346819607923E-01,0.3038867186919715E-01,
     *  0.2033992600067372E-01,0.1003458654783682E-01,
     *  0.2588798610330578E-02/
c 
        data ws20/
     1  0.4059935207530622E-02,0.1696534814307852E-01,
     2  0.3773670871251752E-01,0.6276728896060867E-01,
     3  0.8811560658127433E-01,0.1103266990954048E+00,
     4  0.1267408982954409E+00,0.1357448932667944E+00,
     5  0.1370597116802531E+00,0.1319147187112736E+00,
     6  0.1227303826775500E+00,0.1121260057281790E+00,
     7  0.1018983669320483E+00,0.9283806986061852E-01,
     8  0.8523176935699935E-01,0.7933350094570478E-01,
     9  0.7545644099482743E-01,0.7364615109527947E-01,
     *  0.7314786487803730E-01,0.7229684663331534E-01,
     *  0.6917522746612383E-01,0.6247026053865744E-01,
     *  0.5193803452866878E-01,0.3846907197094767E-01,
     *  0.2396881672034206E-01,0.1111532925500303E-01,
     *  0.2726051737813389E-02/
c 
        data ws21/
     1  0.2324140514541588E-02,0.1014618830501652E-01,
     2  0.2354058152821890E-01,0.4059999340227689E-01,
     3  0.5883421598197862E-01,0.7598401271606869E-01,
     4  0.9040128122817617E-01,0.1011920027066544E+00,
     5  0.1081951174598616E+00,0.1118371321993111E+00,
     6  0.1129266946963980E+00,0.1124552161828050E+00,
     7  0.1114226258535849E+00,0.1106501765145693E+00,
     8  0.1105462411120340E+00,0.1109003217600016E+00,
     9  0.1108837204692504E+00,0.1093271699512287E+00,
     *  0.1051177327513327E+00,0.9751407435697319E-01,
     *  0.8631388786201605E-01,0.7191096088978267E-01,
     *  0.5529194155988725E-01,0.3799442816201221E-01,
     *  0.2200560203419779E-01,0.9500930944471324E-02,
     *  0.2183608857232803E-02/
c 
c       return to the user the nodes and weights of all 27 quadratures
c 
        do 1400 i=1,21
        do 1200 j=1,27
c 
        xsout(j,i)=xz(j,i)
        wsout(j,i)=wz(j,i)
 1200 continue
 1400 continue
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine o3rnear7(xsout,wsout)
        implicit real *8 (a-h,o-z)
        save
        dimension xsout(31,21),wsout(31,21)
c 
        dimension xs1(31),xs2(31),xs3(31),xs4(31),xs5(31),
     1      xs6(31),xs7(31),xs8(31),xs9(31),xs10(31),xs11(31),
     2      xs12(31),xs13(31),xs14(31),xs15(31),xs16(31),xs17(31),
     3      xs18(31),xs19(31),xs20(31),xs21(31),xz(31,21)
c 
        equivalence (xz(1,1),xs1(1)),(xz(1,2),xs2(1)),
     1      (xz(1,3),xs3(1)),(xz(1,4),xs4(1)),(xz(1,5),xs5(1)),
     2      (xz(1,6),xs6(1)),(xz(1,7),xs7(1)),(xz(1,8),xs8(1)),
     3      (xz(1,9),xs9(1)),(xz(1,10),xs10(1)),(xz(1,11),xs11(1))
  
        equivalence (xz(1,12),xs12(1)),
     1      (xz(1,13),xs13(1)),(xz(1,14),xs14(1)),(xz(1,15),xs15(1)),
     2      (xz(1,16),xs16(1)),(xz(1,17),xs17(1)),(xz(1,18),xs18(1)),
     3       (xz(1,19),xs19(1)),(xz(1,20),xs20(1)),(xz(1,21),xs21(1))
        dimension ws1(31),ws2(31),ws3(31),ws4(31),ws5(31),
     1      ws6(31),ws7(31),ws8(31),ws9(31),ws10(31),ws11(31),
     2      ws12(31),ws13(31),ws14(31),ws15(31),ws16(31),ws17(31),
     3      ws18(31),ws19(31),ws20(31),ws21(31),wz(31,21)
c 
        equivalence (wz(1,1),ws1(1)),(wz(1,2),ws2(1)),
     1      (wz(1,3),ws3(1)),(wz(1,4),ws4(1)),(wz(1,5),ws5(1)),
     2      (wz(1,6),ws6(1)),(wz(1,7),ws7(1)),(wz(1,8),ws8(1)),
     3      (wz(1,9),ws9(1)),(wz(1,10),ws10(1)),(wz(1,11),ws11(1))
  
        equivalence (wz(1,12),ws12(1)),
     1      (wz(1,13),ws13(1)),(wz(1,14),ws14(1)),(wz(1,15),ws15(1)),
     2      (wz(1,16),ws16(1)),(wz(1,17),ws17(1)),(wz(1,18),ws18(1)),
     3       (wz(1,19),ws19(1)),(wz(1,20),ws20(1)),(wz(1,21),ws21(1))
c 
        data xs1/
     1  -.9944655075872321E+00,-.9620795332146547E+00,
     2  -.8884263566721751E+00,-.7689920084820934E+00,
     3  -.6063586289832279E+00,-.4081957352983611E+00,
     4  -.1857802165608239E+00,0.4734647000247949E-01,
     5  0.2767335142061920E+00,0.4883972794785234E+00,
     6  0.6701950852457215E+00,0.8131673719007346E+00,
     7  0.9128448901500306E+00,0.9705569233088087E+00,
     8  0.9946095488185065E+00,0.9997292903760449E+00,
     9  0.9999890167242436E+00,0.9999961042176668E+00,
     *  0.9999973140608358E+00,0.9999977825660253E+00,
     *  0.9999980550798229E+00,0.9999982578684359E+00,
     *  0.9999984288596897E+00,0.9999985803693388E+00,
     *  0.9999987173163174E+00,0.9999988432177971E+00,
     *  0.9999989633520243E+00,0.9999990888680785E+00,
     *  0.9999992429745499E+00,0.9999994691136977E+00,
     *  0.9999998143123665E+00/
c 
        data xs2/
     1  -.9944841890971526E+00,-.9621874961153871E+00,
     2  -.8887018601205613E+00,-.7694978715389201E+00,
     3  -.6071335679975014E+00,-.4092505734654041E+00,
     4  -.1870972851328216E+00,0.4581178721926271E-01,
     5  0.2750497350361298E+00,0.4866529723324513E+00,
     6  0.6684945745072830E+00,0.8116257415561834E+00,
     7  0.9115806460715652E+00,0.9696768656195884E+00,
     8  0.9941672479521795E+00,0.9996268068861132E+00,
     9  0.9999781774628427E+00,0.9999921686519018E+00,
     *  0.9999946146094321E+00,0.9999955585256532E+00,
     *  0.9999961057479876E+00,0.9999965121262242E+00,
     *  0.9999968544504703E+00,0.9999971576425222E+00,
     *  0.9999974316285565E+00,0.9999976834463296E+00,
     *  0.9999979235674537E+00,0.9999981740570774E+00,
     *  0.9999984809467486E+00,0.9999989309285779E+00,
     *  0.9999996220682605E+00/
  
c 
        data xs3/
     1  -.9945196866471339E+00,-.9623854572052631E+00,
     2  -.8891918399417874E+00,-.7703740076901996E+00,
     3  -.6084438018780528E+00,-.4109938562208172E+00,
     4  -.1892260142452760E+00,0.4338602140491821E-01,
     5  0.2724484052753173E+00,0.4840219598205331E+00,
     6  0.6659949409642290E+00,0.8094230260621806E+00,
     7  0.9098299436345472E+00,0.9684969933860410E+00,
     8  0.9935862511949611E+00,0.9994820255405512E+00,
     9  0.9999581259893078E+00,0.9999844313813560E+00,
     *  0.9999892522798983E+00,0.9999911275596787E+00,
     *  0.9999922179945366E+00,0.9999930289110119E+00,
     *  0.9999937123715041E+00,0.9999943178183117E+00,
     *  0.9999948650086359E+00,0.9999953680006444E+00,
     *  0.9999958477297508E+00,0.9999963482829639E+00,
     *  0.9999969616088524E+00,0.9999978609739255E+00,
     *  0.9999992430626075E+00/
c 
        data xs4/
     1  -.9945599905058116E+00,-.9626085423456778E+00,
     2  -.8897405138728066E+00,-.7713502037118924E+00,
     3  -.6098980554462048E+00,-.4129231655094896E+00,
     4  -.1915770862407647E+00,0.4071018657406955E-01,
     5  0.2695799636576759E+00,0.4811188214023694E+00,
     6  0.6632308775550617E+00,0.8069764544756772E+00,
     7  0.9078688651498799E+00,0.9671535135849965E+00,
     8  0.9929004412184974E+00,0.9992893654397061E+00,
     9  0.9999228166102230E+00,0.9999694538798944E+00,
     *  0.9999786684387935E+00,0.9999823326822206E+00,
     *  0.9999844861102740E+00,0.9999860963456124E+00,
     *  0.9999874565704969E+00,0.9999886625637713E+00,
     *  0.9999897530159221E+00,0.9999907559663803E+00,
     *  0.9999917136689935E+00,0.9999927151818802E+00,
     *  0.9999939455151903E+00,0.9999957510796305E+00,
     *  0.9999985098123481E+00/
  
c 
        data xs5/
     1  -.9945880022540564E+00,-.9627638090548688E+00,
     2  -.8901242172693021E+00,-.7720382971109934E+00,
     3  -.6109341444951569E+00,-.4143162136399279E+00,
     4  -.1933022258509584E+00,0.3870879939901939E-01,
     5  0.2673849738713272E+00,0.4788349620566625E+00,
     6  0.6609803733984305E+00,0.8048948717208850E+00,
     7  0.9061010562534488E+00,0.9658453304859349E+00,
     8  0.9921588799153720E+00,0.9990408374709223E+00,
     9  0.9998622075917143E+00,0.9999407672116709E+00,
     *  0.9999578836098491E+00,0.9999649313138999E+00,
     *  0.9999691478285988E+00,0.9999723303732913E+00,
     *  0.9999750292289984E+00,0.9999774255418548E+00,
     *  0.9999795939033456E+00,0.9999815902399896E+00,
     *  0.9999835007540661E+00,0.9999855074303623E+00,
     *  0.9999879851364115E+00,0.9999916257035502E+00,
     *  0.9999971180224503E+00/
  
c 
        data xs6/
     1  -.9945777320117367E+00,-.9627075219410635E+00,
     2  -.8899918289690512E+00,-.7718201755839157E+00,
     3  -.6106435987590724E+00,-.4139870930133197E+00,
     4  -.1929837227227353E+00,0.3895858954308329E-01,
     5  0.2675058787618128E+00,0.4787718581964409E+00,
     6  0.6606919255886659E+00,0.8043670357082557E+00,
     7  0.9053708466406312E+00,0.9650397041417318E+00,
     8  0.9915167628985919E+00,0.9987439493111030E+00,
     9  0.9997617026389072E+00,0.9998866752303208E+00,
     *  0.9999173870445450E+00,0.9999306783379600E+00,
     *  0.9999388478725999E+00,0.9999451022784560E+00,
     *  0.9999504375055468E+00,0.9999551852374673E+00,
     *  0.9999594864477983E+00,0.9999634530325703E+00,
     *  0.9999672640737045E+00,0.9999712985139044E+00,
     *  0.9999763244748130E+00,0.9999837183261298E+00,
     *  0.9999945934256018E+00/
  
c 
        data xs7/
     1  -.9945517639878102E+00,-.9625621159538174E+00,
     2  -.8896346898868074E+00,-.7711917181624676E+00,
     3  -.6097250184222318E+00,-.4128009074435299E+00,
     4  -.1915898399541843E+00,0.4046930019512796E-01,
     5  0.2690193891816222E+00,0.4801594866292520E+00,
     6  0.6618205081291240E+00,0.8051141961652186E+00,
     7  0.9056538869129267E+00,0.9648756425247582E+00,
     8  0.9911145160106922E+00,0.9984344470896427E+00,
     9  0.9996052131845350E+00,0.9997874680837480E+00,
     *  0.9998395144151742E+00,0.9998638148482653E+00,
     *  0.9998794628300482E+00,0.9998917513142105E+00,
     *  0.9999023330785963E+00,0.9999117640403805E+00,
     *  0.9999202974611834E+00,0.9999281612117727E+00,
     *  0.9999357498001247E+00,0.9999439053051766E+00,
     *  0.9999542782109800E+00,0.9999696041390320E+00,
     *  0.9999907810857571E+00/
c 
        data xs8/
     1  -.9947834443789991E+00,-.9638287599587008E+00,
     2  -.8926933132370045E+00,-.7765153552881461E+00,
     3  -.6174618717574581E+00,-.4227872490815772E+00,
     4  -.2033934395300107E+00,0.2749045424576879E-01,
     5  0.2556514116080726E+00,0.4672639444188678E+00,
     6  0.6502595082902327E+00,0.7956573048190476E+00,
     7  0.8988470608364960E+00,0.9608470995791951E+00,
     8  0.9893404712694527E+00,0.9978204275824389E+00,
     9  0.9993279459008986E+00,0.9996013444548066E+00,
     *  0.9996890773481004E+00,0.9997331888220020E+00,
     *  0.9997630347193798E+00,0.9997870920856022E+00,
     *  0.9998080190070016E+00,0.9998267325451724E+00,
     *  0.9998436955701380E+00,0.9998593926497072E+00,
     *  0.9998747343707162E+00,0.9998916804235829E+00,
     *  0.9999140004235063E+00,0.9999475990991029E+00,
     *  0.9999885403469571E+00/
  
c 
        data xs9/
     1  -.9949858232208237E+00,-.9649563799279117E+00,
     2  -.8954491113387131E+00,-.7813561537765803E+00,
     3  -.6245587820373248E+00,-.4320359831615742E+00,
     4  -.2144475366532240E+00,0.1517611666483432E-01,
     5  0.2427704826836006E+00,0.4546022996215044E+00,
     6  0.6386283468528537E+00,0.7858079169790262E+00,
     7  0.8913566321911099E+00,0.9559608089589819E+00,
     8  0.9867724976497247E+00,0.9967091320546630E+00,
     9  0.9987850290857325E+00,0.9992350197395032E+00,
     *  0.9993926017641184E+00,0.9994751480207702E+00,
     *  0.9995322673505780E+00,0.9995787610770830E+00,
     *  0.9996193646759484E+00,0.9996557690029875E+00,
     *  0.9996888699408526E+00,0.9997196450488888E+00,
     *  0.9997499437451656E+00,0.9997837061886161E+00,
     *  0.9998284299854805E+00,0.9998954339835335E+00,
     *  0.9999760665859727E+00/
  
c 
        data xs10/
     1  -.9951937779900614E+00,-.9661261291997255E+00,
     2  -.8983227357288680E+00,-.7864240121187505E+00,
     3  -.6320201777786687E+00,-.4418100262106586E+00,
     4  -.2262030088734523E+00,0.1982192193836330E-02,
     5  0.2288459231352840E+00,0.4407628766654814E+00,
     6  0.6257269146977760E+00,0.7746455679277851E+00,
     7  0.8825714166207144E+00,0.9498832138356683E+00,
     8  0.9832331918057353E+00,0.9949226518866134E+00,
     9  0.9977907209089778E+00,0.9985280500805076E+00,
     *  0.9988106235347387E+00,0.9989656196981460E+00,
     *  0.9990756565536585E+00,0.9991662538866475E+00,
     *  0.9992457196341971E+00,0.9993171185047551E+00,
     *  0.9993821727970085E+00,0.9994428644035391E+00,
     *  0.9995029872193408E+00,0.9995705238939239E+00,
     *  0.9996603239560111E+00,0.9997935211728011E+00,
     *  0.9999508140854210E+00/
  
c 
        data xs11/
     1  -.9953581864711351E+00,-.9670675083626927E+00,
     2  -.9006700461830118E+00,-.7906228767412440E+00,
     3  -.6382939887517513E+00,-.4501607211907496E+00,
     4  -.2364255607823527E+00,-.9722439498209834E-02,
     5  0.2162025301733730E+00,0.4278361447853942E+00,
     6  0.6132320810994245E+00,0.7632963788403733E+00,
     7  0.8730125917864742E+00,0.9426049855319811E+00,
     8  0.9783970027079540E+00,0.9920605791697004E+00,
     9  0.9959868905411959E+00,0.9971744198204784E+00,
     *  0.9976748124325086E+00,0.9979646722546858E+00,
     *  0.9981769948209691E+00,0.9983542185521330E+00,
     *  0.9985103858858484E+00,0.9986509037639342E+00,
     *  0.9987790880085807E+00,0.9988990595101541E+00,
     *  0.9990188356220777E+00,0.9991549142275591E+00,
     *  0.9993368674336878E+00,0.9996033235474742E+00,
     *  0.9999061377462735E+00/
  
c 
        data xs12/
     1  -.9953195248249636E+00,-.9668279622057760E+00,
     2  -.9000275257814540E+00,-.7893969906724936E+00,
     3  -.6363578848410634E+00,-.4474630510055759E+00,
     4  -.2330064046816424E+00,-.5722608243811877E-02,
     5  0.2205394735201280E+00,0.4321688059084148E+00,
     6  0.6171330110019564E+00,0.7662730651803253E+00,
     7  0.8745525604017443E+00,0.9423034256229266E+00,
     8  0.9763147353535830E+00,0.9891662135137977E+00,
     9  0.9932861973848259E+00,0.9948251667237822E+00,
     *  0.9955937252734459E+00,0.9960932195082983E+00,
     *  0.9964840403774147E+00,0.9968198361322238E+00,
     *  0.9971191488669131E+00,0.9973900761096073E+00,
     *  0.9976389812572504E+00,0.9978756152399963E+00,
     *  0.9981199637752005E+00,0.9984110264875385E+00,
     *  0.9988105543210768E+00,0.9993636682189759E+00,
     *  0.9998727194810002E+00/
c 
        data xs13/
     1  -.9958946969346982E+00,-.9701054191633891E+00,
     2  -.9080803067795163E+00,-.8035285515969093E+00,
     3  -.6570240203879308E+00,-.4743332026610768E+00,
     4  -.2650589557915891E+00,-.4136572960649232E-01,
     5  0.1833097209567385E+00,0.3955683741849827E+00,
     6  0.5833506314584378E+00,0.7372017071327972E+00,
     7  0.8514938759077646E+00,0.9256267443612525E+00,
     8  0.9651358206447203E+00,0.9815962274200518E+00,
     9  0.9875997386748665E+00,0.9901035804900138E+00,
     *  0.9914516491605361E+00,0.9923721967754005E+00,
     *  0.9931113272797540E+00,0.9937537446754228E+00,
     *  0.9943300937276964E+00,0.9948549613253672E+00,
     *  0.9953411934725684E+00,0.9958096525131808E+00,
     *  0.9963027395731636E+00,0.9969003837221804E+00,
     *  0.9977231060173943E+00,0.9988312262132002E+00,
     *  0.9997813282007444E+00/
c 
        data xs14/
     1  -.9959203623822912E+00,-.9703020797425616E+00,
     2  -.9087072769980732E+00,-.8049036914761275E+00,
     3  -.6594727526709405E+00,-.4781492634315502E+00,
     4  -.2704749953000685E+00,-.4853282971947970E-01,
     5  0.1743336023822226E+00,0.3848212509914573E+00,
     6  0.5709626953473206E+00,0.7233939032248400E+00,
     7  0.8366124076205868E+00,0.9103061888054063E+00,
     8  0.9505931383508084E+00,0.9690958743666387E+00,
     9  0.9771423550614959E+00,0.9810756930525635E+00,
     *  0.9834377283861000E+00,0.9851697441702756E+00,
     *  0.9866110428720056E+00,0.9878787715058857E+00,
     *  0.9890174861848304E+00,0.9900524149935639E+00,
     *  0.9910113380529195E+00,0.9919443155276597E+00,
     *  0.9929500759250327E+00,0.9941914243095329E+00,
     *  0.9958579344176364E+00,0.9979346472329702E+00,
     *  0.9996100881113690E+00/
c 
        data xs15/
     1  -.9961698243216317E+00,-.9717917234292071E+00,
     2  -.9125155870239987E+00,-.8118433061112071E+00,
     3  -.6700205740858590E+00,-.4924426540445149E+00,
     4  -.2883279583866581E+00,-.6947305844756568E-01,
     5  0.1510028012890958E+00,0.3599357191734438E+00,
     6  0.5453952152685171E+00,0.6979449759050886E+00,
     7  0.8119137947939530E+00,0.8868042600968422E+00,
     8  0.9287126200816678E+00,0.9492428394115695E+00,
     9  0.9592970668658096E+00,0.9649340069914740E+00,
     *  0.9687452124905414E+00,0.9717582258122053E+00,
     *  0.9743567783688103E+00,0.9766823316235838E+00,
     *  0.9787971683477614E+00,0.9807467341861720E+00,
     *  0.9825940721328390E+00,0.9844608040088927E+00,
     *  0.9865764083263389E+00,0.9892830651958080E+00,
     *  0.9928598664795454E+00,0.9968339520687579E+00,
     *  0.9994718632346848E+00/
  
c 
        data xs16/
     1  -.9964668572577557E+00,-.9736171411920959E+00,
     2  -.9173262594011719E+00,-.8208760340687278E+00,
     3  -.6841569556708590E+00,-.5121591915789122E+00,
     4  -.3136822257877192E+00,-.1001262224847987E+00,
     5  0.1157263052120174E+00,0.3209569037858593E+00,
     6  0.5037707338330486E+00,0.6547648280078453E+00,
     7  0.7682733347061466E+00,0.8439683772580244E+00,
     8  0.8882321933887385E+00,0.9121599926182984E+00,
     9  0.9255955388074875E+00,0.9342376298886286E+00,
     *  0.9407472215067074E+00,0.9462143770380120E+00,
     *  0.9510491897414171E+00,0.9554150829470276E+00,
     *  0.9594023539787027E+00,0.9631022063381163E+00,
     *  0.9666676407251242E+00,0.9703968634436939E+00,
     *  0.9748015460635597E+00,0.9805184372638562E+00,
     *  0.9877606833572548E+00,0.9949991915144323E+00,
     *  0.9992271430213111E+00/
c 
        data xs17/
     1  -.9967405376337557E+00,-.9753926181981355E+00,
     2  -.9222194126682073E+00,-.8304035034531652E+00,
     3  -.6995418416635963E+00,-.5342520655587543E+00,
     4  -.3429378325120463E+00,-.1366277500731720E+00,
     5  0.7219807072751494E-01,0.2708185955148712E+00,
     6  0.4475273702869449E+00,0.5929989217243805E+00,
     7  0.7019729466963911E+00,0.7752695967650627E+00,
     8  0.8203426171164205E+00,0.8475914723615943E+00,
     9  0.8654001879919984E+00,0.8786505686683186E+00,
     *  0.8896446560539997E+00,0.8993245636488112E+00,
     *  0.9080697304427647E+00,0.9160691802342575E+00,
     *  0.9234795705206629E+00,0.9305273598899289E+00,
     *  0.9376356536814832E+00,0.9455556529500460E+00,
     *  0.9553048190922588E+00,0.9675824448858004E+00,
     *  0.9813720150127796E+00,0.9930292931861895E+00,
     *  0.9989831448672581E+00/
  
c 
        data xs18/
     1  -.9971768143599080E+00,-.9783121557605008E+00,
     2  -.9304533987232991E+00,-.8466929098946688E+00,
     3  -.7261420641093478E+00,-.5727712654429050E+00,
     4  -.3942934069603821E+00,-.2010915323635550E+00,
     5  -.5120431957343907E-02,0.1812564430908421E+00,
     6  0.3465669449582218E+00,0.4818886325812204E+00,
     7  0.5831256111710610E+00,0.6529786059407592E+00,
     8  0.6996838389256239E+00,0.7324025994814276E+00,
     9  0.7576956338451355E+00,0.7790421351325876E+00,
     *  0.7979619602854304E+00,0.8150960590713974E+00,
     *  0.8307734475772028E+00,0.8452643316909285E+00,
     *  0.8589444886792459E+00,0.8724977464369679E+00,
     *  0.8871381940374026E+00,0.9045459471901871E+00,
     *  0.9260433138712574E+00,0.9508058855619583E+00,
     *  0.9744706471999329E+00,0.9912994561433637E+00,
     *  0.9988141537823880E+00/
c 
        data xs19/
     1  -.9976536528087926E+00,-.9816777582614181E+00,
     2  -.9404414784738074E+00,-.8673931199927087E+00,
     3  -.7614488126652300E+00,-.6260946201429638E+00,
     4  -.4684163579206266E+00,-.2981311778937596E+00,
     5  -.1265410154097492E+00,0.3472899468891777E-01,
     6  0.1754802721045550E+00,0.2893942298306225E+00,
     7  0.3762743517191440E+00,0.4416064082959164E+00,
     8  0.4928519859475899E+00,0.5358079950853647E+00,
     9  0.5736889236465526E+00,0.6079961901793089E+00,
     *  0.6394468080874416E+00,0.6684865966229452E+00,
     *  0.6955608448733308E+00,0.7213488413542519E+00,
     *  0.7470755533603697E+00,0.7747844915166053E+00,
     *  0.8070511052104565E+00,0.8456230674294232E+00,
     *  0.8892337760413731E+00,0.9323834894119412E+00,
     *  0.9675729088729634E+00,0.9896281363079520E+00,
     *  0.9986494209416814E+00/
  
c 
        data xs20/
     1  -.9982553254713196E+00,-.9861314978864718E+00,
     2  -.9542662670164700E+00,-.8972301619182923E+00,
     3  -.8143165870634067E+00,-.7089857106488910E+00,
     4  -.5880218559329923E+00,-.4604129714319145E+00,
     5  -.3356776236444065E+00,-.2215171191618026E+00,
     6  -.1215714801447458E+00,-.3506338172369516E-01,
     7  0.4110192037798943E-01,0.1098936387155567E+00,
     8  0.1732182426895841E+00,0.2321173436567132E+00,
     9  0.2871873654681385E+00,0.3388877029268875E+00,
     *  0.3877669852940040E+00,0.4347138872788996E+00,
     *  0.4812990082149182E+00,0.5300677276595822E+00,
     *  0.5842597749821736E+00,0.6464820591838241E+00,
     *  0.7167251763313514E+00,0.7909734706318936E+00,
     *  0.8619101703010505E+00,0.9216483769269148E+00,
     *  0.9645846858112055E+00,0.9891785203144901E+00,
     *  0.9986365817337306E+00/
  
c 
        data xs21/
     1  -.9992631127661414E+00,-.9938744126539175E+00,
     2  -.9789616950083361E+00,-.9512269532328911E+00,
     3  -.9098199664961912E+00,-.8559960124090016E+00,
     4  -.7922873473186982E+00,-.7216475360606280E+00,
     5  -.6468283122732528E+00,-.5700639594562795E+00,
     6  -.4929984884959127E+00,-.4167417982294879E+00,
     7  -.3419588205431744E+00,-.2689313196658408E+00,
     8  -.1975583217709415E+00,-.1272779633158709E+00,
     9  -.5692341548067015E-01,0.1539551892499835E-01,
     *  0.9218058136267615E-01,0.1760141501506842E+00,
     *  0.2686830134635829E+00,0.3702725774297857E+00,
     *  0.4785642546520885E+00,0.5890317183624278E+00,
     *  0.6955572415210118E+00,0.7916252438094795E+00,
     *  0.8715827817468513E+00,0.9317060895985667E+00,
     *  0.9709979052306097E+00,0.9916637784911320E+00,
     *  0.9990062348760766E+00/
c 
c 
c 
c 
c 
        data ws1/
     1  0.1567305174138062E-01,0.5154171523620986E-01,
     2  0.9646355484800159E-01,0.1419425739770370E+00,
     3  0.1820221071074628E+00,0.2124018996829053E+00,
     4  0.2301494866822719E+00,0.2336725842493603E+00,
     5  0.2227480125744991E+00,0.1985333449512740E+00,
     6  0.1635473103107129E+00,0.1216306311337212E+00,
     7  0.7792398803660082E-01,0.3888021304592759E-01,
     8  0.1178977117567845E-01,0.1054352416957817E-02,
     9  0.1993261645135999E-04,0.2173150203852364E-05,
     *  0.6686111220941013E-06,0.3354359675118412E-06,
     *  0.2274271900527882E-06,0.1836673570432841E-06,
     *  0.1601171742973237E-06,0.1436404135399559E-06,
     *  0.1307950997489170E-06,0.1218223328384269E-06,
     *  0.1201971877508939E-06,0.1345751548531195E-06,
     *  0.1807006634712093E-06,0.2829781799760282E-06,
     *  0.3740816050603530E-06/
  
c 
        data ws2/
     1  0.1562330573891008E-01,0.5141189028770778E-01,
     2  0.9626110046863570E-01,0.1416884041284469E+00,
     3  0.1817428427586859E+00,0.2121261637795801E+00,
     4  0.2299052778039175E+00,0.2334855369123402E+00,
     5  0.2226402051179733E+00,0.1985227351006142E+00,
     6  0.1636473447776532E+00,0.1218490343657407E+00,
     7  0.7825877262352991E-01,0.3930603667419522E-01,
     8  0.1221415497831949E-01,0.1268301149019484E-02,
     9  0.3813290413274948E-04,0.4394512930694111E-05,
     *  0.1349347062350317E-05,0.6744713882611632E-06,
     *  0.4560981809924142E-06,0.3678280157198111E-06,
     *  0.3204669299651306E-06,0.2874068885263719E-06,
     *  0.2616489968229482E-06,0.2436025870872701E-06,
     *  0.2401035335380554E-06,0.2682886542642883E-06,
     *  0.3595546118793897E-06,0.5635610730838656E-06,
     *  0.7553270620306970E-06/
  
c 
        data ws3/
     1  0.1553000609120585E-01,0.5118041567900034E-01,
     2  0.9591485977946734E-01,0.1412701464328399E+00,
     3  0.1813011746863440E+00,0.2117095634722154E+00,
     4  0.2295578225103275E+00,0.2332444218110623E+00,
     5  0.2225342583642213E+00,0.1985718685963358E+00,
     6  0.1638617835477637E+00,0.1222268508990968E+00,
     7  0.7877936289682382E-01,0.3991298149523668E-01,
     8  0.1277237819784901E-01,0.1542163157516957E-02,
     9  0.6887780211131587E-04,0.8606164663138482E-05,
     *  0.2675325763611053E-05,0.1342609806592681E-05,
     *  0.9097005993992560E-06,0.7342577634590316E-06,
     *  0.6399000394570517E-06,0.5739601759699536E-06,
     *  0.5225860762019549E-06,0.4866293951976158E-06,
     *  0.4797502611841715E-06,0.5361620310285726E-06,
     *  0.7185931662952173E-06,0.1126474963881344E-05,
     *  0.1511566775766972E-05/
  
c 
        data ws4/
     1  0.1542436116963546E-01,0.5092113602975765E-01,
     2  0.9553023951254758E-01,0.1408084053191708E+00,
     3  0.1808158293548242E+00,0.2112532477929810E+00,
     4  0.2291779304160337E+00,0.2329806378389571E+00,
     5  0.2224171175251717E+00,0.1986222797457441E+00,
     6  0.1640903593844443E+00,0.1226316709242673E+00,
     7  0.7933998519041138E-01,0.4057341490162763E-01,
     8  0.1339266084013486E-01,0.1863909270902019E-02,
     9  0.1161610682172563E-03,0.1620812776851309E-04,
     *  0.5193288295977721E-05,0.2641196265446030E-05,
     *  0.1802924838094605E-05,0.1460237160945363E-05,
     *  0.1274243008590799E-05,0.1143559935252796E-05,
     *  0.1041665926993732E-05,0.9707594276301853E-06,
     *  0.9586396726632388E-06,0.1074198201741971E-05,
     *  0.1442767646403168E-05,0.2259499211555906E-05,
     *  0.2995075328442739E-05/
  
c 
        data ws5/
     1  0.1535093642282632E-01,0.5074019898772410E-01,
     2  0.9525880794458680E-01,0.1404765137270513E+00,
     3  0.1804575675991560E+00,0.2109032727647236E+00,
     4  0.2288689007548153E+00,0.2327416203504053E+00,
     5  0.2222725962631561E+00,0.1985919727398070E+00,
     6  0.1641896017930458E+00,0.1228718547641533E+00,
     7  0.7972787781317124E-01,0.4110013449557780E-01,
     8  0.1397049083361358E-01,0.2215761917426112E-02,
     9  0.1845897390995648E-03,0.2944356953387587E-04,
     *  0.9880615296891757E-05,0.5137621372860285E-05,
     *  0.3551697552054604E-05,0.2893610043665112E-05,
     *  0.2530670003344495E-05,0.2273190371130807E-05,
     *  0.2072158276222024E-05,0.1933834110385130E-05,
     *  0.1915841783060099E-05,0.2158014013498593E-05,
     *  0.2910340971269644E-05,0.4546116494699884E-05,
     *  0.5869925199412618E-05/
  
c 
        data ws6/
     1  0.1537791038544638E-01,0.5080405808373878E-01,
     2  0.9534355599703651E-01,0.1405593399012641E+00,
     3  0.1805161269089272E+00,0.2109192887131784E+00,
     4  0.2288301395176530E+00,0.2326425704252604E+00,
     5  0.2221146743156525E+00,0.1983841104816320E+00,
     6  0.1639511711893792E+00,0.1226397320813879E+00,
     7  0.7957016262395203E-01,0.4112846110971691E-01,
     8  0.1426655373476672E-01,0.2533250782271223E-02,
     9  0.2745636846428049E-03,0.5123260916177435E-04,
     *  0.1833339560226026E-04,0.9854522558749695E-05,
     *  0.6945096314446967E-05,0.5709114611268205E-05,
     *  0.5010075068135659E-05,0.4506679582548975E-05,
     *  0.4113011450306575E-05,0.3847855867227389E-05,
     *  0.3834078878412775E-05,0.4359173594873843E-05,
     *  0.5919572246297566E-05,0.9189415699383861E-05,
     *  0.1127767426586042E-04/
  
c 
        data ws7/
     1  0.1544639929888745E-01,0.5097373534960820E-01,
     2  0.9559154250963657E-01,0.1408470902407098E+00,
     3  0.1808016641859134E+00,0.2111626882113437E+00,
     4  0.2289969216212691E+00,0.2327055351854815E+00,
     5  0.2220548203674007E+00,0.1981913898550024E+00,
     6  0.1636274086620590E+00,0.1222067720497855E+00,
     7  0.7909062246869616E-01,0.4074721317307598E-01,
     8  0.1420280226936921E-01,0.2720695568339052E-02,
     9  0.3706372490621777E-03,0.8310011583635531E-04,
     *  0.3263423384354126E-04,0.1853762987111605E-04,
     *  0.1352467977283689E-04,0.1129301003746884E-04,
     *  0.9952263573263064E-05,0.8948107870556169E-05,
     *  0.8153919284907462E-05,0.7634649200737944E-05,
     *  0.7676414453423037E-05,0.8904036905773481E-05,
     *  0.1230460479105474E-04,0.1883472560724852E-04,
     *  0.2034087825226960E-04/
  
c 
        data ws8/
     1  0.1484102105530438E-01,0.4952169164989217E-01,
     2  0.9350984336192242E-01,0.1384536245747170E+00,
     3  0.1784212142793364E+00,0.2090891813855068E+00,
     4  0.2274710335374408E+00,0.2319049953187465E+00,
     5  0.2220897767784499E+00,0.1991016593875283E+00,
     6  0.1653714048552531E+00,0.1246347937148981E+00,
     7  0.8189439149418057E-01,0.4338814250724484E-01,
     8  0.1595559682677163E-01,0.3347552912253121E-02,
     9  0.5260921014523470E-03,0.1347358285609872E-03,
     *  0.5757710233255557E-04,0.3468570517875379E-04,
     *  0.2623902092788501E-04,0.2226216921609174E-04,
     *  0.1972697632306731E-04,0.1777079091539804E-04,
     *  0.1623002398598285E-04,0.1530522720534978E-04,
     *  0.1569157533689590E-04,0.1881603214658581E-04,
     *  0.2686534218389577E-04,0.4081968568478407E-04,
     *  0.3109044419246348E-04/
  
c 
        data ws9/
     1  0.1430796369575916E-01,0.4821230651903280E-01,
     2  0.9160654094933897E-01,0.1362352422284918E+00,
     3  0.1761736500892067E+00,0.2070736372469766E+00,
     4  0.2259091169920420E+00,0.2309712798256046E+00,
     5  0.2219058614302519E+00,0.1997286494326303E+00,
     6  0.1667958609417143E+00,0.1267455806311228E+00,
     7  0.8444717750013346E-01,0.4595400775128233E-01,
     8  0.1791761771584079E-01,0.4292092834625280E-02,
     9  0.8210589703508053E-03,0.2359990683800498E-03,
     *  0.1061844781944135E-03,0.6585017527253484E-04,
     *  0.5054831345810571E-04,0.4313085757020769E-04,
     *  0.3832734851102279E-04,0.3461859652293028E-04,
     *  0.3173437121889809E-04,0.3010322479229581E-04,
     *  0.3112494357385435E-04,0.3762751287381164E-04,
     *  0.5384763728283849E-04,0.8078892971081582E-04,
     *  0.6235014323418083E-04/
  
c 
        data ws10/
     1  0.1375787534267694E-01,0.4684609792855694E-01,
     2  0.8960972380606975E-01,0.1338928001613681E+00,
     3  0.1737760048362875E+00,0.2048878531755774E+00,
     4  0.2241675083722450E+00,0.2298636072836044E+00,
     5  0.2215704925225051E+00,0.2002407662451153E+00,
     6  0.1681514572224055E+00,0.1288414855220875E+00,
     7  0.8705521734036371E-01,0.4867805713048091E-01,
     8  0.2015690117393390E-01,0.5520451218639473E-02,
     9  0.1274928296967606E-02,0.4116235509249301E-03,
     *  0.1960673415605392E-03,0.1256867080016886E-03,
     *  0.9812554488571406E-04,0.8428646582715163E-04,
     *  0.7510230916186938E-04,0.6796034602229947E-04,
     *  0.6245504235756126E-04,0.5951941806372114E-04,
     *  0.6200845613883258E-04,0.7550783135839255E-04,
     *  0.1079603214812445E-03,0.1588006591675126E-03,
     *  0.1235830193361658E-03/
  
c 
        data ws11/
     1  0.1332000643003517E-01,0.4573158501922663E-01,
     2  0.8794512434036367E-01,0.1318918818697604E+00,
     3  0.1716641340317318E+00,0.2028805795640933E+00,
     4  0.2224620277300043E+00,0.2286306673129014E+00,
     5  0.2209447661028835E+00,0.2003099396254434E+00,
     6  0.1689459185857347E+00,0.1303269757137067E+00,
     7  0.8912592688313775E-01,0.5111760190728069E-01,
     8  0.2248983264527838E-01,0.7041399028069553E-02,
     9  0.1944669053041151E-02,0.7062526922008552E-03,
     *  0.3591911748428529E-03,0.2397838652659448E-03,
     *  0.1910981651335234E-03,0.1654118373679977E-03,
     *  0.1477343003274534E-03,0.1338137406640749E-03,
     *  0.1231938578996028E-03,0.1180045760132805E-03,
     *  0.1242137095790151E-03,0.1528356784343151E-03,
     *  0.2183817581889558E-03,0.3127486219616974E-03,
     *  0.2342276909841958E-03/
c 
        data ws12/
     1  0.1342584431075245E-01,0.4603360372790202E-01,
     2  0.8844501584890105E-01,0.1325496510403276E+00,
     3  0.1724141247252471E+00,0.2036383650004701E+00,
     4  0.2231302735927519E+00,0.2291064956517281E+00,
     5  0.2211267267666516E+00,0.2001049441294476E+00,
     6  0.1682763031092044E+00,0.1291436718784255E+00,
     7  0.8745000584292340E-01,0.4918169654164778E-01,
     8  0.2104247721391757E-01,0.6879615853108142E-02,
     9  0.2311877388729059E-02,0.1017023922799583E-02,
     *  0.5917585892921103E-03,0.4309903787485237E-03,
     *  0.3585703241095415E-03,0.3157632549301945E-03,
     *  0.2840920893937723E-03,0.2587494805618271E-03,
     *  0.2405746276460130E-03,0.2358592405608865E-03,
     *  0.2592441736263394E-03,0.3334116577105946E-03,
     *  0.4768640385832543E-03,0.6011551050786283E-03,
     *  0.3351461784525817E-03/
  
c 
        data ws13/
     1  0.1189323308768094E-01,0.4218941115976453E-01,
     2  0.8287526068309464E-01,0.1261049460705425E+00,
     3  0.1659210913578427E+00,0.1978389651559882E+00,
     4  0.2186606116599262E+00,0.2264644540793782E+00,
     5  0.2206315479358666E+00,0.2018547144336691E+00,
     6  0.1721136999677289E+00,0.1346305285757694E+00,
     7  0.9383397508308533E-01,0.5534496195560899E-01,
     8  0.2568274985379933E-01,0.9499405057061601E-02,
     9  0.3612288008747397E-02,0.1732589427833543E-02,
     *  0.1070097756576928E-02,0.8077149439646942E-03,
     *  0.6831117354101597E-03,0.6062920402994790E-03,
     *  0.5486382133336372E-03,0.5030949846730595E-03,
     *  0.4726215557271604E-03,0.4710593346593559E-03,
     *  0.5283021820418330E-03,0.6878115100440383E-03,
     *  0.9760832433539909E-03,0.1169828910864041E-02,
     *  0.5908764555134290E-03/
  
c 
        data ws14/
     1  0.1181746483370261E-01,0.4190337283570980E-01,
     2  0.8229110880329142E-01,0.1251916996058635E+00,
     3  0.1646925980124072E+00,0.1963431320657458E+00,
     4  0.2169702292024515E+00,0.2246683733445255E+00,
     5  0.2188256977302282E+00,0.2001339249195034E+00,
     6  0.1705675417740692E+00,0.1333558699834242E+00,
     7  0.9300490004035257E-01,0.5540338099219928E-01,
     8  0.2722475373735015E-01,0.1176494487165089E-01,
     9  0.5343291366291650E-02,0.2907923146215867E-02,
     *  0.1958804867219297E-02,0.1556130256269939E-02,
     *  0.1343779612551222E-02,0.1198361644127680E-02,
     *  0.1082968750420829E-02,0.9912167760326887E-03,
     *  0.9344915312428863E-03,0.9475559539003684E-03,
     *  0.1091549981330977E-02,0.1425437898281398E-02,
     *  0.1917794562312057E-02,0.2098694641196311E-02,
     *  0.1042968496355189E-02/
c 
        data ws15/
     1  0.1114079615887291E-01,0.4009306825657101E-01,
     2  0.7951218449465387E-01,0.1217648045748375E+00,
     3  0.1609599050596049E+00,0.1926389469602890E+00,
     4  0.2136034228971680E+00,0.2218984338589658E+00,
     5  0.2168379486173962E+00,0.1990184961697697E+00,
     6  0.1703069073766412E+00,0.1338244557721661E+00,
     7  0.9400217231267607E-01,0.5679619077530952E-01,
     8  0.2908727680620019E-01,0.1386804107633637E-01,
     9  0.7207382293592747E-02,0.4466211427724650E-02,
     *  0.3313659688124056E-02,0.2770165431516991E-02,
     *  0.2448119369226408E-02,0.2212488373663621E-02,
     *  0.2024262786032576E-02,0.1884714465794771E-02,
     *  0.1829277123232813E-02,0.1942110414378954E-02,
     *  0.2348168200346032E-02,0.3122590662510006E-02,
     *  0.3972506289053378E-02,0.3640492374892659E-02,
     *  0.1464763546768443E-02/
c 
        data ws16/
     1  0.1032755778844759E-01,0.3781974466717748E-01,
     2  0.7585597577315439E-01,0.1170372122382980E+00,
     3  0.1555482118302953E+00,0.1869607812611828E+00,
     4  0.2080728392582835E+00,0.2168879149926259E+00,
     5  0.2126430479902612E+00,0.1958313270533018E+00,
     6  0.1682072844086036E+00,0.1328101109088430E+00,
     7  0.9412110332025989E-01,0.5835633325499472E-01,
     8  0.3215089890594262E-01,0.1739026331700676E-01,
     9  0.1040054291525974E-01,0.7302242796914201E-02,
     *  0.5884109385292323E-02,0.5111864564570415E-02,
     *  0.4582494681258513E-02,0.4163293847352906E-02,
     *  0.3825529035126461E-02,0.3598226635139673E-02,
     *  0.3580827233344143E-02,0.3963872280956397E-02,
     *  0.4958875151249577E-02,0.6536552047558460E-02,
     *  0.7704375411698520E-02,0.6175901041974564E-02,
     *  0.2190672765294318E-02/
  
c 
        data ws17/
     1  0.9563510583126924E-02,0.3552009565602860E-01,
     2  0.7194146108649255E-01,0.1117336352682337E+00,
     3  0.1491998741372322E+00,0.1799599126947727E+00,
     4  0.2008109283963105E+00,0.2097090294992136E+00,
     5  0.2058016261466201E+00,0.1894660880191735E+00,
     6  0.1623698956564585E+00,0.1276594863123827E+00,
     7  0.9040428906757651E-01,0.5750197610027940E-01,
     8  0.3447715208440477E-01,0.2143517492169124E-01,
     9  0.1497548834818709E-01,0.1188719676422242E-01,
     *  0.1024491920543260E-01,0.9172538261652746E-02,
     *  0.8347027633687488E-02,0.7676471619449177E-02,
     *  0.7179949481581323E-02,0.6983445796258300E-02,
     *  0.7358582429047456E-02,0.8659637467846537E-02,
     *  0.1097998657221848E-01,0.1343347191118180E-01,
     *  0.1348853439396635E-01,0.9137078791714479E-02,
     *  0.2921523168326093E-02/
  
c 
        data ws18/
     1  0.8331597390959877E-02,0.3165874769890607E-01,
     2  0.6519729945650354E-01,0.1024512900126256E+00,
     3  0.1379702668590188E+00,0.1674681536183224E+00,
     4  0.1877388085701287E+00,0.1966532096345358E+00,
     5  0.1932072416732534E+00,0.1776078900617200E+00,
     6  0.1514870286012620E+00,0.1184313435231424E+00,
     7  0.8451114085768412E-01,0.5666070215141439E-01,
     8  0.3833911565262987E-01,0.2819502346913735E-01,
     9  0.2294767038560412E-01,0.1997995115046362E-01,
     *  0.1795781813827509E-01,0.1636153151910948E-01,
     *  0.1503584629656272E-01,0.1400522405080995E-01,
     *  0.1346460951696656E-01,0.1384395203947997E-01,
     *  0.1573000507107607E-01,0.1933737654787007E-01,
     *  0.2355525784298679E-01,0.2520344625175311E-01,
     *  0.2107720980460687E-01,0.1212385436703839E-01,
     *  0.3467382039328983E-02/
  
c 
        data ws19/
     1  0.6961603888043030E-02,0.2702851270489275E-01,
     2  0.5653143242721934E-01,0.8972940204565981E-01,
     3  0.1215620895744487E+00,0.1479493840036150E+00,
     4  0.1657704451846495E+00,0.1728958885393503E+00,
     5  0.1683187174868164E+00,0.1524791677988072E+00,
     6  0.1279251833753607E+00,0.9987423683363809E-01,
     7  0.7487299385031680E-01,0.5710328613101273E-01,
     8  0.4636504938799323E-01,0.4006330004980829E-01,
     9  0.3593273945689146E-01,0.3279336996722267E-01,
     *  0.3017815582749940E-01,0.2797203787772416E-01,
     *  0.2628434576910827E-01,0.2548986225620114E-01,
     *  0.2630726031118107E-01,0.2956301866481958E-01,
     *  0.3530154433012033E-01,0.4166117607485074E-01,
     *  0.4460164121176238E-01,0.4036488142001461E-01,
     *  0.2910091504893930E-01,0.1502191936269661E-01,
     *  0.3996436993304447E-02/
  
c 
        data ws20/
     1  0.5206369340927443E-02,0.2069056107795850E-01,
     2  0.4394953512657070E-01,0.7022947845230525E-01,
     3  0.9500712540518698E-01,0.1145059949749545E+00,
     4  0.1258964767871100E+00,0.1276941763869927E+00,
     5  0.1204502082643649E+00,0.1072556271627058E+00,
     6  0.9280718879393268E-01,0.8077094011084300E-01,
     7  0.7206800037781318E-01,0.6582733445013042E-01,
     8  0.6098789743417475E-01,0.5690475187473606E-01,
     9  0.5330836535290272E-01,0.5017978331216485E-01,
     *  0.4772099645638130E-01,0.4642830478836298E-01,
     *  0.4716916307787312E-01,0.5093147421790344E-01,
     *  0.5792314119662898E-01,0.6654761733092769E-01,
     *  0.7325785599491332E-01,0.7396743307565085E-01,
     *  0.6654110346195931E-01,0.5197251512817858E-01,
     *  0.3363147555769987E-01,0.1609979808681240E-01,
     *  0.4069306483386881E-02/
  
c 
        data ws21/
     1  0.2229144594125731E-02,0.9411669049797287E-02,
     2  0.2097554838483192E-01,0.3462821547352111E-01,
     3  0.4796168426167274E-01,0.5924929324898813E-01,
     4  0.6766523019134946E-01,0.7315070747541140E-01,
     5  0.7611444457112218E-01,0.7714142857805655E-01,
     6  0.7680615391957664E-01,0.7559751034773801E-01,
     7  0.7392153597336114E-01,0.7214898342769580E-01,
     8  0.7068806503131868E-01,0.7006461810915934E-01,
     9  0.7096077801899195E-01,0.7410308188366332E-01,
     *  0.7991294680711099E-01,0.8806784748641470E-01,
     *  0.9729491411653576E-01,0.1055149429537117E+00,
     *  0.1103013864258157E+00,0.1095928851786440E+00,
     *  0.1023506066724376E+00,0.8882536075759531E-01,
     *  0.7045576174655725E-01,0.4961621168923225E-01,
     *  0.2933506546357618E-01,0.1290104603964836E-01,
     *  0.3012932118861261E-02/
c 
c       return to the user the nodes and weights of all 31 quadratures
c 
        do 1400 i=1,21
        do 1200 j=1,31
c 
        xsout(j,i)=xz(j,i)
        wsout(j,i)=wz(j,i)
 1200 continue
 1400 continue
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine o4rnear7(xsout,wsout)
        implicit real *8 (a-h,o-z)
        save
        dimension xsout(38,21),wsout(38,21)
c 
        dimension xs1(38),xs2(38),xs3(38),xs4(38),xs5(38),
     1      xs6(38),xs7(38),xs8(38),xs9(38),xs10(38),xs11(38),
     2      xs12(38),xs13(38),xs14(38),xs15(38),xs16(38),xs17(38),
     3      xs18(38),xs19(38),xs20(38),xs21(38),xz(38,21)
c 
        equivalence (xz(1,1),xs1(1)),(xz(1,2),xs2(1)),
     1      (xz(1,3),xs3(1)),(xz(1,4),xs4(1)),(xz(1,5),xs5(1)),
     2      (xz(1,6),xs6(1)),(xz(1,7),xs7(1)),(xz(1,8),xs8(1)),
     3      (xz(1,9),xs9(1)),(xz(1,10),xs10(1)),(xz(1,11),xs11(1))
  
        equivalence (xz(1,12),xs12(1)),
     1      (xz(1,13),xs13(1)),(xz(1,14),xs14(1)),(xz(1,15),xs15(1)),
     2      (xz(1,16),xs16(1)),(xz(1,17),xs17(1)),(xz(1,18),xs18(1)),
     3       (xz(1,19),xs19(1)),(xz(1,20),xs20(1)),(xz(1,21),xs21(1))
        dimension ws1(38),ws2(38),ws3(38),ws4(38),ws5(38),
     1      ws6(38),ws7(38),ws8(38),ws9(38),ws10(38),ws11(38),
     2      ws12(38),ws13(38),ws14(38),ws15(38),ws16(38),ws17(38),
     3      ws18(38),ws19(38),ws20(38),ws21(38),wz(38,21)
c 
        equivalence (wz(1,1),ws1(1)),(wz(1,2),ws2(1)),
     1      (wz(1,3),ws3(1)),(wz(1,4),ws4(1)),(wz(1,5),ws5(1)),
     2      (wz(1,6),ws6(1)),(wz(1,7),ws7(1)),(wz(1,8),ws8(1)),
     3      (wz(1,9),ws9(1)),(wz(1,10),ws10(1)),(wz(1,11),ws11(1))
  
        equivalence (wz(1,12),ws12(1)),
     1      (wz(1,13),ws13(1)),(wz(1,14),ws14(1)),(wz(1,15),ws15(1)),
     2      (wz(1,16),ws16(1)),(wz(1,17),ws17(1)),(wz(1,18),ws18(1)),
     3       (wz(1,19),ws19(1)),(wz(1,20),ws20(1)),(wz(1,21),ws21(1))
c 
        data xs1/
     1  -.9936276387206577E+00,-.9577760627534585E+00,
     2  -.8784019728422140E+00,-.7518955413717283E+00,
     3  -.5818084130240981E+00,-.3767632119770460E+00,
     4  -.1488898526924748E+00,0.8759693162605653E-01,
     5  0.3178026002135494E+00,0.5275532658956055E+00,
     6  0.7047665593395038E+00,0.8408050969968233E+00,
     7  0.9318686064230692E+00,0.9805671718763300E+00,
     8  0.9975345390745345E+00,0.9999234647643778E+00,
     9  0.9999945182631860E+00,0.9999971198886042E+00,
     *  0.9999976577925320E+00,0.9999978848890243E+00,
     *  0.9999980241957362E+00,0.9999981329822452E+00,
     *  0.9999982292673497E+00,0.9999983187027539E+00,
     *  0.9999984031128324E+00,0.9999984831825173E+00,
     *  0.9999985592500948E+00,0.9999986315525253E+00,
     *  0.9999987003051216E+00,0.9999987657278877E+00,
     *  0.9999988280964638E+00,0.9999988879132662E+00,
     *  0.9999989464204362E+00,0.9999990070812965E+00,
     *  0.9999990791745548E+00,0.9999991854816008E+00,
     *  0.9999993837275762E+00,0.9999997805641019E+00/
c 
        data xs2/
     1  -.9937265285864976E+00,-.9582578965820439E+00,
     2  -.8794845090876279E+00,-.7536990159595508E+00,
     3  -.5843603208372796E+00,-.3800023520391684E+00,
     4  -.1526739972993156E+00,0.8347547497294627E-01,
     5  0.3135951169120368E+00,0.5235280687699670E+00,
     6  0.7011824173862848E+00,0.8378867352542503E+00,
     7  0.9297816184482763E+00,0.9793777060299740E+00,
     8  0.9971205602635098E+00,0.9998771766347572E+00,
     9  0.9999891912285455E+00,0.9999942453367361E+00,
     *  0.9999953158225403E+00,0.9999957696043305E+00,
     *  0.9999960482319259E+00,0.9999962658780905E+00,
     *  0.9999964585452574E+00,0.9999966375262499E+00,
     *  0.9999968064567486E+00,0.9999969666962327E+00,
     *  0.9999971189183920E+00,0.9999972635977550E+00,
     *  0.9999974011616688E+00,0.9999975320453185E+00,
     *  0.9999976567982858E+00,0.9999977764276364E+00,
     *  0.9999978934223106E+00,0.9999980147183412E+00,
     *  0.9999981588600646E+00,0.9999983712006291E+00,
     *  0.9999987659800440E+00,0.9999995544175486E+00/
  
c 
        data xs3/
     1  -.9937921820877311E+00,-.9585801723651681E+00,
     2  -.8802129845188185E+00,-.7549185693734245E+00,
     3  -.5860928164981944E+00,-.3822086043884352E+00,
     4  -.1552589471987130E+00,0.8065200815531393E-01,
     5  0.3107034617505422E+00,0.5207508680840618E+00,
     6  0.6986970461372788E+00,0.8358490109301442E+00,
     7  0.9283084398039709E+00,0.9785183372848445E+00,
     8  0.9967960362031651E+00,0.9998229160298158E+00,
     9  0.9999792154850078E+00,0.9999885672959035E+00,
     *  0.9999906538707728E+00,0.9999915501017451E+00,
     *  0.9999921037928645E+00,0.9999925377050713E+00,
     *  0.9999929223454458E+00,0.9999932798525641E+00,
     *  0.9999936173471820E+00,0.9999939374890503E+00,
     *  0.9999942416232022E+00,0.9999945306992426E+00,
     *  0.9999948055575307E+00,0.9999950670607194E+00,
     *  0.9999953163091715E+00,0.9999955553269186E+00,
     *  0.9999957891337643E+00,0.9999960316920571E+00,
     *  0.9999963202268420E+00,0.9999967454357768E+00,
     *  0.9999975343190062E+00,0.9999991017706977E+00/
c 
        data xs4/
     1  -.9938631753728774E+00,-.9589370910951307E+00,
     2  -.8810355844232566E+00,-.7563168006478338E+00,
     3  -.5881032981082343E+00,-.3847942683974247E+00,
     4  -.1583141930938392E+00,0.7728899634638919E-01,
     5  0.3072333040101938E+00,0.5173923592961265E+00,
     6  0.6956668104441697E+00,0.8333422019638379E+00,
     7  0.9264769446721473E+00,0.9774334045578457E+00,
     8  0.9963705860451841E+00,0.9997402555573445E+00,
     9  0.9999602910261882E+00,0.9999773038487387E+00,
     *  0.9999813398570294E+00,0.9999831074948349E+00,
     *  0.9999842120906755E+00,0.9999850847968592E+00,
     *  0.9999858612236410E+00,0.9999865832035872E+00,
     *  0.9999872640856313E+00,0.9999879089654680E+00,
     *  0.9999885205271808E+00,0.9999891007047137E+00,
     *  0.9999896511894635E+00,0.9999901736925287E+00,
     *  0.9999906703998827E+00,0.9999911453406247E+00,
     *  0.9999916086680278E+00,0.9999920889333249E+00,
     *  0.9999926620151922E+00,0.9999935106293376E+00,
     *  0.9999950871871219E+00,0.9999982048734106E+00/
  
c 
        data xs5/
     1  -.9939311121944537E+00,-.9592814146909195E+00,
     2  -.8818334468147104E+00,-.7576780568449674E+00,
     3  -.5900662404708097E+00,-.3873251205203887E+00,
     4  -.1613122847120207E+00,0.7397924455300586E-01,
     5  0.3038056351481369E+00,0.5140591989824403E+00,
     6  0.6926402878831932E+00,0.8308161410950181E+00,
     7  0.9246060743474787E+00,0.9762967779615801E+00,
     8  0.9958943233358668E+00,0.9996246883710900E+00,
     9  0.9999255882372780E+00,0.9999553045994443E+00,
     *  0.9999629114778381E+00,0.9999663350410171E+00,
     *  0.9999685074080218E+00,0.9999702380193576E+00,
     *  0.9999717828539318E+00,0.9999732210304255E+00,
     *  0.9999745779389382E+00,0.9999758633808820E+00,
     *  0.9999770825926908E+00,0.9999782393881160E+00,
     *  0.9999793371244039E+00,0.9999803792197623E+00,
     *  0.9999813700868039E+00,0.9999823179543648E+00,
     *  0.9999832436778171E+00,0.9999842058705610E+00,
     *  0.9999853596374655E+00,0.9999870767759620E+00,
     *  0.9999902707581809E+00,0.9999964949185318E+00/
  
c 
        data xs6/
     1  -.9940094720205255E+00,-.9596835812233729E+00,
     2  -.8827735812174164E+00,-.7592918732932890E+00,
     3  -.5924036782538411E+00,-.3903491018384060E+00,
     4  -.1649050473733940E+00,0.7000170537510808E-01,
     5  0.2996736535095996E+00,0.5100266641309108E+00,
     6  0.6889625897413916E+00,0.8277289763497702E+00,
     7  0.9223001681668552E+00,0.9748728973345024E+00,
     8  0.9952697944968463E+00,0.9994481588587729E+00,
     9  0.9998612650682348E+00,0.9999122262535858E+00,
     *  0.9999263820259698E+00,0.9999329639674157E+00,
     *  0.9999372196978444E+00,0.9999406443133068E+00,
     *  0.9999437136421161E+00,0.9999465751669305E+00,
     *  0.9999492765624396E+00,0.9999518364923790E+00,
     *  0.9999542650989872E+00,0.9999565698622808E+00,
     *  0.9999587574243731E+00,0.9999608346122366E+00,
     *  0.9999628103595481E+00,0.9999647015564055E+00,
     *  0.9999665512987659E+00,0.9999684806605103E+00,
     *  0.9999708084782630E+00,0.9999742948777082E+00,
     *  0.9999807871027022E+00,0.9999931962216935E+00/
  
  
c 
        data xs7/
     1  -.9940876714281036E+00,-.9600925407969509E+00,
     2  -.8837427590828852E+00,-.7609722699137119E+00,
     3  -.5948565123301858E+00,-.3935431072306426E+00,
     4  -.1687227146065291E+00,0.6574938564675126E-01,
     5  0.2952268862310121E+00,0.5056538970457454E+00,
     6  0.6849382614409014E+00,0.8243116479507027E+00,
     7  0.9197053576528591E+00,0.9732237575100160E+00,
     8  0.9944943738492455E+00,0.9991838371219284E+00,
     9  0.9997436259050705E+00,0.9998284439210827E+00,
     *  0.9998542249385376E+00,0.9998667146564162E+00,
     *  0.9998749928512511E+00,0.9998817424543982E+00,
     *  0.9998878235702182E+00,0.9998935039770771E+00,
     *  0.9998988711437482E+00,0.9999039599411126E+00,
     *  0.9999087897292682E+00,0.9999133750242801E+00,
     *  0.9999177288722323E+00,0.9999218648672272E+00,
     *  0.9999258011506542E+00,0.9999295726846584E+00,
     *  0.9999332694981597E+00,0.9999371445294043E+00,
     *  0.9999418602254688E+00,0.9999489860224848E+00,
     *  0.9999622756338518E+00,0.9999869755146609E+00/
c 
        data xs8/
     1  -.9942108918517062E+00,-.9608198220542924E+00,
     2  -.8856188418816102E+00,-.7644177494483823E+00,
     3  -.6000877586003446E+00,-.4005429815190042E+00,
     4  -.1772486730550641E+00,0.5613103320490098E-01,
     5  0.2850920120739464E+00,0.4956642752468885E+00,
     6  0.6757826959392237E+00,0.8166389609181635E+00,
     7  0.9140319211124259E+00,0.9697821385134235E+00,
     8  0.9929952617247668E+00,0.9987138095989345E+00,
     9  0.9995265868785750E+00,0.9996668025813875E+00,
     *  0.9997126519364632E+00,0.9997359414864495E+00,
     *  0.9997518844850533E+00,0.9997650733464338E+00,
     *  0.9997770021621151E+00,0.9997881549464029E+00,
     *  0.9997987015353034E+00,0.9998087150422466E+00,
     *  0.9998182373948572E+00,0.9998272998537738E+00,
     *  0.9998359311763148E+00,0.9998441627418884E+00,
     *  0.9998520363854007E+00,0.9998596315338242E+00,
     *  0.9998671478873040E+00,0.9998751356836482E+00,
     *  0.9998850318875192E+00,0.9999002921062812E+00,
     *  0.9999291944817752E+00,0.9999791547857935E+00/
  
  
c 
        data xs9/
     1  -.9923550799323729E+00,-.9519829684542604E+00,
     2  -.8661533913502333E+00,-.7325643821184253E+00,
     3  -.5557881053052802E+00,-.3452528021546931E+00,
     4  -.1137418381879663E+00,0.1240539938570764E+00,
     5  0.3529790657917615E+00,0.5588714958981083E+00,
     6  0.7299829949337961E+00,0.8583648082480122E+00,
     7  0.9413146656396554E+00,0.9831117912562826E+00,
     8  0.9965540898280980E+00,0.9988557137501026E+00,
     9  0.9992417151785704E+00,0.9993602987267365E+00,
     *  0.9994170941399239E+00,0.9994546056309790E+00,
     *  0.9994853507478862E+00,0.9995131562827522E+00,
     *  0.9995391440735357E+00,0.9995636678695372E+00,
     *  0.9995868728406359E+00,0.9996088461345456E+00,
     *  0.9996296568217138E+00,0.9996493669609341E+00,
     *  0.9996680359460766E+00,0.9996857261875953E+00,
     *  0.9997025203137491E+00,0.9997185806192286E+00,
     *  0.9997343493959785E+00,0.9997511487697002E+00,
     *  0.9997725749287099E+00,0.9998072718869373E+00,
     *  0.9998750103280831E+00,0.9999791258620829E+00/
  
c 
        data xs10/
     1  -.9943437382831113E+00,-.9614560812130759E+00,
     2  -.8870140321967268E+00,-.7666947409596947E+00,
     3  -.6032708353845828E+00,-.4045771273598358E+00,
     4  -.1820131067870816E+00,0.5080948512004779E-01,
     5  0.2794241350271587E+00,0.4898844781121257E+00,
     6  0.6701392447415406E+00,0.8113908163355918E+00,
     7  0.9094498274312783E+00,0.9661335219143279E+00,
     8  0.9904413905516625E+00,0.9970981258549994E+00,
     9  0.9984132110432851E+00,0.9987513777358334E+00,
     *  0.9988887767639638E+00,0.9989678287610273E+00,
     *  0.9990266291685652E+00,0.9990774820142523E+00,
     *  0.9991243001428229E+00,0.9991683199377879E+00,
     *  0.9992099811464181E+00,0.9992494893519639E+00,
     *  0.9992869779663531E+00,0.9993225566247084E+00,
     *  0.9993563296271504E+00,0.9993884136707605E+00,
     *  0.9994189870026118E+00,0.9994484408166169E+00,
     *  0.9994778256651315E+00,0.9995099939708789E+00,
     *  0.9995518987946417E+00,0.9996191032610409E+00,
     *  0.9997435551891258E+00,0.9999277764207055E+00/
  
  
c 
        data xs11/
     1  -.9945463286291771E+00,-.9625463126138674E+00,
     2  -.8896347285502771E+00,-.7712667245028411E+00,
     3  -.6099587810225164E+00,-.4132909423582156E+00,
     4  -.1924352204257886E+00,0.3917957377546610E-01,
     5  0.2672161570684253E+00,0.4778007942451092E+00,
     6  0.6588977045009168E+00,0.8016603770374786E+00,
     7  0.9017624439076271E+00,0.9607453653149618E+00,
     8  0.9871342606931746E+00,0.9951220113818284E+00,
     9  0.9970022106699985E+00,0.9975604463124257E+00,
     *  0.9978058990408760E+00,0.9979543684464207E+00,
     *  0.9980682969265737E+00,0.9981682260250631E+00,
     *  0.9982606953529003E+00,0.9983477940599141E+00,
     *  0.9984302874799402E+00,0.9985085509343647E+00,
     *  0.9985828393263891E+00,0.9986533672955782E+00,
     *  0.9987203424026145E+00,0.9987840050837336E+00,
     *  0.9988447385968366E+00,0.9989034092473066E+00,
     *  0.9989623644706180E+00,0.9990279272555445E+00,
     *  0.9991152780370764E+00,0.9992580443828869E+00,
     *  0.9995209803119960E+00,0.9998781939827721E+00/
  
c 
        data xs12/
     1  -.9947462251910808E+00,-.9636573021810650E+00,
     2  -.8923670115978068E+00,-.7761167586250149E+00,
     3  -.6171587743424451E+00,-.4228039149936738E+00,
     4  -.2039775455604722E+00,0.2609968670064289E-01,
     5  0.2532509361824677E+00,0.4637101363432570E+00,
     6  0.6454902285942847E+00,0.7897196216477991E+00,
     7  0.8919410437963167E+00,0.9534013977437996E+00,
     8  0.9821115492907369E+00,0.9916896040946601E+00,
     9  0.9943417544896787E+00,0.9952454957040016E+00,
     *  0.9956775685781221E+00,0.9959541910805212E+00,
     *  0.9961738833203813E+00,0.9963694375506362E+00,
     *  0.9965513130040102E+00,0.9967229207694905E+00,
     *  0.9968855715682330E+00,0.9970399495911233E+00,
     *  0.9971865418165082E+00,0.9973257693821704E+00,
     *  0.9974580484311548E+00,0.9975838783160669E+00,
     *  0.9977040906647343E+00,0.9978206162747039E+00,
     *  0.9979387010987864E+00,0.9980722851569723E+00,
     *  0.9982542463994990E+00,0.9985563507801458E+00,
     *  0.9991050484681168E+00,0.9997909537418466E+00/
  
c 
        data xs13/
     1  -.9946925353970199E+00,-.9634057674783317E+00,
     2  -.8918408546814951E+00,-.7753134529052547E+00,
     3  -.6161289332627012E+00,-.4216371344956642E+00,
     4  -.2027920144349129E+00,0.2716523269334514E-01,
     5  0.2540426900533688E+00,0.4640608893591787E+00,
     6  0.6452155148450656E+00,0.7886057661236467E+00,
     7  0.8897242733736608E+00,0.9497828021505381E+00,
     8  0.9771375416588013E+00,0.9863670351941662E+00,
     9  0.9893514950249045E+00,0.9905932567725850E+00,
     *  0.9912942469400521E+00,0.9918058000994279E+00,
     *  0.9922432934747021E+00,0.9926438899220342E+00,
     *  0.9930195386689817E+00,0.9933745887112670E+00,
     *  0.9937111040694533E+00,0.9940304178320343E+00,
     *  0.9943336006694945E+00,0.9946216265052743E+00,
     *  0.9948954796847471E+00,0.9951563469951967E+00,
     *  0.9954061606349975E+00,0.9956492582924925E+00,
     *  0.9958972775320720E+00,0.9961813669768351E+00,
     *  0.9965749925069418E+00,0.9972314317362285E+00,
     *  0.9983596412588343E+00,0.9996335718089377E+00/
  
c 
        data xs14/
     1  -.9952408607826362E+00,-.9664297400588132E+00,
     2  -.8991800925170527E+00,-.7881537464792341E+00,
     3  -.6349224268849937E+00,-.4461393548752602E+00,
     4  -.2321508180554623E+00,-.5703639324382115E-02,
     5  0.2193728738766498E+00,0.4295099231127226E+00,
     6  0.6127286927944104E+00,0.7599150655430973E+00,
     7  0.8660388140680244E+00,0.9314362166186581E+00,
     8  0.9633771385157297E+00,0.9755925828105925E+00,
     9  0.9801373802215199E+00,0.9822205275818690E+00,
     *  0.9834609762780271E+00,0.9843904823613656E+00,
     *  0.9851947755996135E+00,0.9859352240934953E+00,
     *  0.9866316045461008E+00,0.9872909705100530E+00,
     *  0.9879165806331991E+00,0.9885105246264032E+00,
     *  0.9890745024538548E+00,0.9896101031017075E+00,
     *  0.9901190463424867E+00,0.9906037626785551E+00,
     *  0.9910691403898624E+00,0.9915276715111689E+00,
     *  0.9920127842807169E+00,0.9926055032455372E+00,
     *  0.9934767822379718E+00,0.9949467586533070E+00,
     *  0.9972770395700312E+00,0.9994623015384922E+00/
c 
        data xs15/
     1  -.9951485029642562E+00,-.9659395366616788E+00,
     2  -.8980907501615394E+00,-.7865046068943509E+00,
     3  -.6329956821822679E+00,-.4444162294563212E+00,
     4  -.2312523071477680E+00,-.6331082609793426E-02,
     5  0.2164932873169934E+00,0.4236807099793714E+00,
     6  0.6033434453420667E+00,0.7465587115208852E+00,
     7  0.8487208739226384E+00,0.9111669864655627E+00,
     8  0.9425770831087852E+00,0.9562093367053239E+00,
     9  0.9622935207719547E+00,0.9655312205792379E+00,
     *  0.9676857512960093E+00,0.9694197482724014E+00,
     *  0.9709696210879875E+00,0.9724129077231185E+00,
     *  0.9737753416274668E+00,0.9750669919427002E+00,
     *  0.9762932077453887E+00,0.9774578342638769E+00,
     *  0.9785641988622387E+00,0.9796155575886436E+00,
     *  0.9806157362411885E+00,0.9815709097281997E+00,
     *  0.9824949979586688E+00,0.9834250858295759E+00,
     *  0.9844580581073836E+00,0.9858134051488794E+00,
     *  0.9879157392624476E+00,0.9913831421267862E+00,
     *  0.9960340327313129E+00,0.9993532169693127E+00/
  
c 
        data xs16/
     1  -.9957652698367825E+00,-.9694470819817065E+00,
     2  -.9067589659025056E+00,-.8018649143634948E+00,
     3  -.6557386950825898E+00,-.4744421728683804E+00,
     4  -.2677699773364574E+00,-.4798486987001471E-01,
     5  0.1714582000191487E+00,0.3771995176039216E+00,
     6  0.5572415516298177E+00,0.7021936554207311E+00,
     7  0.8066589805349920E+00,0.8712753645981968E+00,
     8  0.9047829471202209E+00,0.9206044465038019E+00,
     9  0.9286623129829161E+00,0.9336064043592189E+00,
     *  0.9373111450522966E+00,0.9404966664952212E+00,
     *  0.9434179719049447E+00,0.9461624333288656E+00,
     *  0.9487626224824830E+00,0.9512336291142443E+00,
     *  0.9535847849368415E+00,0.9558234135463450E+00,
     *  0.9579562362542118E+00,0.9599903796878284E+00,
     *  0.9619352506043666E+00,0.9638076548835695E+00,
     *  0.9656464924962898E+00,0.9675520537670073E+00,
     *  0.9697700292331899E+00,0.9728145245147140E+00,
     *  0.9775818993083512E+00,0.9849605499608842E+00,
     *  0.9935622842681762E+00,0.9989861519593937E+00/
  
c 
        data xs17/
     1  -.9961762757095662E+00,-.9719518404060261E+00,
     2  -.9133168431499452E+00,-.8141156812923794E+00,
     3  -.6748383080114475E+00,-.5010144667103327E+00,
     4  -.3019045358839216E+00,-.8928422903057698E-01,
     5  0.1237668520233284E+00,0.3240767788747763E+00,
     6  0.4995736532778583E+00,0.6405235982839876E+00,
     7  0.7412598723804273E+00,0.8032155082369465E+00,
     8  0.8366235077674007E+00,0.8544023849226878E+00,
     9  0.8650992673403211E+00,0.8728431729319810E+00,
     *  0.8793663251673072E+00,0.8853029359389051E+00,
     *  0.8908663014406355E+00,0.8961310347944500E+00,
     *  0.9011287756132500E+00,0.9058777840490047E+00,
     *  0.9103919206142581E+00,0.9146835512021474E+00,
     *  0.9187650434331374E+00,0.9226511387435072E+00,
     *  0.9263655793751797E+00,0.9299607005407486E+00,
     *  0.9335713277274935E+00,0.9375366940547329E+00,
     *  0.9425879500568732E+00,0.9500054591728873E+00,
     *  0.9613746872083659E+00,0.9765552560469426E+00,
     *  0.9908865879932709E+00,0.9986359718645417E+00/
  
  
c 
        data xs18/
     1  -.9964532968639627E+00,-.9738190655193993E+00,
     2  -.9187299959218205E+00,-.8252229621914907E+00,
     3  -.6937084011907535E+00,-.5294604920743756E+00,
     4  -.3413769112769710E+00,-.1408369847357894E+00,
     5  0.5942473142140758E-01,0.2464434336713410E+00,
     6  0.4082069964897235E+00,0.5352591002431355E+00,
     7  0.6236003612521923E+00,0.6780650782922568E+00,
     8  0.7101614869879078E+00,0.7305557170065495E+00,
     9  0.7456528545235599E+00,0.7584571286664115E+00,
     *  0.7701317451513190E+00,0.7810821466051991E+00,
     *  0.7914538124339356E+00,0.8013090118692716E+00,
     *  0.8106832291174453E+00,0.8196028622394420E+00,
     *  0.8280911634429216E+00,0.8361711539459526E+00,
     *  0.8438695096564794E+00,0.8512266420456396E+00,
     *  0.8583259391646507E+00,0.8653745771879695E+00,
     *  0.8728937071626135E+00,0.8820336963658173E+00,
     *  0.8948404644756860E+00,0.9139311336638902E+00,
     *  0.9399358912989042E+00,0.9676753508836247E+00,
     *  0.9886579740755739E+00,0.9984132213124013E+00/
c 
        data xs19/
     1  -.9969451346602584E+00,-.9770547239831732E+00,
     2  -.9278428145392844E+00,-.8434629363351275E+00,
     3  -.7241593333702518E+00,-.5749661765001927E+00,
     4  -.4046020337654973E+00,-.2244771183473364E+00,
     5  -.4769907909962423E-01,0.1121656728308584E+00,
     6  0.2432917716269269E+00,0.3396431311319670E+00,
     7  0.4045462846494412E+00,0.4479247735067176E+00,
     8  0.4795459606535646E+00,0.5055799287333367E+00,
     9  0.5289173281529673E+00,0.5506746606751576E+00,
     *  0.5712595663448319E+00,0.5908343810630909E+00,
     *  0.6094800972488794E+00,0.6272505544774494E+00,
     *  0.6441902443707064E+00,0.6603405767737034E+00,
     *  0.6757437555284533E+00,0.6904502209531762E+00,
     *  0.7045395410037567E+00,0.7181789453821585E+00,
     *  0.7317755573438753E+00,0.7463068201584029E+00,
     *  0.7637998400584640E+00,0.7875719602184501E+00,
     *  0.8213881920397073E+00,0.8659956258081827E+00,
     *  0.9150815877584890E+00,0.9578647391128585E+00,
     *  0.9860372460444260E+00,0.9981141664118549E+00/
c 
        data xs20/
     1  -.9976818483891922E+00,-.9821870576503572E+00,
     2  -.9430653646289656E+00,-.8754274617560973E+00,
     3  -.7801597281510972E+00,-.6630530054286641E+00,
     4  -.5338256092202240E+00,-.4047341689738038E+00,
     5  -.2879934438824555E+00,-.1916932604590108E+00,
     6  -.1165937293079858E+00,-.5747648129548355E-01,
     7  -.7885622450616081E-02,0.3657250266821488E-01,
     8  0.7803637500045357E-01,0.1173854111815680E+00,
     9  0.1549704908967926E+00,0.1909465160688846E+00,
     *  0.2253998637254528E+00,0.2583939425989337E+00,
     *  0.2899847365667620E+00,0.3202267669554510E+00,
     *  0.3491784369028732E+00,0.3769129044040864E+00,
     *  0.4035454929765334E+00,0.4293057869506920E+00,
     *  0.4547185972660663E+00,0.4809958581145886E+00,
     *  0.5106403999410074E+00,0.5478167184168019E+00,
     *  0.5976268006995429E+00,0.6633484684744027E+00,
     *  0.7421260144654114E+00,0.8241393770178469E+00,
     *  0.8972864155453351E+00,0.9522036557241416E+00,
     *  0.9849352341055069E+00,0.9980415331510295E+00/
c 
        data xs21/
     1  -.9992852891508569E+00,-.9942324147033921E+00,
     2  -.9808907178887876E+00,-.9573610497647588E+00,
     3  -.9241053764472341E+00,-.8830803672745517E+00,
     4  -.8366942010643945E+00,-.7870952922912594E+00,
     5  -.7358903736863722E+00,-.6841534469419811E+00,
     6  -.6325532299919769E+00,-.5814921413043821E+00,
     7  -.5312149844744460E+00,-.4818799216322492E+00,
     8  -.4336013713610688E+00,-.3864701984727568E+00,
     9  -.3405565263607053E+00,-.2958982882118089E+00,
     *  -.2524868333934854E+00,-.2102521934698133E+00,
     *  -.1690363854021166E+00,-.1285191812889073E+00,
     *  -.8805232853035005E-01,-.4636344555768934E-01,
     *  -.1179946085688277E-02,0.5096383080446731E-01,
     *  0.1143556270251251E+00,0.1929946249670507E+00,
     *  0.2887811159777554E+00,0.3996912539843536E+00,
     *  0.5194899867478533E+00,0.6393702199073610E+00,
     *  0.7501707601198881E+00,0.8441114317167844E+00,
     *  0.9159432391467803E+00,0.9637326340439623E+00,
     *  0.9893889789097698E+00,0.9987131139291487E+00/
c 
c 
c 
c 
c 
        data ws1/
     1  0.1780554089473873E-01,0.5626793826380155E-01,
     2  0.1030203947442069E+00,0.1493680271295429E+00,
     3  0.1893438353654779E+00,0.2187011602093503E+00,
     4  0.2346527909220476E+00,0.2358145757554336E+00,
     5  0.2222115326376198E+00,0.1952537299436553E+00,
     6  0.1577105504278984E+00,0.1137073796362089E+00,
     7  0.6882367753958385E-01,0.3034154580787351E-01,
     8  0.6670473912774839E-02,0.2963322494724191E-03,
     9  0.6712363965049403E-05,0.9181719588695246E-06,
     *  0.3150213121151208E-06,0.1676668826499471E-06,
     *  0.1192080037244396E-06,0.1010092105976836E-06,
     *  0.9238196592590891E-07,0.8675630064134275E-07,
     *  0.8216738774790661E-07,0.7802495508802849E-07,
     *  0.7414855200819412E-07,0.7049194235137905E-07,
     *  0.6704925940127632E-07,0.6383975692285556E-07,
     *  0.6097401025159880E-07,0.5884476091929991E-07,
     *  0.5869026055353415E-07,0.6406351921037089E-07,
     *  0.8363750740724121E-07,0.1376520069125081E-06,
     *  0.2815032683675823E-06,0.4649196529438645E-06/
  
c 
        data ws2/
     1  0.1755903587295352E-01,0.5576038474423054E-01,
     2  0.1023435801067762E+00,0.1486182487817748E+00,
     3  0.1886115937131891E+00,0.2180725890046490E+00,
     4  0.2342023499051840E+00,0.2355978083700877E+00,
     5  0.2222590842489799E+00,0.1955691096489066E+00,
     6  0.1582715829969910E+00,0.1144682450800899E+00,
     7  0.6970977822080610E-01,0.3122068774629705E-01,
     8  0.7288032850928180E-02,0.4278785311815740E-03,
     9  0.1277234937034840E-04,0.1820724804417856E-05,
     *  0.6289021128941663E-06,0.3352695399257525E-06,
     *  0.2384722427010284E-06,0.2021051951919161E-06,
     *  0.1848695104080232E-06,0.1736250657346714E-06,
     *  0.1644411736921891E-06,0.1561433230315552E-06,
     *  0.1483779910035953E-06,0.1410513677233878E-06,
     *  0.1341475907682709E-06,0.1277063782762751E-06,
     *  0.1219532148081083E-06,0.1176765082149130E-06,
     *  0.1173559372472225E-06,0.1281004303438908E-06,
     *  0.1671859270522619E-06,0.2746929244027579E-06,
     *  0.5594207941153619E-06,0.9293866672544108E-06/
  
c 
        data ws3/
     1  0.1739493254980569E-01,0.5541867976266221E-01,
     2  0.1018841917203069E+00,0.1481058890989506E+00,
     3  0.1881083102650603E+00,0.2176377858161460E+00,
     4  0.2338873441152663E+00,0.2354408925407330E+00,
     5  0.2222816636149487E+00,0.1957745265219875E+00,
     6  0.1586462276822697E+00,0.1149825598651854E+00,
     7  0.7031382358621541E-01,0.3182376219479303E-01,
     8  0.7718911094142454E-02,0.5438554441372797E-03,
     9  0.2270662969181585E-04,0.3515589586806957E-05,
     *  0.1237289305849427E-05,0.6647389274448243E-06,
     *  0.4748813546184736E-06,0.4032962718301273E-06,
     *  0.3692015627816323E-06,0.3468553616847725E-06,
     *  0.3285344633164790E-06,0.3119612602834876E-06,
     *  0.2964613428205426E-06,0.2818290632565819E-06,
     *  0.2680285088083929E-06,0.2551515675005835E-06,
     *  0.2436536595815826E-06,0.2351308556942595E-06,
     *  0.2345834685321605E-06,0.2562871563044945E-06,
     *  0.3347974746890245E-06,0.5498763646382468E-06,
     *  0.1115465426779510E-05,0.1850944392744276E-05/
  
c 
        data ws4/
     1  0.1721589469571120E-01,0.5503236245616672E-01,
     2  0.1013513924146218E+00,0.1474996002606014E+00,
     3  0.1875024824939856E+00,0.2171048137418702E+00,
     4  0.2334906035477717E+00,0.2352286672237584E+00,
     5  0.2222826142979411E+00,0.1959964776922170E+00,
     6  0.1590772591492479E+00,0.1155913861419543E+00,
     7  0.7104287459072084E-01,0.3256247010810779E-01,
     8  0.8254759455060500E-02,0.6996799064521617E-03,
     9  0.3955100104607514E-04,0.6713341378471961E-05,
     *  0.2425243620696782E-05,0.1319765259243961E-05,
     *  0.9520460601832518E-06,0.8132073468291260E-06,
     *  0.7456977073354130E-06,0.7002125875655045E-06,
     *  0.6623228610388351E-06,0.6278620440667106E-06,
     *  0.5955769062924691E-06,0.5650564983178180E-06,
     *  0.5361919807804738E-06,0.5091596626825086E-06,
     *  0.4848664647879775E-06,0.4665258518442602E-06,
     *  0.4644253293147899E-06,0.5078666316037087E-06,
     *  0.6665879917014596E-06,0.1099017566955159E-05,
     *  0.2226905618385572E-05,0.3675474308287722E-05/
c 
        data ws5/
     1  0.1704400926124358E-01,0.5465736113886843E-01,
     2  0.1008311277801279E+00,0.1469049761186202E+00,
     3  0.1869058538674278E+00,0.2165769095480494E+00,
     4  0.2330932650101198E+00,0.2350091330226820E+00,
     5  0.2222694128375653E+00,0.1961986124708244E+00,
     6  0.1594853685032435E+00,0.1161779032507230E+00,
     7  0.7175481570235444E-01,0.3329455231374652E-01,
     8  0.8799750810843658E-02,0.8787656935998126E-03,
     9  0.6567440865460299E-04,0.1243481392918593E-04,
     *  0.4655132007530153E-05,0.2580311011729228E-05,
     *  0.1882419532574043E-05,0.1616263095680366E-05,
     *  0.1484856531311648E-05,0.1395210846350548E-05,
     *  0.1320092611298665E-05,0.1251622163756878E-05,
     *  0.1187424138181181E-05,0.1126719234865374E-05,
     *  0.1069310316595316E-05,0.1015576997987754E-05,
     *  0.9674047294040245E-06,0.9314340051613819E-06,
     *  0.9288411676220934E-06,0.1019612078660681E-05,
     *  0.1345510271300212E-05,0.2227449280367118E-05,
     *  0.4501234908543689E-05,0.7226945198243443E-05/
c 
        data ws6/
     1  0.1684477149640810E-01,0.5421504826450484E-01,
     2  0.1002112951004236E+00,0.1461917427219979E+00,
     3  0.1861864897505282E+00,0.2159367943062183E+00,
     4  0.2326069544784523E+00,0.2347337022334900E+00,
     5  0.2222398680320935E+00,0.1964272658048288E+00,
     6  0.1599631439764297E+00,0.1168739772218132E+00,
     7  0.7260709342751738E-01,0.3417689952669246E-01,
     8  0.9464071201914158E-02,0.1112345865864462E-02,
     9  0.1069217425172609E-03,0.2267502379917819E-04,
     *  0.8851184488275595E-05,0.5018470395703330E-05,
     *  0.3711855300882472E-05,0.3207029399783208E-05,
     *  0.2952986773020331E-05,0.2777048320854157E-05,
     *  0.2628579683409580E-05,0.2492892468772692E-05,
     *  0.2365542923252119E-05,0.2245076823128908E-05,
     *  0.2131158974317071E-05,0.2024625201118118E-05,
     *  0.1929427596450596E-05,0.1859379693599250E-05,
     *  0.1858359932117491E-05,0.2049900315408983E-05,
     *  0.2723538777832011E-05,0.4531485491650275E-05,
     *  0.9119103457669171E-05,0.1413835375988345E-04/
  
c 
        data ws7/
     1  0.1664447863431404E-01,0.5375860060234428E-01,
     2  0.9956114008497639E-01,0.1454344754277920E+00,
     3  0.1854144603639509E+00,0.2152410193279825E+00,
     4  0.2320674261177108E+00,0.2344125431959896E+00,
     5  0.2221768262762707E+00,0.1966392777160401E+00,
     6  0.1604454195490522E+00,0.1175988417935685E+00,
     7  0.7351277838428896E-01,0.3513127923905249E-01,
     8  0.1020095884828047E-01,0.1396250988356474E-02,
     9  0.1685765177967673E-03,0.4029219560530727E-04,
     *  0.1655366794616147E-04,0.9668830419223072E-05,
     *  0.7282153790354300E-05,0.6343005988571067E-05,
     *  0.5858049021340586E-05,0.5515609644735886E-05,
     *  0.5224039476117546E-05,0.4956639520335998E-05,
     *  0.4705310811788888E-05,0.4467426190158249E-05,
     *  0.4242480844399329E-05,0.4032363766653391E-05,
     *  0.3845491134144111E-05,0.3710974379694546E-05,
     *  0.3720834293197084E-05,0.4132299312275364E-05,
     *  0.5542500942746327E-05,0.9288161067783519E-05,
     *  0.1857457232353100E-04,0.2738675166089225E-04/
  
c 
        data ws8/
     1  0.1631354749749151E-01,0.5287208052502311E-01,
     2  0.9817255903394821E-01,0.1437190397959946E+00,
     3  0.1835979014341372E+00,0.2135583587493708E+00,
     4  0.2307300161466432E+00,0.2335890339582588E+00,
     5  0.2219816545954789E+00,0.1971295544144253E+00,
     6  0.1616169693842617E+00,0.1193719153068032E+00,
     7  0.7569261143791704E-01,0.3732847711727998E-01,
     8  0.1176582266720547E-01,0.1917092176853284E-02,
     9  0.2684389582973397E-03,0.6995259980859289E-04,
     *  0.3029100870926738E-04,0.1839558249346512E-04,
     *  0.1416674980853001E-04,0.1243051955640330E-04,
     *  0.1149798605223795E-04,0.1083293127020899E-04,
     *  0.1027153773461855E-04,0.9762187893736222E-05,
     *  0.9287606037094150E-05,0.8842024799186338E-05,
     *  0.8425762685523442E-05,0.8043955120324989E-05,
     *  0.7715395991910207E-05,0.7504773154558133E-05,
     *  0.7611063305055074E-05,0.8587640738980822E-05,
     *  0.1173762950279136E-04,0.2008143901301922E-04,
     *  0.4024970508091251E-04,0.4896784053874193E-04/
  
  
c 
        data ws9/
     1  0.2090228328435429E-01,0.6198627950143815E-01,
     2  0.1099850687693039E+00,0.1564016843246290E+00,
     3  0.1955545553845480E+00,0.2233608478996564E+00,
     4  0.2371897879003910E+00,0.2358523957226636E+00,
     5  0.2196143538455974E+00,0.1901891522037731E+00,
     6  0.1506852404275561E+00,0.1055992846919649E+00,
     7  0.6097817883470300E-01,0.2478161278720300E-01,
     8  0.5293537167450298E-02,0.7558533802948669E-03,
     9  0.1856913610580952E-03,0.7567756069415749E-04,
     *  0.4378129188987074E-04,0.3307598552144631E-04,
     *  0.2896674497418608E-04,0.2680455737353023E-04,
     *  0.2522258886290712E-04,0.2384723854984054E-04,
     *  0.2257670875668945E-04,0.2138124553545302E-04,
     *  0.2025041345153598E-04,0.1917968742193576E-04,
     *  0.1816845481756637E-04,0.1722500125106060E-04,
     *  0.1638754238490471E-04,0.1579793381350248E-04,
     *  0.1593724597629532E-04,0.1822031051314459E-04,
     *  0.2597091999793650E-04,0.4660786316602982E-04,
     *  0.9381518462765515E-04,0.7534551005782245E-04/
  
c 
        data ws10/
     1  0.1598384285765025E-01,0.5221628807035989E-01,
     2  0.9733338803707886E-01,0.1428106416049307E+00,
     3  0.1827073100286661E+00,0.2127578097195338E+00,
     4  0.2300788251258254E+00,0.2331323128341695E+00,
     5  0.2217495878282840E+00,0.1971399436556014E+00,
     6  0.1618809291039636E+00,0.1199006582852673E+00,
     7  0.7649670374874748E-01,0.3837695850609888E-01,
     8  0.1284923477659477E-01,0.2653394972784571E-02,
     9  0.5722226724542694E-03,0.1955116312336867E-03,
     *  0.9796703882431916E-04,0.6565146724620650E-04,
     *  0.5376885361462394E-04,0.4850655890215419E-04,
     *  0.4530978287946302E-04,0.4279552685284596E-04,
     *  0.4055899301021668E-04,0.3847900229704718E-04,
     *  0.3651623753012368E-04,0.3465823962022947E-04,
     *  0.3290652795797357E-04,0.3128914869683059E-04,
     *  0.2991627150730285E-04,0.2914869973055332E-04,
     *  0.3005873407675912E-04,0.3537260203675903E-04,
     *  0.5086374990305113E-04,0.8896184171460867E-04,
     *  0.1652633880378697E-03,0.1662380054552571E-03/
  
c 
        data ws11/
     1  0.1545843374848308E-01,0.5097690366334015E-01,
     2  0.9555106282714211E-01,0.1407334998073937E+00,
     3  0.1805954567481908E+00,0.2108561478561429E+00,
     4  0.2285945335611007E+00,0.2322225881260207E+00,
     5  0.2215152803355238E+00,0.1976255088707507E+00,
     6  0.1630718711611800E+00,0.1217099156276383E+00,
     7  0.7872892676905761E-01,0.4066335267290381E-01,
     8  0.1462338867942650E-01,0.3538747717972821E-02,
     9  0.9031912226281995E-03,0.3398042471755166E-03,
     *  0.1803631042107984E-03,0.1257324343215009E-03,
     *  0.1051633281376006E-03,0.9564700817408749E-04,
     *  0.8959674255430824E-04,0.8471565096092635E-04,
     *  0.8033013740950147E-04,0.7623818597774735E-04,
     *  0.7237391162535435E-04,0.6871607948096995E-04,
     *  0.6527268852839662E-04,0.6211178797523216E-04,
     *  0.5948512846707486E-04,0.5820825806262874E-04,
     *  0.6067618925576341E-04,0.7284302946114478E-04,
     *  0.1071150324247881E-03,0.1898542779501641E-03,
     *  0.3419820053524061E-03,0.2947818787203991E-03/
  
c 
        data ws12/
     1  0.1493330461908448E-01,0.4968376112861324E-01,
     2  0.9363964135970539E-01,0.1384540732167444E+00,
     3  0.1782202993145595E+00,0.2086467207538004E+00,
     4  0.2267801233871038E+00,0.2309886935465059E+00,
     5  0.2209965642370670E+00,0.1979016917493623E+00,
     6  0.1641583761959016E+00,0.1235359910042368E+00,
     7  0.8109382845539174E-01,0.4316615552859600E-01,
     8  0.1665553992562550E-01,0.4670770965769521E-02,
     9  0.1396345184829289E-02,0.5800010444284072E-03,
     *  0.3284309572911447E-03,0.2394218509083964E-03,
     *  0.2048085160621534E-03,0.1878235586747034E-03,
     *  0.1764278992874788E-03,0.1669862340464588E-03,
     *  0.1584236313514954E-03,0.1504118094944754E-03,
     *  0.1428419972985902E-03,0.1356812753924409E-03,
     *  0.1289566068625927E-03,0.1228309951700433E-03,
     *  0.1178830688166962E-03,0.1159571434088671E-03,
     *  0.1223691331900766E-03,0.1499903046812143E-03,
     *  0.2252111578152748E-03,0.4024784158203339E-03,
     *  0.6956258404526908E-03,0.5254129562646749E-03/
  
  
c 
        data ws13/
     1  0.1506604662966895E-01,0.4993405875680113E-01,
     2  0.9392590724789091E-01,0.1387134199848888E+00,
     3  0.1784075841071873E+00,0.2087285963448760E+00,
     4  0.2267323761675518E+00,0.2307940080499276E+00,
     5  0.2206415050956907E+00,0.1973720042206668E+00,
     6  0.1634324575093457E+00,0.1225751924822109E+00,
     7  0.7983868193483851E-01,0.4165634700301924E-01,
     8  0.1568022650470914E-01,0.4842806018618718E-02,
     9  0.1772572150906375E-02,0.8775196125029974E-03,
     *  0.5756869941594173E-03,0.4645244451928664E-03,
     *  0.4158885904786121E-03,0.3870517383350155E-03,
     *  0.3648988662263037E-03,0.3455243333101787E-03,
     *  0.3277221990477836E-03,0.3110814595979607E-03,
     *  0.2954450284040133E-03,0.2807680157614659E-03,
     *  0.2671285904612822E-03,0.2549048656290567E-03,
     *  0.2453745632666915E-03,0.2425669139212985E-03,
     *  0.2582832444234444E-03,0.3215399640029041E-03,
     *  0.4904344072513817E-03,0.8654073649432460E-03,
     *  0.1359568056578127E-02,0.9347540859684533E-03/
  
c 
        data ws14/
     1  0.1362760778087279E-01,0.4644880765847716E-01,
     2  0.8889637980654196E-01,0.1328646717238126E+00,
     3  0.1724668574417347E+00,0.2033535121111844E+00,
     4  0.2224770000515529E+00,0.2280899350457065E+00,
     5  0.2197793570542149E+00,0.1984769974547619E+00,
     6  0.1664214454010544E+00,0.1270924226287796E+00,
     7  0.8518411334415840E-01,0.4680056150201374E-01,
     8  0.1954722689351132E-01,0.6993116524475959E-02,
     9  0.2875055433456974E-02,0.1523834493984680E-02,
     *  0.1036371494071272E-02,0.8504594240793408E-03,
     *  0.7671287327999090E-03,0.7166593273271354E-03,
     *  0.6771558858971117E-03,0.6420838359745401E-03,
     *  0.6094742188370945E-03,0.5786933048793533E-03,
     *  0.5495238829521852E-03,0.5219559996314011E-03,
     *  0.4963110989325370E-03,0.4738621366054729E-03,
     *  0.4587698466257771E-03,0.4634154226184871E-03,
     *  0.5196450216726312E-03,0.6931194803422810E-03,
     *  0.1102915746039297E-02,0.1905320573771650E-02,
     *  0.2594053348942965E-02,0.1424060204705048E-02/
c 
        data ws15/
     1  0.1386841552626945E-01,0.4698301103871714E-01,
     2  0.8951837736181688E-01,0.1333210918406974E+00,
     3  0.1725322507476907E+00,0.2028576150549576E+00,
     4  0.2213103635192672E+00,0.2262006023023152E+00,
     5  0.2171690679783378E+00,0.1952026930434567E+00,
     6  0.1626137810948065E+00,0.1230225879425815E+00,
     7  0.8149125891497690E-01,0.4487994713428895E-01,
     8  0.2032248621633931E-01,0.8677982814592047E-02,
     9  0.4213138002554757E-02,0.2528047906090857E-02,
     *  0.1882098910829414E-02,0.1621199416609128E-02,
     *  0.1489894084096102E-02,0.1400467350869191E-02,
     *  0.1325930686995305E-02,0.1258218195972308E-02,
     *  0.1194841010620315E-02,0.1134960579295906E-02,
     *  0.1078302431654196E-02,0.1025027948251094E-02,
     *  0.9762825330740149E-03,0.9361579065797327E-03,
     *  0.9176113193491441E-03,0.9576233667852586E-03,
     *  0.1143555212996233E-02,0.1636323078420324E-02,
     *  0.2683100926301632E-02,0.4272941682869182E-02,
     *  0.4530978802661598E-02,0.1821656703193376E-02/
c 
        data ws16/
     1  0.1222878729390068E-01,0.4285933899253371E-01,
     2  0.8345230042332787E-01,0.1261452290623626E+00,
     3  0.1650750649444811E+00,0.1958646091677884E+00,
     4  0.2154146938167388E+00,0.2218914992564782E+00,
     5  0.2147469193289745E+00,0.1947156428404839E+00,
     6  0.1637815314410169E+00,0.1251906024186955E+00,
     7  0.8379581909199189E-01,0.4693958863336365E-01,
     8  0.2250303601849073E-01,0.1080138870365326E-01,
     9  0.6035337391867486E-02,0.4140806390418536E-02,
     *  0.3378309080463052E-02,0.3030692159214580E-02,
     *  0.2824859444962189E-02,0.2668974457223337E-02,
     *  0.2533757634127642E-02,0.2409759065667808E-02,
     *  0.2293753397659249E-02,0.2184609645422979E-02,
     *  0.2082196927588894E-02,0.1987585280932136E-02,
     *  0.1904746773812646E-02,0.1845953783716174E-02,
     *  0.1847032681890030E-02,0.2002883980245216E-02,
     *  0.2516450960304740E-02,0.3723265114361740E-02,
     *  0.6000179919306654E-02,0.8573267475646098E-02,
     *  0.7743614715571459E-02,0.2865887863529567E-02/
c 
        data ws17/
     1  0.1110713906655937E-01,0.3975843079565211E-01,
     2  0.7850613238519290E-01,0.1197878177344086E+00,
     3  0.1578240779063885E+00,0.1882601015033773E+00,
     4  0.2079769172130297E+00,0.2150634673283194E+00,
     5  0.2088235544499896E+00,0.1897645063012507E+00,
     6  0.1595811045464016E+00,0.1213219530616173E+00,
     7  0.8033079095444842E-01,0.4541370058581692E-01,
     8  0.2370137559439200E-01,0.1325682053399752E-01,
     9  0.8788129778634617E-02,0.6967619109560174E-02,
     *  0.6173571035907297E-02,0.5731057613946892E-02,
     *  0.5406778788864156E-02,0.5127557678985921E-02,
     *  0.4870845408278701E-02,0.4629442101806798E-02,
     *  0.4400871607122594E-02,0.4184418065960579E-02,
     *  0.3980930619576423E-02,0.3794937930582937E-02,
     *  0.3641864093751434E-02,0.3568860549288093E-02,
     *  0.3705525905894905E-02,0.4344100133933261E-02,
     *  0.5972176183926411E-02,0.9157470029434818E-02,
     *  0.1363285047977695E-01,0.1588661879117484E-01,
     *  0.1166562875063144E-01,0.3890811828509624E-02/
  
c 
        data ws18/
     1  0.1032641937356067E-01,0.3725217089059854E-01,
     2  0.7389140613918284E-01,0.1130290142838400E+00,
     3  0.1491009271818227E+00,0.1778930652857377E+00,
     4  0.1963587524293174E+00,0.2025799773112359E+00,
     5  0.1957635912887711E+00,0.1762475690567011E+00,
     6  0.1456566923908020E+00,0.1077050403765863E+00,
     7  0.6981225165616568E-01,0.4120155373757348E-01,
     8  0.2481638752193286E-01,0.1701614007225658E-01,
     9  0.1365079926960017E-01,0.1213350252944680E-01,
     *  0.1127584989110230E-01,0.1064671368698593E-01,
     *  0.1010621937839008E-01,0.9609850256247266E-02,
     *  0.9142901511706828E-02,0.8700199224219550E-02,
     *  0.8280185589703131E-02,0.7884099618800119E-02,
     *  0.7518930073440852E-02,0.7208039776857188E-02,
     *  0.7022045285307819E-02,0.7156398777034097E-02,
     *  0.8068282245601584E-02,0.1055402585461373E-01,
     *  0.1553364190628114E-01,0.2284826667939268E-01,
     *  0.2824905985301587E-01,0.2566425908223427E-01,
     *  0.1549316316325543E-01,0.4602576978640157E-02/
  
  
c 
        data ws19/
     1  0.8948832007955826E-02,0.3300209097371221E-01,
     2  0.6637533469416523E-01,0.1023214760231847E+00,
     3  0.1354286024869706E+00,0.1614882856939804E+00,
     4  0.1773221422694462E+00,0.1807240025999083E+00,
     5  0.1705377849168309E+00,0.1471503963175719E+00,
     6  0.1140134026886696E+00,0.7929758768959450E-01,
     7  0.5230088575372565E-01,0.3613768322397964E-01,
     8  0.2812068236640197E-01,0.2439978449914841E-01,
     9  0.2244363120279341E-01,0.2113213538487230E-01,
     *  0.2006241581764325E-01,0.1909995596263961E-01,
     *  0.1820020137170214E-01,0.1734804225615931E-01,
     *  0.1653818688716414E-01,0.1576940808860308E-01,
     *  0.1504505814827468E-01,0.1438028453177072E-01,
     *  0.1382393368588840E-01,0.1351731483266009E-01,
     *  0.1382983239246585E-01,0.1556487153468846E-01,
     *  0.1998351783177846E-01,0.2824549698937852E-01,
     *  0.3960957370327885E-01,0.4849843597472735E-01,
     *  0.4773829408756206E-01,0.3637388645760142E-01,
     *  0.1970953392722439E-01,0.5517000438710526E-02/
  
c 
        data ws20/
     1  0.6844546598309621E-02,0.2597882773461111E-01,
     2  0.5307248945145074E-01,0.8204002097723424E-01,
     3  0.1075148379709101E+00,0.1250577756165029E+00,
     4  0.1313132625685440E+00,0.1247602254515481E+00,
     5  0.1073064266463237E+00,0.8519349488396223E-01,
     6  0.6596982576988769E-01,0.5339949415170158E-01,
     7  0.4651055569970410E-01,0.4274215074609144E-01,
     8  0.4032021806998282E-01,0.3843134307384658E-01,
     9  0.3676278324796054E-01,0.3520301164948501E-01,
     *  0.3371403577175331E-01,0.3228337271447865E-01,
     *  0.3090722535552069E-01,0.2958640963159920E-01,
     *  0.2832872517865545E-01,0.2715837281997220E-01,
     *  0.2614232962542636E-01,0.2545938407036311E-01,
     *  0.2555715501837472E-01,0.2740591598601404E-01,
     *  0.3258732498931086E-01,0.4265944041623362E-01,
     *  0.5759528618165986E-01,0.7339474319499631E-01,
     *  0.8242534130307047E-01,0.7949562696209854E-01,
     *  0.6519726404589393E-01,0.4396015623295171E-01,
     *  0.2193078547795686E-01,0.5789811587498847E-02/
c 
  
        data ws21/
     1  0.2144890842735288E-02,0.8667248941989610E-02,
     2  0.1833421085850688E-01,0.2863134065480143E-01,
     3  0.3753531368394887E-01,0.4410133349957096E-01,
     4  0.4830865514872050E-01,0.5062196794779610E-01,
     5  0.5161124506583468E-01,0.5175354964255816E-01,
     6  0.5138164852490184E-01,0.5070086477440726E-01,
     7  0.4982780905905044E-01,0.4882355304605747E-01,
     8  0.4771846636634869E-01,0.4653194821231353E-01,
     9  0.4528884266423860E-01,0.4402883551795351E-01,
     *  0.4280540002336474E-01,0.4168906318276632E-01,
     *  0.4079329578439429E-01,0.4034416107335378E-01,
     *  0.4079471266008781E-01,0.4295745811684897E-01,
     *  0.4799301101364086E-01,0.5702766034004356E-01,
     *  0.7045371107051734E-01,0.8717649111975151E-01,
     *  0.1040652482266681E+00,0.1166839564907722E+00,
     *  0.1214060528980976E+00,0.1167987742974570E+00,
     *  0.1034871029206133E+00,0.8350664666934067E-01,
     *  0.5982555602884088E-01,0.3606861971349424E-01,
     *  0.1622893731050887E-01,0.3882416576075897E-02/
  
c 
c       return to the user the nodes and weights of all 38 quadratures
c 
        do 1400 i=1,21
        do 1200 j=1,38
c 
        xsout(j,i)=xz(j,i)
        wsout(j,i)=wz(j,i)
 1200 continue
 1400 continue
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine o5rnear7(xsout,wsout)
        implicit real *8 (a-h,o-z)
        save
        dimension xsout(47,21),wsout(47,21)
c 
        dimension xs1(47),xs2(47),xs3(47),xs4(47),xs5(47),
     1      xs6(47),xs7(47),xs8(47),xs9(47),xs10(47),xs11(47),
     2      xs12(47),xs13(47),xs14(47),xs15(47),xs16(47),xs17(47),
     3      xs18(47),xs19(47),xs20(47),xs21(47),xz(47,42)
c 
        equivalence (xz(1,1),xs1(1)),(xz(1,2),xs2(1)),
     1      (xz(1,3),xs3(1)),(xz(1,4),xs4(1)),(xz(1,5),xs5(1)),
     2      (xz(1,6),xs6(1)),(xz(1,7),xs7(1)),(xz(1,8),xs8(1)),
     3      (xz(1,9),xs9(1)),(xz(1,10),xs10(1)),(xz(1,11),xs11(1))
  
        equivalence (xz(1,12),xs12(1)),
     1      (xz(1,13),xs13(1)),(xz(1,14),xs14(1)),(xz(1,15),xs15(1)),
     2      (xz(1,16),xs16(1)),(xz(1,17),xs17(1)),(xz(1,18),xs18(1)),
     3       (xz(1,19),xs19(1)),(xz(1,20),xs20(1)),(xz(1,21),xs21(1))
        dimension ws1(47),ws2(47),ws3(47),ws4(47),ws5(47),
     1      ws6(47),ws7(47),ws8(47),ws9(47),ws10(47),ws11(47),
     2      ws12(47),ws13(47),ws14(47),ws15(47),ws16(47),ws17(47),
     3      ws18(47),ws19(47),ws20(47),ws21(47),wz(47,42)
c 
        equivalence (wz(1,1),ws1(1)),(wz(1,2),ws2(1)),
     1      (wz(1,3),ws3(1)),(wz(1,4),ws4(1)),(wz(1,5),ws5(1)),
     2      (wz(1,6),ws6(1)),(wz(1,7),ws7(1)),(wz(1,8),ws8(1)),
     3      (wz(1,9),ws9(1)),(wz(1,10),ws10(1)),(wz(1,11),ws11(1))
  
        equivalence (wz(1,12),ws12(1)),
     1      (wz(1,13),ws13(1)),(wz(1,14),ws14(1)),(wz(1,15),ws15(1)),
     2      (wz(1,16),ws16(1)),(wz(1,17),ws17(1)),(wz(1,18),ws18(1)),
     3       (wz(1,19),ws19(1)),(wz(1,20),ws20(1)),(wz(1,21),ws21(1))
  
  
  
        dimension xs22(47),xs23(47),xs24(47),xs25(47),
     1      xs26(47),xs27(47),xs28(47),xs29(47),xs30(47),xs31(47),
     2      xs32(47),xs33(47),xs34(47),xs35(47),xs36(47),xs37(47),
     3      xs38(47),xs39(47),xs40(47),xs41(47),xs42(47)
c 
        equivalence (xz(1,22),xs22(1)),
     1      (xz(1,23),xs23(1)),(xz(1,24),xs24(1)),(xz(1,25),xs25(1)),
     2      (xz(1,26),xs26(1)),(xz(1,27),xs27(1)),(xz(1,28),xs28(1)),
     3       (xz(1,29),xs29(1)),(xz(1,30),xs30(1)),(xz(1,31),xs31(1))
c 
        equivalence (xz(1,32),xs32(1)),(xz(1,42),xs42(1)),
     1      (xz(1,33),xs33(1)),(xz(1,34),xs34(1)),(xz(1,35),xs35(1)),
     2      (xz(1,36),xs36(1)),(xz(1,37),xs37(1)),(xz(1,38),xs38(1)),
     3       (xz(1,39),xs39(1)),(xz(1,40),xs40(1)),(xz(1,41),xs41(1))
  
  
  
  
        dimension ws22(47),ws23(47),ws24(47),ws25(47),
     1      ws26(47),ws27(47),ws28(47),ws29(47),ws30(47),ws31(47),
     2      ws32(47),ws33(47),ws34(47),ws35(47),ws36(47),ws37(47),
     3      ws38(47),ws39(47),ws40(47),ws41(47),ws42(47)
  
  
c 
        equivalence (wz(1,22),ws22(1)),
     1      (wz(1,23),ws23(1)),(wz(1,24),ws24(1)),(wz(1,25),ws25(1)),
     2      (wz(1,26),ws26(1)),(wz(1,27),ws27(1)),(wz(1,28),ws28(1)),
     3      (wz(1,29),ws29(1)),(wz(1,30),ws30(1)),(wz(1,31),ws31(1))
  
c 
        equivalence (wz(1,32),ws32(1)),(wz(1,42),ws42(1)),
     1      (wz(1,33),ws33(1)),(wz(1,34),ws34(1)),(wz(1,35),ws35(1)),
     2      (wz(1,36),ws36(1)),(wz(1,37),ws37(1)),(wz(1,38),ws38(1)),
     3      (wz(1,39),ws39(1)),(wz(1,40),ws40(1)),(wz(1,41),ws41(1))
  
c 
        data xs1/
     1  -.9929349694864852E+00,-.9546335944914043E+00,
     2  -.8717852572260915E+00,-.7414887036898580E+00,
     3  -.5678291872175965E+00,-.3598465064367627E+00,
     4  -.1300030426973473E+00,0.1072530794531561E+00,
     5  0.3369164800396298E+00,0.5448524364662104E+00,
     6  0.7191778290541933E+00,0.8516039621127764E+00,
     7  0.9388086195832428E+00,0.9840021106893595E+00,
     8  0.9985144799124259E+00,0.9999874974352061E+00,
     9  0.9999970736640207E+00,0.9999976982761511E+00,
     *  0.9999978944165662E+00,0.9999980011525864E+00,
     *  0.9999980816660963E+00,0.9999981529547282E+00,
     *  0.9999982198790545E+00,0.9999982838741142E+00,
     *  0.9999983454171930E+00,0.9999984047141484E+00,
     *  0.9999984618894662E+00,0.9999985170416734E+00,
     *  0.9999985702575656E+00,0.9999986216158261E+00,
     *  0.9999986711897382E+00,0.9999987190514874E+00,
     *  0.9999987652733751E+00,0.9999988099304704E+00,
     *  0.9999988531177448E+00,0.9999988949998388E+00,
     *  0.9999989359757959E+00,0.9999989772402327E+00,
     *  0.9999990226766755E+00,0.9999990846314435E+00,
     *  0.9999992031797764E+00,0.9999995319421162E+00,
     a   0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data ws1/
     1  0.1948891020747820E-01,0.5935798536797277E-01,
     2  0.1067557823541683E+00,0.1531261674930292E+00,
     3  0.1926612754996361E+00,0.2212030087289195E+00,
     4  0.2360522847872286E+00,0.2359349248060096E+00,
     5  0.2210130564209895E+00,0.1928560305103248E+00,
     6  0.1543908850179433E+00,0.1098820038509763E+00,
     7  0.6503118307087710E-01,0.2724326810452467E-01,
     8  0.4950408854380404E-02,0.4904351394787876E-04,
     9  0.1284125323338934E-05,0.2967767900181033E-06,
     *  0.1331723262828402E-06,0.8889677194670255E-07,
     *  0.7455408332605487E-07,0.6871674154354716E-07,
     *  0.6533813089706866E-07,0.6272252827280714E-07,
     *  0.6039541787333059E-07,0.5821846975957651E-07,
     *  0.5614845932174779E-07,0.5417032452662999E-07,
     *  0.5227456975967282E-07,0.5045405567174443E-07,
     *  0.4870579515214341E-07,0.4702975763691292E-07,
     *  0.4542622686269433E-07,0.4390315839801963E-07,
     *  0.4249657019779404E-07,0.4132665082879076E-07,
     *  0.4079819347177795E-07,0.4229499683861630E-07,
     *  0.5039682930736324E-07,0.7917825499569254E-07,
     *  0.1787556389073350E-06,0.5590109896915794E-06,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data xs2/
     1  -.9948308629309646E+00,-.9640578203611618E+00,
     2  -.8931454117897650E+00,-.7770806656039737E+00,
     3  -.6178973583188261E+00,-.4227656852492337E+00,
     4  -.2025649145622004E+00,0.2943779598608720E-01,
     5  0.2589371632957183E+00,0.4719666301835863E+00,
     6  0.6562661423099163E+00,0.8026018249432969E+00,
     7  0.9060349269236057E+00,0.9672031604220552E+00,
     8  0.9936091646093450E+00,0.9996247942850531E+00,
     9  0.9999884663660503E+00,0.9999958371310651E+00,
     *  0.9999966939190366E+00,0.9999969796435317E+00,
     *  0.9999971322212762E+00,0.9999972424166643E+00,
     *  0.9999973370713082E+00,0.9999974247724599E+00,
     *  0.9999975083199291E+00,0.9999975887264407E+00,
     *  0.9999976664049783E+00,0.9999977415369630E+00,
     *  0.9999978142380110E+00,0.9999978846251683E+00,
     *  0.9999979527976808E+00,0.9999980188249815E+00,
     *  0.9999980827728689E+00,0.9999981447188169E+00,
     *  0.9999982047448253E+00,0.9999982629225823E+00,
     *  0.9999983193172548E+00,0.9999983740232550E+00,
     *  0.9999984272036400E+00,0.9999984792110175E+00,
     *  0.9999985309741208E+00,0.9999985851326314E+00,
     *  0.9999986492355142E+00,0.9999987443075713E+00,
     *  0.9999989335925765E+00,0.9999994305660074E+00,
     a   0.0d0/
c 
        data ws2/
     1  0.1472129550804123E-01,0.4929528778827950E-01,
     2  0.9331656787742840E-01,0.1384431179738539E+00,
     3  0.1787055535594314E+00,0.2097226773236341E+00,
     4  0.2284451108522742E+00,0.2331532963579966E+00,
     5  0.2234949222734104E+00,0.2004986049856327E+00,
     6  0.1665403915725978E+00,0.1252867396298554E+00,
     7  0.8166278459273880E-01,0.4191118958601264E-01,
     8  0.1338819376933992E-01,0.1383797982651257E-02,
     9  0.2525392181001092E-04,0.1666108540084189E-05,
     *  0.4290903609086513E-06,0.1935221627712151E-06,
     *  0.1241927349299030E-06,0.1001254187384609E-06,
     *  0.9044168283068640E-07,0.8537132414073449E-07,
     *  0.8187672302130234E-07,0.7900017571544875E-07,
     *  0.7638448268872722E-07,0.7389753069624654E-07,
     *  0.7152479132255088E-07,0.6926629052895061E-07,
     *  0.6708968641730531E-07,0.6497564006494854E-07,
     *  0.6293318009790270E-07,0.6097270353436628E-07,
     *  0.5909130118959799E-07,0.5727424518077745E-07,
     *  0.5553045579762284E-07,0.5390705942434073E-07,
     *  0.5250693614476564E-07,0.5164664674288179E-07,
     *  0.5227223635490646E-07,0.5718156540430345E-07,
     *  0.7416961030749466E-07,0.1250419822163543E-06,
     *  0.2845308165578807E-06,0.7895648325186058E-06,
     a  0.0d0/
  
c 
        data xs3/
     1  -.9948620670880578E+00,-.9642274658373227E+00,
     2  -.8935470942640659E+00,-.7777639324670259E+00,
     3  -.6188685035252204E+00,-.4239939600632840E+00,
     4  -.2039904211631161E+00,0.2789594772931214E-01,
     5  0.2573726635704534E+00,0.4704778163693797E+00,
     6  0.6549471605213413E+00,0.8015323797383098E+00,
     7  0.9052678296221153E+00,0.9667490303799366E+00,
     8  0.9934192749048574E+00,0.9995868701163781E+00,
     9  0.9999845940952687E+00,0.9999941527246623E+00,
     *  0.9999953333795091E+00,0.9999957323065329E+00,
     *  0.9999959466919108E+00,0.9999961021140475E+00,
     *  0.9999962358744809E+00,0.9999963598928951E+00,
     *  0.9999964780319222E+00,0.9999965916904410E+00,
     *  0.9999967014574809E+00,0.9999968076053758E+00,
     *  0.9999969103109829E+00,0.9999970097447988E+00,
     *  0.9999971060489663E+00,0.9999971993320884E+00,
     *  0.9999972896887430E+00,0.9999973772175741E+00,
     *  0.9999974620262224E+00,0.9999975442132962E+00,
     *  0.9999976238807908E+00,0.9999977011748744E+00,
     *  0.9999977763339496E+00,0.9999978498653430E+00,
     *  0.9999979230998468E+00,0.9999979998079607E+00,
     *  0.9999980907744948E+00,0.9999982260554743E+00,
     *  0.9999984960481552E+00,0.9999992034122889E+00,
     a  0.0d0/
c 
        data ws3/
     1  0.1463965752761892E-01,0.4910332451517414E-01,
     2  0.9305217769278948E-01,0.1381515906886437E+00,
     3  0.1784275151113734E+00,0.2094913220273606E+00,
     4  0.2282855272226948E+00,0.2330823619683943E+00,
     5  0.2235213737195730E+00,0.2006228574036873E+00,
     6  0.1667534318645852E+00,0.1255683317163964E+00,
     7  0.8197890757977169E-01,0.4221098591675810E-01,
     8  0.1360525082294536E-01,0.1466625802235010E-02,
     9  0.3154022511317785E-04,0.2275879442343147E-05,
     *  0.5968642851388984E-06,0.2712780309866931E-06,
     *  0.1749100004705489E-06,0.1413952620066618E-06,
     *  0.1278709122505578E-06,0.1207312001221355E-06,
     *  0.1157582120324875E-06,0.1116505154139568E-06,
     *  0.1079253845937768E-06,0.1043982734427909E-06,
     *  0.1010423123428573E-06,0.9784894494787608E-07,
     *  0.9477712146101555E-07,0.9180458244128715E-07,
     *  0.8892497843984050E-07,0.8615114762456150E-07,
     *  0.8348236742221631E-07,0.8090790555194595E-07,
     *  0.7845101546993312E-07,0.7617445889165677E-07,
     *  0.7422020928935701E-07,0.7304177936865034E-07,
     *  0.7398795728315347E-07,0.8105214773365589E-07,
     *  0.1053767666550596E-06,0.1781528429116187E-06,
     *  0.4060412312312386E-06,0.1118310075022598E-05,
     a  0.0d0/
c 
        data xs4/
     1  -.9948942257696641E+00,-.9643998840110680E+00,
     2  -.8939543657234619E+00,-.7784573275728482E+00,
     3  -.6198552507733976E+00,-.4252430820092907E+00,
     4  -.2054409662436956E+00,0.2632638033575652E-01,
     5  0.2557793684953877E+00,0.4689603463929206E+00,
     6  0.6536001398611095E+00,0.8004351624349870E+00,
     7  0.9044726177129392E+00,0.9662676980829425E+00,
     8  0.9932093556059851E+00,0.9995419724774943E+00,
     9  0.9999794118336425E+00,0.9999917945803458E+00,
     *  0.9999934160618884E+00,0.9999939715208741E+00,
     *  0.9999942721773596E+00,0.9999944910742715E+00,
     *  0.9999946798298895E+00,0.9999948549840418E+00,
     *  0.9999950219222775E+00,0.9999951825658813E+00,
     *  0.9999953377052340E+00,0.9999954877033830E+00,
     *  0.9999956328113355E+00,0.9999957732890947E+00,
     *  0.9999959093517647E+00,0.9999960411503239E+00,
     *  0.9999961688086794E+00,0.9999962924607513E+00,
     *  0.9999964122640872E+00,0.9999965283677028E+00,
     *  0.9999966409247598E+00,0.9999967501458242E+00,
     *  0.9999968563661439E+00,0.9999969603016750E+00,
     *  0.9999970638458611E+00,0.9999971723989891E+00,
     *  0.9999973013793770E+00,0.9999974937215435E+00,
     *  0.9999978785827668E+00,0.9999988852902921E+00,
     a  0.0d0/
c 
        data ws4/
     1  0.1455613800565640E-01,0.4890935782127279E-01,
     2  0.9278408670736497E-01,0.1378548827255364E+00,
     3  0.1781439429507302E+00,0.2092551476297313E+00,
     4  0.2281224294746976E+00,0.2330095423596237E+00,
     5  0.2235474461411787E+00,0.2007476314965190E+00,
     6  0.1669670413657897E+00,0.1258498252036376E+00,
     7  0.8229472077514921E-01,0.4251385546148265E-01,
     8  0.1383392881404580E-01,0.1560575266863172E-02,
     9  0.3941970917232509E-04,0.3097817760086794E-05,
     *  0.8276961267050803E-06,0.3794124343589124E-06,
     *  0.2459544412561049E-06,0.1993988849593675E-06,
     *  0.1805355779846744E-06,0.1705632284175428E-06,
     *  0.1636017185788695E-06,0.1578083721873947E-06,
     *  0.1525254753723679E-06,0.1475105663147389E-06,
     *  0.1427504343792009E-06,0.1382410905297209E-06,
     *  0.1339086735035987E-06,0.1297084698393308E-06,
     *  0.1256301119063530E-06,0.1217012486734430E-06,
     *  0.1179299868266278E-06,0.1143016098453894E-06,
     *  0.1108470922722882E-06,0.1076472085681955E-06,
     *  0.1049014872923203E-06,0.1032520836855205E-06,
     *  0.1046403685876232E-06,0.1147879978317819E-06,
     *  0.1495935296580761E-06,0.2536238443006931E-06,
     *  0.5791623634255399E-06,0.1584108617324561E-05,
     a  0.0d0/
c 
        data xs5/
     1  -.9949114928470694E+00,-.9644956273997801E+00,
     2  -.8941901223161488E+00,-.7788766164364096E+00,
     3  -.6204789086147513E+00,-.4260687467266226E+00,
     4  -.2064443896797561E+00,0.2518893979116776E-01,
     5  0.2545678426947973E+00,0.4677466266848149E+00,
     6  0.6524626595763410E+00,0.7994517829541863E+00,
     7  0.9037110836628279E+00,0.9657717926763876E+00,
     8  0.9929758769286541E+00,0.9994878650468919E+00,
     9  0.9999724118922244E+00,0.9999884858351262E+00,
     *  0.9999907115747073E+00,0.9999914847871980E+00,
     *  0.9999919062035059E+00,0.9999922140632918E+00,
     *  0.9999924798791420E+00,0.9999927267772997E+00,
     *  0.9999929623410238E+00,0.9999931892011746E+00,
     *  0.9999934083968274E+00,0.9999936204440683E+00,
     *  0.9999938256844747E+00,0.9999940243969163E+00,
     *  0.9999942168270621E+00,0.9999944031976473E+00,
     *  0.9999945837197881E+00,0.9999947586034179E+00,
     *  0.9999949280544049E+00,0.9999950922781000E+00,
     *  0.9999952515009829E+00,0.9999954060107623E+00,
     *  0.9999955562758143E+00,0.9999957033113678E+00,
     *  0.9999958497928385E+00,0.9999960034030084E+00,
     *  0.9999961861838642E+00,0.9999964595275419E+00,
     *  0.9999970080974896E+00,0.9999984409287574E+00,
     a  0.0d0/
c 
        data ws5/
     1  0.1451085698452428E-01,0.4879791500834598E-01,
     2  0.9261872574455886E-01,0.1376570269858090E+00,
     3  0.1779369592482559E+00,0.2090617826933348E+00,
     4  0.2279635265269891E+00,0.2329030852534003E+00,
     5  0.2235077209635528E+00,0.2007840886561713E+00,
     6  0.1670830100041986E+00,0.1260404424066818E+00,
     7  0.8254394631141011E-01,0.4278831479695321E-01,
     8  0.1407048020759697E-01,0.1668685207360805E-02,
     9  0.4945322343160273E-04,0.4213953808174894E-05,
     *  0.1147413020557955E-05,0.5304961514723998E-06,
     *  0.3455267882978840E-06,0.2806740226187982E-06,
     *  0.2543553834578468E-06,0.2405573704263803E-06,
     *  0.2309622140700937E-06,0.2229164094518161E-06,
     *  0.2155560422040854E-06,0.2085937896540875E-06,
     *  0.2019331501820264E-06,0.1955325436889254E-06,
     *  0.1893646947412183E-06,0.1834113989862599E-06,
     *  0.1776680213776754E-06,0.1721335140867419E-06,
     *  0.1668022657255077E-06,0.1616823752541377E-06,
     *  0.1568096595733633E-06,0.1522828378723662E-06,
     *  0.1484008733444485E-06,0.1460687332372103E-06,
     *  0.1480361461674811E-06,0.1625049631471485E-06,
     *  0.2122369996390567E-06,0.3609565558090389E-06,
     *  0.8261846959454907E-06,0.2243536618607448E-05,
     a  0.0d0/
c 
        data xs6/
     1  -.9949456278622959E+00,-.9646808683374459E+00,
     2  -.8946319613581606E+00,-.7796361800896779E+00,
     3  -.6215712634658639E+00,-.4274677764906883E+00,
     4  -.2080899144180757E+00,0.2338345476058223E-01,
     5  0.2527071181473086E+00,0.4659448345805328E+00,
     6  0.6508334933158682E+00,0.7980966837506429E+00,
     7  0.9027050055901888E+00,0.9651455864834935E+00,
     8  0.9926931597173301E+00,0.9994225072816068E+00,
     9  0.9999629806066801E+00,0.9999838470451208E+00,
     *  0.9999868977874251E+00,0.9999879729227892E+00,
     *  0.9999885634524761E+00,0.9999889968995426E+00,
     *  0.9999893720099702E+00,0.9999897207316111E+00,
     *  0.9999900535309723E+00,0.9999903740420020E+00,
     *  0.9999906837173231E+00,0.9999909832932696E+00,
     *  0.9999912732600569E+00,0.9999915540141204E+00,
     *  0.9999918259003169E+00,0.9999920892305303E+00,
     *  0.9999923443038450E+00,0.9999925914131345E+00,
     *  0.9999928308460866E+00,0.9999930628945741E+00,
     *  0.9999932878811222E+00,0.9999935062153134E+00,
     *  0.9999937185654100E+00,0.9999939263823017E+00,
     *  0.9999941334970261E+00,0.9999943509197295E+00,
     *  0.9999946102086741E+00,0.9999949992015567E+00,
     *  0.9999957820702263E+00,0.9999978222204848E+00,
     a  0.0d0/
c 
        data ws6/
     1  0.1442178284584919E-01,0.4858757958254147E-01,
     2  0.9232368889265644E-01,0.1373243796850414E+00,
     3  0.1776109274201818E+00,0.2087801531418961E+00,
     4  0.2277565685415041E+00,0.2327931508638264E+00,
     5  0.2235087708259406E+00,0.2009008567368143E+00,
     6  0.1673094946470465E+00,0.1263577297339081E+00,
     7  0.8291744529418828E-01,0.4316344721261945E-01,
     8  0.1436651554146949E-01,0.1795946352463310E-02,
     9  0.6210245902401683E-04,0.5723219609429169E-05,
     *  0.1588531855180321E-05,0.7411836619163714E-06,
     *  0.4856109004184605E-06,0.3957628946953183E-06,
     *  0.3591519222023706E-06,0.3398316116549019E-06,
     *  0.3263085887527572E-06,0.3149353389670983E-06,
     *  0.3045314880062984E-06,0.2946999397417463E-06,
     *  0.2852992063531304E-06,0.2762661440110816E-06,
     *  0.2675578107114644E-06,0.2591524031617738E-06,
     *  0.2510430413578515E-06,0.2432233729158865E-06,
     *  0.2356906959656997E-06,0.2284593216730894E-06,
     *  0.2215796095785376E-06,0.2151928331961363E-06,
     *  0.2097268051548021E-06,0.2064774172173250E-06,
     *  0.2093889372119774E-06,0.2302123161282602E-06,
     *  0.3014981197564166E-06,0.5144228728759304E-06,
     *  0.1179775488182915E-05,0.3176353998551425E-05,
     a  0.0d0/
c 
        data xs7/
     1  -.9950240921807282E+00,-.9651016767321898E+00,
     2  -.8956225838211471E+00,-.7813163159636178E+00,
     3  -.6239546053183669E+00,-.4304777156086835E+00,
     4  -.2115790442160642E+00,0.1961365633030221E-01,
     5  0.2488864802872095E+00,0.4623139891761363E+00,
     6  0.6476211305581343E+00,0.7954933179207520E+00,
     7  0.9008328191910178E+00,0.9640248526649480E+00,
     8  0.9922080594756399E+00,0.9993118582248073E+00,
     9  0.9999471154148403E+00,0.9999766740980709E+00,
     *  0.9999812264206824E+00,0.9999828330357379E+00,
     *  0.9999836981137125E+00,0.9999843167621039E+00,
     *  0.9999848428572036E+00,0.9999853282390434E+00,
     *  0.9999857903888685E+00,0.9999862353921289E+00,
     *  0.9999866656085865E+00,0.9999870821651746E+00,
     *  0.9999874857562125E+00,0.9999878769106096E+00,
     *  0.9999882560828177E+00,0.9999886236906950E+00,
     *  0.9999889801322976E+00,0.9999893257798744E+00,
     *  0.9999896609843157E+00,0.9999899861000243E+00,
     *  0.9999903014952054E+00,0.9999906076018423E+00,
     *  0.9999909050238849E+00,0.9999911948384168E+00,
     *  0.9999914794903548E+00,0.9999917653627126E+00,
     *  0.9999920703165810E+00,0.9999924435708791E+00,
     *  0.9999930172839340E+00,0.9999941767891507E+00,
     *  0.9999970982040363E+00/
c 
        data ws7/
     1  0.1421780801460081E-01,0.4811513033844099E-01,
     2  0.9167545803562738E-01,0.1366114309663752E+00,
     3  0.1769325769808265E+00,0.2082171137608964E+00,
     4  0.2273697572269029E+00,0.2326238390182035E+00,
     5  0.2235782442520796E+00,0.2012090388850110E+00,
     6  0.1678319847166894E+00,0.1270419016180320E+00,
     7  0.8367868890608676E-01,0.4388218722121314E-01,
     8  0.1489284760466838E-01,0.2008894246821693E-02,
     9  0.8482272767243177E-04,0.8478291093831059E-05,
     *  0.2382363562876125E-05,0.1098249593641101E-05,
     *  0.7013089042067874E-06,0.5587316112123132E-06,
     *  0.5011488033388054E-06,0.4722024293597779E-06,
     *  0.4530039052870719E-06,0.4373611962935806E-06,
     *  0.4232485002034267E-06,0.4099765601593008E-06,
     *  0.3972924488670385E-06,0.3850914880863359E-06,
     *  0.3733222540088035E-06,0.3619600648183291E-06,
     *  0.3509851530573553E-06,0.3403678594072586E-06,
     *  0.3301002110538532E-06,0.3201916259647683E-06,
     *  0.3106692210976703E-06,0.3016406814761008E-06,
     *  0.2933696886691431E-06,0.2866313243069364E-06,
     *  0.2836272517643653E-06,0.2907948647827421E-06,
     *  0.3265953702225743E-06,0.4397535156091245E-06,
     *  0.7635564566258294E-06,0.1736130191861373E-05,
     *  0.4409975429807926E-05/
  
c 
        data xs8/
     1  -.9950659424626611E+00,-.9653305250201182E+00,
     2  -.8961709758517893E+00,-.7822626203689349E+00,
     3  -.6253208289744983E+00,-.4322352157006627E+00,
     4  -.2136564971305051E+00,0.1732167559257136E-01,
     5  0.2465100793460315E+00,0.4599973657377241E+00,
     6  0.6455100708714492E+00,0.7937202842990123E+00,
     7  0.8994991704028761E+00,0.9631787714749178E+00,
     8  0.9918133944947627E+00,0.9992116836663942E+00,
     9  0.9999290486633874E+00,0.9999673362398629E+00,
     *  0.9999735429495047E+00,0.9999757664390565E+00,
     *  0.9999769744009031E+00,0.9999778435486415E+00,
     *  0.9999785850472176E+00,0.9999792700229500E+00,
     *  0.9999799224878376E+00,0.9999805508435800E+00,
     *  0.9999811583620069E+00,0.9999817466025205E+00,
     *  0.9999823165336708E+00,0.9999828688993935E+00,
     *  0.9999834043335607E+00,0.9999839234300081E+00,
     *  0.9999844267666500E+00,0.9999849148800268E+00,
     *  0.9999853882622158E+00,0.9999858474022156E+00,
     *  0.9999862928143046E+00,0.9999867251106807E+00,
     *  0.9999871451546232E+00,0.9999875544905910E+00,
     *  0.9999879566193497E+00,0.9999883606858735E+00,
     *  0.9999887922780086E+00,0.9999893217906461E+00,
     *  0.9999901380781289E+00,0.9999917917006358E+00,
     *  0.9999959471495492E+00/
c 
        data ws8/
     1  0.1410823749585400E-01,0.4785396459811838E-01,
     2  0.9130709554531141E-01,0.1361934675792144E+00,
     3  0.1765191802272488E+00,0.2078550826583606E+00,
     4  0.2270975825156295E+00,0.2324710760231148E+00,
     5  0.2235645457145783E+00,0.2013426001652986E+00,
     6  0.1681073797300681E+00,0.1274374375477875E+00,
     7  0.8415320676889758E-01,0.4436898572686881E-01,
     8  0.1528774826613207E-01,0.2186599778586870E-02,
     9  0.1064993187953012E-03,0.1145185617329062E-04,
     *  0.3281632458443737E-05,0.1528094111490013E-05,
     *  0.9829116956522648E-06,0.7866057458890446E-06,
     *  0.7069397611645627E-06,0.6665668858215487E-06,
     *  0.6396164985771512E-06,0.6175926737531941E-06,
     *  0.5976883581315997E-06,0.5789488687030850E-06,
     *  0.5610361416137761E-06,0.5437995611434986E-06,
     *  0.5271673602735776E-06,0.5111225076707161E-06,
     *  0.4956403308556541E-06,0.4806674471107442E-06,
     *  0.4661782785828919E-06,0.4521861240828533E-06,
     *  0.4387379574767887E-06,0.4259928787318629E-06,
     *  0.4143336528048276E-06,0.4048707192773638E-06,
     *  0.4007543712935318E-06,0.4112159828392915E-06,
     *  0.4626734253445627E-06,0.6246983843703721E-06,
     *  0.1087756685176471E-05,0.2477042916364211E-05,
     *  0.6235931671460165E-05/
  
c 
        data xs9/
     1  -.9937592989546193E+00,-.9584957156004691E+00,
     2  -.8801808849524119E+00,-.7550887873388992E+00,
     3  -.5865980152829318E+00,-.3831303951458523E+00,
     4  -.1566142892338858E+00,0.7891145422334733E-01,
     5  0.3086854142212521E+00,0.5186115727750924E+00,
     6  0.6966263363020264E+00,0.8340483873690692E+00,
     7  0.9269556881327394E+00,0.9777077553040472E+00,
     8  0.9964591046335487E+00,0.9997279859170343E+00,
     9  0.9999413900682463E+00,0.9999603177253333E+00,
     *  0.9999650442612407E+00,0.9999671759983559E+00,
     *  0.9999685807526117E+00,0.9999697427165413E+00,
     *  0.9999708062504942E+00,0.9999718152166615E+00,
     *  0.9999727837527532E+00,0.9999737171893297E+00,
     *  0.9999746181387632E+00,0.9999754882907470E+00,
     *  0.9999763289734551E+00,0.9999771413685830E+00,
     *  0.9999779265795510E+00,0.9999786856391952E+00,
     *  0.9999794195345346E+00,0.9999801292520627E+00,
     *  0.9999808158408907E+00,0.9999814805545579E+00,
     *  0.9999821253421308E+00,0.9999827543921700E+00,
     *  0.9999833791803696E+00,0.9999840347600126E+00,
     *  0.9999848282802990E+00,0.9999860761980850E+00,
     *  0.9999887716064833E+00,0.9999954774606028E+00,
     a  0.0d0,0.0d0,0.0d0/
c 
  
        data ws9/
     1  0.1746321581492193E-01,0.5543245224257121E-01,
     2  0.1017581763115299E+00,0.1478307262956063E+00,
     3  0.1877223082227413E+00,0.2172016158274593E+00,
     4  0.2334673745490377E+00,0.2351004457408108E+00,
     5  0.2220752764498276E+00,0.1957441212366385E+00,
     6  0.1588159625398695E+00,0.1153494621611288E+00,
     7  0.7082928838728603E-01,0.3236471844727352E-01,
     8  0.8093010109667785E-02,0.6647429993275806E-03,
     9  0.4231730277644703E-04,0.7793866737035458E-05,
     *  0.2874120786145756E-05,0.1634832889431528E-05,
     *  0.1243934951822335E-05,0.1100863523875995E-05,
     *  0.1032549104224789E-05,0.9874579096400172E-06,
     *  0.9504208625997126E-06,0.9168634036556762E-06,
     *  0.8853094530559405E-06,0.8552124671317667E-06,
     *  0.8263489731078218E-06,0.7986253133004467E-06,
     *  0.7719682504109985E-06,0.7463148665106903E-06,
     *  0.7216391289770802E-06,0.6979692305989549E-06,
     *  0.6754102337256468E-06,0.6543210048427801E-06,
     *  0.6358830771555545E-06,0.6239037882757948E-06,
     *  0.6308077424624275E-06,0.6962990380483741E-06,
     *  0.9375034790736982E-06,0.1701128224051230E-05,
     *  0.4158859292363440E-05,0.8874875902421340E-05,
     a  0.0d0,0.0d0,0.0d0/
c 
        data xs10/
     1  -.9951255682402216E+00,-.9656678290032594E+00,
     2  -.8970042711256795E+00,-.7837417581093840E+00,
     3  -.6275150240204389E+00,-.4351338847821891E+00,
     4  -.2171746392506158E+00,0.1333552069643917E-01,
     5  0.2422627244360023E+00,0.4557362946535617E+00,
     6  0.6415044063607771E+00,0.7902368553889421E+00,
     7  0.8967722254985340E+00,0.9613670097498397E+00,
     8  0.9909226288125701E+00,0.9989680608248899E+00,
     9  0.9998766382443458E+00,0.9999373410891073E+00,
     *  0.9999481238483221E+00,0.9999521452971212E+00,
     *  0.9999544246361820E+00,0.9999561301572346E+00,
     *  0.9999576196786213E+00,0.9999590090705852E+00,
     *  0.9999603364052503E+00,0.9999616150424084E+00,
     *  0.9999628504378704E+00,0.9999640454489445E+00,
     *  0.9999652020150711E+00,0.9999663216943642E+00,
     *  0.9999674058648236E+00,0.9999684558126458E+00,
     *  0.9999694727443560E+00,0.9999704578151651E+00,
     *  0.9999714121461072E+00,0.9999723368623523E+00,
     *  0.9999732332505396E+00,0.9999741030119265E+00,
     *  0.9999749489504566E+00,0.9999757771818907E+00,
     *  0.9999766041333553E+00,0.9999774772534984E+00,
     *  0.9999785312835843E+00,0.9999801376567161E+00,
     *  0.9999834058620503E+00,0.9999917890577492E+00,
     a  0.0d0/
c 
        data ws10/
     1  0.1395018386399575E-01,0.4745817413266162E-01,
     2  0.9072308261948257E-01,0.1354994192410846E+00,
     3  0.1757964875429762E+00,0.2071812753993206E+00,
     4  0.2265435693891833E+00,0.2320983139828617E+00,
     5  0.2234217481510049E+00,0.2014616600616215E+00,
     6  0.1684988575701339E+00,0.1280851793683064E+00,
     7  0.8500635936828491E-01,0.4532191737365363E-01,
     8  0.1612866740822501E-01,0.2591226342468565E-02,
     9  0.1601441499759759E-03,0.1951847668250808E-04,
     *  0.5838343531383873E-05,0.2826370324455542E-05,
     *  0.1896944700940157E-05,0.1566636550139228E-05,
     *  0.1429424321682945E-05,0.1354966544799175E-05,
     *  0.1301697676915743E-05,0.1256412381222969E-05,
     *  0.1214831571971866E-05,0.1175506432043135E-05,
     *  0.1137882649579798E-05,0.1101703950817304E-05,
     *  0.1066851642081826E-05,0.1033244306927783E-05,
     *  0.1000811476230592E-05,0.9695167811451501E-06,
     *  0.9393294972430140E-06,0.9103142301025251E-06,
     *  0.8827360928958374E-06,0.8572239594168746E-06,
     *  0.8355736741074822E-06,0.8233305207478224E-06,
     *  0.8376422556610332E-06,0.9289619115090326E-06,
     *  0.1234649623292909E-05,0.2137966822260639E-05,
     *  0.4929735780211138E-05,0.1268347070895575E-04,
     a  0.0d0/
c 
        data xs11/
     1  -.9952077697360574E+00,-.9661187007515741E+00,
     2  -.8980848992586246E+00,-.7856053664580702E+00,
     3  -.6302040089055656E+00,-.4385916354688519E+00,
     4  -.2212607692533157E+00,0.8828925202720495E-02,
     5  0.2375927801318287E+00,0.4511879514940729E+00,
     6  0.6373648185361263E+00,0.7867651148066388E+00,
     7  0.8941642568455510E+00,0.9597128692982324E+00,
     8  0.9901460655303596E+00,0.9987586375057309E+00,
     9  0.9998307546087275E+00,0.9999117557751858E+00,
     *  0.9999267295302545E+00,0.9999323717751847E+00,
     *  0.9999355848630056E+00,0.9999379959696892E+00,
     *  0.9999401048328669E+00,0.9999420730854987E+00,
     *  0.9999439537415267E+00,0.9999457653784274E+00,
     *  0.9999475155963629E+00,0.9999492084068655E+00,
     *  0.9999508465656288E+00,0.9999524322812332E+00,
     *  0.9999539674810753E+00,0.9999554539711173E+00,
     *  0.9999568934640914E+00,0.9999582875887985E+00,
     *  0.9999596379167193E+00,0.9999609460046675E+00,
     *  0.9999622136277158E+00,0.9999634431873567E+00,
     *  0.9999646387076411E+00,0.9999658091121806E+00,
     *  0.9999669785950427E+00,0.9999682171312620E+00,
     *  0.9999697218807353E+00,0.9999720301929242E+00,
     *  0.9999767300328606E+00,0.9999886481454527E+00,
     a  0.0d0/
c 
        data ws11/
     1  0.1373459456155187E-01,0.4694309646767363E-01,
     2  0.8999766478895242E-01,0.1346775396168795E+00,
     3  0.1749841626576939E+00,0.2064700497721745E+00,
     4  0.2260093302820024E+00,0.2317998490983233E+00,
     5  0.2233982233290673E+00,0.2017288371194961E+00,
     6  0.1690445950334200E+00,0.1288647148872162E+00,
     7  0.8593701091202434E-01,0.4627168443622228E-01,
     8  0.1689273174091132E-01,0.2934249314105930E-02,
     9  0.2081223043952888E-03,0.2690849172420329E-04,
     *  0.8167909604421820E-05,0.3977203232033945E-05,
     *  0.2678723048300635E-05,0.2216876554157773E-05,
     *  0.2024579675764907E-05,0.1919730038793634E-05,
     *  0.1844354809673390E-05,0.1780083449894511E-05,
     *  0.1720988609506862E-05,0.1665082301222302E-05,
     *  0.1611599652818946E-05,0.1560149630258688E-05,
     *  0.1510551156828290E-05,0.1462714245013457E-05,
     *  0.1416542512021025E-05,0.1371968898740031E-05,
     *  0.1328942019990744E-05,0.1287523786840430E-05,
     *  0.1248109380681801E-05,0.1211635584286150E-05,
     *  0.1180746166139377E-05,0.1163656145137173E-05,
     *  0.1185750770163329E-05,0.1321203381498364E-05,
     *  0.1768799313120938E-05,0.3077331003760257E-05,
     *  0.7074114365814196E-05,0.1782474848733653E-04,
     a  0.0d0/
  
c 
        data xs12/
     1  -.9932764815535604E+00,-.9562871421531310E+00,
     2  -.8755064953720098E+00,-.7477250128869437E+00,
     3  -.5767284963449109E+00,-.3712683983686061E+00,
     4  -.1435321517898747E+00,0.9227359444042760E-01,
     5  0.3213297365858179E+00,0.5296099190312176E+00,
     6  0.7052495164319776E+00,0.8398889179969391E+00,
     7  0.9300500354092258E+00,0.9786605853744712E+00,
     8  0.9963981883642075E+00,0.9995813987944435E+00,
     9  0.9998549255107585E+00,0.9998920680026611E+00,
     *  0.9999029951536615E+00,0.9999083962255153E+00,
     *  0.9999122047849758E+00,0.9999154725730739E+00,
     *  0.9999185051179187E+00,0.9999213936381750E+00,
     *  0.9999241682133160E+00,0.9999268411193110E+00,
     *  0.9999294190913168E+00,0.9999319069277883E+00,
     *  0.9999343086178085E+00,0.9999366277195212E+00,
     *  0.9999388675336518E+00,0.9999410312414423E+00,
     *  0.9999431219890702E+00,0.9999451429959059E+00,
     *  0.9999470979881140E+00,0.9999489923844807E+00,
     *  0.9999508369726207E+00,0.9999526603696696E+00,
     *  0.9999545509947529E+00,0.9999567904997804E+00,
     *  0.9999602445601452E+00,0.9999676608399710E+00,
     *  0.9999861649272348E+00,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data ws12/
     1  0.1864021078115810E-01,0.5762003452050739E-01,
     2  0.1044125897270160E+00,0.1504885236588368E+00,
     3  0.1900209616300190E+00,0.2188442827013445E+00,
     4  0.2342371377684144E+00,0.2348810568703692E+00,
     5  0.2208708012876840E+00,0.1936898154899376E+00,
     6  0.1561734650671833E+00,0.1125002870888298E+00,
     7  0.6828206139750334E-01,0.3072940041190430E-01,
     8  0.7675097928502930E-02,0.7414060744501512E-03,
     9  0.7511566584102510E-04,0.1721532472213776E-04,
     *  0.7019241706651689E-05,0.4313546035395801E-05,
     *  0.3453908279725035E-05,0.3125447316017041E-05,
     *  0.2952793915386297E-05,0.2828684334536528E-05,
     *  0.2722348576670642E-05,0.2624536952745005E-05,
     *  0.2532191018165156E-05,0.2444141374918596E-05,
     *  0.2359829330048425E-05,0.2278921639502726E-05,
     *  0.2201235775404022E-05,0.2126703240792345E-05,
     *  0.2055322068188879E-05,0.1987295881723094E-05,
     *  0.1923542719224519E-05,0.1866878364221111E-05,
     *  0.1826498175298796E-05,0.1833132004715111E-05,
     *  0.1989543906651232E-05,0.2616992129258349E-05,
     *  0.4687217467380602E-05,0.1142977940222160E-04,
     *  0.2513246157420696E-04,0.0d0,0.0d0,0.0d0,0.0d0/
  
  
  
  
  
        data xs13/
     1  -.9942892143104046E+00,-.9612220582454930E+00,
     2  -.8865497325024589E+00,-.7659883061334046E+00,
     3  -.6023113654919821E+00,-.4033413338262190E+00,
     4  -.1804662973385407E+00,0.5270604689584490E-01,
     5  0.2817043400834260E+00,0.4925674309130303E+00,
     6  0.6732160246284491E+00,0.8148022545879732E+00,
     7  0.9130480607708280E+00,0.9696264453714039E+00,
     8  0.9933723781033912E+00,0.9990965374217491E+00,
     9  0.9997535410765566E+00,0.9998353891922502E+00,
     *  0.9998572116434352E+00,0.9998667165056656E+00,
     *  0.9998726127882698E+00,0.9998772735760104E+00,
     *  0.9998814494295667E+00,0.9998853810765422E+00,
     *  0.9998891476776044E+00,0.9998927780264337E+00,
     *  0.9998962849348905E+00,0.9998996759447112E+00,
     *  0.9999029565955757E+00,0.9999061314209228E+00,
     *  0.9999092044683247E+00,0.9999121795500618E+00,
     *  0.9999150602944382E+00,0.9999178502402282E+00,
     *  0.9999205528686721E+00,0.9999231718648265E+00,
     *  0.9999257118043512E+00,0.9999281795080572E+00,
     *  0.9999305882871481E+00,0.9999329716466705E+00,
     *  0.9999354262475713E+00,0.9999382354773154E+00,
     *  0.9999421901313457E+00,0.9999494954476229E+00,
     *  0.9999662638389255E+00,0.9999964555571932E+00,
     a  0.0d0/
c 
        data ws13/
     1  0.1611401079789932E-01,0.5243176614272833E-01,
     2  0.9757198720922801E-01,0.1430566289525430E+00,
     3  0.1829696528987093E+00,0.2130500885986832E+00,
     4  0.2304092174407026E+00,0.2335006034909106E+00,
     5  0.2221461707557784E+00,0.1975443463734119E+00,
     6  0.1622560937229564E+00,0.1201804082980457E+00,
     7  0.7656634086096707E-01,0.3806558730134504E-01,
     8  0.1203889395944485E-01,0.1738747468881052E-02,
     9  0.1717070294715531E-03,0.3589153209548273E-04,
     *  0.1311700376394722E-04,0.7057331876018628E-05,
     *  0.5077886927455506E-05,0.4353942087964659E-05,
     *  0.4032996542926421E-05,0.3841956553037945E-05,
     *  0.3695574573653049E-05,0.3567109690156211E-05,
     *  0.3447919896396301E-05,0.3335010414685136E-05,
     *  0.3227037955122892E-05,0.3123283120319233E-05,
     *  0.3023447342926611E-05,0.2927318268940680E-05,
     *  0.2834760661134382E-05,0.2745708018930747E-05,
     *  0.2660150203683463E-05,0.2578587845579162E-05,
     *  0.2502368805207544E-05,0.2435069880082298E-05,
     *  0.2387466006158619E-05,0.2393552164191738E-05,
     *  0.2558135255780179E-05,0.3179884323353273E-05,
     *  0.5056942926955942E-05,0.1050315973363230E-04,
     *  0.2503612595298318E-04,0.2090338508071873E-04,
     a  0.0d0/
  
  
  
c 
        data xs14/
     1  -.9943391303682165E+00,-.9614788003820388E+00,
     2  -.8871448669253550E+00,-.7669953775543426E+00,
     3  -.6037445778176899E+00,-.4051598936659937E+00,
     4  -.1825835241790237E+00,0.5041024077154166E-01,
     5  0.2793694170193652E+00,0.4903381136182261E+00,
     6  0.6712289386791219E+00,0.8131721991077703E+00,
     7  0.9118512698365149E+00,0.9688779788286703E+00,
     8  0.9930003226979434E+00,0.9989414527261056E+00,
     9  0.9996690952163749E+00,0.9997699321811926E+00,
     *  0.9997987490493243E+00,0.9998117643325106E+00,
     *  0.9998199764431051E+00,0.9998265138232092E+00,
     *  0.9998323863713778E+00,0.9998379212234560E+00,
     *  0.9998432265058166E+00,0.9998483415872175E+00,
     *  0.9998532839922032E+00,0.9998580639405330E+00,
     *  0.9998626888712885E+00,0.9998671650070863E+00,
     *  0.9998714980551506E+00,0.9998756933957362E+00,
     *  0.9998797560783117E+00,0.9998836909612385E+00,
     *  0.9998875028993206E+00,0.9998911970345365E+00,
     *  0.9998947795118112E+00,0.9998982593526051E+00,
     *  0.9999016542635100E+00,0.9999050085156207E+00,
     *  0.9999084485138390E+00,0.9999123429060492E+00,
     *  0.9999177156895804E+00,0.9999273930521837E+00,
     *  0.9999492890623149E+00,0.9999904705712346E+00,
     a  0.0d0/
c 
        data ws14/
     1  0.1598676400628937E-01,0.5215037202965101E-01,
     2  0.9718676887233320E-01,0.1426279065242131E+00,
     3  0.1825553525814528E+00,0.2127014038494772E+00,
     4  0.2301661540741270E+00,0.2333897295997255E+00,
     5  0.2221796686138953E+00,0.1977207145254158E+00,
     6  0.1625606023766329E+00,0.1205836406929136E+00,
     7  0.7701954236405138E-01,0.3849429408315635E-01,
     8  0.1234658476483797E-01,0.1866262100226456E-02,
     9  0.2045899534086914E-03,0.4636947956307494E-04,
     *  0.1775139222856338E-04,0.9773569226057931E-05,
     *  0.7105923446944902E-05,0.6117688112972433E-05,
     *  0.5675376082358247E-05,0.5410268256027016E-05,
     *  0.5206213735546044E-05,0.5026670032165686E-05,
     *  0.4859781240155371E-05,0.4701333544239692E-05,
     *  0.4549553339212040E-05,0.4403667165208445E-05,
     *  0.4263326882810885E-05,0.4128192953509400E-05,
     *  0.3997978481913875E-05,0.3872593326303312E-05,
     *  0.3752127826079627E-05,0.3637142495717436E-05,
     *  0.3529232536447035E-05,0.3433155277545725E-05,
     *  0.3363261866717796E-05,0.3363964901859570E-05,
     *  0.3571314997232323E-05,0.4372637205834408E-05,
     *  0.6793392440202207E-05,0.1377494598402248E-04,
     *  0.3278772630083665E-04,0.3380428461727960E-04,
     a  0.0d0/
c 
        data xs15/
     1  -.9944527554024211E+00,-.9620792293323837E+00,
     2  -.8885650957725816E+00,-.7694359126799850E+00,
     3  -.6072613089176417E+00,-.4096704176790649E+00,
     4  -.1878866768463257E+00,0.4460653358845534E-01,
     5  0.2734163258253048E+00,0.4846134943832139E+00,
     6  0.6661052109766126E+00,0.8089804858583520E+00,
     7  0.9088323488622694E+00,0.9671034028047145E+00,
     8  0.9922538775961456E+00,0.9987071371365035E+00,
     9  0.9995539144790007E+00,0.9996795514382215E+00,
     *  0.9997170589717721E+00,0.9997345969905926E+00,
     *  0.9997459207384476E+00,0.9997550439590999E+00,
     *  0.9997632799595284E+00,0.9997710567915323E+00,
     *  0.9997785168645936E+00,0.9997857126354828E+00,
     *  0.9997926677807820E+00,0.9997993960027323E+00,
     *  0.9998059072077690E+00,0.9998122099915139E+00,
     *  0.9998183122334541E+00,0.9998242212737132E+00,
     *  0.9998299441014891E+00,0.9998354875314319E+00,
     *  0.9998408584484979E+00,0.9998460640568509E+00,
     *  0.9998511127738208E+00,0.9998560169857833E+00,
     *  0.9998608008238576E+00,0.9998655245608932E+00,
     *  0.9998703595283362E+00,0.9998758030687451E+00,
     *  0.9998832329577967E+00,0.9998964388648905E+00,
     *  0.9999261447896413E+00,0.9999835305779745E+00,
     a  0.0d0/
c 
        data ws15/
     1  0.1569405098667137E-01,0.5147810056958881E-01,
     2  0.9624277029151886E-01,0.1415555964719376E+00,
     3  0.1814985624565698E+00,0.2117905503178576E+00,
     4  0.2295069831397112E+00,0.2330576618830067E+00,
     5  0.2222174107264295E+00,0.1981390512427457E+00,
     6  0.1633373633766419E+00,0.1216563249432230E+00,
     7  0.7826509657743194E-01,0.3968908957508292E-01,
     8  0.1314742148779711E-01,0.2115777694171182E-02,
     9  0.2501412624202026E-03,0.5936767001689306E-04,
     *  0.2359930244793368E-04,0.1335989146292239E-04,
     *  0.9874826223149834E-05,0.8565598721212880E-05,
     *  0.7969225652932500E-05,0.7605462891162414E-05,
     *  0.7322583372160934E-05,0.7072649233274218E-05,
     *  0.6839839807618621E-05,0.6618213902714525E-05,
     *  0.6405624671466365E-05,0.6201253009434937E-05,
     *  0.6004450781298510E-05,0.5814792453022625E-05,
     *  0.5631993448073646E-05,0.5456012073702955E-05,
     *  0.5287010371306709E-05,0.5125568046245252E-05,
     *  0.4973823978800915E-05,0.4838283344704153E-05,
     *  0.4738380097101526E-05,0.4734507910081520E-05,
     *  0.5010027288999474E-05,0.6086020463593940E-05,
     *  0.9337510447072140E-05,0.1870558641795381E-04,
     *  0.4460393830331468E-04,0.5066754035946578E-04,
     a  0.0d0/
c 
        data xs16/
     1  -.9942738710669105E+00,-.9611601796418081E+00,
     2  -.8864419363222030E+00,-.7658619645020912E+00,
     3  -.6022105357918845E+00,-.4033184948145269E+00,
     4  -.1805744222971716E+00,0.5242119169916371E-01,
     5  0.2812109803208613E+00,0.4918534912605703E+00,
     6  0.6722939753235808E+00,0.8137146182383239E+00,
     7  0.9118763001223043E+00,0.9685058283548127E+00,
     8  0.9924935082984350E+00,0.9985636493775176E+00,
     9  0.9994088554349398E+00,0.9995475972391128E+00,
     *  0.9995914540471973E+00,0.9996136439856249E+00,
     *  0.9996292332997112E+00,0.9996424995353186E+00,
     *  0.9996547623852656E+00,0.9996664346773219E+00,
     *  0.9996776545865815E+00,0.9996884766262177E+00,
     *  0.9996989287229814E+00,0.9997090295817979E+00,
     *  0.9997187942028702E+00,0.9997282357911770E+00,
     *  0.9997373664948797E+00,0.9997461977595816E+00,
     *  0.9997547405376013E+00,0.9997630054511137E+00,
     *  0.9997710030208402E+00,0.9997787440670283E+00,
     *  0.9997862406595226E+00,0.9997935087931222E+00,
     *  0.9998005764735776E+00,0.9998075092686391E+00,
     *  0.9998144922705698E+00,0.9998220851713351E+00,
     *  0.9998319364822508E+00,0.9998488069701454E+00,
     *  0.9998870076154076E+00,0.9999668868575016E+00,
     a  0.0d0/
c 
        data ws16/
     1  0.1615004789119528E-01,0.5248325975074779E-01,
     2  0.9760777642490745E-01,0.1430552046167296E+00,
     3  0.1829186530908038E+00,0.2129450328564640E+00,
     4  0.2302535748586234E+00,0.2333051337521536E+00,
     5  0.2219279307654655E+00,0.1973255075417238E+00,
     6  0.1620636678371875E+00,0.1200480295854493E+00,
     7  0.7653941539051486E-01,0.3820593362011752E-01,
     8  0.1237326707619097E-01,0.2033668488842276E-02,
     9  0.2694142996814896E-03,0.6769981867167573E-04,
     *  0.2871613694173354E-04,0.1772830418841314E-04,
     *  0.1407043781729702E-04,0.1265390194199328E-04,
     *  0.1193191969168123E-04,0.1143310273547255E-04,
     *  0.1101502673091466E-04,0.1063346937474172E-04,
     *  0.1027376500990355E-04,0.9930425185966801E-05,
     *  0.9601007172663466E-05,0.9284188726263817E-05,
     *  0.8979123184673414E-05,0.8685229430585772E-05,
     *  0.8402093275558332E-05,0.8129481232983682E-05,
     *  0.7867447953685868E-05,0.7616618184453401E-05,
     *  0.7379120696597895E-05,0.7161462677015895E-05,
     *  0.6983697629712403E-05,0.6909128853796828E-05,
     *  0.7139576982320536E-05,0.8293658734125402E-05,
     *  0.1211113031161723E-04,0.2376246624947081E-04,
     *  0.5830420799412502E-04,0.8154263557507857E-04,
     a  0.0d0/
c 
        data xs17/
     1  -.9947980938884957E+00,-.9639314882755070E+00,
     2  -.8929734072036359E+00,-.7770234188775325E+00,
     3  -.6181897610585919E+00,-.4236711088217255E+00,
     4  -.2043283313729067E+00,0.2663187325945926E-01,
     5  0.2549994935204163E+00,0.4669322030161671E+00,
     6  0.6503222508292241E+00,0.7961216491175962E+00,
     7  0.8996188771801386E+00,0.9617018053356794E+00,
     8  0.9899477496221689E+00,0.9979447797633170E+00,
     9  0.9991692079719052E+00,0.9993755158821243E+00,
     *  0.9994407723776229E+00,0.9994728268554536E+00,
     *  0.9994943634762818E+00,0.9995121138084326E+00,
     *  0.9995282936356151E+00,0.9995436292496307E+00,
     *  0.9995583640982999E+00,0.9995725889236264E+00,
     *  0.9995863456593508E+00,0.9995996593455627E+00,
     *  0.9996125484945739E+00,0.9996250291634816E+00,
     *  0.9996371163081665E+00,0.9996488240901853E+00,
     *  0.9996601660996360E+00,0.9996711556664095E+00,
     *  0.9996818060363859E+00,0.9996921310101194E+00,
     *  0.9997021469610070E+00,0.9997118780046683E+00,
     *  0.9997213710575583E+00,0.9997307417088452E+00,
     *  0.9997403142857411E+00,0.9997510252668813E+00,
     *  0.9997654653095025E+00,0.9997907678439613E+00,
     *  0.9998473679120106E+00,0.9999594700273391E+00,
     a  0.0d0/
c 
c 
        data ws17/
     1  0.1479834063189082E-01,0.4938580272221483E-01,
     2  0.9329805673002138E-01,0.1382196452541100E+00,
     3  0.1782250160674166E+00,0.2089803589194528E+00,
     4  0.2274815601594297E+00,0.2320476038036737E+00,
     5  0.2223577229416183E+00,0.1994676216820400E+00,
     6  0.1657833286929679E+00,0.1250091955785070E+00,
     7  0.8211250020817376E-01,0.4331538449431646E-01,
     8  0.1554362100955557E-01,0.2888306031743207E-02,
     9  0.3985963890579271E-03,0.1011557605758503E-03,
     *  0.4222873091413611E-04,0.2500502391073524E-04,
     *  0.1905599288758689E-04,0.1677212044934018E-04,
     *  0.1569430478082904E-04,0.1501290638027758E-04,
     *  0.1447048268593719E-04,0.1398572317604512E-04,
     *  0.1353177852166340E-04,0.1309861528667141E-04,
     *  0.1268234111564460E-04,0.1228148626663485E-04,
     *  0.1189516027901054E-04,0.1152266187121407E-04,
     *  0.1116357784741614E-04,0.1081775487013771E-04,
     *  0.1048525371776807E-04,0.1016735025672174E-04,
     *  0.9868359746138565E-05,0.9600908441200321E-05,
     *  0.9402491340183430E-05,0.9386967919351321E-05,
     *  0.9898688365553282E-05,0.1191622199455437E-04,
     *  0.1802353731259170E-04,0.3567261494845622E-04,
     *  0.8511791833744753E-04,0.1072763676540649E-03,
     a  0.0d0/
  
  
        data xs18/
     1  -.9950035009213841E+00,-.9650474124123359E+00,
     2  -.8956388943101541E+00,-.7816096274887228E+00,
     3  -.6247848860541989E+00,-.4321073376957841E+00,
     4  -.2142253795469082E+00,0.1581755554959297E-01,
     5  0.2439220165027369E+00,0.4562997789665633E+00,
     6  0.6408322959468429E+00,0.7883804719858096E+00,
     7  0.8940372207954317E+00,0.9583558224527990E+00,
     8  0.9884204245687125E+00,0.9973688596539452E+00,
     9  0.9988570609923857E+00,0.9991257441643602E+00,
     *  0.9992131134371816E+00,0.9992568649544147E+00,
     *  0.9992867211468155E+00,0.9993115497072482E+00,
     *  0.9993342646761245E+00,0.9993558225980332E+00,
     *  0.9993765466704863E+00,0.9993965585465712E+00,
     *  0.9994159146467405E+00,0.9994346494088866E+00,
     *  0.9994527889425222E+00,0.9994703555919622E+00,
     *  0.9994873699918274E+00,0.9995038517243542E+00,
     *  0.9995198199150699E+00,0.9995352932304391E+00,
     *  0.9995502899998314E+00,0.9995648297716642E+00,
     *  0.9995789358895939E+00,0.9995926423368164E+00,
     *  0.9996060149650967E+00,0.9996192167019415E+00,
     *  0.9996327031776899E+00,0.9996477847805110E+00,
     *  0.9996680805507698E+00,0.9997035521909643E+00,
     *  0.9997826707145174E+00,0.9999393578149547E+00,
     a  0.0d0/
c 
        data ws18/
     1  0.1426220181383326E-01,0.4811708719110026E-01,
     2  0.9151489246197726E-01,0.1362090806437991E+00,
     3  0.1762606846347774E+00,0.2072979870213199E+00,
     4  0.2262696059869694E+00,0.2314435077516763E+00,
     5  0.2224446755752144E+00,0.2002692182367312E+00,
     6  0.1672518719976075E+00,0.1270061873117333E+00,
     7  0.8437946145823837E-01,0.4543563412099401E-01,
     8  0.1697180886050652E-01,0.3398798719035943E-02,
     9  0.5106208493062753E-03,0.1342408754629356E-03,
     *  0.5716807259447701E-04,0.3444453271232865E-04,
     *  0.2656937597906280E-04,0.2351805658247176E-04,
     *  0.2205257592956809E-04,0.2111113839701275E-04,
     *  0.2035540527683267E-04,0.1967714752610001E-04,
     *  0.1904060556226484E-04,0.1843320856833502E-04,
     *  0.1784954620130195E-04,0.1728718837361232E-04,
     *  0.1674486035099345E-04,0.1622179434116388E-04,
     *  0.1571770268296767E-04,0.1523195042079184E-04,
     *  0.1476481298485372E-04,0.1431852541850519E-04,
     *  0.1389908106572461E-04,0.1352386771920470E-04,
     *  0.1324579273795986E-04,0.1322540760317839E-04,
     *  0.1394421350945048E-04,0.1676794058237126E-04,
     *  0.2530340177619964E-04,0.4995516958025669E-04,
     *  0.1187427532750262E-03,0.1532197045985854E-03,
     a  0.0d0/
c 
        data xs19/
     1  -.9951966339090805E+00,-.9661088824850906E+00,
     2  -.8981901512004211E+00,-.7860192773484290E+00,
     3  -.6311552632098202E+00,-.4402996074298553E+00,
     4  -.2238962291371619E+00,0.5175133258910725E-02,
     5  0.2329328413650364E+00,0.4456535356913818E+00,
     6  0.6312193947256772E+00,0.7804127648850141E+00,
     7  0.8881487197666094E+00,0.9546738601223848E+00,
     8  0.9866077537774733E+00,0.9966082766424800E+00,
     9  0.9984217420149313E+00,0.9987739386333176E+00,
     *  0.9988918876354582E+00,0.9989520227153412E+00,
     *  0.9989936110724411E+00,0.9990284546846930E+00,
     *  0.9990604252183259E+00,0.9990907947391177E+00,
     *  0.9991199962305748E+00,0.9991481947139094E+00,
     *  0.9991754681372595E+00,0.9992018649894268E+00,
     *  0.9992274226368152E+00,0.9992521732971806E+00,
     *  0.9992761461375181E+00,0.9992993684830189E+00,
     *  0.9993218667879198E+00,0.9993436668092048E+00,
     *  0.9993647944224529E+00,0.9993852775411904E+00,
     *  0.9994051493944760E+00,0.9994244584233976E+00,
     *  0.9994432990463529E+00,0.9994619049523612E+00,
     *  0.9994809289846933E+00,0.9995022414593817E+00,
     *  0.9995309819671954E+00,0.9995812377314126E+00,
     *  0.9996928626090416E+00,0.9999117733695635E+00,
     a  0.0d0/
c 
        data ws19/
     1  0.1375546256792406E-01,0.4690154595015853E-01,
     2  0.8979548861369277E-01,0.1342565852394687E+00,
     3  0.1743320881253230E+00,0.2056169478252642E+00,
     4  0.2250208944410695E+00,0.2307688236165706E+00,
     5  0.2224368624478499E+00,0.2009632601240992E+00,
     6  0.1686117358011371E+00,0.1289075388002884E+00,
     7  0.8658613586282654E-01,0.4756367118072332E-01,
     8  0.1848975373976336E-01,0.4000863485530274E-02,
     9  0.6573546130799527E-03,0.1795947611137290E-03,
     *  0.7800396057744634E-04,0.4772143912315707E-04,
     *  0.3718826533813532E-04,0.3307018677763055E-04,
     *  0.3105806856706401E-04,0.2974521771797476E-04,
     *  0.2868276983424356E-04,0.2772638657117467E-04,
     *  0.2682828931772669E-04,0.2597155827799417E-04,
     *  0.2514907745708943E-04,0.2435708711523712E-04,
     *  0.2359312641619113E-04,0.2285597747894196E-04,
     *  0.2214491257474729E-04,0.2145940698688033E-04,
     *  0.2080043082314714E-04,0.2017120669501110E-04,
     *  0.1958014240703845E-04,0.1905234055411293E-04,
     *  0.1866390611377620E-04,0.1864501848863166E-04,
     *  0.1968409855480868E-04,0.2372101896334557E-04,
     *  0.3585479450344335E-04,0.7072912453653135E-04,
     *  0.1667873788658514E-03,0.2160875588462321E-03,
     a  0.0d0/
c 
        data xs20/
     1  -.9926438775996931E+00,-.9533105320260309E+00,
     2  -.8689946197718250E+00,-.7371043066652959E+00,
     3  -.5619738232143168E+00,-.3528454644244927E+00,
     4  -.1223521581179309E+00,0.1149297598665073E+00,
     5  0.3439191126075318E+00,0.5504767819433227E+00,
     6  0.7228077318082580E+00,0.8528289274832785E+00,
     7  0.9375974819353456E+00,0.9810036159326900E+00,
     8  0.9953059071265902E+00,0.9978329670223396E+00,
     9  0.9983004462690352E+00,0.9984568791709991E+00,
     *  0.9985393795855064E+00,0.9985994500784978E+00,
     *  0.9986517361330272E+00,0.9987005412397021E+00,
     *  0.9987471336985017E+00,0.9987919237855460E+00,
     *  0.9988350786470068E+00,0.9988766932566326E+00,
     *  0.9989168386869340E+00,0.9989555763577951E+00,
     *  0.9989929628611443E+00,0.9990290518194671E+00,
     *  0.9990638950784210E+00,0.9990975443454702E+00,
     *  0.9991300543347643E+00,0.9991614906704618E+00,
     *  0.9991919551147894E+00,0.9992216697875189E+00,
     *  0.9992512661605457E+00,0.9992827609938535E+00,
     *  0.9993225445042400E+00,0.9993900012878221E+00,
     *  0.9995447851390442E+00,0.9998706821928496E+00,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data ws20/
     1  0.2019707399579139E-01,0.6066126994846471E-01,
     2  0.1083335847167645E+00,0.1546949970127301E+00,
     3  0.1940005237492320E+00,0.2221262343866491E+00,
     4  0.2364079653172650E+00,0.2356188059961511E+00,
     5  0.2199805093422290E+00,0.1911460686580741E+00,
     6  0.1521451485007884E+00,0.1073786454620843E+00,
     7  0.6277310008003004E-01,0.2610289802773415E-01,
     8  0.5705558638987483E-02,0.8810649880819061E-03,
     9  0.2368572818937826E-03,0.1048739520223694E-03,
     *  0.6726668681556300E-04,0.5498903576035382E-04,
     *  0.5019966185785490E-04,0.4759264548296091E-04,
     *  0.4565251216117395E-04,0.4395325567944056E-04,
     *  0.4237165387184646E-04,0.4086921868589956E-04,
     *  0.3943180856059833E-04,0.3805294084418790E-04,
     *  0.3672898271909916E-04,0.3545755060783177E-04,
     *  0.3423753452664730E-04,0.3307001213227317E-04,
     *  0.3196055815448449E-04,0.3092794196497831E-04,
     *  0.3003408200487918E-04,0.2948713717244521E-04,
     *  0.2999833818280665E-04,0.3393740466784467E-04,
     *  0.4847873631570697E-04,0.9528961253772199E-04,
     *  0.2381833408962277E-03,0.3267005842922957E-03,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data xs21/
     1  -.9951040511510197E+00,-.9656687592660109E+00,
     2  -.8973097388184573E+00,-.7848250643235506E+00,
     3  -.6299406599654161E+00,-.4394545028209644E+00,
     4  -.2238357796556922E+00,0.4075320722436881E-02,
     5  0.2303953648240855E+00,0.4415559122965564E+00,
     6  0.6256563106164483E+00,0.7737691208810542E+00,
     7  0.8811730156732283E+00,0.9484717033034527E+00,
     8  0.9822608757496651E+00,0.9941814629658989E+00,
     9  0.9969207494949348E+00,0.9975644750886019E+00,
     *  0.9977925924051301E+00,0.9979110502910251E+00,
     *  0.9979939863300389E+00,0.9980640053272198E+00,
     *  0.9981284412202385E+00,0.9981896798342951E+00,
     *  0.9982485373941228E+00,0.9983053347135067E+00,
     *  0.9983602310749642E+00,0.9984133284914605E+00,
     *  0.9984647060536624E+00,0.9985144318103376E+00,
     *  0.9985625678390844E+00,0.9986091727719113E+00,
     *  0.9986543026808690E+00,0.9986980120791880E+00,
     *  0.9987403553162364E+00,0.9987813907008096E+00,
     *  0.9988211909122453E+00,0.9988598706984211E+00,
     *  0.9988976720700578E+00,0.9989352244357381E+00,
     *  0.9989743193951808E+00,0.9990199912599722E+00,
     *  0.9990853457949855E+00,0.9992028858512292E+00,
     *  0.9994451227052186E+00,0.9998458089087612E+00,
     a  0.0d0/
c 
        data ws21/
     1  0.1398737522286848E-01,0.4733386115356048E-01,
     2  0.9020578557212179E-01,0.1344458561535895E+00,
     3  0.1741675356939313E+00,0.2050384297314148E+00,
     4  0.2240367917854748E+00,0.2294487708974988E+00,
     5  0.2209072582843205E+00,0.1994085578852238E+00,
     6  0.1672837753317238E+00,0.1281358236566034E+00,
     7  0.8675862433336919E-01,0.4894952105897607E-01,
     8  0.2063258093966296E-01,0.5473062646716955E-02,
     9  0.1148836058741209E-02,0.3429597375636260E-03,
     *  0.1526683656158005E-03,0.9467027883595697E-04,
     *  0.7452140778476928E-04,0.6659472868679765E-04,
     *  0.6262739106118034E-04,0.5997107258753510E-04,
     *  0.5779292624759017E-04,0.5582683073545037E-04,
     *  0.5398240561587846E-04,0.5222544427328609E-04,
     *  0.5054094241048406E-04,0.4892087628117273E-04,
     *  0.4736094427533256E-04,0.4585824862271444E-04,
     *  0.4441063895470279E-04,0.4301715517037432E-04,
     *  0.4167887961354003E-04,0.4040362542928915E-04,
     *  0.3921472765399737E-04,0.3818153743765091E-04,
     *  0.3751504170034464E-04,0.3785676978497034E-04,
     *  0.4109839584273283E-04,0.5226270512682273E-04,
     *  0.8346186814821988E-04,0.1635680951002824E-03,
     *  0.3379724749337096E-03,0.3703127159730763E-03,
     a  0.0d0/
c 
        data xs22/
     1  -.9947673264185554E+00,-.9637835336745391E+00,
     2  -.8926742437726787E+00,-.7766149088251633E+00,
     3  -.6177763841640456E+00,-.4233996206978433E+00,
     4  -.2043657991927640E+00,0.2611850217173384E-01,
     5  0.2538611152363313E+00,0.4650522088361875E+00,
     6  0.6476302409809906E+00,0.7926105200266851E+00,
     7  0.8953735798571489E+00,0.9569558981910677E+00,
     8  0.9851854759619925E+00,0.9939296554269335E+00,
     9  0.9960408407671286E+00,0.9966665094453069E+00,
     *  0.9969266861520381E+00,0.9970758662350383E+00,
     *  0.9971868222132341E+00,0.9972832799500480E+00,
     *  0.9973730553074496E+00,0.9974587184191475E+00,
     *  0.9975411709696642E+00,0.9976207842734863E+00,
     *  0.9976977543112233E+00,0.9977722141019816E+00,
     *  0.9978442707500954E+00,0.9979140185748256E+00,
     *  0.9979815443946471E+00,0.9980469303248156E+00,
     *  0.9981102555360548E+00,0.9981715978472623E+00,
     *  0.9982310366423122E+00,0.9982886595368530E+00,
     *  0.9983445800591604E+00,0.9983989891058116E+00,
     *  0.9984523104061380E+00,0.9985056752484551E+00,
     *  0.9985623331009090E+00,0.9986314054343903E+00,
     *  0.9987369765052867E+00,0.9989405473597228E+00,
     *  0.9993550543353598E+00,0.9998606637041936E+00,
     a  0.0d0/
c 
        data ws22/
     1  0.1487571710366096E-01,0.4953286978688549E-01,
     2  0.9344025917998851E-01,0.1382856185891254E+00,
     3  0.1781617114239505E+00,0.2087565308429891E+00,
     4  0.2270875072379329E+00,0.2314928367329240E+00,
     5  0.2216678794524420E+00,0.1986819005176828E+00,
     6  0.1649554091533384E+00,0.1242137908238595E+00,
     7  0.8146349543637163E-01,0.4300433589546448E-01,
     8  0.1588437904630543E-01,0.3927872801255707E-02,
     9  0.1024854645872442E-02,0.3701879078888318E-03,
     *  0.1851154354146285E-03,0.1237996590059079E-03,
     *  0.1016271108251673E-03,0.9243868772938331E-04,
     *  0.8749026488742106E-04,0.8397055575318410E-04,
     *  0.8099154953125317E-04,0.7826616414859795E-04,
     *  0.7569556088760252E-04,0.7324168249809118E-04,
     *  0.7088725712738407E-04,0.6862281764103395E-04,
     *  0.6644248117145325E-04,0.6434255244018140E-04,
     *  0.6232079589417218E-04,0.6037696011997482E-04,
     *  0.5851498190596050E-04,0.5674900603912158E-04,
     *  0.5512168291601978E-04,0.5376038776890740E-04,
     *  0.5305054649113726E-04,0.5416378569935225E-04,
     *  0.6053846234160898E-04,0.8127795072183745E-04,
     *  0.1393925086325846E-03,0.2893816348301497E-03,
     *  0.5349976746692618E-03,0.3600943980974334E-03,
     a  0.0d0/
c 
        data xs23/
     1  -.9949374874662413E+00,-.9647155916033975E+00,
     2  -.8949242148125345E+00,-.7805301546480525E+00,
     3  -.6234723044071621E+00,-.4307737001705397E+00,
     4  -.2131288610000641E+00,0.1640374312972170E-01,
     5  0.2437371103234487E+00,0.4551196951370600E+00,
     6  0.6384938049515028E+00,0.7848101926990658E+00,
     7  0.8892943560853712E+00,0.9527066788076717E+00,
     8  0.9824692340293371E+00,0.9921076183097630E+00,
     9  0.9945789674152377E+00,0.9953490746879472E+00,
     *  0.9956846022119599E+00,0.9958848168126266E+00,
     *  0.9960376618121569E+00,0.9961721755629513E+00,
     *  0.9962979541358641E+00,0.9964181725693579E+00,
     *  0.9965339598181140E+00,0.9966457917496461E+00,
     *  0.9967539271964769E+00,0.9968585464749877E+00,
     *  0.9969597970787025E+00,0.9970578099301977E+00,
     *  0.9971527060203846E+00,0.9972445999434019E+00,
     *  0.9973336023083503E+00,0.9974198220597756E+00,
     *  0.9975033705190763E+00,0.9975843710089956E+00,
     *  0.9976629845781477E+00,0.9977394848730456E+00,
     *  0.9978144843742354E+00,0.9978896210510961E+00,
     *  0.9979695752295958E+00,0.9980673335384616E+00,
     *  0.9982166291726551E+00,0.9985018487828368E+00,
     *  0.9990772800051453E+00,0.9997948891698736E+00,
     a  0.0d0/
c 
        data ws23/
     1  0.1443054869920786E-01,0.4846406977802826E-01,
     2  0.9191014054248410E-01,0.1365243499262420E+00,
     3  0.1763982321340951E+00,0.2071949131226592E+00,
     4  0.2258959358260022E+00,0.2307985957144837E+00,
     5  0.2215534810480780E+00,0.1991801492008999E+00,
     6  0.1660389169019041E+00,0.1257771826677450E+00,
     7  0.8329682209697400E-01,0.4476099449066401E-01,
     8  0.1712681661615264E-01,0.4496977933453570E-02,
     9  0.1237307845120950E-02,0.4680743177038598E-03,
     *  0.2443956862509965E-03,0.1688685360397489E-03,
     *  0.1411288963241849E-03,0.1293114315219923E-03,
     *  0.1227136380830186E-03,0.1178921770019297E-03,
     *  0.1137556282262328E-03,0.1099494513414499E-03,
     *  0.1063508544521356E-03,0.1029120515718834E-03,
     *  0.9961087959550733E-04,0.9643493350139446E-04,
     *  0.9337630868550136E-04,0.9042993795409972E-04,
     *  0.8759289031830349E-04,0.8486502513194128E-04,
     *  0.8225208166966934E-04,0.7977469051905090E-04,
     *  0.7749488883016191E-04,0.7559814386906812E-04,
     *  0.7464479720680494E-04,0.7632821329499292E-04,
     *  0.8556060206716752E-04,0.1151182657737043E-03,
     *  0.1965999602072919E-03,0.4030750404132988E-03,
     *  0.7445646873644774E-03,0.5257452551578952E-03,
     a  0.0d0/
c 
        data xs24/
     1  -.9950929276131520E+00,-.9655790529600533E+00,
     2  -.8970282844419225E+00,-.7842201940874076E+00,
     3  -.6288829979035566E+00,-.4378403423687028E+00,
     4  -.2216122267268278E+00,0.6886666622879734E-02,
     5  0.2336771141086140E+00,0.4450732063546738E+00,
     6  0.6290312851453743E+00,0.7764514414874029E+00,
     7  0.8824284747459601E+00,0.9475004780835832E+00,
     8  0.9787945937603115E+00,0.9895216555013828E+00,
     9  0.9925090893450530E+00,0.9934908215991765E+00,
     *  0.9939338454536046E+00,0.9942061661516017E+00,
     *  0.9944182652705333E+00,0.9946066907083702E+00,
     *  0.9947834867609615E+00,0.9949526565059541E+00,
     *  0.9951156452222840E+00,0.9952730785226118E+00,
     *  0.9954253075650912E+00,0.9955725825627979E+00,
     *  0.9957151100006014E+00,0.9958530738656339E+00,
     *  0.9959866448101575E+00,0.9961159847999426E+00,
     *  0.9962412504348330E+00,0.9963625961253395E+00,
     *  0.9964801798572394E+00,0.9965941777696208E+00,
     *  0.9967048237657767E+00,0.9968125249849187E+00,
     *  0.9969182084550661E+00,0.9970243697686345E+00,
     *  0.9971381130541910E+00,0.9972788725423536E+00,
     *  0.9974960622216587E+00,0.9979085844072360E+00,
     *  0.9987178150868286E+00,0.9997122515739980E+00,
     a  0.0d0/
c 
        data ws24/
     1  0.1402150239203374E-01,0.4746424609127411E-01,
     2  0.9046194979964656E-01,0.1348361880359467E+00,
     3  0.1746785275081920E+00,0.2056318976897863E+00,
     4  0.2246495792113052E+00,0.2299955846385295E+00,
     5  0.2212817594965038E+00,0.1994814430747643E+00,
     6  0.1668973791073950E+00,0.1271047120065679E+00,
     7  0.8491783066541712E-01,0.4640867108159857E-01,
     8  0.1849014440295712E-01,0.5276079369089038E-02,
     9  0.1551029895694315E-02,0.6092026001810326E-03,
     *  0.3282682355628004E-03,0.2325532717682409E-03,
     *  0.1970611521058886E-03,0.1815617998595978E-03,
     *  0.1726209586984278E-03,0.1659354616679799E-03,
     *  0.1601387190892622E-03,0.1547839996720963E-03,
     *  0.1497149914102665E-03,0.1448691081587209E-03,
     *  0.1402162950897415E-03,0.1357398080080069E-03,
     *  0.1314290219779013E-03,0.1272770098860868E-03,
     *  0.1232799129864851E-03,0.1194375836170509E-03,
     *  0.1157586766035021E-03,0.1122745046127468E-03,
     *  0.1090804975762903E-03,0.1064616870146367E-03,
     *  0.1052770600233365E-03,0.1081105102803362E-03,
     *  0.1223520914196917E-03,0.1667690194756243E-03,
     *  0.2863509436010686E-03,0.5774529998739252E-03,
     *  0.1032304770697933E-02,0.7350707447859947E-03,
     a  0.0d0/
c 
c 
        data xs25/
     1  -.9952103811404152E+00,-.9662438087556234E+00,
     2  -.8986745685208340E+00,-.7871529773594461E+00,
     3  -.6332539122702368E+00,-.4436496357085439E+00,
     4  -.2287203298392633E+00,-.1257822465520531E-02,
     5  0.2248583075964359E+00,0.4360111868454903E+00,
     6  0.6201861390304559E+00,0.7682615952865568E+00,
     7  0.8752514405547389E+00,0.9415594669921663E+00,
     8  0.9741663288473799E+00,0.9860273764589380E+00,
     9  0.9896344142887756E+00,0.9908881525039008E+00,
     *  0.9914745530122273E+00,0.9918459150549157E+00,
     *  0.9921409385396949E+00,0.9924054229669861E+00,
     *  0.9926543877043639E+00,0.9928928504461585E+00,
     *  0.9931226582147494E+00,0.9933446365740464E+00,
     *  0.9935592637873116E+00,0.9937668880789836E+00,
     *  0.9939678008936070E+00,0.9941622635323748E+00,
     *  0.9943505184458867E+00,0.9945327957014474E+00,
     *  0.9947093181999644E+00,0.9948803068432903E+00,
     *  0.9950459900468973E+00,0.9952066284055623E+00,
     *  0.9953625812695692E+00,0.9955144992177349E+00,
     *  0.9956639062046602E+00,0.9958149581753488E+00,
     *  0.9959794747476565E+00,0.9961893806584938E+00,
     *  0.9965243403257644E+00,0.9971625786067164E+00,
     *  0.9983311149616552E+00,0.9996376073712174E+00,
     a  0.0d0/
c 
        data ws25/
     1  0.1371023132607394E-01,0.4668312724950054E-01,
     2  0.8930305192824227E-01,0.1334480007869598E+00,
     3  0.1732158242310376E+00,0.2042412339980767E+00,
     4  0.2234632236793640E+00,0.2291258932563945E+00,
     5  0.2208144848003660E+00,0.1994673014676816E+00,
     6  0.1673418197186921E+00,0.1279573740531807E+00,
     7  0.8606759399630443E-01,0.4770640654573670E-01,
     8  0.1978395012042153E-01,0.6167564757783217E-02,
     9  0.1946057653364518E-02,0.7946578734183881E-03,
     *  0.4420350228334954E-03,0.3210602568959724E-03,
     *  0.2757611554069425E-03,0.2554151774711137E-03,
     *  0.2432556582775661E-03,0.2339480966424542E-03,
     *  0.2257957361004434E-03,0.2182376655665978E-03,
     *  0.2110737644282216E-03,0.2042230133047291E-03,
     *  0.1976459994951277E-03,0.1913195865825224E-03,
     *  0.1852284979239722E-03,0.1793630840832828E-03,
     *  0.1737185931190801E-03,0.1682964010353653E-03,
     *  0.1631124850386257E-03,0.1582214345550831E-03,
     *  0.1537845767699324E-03,0.1502805150476088E-03,
     *  0.1491524611125059E-03,0.1547282209126265E-03,
     *  0.1792094283838350E-03,0.2529353922962755E-03,
     *  0.4461090063192698E-03,0.8780745989784685E-03,
     *  0.1416598482849384E-02,0.9335142105354574E-03,
     a  0.0d0/
  
c 
        data xs26/
     1  -.9953291967699994E+00,-.9669238785392430E+00,
     2  -.9003743012902230E+00,-.7902070627261963E+00,
     3  -.6378450429968741E+00,-.4498060591922696E+00,
     4  -.2363226171679226E+00,-.1005186568882110E-01,
     5  0.2152404721149954E+00,0.4260213503654653E+00,
     6  0.6103191404105166E+00,0.7590032996796428E+00,
     7  0.8670161413644162E+00,0.9346200405775575E+00,
     8  0.9685805969360086E+00,0.9815269199830359E+00,
     9  0.9857319788433723E+00,0.9872784972603286E+00,
     *  0.9880383185573880E+00,0.9885395981359388E+00,
     *  0.9889479677889919E+00,0.9893181074837831E+00,
     *  0.9896679093219367E+00,0.9900034064936356E+00,
     *  0.9903268756680366E+00,0.9906393756779779E+00,
     *  0.9909415450209900E+00,0.9912338597192709E+00,
     *  0.9915167225811263E+00,0.9917904981054961E+00,
     *  0.9920555267713431E+00,0.9923121329417193E+00,
     *  0.9925606318455553E+00,0.9928013380923344E+00,
     *  0.9930345825296245E+00,0.9932607541550154E+00,
     *  0.9934804150532117E+00,0.9936946374219000E+00,
     *  0.9939060218228550E+00,0.9941217661423281E+00,
     *  0.9943623331012417E+00,0.9946825595909601E+00,
     *  0.9952189641734681E+00,0.9962556754284539E+00,
     *  0.9979799741445605E+00,0.9995964127215438E+00,
     a  0.0d0/
c 
c 
        data ws26/
     1  0.1339397218831067E-01,0.4587725400495118E-01,
     2  0.8809171889796823E-01,0.1319763124742280E+00,
     3  0.1716392909341415E+00,0.2027118986989270E+00,
     4  0.2221234302011706E+00,0.2281013038564317E+00,
     5  0.2202056438357740E+00,0.1993395014060147E+00,
     6  0.1677139762748210E+00,0.1287908335854329E+00,
     7  0.8725651840598900E-01,0.4907406842229829E-01,
     8  0.2106834141060372E-01,0.7004768100283533E-02,
     9  0.2348671618350209E-02,0.1008334501978650E-02,
     *  0.5866582907666720E-03,0.4403109931640915E-03,
     *  0.3845028441147800E-03,0.3584053556525387E-03,
     *  0.3420916433350655E-03,0.3292472405727792E-03,
     *  0.3178565704392074E-03,0.3072462784890584E-03,
     *  0.2971704790678185E-03,0.2875255273636920E-03,
     *  0.2782607416999346E-03,0.2693469275054900E-03,
     *  0.2607643463992027E-03,0.2525004025033518E-03,
     *  0.2445495898784776E-03,0.2369174372048276E-03,
     *  0.2296346456995486E-03,0.2227975623759407E-03,
     *  0.2166886879403911E-03,0.2121461274923863E-03,
     *  0.2116943504070452E-03,0.2228783425810255E-03,
     *  0.2666802856052077E-03,0.3949554014598967E-03,
     *  0.7266150483034323E-03,0.1399016087689895E-02,
     *  0.1914973987461759E-02,0.1062637880496613E-02,
     a  0.0d0/
  
  
c 
        data xs27/
     1  -.9955104390996922E+00,-.9679463381194869E+00,
     2  -.9028844710645931E+00,-.7946329862477439E+00,
     3  -.6443721930532729E+00,-.4583905632964291E+00,
     4  -.2467145575688032E+00,-.2182542457406944E-01,
     5  0.2026464883478112E+00,0.4132497979148055E+00,
     6  0.5980252769498773E+00,0.7477695180838888E+00,
     7  0.8572509247359053E+00,0.9264479771075533E+00,
     8  0.9617824695307733E+00,0.9756632905334361E+00,
     9  0.9804062138424350E+00,0.9822644950551020E+00,
     *  0.9832361183963143E+00,0.9839093047507086E+00,
     *  0.9844731560607221E+00,0.9849901721323335E+00,
     *  0.9854808239564623E+00,0.9859521124919164E+00,
     *  0.9864067638836591E+00,0.9868461067969412E+00,
     *  0.9872709794592223E+00,0.9876820257646231E+00,
     *  0.9880797983851151E+00,0.9884647987333266E+00,
     *  0.9888374960224384E+00,0.9891983391121220E+00,
     *  0.9895477664435057E+00,0.9898862187369315E+00,
     *  0.9902141665802600E+00,0.9905321826786870E+00,
     *  0.9908411438586637E+00,0.9911428273718283E+00,
     *  0.9914417128479107E+00,0.9917503740361897E+00,
     *  0.9921043458492249E+00,0.9925976475661735E+00,
     *  0.9934627010252346E+00,0.9951417862725489E+00,
     *  0.9976516007848020E+00,0.9995782222454819E+00,
     a  0.0d0/
c 
        data ws27/
     1  0.1291356287483142E-01,0.4468283390847467E-01,
     2  0.8634949589596620E-01,0.1299279465546809E+00,
     3  0.1695235467682363E+00,0.2007473915119682E+00,
     4  0.2205025886631954E+00,0.2269817057035191E+00,
     5  0.2196982800595600E+00,0.1994940600933744E+00,
     6  0.1685029063871334E+00,0.1300933170498692E+00,
     7  0.8884336252137762E-01,0.5061273640531937E-01,
     8  0.2223741902929324E-01,0.7715245129388052E-02,
     9  0.2745761231317028E-02,0.1255624972693383E-02,
     *  0.7724332413063580E-03,0.6019098855499914E-03,
     *  0.5350097017862664E-03,0.5020314351266698E-03,
     *  0.4803120315359437E-03,0.4626798014563763E-03,
     *  0.4468319036510730E-03,0.4319894796409390E-03,
     *  0.4178616868210019E-03,0.4043224472800318E-03,
     *  0.3913060623018682E-03,0.3787726460638277E-03,
     *  0.3666965783817861E-03,0.3550624918326866E-03,
     *  0.3438652279459324E-03,0.3331167517694638E-03,
     *  0.3228715367945171E-03,0.3132984424757620E-03,
     *  0.3048948836062162E-03,0.2991421997888699E-03,
     *  0.3005021478136083E-03,0.3221982111408917E-03,
     *  0.4001904952431944E-03,0.6226497579301408E-03,
     *  0.1187636428611275E-02,0.2210068254904190E-02,
     *  0.2523510323021894E-02,0.1140994328898936E-02,
     a  0.0d0/
  
c 
        data xs28/
     1  -.9953531170792721E+00,-.9670962571717365E+00,
     2  -.9008997917624749E+00,-.7913275934891781E+00,
     3  -.6398042544724982E+00,-.4528185489591424E+00,
     4  -.2405487024603478E+00,-.1557997786391530E-01,
     5  0.2084020904516022E+00,0.4179421453952525E+00,
     6  0.6011297202613015E+00,0.7488570664921805E+00,
     7  0.8560175657517546E+00,0.9227574839023862E+00,
     8  0.9559698886242141E+00,0.9688630164827434E+00,
     9  0.9735157328870478E+00,0.9755358629679971E+00,
     *  0.9767148855439765E+00,0.9776055666102386E+00,
     *  0.9783872656317998E+00,0.9791175872576059E+00,
     *  0.9798151285486449E+00,0.9804864197608345E+00,
     *  0.9811342244501695E+00,0.9817600686107364E+00,
     *  0.9823650279786960E+00,0.9829499864172924E+00,
     *  0.9835157299325605E+00,0.9840629867085067E+00,
     *  0.9845924476137140E+00,0.9851047816803609E+00,
     *  0.9856006561610215E+00,0.9860807743431828E+00,
     *  0.9865459711718867E+00,0.9869974863842337E+00,
     *  0.9874377905301053E+00,0.9878731562972638E+00,
     *  0.9883216019326573E+00,0.9888354289528485E+00,
     *  0.9895553828686257E+00,0.9908242454287390E+00,
     *  0.9932378112784477E+00,0.9967397861132132E+00,
     *  0.9994203669336863E+00,0.0d0,0.0d0/
c 
        data ws28/
     1  0.1332484527919674E-01,0.4563593386467120E-01,
     2  0.8762072211076506E-01,0.1312570175048774E+00,
     3  0.1706861546302794E+00,0.2015676577434214E+00,
     4  0.2208526240847442E+00,0.2267815985441695E+00,
     5  0.2189179774527111E+00,0.1981560580324014E+00,
     6  0.1666813810970625E+00,0.1279021644319423E+00,
     7  0.8641716098183463E-01,0.4819260052679081E-01,
     8  0.2059168820079687E-01,0.7316139001690820E-02,
     9  0.2846179246569441E-02,0.1451750641095387E-02,
     *  0.9870839018230246E-03,0.8207334265143718E-03,
     *  0.7510886976162548E-03,0.7122812157703858E-03,
     *  0.6837764891186578E-03,0.6592319656119930E-03,
     *  0.6366208894762251E-03,0.6152430831677299E-03,
     *  0.5948215932513228E-03,0.5752256773295741E-03,
     *  0.5563824128349095E-03,0.5382460617305882E-03,
     *  0.5207869606869138E-03,0.5039918507370552E-03,
     *  0.4878729761311804E-03,0.4724987144931455E-03,
     *  0.4580900915567569E-03,0.4453164315795208E-03,
     *  0.4362256291787169E-03,0.4371741425282539E-03,
     *  0.4675948707003631E-03,0.5816986147091296E-03,
     *  0.9123137420052798E-03,0.1736682085800555E-02,
     *  0.3121297354431232E-02,0.3509737635981101E-02,
     *  0.1578246057646779E-02,0.0d0,0.0d0/
c 
        data xs29/
     1  -.9958404486938160E+00,-.9698455005830877E+00,
     2  -.9076123078553266E+00,-.8030636185527186E+00,
     3  -.6569356275399555E+00,-.4750886220721297E+00,
     4  -.2671576805167347E+00,-.4528842999421012E-01,
     5  0.1771417309598526E+00,0.3868202512089498E+00,
     6  0.5717802256992930E+00,0.7226343552439323E+00,
     7  0.8337708319014392E+00,0.9046518936538390E+00,
     8  0.9414338341270128E+00,0.9567277082131853E+00,
     9  0.9626734746298207E+00,0.9653903525594709E+00,
     *  0.9670117966085645E+00,0.9682412353982892E+00,
     *  0.9693185599004954E+00,0.9703241335760184E+00,
     *  0.9712847062900178E+00,0.9722097748163783E+00,
     *  0.9731032944255786E+00,0.9739673738644659E+00,
     *  0.9748034557346754E+00,0.9756127080666579E+00,
     *  0.9763961635526157E+00,0.9771547752999479E+00,
     *  0.9778894441474903E+00,0.9786010366275228E+00,
     *  0.9792904035018082E+00,0.9799584102726084E+00,
     *  0.9806060056594368E+00,0.9812344040031826E+00,
     *  0.9818456102680384E+00,0.9824439857739569E+00,
     *  0.9830409696217008E+00,0.9836689297852235E+00,
     *  0.9844173744641862E+00,0.9855142658024096E+00,
     *  0.9874865636997161E+00,0.9911331952494036E+00,
     *  0.9960120487683937E+00,0.9993434413901838E+00,
     a  0.0d0/
c 
        data ws29/
     1  0.1203164177616014E-01,0.4243331841088184E-01,
     2  0.8301042264201379E-01,0.1259364219655384E+00,
     3  0.1653210646093062E+00,0.1967461780183616E+00,
     4  0.2170700952254632E+00,0.2244171306941097E+00,
     5  0.2182030407690249E+00,0.1991389423337117E+00,
     6  0.1691967771149833E+00,0.1315555296410502E+00,
     7  0.9059924736711072E-01,0.5216848013754961E-01,
     8  0.2364191581554106E-01,0.9084136754449001E-02,
     9  0.3762247335174795E-02,0.1984452700213945E-02,
     *  0.1362542780835627E-02,0.1132087883888316E-02,
     *  0.1034422450851651E-02,0.9806715136188471E-03,
     *  0.9418947178623309E-03,0.9088505406957050E-03,
     *  0.8785237392399330E-03,0.8498700596945763E-03,
     *  0.8224863509659829E-03,0.7961896090085039E-03,
     *  0.7708798018688795E-03,0.7464935647370516E-03,
     *  0.7229882774257586E-03,0.7003380260678326E-03,
     *  0.6785389841040470E-03,0.6576300798559494E-03,
     *  0.6377556827313395E-03,0.6193532098516070E-03,
     *  0.6037198520846221E-03,0.5947384701683890E-03,
     *  0.6040696532453822E-03,0.6655164327561068E-03,
     *  0.8666579974816876E-03,0.1411591933074453E-02,
     *  0.2689690410506988E-02,0.4578598530484366E-02,
     *  0.4613242612399429E-02,0.1829584054061317E-02,
     a  0.0d0/
  
  
c 
        data xs30/
     1  -.9961835958615686E+00,-.9718795942596589E+00,
     2  -.9127684068384633E+00,-.8123552037859059E+00,
     3  -.6708541651840542E+00,-.4936021193809852E+00,
     4  -.2897475382516094E+00,-.7102239805770039E-01,
     5  0.1494979186319092E+00,0.3586559978995192E+00,
     6  0.5444649659047174E+00,0.6973181769836305E+00,
     7  0.8111618635909716E+00,0.8848232572920416E+00,
     8  0.9238953554394500E+00,0.9408342399062666E+00,
     9  0.9478725910120163E+00,0.9513118958999003E+00,
     *  0.9534613198439714E+00,0.9551282696384140E+00,
     *  0.9566008485483222E+00,0.9579791528984376E+00,
     *  0.9592975786139208E+00,0.9605686137297786E+00,
     *  0.9617974941987852E+00,0.9629869899223251E+00,
     *  0.9641389774181013E+00,0.9652549711533640E+00,
     *  0.9663363129878673E+00,0.9673842464052766E+00,
     *  0.9683999507959384E+00,0.9693845622381270E+00,
     *  0.9703391924907729E+00,0.9712649564291698E+00,
     *  0.9721630288984144E+00,0.9730347868395575E+00,
     *  0.9738821959054764E+00,0.9747089124986459E+00,
     *  0.9755234886144744E+00,0.9763486568282760E+00,
     *  0.9772467640614251E+00,0.9783799777454036E+00,
     *  0.9801370895369654E+00,0.9833391090073491E+00,
     *  0.9888790496833069E+00,0.9953820936280999E+00,
     *  0.9992885137887741E+00/
  
  
        data ws30/
     1  0.1110318930235963E-01,0.3997690396516652E-01,
     2  0.7929777651282971E-01,0.1214664713703592E+00,
     3  0.1606248582304631E+00,0.1923341934089972E+00,
     4  0.2133989636376421E+00,0.2218509657231246E+00,
     5  0.2169758061353805E+00,0.1993223528605134E+00,
     6  0.1706752403645005E+00,0.1340016180947645E+00,
     7  0.9345760833649173E-01,0.5476790654306071E-01,
     8  0.2560346564830215E-01,0.1042873913346403E-01,
     9  0.4635670740877079E-02,0.2586163969888131E-02,
     *  0.1833650140360112E-02,0.1543552673337895E-02,
     *  0.1416493468753922E-02,0.1345223868693455E-02,
     *  0.1293503730256550E-02,0.1249372396561664E-02,
     *  0.1208828171883288E-02,0.1170468523085827E-02,
     *  0.1133756163315381E-02,0.1098453866416893E-02,
     *  0.1064436534612588E-02,0.1031626546444343E-02,
     *  0.9999712847524509E-03,0.9694364738680679E-03,
     *  0.9400090045214994E-03,0.9117126141073129E-03,
     *  0.8846562937550135E-03,0.8591715596459903E-03,
     *  0.8362105703993791E-03,0.8185152882299108E-03,
     *  0.8140485750968744E-03,0.8457344260734337E-03,
     *  0.9755642985509003E-03,0.1351576534098974E-02,
     *  0.2298941332975711E-02,0.4301361671817030E-02,
     *  0.6573125909031834E-02,0.5705039781911453E-02,
     *  0.2017638396973912E-02/
  
c 
        data xs31/
     1  -.9963207644720740E+00,-.9727550713801462E+00,
     2  -.9151536097556596E+00,-.8169643265950007E+00,
     3  -.6782526000546915E+00,-.5041589595190312E+00,
     4  -.3036081976403223E+00,-.8810012430111835E-01,
     5  0.1295112085029277E+00,0.3362616940790109E+00,
     6  0.5203176810088461E+00,0.6721804545539514E+00,
     7  0.7858643325865151E+00,0.8602339338077858E+00,
     8  0.9007577273595242E+00,0.9192957762711054E+00,
     9  0.9275520852361587E+00,0.9318606812372741E+00,
     *  0.9347016188290829E+00,0.9369822363398423E+00,
     *  0.9390297606081623E+00,0.9409578478874221E+00,
     *  0.9428059821998578E+00,0.9445889249469819E+00,
     *  0.9463131729256011E+00,0.9479823503549367E+00,
     *  0.9495990009971771E+00,0.9511652063875238E+00,
     *  0.9526828143888220E+00,0.9541535352103190E+00,
     *  0.9555789887421962E+00,0.9569607340217205E+00,
     *  0.9583002973585948E+00,0.9595992145361041E+00,
     *  0.9608591209819854E+00,0.9620819807834746E+00,
     *  0.9632707201539128E+00,0.9644310448513109E+00,
     *  0.9655767112610998E+00,0.9667445845604363E+00,
     *  0.9680344500502309E+00,0.9696979970667714E+00,
     *  0.9723150683870265E+00,0.9770274836675890E+00,
     *  0.9848636222410239E+00,0.9937503244321670E+00,
     *  0.9990395650212009E+00/
c 
        data ws31/
     1  0.1072258486818811E-01,0.3885525126568586E-01,
     2  0.7741150011012347E-01,0.1189303111507443E+00,
     3  0.1576153996381504E+00,0.1890644423237337E+00,
     4  0.2101000933423899E+00,0.2187528748113498E+00,
     5  0.2142883342844222E+00,0.1972202683876735E+00,
     6  0.1692892832508026E+00,0.1334177010313564E+00,
     7  0.9372914563062615E-01,0.5590165532012459E-01,
     8  0.2726166256769025E-01,0.1186922026787726E-01,
     9  0.5648437602908171E-02,0.3342871466111524E-02,
     *  0.2476062920246847E-02,0.2134344160832068E-02,
     *  0.1977603167676653E-02,0.1884450322963862E-02,
     *  0.1814048330447787E-02,0.1752844634418832E-02,
     *  0.1696230173299777E-02,0.1642539864018262E-02,
     *  0.1591105048235866E-02,0.1541612279721559E-02,
     *  0.1493888041000646E-02,0.1447823138778191E-02,
     *  0.1403343328517988E-02,0.1360400986711245E-02,
     *  0.1318980325763931E-02,0.1279123309376825E-02,
     *  0.1241007910404752E-02,0.1205173607984333E-02,
     *  0.1173181201859904E-02,0.1149544091621617E-02,
     *  0.1147332404407553E-02,0.1203629955561815E-02,
     *  0.1415304918698834E-02,0.2003161179540438E-02,
     *  0.3424454479583007E-02,0.6236413679855494E-02,
     *  0.9104382376912462E-02,0.7733273451516578E-02,
     *  0.2727671072473307E-02/
c 
        data xs32/
     1  -.9960850152568159E+00,-.9713642947026677E+00,
     2  -.9116838300290891E+00,-.8108702673035235E+00,
     3  -.6694490553599235E+00,-.4930180177539323E+00,
     4  -.2909245574299754E+00,-.7503048510309318E-01,
     5  0.1415284680098026E+00,0.3455940323824485E+00,
     6  0.5252223141799726E+00,0.6709071084767939E+00,
     7  0.7769324252449614E+00,0.8432991054772763E+00,
     8  0.8778554497246858E+00,0.8937874127876742E+00,
     9  0.9015795888028562E+00,0.9062593200431857E+00,
     *  0.9097680477615029E+00,0.9128159497103945E+00,
     *  0.9156497431593531E+00,0.9183529157045081E+00,
     *  0.9209550996137763E+00,0.9234684220674307E+00,
     *  0.9258991617562271E+00,0.9282515011217225E+00,
     *  0.9305287704182179E+00,0.9327338859711658E+00,
     *  0.9348695215590728E+00,0.9369381883957766E+00,
     *  0.9389422850941913E+00,0.9408841439788582E+00,
     *  0.9427660986494397E+00,0.9445906240561244E+00,
     *  0.9463606912601743E+00,0.9480807559287296E+00,
     *  0.9497596381591609E+00,0.9514190571793522E+00,
     *  0.9531185780239442E+00,0.9550225185345358E+00,
     *  0.9575492394748688E+00,0.9616516264578525E+00,
     *  0.9690099650470818E+00,0.9804030689629881E+00,
     *  0.9921864610868736E+00,0.9988213369848569E+00,
     a  0.0d0/
c 
        data ws32/
     1  0.1136060888886296E-01,0.4051989948976139E-01,
     2  0.7983757663049135E-01,0.1216749869854466E+00,
     3  0.1602130413816060E+00,0.1910706989617200E+00,
     4  0.2111185565843142E+00,0.2184593469527335E+00,
     5  0.2124448770605909E+00,0.1936754561179633E+00,
     6  0.1639753242273313E+00,0.1264109305248668E+00,
     7  0.8559193192729549E-01,0.4847480136960413E-01,
     8  0.2302081333853095E-01,0.1062356690820499E-01,
     9  0.5742788018464034E-02,0.3911076278722807E-02,
     *  0.3213491179549988E-02,0.2918973653169005E-02,
     *  0.2760965365378188E-02,0.2649811787752082E-02,
     *  0.2556415678115419E-02,0.2471223490474273E-02,
     *  0.2390935479364350E-02,0.2314292006507753E-02,
     *  0.2240730106047790E-02,0.2169945401162598E-02,
     *  0.2101743517247904E-02,0.2035989204041902E-02,
     *  0.1972591896648853E-02,0.1911512769248358E-02,
     *  0.1852803990740945E-02,0.1796727998294654E-02,
     *  0.1744108432938396E-02,0.1697370552960580E-02,
     *  0.1663656766141572E-02,0.1664099191784153E-02,
     *  0.1759947457969427E-02,0.2113099799786991E-02,
     *  0.3092134659760037E-02,0.5411570347218005E-02,
     *  0.9525892666385807E-02,0.1260757865505377E-01,
     *  0.9874825372985530E-02,0.3361247900705708E-02,
     a  0.0d0/
c 
        data xs33/
     1  -.9962287983460283E+00,-.9722470501639205E+00,
     2  -.9140223821098656E+00,-.8152956814796113E+00,
     3  -.6764423475183157E+00,-.5028902072215747E+00,
     4  -.3038180229923926E+00,-.9094359890254611E-01,
     5  0.1226864871255460E+00,0.3239096660144243E+00,
     6  0.5006391668033809E+00,0.6430704491790502E+00,
     7  0.7452349868039272E+00,0.8077061791600020E+00,
     8  0.8399508775661020E+00,0.8555569644731384E+00,
     9  0.8640147503727427E+00,0.8696854390249387E+00,
     *  0.8742800334451624E+00,0.8784244049647170E+00,
     *  0.8823340792827823E+00,0.8860830563670397E+00,
     *  0.8896989045207626E+00,0.8931939297100984E+00,
     *  0.8965752095396113E+00,0.8998478744048097E+00,
     *  0.9030162262537689E+00,0.9060841492186262E+00,
     *  0.9090552793818185E+00,0.9119330904058205E+00,
     *  0.9147209523286006E+00,0.9174221943534049E+00,
     *  0.9200402140562041E+00,0.9225787365996467E+00,
     *  0.9250425142654659E+00,0.9274393215034609E+00,
     *  0.9297858059781641E+00,0.9321247194433671E+00,
     *  0.9345738655020996E+00,0.9374479778959324E+00,
     *  0.9415036288601550E+00,0.9483114379683826E+00,
     *  0.9599879450950286E+00,0.9760802203505734E+00,
     *  0.9908959725986519E+00,0.9986611172100199E+00,
     a  0.0d0/
c 
        data ws33/
     1  0.1096766187539867E-01,0.3941788248259565E-01,
     2  0.7804432483646828E-01,0.1193201369865604E+00,
     3  0.1574607685385605E+00,0.1880936719475171E+00,
     4  0.2080775417050222E+00,0.2154765845191520E+00,
     5  0.2095699916521430E+00,0.1908423752057071E+00,
     6  0.1609590565538101E+00,0.1228643651421676E+00,
     7  0.8149558312038253E-01,0.4511001922661589E-01,
     8  0.2182680872859787E-01,0.1095971096073970E-01,
     9  0.6632113544083238E-02,0.4970425537704509E-02,
     *  0.4312863901691169E-02,0.4007837953002052E-02,
     *  0.3822499169344290E-02,0.3679616040406393E-02,
     *  0.3553991813879168E-02,0.3437192223057865E-02,
     *  0.3326206606249994E-02,0.3219834140941916E-02,
     *  0.3117514372458848E-02,0.3018936241238288E-02,
     *  0.2923902098418026E-02,0.2832280646391822E-02,
     *  0.2743996431763700E-02,0.2659050106576207E-02,
     *  0.2577600134530985E-02,0.2500208732402935E-02,
     *  0.2428559741196422E-02,0.2367584998082709E-02,
     *  0.2331796743787005E-02,0.2363761763452415E-02,
     *  0.2583019651683875E-02,0.3284410128248637E-02,
     *  0.5085092120989761E-02,0.8937955684552593E-02,
     *  0.1442061812624145E-01,0.1669998389701246E-01,
     *  0.1183629613884086E-01,0.3838337775741875E-02,
     a  0.0d0/
c 
        data xs34/
     1  -.9962210550147963E+00,-.9723311174140808E+00,
     2  -.9145979419798353E+00,-.8170164070393446E+00,
     3  -.6800955705421852E+00,-.5092884722196467E+00,
     4  -.3137170731493151E+00,-.1049899298957913E+00,
     5  0.1039732274546312E+00,0.3001064273333476E+00,
     6  0.4713269115806065E+00,0.6077187507885136E+00,
     7  0.7033794173028844E+00,0.7600346845471160E+00,
     8  0.7890325986626547E+00,0.8040058546542281E+00,
     9  0.8132292264038192E+00,0.8202542563608179E+00,
     *  0.8264093933287278E+00,0.8321510264262979E+00,
     *  0.8376319249752386E+00,0.8429065845402035E+00,
     *  0.8479979829799087E+00,0.8529185898222285E+00,
     *  0.8576770452551556E+00,0.8622803661135277E+00,
     *  0.8667347300072001E+00,0.8710457938290876E+00,
     *  0.8752188531872642E+00,0.8792589532193396E+00,
     *  0.8831710098031738E+00,0.8869600075740880E+00,
     *  0.8906314355494467E+00,0.8941924395141039E+00,
     *  0.8976551385562792E+00,0.9010465235204630E+00,
     *  0.9044381873678828E+00,0.9080318302896748E+00,
     *  0.9123703122537157E+00,0.9187488560378858E+00,
     *  0.9297059578369802E+00,0.9477845571826928E+00,
     *  0.9703975535803278E+00,0.9892363420040846E+00,
     *  0.9984605101813473E+00,0.0d0,0.0d0/
c 
        data ws34/
     1  0.1096963075621450E-01,0.3917778418753773E-01,
     2  0.7725839454541951E-01,0.1177937467090502E+00,
     3  0.1551184225420143E+00,0.1849558841370124E+00,
     4  0.2042325114398812E+00,0.2110489470475969E+00,
     5  0.2046809514774586E+00,0.1855530442833493E+00,
     6  0.1552086052955728E+00,0.1165284701463707E+00,
     7  0.7500936938876110E-01,0.4035131907019439E-01,
     8  0.2007227882318971E-01,0.1121292978155313E-01,
     9  0.7781038369171537E-02,0.6469544960204133E-02,
     *  0.5908292797531340E-02,0.5597522042558075E-02,
     *  0.5372433987715412E-02,0.5180425411593971E-02,
     *  0.5004357303771848E-02,0.4838265571161535E-02,
     *  0.4679803612311525E-02,0.4527863013100261E-02,
     *  0.4381805135304111E-02,0.4241203336595208E-02,
     *  0.4105754277504806E-02,0.3975262032902643E-02,
     *  0.3849675529930400E-02,0.3729218412418321E-02,
     *  0.3614781812654976E-02,0.3509107185068921E-02,
     *  0.3420370252409308E-02,0.3373086721940642E-02,
     *  0.3440536478341856E-02,0.3830490686459511E-02,
     *  0.5051457984895494E-02,0.8135508402703785E-02,
     *  0.1431447166172885E-01,0.2144508194994752E-01,
     *  0.2221084353712943E-01,0.1439852638377590E-01,
     *  0.4440947535721182E-02,0.0d0,0.0d0/
c 
        data xs35/
     1  -.9966510750991830E+00,-.9749720241592498E+00,
     2  -.9215329278657105E+00,-.8299609886385974E+00,
     3  -.7002138866656071E+00,-.5371616972498704E+00,
     4  -.3493781814688475E+00,-.1480234552407970E+00,
     5  0.5425534750287642E-01,0.2444077466153785E+00,
     6  0.4100582482798061E+00,0.5408319712769052E+00,
     7  0.6309766092578797E+00,0.6839145183447968E+00,
     8  0.7122129908173587E+00,0.7283881396956757E+00,
     9  0.7394772740825155E+00,0.7485283603453863E+00,
     *  0.7567071607287845E+00,0.7644257930488514E+00,
     *  0.7718288092903780E+00,0.7789710639380831E+00,
     *  0.7858773082536592E+00,0.7925617309850295E+00,
     *  0.7990344609945628E+00,0.8053038347077053E+00,
     *  0.8113772643290438E+00,0.8172616203915311E+00,
     *  0.8229634359883509E+00,0.8284890491835000E+00,
     *  0.8338447355659797E+00,0.8390368959128114E+00,
     *  0.8440724543232754E+00,0.8489598999762000E+00,
     *  0.8537122105827861E+00,0.8583552807653678E+00,
     *  0.8629523833457415E+00,0.8676735054528577E+00,
     *  0.8729725847666965E+00,0.8799514604011378E+00,
     *  0.8909251957302293E+00,0.9094229046337805E+00,
     *  0.9365790570641110E+00,0.9660610218658774E+00,
     *  0.9881843982368550E+00,0.9983575356720777E+00,
     a  0.0d0/
c 
        data ws35/
     1  0.9792123375687179E-02,0.3589433205678722E-01,
     2  0.7202198815784074E-01,0.1111049284841143E+00,
     3  0.1475578663296468E+00,0.1770990825892073E+00,
     4  0.1965895773589052E+00,0.2039939216048679E+00,
     5  0.1983709824420185E+00,0.1798514203451931E+00,
     6  0.1496819556739958E+00,0.1108162021403531E+00,
     7  0.7001140908512356E-01,0.3817306184529694E-01,
     8  0.2059802911916622E-01,0.1287937160907447E-01,
     9  0.9771193538233473E-02,0.8507515492153301E-02,
     *  0.7911500667708642E-02,0.7547319571834854E-02,
     *  0.7267021022013562E-02,0.7021340394010978E-02,
     *  0.6793421530322587E-02,0.6577079651773303E-02,
     *  0.6369756961287036E-02,0.6170219981932728E-02,
     *  0.5977781149707290E-02,0.5792017507117309E-02,
     *  0.5612667958532916E-02,0.5439601543475065E-02,
     *  0.5272832267011505E-02,0.5112627517737269E-02,
     *  0.4959862584684232E-02,0.4817072106614292E-02,
     *  0.4691517023508293E-02,0.4604164378474769E-02,
     *  0.4615646482278865E-02,0.4896273799223824E-02,
     *  0.5875602244226301E-02,0.8459010998958622E-02,
     *  0.1413198175579811E-01,0.2318698231431191E-01,
     *  0.2997458121772731E-01,0.2720701304855322E-01,
     *  0.1622753785041729E-01,0.4772588050534671E-02,
     a  0.0d0/
c 
        data xs36/
     1  -.9967759053483217E+00,-.9758619275525324E+00,
     2  -.9242357943139140E+00,-.8357345725370844E+00,
     3  -.7103921528339198E+00,-.5530662787867442E+00,
     4  -.3722706513980436E+00,-.1791115448287242E+00,
     5  0.1373804074469098E-01,0.1930464774459577E+00,
     6  0.3461444098178006E+00,0.4628656425911428E+00,
     7  0.5397959564246526E+00,0.5844292240509146E+00,
     8  0.6101586452814443E+00,0.6271187338380999E+00,
     9  0.6404173102551021E+00,0.6521730945678712E+00,
     *  0.6631683342826318E+00,0.6736798491056675E+00,
     *  0.6838100200558603E+00,0.6936022263408166E+00,
     *  0.7030793410250225E+00,0.7122566256095817E+00,
     *  0.7211462213063824E+00,0.7297588246413433E+00,
     *  0.7381043693068954E+00,0.7461923352368447E+00,
     *  0.7540319039745378E+00,0.7616320642359251E+00,
     *  0.7690017547851395E+00,0.7761501903881030E+00,
     *  0.7830877311301752E+00,0.7898282812439281E+00,
     *  0.7963960671744383E+00,0.8028450320185729E+00,
     *  0.8093140156730898E+00,0.8161752979715289E+00,
     *  0.8243723737121491E+00,0.8360049050233451E+00,
     *  0.8548804195356612E+00,0.8848093710583356E+00,
     *  0.9234663448391436E+00,0.9607593313023970E+00,
     *  0.9867191591345246E+00,0.9981828848398375E+00,
     a  0.0d0/
c 
        data ws36/
     1  0.9433178496235343E-02,0.3465442192187947E-01,
     2  0.6960218026025025E-01,0.1073712896161231E+00,
     3  0.1424848370003338E+00,0.1707303734131353E+00,
     4  0.1889927250300206E+00,0.1951966664535122E+00,
     5  0.1882798405429850E+00,0.1681978432825040E+00,
     6  0.1362331324128612E+00,0.9653397425727088E-01,
     7  0.5867446756724163E-01,0.3302203741150333E-01,
     8  0.2013151340085500E-01,0.1459781950045566E-01,
     9  0.1232501143919730E-01,0.1130403500227755E-01,
     *  0.1072776390901031E-01,0.1031069898667050E-01,
     *  0.9956333245786275E-02,0.9631710206454947E-02,
     *  0.9324998828098497E-02,0.9031569598393789E-02,
     *  0.8749394523623725E-02,0.8477463084893378E-02,
     *  0.8215203899747957E-02,0.7962256642219083E-02,
     *  0.7718376851731237E-02,0.7483428742525018E-02,
     *  0.7257475573835915E-02,0.7041087957148326E-02,
     *  0.6836214256244549E-02,0.6648647451151428E-02,
     *  0.6495122233943463E-02,0.6423735284828474E-02,
     *  0.6570832528018372E-02,0.7299082936838585E-02,
     *  0.9431892844740115E-02,0.1448181598578451E-01,
     *  0.2403501628964841E-01,0.3548774161185994E-01,
     *  0.3996557733446365E-01,0.3284570006021933E-01,
     *  0.1852651834832715E-01,0.5298976186552273E-02,
     a  0.0d0/
c 
        data xs37/
     1  -.9971843318345427E+00,-.9785667079620213E+00,
     2  -.9318223164751303E+00,-.8507325182584013E+00,
     3  -.7349763385052072E+00,-.5889727800251370E+00,
     4  -.4208339975053584E+00,-.2414099750497374E+00,
     5  -.6335507314425270E-01,0.9989284061643566E-01,
     6  0.2356813497505110E+00,0.3354960637274351E+00,
     7  0.4000446539258955E+00,0.4394188166114923E+00,
     8  0.4651235404834384E+00,0.4845395083744559E+00,
     9  0.5012183492437925E+00,0.5166054001119684E+00,
     *  0.5312397980355608E+00,0.5453226735424592E+00,
     *  0.5589357650160105E+00,0.5721184491494773E+00,
     *  0.5848943795135231E+00,0.5972809153946876E+00,
     *  0.6092926479781313E+00,0.6209428257970826E+00,
     *  0.6322439722514468E+00,0.6432081547319601E+00,
     *  0.6538471014026686E+00,0.6641722832856813E+00,
     *  0.6741950786397066E+00,0.6839271945347740E+00,
     *  0.6933817875191011E+00,0.7025764859612280E+00,
     *  0.7115415607480702E+00,0.7203422902890733E+00,
     *  0.7291401390433296E+00,0.7383521649915428E+00,
     *  0.7490061791151068E+00,0.7633328490370239E+00,
     *  0.7853248025582564E+00,0.8195632310328624E+00,
     *  0.8659800688107202E+00,0.9163014309698652E+00,
     *  0.9590798922178903E+00,0.9866248547881935E+00,
     *  0.9982152628490837E+00/
c 
        data ws37/
     1  0.8284709381022180E-02,0.3109889248576235E-01,
     2  0.6341116831944685E-01,0.9880172263535576E-01,
     3  0.1319569063979251E+00,0.1586824149588372E+00,
     4  0.1757714091886660E+00,0.1809513737163130E+00,
     5  0.1728957988988310E+00,0.1514392024424430E+00,
     6  0.1186153928905716E+00,0.8110830918203328E-01,
     7  0.4984947407962845E-01,0.3091106974475082E-01,
     8  0.2170844730178230E-01,0.1768354970112336E-01,
     9  0.1589519741757605E-01,0.1496013197393411E-01,
     *  0.1433902291534883E-01,0.1383925588543561E-01,
     *  0.1339311168284412E-01,0.1297606359366572E-01,
     *  0.1257863881092259E-01,0.1219684335096551E-01,
     *  0.1182882030216807E-01,0.1147361826821408E-01,
     *  0.1113068397428138E-01,0.1079963439379151E-01,
     *  0.1048016821187948E-01,0.1017208562421110E-01,
     *  0.9875439919372474E-02,0.9590940617299282E-02,
     *  0.9321065257039072E-02,0.9073054506732325E-02,
     *  0.8867086656245129E-02,0.8758894812884352E-02,
     *  0.8900571245168009E-02,0.9684560478674581E-02,
     *  0.1198304611643759E-01,0.1734387532912172E-01,
     *  0.2750548480196370E-01,0.4105043846420040E-01,
     *  0.5028955363983602E-01,0.4831684914242589E-01,
     *  0.3591081435159934E-01,0.1905533076892946E-01,
     *  0.5239867842556202E-02/
c 
        data xs38/
     1  -.9971305181790984E+00,-.9783921344784753E+00,
     2  -.9319245926967599E+00,-.8521986566622971E+00,
     3  -.7396322913154771E+00,-.5994475346168715E+00,
     4  -.4407176847427319E+00,-.2755472119990952E+00,
     5  -.1181917798523121E+00,0.1660096154628321E-01,
     6  0.1178996760259059E+00,0.1853194886561091E+00,
     7  0.2289633118255210E+00,0.2601053498188963E+00,
     8  0.2856464288134568E+00,0.3086881433643791E+00,
     9  0.3304081726213586E+00,0.3512317171357235E+00,
     *  0.3713208337956908E+00,0.3907479204134978E+00,
     *  0.4095539023547934E+00,0.4277675669061924E+00,
     *  0.4454122303139385E+00,0.4625083460563693E+00,
     *  0.4790748069841002E+00,0.4951297479765810E+00,
     *  0.5106910124437514E+00,0.5257763790367564E+00,
     *  0.5404037772436654E+00,0.5545919039506362E+00,
     *  0.5683621427723680E+00,0.5817441277884678E+00,
     *  0.5947914879523348E+00,0.6076264951180168E+00,
     *  0.6205651136806916E+00,0.6344418023510843E+00,
     *  0.6512887839992918E+00,0.6752873561896153E+00,
     *  0.7128118421726632E+00,0.7678567006432964E+00,
     *  0.8347550750247650E+00,0.9003389819137797E+00,
     *  0.9524266610686798E+00,0.9846650667701635E+00,
     *  0.9979654340891015E+00,0.0d0,0.0d0/
c 
        data ws38/
     1  0.8413121127381169E-02,0.3112776594805041E-01,
     2  0.6271131201148959E-01,0.9665117745811480E-01,
     3  0.1275955832149486E+00,0.1512483557513484E+00,
     4  0.1641805329679595E+00,0.1637616672771547E+00,
     5  0.1484248297279825E+00,0.1191980601770978E+00,
     6  0.8331693381251276E-01,0.5341870092991467E-01,
     7  0.3583543435216240E-01,0.2757494177651255E-01,
     8  0.2398905816900164E-01,0.2227193576433215E-01,
     9  0.2123181517618692E-01,0.2043982186103727E-01,
     *  0.1974966080451101E-01,0.1911106235063615E-01,
     *  0.1850560492616176E-01,0.1792557206395804E-01,
     *  0.1736714353601946E-01,0.1682822360117308E-01,
     *  0.1630771324484212E-01,0.1580513954686974E-01,
     *  0.1532034996240690E-01,0.1485336319745225E-01,
     *  0.1440452229349093E-01,0.1397521539650413E-01,
     *  0.1356997637147541E-01,0.1320229363808050E-01,
     *  0.1291099622474517E-01,0.1280675250142889E-01,
     *  0.1319822779422890E-01,0.1487635920341424E-01,
     *  0.1950325854569649E-01,0.2962604307738218E-01,
     *  0.4626379890493871E-01,0.6282586656737886E-01,
     *  0.6860275649158744E-01,0.6046198716369746E-01,
     *  0.4263577947644679E-01,0.2199295241327081E-01,
     *  0.5978316599912607E-02,0.0d0,0.0d0/
c 
        data xs39/
     1  -.9971719347823668E+00,-.9789739465013674E+00,
     2  -.9345377733647775E+00,-.8594701956557286E+00,
     3  -.7554158630894393E+00,-.6291284125975649E+00,
     4  -.4917606606130811E+00,-.3578740431758339E+00,
     5  -.2426813373618336E+00,-.1556051191238306E+00,
     6  -.9441210466719352E-01,-.4949052846522584E-01,
     7  -.1250741667718502E-01,0.2086069235819865E-01,
     8  0.5232667368861987E-01,0.8250728505545246E-01,
     9  0.1116321842934440E+00,0.1398000070442642E+00,
     *  0.1670640210734491E+00,0.1934611084427171E+00,
     *  0.2190222366629226E+00,0.2437762289454808E+00,
     *  0.2677506968651016E+00,0.2909719824508459E+00,
     *  0.3134652533615201E+00,0.3352552304573040E+00,
     *  0.3563675111551148E+00,0.3768306465564379E+00,
     *  0.3966805789998391E+00,0.4159730008942373E+00,
     *  0.4348200612034121E+00,0.4534984032881146E+00,
     *  0.4727525173522811E+00,0.4945262075186994E+00,
     *  0.5231869428350246E+00,0.5662186358312506E+00,
     *  0.6306743882432309E+00,0.7143274500433725E+00,
     *  0.8042224427139983E+00,0.8851672291428037E+00,
     *  0.9462101146844832E+00,0.9828563863345936E+00,
     *  0.9977351654826454E+00,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data ws39/
     1  0.8257165586526029E-02,0.3002843966604137E-01,
     2  0.5955867522031813E-01,0.9028661108467821E-01,
     3  0.1166856374740235E+00,0.1340171120854554E+00,
     4  0.1382528380017460E+00,0.1268719163051007E+00,
     5  0.1018043774733280E+00,0.7281570497482934E-01,
     6  0.5135016021692557E-01,0.3991491315161756E-01,
     7  0.3473756256212573E-01,0.3225817265651511E-01,
     8  0.3076600899612992E-01,0.2962985471477475E-01,
     9  0.2863518450642250E-01,0.2770893007546050E-01,
     *  0.2682509408285224E-01,0.2597418848511909E-01,
     *  0.2515286276257006E-01,0.2435972964826949E-01,
     *  0.2359359792928386E-01,0.2285315669081153E-01,
     *  0.2213748395899635E-01,0.2144670529119090E-01,
     *  0.2078255172865811E-01,0.2014953605850141E-01,
     *  0.1955909293248033E-01,0.1904290545737768E-01,
     *  0.1869288287590050E-01,0.1877624562716304E-01,
     *  0.2003370841431706E-01,0.2422182252074601E-01,
     *  0.3439045603845218E-01,0.5304496708103107E-01,
     *  0.7552755793237381E-01,0.8943141920622641E-01,
     *  0.8772108253550891E-01,0.7231473798697845E-01,
     *  0.4901199947314451E-01,0.2471940514329627E-01,
     *  0.6657535140990148E-02,0.0d0,0.0d0,0.0d0,0.0d0/
c 
  
  
        data xs40/
     1  -.9978244420523276E+00,-.9835221754200973E+00,
     2  -.9482952548953438E+00,-.8895874057257239E+00,
     3  -.8114661984829510E+00,-.7235321523308536E+00,
     4  -.6380859923650637E+00,-.5644673386277883E+00,
     5  -.5042405250588234E+00,-.4533301788287747E+00,
     6  -.4074435340280291E+00,-.3642103095882975E+00,
     7  -.3226418180900480E+00,-.2823658446745625E+00,
     8  -.2432422313254129E+00,-.2052123239547679E+00,
     9  -.1682452956641383E+00,-.1323197051623446E+00,
     *  -.9741649371817862E-01,-.6351602732433059E-01,
     *  -.3059754137466630E-01,0.1360363648836849E-02,
     *  0.3237979426699078E-01,0.6248462280579159E-01,
     *  0.9170184538868731E-01,0.1200633865797773E+00,
     *  0.1476113355351869E+00,0.1744147700762144E+00,
     *  0.2006192627539403E+00,0.2265878263903367E+00,
     *  0.2532856309225056E+00,0.2831783157743941E+00,
     *  0.3216304806362389E+00,0.3773678864101457E+00,
     *  0.4581491416992061E+00,0.5619170197180230E+00,
     *  0.6759751941276166E+00,0.7849312443170648E+00,
     *  0.8763634077983538E+00,0.9426847399494902E+00,
     *  0.9817788017052542E+00,0.9975857023671646E+00,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
  
c 
        data ws40/
     1  0.6397179033967202E-02,0.2376331155442371E-01,
     2  0.4711413880058928E-01,0.6960724482080320E-01,
     3  0.8499316712650898E-01,0.8872285642496532E-01,
     4  0.8050767880357926E-01,0.6651943541089826E-01,
     5  0.5472515461454806E-01,0.4784046351614240E-01,
     6  0.4431190622155098E-01,0.4230604656567118E-01,
     7  0.4088712379406020E-01,0.3968599025438293E-01,
     8  0.3857024954933507E-01,0.3749444689568542E-01,
     9  0.3644308062172245E-01,0.3541122810565634E-01,
     *  0.3439845936237367E-01,0.3340594536169465E-01,
     *  0.3243459145649729E-01,0.3148488227089260E-01,
     *  0.3055793518762882E-01,0.2965619908271244E-01,
     *  0.2878349847262395E-01,0.2794636198863232E-01,
     *  0.2716034303618445E-01,0.2646876810384896E-01,
     *  0.2599504707358951E-01,0.2608848673622935E-01,
     *  0.2769184910040160E-01,0.3297025325406417E-01,
     *  0.4544947814744619E-01,0.6751621672029912E-01,
     *  0.9371601836500958E-01,0.1115308228909753E+00,
     *  0.1139292597238274E+00,0.1018889484241901E+00,
     *  0.7970266354509672E-01,0.5257712779248135E-01,
     *  0.2625943291069229E-01,0.7086701701652775E-02,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data xs41/
     1  -.9984889570890517E+00,-.9891474146800662E+00,
     2  -.9679563824614772E+00,-.9354554196311098E+00,
     3  -.8946827024018078E+00,-.8489277890465374E+00,
     4  -.8006526690119833E+00,-.7513810506001187E+00,
     5  -.7019758326354778E+00,-.6529191852069010E+00,
     6  -.6044902029553474E+00,-.5568594542001148E+00,
     7  -.5101362487195589E+00,-.4643924793162697E+00,
     8  -.4196753700687634E+00,-.3760143977529493E+00,
     9  -.3334252487972508E+00,-.2919121516894750E+00,
     *  -.2514672642989827E+00,-.2120644406558948E+00,
     *  -.1736403109995601E+00,-.1360373948400545E+00,
     *  -.9883517937963016E-01,-.6087804368179051E-01,
     *  -.1915334699161186E-01,0.3287349815451845E-01,
     *  0.1052737038908849E+00,0.2062611072837826E+00,
     *  0.3340084617288603E+00,0.4768534193727703E+00,
     *  0.6200020684965595E+00,0.7500137667878022E+00,
     *  0.8565429474036721E+00,0.9332274984549707E+00,
     *  0.9785773749467332E+00,0.9971246091715553E+00,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     b  0.0d0,0.0d0,0.0d0/
  
  
        data ws41/
     1  0.4376265245379458E-02,0.1501296967541977E-01,
     2  0.2723450295699060E-01,0.3722119187038531E-01,
     3  0.4375887545570258E-01,0.4734040616276701E-01,
     4  0.4896142360838269E-01,0.4944283231745544E-01,
     5  0.4929011999263463E-01,0.4877813931405894E-01,
     6  0.4805212323437792E-01,0.4719152376310962E-01,
     7  0.4624311228079219E-01,0.4523667082964065E-01,
     8  0.4419269866909216E-01,0.4312667251184408E-01,
     9  0.4205091147435884E-01,0.4097649287948939E-01,
     *  0.3991737222096629E-01,0.3889803884916698E-01,
     *  0.3797405957889553E-01,0.3729493330097350E-01,
     *  0.3728409272785383E-01,0.3909940147312107E-01,
     *  0.4542748270784752E-01,0.6039327261549629E-01,
     *  0.8595931253183925E-01,0.1156937920181903E+00,
     *  0.1377171398610726E+00,0.1454366708455651E+00,
     *  0.1385965989938824E+00,0.1196893378428169E+00,
     *  0.9229445270090091E-01,0.6080670749683954E-01,
     *  0.3062795243757721E-01,0.8402436227337867E-02,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     b  0.0d0,0.0d0,0.0d0/
c 
        data xs42/
     1  -.9985556909909618E+00,-.9895284360387931E+00,
     2  -.9688727557483384E+00,-.9370447156681621E+00,
     3  -.8970238376973548E+00,-.8520483201605287E+00,
     4  -.8045344079768811E+00,-.7559697834585407E+00,
     5  -.7071794075344974E+00,-.6585784523448366E+00,
     6  -.6102904511144012E+00,-.5620884394523834E+00,
     7  -.5130428094598277E+00,-.4606321091514481E+00,
     8  -.3993653865139211E+00,-.3204441296507803E+00,
     9  -.2148957715288459E+00,-.7944536339151720E-01,
     *  0.8083837278352091E-01,0.2552850771559281E+00,
     *  0.4312491271990662E+00,0.5963530929789125E+00,
     *  0.7398389395173328E+00,0.8536177350086949E+00,
     *  0.9332649238515260E+00,0.9790850210926054E+00,
     *  0.9972643737740410E+00,0.0d0,0.0d0,0.0d0,0.0d0,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     b  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
c 
        data ws42/
     1  0.4196947679652322E-02,0.1457255404392639E-01,
     2  0.2661887382766595E-01,0.3650018913878102E-01,
     3  0.4298483001495264E-01,0.4656341141722103E-01,
     4  0.4822179381394116E-01,0.4877437215211892E-01,
     5  0.4873919440383277E-01,0.4844138570508708E-01,
     6  0.4816723621210245E-01,0.4837744561817908E-01,
     7  0.5010757331863923E-01,0.5561557702074817E-01,
     8  0.6843588849735691E-01,0.9101983337658533E-01,
     9  0.1206741503155796E+00,0.1493186392796292E+00,
     *  0.1694060708021879E+00,0.1773446963740976E+00,
     *  0.1724993568791149E+00,0.1558985413227069E+00,
     *  0.1297086739161130E+00,0.9709397137804725E-01,
     *  0.6223090657459907E-01,0.3042633072530855E-01,
     *  0.8061551988357710E-02,0.0d0,0.0d0,0.0d0,0.0d0,
     a  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     b  0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
c 
c 
c 
c 
c 
  
  
  
  
c 
c       return to the user the nodes and weights of all 42 quadratures
c 
        do 1400 i=1,42
        do 1200 j=1,47
c 
        xsout(j,i)=xz(j,i)
        wsout(j,i)=wz(j,i)
 1200 continue
 1400 continue
  
        return
        end
