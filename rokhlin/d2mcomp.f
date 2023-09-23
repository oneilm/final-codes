        implicit real *8 (a-h,o-z)
        real *8 udad(1000 000),vdad(1000 000),
     1    w(4 000 000)
        external funmat
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
         PRINT *, 'ENTER nrowmax'
         READ *,nrowmax
         CALL PRINF('nrowmax=*',nrowmax,1 )
c 
         PRINT *, 'ENTER itypein'
         READ *,itypein
         CALL PRINF('itypein=*',itypein,1 )
c 
c        initialize the matrix computer
c 
        call funmatin(n,m)
c 
c       compress the matrix recursively
c 
        eps=1.0d-12
c 
        nsmall=nrowmax
        maxcols=40
        ltot=2 000 000
ccc        ltot=60000
cccc        itypein=1
        call d2mcomp(ier,funmat,n,m,eps,nsmall,
     1    maxcols,itypein,udad,vdad,ncols,itype,w,ltot,lused,
     2    w1,w2,w3,w4,w5)
c 
          call prinf('after d2mcomp, ier=*',ier,1)
          call prinf('after d2mcomp, lused=*',lused,1)
          call prinf('after d2mcomp, ncols=*',ncols,1)
ccc          call prin2('after d2mcomp, udad=*',udad,ncols*n,1)
ccc          call prin2('after d2mcomp, vdad=*',vdad,ncols*m,1)
c 
c       check the compression of the matrix
c 
        call exptes(udad,vdad,n,m,ncols,itype,funmat,
     1      'checking compressed expansion after d2comp0*')
  
        stop
  
        end
c 
c 
c 
c 
c 
        subroutine actplot(iw,n,m,boxes,nboxes,
     1      actions,nact,mes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(12,1),actions(4,1)
        real *8 xx(5),yy(5)
        character *1 mes(1)
c 
c        plot the title
c 
        IF(IW.EQ.0) RETURN
        call titlpr(MES,Iw)
C 
C 
C . . . FIND THE LIMITS TO BE SET
C 
        dMIN=0.5
        dMAX=n+0.5
c 
 1200 FORMAT(1X,37H#SET DEVICE postscript FILE 'plot.ps',/,
cccc     2  1X,'SET AXIS TOP OFF',/,
cccc     3  1X,'SET AXIS RIGHT OFF',/,
     4  1X,'SET LIMITS X ', 2(1X,E11.5),' Y ',2(1X,E11.5),/,
     5  1X,'SET SCALE EQUAL')
c 
          WRITE(IW,1200) DMIN,DMAX,DMIN,DMAX
 1400 FORMAT(2(2X,E15.7))
c 
c       plot all boxes
c 
        do 3000 i=1,nboxes
        yy(1)=boxes(2,i)-0.5
        yy(2)=yy(1)
        yy(3)=boxes(3,i)+0.5
        yy(4)=yy(3)
        yy(5)=yy(1)
c 
        xx(1)=boxes(4,i)-0.5
        xx(2)=boxes(5,i)+0.5
        xx(3)=xx(2)
        xx(4)=xx(1)
        xx(5)=xx(1)
c 
 1410 FORMAT(2(2X,E15.7))
        write(iw,1410) (xx(j),yy(j),j=1,5)
        write(iw,1450)
 3000 continue
c 
c       now, plot the actions
c 
        do 5000 i=1,nact
c 
c       if this is a childless box being compressed
c       - plot this fact
c 
        if(actions(1,i) .ne. 1) goto 4000
c 
        k=actions(2,i)
        xxx=boxes(4,k)+boxes(5,k)
        xxx=xxx/2
        yyy=boxes(2,k)+boxes(3,k)
        yyy=yyy/2
c 
        write(iw,1410) xxx,yyy
 3200 format(1x,'plot ''*'' ')
        write(iw,3200)
        goto 5000
 4000 continue
c 
c      this is two boxes being combined. plot this fact
c 
        if(actions(1,i) .ne. 2) goto 5000
        k=actions(2,i)
c 
        y1=boxes(2,k)
        y2=boxes(3,k)
        x1=boxes(4,k)
        x2=boxes(5,k)
c 
        write(iw,1410) x1,y1
        write(iw,1410) x2,y2
        write(iw,1450)
c 
        write(iw,1410) x1,y2
        write(iw,1410) x2,y1
        write(iw,1450)
 5000 continue
c 
c 
c 
C1600 FORMAT(1X,'FRAME')
C1800 FORMAT(1X,'EXIT')
 1450 FORMAT(1X,'join ')
c1500 FORMAT(1X,'PLOT')
 1600 FORMAT(1X,'#FRAME')
 1800 FORMAT(1X,'#EXIT')
cccc        WRITE(IW,1500)
        WRITE(IW,1600)
        WRITE(IW,1800)
        return
        end
c 
c 
c 
c 
c 
        subroutine allplt(iw,n,m,boxes,nboxes,mes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(12,1)
ccc        real *8 z(2,1),xx(5),yy(5),center0(2),center(2)
        real *8 xx(5),yy(5)
        character *1 mes(1)
c 
c        plot the title
c 
        IF(IW.EQ.0) RETURN
        call titlpr(MES,Iw)
C 
C 
C . . . FIND THE LIMITS TO BE SET
C 
        dMIN=0.5
        dMAX=n+0.5
c 
 1200 FORMAT(1X,37H#SET DEVICE postscript FILE 'plot.ps',/,
cccc     2  1X,'SET AXIS TOP OFF',/,
cccc     3  1X,'SET AXIS RIGHT OFF',/,
     4  1X,'SET LIMITS X ', 2(1X,E11.5),' Y ',2(1X,E11.5),/,
     5  1X,'SET SCALE EQUAL')
c 
          WRITE(IW,1200) DMIN,DMAX,DMIN,DMAX
 1400 FORMAT(2(2X,E15.7))
c 
c       plot all boxes
c 
        do 3000 i=1,nboxes
        yy(1)=boxes(2,i)-0.5
        yy(2)=yy(1)
        yy(3)=boxes(3,i)+0.5
        yy(4)=yy(3)
        yy(5)=yy(1)
c 
        xx(1)=boxes(4,i)-0.5
        xx(2)=boxes(5,i)+0.5
        xx(3)=xx(2)
        xx(4)=xx(1)
        xx(5)=xx(1)
c 
 1410 FORMAT(2(2X,E15.7))
        write(iw,1410) (xx(j),yy(j),j=1,5)
        write(iw,1450)
 3000 continue
C1600 FORMAT(1X,'FRAME')
C1800 FORMAT(1X,'EXIT')
 1450 FORMAT(1X,'join ')
c1500 FORMAT(1X,'PLOT')
 1600 FORMAT(1X,'#FRAME')
 1800 FORMAT(1X,'#EXIT')
cccc        WRITE(IW,1500)
        WRITE(IW,1600)
        WRITE(IW,1800)
        return
        end
c 
c 
c 
c 
c 
        subroutine orttes(u,n,m,w,mes)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,m),w(m,m)
        character *1 mes(1)
c 
c       calculate the inner products between all columns of u
c 
        error=0
        do 1400 i=1,m
        do 1200 j=1,m
        call scapro(u(1,i),u(1,j),n,dd )
        if(i .eq. j) error=error+(dd-1)**2
        if(i .ne. j) error=error+dd**2
 1200 continue
 1400 continue
        call prin2(mes,dsqrt(error),1)
        return
        end
c 
c 
c 
c 
c 
        subroutine boxesprt(boxes,nboxes)
        save
        integer *4 boxes(12,1)
c 
        call prinf('printing boxes; of them there are*',nboxes,1)
 1200 format(12(1x,i5))
        do 2000 i=1,nboxes
        write(6,1200)(boxes(j,i),j=1,12)
        write(13,1200)(boxes(j,i),j=1,12)
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine submcomp(a,j1,j2,i1,i2)
        implicit real *8 (a-h,o-z)
        save
        dimension a(j2-j1+1,1)
c 
        do 1400 i=1,i2-i1+1
        do 1200 j=1,j2-j1+1
        call funmat(j+j1-1,i+i1-1,a(j,i),w1,w2,w3,w4,w5)
 1200 continue
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine mulmata0(b,v,n,m,ncols,i,j,cij)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,ncols),v(m,ncols)
c 
c       this subroutine multiplies the matrix b by the adjoint of v,
c       obtaining the matrix c
c 
c                  input parameters:
c 
c  b,v - the matrices to be multiplied
c  n,m,ncol - the dimensionalities of the matrices to be multiplied,
c       and of the result (see array declarations above)
c 
c                  output parameters:
c 
c  c - the result of multiplication
c 
c        multiply the matrices a and b ^*  obtaining the matrix c
c 
        cij=0
        do 2200 k=1,ncols
        cij=cij+b(i,k)*v(j,k)
 2200 continue
cccc         call prin2('in mulmata0, cij=*',cij,1)
        return
        end
  
c 
c 
c 
c 
c 
        subroutine exptes(u,v,n,m,ncols,itype,funmat,mes)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,1),v(m,1)
        character *1 mes(1)
c 
c       print the message
c 
        call prin2('   *',a,0)
        call prin2('   *',a,0)
        call prin2('   *',a,0)
        call prin2(mes,u,0)
c 
c        depending on the type of the expansion, check the
c        orthogonality of either u or v
c 
        if(itype .ne. 1) goto 1200
        call orttes(v,m,ncols,w,
     1      'since this matrix will be compressed horizontally, checking
     2 orthogonality of v*')
 1200 continue
c 
        if(itype .ne. 2) goto 1400
        call orttes(u,n,ncols,w,
     1      'since this matrix will be compressed vertically, checking
     2orthogonality of u*')
 1400 continue
c 
c       now, check the accuracy of the expansion
c 
ccccc        call mulmatad(u,v,n,m,ncols,w)
c 
        d=0
        dd=0
        do 3800 i=1,m
        do 3600 j=1,n
        call mulmata0(u,v,n,m,ncols,j,i,cji)
        call funmat(j,i,aji,w1,w2,w3,w4,w5)
ccc        d=d+(a(j,i)-w(j,i))**2
cccc        d=d+(a(j,i)-cji)**2
        d=d+(aji-cji)**2
cccc        dd=dd+a(j,i)**2
        dd=dd+aji**2
 3600 continue
 3800 continue
        d=dsqrt(d)
        dd=dsqrt(dd)
cccc         call prin2('in exptes, a=*',a,n*m)
cccc         call prin2('in exptes, w=*',w,n*m)
c 
        call prin2('cumulative error of decomposition is*',
     1      d,1)
        call prin2('and relative error is*',d/dd,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine funmatin(n,m)
        implicit real *8 (a-h,o-z)
        save
        dimension x(10000),y(10000)
c 
        done=1
        hm=done/m
        hn=done/n
c 
        do 1200 i=1,n
        x(i)=i*hn
 1200 continue
c 
        do 1400 i=1,m
        y(i)=i*hm+2
 1400 continue
        return
c 
c 
c 
c 
        entry funmat(jj,ii,aji,w1,w2,w3,w4,w5)
c 
        d=x(jj)-y(ii)
        d=d**2
cccc        aji=dlog(d)*dsin((x(jj)-y(ii))*10)
        aji=dlog(d)
        return
        end
  
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c           this is the end of the debugging code and the beginning
c           of the actual matrix-handling subroutines
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine d2mcomp(ier,funmat,n,m,eps,nsmall,
     1    maxcols,itypein,u,v,ncols,itypeout,w,ltot,lused,
     2    w1,w2,w3,w4,w5)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,1),v(m,1),w(1),
     1      w1(1),w2(1),w3(1),w4(1),w5(1)
c 
c 
c       this subroutine produces the compressed form of a
c       user-supplied matrix amatr. the user supplies the matrix
c       in the form of the subroutine funmat, evaluating
c       the matrix elements as a function of the indices.
c       (see the description of the subroutine funmat below).
c       the output of the subroutine is in the form of two
c       matrices u(n,ncols), v(n,ncols), such that
c 
c       amatr=u \times v^T,
c 
c       and ncols is the (numerical) rank of the user-supplied
c       matrix a.
c 
c                 input parameters:
c 
c  funmat - the user-provided subroutinbe for the evaluation of the
c        elements of the matrix to be compressed. its calling sequence
c        must be:
c 
c         funmat(j,i,aji,w1,w2,w3,w4,w5)
c 
c        with the following input parameters:
c 
c        (j,i) - the indices of the element to be crated by this
c                call
c        w1,w2,w3,w4,w5 - whatever parameters the user's subroutine funmat
c                needs (real, integer, array, variable - whatever),
c 
c        and the following output parameters:
c 
c        aji - the element of the matrix with the indices (j,i)
c 
c  n - the number of rows in the matrix whose compressed form
c        is to be constructed
c  m - the number of columns in the matrix whose compressed form
c        is to be constructed
c  eps - the precision of computations
c  nsmall - the maximum size of the box on the finest level
c        of the tree structure into which the matrix will be
c        subdivided. recommended value: 200
c  maxcols - the maximum permitted number of columns in the matrices
c        u,v. note that if for the specified precision eps, the
c        matrices u,v need more columns than maxcols, the error code
c        ier=1024 will be generated, and the execution of the
c        subroutine terminated.
c  itypein - the type of the compressed expansion to be produced
c        by the subroutine. explanation:
c        itypein=1 means that the columns of v will be orthonormal
c        itypein=2 means that the columns of u will be orthonormal
c        itypein=0 will permit the subroutine to produce the
c           expansion of the type it deems appropriate, informing
c           the user of its decision via the parameter itypeout (see
c           below).
c  ltot - the amount of space provided by the user in the array w,
c        in terms of real *8 words.
c        used by the subroutine to avoid exceeding the array w. if
c        ltot is too small for the atsk at hand, the error return code
c        will be set ier=256, and the execution of the subroutine
c        terminated.
c  w1,w2,w3,w4,w5 - whatever parameters the user's subroutine
c        funmat needs (real, integer, array, variable - whatever).
c 
c                 output parameters:
c 
c  ier - error return code.
c         ier=0 means successful conclusion
c         ier=256 means that the amount of storage provided in array
c                   arr has been insufficient. this is a fatal error.
c         ier=1024 means that the permitted number maxcol of columns in
c                    the matrices udad, vdad has been exceeded.
c                    this is a fatal error.
c         ier=32000 means that something has gone terribly wrong,
c                    so much so that the subroutine can not come up
c                    with a diagnosis. possibly, the array w is shorter
c                    than ltot specified by the user, so that the
c                    subroutine has exceeded the array w without knowing
c                    it. possibly, the user has managed to hit an
c                    internal bug in the subroutine. in any case, this
c                    is bigtime trouble.
c  u - the matrix u in the compressed form of the user-supplied
c         matrix; must be dimensioned at least max(n*maxcols,nsmall**2)
c  v - the matrix v in the compressed form of the user-supplied
c         matrix; must be dimensioned at least max(m*maxcols,nsmall**2)
c  ncols - the number of columns in the arrays u, v
c  itypeout - the type of the decomposition u * v^* of the user's
c         matrix. explanation:
c              itypeout=1 means that the columns of the matrix v are
c                          orthonormal
c              itypeout=2 means that the columns of the matrix u are
c                          orthonormal
c  lused - the amount of memory in array w actually used by this call
c           to the subroutine
c 
c                         work arrays:
c 
c  w - must be long; at least ????????
c 
        ninire=2
        ier=0
c 
c       check if the amount of memory provided by the user
c       is obviously insufficient
c 
        ll=nsmall**2+50+maxcols*n+maxcols*m
        if(ltot .gt. ll) goto 1200
        ier=256
        return
 1200 continue
c 
c 
c       construct the quad-tree structure
c 
        iboxes=1
        call boxesbld(n,m,nsmall,w(iboxes),nboxes)
        call boxesprt(w(iboxes),nboxes)
c 
c       allocate memory for the action tree,
c       and construct the action tree
c 
        lboxes=(nboxes*12+12)/ninire
        iactions=iboxes+lboxes
c 
        call acttree(w(iboxes),nboxes,w(iactions),nact)
c 
c        allocate memory for the main compression routine
c 
        lactions=(nact*4+4)/ninire
c 
        imap=iactions+lactions
        lmap=100*10/ninire
c 
        irnorms=imap+lmap
        lrnorms=m
        if(n .gt. m) lrnorms=n
        lrnorms=lrnorms+10
c 
        iw=irnorms+lrnorms
c 
        nm=n
        if(m .gt. nm) nm=m
        lw=2*nm*(maxcols+1)+16*maxcols**2+2*maxcols+50
        lw7=nsmall**2+100
        if(lw7 .gt. lw) lw=lw7
c 
ccccc           lw=lw*2
c 
        iarr=iw+lw
        lleft=ltot-iarr
         call prinf('in d2mcomp, lleft=*',lleft,1)
c 
        if(lleft .gt. 10000) goto 1400
        ier=256
        return
 1400 continue
c 
c        compress the user-supplied matrix
c 
        call d2mcomp0(jer,w(iboxes),w(iactions),nact,
     1    funmat,eps,w1,w2,w3,w4,w5,maxcols,lleft,
     2    w(iarr),w(imap),u,v,w(irnorms),
     3    w(iw),ncols,itypedad,larrmax)
           call prinf('after d2mcomp0, larrmax=*',larrmax,1)
          lused=iarr+larrmax+10
c 
c        if the compression failed for whatever reason - bomb out
c 
        if(jer .eq. 0) goto 1600
        ier=jer
        return
 1600 continue
c 
c       if the user has requested a particular type of expansion
c       - enforce this type
c 
        itypeout=itypedad
c 
        if(itypein .eq. 0) return
        if(itypein .eq. itypeout) return
c 
        iu=1
        lu=ncols*n+10
c 
        iv=iu+lu
        lv=ncols*m+10
c 
        iw=iv+lv
c 
        if(itypeout .eq. 1)
     1  call prinf('changing orthogonality, itypeout=*',itypeout,1)
        if(itypeout .eq. 2)
     1  call prinf('changing orthogonality, itypeout=*',itypeout,1)
  
        if(itypein .eq. 1)
     1      call makevort(u,v,n,m,ncols,
     2      w(iu),w(iv),w(iw))
  
        if(itypein .eq. 2)
     1      call makeuort(u,v,n,m,ncols,
     2      w(iu),w(iv),w(iw))
c 
         call arrcopy(w(iu),u,n*ncols)
         call arrcopy(w(iv),v,m*ncols)
         itypeout=itypein
  
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mcomp0(ier,boxes,actions,nact,
     1    funmat,eps,w1,w2,w3,w4,w5,maxcols,lenarr,
     2    arr,map,udad,vdad,rnorms,
     3    w,ncolsdad,itypedad,larrmax)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(12,1),actions(4,1),map(10,1)
        real *8 udad(1),w(1),arr(1),vdad(1)
c 
c       this subroutine performs the actual compression
c       of the user-specified matrix. it assumes that both
c       the box structure and the action tree have already
c       been constructed.
c 
c 
c                 input parameters:
c 
c  boxes - integer array dimensioned boxes(12,nboxes) describing
c        the tree structure (the subdivision of the matrix); normally,
c        it is constructed by a call to the subroutine boxesbld (see)
c  actions - integer *4 array dimensioned (4,nact). in fact, this is an
c        action table for the  actual compression of a matrix, with each
c        column describing one action normally, it is constructed
c        by a call to the subroutine acttree (see).
c  nact - the number of columns in array actions
c  funmat - the user-provided subroutinbe for the evaluation of the
c        elements of the matrix to be compressed. its calling sequence
c        must be:
c 
c         funmat(j,i,aji,w1,w2,w3,w4,w5)
c 
c        with the following input parameters:
c 
c        (j,i) - the indices of the element to be crated by this
c                call
c        w1,w2,w3,w4,w5 - whatever parameters the user's subroutine funmat
c                needs (real, integer, array, variable - whatever),
c 
c        and the following output parameters:
c 
c        aji - the element of the matrix with the indices (j,i)
c 
c  eps - the precision of computations
c 
c  w1,w2,w3,w4,w5 - whatever parameters the user's subroutine
c        funmat needs (real, integer, array, variable - whatever).
c 
c                 output parameters:
c 
c  ier - error return code.
c         ier=0 means successful conclusion
c         ier=256 means that the amount of storage provided in array
c                   arr has been insufficient. this is a fatal error.
c         ier=1024 means that the permitted number maxcol of columns in
c                    the matrices udad, vdad has been exceeded.
c                    this is a fatal error.
c  udad - the matrix u in the compressed form of the user-supplied
c         matrix
c  vdad - the matrix v in the compressed form of the user-supplied
c         matrix
c  ncolsdad - the number of columns in the arrays udad, vdad
c  itypedad - the type of the decomposition u * v^* of the user's
c         matrix. explanation:
c              itypedad=1 means that the columns of the matrix u are
c                          orthonormal
c              itypedad=2 means that the columns of the matrix v are
c                          orthonormal
c  larrmax - the amount of memory in array arr actually used by this call
c           to the subroutine
c 
c                         work arrays:
c 
c  arr - should be long, since it is the main storage array for
c          the intermediate compressions created by this subroutine
c  map - must be at least  100*10/ninire real *8 elements long
c  rnorms - must be at least max(m,n) +1 real *8 elements long
c  w - must be at least
c 
c     max(2*nm*(maxcols+1)+16*maxcols**2+2*maxcols+50, nsmall**2+100)
c 
c          real *8 elements long
c 
c 
c         . . . execute the instructions in the array actions
c               one after another
c 
          call prin2('in d2comp0, eps=*',eps,1)
          call prinf('in d2comp0, maxcols=*',maxcols,1)
          call prinf('in d2comp0, lenarr=*',lenarr,1)
        ier=0
        nmap=0
        larrmax=0
          icount=1
        do 4000 iact=1,nact
c 
         icount=icount+1
         if(icount .eq. 10) icount=0
         if(icount .eq. 0) call prinf('iact=*',iact,1)
  
  
c 
c        if this is a compression of a new lowest level chunk
c        - generate and compress it
c 
        if(actions(1,iact) .ne. 1) goto 2000
c 
        k=actions(2,iact)
        j1=boxes(2,k)
        j2=boxes(3,k)
        i1=boxes(4,k)
        i2=boxes(5,k)
        mi=i2-i1+1
        nj=j2-j1+1
c 
        itype=boxes(11,k)
c 
c       . . . create the lowest level submatrix and compress it
c 
        call crecomlo(funmat,j1,j2,i1,i2,w,udad,vdad,
     1      ncols,eps,itype,rnorms,w1,w2,w3,w4,w5)
c 
cccccc          call prinf('after crecomlo, ncols=*',ncols,1)
        if(ncols .le. maxcols) goto 1600
        ier=1024
        return
 1600 continue
c 
c       store the newly created compressed form of the submatrix
c 
        call chnkstor(jer,map,arr,nmap,
     1      k,nj,mi,ncols,itype,udad,vdad,larr,lenarr)
        if(larrmax .lt. larr) larrmax=larr
        if(jer .eq. 0) goto 1800
        if(jer. eq. 4) ier=32000
        if(jer .eq. 1024) ier=256
        return
 1800 continue
cccccc         call prinf('after chnkstor, larr=*',larr,1)
c 
        goto 4000
 2000 continue
c 
c       this action is the merging of two boxes. extract the first box
c       from array arr
c 
        kson1=actions(3,iact)
        kson2=actions(4,iact)
        kdad=actions(2,iact)
        ib11=boxes(11,kson1)
c 
c        extract the address of the expansion for the
c        first son in array arr
c 
        call chnkreti(jer1,map,nmap,arr,kson1,
     1      nson1,mson1,ncolson1,itype1,iuson1,ivson1)
c 
c        extract the address of the expansion for the
c        second son in array arr
c 
        call chnkreti(ier2,map,nmap,arr,kson2,
     1      nson2,mson2,ncolson2,itype2,iuson2,ivson2)
c 
c        merge the compressed forms of the sonnies to obtain
c        the daddy's compressed form
c 
        if( (jer1 .eq. 0) .and. (jer2 .eq. 0) ) goto 2200
        ier=32000
        return
 2200 continue
c 
        if(ib11 .eq. 1)
     1      call combhor(arr(iuson1),arr(ivson1),
     2      arr(iuson2),arr(ivson2),
     2              nson1,mson1,mson2,ncolson1,ncolson2,
     3              udad,vdad,ncolsdad,eps,w)
c 
        if(ib11 .eq. 2)
     1      call combvert(arr(iuson1),arr(ivson1),
     2      arr(iuson2),arr(ivson2),
     2              nson1,nson2,mson1,ncolson1,ncolson2,
     3              udad,vdad,ncolsdad,eps,w)
c 
ccccc          call prinf('after combvert, ncolsdad=*',ncolsdad,1)
        if(ncolsdad .le. maxcols) goto 2400
        ier=1024
        return
 2400 continue
c 
c         remove the sonnies' corpses from the array arr
c 
        call chnkrem(jer1,map,arr,nmap,kson1,larr)
        call chnkrem(jer2,map,arr,nmap,kson2,larr)
c 
        if( (jer1 .eq. 0) .and. (jer2 .eq. 0) ) goto 2600
        ier=32000
        return
 2600 continue
c 
c 
c        store daddy in array arr
c 
c 
        ndad=boxes(3,kdad)-boxes(2,kdad)+1
        mdad=boxes(5,kdad)-boxes(4,kdad)+1
c 
        itypedad=1
        if(ib11 .eq. 1) itypedad=2
c 
c        if the type of the daddy's expansion does not
c        correspond to the type of compression it will
c        undergo - change the type
c 
        ib11dad=boxes(11,kdad)
        if(ib11dad .eq. itypedad) goto 3600
c 
c       store in the array map the information about
c       u, v to be stored
c 
        call chnkstoi(jer,map,arr,nmap,kdad,ndad,
     1       mdad,ncolsdad,ib11dad,iu7,iv7,larr,lenarr)
        if(larrmax .lt. larr) larrmax=larr
        if(jer .eq. 0) goto 3400
        if(jer. eq. 4) ier=32000
        if(jer .eq. 1024) ier=256
        return
 3400 continue
c        change the type of the expansion and store
c        the result in array arr
c 
        if(ib11dad .eq. 1)
     1      call makevort(udad,vdad,ndad,mdad,ncolsdad,
     2      arr(iu7),arr(iv7),w)
c 
        if(ib11dad .eq. 2)
     1      call makeuort(udad,vdad,ndad,mdad,ncolsdad,
     2      arr(iu7),arr(iv7),w)
        goto 4000
 3600 continue
c 
c        the type of the daddy's expansion does correspond
c        to the type of compression it will undergo
c        store it in array arr as is
c 
        call chnkstor(jer,map,arr,nmap,kdad,ndad,
     1      mdad,ncolsdad,ib11dad,udad,vdad,larr,lenarr)
        if(larrmax .lt. larr) larrmax=larr
        if(jer .eq. 0) goto 3800
        if(jer. eq. 4) ier=32000
        if(jer .eq. 1024) ier=256
        return
 3800 continue
c 
 4000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arrcopy(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine crecomlo(funmat,j1,j2,i1,i2,a,u,v,
     1      ncols,eps,itype,rnorms,w1,w2,w3,w4,w5)
        implicit real *8 (a-h,o-z)
        save
        dimension a(j2-j1+1,1),u(j2-j1+1,1),v(i2-i2+1,1),
     1      w1(1),w2(1),w3(1),w4(1),w5(1),rnorms(1)
c 
c        this subroutine constructs the submatrix on the final
c        level and compresses  it.
c 
c 
c                 input parameters:
c 
c  funmat - the user-provided subroutine for the evaluation of the
c        elements of the matrix to be compressed.
c  j1 - the number in the big matrix of the first row
c        of the submatrix to be created
c  j2 - the number in the big matrix of the last row
c        of the submatrix to be created
c  i1 - the number in the big matrix of the first column
c        of the submatrix to be created
c  i2 - the number in the big matrix of the last column
c        of the submatrix to be created
c  eps - the accuracy to which the computations are to
c        be conducted
c  itype - the type of compression to be produced:
c        itype=1 means that the columns of u will be orthonormal
c        itype=2 means that the columns of v will be orthonormal
c 
c  w1,w2,w3,w4,w5 - whatever parameters the user's subroutine funmat
c                needs (real, integer, array, variable - whatever),
c 
c                 output parameters:
c 
c  u - the matrix u in the compressed form of the submatrix;
c          must be dimensioned at least (j2-j1+1,i2-i1+1)
c  v - the matrix v in the compressed form of the submatrix;
c          must be dimensioned at least (j2-j1+1,i2-i1+1)
c  ncols - the number of columns in the arrays u, v actually
c          generated (hopefully, less much less than i2-i2+1)
c 
c                 work arrays:
c 
c  a - must be at least (j2-j1+1)*(i2-i1+1) real *8 elements long
c  rnorms - must be at least max(j2-j1+1,i2-i1+1)+2 real *8 elements long
c 
c 
c       . . . construct the submatrix to be compressed
c 
        m=i2-i1+1
        n=j2-j1+1
        do 2000 i=1,m
        do 1800 j=1,n
        call funmat(j+j1-1,i+i1-1,a(j,i),w1,w2,w3,w4,w5)
 1800 continue
 2000 continue
c 
c       compress the newly constructed submatrix
c 
        ifpivot=1
c 
c        if this box is to be merged vertically - compress
c        it in such a fashion that the columns of the
c        matrix u are orthonormal
c 
         if(itype. eq. 2) call grmright(a,n,m,u,v,
     1      ncols,rnorms,eps,ifpivot)
c 
c        if this box is to be merged horizontally - compress
c        it in such a fashion that the columns of the
c        matrix v are orthonormal
c 
         if(itype. eq. 1) call grmleft(a,n,m,u,v,
     1      ncols,rnorms,eps,ifpivot)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chnkstoi(ier,map,arr,nmap,
     1      k,n,m,ncols,itype,iu7,iv7,larr,lenarr)
        implicit real *8 (a-h,o-z)
        save
        dimension arr(1),u(1),v(1)
        integer *4 map(10,1)
c 
c        this subroutine stores in the array map the information
c        about an area in the array arr to contain
c        the standard compression of the matrix u * v^*.
c        it is assumed that the user will store the actual
c        compressed form in the array arr directly, starting
c        with positions iu7, iv7 provided by this subroutine.
c        once stored, the information can be retrieved by the
c        entry chnkretr (see below), or deleted by the entry
c        chnkrem (see below).
c 
c                   input parameters:
c 
c  map - the map of the storage area BEFORE this subroutine
c        altered it
c 
c        IMPORTANT: description of array map - the map of the storage
c        area arr
c 
c       map(1,i) - the address in the array boxes of the box
c          whose storage is described by this entry in array map (k)
c       map(2,i) - the number of rows in this submatrix (n)
c       map(3,i) - the number of columns in this submatrix (m)
c       map(4,i) - the rank of this submatrix (ncols)
c       map(5,i) - the location in array arr of the left part (u)
c          of the compressed form of this submatrix
c       map(6,i) - the length in array arr of the left part (u)
c          of the compressed form of this submatrix
c       map(7,i) - the location in array arr of the right part (v)
c          of the compressed form of this submatrix
c       map(8,i) - the length in array arr of the right part (v)
c          of the compressed form of this submatrix
c       map(9,i) - the type of the expansion being stored:
c        map(9,i)=1 means that the columns of u are orthonormal
c        map(9,i)=2 means that the columns of v are orthonormal
c  arr - the main storage area BEFORE this subroutine altered it
c  nmap - the number of pairs (u,v) stored in the arrays map, arr
c          BEFORE this call to this subroutine; this is also the
c          number of columns of the array map in use before this call
c  k - the label of this storage area. conceptually, the sequence number
c          in the array BOXES of the submatrix whose compression is stored
c          by this call.
c  n - the number of rows in the matrix u
c  m - the number of rows in the matrix v
c  ncols - the number of columns in each of the matrices u, v
c  itype - the type of the compression used on this box:
c        itype=1 means that the columns of u are orthonormal
c        itype=2 means that the columns of v are orthonormal
c  u - the n*ncols matrix to be stored by this call
c  v - the m*ncols matrix to be stored by this call
c  lenarr - the total amount of storage available in array arr
c          (in terms of real *8 words)
c 
c                     output parameters:
c 
c  ier - error return code.
c          ier=0 means successful conclusion
c          ier=4 means that an attempt has been made to store
c                arrays u, v with a duplicate label. in this
c                case, the request has been ignored.
c          ier=1024 means that the storage allocation (given by
c                the input parameter lenarr) would be exceeded  if
c                this storage were executed. the request ignored.
c  arr - the main storage area AFTER this subroutine altered it
c  map - the map of the storage area AFTER this subroutine
c        altered it
c  larr - the number of real *8 elements in array arr that are in
c        use on exit from this subroutine
c 
c        . . . if an attempt is being made to store arrays
c              with a duplicate label - bomb out
c 
        ier=0
        if(nmap .eq. 0) goto 300
        do 200 i=1,nmap
        if(k .eq. map(1,i)) ier=4
  200 continue
        if(ier .ne. 0) return
  300 continue
c 
c        . . . store the matrices u, v in the array arr
c 
        nm=nmap+1
        if(nmap .ne. 0) iu=map(7,nmap)+map(8,nmap)+1
        if(nmap .eq. 0) iu=1
        nu=ncols*n
        iu7=iu
c 
        iv=iu+nu
        nv=ncols*m
        iv7=iv
c 
        if(iv+nv .lt. lenarr) goto 400
        ier=1024
        return
  400 continue
c 
c       enter the information about the last event in the
c       array map
c 
        map(1,nm)=k
        map(2,nm)=n
        map(3,nm)=m
        map(4,nm)=ncols
        map(5,nm)=iu
        map(6,nm)=nu
        map(7,nm)=iv
        map(8,nm)=nv
        map(9,nm)=itype
c 
        nmap=nm
        larr=map(7,nm)+map(8,nm)
c 
        return
c 
c 
c 
c 
        entry chnkstor(ier,map,arr,nmap,
     1      k,n,m,ncols,itype,u,v,larr,lenarr)
c 
c        this subroutine stores in the array arr the matrices u, v
c        of the standard compression of the matrix u * v^*.
c        it also stores in the array map the map of the storage
c        areas in arr. once stored by this entry, the information
c        can be retrieved by the entry chnkretr (see below), or
c        deleted by the entry chnkrem (see below).
c 
c                   input parameters:
c 
c  map - the map of the storage area BEFORE this subroutine
c        altered it
c 
c        IMPORTANT: description of array map - the map of the storage
c        area arr
c 
c       map(1,i) - the address in the array boxes of the box
c          whose storage is described by this entry in array map (k)
c       map(2,i) - the number of rows in this submatrix (n)
c       map(3,i) - the number of columns in this submatrix (m)
c       map(4,i) - the rank of this submatrix (ncols)
c       map(5,i) - the location in array arr of the left part (u)
c          of the compressed form of this submatrix
c       map(6,i) - the length in array arr of the left part (u)
c          of the compressed form of this submatrix
c       map(7,i) - the location in array arr of the right part (v)
c          of the compressed form of this submatrix
c       map(8,i) - the length in array arr of the right part (v)
c          of the compressed form of this submatrix
c 
c  arr - the main storage area BEFORE this subroutine altered it
c  nmap - the number of pairs (u,v) stored in the arrays map, arr
c          BEFORE this call to this subroutine; this is also the
c          number of columns of the array map in use before this call
c  k - the label of this storage area. conceptually, the sequence number
c          in the array BOXES of the submatrix whose compression is stored
c          by this call.
c  n - the number of rows in the matrix u
c  m - the number of rows in the matrix v
c  ncols - the number of columns in each of the matrices u, v
c  itype - the type of the compression used on this box:
c        itype=1 means that the columns of u are orthonormal
c        itype=2 means that the columns of v are orthonormal
c  u - the n*ncols matrix to be stored by this call
c  v - the m*ncols matrix to be stored by this call
c  lenarr - the total amount of storage available in array arr
c          (in terms of real *8 words)
c 
c                     output parameters:
c 
c  ier - error return code.
c          ier=0 means successful conclusion
c          ier=4 means that an attempt has been made to store
c                arrays u, v with a duplicate label. in this
c                case, the request has been ignored.
c          ier=1024 means that the storage allocation (given by
c                the input parameter lenarr) would be exceeded  if
c                this storage were executed. the request ignored.
c  arr - the main storage area AFTER this subroutine altered it
c  map - the map of the storage area AFTER this subroutine
c        altered it
c  larr - the number of real *8 elements in array arr that are in
c        use on exit from this subroutine
c 
c        . . . if an attempt is being made to store arrays
c              with a duplicate label - bomb out
c 
        ier=0
        if(nmap .eq. 0) goto 1300
        do 1200 i=1,nmap
        if(k .eq. map(1,i)) ier=4
 1200 continue
        if(ier .ne. 0) return
 1300 continue
c 
c        . . . store the matrices u, v in the array arr
c 
        nm=nmap+1
        if(nmap .ne. 0) iu=map(7,nmap)+map(8,nmap)+1
        if(nmap .eq. 0) iu=1
        nu=ncols*n
c 
        if(iu+nu .lt. lenarr) goto 1350
        ier=1024
        return
 1350 continue
c 
        do 1400 i=1,nu
        arr(iu+i-1)=u(i)
 1400 continue
c 
        iv=iu+nu
        nv=ncols*m
c 
        if(iv+nv .lt. lenarr) goto 1500
        ier=1024
        return
 1500 continue
c 
        do 1600 i=1,nv
        arr(iv+i-1)=v(i)
 1600 continue
c 
c       enter the information about the last event in the
c       array map
c 
        map(1,nm)=k
        map(2,nm)=n
        map(3,nm)=m
        map(4,nm)=ncols
        map(5,nm)=iu
        map(6,nm)=nu
        map(7,nm)=iv
        map(8,nm)=nv
        map(9,nm)=itype
c 
        nmap=nm
        larr=map(7,nm)+map(8,nm)
c 
        return
c 
c 
c 
c 
        entry chnkrem(ier,map,arr,nmap,k,larr)
c 
c        this entry deletes from array arr the pair of
c        matrices u, v that has been previously stored
c        there by the entry chnkstor (see above). it also
c        deletes from the array map the information about
c        the area deleted from the array arr.
c 
c 
c                   input parameters:
c 
c  map - the map of the storage area BEFORE this subroutine
c        altered it
c 
c  arr - the main storage area BEFORE this subroutine altered it
c  nmap - the number of pairs (u,v) stored in the arrays map, arr
c          BEFORE this call to this subroutine; this is also the
c          number of columns of the array map in use before this call
c  k - the label of the storage area to be deleted. conceptually,
c           the sequence number in the array BOXES of the submatrix
c           whose compression is to be deleted by this call.
c 
c                     output parameters:
c 
c  arr - the main storage area AFTER this subroutine altered it
c  map - the map of the storage area AFTER this subroutine
c        altered it
c  larr - the number of real *8 elements in array arr that are in
c        use on exit from this subroutine
c 
c       . . . find the location in the array map of the map of
c             the submatrix to be removed
c 
        ier=0
        do 2200 i=1,nmap
        ii=i
        if(map(1,i) .eq. k) goto 2400
 2200 continue
        ier=4
        return
 2400 continue
c 
c       if there is only one box stored in array arr - do
c       not compress anything
c 
        if(nmap .ne. 1) goto 2500
        nmap=0
        return
 2500 continue
c 
c       delete in array arr the storage areas corresponding to
c       the submatrix whose address in array boxes is k
c 
        iu=map(5,ii)
        iv=map(7,ii)
        nu=map(6,ii)
        nv=map(8,ii)
        ishift=nu+nv
        nn=map(7,nmap)+map(8,nmap)-iu+10
c 
        do 2600 i=1,nn
        arr(iu+i-1)=arr(iu+i-1+ishift)
 2600 continue
c 
c       adjust all addresses in array map that are subsequent
c       to the one that has been changed
c 
        do 2800 i=ii+1,nmap
        map(5,i)=map(5,i)-ishift
        map(7,i)=map(7,i)-ishift
 2800 continue
c 
c       eliminate column number ii from the array map
c 
        do 3200 i=ii,nmap
        do 3000 j=1,10
        map(j,i)=map(j,i+1)
 3000 continue
 3200 continue
        nmap=nmap-1
        larr=map(7,nmap)+map(8,nmap)
        return
c 
c 
c 
c 
        entry chnkretr(ier,map,nmap,arr,k,
     1      n,m,ncols,itype,u,v)
c 
c        this entry retrieves from the array arr the matrices u, v
c        of the standard compression of the matrix u * v^*.
c        it also retrieves from array map the description of
c        the matrices u, v. it is assumed that this information
c        has been previously stored by the entry chnkstor
c        (see above).
c 
c                   input parameters:
c 
c  map - the map of the storage area.
c  arr - the main storage area
c  nmap - the number of pairs (u,v) stored in the arrays map, arr.
c  k - the label of the matrices u, v to be retrieved.
c          conceptually, the sequence number in the array BOXES
c          of the submatrix whose compression is being retrieved
c          by this call.
c 
c                     output parameters:
c 
c  n - the number of rows in the matrix u
c  m - the number of rows in the matrix v
c  ncols - the number of columns in each of the matrices u, v
c  itype - the type of the compression used on this box:
c        itype=1 means that the columns of u are orthonormal
c        itype=2 means that the columns of v are orthonormal
c  u - the n*ncols matrix retrieved by this call
c  v - the m*ncols matrix retrieved by this call
c 
c       . . . find in the array map the column describing the submatrix
c             number k
c 
        ier=0
        do 4200 i=1,nmap
        ii=i
        if(map(1,i) .eq. k) goto 4400
 4200 continue
        ier=4
        return
 4400 continue
c 
c       retrieve from array arr the matrices u,v corresponding
c       to the chunk k
c 
        iu=map(5,ii)
        iv=map(7,ii)
        nu=map(6,ii)
        nv=map(8,ii)
        itype=map(9,ii)
c 
        do 4600 i=1,nu
        u(i)=arr(iu+i-1)
 4600 continue
c 
        do 4800 i=1,nv
        v(i)=arr(iv+i-1)
 4800 continue
c 
        n=map(2,ii)
        m=map(3,ii)
        ncols=map(4,ii)
        return
c 
c 
c 
c 
        entry chnkreti(ier,map,nmap,arr,k,
     1      n,m,ncols,itype,iu7,iv7)
c 
c        this entry retrieves from the map the location
c        in array arr of the the matrices u, v
c        of the standard compression of the matrix u * v^*.
c        it also retrieves from array map the description of
c        the matrices u, v. it is assumed that this information
c        has been previously stored by the entry chnkstor
c        (see above).
c 
c                   input parameters:
c 
c  map - the map of the storage area.
c  arr - the main storage area
c  nmap - the number of pairs (u,v) stored in the arrays map, arr.
c  k - the label of the matrices u, v to be retrieved.
c          conceptually, the sequence number in the array BOXES
c          of the submatrix whose compression is being retrieved
c          by this call.
c 
c                     output parameters:
c 
c  n - the number of rows in the matrix u
c  m - the number of rows in the matrix v
c  ncols - the number of columns in each of the matrices u, v
c  itype - the type of the compression used on this box:
c        itype=1 means that the columns of u are orthonormal
c        itype=2 means that the columns of v are orthonormal
c  iu7 - the location in array arr of the matrix u
c  iv7 - the location in array arr of the matrix u
c 
c       . . . find in the array map the column describing the submatrix
c             number k
c 
        ier=0
        do 5200 i=1,nmap
        ii=i
        if(map(1,i) .eq. k) goto 5400
 5200 continue
        ier=4
        return
 5400 continue
c 
c       retrieve from array arr the matrices u,v corresponding
c       to the chunk k
c 
        iu=map(5,ii)
        iv=map(7,ii)
        nu=map(6,ii)
        nv=map(8,ii)
        itype=map(9,ii)
c 
        iu7=iu
        iv7=iv
c 
        n=map(2,ii)
        m=map(3,ii)
        ncols=map(4,ii)
        return
        end
  
c 
c 
c 
c 
c 
        subroutine boxesbld(n,m,nrowmax,boxes,nboxes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(12,1)
c 
c        this subroutine subdivides the user-supplied matrix
c        (actually, only dimensionalities n,m are used) into
c        a quad-tree structure, the leaves of the tree being
c        submatrices of a size no greater that nrowmax, and
c        as close to nrowmax as possible.
c 
c                 input parameters:
c 
c  n - the number of rows of the matrix to be subdivided
c  m - the number of columns of the matrix to be subdivided
c  nrowmax - the maximum number of rows (and columns) that
c        will be permitted at the finest level of subdivision
c 
c                 output parameters:
c 
c  boxes - integer array dimensioned boxes(12,nboxes) describing
c        the tree structure (the subdivision of the matrix). each
c        column of the array boxes has 10 entries, as follows
c 
c       boxes(1,i) - the location of the i-th box in the array boxes
c       boxes(2,i) - the location of the first row of the i-th
c                      box in the big matrix
c       boxes(3,i) - the location of the last row of the i-th
c                      box in the big matrix
c       boxes(4,i) - the location of the first column of the i-th
c                      box in the big matrix
c       boxes(5,i) - the location of the last column of the i-th
c                      box in the big matrix
c       boxes(6,i) - the location of the daddy of this box
c                      in array boxes
c       boxes(7,i) - the location of the first son of this box
c                      in array boxes
c       boxes(8,i) - the location of the second son of this box
c                      in array boxes
c       boxes(9,i) - the slot to be used by the subroutine acttree
c                      to keep track of thether this box has been
c                      processed
c       boxes(10,i) - the level of subdivision of this box
c       boxes(11,i) - the type of merging this box will undergo:
c                      boxes(11,k)=1 means that the box will be merged
c                           horizontally when its turn comes
c                      boxes(11,k)=2 means that the box will be merged
c                           vertically when its turn comes
c       boxes(12,i) - not used in this version
c 
c  nboxes - the number of columns in array boxes, i.e. the total
c       number of submatrices created on all levels.
c 
c        note that this subroutine only creates the entries 1-8 for
c        each box. the entries 9, 10 will be created later.
c 
c       . . . construct the main box
c 
        boxes(1,1)=1
        boxes(2,1)=1
        boxes(3,1)=n
        boxes(4,1)=1
        boxes(5,1)=m
        boxes(6,1)=0
        boxes(7,1)=2
        boxes(8,1)=3
        boxes(9,1)=0
        boxes(10,1)=0
        boxes(11,1)=0
        boxes(12,1)=0
c 
c         subdivide all boxes recursively until none are
c         left with dimensions
c 
        nboxes=1
        nbold=1
        icurr=2
        iw=30
        do 2000 i=1,100000
c 
        do 1800 j=nbold,nboxes
c 
        do 1200 l=1,12
        boxes(l,icurr)=0
        boxes(l,icurr+1)=0
 1200 continue
        iw=iw+1
c 
c       attempt to subdivide this box
c 
        call boxbld(boxes(1,j),nrowmax,boxes(1,icurr),
     1      boxes(1,icurr+1),ifdiv)
c 
c       if this box has been subdivided - enter the results in the
c       array boxes
c 
        if(ifdiv .eq. 0) goto 1400
        boxes(7,j)=icurr
        icurr=icurr+1
        boxes(8,j)=icurr
        icurr=icurr+1
c 
        boxes(6,icurr-2)=j
        boxes(6,icurr-1)=j
c 
        goto 1800
c 
 1400 continue
c 
c       this box has not been subdivided. enter in the array
c       boxes the information that it is childless
c 
        boxes(7,j)=-7
 1800 continue
        nbold=nboxes+1
        nboxes=icurr
         if(nboxes .le. nbold) goto 2200
 2000 continue
c 
 2200 continue
        nboxes=nboxes-1
c 
        do 2400 i=1,nboxes
        boxes(1,i)=i
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine boxbld(box,nrowmax,son1,son2,ifdiv)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(1),son1(1),son2(1)
c 
c       box(1) - the location of the this box in the array boxes
c       box(2) - the location of the first row of the this
c                      box in the big matrix
c       box(3) - the location of the last row of the this
c                      box in the big matrix
c       box(4) - the location of the first column of the this
c                      box in the big matrix
c       box(5) - the location of the last column of the this
c                      box in the big matrix
c       box(6) - the location of the daddy of this box
c                      in array boxes
c       box(7) - the location of the first son of this box
c                      in array boxes
c       box(8) - the location of the second son of this box
c                      in array boxes
c       box(9) - the slot to be used by the subroutine acttree
c                      to keep track of thether this box has been
c                      processed
c       box(10) - the level of subdivision of this box
c       box(11) - the type of merging this box will undergo:
c                      box(11)=1 means that the box will be merged
c                           horizontally when its turn comes
c                      box(11)=2 means that the box will be merged
c                           vertically when its turn comes
c       box(11) - not used in this version
c 
c 
c       . . . determine if this box should be subdivided
c 
        ifdiv=0
        nrows=box(3)-box(2)+1
        ncols=box(5)-box(4)+1
c 
        if(nrows. gt. nrowmax) ifdiv=1
        if(ncols. gt. nrowmax) ifdiv=1
        if(ifdiv .eq. 0) box(7)=-7
        if(ifdiv .eq. 0) return
c 
c       this box should be subdivided. if it should be subdivided
c       vertically - do so.
c 
        if(ncols .gt. nrows) goto 2000
        nn=nrows/2
        son1(2)=box(2)
        son1(3)=son1(2)+nn-1
c 
        son2(2)=son1(3)+1
        son2(3)=box(3)
c 
        son1(4)=box(4)
        son1(5)=box(5)
        son2(4)=box(4)
        son2(5)=box(5)
c 
        son1(7)=0
        son1(8)=0
        son2(7)=0
        son2(8)=0
c 
        son1(10)=box(10)+1
        son2(10)=son1(10)
c 
        son1(11)=2
        son2(11)=2
        return
 2000 continue
c 
c       this box should be subdivided
c       horizontally - do so.
c 
        nn=ncols/2
        son1(2)=box(2)
        son1(3)=box(3)
c 
        son2(2)=box(2)
        son2(3)=box(3)
c 
        son1(4)=box(4)
        son1(5)=son1(4)+nn-1
        son2(4)=son1(5)+1
        son2(5)=box(5)
c 
        son1(7)=0
        son1(8)=0
        son2(7)=0
        son2(8)=0
c 
        son1(10)=box(10)+1
        son2(10)=son1(10)
c 
        son1(11)=1
        son2(11)=1
        return
        end
c 
c 
c 
c 
c 
        subroutine acttree(boxes,nboxes,actions,nact)
        save
        integer *4 boxes(12,1),actions(4,1)
c 
c        this subroutine constructs the action table controlling
c        the recursive compression of a matrix.
c 
c                     input parameters:
c 
c  boxes - the table of boxes as produced by the subroutine boxesbld (see)
c  nboxes - the number of columns in the array boxes - also produced
c       by boxesbld
c 
c                     output parameters:
c 
c  actions - integer *4 array dimensioned (4,nact). in fact, this is a
c       table for the  actual compression of a matrix, with each column
c       describing one action, in the following manner:
c 
c       actions(1,i) describes the type of action. action(1,i)=1
c                    means that a lowest level submatrix will be
c                    constructed and  QR - decomposed directly.
c                    action(1,i)=2 means that two submatrices in
c                    the QR - form will be merged.
c       actions(2,i) is the address in the array boxes of the
c                    submatrix whose QR-decomposition is being created.
c       actions(3,i),actions(4,i) are addresses in array boxes of the
c                    submatrices whose QR - forms are being merged.
c                    these two entries are only used if actions(1,i)=1
c  nact - the number of columns in array actions
c 
c       . . . process one box after another, going up and
c             down the array boxes
c 
        k=1
        nact=0
        icount=0
        do 4000 i=1,100 000
        icount=icount+1
        if(icount .eq. 10) icount=0
c 
c       if this box has been compressed - go to his daddy
c 
        if(boxes(9,k) .le. 0) goto 1400
        k=boxes(6,k)
        goto 3000
 1400 continue
        if(icount .eq. 0) call prinf('i=*',i,1)
c 
c       this box has not been compressed. if it is childless,
c       generate and compress it directly; after that, go to
c       his daddy
c 
        if (boxes(7,k) .gt. 0) goto 1600
        nact=nact+1
        actions(1,nact)=1
        actions(2,nact)=k
c 
        boxes(9,k)=7
c 
        k=boxes(6,k)
        goto 3000
 1600 continue
c 
c       this box has children. if both of them have been
c       compressed, combine them to obtain the compression
c       of this box. after that, go to his daddy
c 
        kson1=boxes(7,k)
        kson2=boxes(8,k)
        if(boxes(9,kson1) .le. 0) goto 1800
        if(boxes(9,kson2) .le. 0) goto 1800
c 
        nact=nact+1
        actions(1,nact)=2
        actions(2,nact)=k
        actions(3,nact)=kson1
        actions(4,nact)=kson2
c 
        boxes(9,k)=7
        k=boxes(6,k)
        goto 3000
 1800 continue
c 
c       at least one of the sons has not been compressed.
c       goto it
c 
        if(boxes(9,kson2) .le. 0) k=kson2
        if(boxes(9,kson1) .le. 0) k=kson1
 3000 continue
c 
c        if the whole matrix has been compressed - exit
c        the subroutine
c 
         if (boxes(9,1) .gt. 0)  return
 4000 continue
       return
        end
  
c 
c 
c 
c 
c 
  
        subroutine combvert(u1,v1,u2,v2,n1,n2,m,
     1    ncols1,ncols2,u3,v3,ncols3,eps,w)
        implicit real *8 (a-h,o-z)
        save
        dimension u1(n1,ncols1),v1(m,ncols1),
     1      u2(n2,ncols2),v2(m,ncols2),
     2      u3(n1+n2,1),v3(m,1),w(1)
c 
c       this subroutine combines two matrices in compressed
c       form that have been put side by side vertically.
c       more specifically, the user supplies two matrices
c 
c         a1= u1 * v1^*,                             (1)
c         a2= u2 * v2^*,                             (2)
c 
c       and this subroutine obtains the matrix
c 
c             a1
c         a3=     =u3 * v3^*.                       (3)
c             a2
c 
c       the point of this subroutine is that the numerical
c       ranks of a1, a2, a3 are expected to be much smaller
c       than their dimensionalities, so that  (1), (2), (3)
c       are much compressed forms of a1, a2, a3. in fact, this
c       is simply a memory management routine; all the actual
c       work is done by the subroutine combhor0 (see)
c 
c  PLEASE NOTE: the subroutine assumes that in the representations
c  (1), (2), the columns in each of the matrices u1, u2, are
c  orthonormal. However, on exit the columns of the matrix v3
c  are orthonormal, not those of u3.
c 
c                     input parameters:
c 
c  u1 - n1 * ncols1 matrix in (1)
c  v1 - m * ncols1 matrix in (1)
c  u2 - n2 * ncols2 matrix in (2)
c  v2 - m * ncols2 matrix in (2)
c  n1 - the number of rows in the matrices a1,u1
c  n2 - the number of rows in the matrices a2,u2
c  m - the number of columns in the matrics a1, a2, a3; also the
c       number of rows in the matrices v1,v2,v3
c  ncols1 - the rank of the matrix a1; also the number of columns in
c       the matrices u1,v1
c  ncols2 - the rank of the matrix a2; also the number of columns in
c       the matrices u2,v2
c  eps - the precision of computations
c 
c                     output parameters:
c 
c  u3 - the (n1+n2)*ncols3 matrix in (3)
c  v3 - the  m*ncols3 matrix in (3)
c  ncols3 - the numerical rank of the matrix  a3 in (3);
c        also the number of columns in the matrices u3, v3
c 
c                     work arrays:
c 
c  w - must be at least (ncols1+ncols2)*m +2*(ncols1+ncols2)**2+
c        (ncols1+ncols2)+50
c 
c      . . . transpose the problem and use the horizontal
c            compressor combhor
c 
        call combhor(v1,u1,v2,u2,m,n1,n2,
     1    ncols1,ncols2,v3,u3,ncols3,eps,w)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine combhor(u1,v1,u2,v2,n,m1,m2,
     1    ncols1,ncols2,u3,v3,ncols3,eps,w)
        implicit real *8 (a-h,o-z)
        save
        dimension u1(n,ncols1),v1(m1,ncols1),
     1      u2(n,ncols2),v2(m2,ncols2),
     2      u3(n,1),v3(m1+m2,1),w(1)
c 
c       this subroutine combines two matrices in compressed
c       form that have been put side by side horizontally.
c       more specifically, the user supplies two matrices
c 
c         a1= u1 * v1^*,                             (1)
c         a2= u2 * v2^*,                             (2)
c 
c       and this subroutine obtains the matrix
c 
c         a3=(a1a2)=u3 * v3^*.                       (3)
c 
c       the point of this subroutine is that the numerical
c       ranks of a1, a2, a3 are expected to be much smaller
c       than their dimensionalities, so that  (1), (2), (3)
c       are much compressed forms of a1, a2, a3. in fact, this
c       is simply a memory management routine; all the actual
c       work is done by the subroutine combhor0 (see)
c 
c  PLEASE NOTE: the subroutine assumes that in the representations
c  (1), (2), the columns in each of the matrices v1, v2, are
c  orthonormal. However, on exit the columns of the matrix u3
c  are orthonormal, not those of v3.
c 
c                     input parameters:
c 
c  u1 - n * ncols1 matrix in (1)
c  v1 - m1 * ncols1 matrix in (1)
c  u2 - n * ncols2 matrix in (2)
c  v2 - m2 * ncols2 matrix in (2)
c  n - the number of rows in the matrices a1,a2,u1,u2,a3,u3
c  m1 - the number of columns in the matrix a1; also the number
c       of rows in the matrix v1
c  m2 - the number of columns in the matrix a2; also the number
c       of rows in the matrix v2
c  ncols1 - the rank of the matrix a1; also the number of columns in
c       the matrices u1,v1
c  ncols2 - the rank of the matrix a2; also the number of columns in
c       the matrices u2,v2
c  eps - the precision of computations
c 
c                     output parameters:
c 
c  u3 - the n*ncols3 matrix in (3)
c  v3 - the (m1+m2)*ncols3 matrix in (3)
c  ncols3 - the numerical rank of the matrix  a3 in (3);
c        also the number of columns in the matrices u3, v3
c 
c                     work arrays:
c 
c  w - must be at least (ncols1+ncols2)*n +2*(ncols1+ncols2)**2+
c        (ncols1+ncols2)+50
c 
c      . . . allocate memory for the subroutine combhor0
c 
cccc        call orttes(v1,m1,ncols1,u3,'in combhor, error in v1*')
cccc        call orttes(v2,m2,ncols2,u3,'in combhor, error in v2*')
        iw=1
        lw=(ncols1+ncols2)*n+10
c 
        ialpha=iw+lw
        lalpha=(ncols1+ncols2)*ncols1 +10
c 
        ibeta=ialpha+lalpha
        lbeta=(ncols1+ncols2)*ncols2 +10
c 
        ialpha0=ibeta+lbeta
        lalpha0=(ncols1+ncols2)**2+10
c 
        irnorms=ialpha0+lalpha0
        lrnorms=ncols1+ncols2+2
c 
c       combine the two matrices
c 
        call combhor0(u1,v1,u2,v2,n,m1,m2,
     1    ncols1,ncols2,u3,v3,ncols3,w(iw),w(ialpha),
     2    w(ibeta),w(ialpha0),eps,w(irnorms) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine combhor0(u1,v1,u2,v2,n,m1,m2,
     1    ncols1,ncols2,u3,v3,ncols3,w,alpha,beta,
     2    alpha0,eps,rnorms)
        implicit real *8 (a-h,o-z)
        save
        dimension u1(n,ncols1),v1(m1,ncols1),
     1      u2(n,ncols2),v2(m2,ncols2),w(n,1),
     2      alpha(ncols1,1),beta(ncols2,1),
     3      alpha0(ncols1+ncols2,1),u3(n,1),
     4       v3(m1+m2,1),rnorms(1)
c 
c       this subroutine combines two matrices in compressed
c       form that have been put side by side horizontally.
c       more specifically, the user supplies two matrices
c 
c         a1= u1 * v1^*,                             (1)
c         a2= u2 * v2^*,                             (2)
c 
c       and this subroutine obtains the matrix
c 
c         a3=(a1a2)=u3 * v3^*.                       (3)
c 
c       the point of this subroutine is that the numerical
c       ranks of a1, a2, a3 are expected to be much smaller
c       than their dimensionalities, so that  (1), (2), (3)
c       are much compressed forms of a1, a2, a3.
c 
c  PLEASE NOTE: the subroutine assumes that in the representations
c  (1), (2), the columns in each of the matrices v1, v2, are
c  orthonormal. However, on exit the columns of the matrix u3
c  are orthonormal, not those of v3.
c 
c                     input parameters:
c 
c  u1 - n * ncols1 matrix in (1)
c  v1 - m1 * ncols1 matrix in (1)
c  u2 - n * ncols2 matrix in (2)
c  v2 - m2 * ncols2 matrix in (2)
c  n - the number of rows in the matrices a1,a2,u1,u2,a3,u3
c  m1 - the number of columns in the matrics a1; also the number
c       of rows in the matrix v1
c  m2 - the number of columns in the matrics a2; also the number
c       of rows in the matrix v2
c  ncols1 - the rank of the matrix a1; also the number of columns in
c       the matrices u1,v1
c  ncols2 - the rank of the matrix a2; also the number of columns in
c       the matrices u2,v2
c  eps - the precision of computationsc
c 
c                     output parameters:
c 
c  u3 - the n*ncols3 matrix in (3)
c  v3 - the (m1+m2)*ncols3 matrix in (3)
c  ncols3 - the numerical rank of the matrix  a3 in (3);
c        also the number of columns in the matrices u3, v3
c 
c                     work arrays:
c 
c  w - must be at least (m1+m2)*n +10 real *8 elements long
c  alpha - must be at least (ncols1+ncols2)*ncols1+10
c          real *8 elements long
c  beta -  must be at least (ncols1+ncols2)*ncols2+10
c          real *8 elements long
c  alpha0 - must be at least (ncols1+ncols2)**2+10
c          real *8 elements long
c  rnorms - must be at least ncols1+ncols2+2 real *8 elements long
c 
c        . . . combine the matrices u1, u2 horizontally
c              into a single matrix; it is assumed that the
c              matrices v1, v2 are orthogonal
c 
        do 1400 i=1,ncols1
        do 1200 j=1,n
        w(j,i)=u1(j,i)
 1200 continue
 1400 continue
c 
        do 2400 i=1,ncols2
        do 2200 j=1,n
        w(j,i+ncols1)=u2(j,i)
 2200 continue
 2400 continue
c 
c       gram-schmidt compress the resulting combined matrix
c 
        ifpivot=1
c 
        call grmpivot(w,n,ncols1+ncols2,rnorms,eps,ncols3,ifpivot)
c 
c       project the columns of matrices u1, u2 on the
c       (now) orthogonalized columns of the matrix w
c 
        do 2460 i=1,ncols3
        do 2440 j=1,ncols1
        call scapro(w(1,i),u1(1,j),n,prod)
        alpha(j,i)=prod
 2440 continue
 2460 continue
c 
        do 2500 i=1,ncols3
        do 2480 j=1,ncols2
        call scapro(w(1,i),u2(1,j),n,prod)
        beta(j,i)=prod
 2480 continue
 2500 continue
c 
c       copy the appropriate piece of array w into u3
c 
        call arrcopy(w,u3,ncols3*n)
c 
cccc         call prin2('alpha as constructed*',alpha,ncols3*ncols1)
cccc         call prin2('beta as constructed*',beta,ncols3*ncols2)
c 
c       now, construct the matrix v3
c 
        do 3600 i=1,m1
        do 3400 j=1,ncols3
        d=0
        do 3300 k=1,ncols1
        d=d+v1(i,k)*alpha(k,j)
 3300 continue
        v3(i,j)=d
 3400 continue
 3600 continue
c 
        do 4200 i=1,m2
        do 4000 j=1,ncols3
        d=0
        do 3800 k=1,ncols2
        d=d+v2(i,k)*beta(k,j)
 3800 continue
        v3(i+m1,j)=d
 4000 continue
 4200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine makeuort(uin,vin,n,m,ncols,uout,vout,
     1      alpha)
        implicit real *8 (a-h,o-z)
        save
        dimension uin(n,ncols),uout(n,ncols),vin(m,ncols),
     1      vout(m,ncols),alpha(ncols,ncols)
c 
c        this subroutine converts an expansion of the
c        form
c 
c        a=uin * vin ^*                                      (1)
c 
c        with orthogonal vin into an expansion of the form
c 
c        a=uout * vout ^*                                    (2)
c 
c        with orthogonal uout.
c 
c                        input parameters:
c 
c  uin - the matrix in (1).
c  vin - the matrix in (1).
c  n - the number of rows in the matrices uin, uout, a
c  m - the number of rown in the matrices vin, vout. also the number
c        of columns in the matrix a.
c  ncols - the number of columns in the matrices uin,uout,vin,vout.
c        conceptually, it is much smaller than n,m
c 
c                        output parameters:
c 
c  uout - the matrix in (2).
c  vout - the matrix in (2).
c 
c                        work arrays:
c 
c  alpha - must be at least ncols**2+2*ncols+10 real *8 elements long
c 
c        . . . apply the gram-schmidt procedure to the matrix uin
c 
        ifpivot=1
        call grmright(uin,n,ncols,uout,alpha,ncols2,
     1       alpha(1,ncols+1),eps,ifpivot)
c 
c       post-multiply the matrix v by alpha, obtaining vout
c 
        do 2000 i=1,ncols
        do 1800 j=1,m
        d=0
        do 1600 k=1,ncols
        d=d+vin(j,k)*alpha(k,i)
 1600 continue
        vout(j,i)=d
 1800 continue
 2000 continue
        return
c 
c 
c 
c 
        entry makevort(uin,vin,n,m,ncols,uout,vout,
     1      alpha)
c 
c        this subroutine converts an expansion of the
c        form
c 
c        a=uin * vin ^*                                      (1)
c 
c        with orthogonal uin into an expansion of the form
c 
c        a=uout * vout ^*                                    (2)
c 
c        with orthogonal vout.
c 
c                        input parameters:
c 
c  uin - the matrix in (1).
c  vin - the matrix in (1).
c  n - the number of rows in the matrices uin, uout, a
c  m - the number of rown in the matrices vin, vout. also the number
c        of columns in the matrix a.
c  ncols - the number of columns in the matrices uin,uout,vin,vout.
c        conceptually, it is much smaller than n,m
c 
c                        output parameters:
c 
c  uout - the matrix in (2).
c  vout - the matrix in (2).
c 
c                        work arrays:
c 
c  alpha - must be at least ncols**2+2*ncols+10 real *8 elements long
c 
c      . . .  apply gram-schmidt procedure to the matrix vin
c 
        ifpivot=1
        call grmright(vin,m,ncols,vout,alpha,ncols2,
     1       alpha(1,ncols+1),eps,ifpivot)
c 
c       post-multiply the matrix u by alpha,obtaining vout
c 
        do 3000 i=1,ncols
        do 2800 j=1,n
        d=0
        do 2600 k=1,ncols
        d=d+uin(j,k)*alpha(k,i)
 2600 continue
        uout(j,i)=d
 2800 continue
 3000 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine grmleft(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m),v(m,1),rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively, the columns
c       of v being orthonormal, and ncols being the rank of
c       the matrix a to the precision eps. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is considerably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least max(n,m) + 1
c        real *8 elements long.
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        v(i,j)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt process (with pivoting) to b
c 
         call grmpivot(v,m,n,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,n
c 
        call scapro2(a,j,v(1,i),n,m,prod)
        b(j,i)=prod
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine scapro2(a,k,x,n,m,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,1),x(1)
c 
c       evaluate the inner product of the k-th row of the
c       matrix a with the vector x
c 
        prod=0
        do 1200 i=1,m
        prod=prod+a(k,i)*x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine grmright(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m),v(m,1),rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively, the columns
c       of b being orthonormal, and ncols being the rank of
c       the matrix a to the precision eps. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is considerably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least max(n,m) + 1
c        real *8 elements long.
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt proces (with pivoting) to b
c 
         call grmpivot(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call scapro(a(1,j),b(1,i),n,prod)
        v(j,i)=prod
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine grmpivot(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=dsqrt(d)
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,m
c 
c       find the pivot
c 
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call scapro(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2785 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call scapro(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if(d .lt. thresh) return
        ncols=i
c 
        d=done/dsqrt(d)
c 
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call scapro(b(1,i),b(1,j),n,d)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine scapro(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine mulmatad(b,v,n,m,ncols,c)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,ncols),v(m,ncols),c(n,m)
c 
c       this subroutine multiplies the matrix b by the adjoint of v,
c       obtaining the matrix c
c 
c                  input parameters:
c 
c  b,v - the matrices to be multiplied
c  n,m,ncol - the dimensionalities of the matrices to be multiplied,
c       and of the result (see array declarations above)
c 
c                  output parameters:
c 
c  c - the result of multiplication
c 
c        multiply the matrices a and b ^*  obtaining the matrix c
c 
        do 2600 i=1,n
        do 2400 j=1,m
        d=0
        do 2200 k=1,ncols
        d=d+b(i,k)*v(j,k)
 2200 continue
        c(i,j)=d
 2400 continue
 2600 continue
        return
        end
  
