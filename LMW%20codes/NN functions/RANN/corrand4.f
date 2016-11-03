c
c Vladimir Rokhlin (Yale University)
c
c        implicit real *8 (a-h,o-z)
c        real *8 xs(40 000 000),ys(40 000 000),zs(40 000 000)
c        DIMENSION vals(1000 000),bins(1000 000),fs(1000 000),
c     1      ixs(1000 000)
c        DIMENSION vals2(1000 000),bins2(1000 000),fs2(1000 000)
cc 
c        call prini(6,13)
cC 
cC       SET ALL PARAMETERS
cC 
c         PRINT *, 'ENTER n'
c         READ *,n
c         CALL PRINF('n=*',n,1)
c
c        call corrand4(N,xs)
c
c        k=30
ccccc        k=200
ccccc        k=100
c
c        iw=31
c        call hystog0(N,xs,k,iw,bins,vals)
c
c        iw=32
c        call zquaplot(iw,xs,n/2,1,'area*')
cc
cc        construct the normal distributions
cc
c        call corrand_norm_4(N,xs,ys)
c        call corrand_norm_4(N,ys,zs)
c
c
c        iw=33
c        call hystog0(N,xs,k,iw,bins,vals)
cc
c        iw=34
c        call hystog0(N,ys,k,iw,bins2,vals2)
cc
cc        construct the distributions analytically and plot
cc        them things together with the experimental ones
cc
c        done=1
c        pi=atan(done)*4
cc
c        h=bins(2)-bins(1)
c        rint=0
c        do 2200 i=1,k
cc
c        fs(i)=exp(-bins(i)**2/2) /sqrt(pi*2)
c        rint=rint+fs(i)
c 2200 continue
cc
c        rint=rint*h
c        call prin2('and rint=*',rint,1)
cc
cc
c        iw2=53
c        call quagraph2(iw2,bins,vals,k-1,3,bins,fs,k-1,3,
c     1      'hystogram for vals vs. analytical result *')
cc
c        h=bins2(2)-bins2(1)
c        rint2=0
c        do 2400 i=1,k
cc
c        fs2(i)=exp(-bins2(i)**2/2) /sqrt(pi*2)
c        rint2=rint2+fs2(i)
c 2400 continue
cc
c        rint2=rint2*h
c        call prin2('and rint2=*',rint2,1)
cc
c        iw2=54
c        call quagraph2(iw2,bins2,vals2,k-1,3,bins2,fs2,k-1,3,
c     1      'hystogram for vals2 vs. analytical result *')
cc
c        call arbscap(xs,xs,n,d)
c
c        call prin2('(xs,xs)=*',d,1)
c
c        call arbscap(xs,ys,n,d2)
c
c        call prin2('(xs,ys)=*',d2,1)
cc
cc        test the integer number generator
cc 
c        call corrand_integer_knuth_4(n,ixs)
c        call prinf('after _integr, ixs=*',ixs,n)
c
c        call corrand_integer_knuth_4(n,ixs)
c        call prinf('after _integr, ixs=*',ixs,n)
cc
cccc        call bubba(N,ixs)
c
cccc        call prinf('after bubba, ixs=*',ixs,n)
c
c
c
c        do 3200 i=1,n
cc
c        fs2(i)=ixs(i)
c 3200 continue
c
c
c        iw=36
c        call zquaplot(iw,fs2,n/2,2,'area*')
cc
cc       test the gold-plated random number generator
cc
c        call corrand_norm_4(N,ys,zs)
c
c        call prin2('ys as created*',ys,20)
c
c        done=1
c        pi=atan(done)*4
c
c        dd=sqrt(2*done)
c
c        do 4200 i=1,n
cc
ccc        zs(i)=exp(-ys(i)**2/2)/sqrt(2*pi)
c
ccc        zs(i)=exp(-ys(i)**2/2)/sqrt(2*pi) *ys(i)
c
c        zs(i)=-exp(-ys(i)**2/2)/sqrt(2*pi) *ys(i) 
c
c        zs(i)=ys(i)/abs(zs(i))
c
c
c       
c        call qerrfun(ys(i)/dd,zs(i))
c
c        zs(i)=zs(i)/2
c
ccccc        zs(i)=1/zs(i)
c 4200 continue
c
c
c        call prin2('zs as created*',zs,20)
c
c        iw=37
c        k=300
c        k=100
c        call hystog0(N,zs,k,iw,bins2,vals2)
c
c
c        n=10
c        nn=100000
c        do 4400 i=1,nn
c
c        call corrand_one_integer_4(n,ixs(i))
c
c 4400 continue
c
c        call prinf('and ixs=*',ixs,nn)
c
c        do 4600 i=1,n
cc
c        vals(i)=0
c        xs(i)=i
c 4600 continue
cc
c        do 4800 i=1,nn
cc
c        do 4700 j=1,n
cc
c        if(ixs(i) .eq. j) vals(j)=vals(j)+1
c 4700 continue
c 4800 continue
c
c        iw=41        
c        call quagraph(iw,xs,vals,n,2,'hystogram*')
cc
c
c
c
c
c
c
c        stop
c        end
cC 
cC 
cC 
cC 
cC 
c        SUBROUTINE bubba(N,ixs)
c        IMPLICIT REAL *8 (A-H,O-Z)
c        save
c        DIMENSION ixs(1)
cc
c        do 1400 i=1,n
c        do 1200 j=1,n-1
c        if(ixs(j) .gt. ixs(j+1) ) then
c            jj=ixs(j)
c            ixs(j)=ixs(j+1)
c            ixs(j+1)=jj
c        endif
cc
c 1200 continue
c 1400 continue
c
c
c        return
c        end
cC 
cC 
cC 
cC 
cC 
c        SUBROUTINE hystog0(N,xs,k,iw,bins,vals)
c        IMPLICIT REAL *8 (A-H,O-Z)
c        save
c        DIMENSION xs(1),vals(1),bins(1)
cc
cc        determine the maximum value of the random variable
cc
c        xmax=-1.0d30
c        xmin=1.0d30
c        do 1200 i=1,n
cc
c        if(xmax .lt. xs(i)) xmax=xs(i)
c        if(xmin .gt. xs(i)) xmin=xs(i)
c 1200 continue
cc
c        call prin2('xmin=*',xmin,1)
c        call prin2('xmax=*',xmax,1)
cc
cc       . . . bin them things
cc
c        h=(xmax-xmin)/k
c
c        call prin2('h=*',h,1)
cc
c        done=1
c        pi=atan(done)*4
cc
c        rint=0
c        do 1400 i=1,k
cc
c        vals(i)=0
c        bins(i)=xmin+h/2+(i-1)*h
cc
c 1400 continue
c
c
c        call prin2('bins=*',bins,k)
cc
c        do 1600 i=1,n
cc
c        d=(xs(i)-xmin)/h
c        j=d+1
c        vals(j)=vals(j)+1
c 1600 continue
cc
c        do 1800 i=1,k
cc
c        vals(i)=vals(i)/n*k/(xmax-xmin)
c 1800 continue
cc
c        call prin2('vals=*',vals,k)
cc
cc       . . . plot
cc
c        call quagraph(iw,bins,vals,k-1,3,'hystogram*')
cc
c        return
c        end
cc
cc 
cc 
cc 
cc 
c        subroutine arbscap(x,y,n,d)
c        implicit real *8 (a-h,o-z)
c        save
c        dimension x(1),y(1)
cc 
c        d=0
c        do 2000 i=1,n
c        d=d+x(i)*y(i)
c 2000 continue
c        return
c        end
  
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
C        This file contains four user-callable subroutines: corrand4,
c        corrand_norm_4, corrand_integer_knuth_4, corrand_one_integer_4.
c        Following is a brief description of the said four subroutines.
c
c   corrand4 - returns to the user a pseudo-random vector distributed 
c        uniformly on the interval [0,1]. The algorithm used by this
c        subroutine is a pompous one: constructs nine pseudo-random 
c        sequences using nine different congruential generators, 
c        averages them, and uses the obvious transformation to bring it 
c        back to the uniform distribution. 
c
c   corrand_norm_4 - returns to the user two random vectors y1, y2, 
c        each of which is distributed normally with the 
c        distribution density
c
c        1/sqrt(2 * \pi) * e^(y^2/2)                         (1)
c   
c   corrand_integer_knuth_4 - for a user-specified n, returns to 
c        the user a random permutation of n integers 1,2,...,n
c
c   corrand_one_integer_4 - for a user-specified n, returns to 
c        the user a random integer on the interval [1,n]
c
C 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c 
c 
        subroutine corrand_one_integer_4(n,i)
        implicit real *8 (a-h,o-z)
        dimension tt(10)
        save
c       
c        This subroutine returns to the user a random integer
c        number i on the interval [0,n]
c
        call corrand4(1,tt(6))
c
        i=tt(6)*n
        i=i+1
c        
        return
        end
c
c
c
c 
c 
        subroutine corrand_integer_knuth_4(n,ixs)
        implicit real *8 (a-h,o-z)
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
        call corrand4(11,tt)
        call corrand4(11,tt)
        call corrand4(11,tt)
        do 1200 i=1,n
c
        ixs(i)=i
 1200 continue
c
        done=1
        do 1400 i=1,n-1
c
        call corrand4(1,tt)
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
        SUBROUTINE corrand_norm_4(N,Y1,y2)
        IMPLICIT REAL *8 (A-H,O-Z)
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
        call corrand4(N,Y1)
        call corrand4(N,Y2)
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

        SUBROUTINE corrand4(n,y)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        dimension js21(10),js22(10),js23(10),y(1),
     1      ias(10),ixks(10),ds(10)
ccc        data ias/7,17,31,37,43,53,61,67,79,91/
        data ias/97,17,31,37,43,53,61,67,79,91/
        data ifcalled/0/

        done=1
c
c       retrieve the prime numbers and conduct preliminary
c       randomization
c
        if(ifcalled .ne. 0) goto 1350
c
        call corrand_primes_4(js21,js22,js23)
c
        do 1200 i=1,10
c
        ixks(i)=js21(i)
 1200 continue
        ifcalled=1
c
        do 1300 i=1,100
        do 1250 j=1,9    
c
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))

 1250 continue
 1300 continue
c
 1350 continue
c
        do 2000 i=1,n
c
        do 1400 j=1,3
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))
c
        ds(j)=ixks(j)*done/(js23(j)-1)
 1400 continue
c
        d=ds(1)+ds(2)+ds(3)
        call corrand_comp_4(d,dd1)
c
        do 1600 j=4,6
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))
c
        ds(j)=ixks(j)*done/(js23(j)-1)
 1600 continue
c
        d=ds(4)+ds(5)+ds(6)
        call corrand_comp_4(d,dd2)

c
        do 1800 j=7,9
        call corrand_onestep_4(js23(j),js22(j),ias(j),ixks(j))
c
        ds(j)=ixks(j)*done/(js23(j)-1)
 1800 continue
c
        d=ds(7)+ds(8)+ds(9)
        call corrand_comp_4(d,dd3)
c
        d=dd1+dd2+dd3
        call corrand_comp_4(d,y(i))
c
 2000 continue
c
        return

        entry corrand4_restart
        ifcalled = 0
        return


        return
        end
c
c
c
c
c
        SUBROUTINE corrand_comp_4(x,rint)
        IMPLICIT REAL *8 (A-H,O-Z)
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
c
c
c
c
c
        SUBROUTINE corrand_onestep_4(m,ic,ia,ixk)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
c
        jj=ia*ixk+ic  
        j=jj/m
        ixk=jj-j*m
c
        return
        END
c
c
c
c
c
        SUBROUTINE corrand_primes_4(js21,js22,js23)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION is21(10),is22(10),is23(10),
     1      js21(1),js22(1),js23(1)
c
        data is21/9,19,21,55,61,69,105,111,121,129/
        data is22/3,17,27,33,57,87,105,113,117,123/
        data is23/15,21,27,37,61,69,135,147,157,159/
c
        i21=2**21
        i22=i21*2
        i23=i22*2
c
        do 1200 i=1,10
c
        js21(i)=i21-is21(i)
        js22(i)=i22-is22(i)
        js23(i)=i23-is23(i)
 1200 continue
c
        return
        end
