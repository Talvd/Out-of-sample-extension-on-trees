c
c  Vladimir Rokhlin, Andrei Osipov (Yale University)
c
c
c       INT8
c
c
        subroutine peter_heapsort_rea_part(ra,ia,n,k)
        implicit none
        integer *8 n,ia(n),k,i,ii,ii0,m,j,jj,i2,lone
        real *8 ra(1)
c
c        This subroutine uses the standard heapsort scheme to
c        find and sort the first k smallest elements of the 
c        user-supplied real *8 array ra. The sort is returned 
c        in the form of the integer array ia, such that the 
c        values
c
c        ra(ia(1)), ra(ia(2)), ra(ia(3)), . . .  ra(ia(k))        (1)
c
c        are in increasing order, and are the k smallest 
c        elements of array ra
c
c            Input parameters:
c
c  ra - the array to be sorted
c  n - the number of elements in array ra
c  k - the first k smallest elements in array ra will be found and
c        returned in increasing order
c
c            Output parameters:
c
c  ia - integer array containing the permutation that sorts ra 
c        (see (1) above). Please note that while only first k 
c        elements of ia are meaningful, its length has to be at
c        least n
c
c
        do 1200 i=1,n
c
        ia(i)=i
 1200 continue
c
        do 1600 ii=n,2,-2
c
        ii0=ii/2
        call peter_heapit_rea_decr(ra,ia,n,ii0)
 1600 continue
c        call prinl("in phrp, ia = *", ia, k)

c
c       . . . sort
c
cccc        do 2000 i=1,n-1
        do 2000 i=1,k
c
        m=n-i+1
c
        jj=ia(m)
        ia(m)=ia(1)
        ia(1)=jj
c
        m=m-1
        lone = 1
        call peter_heapit_rea_decr(ra,ia,m,lone)
 2000 continue
c        call prinl("again in phrp, ia = *", ia, k)
c
c       move the k smallest values to the beginning of the array
c
        do 2200 i=1,n/2
c
        i2=n-i+1
        j=ia(i2)
        ia(i2)=ia(i)
        ia(i)=j
 2200 continue

        return
        end
c
c
c       INT8
c
c
        subroutine peter_heapit_rea_decr(ra,ia,n,ii)
        implicit none
        integer *8 n,ia(n),ii,i,ijk,ison1,ison2,j0,j1,j2
        real *8 ra(1)
c
c       heapify
c
        i=ii
c        
        do 1400 ijk=1,100
c
        ison1=i*2
        if(ison1 .gt. n) return
c
        ison2=i*2+1
        if(ison2 .gt. n) then
c
            j0=ia(i)
            j1=ia(ison1)
cccc            if(ra(j0) .ge. ra(j1)) return
            if(ra(j0) .le. ra(j1)) return
c
            ia(i)=j1
            ia(ison1)=j0
            return
         endif
c
        j0=ia(i)
        j1=ia(ison1)
        j2=ia(ison2)
c
cccc        if( (ra(j0) .ge. ra(j1)) .and. (ra(j0) .ge. ra(j2))) return
        if( (ra(j0) .le. ra(j1)) .and. (ra(j0) .le. ra(j2))) return
c
cccc        if(ra(j1) .ge. ra(j2)) then 
        if(ra(j1) .lt. ra(j2)) then
           ia(i)=j1
            ia(ison1)=j0
            i=ison1
            goto 1400
        endif
c
            ia(i)=j2
            ia(ison2)=j0
            i=ison2
c
 1400 continue
            return
            end
