        real *8 a(100 000 00)
        integer *4 m, k, numit, ion, l
        integer *8 i, n
        character(17) file_name

c 
c        call prini(6,13)
        ion = 1
        call print_on_off(ion)


        m = 5
        k = 3
        n = k * 2**3
        call prinf("n = *", n, 1)
        call prinf("m = *", m, 1)   

c        write(*,*) "input m"
c        read(*,*) m
c        write(*,*) "input k"
c        read(*,*) k
c        write(*,*) "input l"
c        read(*,*) l
c        n = k*2**l
c        write(*,*) "n=",n
c        read(*,*) n
              

        m = 50
        k = 20
        n = k * 2**12
       
c        m = 25
c        k = 10
c        n = k * 2**10

c        m = 10
c        k = 5
c        n = k * 2**8

         m = 40
         k = 30
         n = k * 2**13


c        write(file_name, 2500) m,n
c 2500   FORMAT ("frt_m", I3.3, "_n", I7.7)
c        write(*, *) file_name
c
c        open(10, file = file_name)
c        do 30 i = 1, n*m
c           read(10, 3500) a(i)
c 30     continue
c 3500   FORMAT (E24.18)
c        close(10)

c
c       CALL RANN
c

        numit = 2
        call one_div_test(n,m,k,a,numit)


        stop
        end
c
c
c
 
c
c
c       INT8
c       n, lenw, i, near_convicts
c
c
        subroutine one_div_test(n,m,k,a,numit)
        implicit none
        save
        integer *4 m,k,numit,isuper,istat
        integer *8 n,lenw,i,near_convicts(80 000 00),nwords
        real *8 a(m,n),w(300 000 00),
     1      dists_convicts(80 000 00),inds_many(10 000 00),
     2      dists_many(10 000 00)
c
        character(23) file_idxs, file_dist

        call corrand_norm_4(n*m,a,w)
        isuper = 2
        istat = 2

        call prinf("first 10 points = *", i, 0)
        do 3104 i = 1, 10
           call prin2("*", a(1, i), m)
 3104   continue


        lenw = 300 000 00
        call get_memory_size(n,m,k,numit,isuper,nwords)
        call prinf("nwords = *", nwords, 1)
        lenw = nwords

        call prinf("nwords = *", nwords, 1)
        lenw = nwords
        if (lenw .GT. 300 000 00) then
           call prinl("error, lenw = *", lenw, 1)
           stop
        end if
        
        w(lenw+1) = 123

        call tree_test(n,m,a,numit,isuper,istat,
     1      near_convicts,dists_convicts,k,w,lenw,
     2      inds_many,dists_many)

c        call prinf("near_convicts = *", near_convicts, n*k)
c        call prin2("dists_convicts = *", dists_convicts, n*k)

        call prin2("avg true dist = *", w(20), 1)
        call prin2("avg susp dist = *", w(21), 1)
        call prin2("avg proportion = *", w(22), 1)

        call prin2("check = *", w(lenw+1), 1)
        call prinf("lenw = *", lenw, 1)

c
c       WRITE INDICES AND DISTANCES INTO FLAT FILES
c

        write(file_idxs, 5000) m, n, k
 5000   FORMAT ("i10d_m", I3.3, "_n", I7.7, "_k", I3.3)
        write (*, *) file_idxs

        open(10, file = file_idxs)
        do 20 i = 1, n*k
           write(10, 7000) near_convicts(i)
 20     continue
 7000   FORMAT (I12)
        close(10)

        write(file_dist, 5500) m, n, k
 5500   FORMAT ("d10d_m", I3.3, "_n", I7.7, "_k", I3.3)
        write (*, *) file_dist

        open(10, file = file_dist)
        do 30 i = 1, n*k
           write(10, 7500) dists_convicts(i)
 30     continue
 7500   FORMAT (E24.18)
        close(10)

        return
        end
c
c
c
