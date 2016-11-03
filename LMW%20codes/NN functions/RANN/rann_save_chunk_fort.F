c
c  Andrei Osipov (Yale University)
c
c
c
        subroutine save_points64(m,n,a)
        implicit none

        integer *4 m
        integer *8 n,ichunk,mn
        real *8 a(m,n),t1,t2,t3
        character(20) file_name

        ichunk = 125 000 000
        write(file_name, "(A)") "rann_pts"

        call cpu_time(t1)
        mn = m*n
        call save_real_by_chunks(file_name,mn,a,ichunk)
        call cpu_time(t2)
        t3=t2-t1
        call prin2("In save_points64, CPU time = *", t3, 1)

        return
        end

c
c
c
        subroutine load_points64(m,n,a)
        implicit none

        integer *4 m
        integer *8 n,ichunk,mn
        real *8 a(m,n),t1,t2,t3
        character(20) file_name

        ichunk = 125 000 000
        write(file_name, "(A)") "rann_pts"

        call cpu_time(t1)
        mn = m*n
        call load_real_by_chunks(file_name,mn,a,ichunk)
        call cpu_time(t2)
        t3=t2-t1
        call prin2("In load_points64, CPU time = *", t3, 1)

        return
        end


c
c
c
        subroutine save_distances64(m,n,a)
        implicit none

        integer *4 m
        integer *8 n,ichunk,mn
        real *8 a(m,n),t1,t2,t3
        character(20) file_name

        ichunk = 125 000 000
        write(file_name, "(A)") "rann_dst"

        call cpu_time(t1)
        mn = m*n
        call save_real_by_chunks(file_name,mn,a,ichunk)
        call cpu_time(t2)
        t3=t2-t1
        call prin2("In save_distances64, CPU time = *", t3, 1)

        return
        end

c
c
c
        subroutine load_distances64(m,n,a)
        implicit none

        integer *4 m
        integer *8 n,ichunk,mn
        real *8 a(m,n),t1,t2,t3
        character(20) file_name

        ichunk = 125 000 000
        write(file_name, "(A)") "rann_dst"

        call cpu_time(t1)
        mn = m*n
        call load_real_by_chunks(file_name,mn,a,ichunk)
        call cpu_time(t2)
        t3=t2-t1
        call prin2("In load_distances64, CPU time = *", t3, 1)

        return
        end

c
c
c
        subroutine save_indices64(m,n,a)
        implicit none

        integer *4 m
        integer *8 n,ichunk,mn
        real *8 a(m,n),t1,t2,t3
        character(20) file_name

        ichunk = 125 000 000
        write(file_name, "(A)") "rann_idx"

        call cpu_time(t1)
        mn = m*n
        call save_real_by_chunks(file_name,mn,a,ichunk)
        call cpu_time(t2)
        t3=t2-t1
        call prin2("In save_indices64, CPU time = *", t3, 1)

        return
        end

c
c
c
        subroutine load_indices64(m,n,a)
        implicit none

        integer *4 m
        integer *8 n,ichunk,mn
        real *8 a(m,n),t1,t2,t3
        character(20) file_name

        ichunk = 125 000 000
        write(file_name, "(A)") "rann_idx"

        call cpu_time(t1)
        mn = m*n
        call load_real_by_chunks(file_name,mn,a,ichunk)
        call cpu_time(t2)
        t3=t2-t1
        call prin2("In load_indices64, CPU time = *", t3, 1)

        return
        end


c
c
c       file_name: character of length 8
c       a(n): real array
c       ichunk: maximal size of a file (in real words)
c
c
        subroutine save_real_by_chunks(file_name,n,a,ichunk)
        implicit none

        character *1 file_name(1)
        character(13) file_name_chunk
        integer *8 n,ichunk,nrem,ncur,ipt
        integer *4 i,j
        real *8 a(n)

        nrem = n
        ipt = 1
        do 10 i = 0, 9999
           call prinf("i = *", i, 1)
           call prinl("nrem = *", nrem, 1)
           call prinl("ipt = *", ipt, 1)
           call prinl("ichunk = *", ichunk, 1)
           if (nrem .LE. 0) goto 20
          
           ncur = ichunk
           if (ncur .GT. nrem) ncur = nrem
           write(file_name_chunk, "(8A1,I4.4)") 
     1          (file_name(j),j=1,8), i
           call prina("file_name_chunk = *", file_name_chunk, 12)
           write(*,*) file_name_chunk

           call save_real_into_file64(file_name_chunk,ncur,a(ipt))

           call prinl("ncur = *", ncur, 1)

           nrem = nrem-ncur
           ipt = ipt+ncur
 10     continue
 20     continue

        return
        end
      
c
c
c
c
c
c
        subroutine save_real_into_file64(file_name,k,a)
        implicit none
        
        integer *8 k,k8,ksmall
        real *8 a(k)
        character(13) file_name

        call prinl("k = *", k, 1)
        k8 = k*8
        call prinl("k8 = *", k8, 1)

        ksmall = k
        if (ksmall .GT. 100) ksmall = 100
        call prin2("a = *", a, ksmall)

        open(10,file=file_name,status="replace",
     1          access="direct",
     2          form="unformatted",recl=k8)
        write(10,rec=1) a
        close(10)

        return
        end
   
c
c
c       file_name: character of length 8
c       a(n): real array
c       ichunk: maximal size of a file (in real words)
c
c
        subroutine load_real_by_chunks(file_name,n,a,ichunk)
        implicit none

        character *1 file_name(1)
        character(13) file_name_chunk
        integer *8 n,ichunk,nrem,ncur,ipt
        integer *4 i,j
        real *8 a(n)

        nrem = n
        ipt = 1
        do 10 i = 0, 9999
           call prinf("i = *", i, 1)
           call prinl("nrem = *", nrem, 1)
           call prinl("ipt = *", ipt, 1)
           call prinl("ichunk = *", ichunk, 1)
           if (nrem .LE. 0) goto 20
          
           ncur = ichunk
           if (ncur .GT. nrem) ncur = nrem
           write(file_name_chunk, "(8A1,I4.4)") 
     1          (file_name(j),j=1,8), i
           call prina("file_name_chunk = *", file_name_chunk, 12)
           write(*,*) file_name_chunk

           call load_real_from_file64(file_name_chunk,ncur,a(ipt))

           call prinl("ncur = *", ncur, 1)

           nrem = nrem-ncur
           ipt = ipt+ncur
 10     continue
 20     continue

        return
        end
      
c
c
c
c
c
c
        subroutine load_real_from_file64(file_name,k,a)
        implicit none
        
        integer *8 k,k8,ksmall
        real *8 a(k)
        character(13) file_name

        call prinl("k = *", k, 1)
        k8 = k*8
        call prinl("k8 = *", k8, 1)

        ksmall = k
        if (ksmall .GT. 100) ksmall = 100
        call prin2("a = *", a, ksmall)

        open(10,file=file_name,
     1          access="direct",
     2          form="unformatted",recl=k8)
        read(10,rec=1) a
        close(10)

        return
        end
