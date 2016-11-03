c
c  Vladimir Rokhlin (Yale University)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c 
c       This is the end of the debugging code, and the beginning of the
c       random transform code proper
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This file contains 6 user-callable subroutines: random_transf, 
c       random_transf_inverse, random_transf_c,random_transf_c_inverse,
c       random_transf_init. Following is a brief description of these 
c       subroutines.
c
c   random_transf_matr_c - applies rapidly a fairly random unitary
c       matrix to the user-supplied matrix a
c
c   random_transf_matr - applies rapidly a fairly random orthogonal 
c       matrix to the user-supplied matrix a
c   
c   random_transf - applies rapidly a fairly random orthogonal matrix to 
c       the user-supplied vector x, using the data in array w stored there 
c       by a preceding call to the subroutine random_transf_init (see)
c   
c   random_transf_inverse - applies rapidly the inverse of the operator 
c       applied by the subroutine random_transf_inverse 
c
c   random_transf_c - same as random_transf, but both the input and the 
c       output vectors are complex. Please note that the random matrix 
c       applied by this subroutine is REAL.
c
c   random_transf_c_inverse - applies rapidly the inverse of the operator 
c       applied by the subroutine random_transf_c (see)
c
c   random_transf_init - initialization subroutine. It prepares and stores 
c        in array w the data used by the subroutines random_transf,
c        random_transf_c (see).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c 
        subroutine random_transf_matr_c(if_inverse,a,n,w,nsteps,
     1      ifinit,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),w(1)
c
c       This subroutine applies a "fast" random unitary
c       transformation to the user-supplied complex matrix a(n,n). 
c       In other words, it replaces the input matrix a with the matrix
c
c                U^{-1} A U,
c
c       where A is the user-supplied matrix, and U is defined by the
c       subroutine random_transf_c (see)
c
c                     Input parameters:
c
c  if_inverse - integer parameter telling the subroutine whether it 
c       should apply the forward (if_inverse=0) or inverse (if_inverse=1)
c       transformation
c  a - the n \times n matrix to be transformed
c  n - the dimensionality of the matrix a
c  w - the array used by the subroutine random_transf (see). Please 
c       note that w is an input array only if ifinit=0; otherwise,
c       it is an output/work array
c  nsteps - the number of steps used by the randomizer; 4 is a good 
c       value
c  ifinit - integer parameter telling the subroutine whether a new 
c       random orthogonal transformation is to be created (ifinit=1),
c       or the old one (stored in the array w) is to be used (ifinit=0). 
c       Please note that for the first call to this subroutine (with a 
c       given n), ifinit must be set to 1.
c
c                      Output parameters:
c
c  a - the transformed matrix
c  w - the array containing data used by the subroutine random_transf 
c       (see) to apply a random orthogonal matrix to a vector. Please 
c       note that w is an output array only if ifinit=1; otherwise,
c       it is an input/work array. 
c  keep - the number of real *8 elements in array w that must be kept
c       unchanged between calls to this subroutine, if the said calls 
c       are to apply the same orthogonal matrix
c  lused - the total number of (real 88) elements of array w actually
c       used by this subroutine. Always equal to keep+2*n+4.
c
c
c
c
c       . . . if the user so requested, initialize the random 
c             matrix subroutine
c
        if(ifinit .ne. 0) call random_transf_init(nsteps,n,w,keep)
c
c       allocate memory
c
        iww=keep+2 
        lww=n+2 +n
        lww=n+2 
c
        iwww=iww+lww
        lwww=n+2 +n
        lwww=n+2 
c
        lused=iwww+lwww
c
c       multiply a from the left by the random matrix
c
        do 1400 i=1,n
c
        if(if_inverse .eq. 0) call random_transf_c(a(1,i),w(iww),w)
        if(if_inverse .ne. 0) 
     1      call random_transf_c_inverse(a(1,i),w(iww),w)
        call random_transf_copy(w(iww),a(1,i),n*2)
 1400 continue
c
c       multiply a from the right by the adjoint of the random
c       matrix
c
        do 2600 i=1,n
c
        do 2200 j=1,n
c
cccc        w(iww+j-1)=a(i,j)
        w(iww+j-1)=conjg(a(i,j))
 2200 continue
c
        if(if_inverse .eq. 0) call random_transf_c(w(iww),w(iwww),w)
        if(if_inverse .ne. 0) 
     1      call random_transf_c_inverse(w(iww),w(iwww),w)
c
        do 2400 j=1,n
c
cccc        a(i,j)=conjg(w(iwww+j-1))
        a(i,j)=conjg(w(iwww+j-1))
 2400 continue
 2600 continue
c
        return
        end
c
c
c
c
c 
        subroutine random_transf_matr(if_inverse,a,n,w,nsteps,
     1      ifinit,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),w(1)
c
c       This subroutine applies a "fast" random orthogonal
c       transformation to the user-supplied matrix a(n,n). In other
c       words, it replaces the input matrix a with the matrix
c
c                U^{-1} A U,
c
c       where A is the user-supplied matrix, and U is defined by the
c       subroutine random_transf (see)
c
c                     Input parameters:
c
c  if_inverse - integer parameter telling the subroutine whether it 
c       should apply the forward (if_inverse=0) or inverse (if_inverse=1)
c       transformation
c  a - the n \times n matrix to be transformed
c  n - the dimensionality of the matrix a
c  w - the array used by the subroutine random_transf (see). Please 
c       note that w is an input array only if ifinit=0; otherwise,
c       it is an output/work array
c  nsteps - the number of steps used by the randomizer; 4 is a good 
c       value
c  ifinit - integer parameter telling the subroutine whether a new 
c       random orthogonal transformation is to be created (ifinit=1),
c       or the old one (stored in the array w) is to be used (ifinit=0). 
c       Please note that for the first call to this subroutine (with a 
c       given n), ifinit must be set to 1.
c
c                      Output parameters:
c
c  a - the transformed matrix
c  w - the array containing data used by the subroutine random_transf 
c       (see) to apply a random orthogonal matrix to a vector. Please 
c       note that w is an output array only if ifinit=1; otherwise,
c       it is an input/work array. 
c  keep - the number of real *8 elements in array w that must be kept
c       unchanged between calls to this subroutine, if the said calls 
c       are to apply the same orthogonal matrix
c  lused - the total number of (real 88) elements of array w actually
c       used by this subroutine. Always equal to keep+2*n+4.
c
c
c       . . . if the user so requested, initialize the random 
c             matrix subroutine
c
        if(ifinit .ne. 0) call random_transf_init(nsteps,n,w,keep)
c
c       allocate memory
c
        iww=keep+2
        lww=n+2
c
        iwww=iww+lww
        lwww=n+2
c
        lused=iwww+lwww
c
c       multiply a from the left by the random matrix
c
        do 1400 i=1,n
c
        if(if_inverse .eq. 0) call random_transf(a(1,i),w(iww),w)
        if(if_inverse .ne. 0) 
     1      call random_transf_inverse(a(1,i),w(iww),w)
        call random_transf_copy(w(iww),a(1,i),n)
 1400 continue
c
c       multiply a from the right by the adjoint of the random
c       matrix
c
        do 2600 i=1,n
c
        do 2200 j=1,n
c
        w(iww+j-1)=a(i,j)
 2200 continue
c
        if(if_inverse .eq. 0) call random_transf(w(iww),w(iwww),w)
        if(if_inverse .ne. 0) 
     1      call random_transf_inverse(w(iww),w(iwww),w)
c
        do 2400 j=1,n
c
        a(i,j)=w(iwww+j-1)
 2400 continue
 2600 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),w(1)
c
c        This subroutine applies rapidly a fairly random 
c        orthogonal matrix to the user-specified real vector x, 
c        using the data in array w stored there by a preceding 
c        call to the subroutine random_transf_init (see)
c
c                Input parameters:
c
c  x - the vector of length n to which the matrix is to be applied
c  w - array containing all the data to be used by this subroutine.
c
c                Output parameters:
c
c  y - the results of applying the said matrix to the said vector
c
c        . . . if n is less than 21 -  act accordingly
c
       n=w(5)       
       if(n .gt. 5) goto 1100
       iu=11
c
       call random_transf_small(n,w(iu),x,y)
c
       return
 1100 continue
c
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
        k=w(6)
        iwsave=w(7)
        iixs2=w(8)
        ialphas2=w(9)
        ibetas2=w(10)
        iixs3=w(11)
c
        call random_transf0(nsteps,x,y,n,w(iww),w(ialbetas),w(iixs))
c
        call random_transf_usefft(n,y,w(iwsave),
     1      w(ialphas2),w(ibetas2),w(iixs3),w(iixs2),w(iww),k)
        return
        end
c
c
c
c
c 
        subroutine random_transf_init(nsteps,n,w,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
c
c        This is the initialization subroutine. It prepares and stores 
c        in array w the data used by the subroutines random_transf,
c        random_transf_c (see) to apply rapidly a fairly random 
c        orthogonal matrix to an arbitrary user-specified vector.
c
c                Input parameters:
c
c  nsteps - the degree of randomness of the operator to be applied
c  n - the dimensionality of the matrix to be applied     
c
c                Output parameters:
c
c  w - the first keep elements of w contain all the data to be used 
c        by the subroutine. Please note that the number of elements 
c        by this subroutine is also equal to keep. The length of 
c        this array should be at least 
c        2*nsteps*n + nsteps*n/ninire +8*n +n/4 + 2*n/ninire + 300.
c
c  keep - the number of elements in w that have been actually used 
c        by this subroutine; also the number that should not be 
c        changed between the call to this subroutine and the subsequent
c        calls to the subroutine random_transf (see).
c
c        . . . if n is less than 21, act accordingly
c
       call prinf('n=*',n,1)
c

        if(n .gt. 5) goto 1100
        iumatr=11
        lumatr=n*n+2
c       
        iw=iumatr+lumatr
        lw=n*(n+2)+4
c
        call random_transf_small_init(n,w(iumatr),w(iw))
        keep=n*n+12
        w(5)=n+0.1
        return
c
 1100 continue
c
c        . . . allocate memory 
c
        ninire=2
c
        ialbetas=20
        lalbetas=2*n*nsteps+10 
c
        iixs=ialbetas+lalbetas
        lixs=n*nsteps/ninire+10
c
c        find the size of the FFT to be used
c
        k=1
        do 1200 i=1,1000
c
        if(k .gt. n) goto 1400
c
        k=k*2
 1200 continue
 1400 continue
c
        k=k/4
c
        iwsave=iixs+lixs
        lwsave=k*4+30
c
        iixs2=iwsave+lwsave
        lixs2=n/ninire+5
c
        ialphas2=iixs2+lixs2
        lalphas2=n+4
c
        ibetas2=ialphas2+lalphas2
        lbetas2=n+4
c
        iixs3=ibetas2+lbetas2
        lixs3=n/ninire+5

        iww=iixs3+lixs3
        lww=2*n+n/4+20
c
        iww2=iww+lww
        lww2=2*n+n/4+20
c
        keep=iww2+lww2
c
        w(1)=ialbetas+0.1
        w(2)=iixs+0.1
        w(3)=nsteps+0.1
        w(4)=iww+0.1        
        w(5)=n+0.1
        w(6)=k+0.1
        w(7)=iwsave+0.1
        w(8)=iixs2+0.1
        w(9)=ialphas2+0.1
        w(10)=ibetas2+0.1
        w(11)=iixs3+0.1
        w(12)=iww2+0.1
c
        call random_transf_init0(nsteps,n,w(ialbetas),w(iixs))
c
        call DCFFTI(k,w(iWSAVE))
        call corrand_integer_knuth_4(n,w(iixs2))
        call random_transf_init00(n,w(ialphas2),w(ibetas2),w(iixs3))
c
        return
        end
c
c
c
c
        subroutine random_transf_inverse(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),w(1),x2(10000)
c
c        This subroutine applies rapidly a fairly random 
c        orthogonal matrix to the user-specified real vector x, 
c        using the data in array w stored there by a preceding 
c        call to the subroutine random_transf_init (see)
c
c        PLEASE NOTE THAT the transformation applied by this 
c        subroutine is the inverse of the transformation applied
c        by the subroutine random_transf (see)
c
c
c                Input parameters:
c
c  x - the vector of length n to which the matrix is to be applied
c  w - array containing all the data to be used by this subroutine.
c
c                Output parameters:
c
c  y - the result of applying the said matrix to the said vector
c
c        . . . if n is less than 21 -  act accordingly
c
c
       n=w(5)
       if(n .gt. 5) goto 1100
       iu=11
c
       call random_transf_small_inv(n,w(iu),x,y)
c
       return
 1100 continue
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
        k=w(6)
        iwsave=w(7)
        iixs2=w(8)
        ialphas2=w(9)
        ibetas2=w(10)
        iixs3=w(11)
        iww2=w(12)
c
c       undo the effects of the FFT block
c
        call random_transf_copy(x,w(iww2),n)
        call random_transf_usefft_inv(n,w(iww2),w(iwsave),
     1      w(ialphas2),w(ibetas2),w(iixs3),w(iixs2),w(iww),k)
c
c       undo the effects of the rotations
c
        call random_transf0_inv(nsteps,w(iww2),y,n,w(iww),
     1      w(ialbetas),w(iixs))
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_c(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1)
        dimension w(1)
c
c        This subroutine applies rapidly a fairly random 
c        orthogonal matrix to the user-specified complex vector x, 
c        using the data in array w stored there by a preceding 
c        call to the subroutine random_transf_init (see)
c
c                Input parameters:
c
c  x - the vector of length n to which the matrix is to be applied
c  w - array containing all the data to be used by this subroutine.
c
c                Output parameters:
c
c  y - the results of applying the said matrix to the said vector
c
c        . . . if n is less than 21 -  act accordingly
c
c
       n=w(5)       
       if(n .gt. 5) goto 1100
       iu=11
c
       call random_transf_small_c(n,w(iu),x,y)
c
       return
 1100 continue
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
c
        k=w(6)
        iwsave=w(7)
        iixs2=w(8)
        ialphas2=w(9)
        ibetas2=w(10)
        iixs3=w(11)
c
        call random_transf0_c(nsteps,x,y,n,w(iww),w(ialbetas),w(iixs))
c
        call random_transf_usefft_c(n,y,w(iwsave),
     1      w(ialphas2),w(ibetas2),w(iixs3),w(iixs2),w(iww),k)
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_c_inverse(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1)
        dimension w(1)
c
c        This subroutine. It applies rapidly a fairly random 
c        orthogonal matrix to the user-specified complex vector x, 
c        using the data in array w stored there by a preceding 
c        call to the subroutine random_transf_init (see). 
c
c        PLEASE NOTE THAT the transformation applied by this 
c        subroutine is the inverse of the transformation applied
c        by the subroutine random_transf_c (see)
c
c                Input parameters:
c
c  x - the vector of length n to which the matrix is to be applied
c  w - array containing all the data to be used by this subroutine.
c
c                Output parameters:
c
c  y - the results of applying the said matrix to the said vector
c
c
c        . . . if n is less than 21 -  act accordingly
c
c
       n=w(5)
       if(n .gt. 5) goto 1100
       iu=11
c
       call random_transf_small_c_inv(n,w(iu),x,y)
c
       return
 1100 continue
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
c
        k=w(6)
        iwsave=w(7)
        iixs2=w(8)
        ialphas2=w(9)
        ibetas2=w(10)
        iixs3=w(11)
        iww2=w(12)
c
        call random_transf_copy(x,w(iww2),n*2)
cccc        call random_transf_usefft_c_inv(n,x,w(iwsave),
        call random_transf_usefft_c_inv(n,w(iww2),w(iwsave),
     1      w(ialphas2),w(ibetas2),w(iixs3),w(iixs2),w(iww),k)
c
cccc        call random_transf0_c_inv(nsteps,x,y,n,w(iww),
        call random_transf0_c_inv(nsteps,w(iww2),y,n,w(iww),
     1      w(ialbetas),w(iixs))
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_usefft(n,x,wsave,alphas,betas,
     1      ixs,ixs2,w2,k)
        implicit real *8 (a-h,o-z)
        save
        dimension x(n),wsave(1),alphas(1),betas(1),
     1      ixs(1),w2(1),ixs2(1)
c
c       perform the extra permutation
c
        do 1200 i=1,n
c
        j=ixs2(i)
        w2(i)=x(j)
 1200 continue
c
c        apply the FFT
c
        call DCFFTF(k,w2,WSAVE)
c  
        done=1
        d=sqrt(done/k)
        do 1600 i=1,k*2
c
        w2(i)=w2(i)*d
 1600 continue
c
c        . . . apply the permutation & combination
c
        call random_transf00(w2,x,n,alphas,betas,ixs)
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_usefft_inv(n,x,wsave,alphas,betas,
     1      ixs,ixs2,w2,k)
        implicit real *8 (a-h,o-z)
        save
        dimension x(n),wsave(1),alphas(1),betas(1),
     1      ixs(1),w2(1),ixs2(1)
c
c       undo the permutation and combination
c
        call random_transf00_inv(x,w2,n,alphas,betas,ixs)
c
c       . . . and the FFT
c
        call DCFFTB(k,w2,WSAVE)
cccc        call DCFFTf(k,w2,WSAVE)
c  
        done=1
        d=sqrt(done/k)
        do 1200 i=1,k*2
c
        w2(i)=w2(i)*d
 1200 continue
c
c       undo the extra permutation
c
        do 1100 i=1,n
c
        j=ixs2(i)
        x(j)=w2(i)
 1100 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_usefft_c(n,x,wsave,alphas,betas,
     1      ixs,ixs2,w2,k)
        implicit real *8 (a-h,o-z)
        save
        dimension wsave(1),alphas(1),betas(1),ixs(1),ixs2(1)
        complex *16 x(n),w2(1)
c
c       perform the extra permutation
c
        do 1200 i=1,n
c
        j=ixs2(i)
        w2(i)=x(j)
 1200 continue
c
c        apply the FFT
c
        call DCFFTf(k,w2,WSAVE)
        call DCFFTf(k,w2(k+1),WSAVE)
c  
        done=1
        d=sqrt(done/k)
        do 1600 i=1,k *2
c
        w2(i)=w2(i)*d
 1600 continue
c
c        . . . apply the permutation & combination
c
        call random_transf00_c(w2,x,n,alphas,betas,ixs)
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_usefft_c_inv(n,x,wsave,alphas,betas,
     1      ixs,ixs2,w2,k)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(n),w2(1)
        dimension wsave(1),alphas(1),betas(1),ixs(1),ixs2(1)
c
c       undo the permutation and combination
c
        call random_transf00_c_inv(x,w2,n,alphas,betas,ixs)
c
c       . . . and the FFT
c
        call DCFFTB(k,w2,WSAVE)
        call DCFFTB(k,w2(k+1),WSAVE)
c  
        done=1
        d=sqrt(done/k)
        do 1200 i=1,k*2
c
        w2(i)=w2(i)*d
 1200 continue
c
c       undo the extra permutation
c
        do 1100 i=1,n
c
        j=ixs2(i)
        x(j)=w2(i)
 1100 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf0_inv(nsteps,x,y,n,w2,albetas,iixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),w2(1),albetas(n,2,1),iixs(n,1)
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=nsteps,1,-1
c
        call random_transf00_inv(w2,y,n,albetas(1,1,ijk),
     1      albetas(1,2,ijk),iixs(1,ijk) )
c
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf00_inv(x,y,n,alphas,betas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),alphas(1),betas(1),ixs(1)
c
c        implement 2 \times 2 matrices
c
        do 1600 i=1,n
        y(i)=x(i)
 1600 continue
c
        do 1800 i=n-1,1,-1
c
        d1=alphas(i)*y(i)-betas(i)*y(i+1)
        d2=betas(i)*y(i)+alphas(i)*y(i+1)
c
        y(i)=d1
        y(i+1)=d2
 1800 continue
c
c       implement the permutation
c
        do 2600 i=1,n
c
        j=ixs(i)
        x(j)=y(i)
 2600 continue
c
        do 2800 i=1,n
c
        y(i)=x(i)
 2800 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf0_c_inv(nsteps,x,y,n,w2,albetas,iixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1),w2(1)
        dimension albetas(n,2,1),iixs(n,1)
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=nsteps,1,-1
c
        call random_transf00_c_inv(w2,y,n,albetas(1,1,ijk),
     1      albetas(1,2,ijk),iixs(1,ijk) )
c
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf00_c_inv(x,y,n,alphas,betas,ixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1),d1,d2
        dimension alphas(1),betas(1),ixs(1)
c
c        implement 2 \times 2 matrices
c
        do 1600 i=1,n
        y(i)=x(i)
 1600 continue
c
        do 1800 i=n-1,1,-1
c
        d1=alphas(i)*y(i)-betas(i)*y(i+1)
        d2=betas(i)*y(i)+alphas(i)*y(i+1)
c
        y(i)=d1
        y(i+1)=d2
 1800 continue
c
c       implement the permutation
c
        do 2600 i=1,n
c
        j=ixs(i)
        x(j)=y(i)
 2600 continue
c
        do 2800 i=1,n
c
        y(i)=x(i)
 2800 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf0(nsteps,x,y,n,w2,albetas,iixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),w2(1),albetas(n,2,1),iixs(n,1)
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=1,nsteps
c
        call random_transf00(w2,y,n,albetas(1,1,ijk),
     1      albetas(1,2,ijk),iixs(1,ijk) )
c
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf00(x,y,n,alphas,betas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),alphas(1),betas(1),ixs(1)
c
c       implement the permutation
c
        do 1600 i=1,n
c
        j=ixs(i)
        y(i)=x(j)
 1600 continue
c
c        implement 2 \times 2 matrices
c
        do 1800 i=1,n-1
c
        d1=alphas(i)*y(i)+betas(i)*y(i+1)
        d2=-betas(i)*y(i)+alphas(i)*y(i+1)
c
        y(i)=d1
        y(i+1)=d2
 1800 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_init0(nsteps,n,albetas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension albetas(n,2,1),ixs(n,1)
c
        do 2000 ijk=1,nsteps
c
        call random_transf_init00(n,albetas(1,1,ijk),
     1      albetas(1,2,ijk),ixs(1,ijk) )
c
 2000 continue
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_init00(n,alphas,betas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension alphas(1),betas(1),ixs(1)
c
c       construct the random permutation
c
        ifrepeat=0
        call corrand_integer_knuth_4(n,ixs)
c
c       construct the random variables
c
        call corrand4(n/2,alphas)
        call corrand4(n/2,betas)
        call corrand4(n,alphas)
        call corrand4(n,betas)
c
c       construct the random 2 \times 2 transformations
c
        do 1400 i=1,n
c
        d=alphas(i)**2+betas(i)**2
        d=1/sqrt(d)
        alphas(i)=alphas(i)*d
        betas(i)=betas(i)*d
 1400 continue
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf0_c(nsteps,x,y,n,w2,albetas,iixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1),w2(1)
        dimension albetas(n,2,1),iixs(n,1)
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=1,nsteps
c
        call random_transf00_c(w2,y,n,albetas(1,1,ijk),
     1      albetas(1,2,ijk),iixs(1,ijk) )
c
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf00_c(x,y,n,alphas,betas,ixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1),d1,d2
        dimension alphas(1),betas(1),ixs(1)
c
c       implement the permutation
c
        do 1600 i=1,n
c
        j=ixs(i)
        y(i)=x(j)
 1600 continue
c
c        implement 2 \times 2 matrices
c
        do 1800 i=1,n-1
c
        d1=alphas(i)*y(i)+betas(i)*y(i+1)
        d2=-betas(i)*y(i)+alphas(i)*y(i+1)
c
        y(i)=d1
        y(i+1)=d2
 1800 continue
c
        return
        end
c
c
c
c
c 
        subroutine random_transf_small_init(n,umatr,w)
        implicit real *8 (a-h,o-z)
        save
        dimension umatr(n,n),w(n,1),ipivots(10000),rnorms(1000)
c
c       generate the random matrix 
c
        call corrand_norm_4(N*n,w,umatr)
        call corrand_norm_4(n*2,w(1,n+1),umatr)

        call prin2('w as created, w=*',w,n*(n+2) )
c
c       run the obtained random matrix through
c       the Gram-Schmidr
c
        eps=0
        call random_transf_piv(w,n,n+2,eps,rnorms,ipivots,ncols)
c
        call prin2('and rnorms=*',rnorms,ncols)
c
        do 2400 i=1,n
        do 2200 j=1,n
c
        umatr(j,i)=w(j,i)
 2200 continue
 2400 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine random_transf_piv(b,n,m,eps,rnorms,ipivots,ncols)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),ipivots(1)
        real *8 b(n,m),cd
c 
c        . . . initialize the array of pivots
c 
        do 1100 i=1,m
        ipivots(i)=i
 1100 continue
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
c 
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
        dtot=dtot+d
 1400 continue
c 
        thresh=dtot*eps**2
  
        thresh=sqrt(thresh)
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,n
c 
c       find the pivot
c 
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
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
c       orthogonalize the i-th column to all preceding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call random_transf_scap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call random_transf_scap(b(1,i),b(1,i),n,cd)
c 
        d=cd
        if(d .lt. thresh**2 ) return
c 
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        if(i .eq. m) goto 3400
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call random_transf_scap(b(1,i),b(1,j),n,cd)
c 
        cd=(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
c 
 4000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine random_transf_scap(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*(y(i))
 1200 continue
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_small_c(n,u,x,y)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n)
        complex *16 x(1),y(1),cd
c
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
c
        cd=cd+u(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
c
        return
c
c
c
c
        entry random_transf_small_c_inv(n,u,x,y)
c
        do 2400 i=1,n
        cd=0
        do 2200 j=1,n
c
        cd=cd+u(j,i)*x(j)
 2200 continue
        y(i)=cd
 2400 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine random_transf_small(n,u,x,y)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n)
        real *8 x(1),y(1),cd
c
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
c
        cd=cd+u(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
c
        return
c
c
c
c
        entry random_transf_small_inv(n,u,x,y)
c
        do 2400 i=1,n
        cd=0
        do 2200 j=1,n
c
        cd=cd+u(j,i)*x(j)
 2200 continue
        y(i)=cd
 2400 continue
c
        return
        end
c
c
c
c
c 
        subroutine random_transf_copy(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1)
c
        do 1400 i=1,n
c
        b(i)=a(i)
 1400 continue
c
        return
        end
