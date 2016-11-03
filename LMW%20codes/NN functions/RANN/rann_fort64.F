c
c  Vladimir Rokhlin, Andrei Osipov (Yale University)
c
c       INT8
c       n,lenw,near_convicts,lused,keep
c
c
        subroutine tree_test(n,m,a,numit,isuper,istat,
     1      near_convicts,dists_convicts,k,w,lenw)
        implicit none
        save
        integer *4 m,numit,k,ier,iwrite,i,isuper,istat,
     1             itimes(15),it1,it2
        integer *8 n,lenw,near_convicts(k,1),lused,keep,ntest
        real *8 a(m,n),w(1),dists_convicts(k,1),
     2      t1,t2,t3
c
c       For each of the n input points in R^m, select k
c       "approximately nearest" points in array a
c
c       
c        call print_on_off(1)
c
c       ANDREI IMPORTANT: OTHERWISE HARD TO COMPARE!!!
c       in the real version, one can delete this line
c
        call corrand4_restart

        call prinf("n = *", n, 1)
        call prinf("m = *", m, 1)         

        call prin2("a = *", a, 12)

        call prina("first 10 points = *", 18)
        do 3104 i = 1, 10
           call prin2("*", a(1, i), m)
 3104   continue

c
c       ANDREI: replace with cpu_time
c        t1=clotatim()
        call cpu_time(t1)

        call prin2("t1 = *", t1, 1)
c
        it1=time()

        call all_iter(ier,n,m,a,k,numit,isuper,
     1      near_convicts,dists_convicts,w,lenw,lused,keep,
     2      itimes)

        it2=time()
        itimes(11)=it2-it1


c
c       ANDREI: replace with cpu_time
c        t2=clotatim()
        call cpu_time(t2)

        call prin2('after all_iter, w=*',w,3)
        call prinl('and lused/1000=*',lused/1000,1)
        call prinl('and keep/1000=*',keep/1000,1)

cccc        stop
c
c       . . . print the report
c
        itimes(12)=0
        if (istat .EQ. 0) goto 3111

        iwrite=13
        
        it1=time()
        call get_stats(n,m,a,numit,
     1      near_convicts,dists_convicts,k,w,lused,ntest,
     2      iwrite,istat)
        it2=time()
        itimes(12)=it2-it1

 3111   continue

        call cpu_time(t3)

        w(20) = w(1)
        w(21) = w(2)
        w(22) = w(3)
        w(23) = t2-t1
        w(24) = t3-t2

        w(30) = 1D0 * itimes(11) + 1D-10
        w(31) = 1D0 * itimes(12) + 1D-10
        w(32) = 1D0 * itimes(4) + 1D-10
        w(33) = 1D0 * itimes(3) + 1D-10
        w(34) = 1D0 * itimes(6) + 1D-10
        w(35) = 1D0 * itimes(7) + 1D-10
        w(36) = 1D0 * itimes(8) + 1D-10
        w(37) = 1D0 * itimes(9) + 1D-10
        w(38) = 1D0 * itimes(10) + 1D-10
        w(39) = 1D0 * itimes(5) + 1D-10
        w(40) = 1D0 * itimes(1) + 1D-10
        w(41) = 1D0 * itimes(2) + 1D-10

 3400 format(10x)

c
 5200 format(6x,'The CPU time in seconds: ',33x,E11.5)
c
c         write(6,5200) t2-t1
c         write(iwrite,5200) t2-t1

        call prinf("TIMES in sec = *", itimes, 10)

        call prinf("Total all_iter0=*",itimes(4),1)
        call prinf("   -> one_iter=*", itimes(3),1)
        call prinf("   -> -> centralize=*",itimes(6),1)
        call prinf("   -> -> rotate=*",itimes(7),1)
        call prinf("   -> -> one_iter0=*",itimes(8),1)
        call prinf("   -> -> -> struct_bld=*",itimes(9),1)
        call prinf("   -> -> -> boxes=*",itimes(10),1)
        call prinf("   -> merge=*",itimes(5),1)
        call prinf("Total second_search=*",itimes(1),1)
        call prinf("   -> copy to dst,susp=*",itimes(2),1)

        call prinf("TIME for all_iter=*",itimes(11),1)
        call prinf("TIME for get_stats=*",itimes(12),1)

        return
        end
c
c
c
c
c
        subroutine get_stats_lookup(rat_ave,ntest,
     1      dnum_belongs_ave,smalls,rat0s,n,m,k,
     2      a,inds_many,dists_many,cpu)
        implicit real *8 (a-h,o-z)
        dimension a(m,n),inds_many(k,1),dists_many(k,1)
c
        ic=0
c
        nwrong=0
        sum_true_all=0
        sum_appr_all=0
c       ANDREI
        sum_true_all2 = 0
        sum_appr_all2 = 0
        smalls=0
        rat0s=0
        do 1800 jkj=1,ntest
c
        ic=ic+1
        if(ic .eq. 100) then
            call prinf('in get_stats_lookup, jkj=*',jkj,1)
            ic=0
        endif
c
        call get_stats_onepnt(a(1,jkj),a,m,n,inds_many(1,jkj),
     1      dists_many(1,jkj),k,rat,prop_belong,num_belong,
     2      nsmall,sum_true,sum_appr,rat0, sum_true2, sum_appr2)
c
        sum_true_all=sum_true_all+sum_true
        sum_appr_all=sum_appr_all+sum_appr
        sum_true_all2 = sum_true_all2 + sum_true2
        sum_appr_all2 = sum_appr_all2 + sum_appr2
        rat0s=rat0s+rat0

cccc        call prinf('and nsmall=*',nsmall,1)
        if(nsmall .gt. 1) nwrong=nwrong+1
        smalls=smalls+nsmall-1

cccc        call prin2('and prop_belong=*',prop_belong,1)
cccc        call prinf('and num_belong=*',num_belong,1)

cccc        call prin2('and rat=*',rat,1)

cccc        rats(jkj)=rat
cccc        num_belongs(jkj)=num_belong
cccc        prop_belongs(jkj)=prop_belong-1

        rat_ave=rat_ave+rat
        prop_belongs_ave=prop_belongs_ave+prop_belong
        dnum_belongs_ave=dnum_belongs_ave+num_belong


 1800 continue

        rat_ave=rat_ave/ntest

        call prin2('and rat_ave=*',rat_ave,1)
        call prin2('and rat_ave-1=*',rat_ave-1,1)


        dnum_belongs_ave=dnum_belongs_ave/ntest
        

        call prin2('and dnum_belongs_ave=*',dnum_belongs_ave,1)

        call prin2('smalls=*',smalls,1)

        rat0s=rat0s/ntest
c
        iwrite=13

 3400 format(10x)

c
 5200 format(6x,'The CPU time in seconds: ',33x,E11.5)
c

 6200 format(25x, 'Report on the Table Look-Up')


         write(6,6200)
         write(iwrite,6200)

         write(6,3400)
         write(iwrite,3400)

         write(6,3400)
         write(iwrite,3400)

 6250 format(6x,'Average proportion of points that are closer'
     1      ' to one',/,6x,'of the first NTEST points than the',
     2      ' nearest neighbor',/,6x,
     3      'as determined by the subroutine: ', 25x,E11.5)
c
        write(6,6250) smalls/n/ntest
        write(iwrite,6250) smalls/n/ntest

         write(6,3400)
         write(iwrite,3400)

c
 6300 format(6x,'The average ratio RAT between the TOTAL',
     1      ' sum of',/,6x, 'distances',
     2      ' of allegedly nearest k points to the same',/,6x,
     3      'sum of distances ','to true nearest k points,',
     4      ' taken',/,6x,'over the first NTEST',
     5      ' points:',30x,E11.5) 
c
         write(6,6300) sum_appr_all/sum_true_all
         write(iwrite,6300) sum_appr_all/sum_true_all

 6311 format(6x,'Appr Dist is ',45x,E11.5)
         write(6,6311) sum_appr_all/ntest
         write(iwrite,6311) sum_appr_all/ntest

c
 6400 format(6x,'And RAT-1 is ',45x,E11.5)
c
         write(6,6400) sum_appr_all/sum_true_all-1
         write(iwrite,6400) sum_appr_all/sum_true_all-1
            
c
         write(6,3400)
         write(iwrite,3400)

 6500 format(6x,'The average RAT_AVER of the individual',
     1      ' ratios between',/,6x,'the distance to the putative ',
     2      'nearest points',
     1      ' to the ',/,6x,
     2      'distances ','to their respective true nearest points,',
     3      /,6x,'taken over',' the first NTEST',
     3      ' points:',24x,E11.5) 
c
         write(6,6500) rat0s
         write(iwrite,6500) rat0s
c
 6550 format(6x,'And RAT_AVER-1 is ',40x,E11.5)
c
         write(6,6550) rat0s-1
         write(iwrite,6550) rat0s-1
            
         write(6,3400)
         write(iwrite,3400)


 6600 format(6x,'The number of nearest points determined ',
     2      'incorrectly',/,6x,'for the first NTEST points:',36x,i6)
c
         write(6,6600) nwrong
         write(iwrite,6600) nwrong
c
 6800 format(6x,'And the proportion of nearest points determined ',
     2      /,6x,'incorrectly:',46x,e11.5)
c
         done=1
         write(6,6800) nwrong*done/ntest
         write(iwrite,6800) nwrong*done/ntest

c
         write(6,3400)
         write(iwrite,3400)
c
c
 7000 format(6x,'On the average, among the k allegedly nearest points ',
     1     /,6x,'found by the algorithm, the number that belong there:',
     2      5     x,E11.5)
c
         write(6,7000) dnum_belongs_ave
         write(iwrite,7000) dnum_belongs_ave
c
         write(6,3400)
         write(iwrite,3400)
c
         write(6,5200) cpu
         write(iwrite,5200) cpu
c

        call prin2("SQ rat_ave = AVG(ratio) =*", rat_ave, 1)
        call prin2("SQ ratio(AVG) =*", sum_appr_all2/sum_true_all2, 1)
        call prin2("avg(true_k) = *", sum_true_all2/(ntest*k), 1)
        call prin2("avg(appr_k) = *", sum_appr_all2/(ntest*k), 1)
 
        return
        end
c
c
c
c
c
c
        subroutine get_stats_onepnt(x,a,m,n,inds,dists,k,
     1      rat,prop_belong,num_belong,nsmall,
     2      sum_true,sum_appr,rat0,sum_true2, sum_appr2)
        implicit real *8 (a-h,o-z)
        dimension a(m,n),dists(1),inds(1),x(1),dds(10 000 000),
     1      ia(10 000 000),dds2(10000),iaa(10000),ib(10000)
c
c        analyze the results and print a summary for the single
c        point x
c
c
c        . . . construct the nearest k neighbors via brute force
c
cccc        call prin2('in _onepnt, dists=*',dists,k)

        nsmall=0
        do 1600 i=1,n
c
        call dist_compute(x,a(1,i),m,d)
c
        dds(i)=d
        if(d .lt. dists(2)) nsmall=nsmall+1


cccc        call prin2('dists(1)=*',dists(1),k)
cccc        call prin2('and d=*',d,1)
 1600 continue


cccc        stop
c
c
c       sort the distances from the point x to all 
c       points in a, and take the first k elements
c
        call peter_heapsort_rea_part(dds,ia,n,k)
c
        d1=0
        d2=0
        sum_true=0
        sum_appr=0
        do 1800 j=1,k
c
        dds2(j)=dds(ia(j))
c
        d1=d1+dds2(j)
        d2=d2+dists(j)
c
c ANDREI DIST
        sum_true=sum_true+sqrt(dds2(j))
        sum_true2 = sum_true + dds2(j)
c ANDREI DIST
        sum_appr=sum_appr+sqrt(dists(j))
        sum_appr2 = sum_appr + dists(j)
 1800 continue
c
        rat0=dists(2)/dds2(2)
        
c       ANDREI
        sum_true2 = d1
        sum_appr2 = d2


cccc        call prin2('and dists=*',dists,6)
cccc        call prin2('and dds2=*',dds2,6)
cccc        call prin2('and rat0=*',rat0,1)



c
cccc        call prin2('and dds2=*',dds2,k)
cccc        call prinf('and ia=*',ia,k)


        rat=d2/d1
cccc        call prin2('and ratio between sums of squares of distances*',
cccc     1      rat,1)

cccc        call prin2('and that same ratio -1*',rat-1,1)

c
c       count the number of discrepancies between the "brute force"
c       and approximate nearest k neighbors for this point
c
        do 2000 j=1,k
c
        dds(j)=dds2(j)
        dds(j+k)=dists(j)
        ia(k+j)=inds(j)
 2000 continue
c
c       . . . sort them things
c
        call peter_heapsort_inte(ia,iaa,2*k)
c
        do 2200 j=1,2*k
c
        ib(j)=ia(iaa(j))
 2200 continue
c
c       squeeze out the doubles
c
        nleft=1
        do 2400 j=2,k*2
        if(ib(j-1) .ne. ib(j)) nleft=nleft+1
 2400 continue
c
        done=1
        num_belong=2*k-nleft
        prop_belong=num_belong*done/k        

cccc        call prinf('and nleft=*',nleft,1)
cccc        call prin2('and prop_belong=*',prop_belong,1)
cccc        call prinf('and num_belong=*',num_belong,1)
c
c
        return
        end
c
c
c
c
c
        subroutine onepoint_all(x,w,a,
     1      near_convicts,dists_convicts,inds,dists,w2)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),inds(1),dists(1),a(1),
     1      near_convicts(1),dists_convicts(1),w(1),w2(1)
c
c        given a structure built by a preceding call to the 
c        subroutine all_iter (see) and a user-supplied point
c        x in R^m, it finds the k points (approximately) 
c        nearest to x. 
c
c            
c                        Input parameters:
c
c  x - the user-specified point for which the nearest neighbors
c        are to be found
c  w - the array of data (hopefully) produced by a preceding call
c        to the subroutine all_iter (see)
c  near_convicts - integer array dimensioned (k,n) produced via a 
c        preceding call to the subroutine all_iter
c     Explanation: the element near_convicts(j,i) contains the 
c        sequence number in array a of the j-th neighbor of
c        the i-th point in a;
c  dists_convicts - real array dimensioned (k,n) produced via a 
c        preceding call to the subroutine all_iter: it contains
c        squares of distances corresponding to the points specified 
c        in array near_convicts 
c     Explanation: the element dists_convicts(j,i) contains the 
c        square of the distance between the point a(j) and 
c        the element whose number is in near_convicts
c
c                        Output parameters:
c
c  inds - integer array containing the sequence numbers in array
c        a of the (allegedly) nearest points in a to the user-supplied
c        point xx
c  dists - dists(j) is the square of the distance of the point 
c        a(inds(j)) to the user-supplied point x
c
c                        Work array:
c
c  w2 - need not be very long ??????????????????????????????????????????
c
c       conduct the numit iterations 
c
        numit=w(1)
        iv=w(2)
        iindsv=w(3)
c
        k=w(4)
        m=w(5)
        levmax=w(6)
c
        iinds2=1
        linds2=numit*k+10 
c
        if(linds2 .lt. k*k+10) linds2=k*k+10
c
        idists2=iinds2+linds2
        ldists2=linds2
c
        ijjs=idists2+ldists2
        ljjs=(levmax+1)*k+10
        ljjs=ljjs*2
c
        isuspects=ijjs+ljjs
        lsuspects=(levmax+1)*k+10
        lsuspects=lsuspects*2
c
        idds=isuspects+lsuspects
        ldds=(levmax+1)*k+10
        ldds=ldds*2
c
        call all_sets_suspects(x,w(iv),w(iindsv),a,numit,k,
     1      inds,dists,w2(iinds2),w2(idists2),w2(ijjs),
     2      w2(isuspects),w2(idds))
c
c       conduct the supercharging for a single point
c
        call onepoint_super(x,inds,dists,n,m,a,k,
     1      near_convicts,dists_convicts,w2(iinds2),w2(idists2))
c
        return
        end
c
c
c
c
c
        subroutine onepoint_super(x,inds,dists,n,m,a,k,
     1      near_convicts,dists_convicts,inds2,dists2)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),inds(1),dists(1),a(m,n),near_convicts(k,1),
     1      dists_convicts(k,1),inds2(1),dists2(1)
c
c       scan the lists of neighbors of all neighbors of x,
c       constructing one big happy neighbor list
c
        do 1200 i=1,k
c
        inds2(i)=inds(i)
 1200 continue
c
        ii=k
        do 2000 i=1,k
        ik=inds(i)
        do 1600 j=1,k
c
        ii=ii+1
        inds2(ii)=near_convicts(j,ik)
 1600 continue
 2000 continue
c
c       evaluate the distances for all suspects
c
        do 2200 i=1,k*k
c
        call dist_compute(x,a(1,inds2(i)),m,dists2(i))
c
 2200 continue
c
        k22=k*k
        call peter_heapsort_rea_up(dists2,inds2,k22,k+2,kk)
c
c       extract from arrays iarr, rarr the non-repeating 
c       elements 
c
        jjj=1
        inds(1)=inds2(k22)
        dists(1)=dists2(k22)
c
        do 2400 j=k22-1,1,-1
c
        if(inds2(j) .eq. inds2(j+1)) goto 2400
c
        jjj=jjj+1
        inds(jjj)=inds2(j)
        dists(jjj)=dists2(j)
c
        if(jjj .eq. k) goto 2400
 2400 continue
c
        return
        end
c
c
c
c
c
        subroutine all_sets_suspects(x,v,indsv,a,numit,k,
     1      inds,dists,inds2,dists2,jjs,suspects,dds)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),v(1),inds(1),dists(1),a(1),indsv(1),
     1      inds2(1),dists2(1),jjs(1),suspects(1),dds(1)
c
c       conduct the numit iterations one after another
c
        ii=indsv(1)
        call one_set_suspects(x,v(ii),a,inds2,dists2,jjs,suspects,dds)
c
        do 1600 i=2,numit
c
c       conduct the iteration
c
        ii=indsv(i)
        call one_set_suspects(x,v(ii),a,inds,dists,jjs,suspects,dds)
c
c       merge the arrays inds, dists with the arrays inds2,dists2
c
        j0=(i-1)*k
        do 1400 j=1,k
c
        inds2(j0+j)=inds(j)
        dists2(j0+j)=dists(j)
 1400 continue
c
 1600 continue
c
c       compress them things
c
c
c       . . . sort
c
        nmk=numit*k
c
        call peter_heapsort_rea_up(dists2,inds2,nmk,k+2,kk)
c
c       extract from arrays iarr, rarr the non-repeating 
c       elements 
c
        jjj=1
        inds(1)=inds2(nmk)
        dists(1)=dists2(nmk)
c
        do 2200 j=nmk-1,1,-1
c
        if(inds2(j) .eq. inds2(j+1)) goto 2200
c
        jjj=jjj+1
        inds(jjj)=inds2(j)
        dists(jjj)=dists2(j)
c
        if(jjj .eq. k) goto 2400
 2200 continue
 2400 continue
c
cccc        call prinf('and final inds=*',inds,k)
cccc        call prin2('and final dists=*',dists,k)
c
        return
        end
c
c
c
c
c
        subroutine one_set_suspects(x,v,a,inds,dists,
     1      jjs,suspects,dds)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),v(1),inds(1),dists(1),a(1),
     1      jjs(1),suspects(1),dds(1)
c
c       retrieve discrete information from array v
c
        call retr_arrs(n,m,k,ilevels,levmax,
     1      iboxes,nbout,icenters,iw_transf,keepw,
     2      iipts,ids,v)
c
c       find the suspects for this values of x, v
c
        call one_set_suspects0(x,n,m,a,k,v(ilevels),levmax,
     1      v(iboxes),v(icenters),v(iw_transf),
     2      v(iipts),v(ids),inds,dists,jjs,suspects,dds)
c
        call index_aid(v(iipts),inds,k)
c
        return
        end
c
c
c
c
c
        subroutine one_set_suspects0(x,n,m,a,k,levels,levmax,
     1      boxes,centers,w_transf,ipts,ds,
     2      inds,dists,jjs,suspects,dds)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),iconts(60),
     1      levels(1),ipts(1),suspects(1)
        dimension a(m,n),centers(1),jjs(1),dds(1),x(1),
     1      inds(1),dists(1)
c       
c       find the box where the user-supplied point lives
c
        call inwhich_box(x,m,ds,centers,levmax,
     1      w_transf,ibox,boxes,suspects,dds)
c
c       find the list of neighbors boxes, and scan the said list, 
c       finding the nearest neighbors in it
c
c
c       . . .  extract the list of this box's neighboring boxes
c
        call one_box_lists(ibox,levels,levmax,iconts)
c
        nbs=levmax+1
c
c       . . . construct the list of all points in all nbs 
c             neighboring boxes
c
        nnj=0
        do 2400 j=1,nbs
c
        jbox=iconts(j)
        j1=boxes(6,jbox)
        nj1=boxes(7,jbox)
        do 2200 jj=1,nj1
c
        nnj=nnj+1
        suspects(nnj)=j1+jj-1
 2200 continue
 2400 continue
c
c       scan the array suspects of points in neighboring boxes, 
c       finding the nearest k points
c
        do 2600 j=1,nnj
c
        jjj=suspects(j)
c
        jj=ipts(jjj)
c
        call dist_compute(x,a(1,jj),m,d)
c
        jjs(j)=jjj
        dds(j)=d
c
 2600 continue
c
c       collect the nearest k suspects
c
        call peter_heapsort_rea_up(dds,jjs,nnj,k+2,kk)
c
        do 2800 i=1,k
        inds(i)=jjs(nnj-i+1)
        dists(i)=dds(nnj-i+1)
 2800 continue
c
        return
        end
c
c
c
c
c
        subroutine inwhich_box(x,m,ds,centers,levmax,
     1      w,ibox,boxes,y,y2)
        implicit real *8 (a-h,o-z)
        save
        DIMENSION x(m),centers(1),y(1),ds(1),y2(1)
        integer *4 boxes(10,1)
c
c       apply to the vector x the (once) random transformation 
c
        do 1200 i=1,m
c
        y2(i)=x(i)-ds(i)
 1200 continue
c
        call random_transf(y2,y,w)
c       
c       find the box where the point x lives
c
        ibox=1
c
        do 1600 i=1,60
c
        if(y(i) .lt. centers(ibox) ) then
            ibox=boxes(4,ibox)
            goto 1400
        endif
c
            ibox=boxes(5,ibox)
 1400 continue
c
        if(i .eq. levmax) return
 1600 continue
c
        return
        end
c
c
c
c
c
        subroutine index_aid(ipts,inds,k)
        implicit real *8 (a-h,o-z)
        save
        integer *4 ipts(1),inds(1)
c
        do 1200 i=1,k
c
        j=ipts(inds(i))
        inds(i)=j
 1200 continue
c
        return
        end
c
c
c
c       INT8
c       n,near_convicts,lenw,lused,keep etc.
c
c    
        subroutine all_iter(ier,n,m,a,k,numit,isuper,
     1      near_convicts,dists_convicts,w,lenw,lused,keep,
     2      itimes)
        implicit none
        save
        integer *4 m,k,numit,jer,ier,isuper,it1,it2,itimes(15)
        integer *8 n,near_convicts(k,1),lenw,lused,keep,
     1             levmax,nsteps,nboxes_max,keepw_max,lenstore,
     2             iindsv,lindsv,iv,lv,iipts,lipts,iib,lib,
     3             irs,lrs,iib2,lib2,irs2,lrs2,iia,lia,
     4             inear_suspects,lnear_suspects,idists_near,
     5             ldists_near,iw,lenw2,iiarr,liarr,irarr,lrarr,
     6             iiarr2,liarr2,irarr2,lrarr2,lused2
        real *8 a(m,n),w(1),dists_convicts(k,1),done,d

c
c        Given a collection of n points in R^m, this subroutine
c        finds for each of the said n points k (approximately)
c        nearest neighbors. It uses a randomized search, and 
c        returns two arrays: the element near_convicts(j,i) contains
c        the sequence number in array a of the j-the neighbor of
c        the i-th point in a; element dists_convicts(j,i) contains
c        the square of the distance between the point a(i) and 
c        the element whose number in in near_convicts(j,i).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     IMPORTANT! ATTENZIONE! ACHTUNG! UNLESS YOU UNDERSTAND THE 
c     OPERATION OF THIS SUBROUTINE WELL, DO NOT DEVIATE FROM THE 
c     RECOMMENDED VALUES OF THE PARAMETERS k, numit. OTHERWISE, 
c     SEVERE DEGRADATION OF PERFORMANCE IS PROBABLE. EVERY STEP 
c     ASIDE WILL BE CONSIDERED AND ACTED UPON AS AN ATTEMPTED ESCAPE!
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      Input parameters:
c
c  a - the points whose neighbors are to be found
c  m - the dimensionality of the space where the n ponts live
c  n - the number of said points
c  k - the algorithm will keep subdividing the group as long
c        as its size is greater than nb_max. 
c     RECOMMENDED VALUE: 30
c  numit - the number of tests the algorithm will perform. 
c     RECOMMENDED VALUE: 10-20
c  lenw - the length of the work array w in real *8 words
c  
c                       Output parameters:
c
c  near_convicts - array dimensioned (k,n). 
c     Explanation: the element near_convicts(j,i) contains the 
c        sequence number in array a of the j-th neighbor of
c        the i-th point in a;
c  dists_convicts - squares of distances corresponding to the 
c        points specified in array near_convicts 
c     Explanation: the element dists_convicts(j,i) contains the 
c        square of the distance between the point a(j) and 
c        the element whose number is in near_convicts
c
c                       Work array:
c
c  w - must be large
c
c        . . . determine the amount of memory to be stored
c


        done=1
        d=log(n*done/k)/log(2*done)
        levmax=d
        call prinl('levmax=*',levmax,1)
        if (levmax .GT. m) levmax = m

        nsteps=6
        nboxes_max=2**levmax*2 
        keepw_max=3*m*nsteps+10*m

        call prinl('keepw_max=*',keepw_max,1)

c
        lenstore=nboxes_max*10+levmax+nboxes_max+keepw_max
     1      +nboxes_max*(levmax+1)+n+m+1000
c
c        lenstore=lenstore*numit
c
        iindsv=100
        lindsv=numit+100
c
        iv=iindsv+lindsv
        lv=lenstore+1
c
        keep=iv+lv 
c
c       ANDREI: MEMORY
c
        call prinl("iindsv = *", iindsv, 1)
        call prinl("iv = *", iv, 1)
        call prinl("iv+lv = keep = *", keep, 1)

c
c        . . . allocate memory
c
        iipts=1 + lenstore+1100
c
c       ANDREI: MEMORY FIX
c       iipts might get less than iv+lv
c       (if numit > 900)
c       in this case w(iipts) overrides w(iv) !!!!
c       solution: make it later
c
c
        iipts = keep + 1000
        lipts=n+2
c
c
c
        call prinl("iipts = *", iipts, 1)

c
        iib=iipts+lipts
        lib=2*k+2
c
        irs=iib+lib
        lrs=2*k+2
c
        iib2=irs+lrs
        lib2=2*k+2
c
        irs2=iib2+lib2
        lrs2=2*k+2
c
        iia=irs2+lrs2
        lia=2*k+2
c
        inear_suspects=iia+lia  
        lnear_suspects=n*k+100
c
        idists_near=inear_suspects+lnear_suspects
        ldists_near=n*k+100
c
        iw=idists_near+ldists_near
        lenw2=lenw-iw-1
c
        call prinl('and lenw2/1000=*',lenw2/1000,1)
        call prinl('while lenw/1000=*',lenw/1000,1)

c AAA
c        stop
c
        if(lenw2 .lt. 2000) then
            ier=32
            call prinl('bombing from all_iter with lenw2=*',
     1          lenw2,1)
            stop
        endif

c
        it1=time()
        call all_iter0(jer,n,m,a,k,numit,nsteps,
     1      near_convicts,w(inear_suspects),w(idists_near),
     2      dists_convicts,w(iw),lenw2,lused2,
     3      w(iipts),w(iia),w(iib),w(irs),w(irs2),w(iib2),w(iv),
     4      w(iindsv),itimes)
        it2=time()
        itimes(4)=it2-it1

c
c        call prin2("at iv = *", w(iv), lv)

c        stop


        w(1)=numit+0.1
        w(2)=iv+0.1
        w(3)=iindsv+0.1
        w(4)=k+0.1
        w(5)=m+0.1
        w(6)=levmax+0.1
        lused=lused2+iw
c
        if(jer .ne. 0) then
            call prinf('bombing; jer from all_iter0 is*',jer,1)
            stop
        endif
c
c        . . . now, supercharge
c
        iiarr=1 + lenstore+1100
c
c       ANDREI: MEMORY FIX
c       SEE iipts ABOVE
c
        iiarr = keep + 1000
        liarr=n+2
        if ( (k*(k+1)) .GT. n ) liarr = k*(k+1)+2
c
        irarr=iiarr+liarr
        lrarr=n+2             
        if ( (k*(k+1)) .GT. n ) lrarr = k*(k+1)+2
c
        iiarr2=irarr+lrarr
        liarr2=n+2
        if ( (k*(k+1)) .GT. n ) liarr2 = k*(k+1)+2
c
        iia=iiarr2+liarr2
        lia=n+2
        if ( (k*(k+1)) .GT. n ) lia = k*(k+1)+2
c
        irarr2=iia+lia
        lrarr2=n+2
        if ( (k*(k+1)) .GT. n ) lrarr2 = k*(k+1)+2

c
        inear_suspects=irarr2+lrarr2
        lnear_suspects=n*k+2
c
        idists_near=inear_suspects+lnear_suspects
        ldists_near=n*k+2

        itimes(1)=0
        if (isuper .EQ. 0) goto 1222
c
c       ANDREI: COMMENT THIS SUPERCHAGING
c

        it1=time()

c
c       OLD SUPERCHARGING ("liberal")
c        
        if (isuper .EQ. 1) then
           call second_search_old(a,m,n,k,near_convicts,
     1         dists_convicts,
     2         w(iiarr),w(irarr),w(iiarr2),w(iia),w(irarr2))
        end if

c
c       NEW SUPERCHARGING ("conservative")
c
        if (isuper .EQ. 2) then
           call second_search(a,m,n,k,near_convicts,
     1         dists_convicts,
     2         w(iiarr),w(irarr),w(iiarr2),w(iia),w(irarr2),
     3         w(inear_suspects),w(idists_near),itimes )
        end if


        it2=time()
        it2=it2-it1
        call prinf("SUPERCHARGING time = *",it2,1)
        itimes(1)=it2


 1222   continue

c
c        call prinl("near_convicts = *", near_convicts, n*k)
c        call prin2("dists_convicts = *", dists_convicts, n*k)

        return
        end
c
c
c       INT8
c
c
c
        subroutine all_iter0(ier,n,m,a,k,numit,nsteps,
     1      near_convicts,near_suspects,dists_near,
     2      dists_convicts,w,lenw2,lused,
     3      ipts,ia,ib,rs,rs2,ib2,v,indsv,itimes)
        implicit none
        save
        integer *4 ier,m,k,numit,ijk,jer,itimes(15),it1,it2
        integer *8 n,nsteps,near_convicts(k,1),ipts(1),indsv(1),
     1             ia(1),ib(1),ib2(1),near_suspects(k,1),
     2             ibb,lbb,lused,lenw2,iw,lenw3,ii,i,j,
     3             keepv,i3,i2,iv,lused3,l2k
        real *8 a(m,n),w(1),rs(1),
     1      rs2(1),dists_near(k,1),dists_convicts(k,1),v(1)
c
        itimes(3)=0
        itimes(5)=0
        itimes(6)=0
        itimes(7)=0
        itimes(8)=0
        itimes(9)=0
        itimes(10)=0

        indsv(1)=1
        do 3000 ijk=1,numit
c
c       . . . perform the iteration 
c
        call prinf('ijk=*',ijk,1)

        ibb=1
        lbb=m*n+100
c
        lused=ibb+lbb+2000
c
c        call prinf("lenw2 = *", lenw2, 1)
c        call prinf("lused = *", lused, 1)
c        stop


        if(lused .gt. lenw2) then
            ier=16
            call prinl('bombing from all_iter0 with lused-lenw2=*',
     1          lused-lenw2,1)
            stop
        endif
c
        iw=ibb+lbb
        lenw3=lenw2-iw-1
        call prinl("lenw3 = *", lenw3, 1)
c        stop
c
        call prinl('and lenw3/1000=*',lenw3/1000,1)
        call prinl('while ibb/1000=*',ibb/1000,1)
        call prinl('and lbb/1000=*',lbb/1000,1)

cccc        stop

        iv=indsv(ijk)
        call prinl("iv = *", iv, 1)

        it1=time()
        call one_iter(jer,n,m,a,k,nsteps,
     1      near_suspects,w(ibb),ipts,w(iw),lenw3,lused3,
     2      dists_near,v(iv),keepv,rs,itimes)
        it2=time()
        itimes(3)=itimes(3)+it2-it1
c
        indsv(ijk+1)=indsv(ijk)+keepv+2

        call prinl('and indsv=*',indsv,ijk+1)

c       ANDREI STOP
c        stop

cccc        stop
c
        call prinl("BEFORE lused = *", lused/1000, 1)
        lused=iw+lused3
        call prinl("AFTER lused = *", lused/1000, 1)
c        stop
c
        if(ijk .ne. 1) goto 1600
c
        do 1400 i=1,n
c
        ii=ipts(i)
        do 1200 j=1,k
c
        near_convicts(j,ii)=ipts(near_suspects(j,i))
        dists_convicts(j,ii)=dists_near(j,i)
 1200 continue
 1400 continue
c ANDREI: PRINT IPTS
c        call prinl("IPTS = *", ipts, n)
c        call prinl("near convicts = *", near_convicts, n*k)
c        call prin2("dists_convicts = *", dists_convicts, n*k)
c
c        stop
        goto 3000
c
 1600 continue
c
c       this is not the first iteration; one column after another,
c       unscramble the array near_suspects, and merge it with
c       the array near_convicts
c
c        stop

        it1=time()

        do 2600 i=1,n
c
        ii=ipts(i)
        do 1800 j=1,k
c
        ib(j)=ipts(near_suspects(j,i))
        ib(k+j)=near_convicts(j,ii)
        rs(j)=dists_near(j,i)
        rs(j+k)=dists_convicts(j,ii)
 1800 continue
c
c        sort the array rs
c
        l2k = 2*k
        call peter_heapsort_rea(rs,ia,l2k)
c
        do 2000 j=1,k*2
c
        rs2(j)=rs(ia(j))
        ib2(j)=ib(ia(j))
 2000   continue
c
c       . . . compress the arrays ib2,rs2, removing duplicates
c
        i2=1
        i3=1
        ib(1)=ib2(1)
        rs(1)=rs2(1)
c
        do 2200 j=1,k*2
c
        i2=i2+1
        if(i2 .gt. k*2) goto 2400
        if(ib2(i2) .eq. ib2(i2-1) ) goto 2200
c
        i3=i3+1
        ib(i3)=ib2(i2)       
        rs(i3)=rs2(i2)       
 2200 continue
c
 2400 continue
c
        do 2500 j=1,k
c
        near_convicts(j,ii)=ib(j)
        dists_convicts(j,ii)=rs(j)
 2500 continue
c
 2600 continue
        it2=time()
        itimes(5)=itimes(5)+it2-it1

 3000 continue

        call prinl('and indsv=*',indsv,numit+1)
c        call prinl("near convicts = *", near_convicts, n*k)
c        call prin2("dists_convicts = *", dists_convicts, n*k)

        return
        end
c
c
c       INT8
c
c
        subroutine one_iter(ier,n,m,a,k,nsteps,
     1      near_suspects,b,ipts,w,lenw3,lused,dists_near,v,
     2      keepv,ds,itimes)
        implicit none
        save
        integer *4 ier,m,k,nsteps_short,it1,it2,itimes(15)
        integer *8 n,nsteps,near_suspects(k,1),ipts(1),lenw3,lused,
     1             keepv,levmax,ndmax,i,j,keep_appr,iww,
     2             iboxes,lboxes,iinboxes,linboxes,iiconts,liconts,
     3             iiwork,liwork,isuspects,lsuspects,irs,lrs,
     4             iinbox,linbox,iia,lia,inear_nums,lnear_nums,
     5             icenters,lcenters,keep,nbout,levels(100)
        real *8 a(m,n),w(1),
     1      b(m,n),dists_near(k,1),ds(1),v(1),done,d
c       
c       perform one iteration
c
        ier=0
c
        done=1
        d=log(n*done/k)/log(2*done)
        levmax=d
        call prinl('levmax=*',levmax,1)
        if (levmax .GT. m) levmax = m
c
        ndmax=levmax+1
c
c       . . . find the center of the collection of points in R^m
c
cccc        do 1900 ijk=1,2

        it1=time()

        do 1800 i=1,m
c
        d=0
        do 1200 j=1,n
c
        d=d+a(i,j)
 1200 continue
c
        d=d/n
        ds(i)=d
        do 1400 j=1,n
c
        b(i,j)=a(i,j)-d

 1400 continue
 1800 continue

        it2=time()
        itimes(6)=itimes(6)+it2-it1
c
 1900 continue
c        call prin2("ds = *", ds, m)
c        call prin2("first 10 b = *", b, m*10)

c
c        . . . apply a random orthogonal transformation to
c              each column of the matrix a
c
c        call random_transf(b(1,i),w(iww),w)
c
c       ANDREI: MEMORY FIX
c       THIS LINE ABOVE IS NOT NEEDED AND DANGEROUS!
c       BECAUSE FOR THE FIRST TIME iww = 0
c       and we write something into w(iww)
c
        keep_appr=3*m*nsteps+10*m

        call prinl('keep_appr=*',keep_appr,1)
        call prinl("lenw3 = *", lenw3, 1)
c
c       ANDREI: MEMORY just to be sure
c
        if((keep_appr+100+m) .gt. lenw3) then
            ier=8
            call prinl('bombing from one_iter with lenw3=*',
     1          lenw3,1)
            stop
        endif

c
c       note, keep and nsteps are integer *8
c       the procedure expects to get integer *4
c
        it1=time()
        keep = 0
        nsteps_short = nsteps
        call random_transf_init(nsteps_short,m,w,keep)

        call prinl('keep=*',keep,1)

c        stop
c
        iww=keep+100
        call prinl("iww = *", iww, 1)
        do 2400 i=1,n
c
c       do w(iww)(1,..,m) = transform( b(1,1),...,b(1,m) ) (based on w)
c
        call random_transf(b(1,i),w(iww),w)
c
        do 2200 j=1,m
c
        b(j,i)=w(iww+j-1)

c ANDREI: UNDO THE RANDOM TRANSFORMATION
c        b(j,i) = a(j,i)
c        w(iww+j-1) = a(j,i)
c ANDREI

 2200 continue
 2400 continue
        it2=time()
        itimes(7)=itimes(7)+it2-it1

        call prin2("after random, b = *", b, 10*m)

c
c       allocate memory 
c
        iboxes=keep+100
c       ANDREI: MEMORY MULTIPLY BY 2, NOT ENOUGH FOR BOXES
        lboxes=2*(2**levmax)*10+100
c
        iinboxes=iboxes+lboxes
        linboxes=n+10
c
        iiconts=iinboxes+linboxes
        liconts=1000
c
        iiwork=iiconts+liconts
        liwork=n*6+200
c
        isuspects=iiwork+liwork
        lsuspects=n+10
c
        irs=isuspects+lsuspects
        lrs=n+10
c
        iinbox=irs+lrs
        linbox=k*10+10
c
        iia=iinbox+linbox
        lia=n+10
c
        inear_nums=iia+lia
        lnear_nums=n+10
c
        icenters=inear_nums+lnear_nums
        lcenters=4*2**levmax+100
c
        lused=2**levmax*(10+ndmax)+7*n+10*k+2000+4*2**levmax+100

c
c       ANDREI: MEMORY
c
        call prinl("lused = *", lused, 1)
        call prinl("icen+lcen=*", icenters+lcenters,1)

        call prinl('and lused/1000=*',lused/1000,1)
c
c       ANDREI: MEMORY
c       WHAT ABOUT icenters+lcenters?
c
        if ((lused .gt. lenw3) .or. 
     1      ((icenters+lcenters) .gt. lenw3)) then
            ier=8
            call prinl('bombing from one_iter with lused-lenw3=*',
     1          lused-lenw3,1)
            stop
        endif
c
        it1=time()
        call one_iter0(n,m,b,k,levmax,
     1      w(iinboxes),near_suspects,
     2      w(iboxes),w(iiconts),levels,ipts,w(iiwork),
     3      w(isuspects),w(iinbox),w(irs),w(iia),
     4      dists_near,w(inear_nums),w(icenters),nbout,
     5      itimes)

        it2=time()
        it2=it2-it1
        itimes(8)=itimes(8)+it2
        call prinf("ONE ITERATION time = *",it2,1)        

c        call prinl("after oi0, near_susp = *", near_suspects, 50)
c
c       store the produced data for future use
c
c        call prini(0, 0)
c        call store_arrs(n,m,k,levels,levmax,
c     1     w(iboxes),nbout,w(icenters),w,keep,
c     2      ipts,ds,v,keepv)
c        call prini(6, 13)
c
        call prinl('keepv=*',keepv,1)

        return
        end
c
c
c       USED IN one_set_suspects ONLY
c
c
        subroutine retr_arrs(n,m,k,ilevels,levmax,
     1      iboxes,nbout,icenters,iw_transf,keepw,
     2      iipts,ids,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)


c
c       retrieve the scalar data and the addresses
c
        n=w(1)
        m=w(2)
        k=w(3)
        levmax=w(4)
        nbout=w(5)
c
        iboxes=w(11)
        ilevels=w(12)
        icenters=w(13)
        iw_transf=w(14)
        iipts=w(16)
        ids=w(17)

        return
        end
c
c
c       USED ONLY IF QUERIES FOR NEW POINTS ARE ENABLED
c
c
        subroutine store_arrs(n,m,k,levels,levmax,
     1      boxes,nbout,centers,w_transf,keepw,
     2      ipts,ds,w,keep)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),levels(1),ipts(1)
        dimension centers(1),ds(1),w(1)
c       
c       store in array w all of the data needed by the 
c       subroutine one_set_suspects
c
c        . . . allocate memory
c
        call prinf('nbout=*',nbout,1)

        nbout_appr=2*2**levmax

        call prinf('nbout_appr=*',nbout_appr,1)
        call prinf('keepw=*',keepw,1)
c
        iboxes=100
        lboxes=nbout*10+10
c
        ilevels=iboxes+lboxes
        llevels=levmax+2
c
        icenters=ilevels+llevels
        lcenters=nbout+2
c
        iw=icenters+lcenters
        lw=keepw+6
c
        iipts=iw+lw
        lipts=n+2
c
        ids=iipts+lipts
        lds=m+2
c
        keep=ids+lds
c
c       . . . copy the arrays
c
        call peter_copy(boxes,w(iboxes),lboxes)
        call peter_copy(levels,w(ilevels),llevels)
        call peter_copy(centers,w(icenters),lcenters)
        call peter_copy(w_transf,w(iw),lw)
        call peter_copy(ipts,w(iipts),lipts)
        call peter_copy(ds,w(ids),lds)
c
c       store the scalar data
c
        w(1)=n+0.1
        w(2)=m+0.1
        w(3)=k+0.1
        w(4)=levmax+0.1
        w(5)=nbout+0.1
c
        w(11)=iboxes+0.1
        w(12)=ilevels+0.1
        w(13)=icenters+0.1
        w(14)=iw+0.1
        w(15)=itrack_coord+0.1
        w(16)=iipts+0.1
        w(17)=ids+0.1
c
        return
        end
c
c
c       INT8
c
c
        subroutine one_iter0(n,m,a,k,levmax,
     1      inboxes,near_suspects,
     2      boxes,iconts,levels,ipts,work,suspects,inbox,
     3      rs,ia,dists_near,near_nums,centers,nbout,itimes)
        implicit none
        save
        integer *4 m,k,itimes(15),it1,it2
        integer *8 n,boxes(10,1),iconts(1),ipts(1),suspects(1),inbox(1),
     1             inboxes(1),ia(1),near_suspects(k,1),near_nums(1),
     2             levmax,i,ibox,nbout,lev,i1,n1,nbs,nnj,jbox,j1,nj1,
     3             iii,jjj,ipt,jj,ipoint,j,np,levels(1),itmp
        real *8 a(m,n),rs(1),dists_near(k,1),centers(1),work(1),d
c       
c       build the structure
c
        it1=time()
        call struct_bld(n,m,a,
     1      boxes,levels,levmax,nbout,ipts,work,centers)
        it2=time()
        itimes(9)=itimes(9)+it2-it1

        call prinl("after sb, n = *", n, 1)
        call prinf("after sb, m = *", m, 1)
        call prin2("after sb, a(1) = *", a, m)
c        call prinl("after sb, boxes = *", boxes, nbout*10)
        call prinl("after sb, levmax = *", levmax, 1)
        call prinl("after sb, levels = *", levels, levmax)
        call prinl("after sb, nbout = *", nbout, 1)
c        call prinl("after sb, ipts = *", ipts, n)
        call prin2("after sb, centers = *", centers, 6)

        itmp=0
c
        do 1400 i=1,n
c
        near_nums(i)=0
 1400 continue
c
c       one box after another, find the list of neighbors, and
c       scan the said list, finding the nearest neighbors in it
c

        it1=time()
        call prinl("nbout = *", nbout, 1)
        do 3000 ibox=1,nbout

        itmp = itmp + 1
        if (itmp .EQ. 10000) then
           itmp = 0
           call prinl("ibox = *", ibox, 1)
        end if
c
        lev=boxes(2,ibox)
        if(lev .ne. levmax) goto 3000
c
        i1=boxes(6,ibox)
        n1=boxes(7,ibox)

c        call prinf("lev = *", lev, 1)
c        call prinf("i1 = *", i1, 1)
c        call prinf("n1 = *", n1, 1)
c        stop

c
        do 2000 i=1,n1
c
        inbox(i)=i1+i-1
 2000 continue
c
c       . . .  extract the list of this box's neighbors
c
        call one_box_lists(ibox,levels,levmax,iconts)
c
        nbs=levmax+1
c
c       construct the list of all points in all nbs 
c       neighboring boxes
c
        nnj=0
        do 2400 j=1,nbs
c
        jbox=iconts(j)
        j1=boxes(6,jbox)
        nj1=boxes(7,jbox)
c
        do 2200 jj=1,nj1
c
        nnj=nnj+1
        suspects(nnj)=j1+jj-1
 2200 continue
 2400 continue
c        call prinl("nbs = *", nbs, 1)
c        call prinl("suspects = *", suspects, nbs)
c        stop
c
        iii=0
c
c       process points in the box ibox one after another
c
        do 2800 i=1,n1
c
        ipoint=inbox(i)
        ipt=ipts(ipoint)
c
c       for the point number ipt, scan the array suspects of 
c       points in neighboring boxes, finding the nearest k points
c
        do 2600 j=1,nnj
c
        jjj=suspects(j)
c
c
        if(jjj .ge. ipoint) goto 2600
c
c
        jj=ipts(jjj)
c
        call dist_compute(a(1,ipt),a(1,jj),m,d)
c        call prinl("jjj = *", jjj, 1)
c        call prinl("jj = *", jj, 1)
c        call prinl("ipt = *", ipt, 1)
c        call prin2("d = *", d, 1)
c
        np=near_nums(ipoint)
c        call prinl("np = *", np, 1)
        if( (np .eq. k) .and. (d .ge. dists_near(1,ipoint)) ) goto 2500
        call store_one(ipoint,jjj,d,near_suspects,dists_near,
     1      near_nums,k)

c        call prinl("ipoint = *", ipoint, 1)
c        call prinl("ns(ipoint) = *", near_suspects(1,ipoint), k)
c
 2500 continue
c
        np=near_nums(jjj)
c        call prinl("np = *", np, 1)
        if( (np .eq. k) .and. (d .ge. dists_near(1,jjj)) ) goto 2600
        
c        call prin2("d = *", d, 1)
c        call prinl("near_suspects = *", near_suspects(1,jjj), k)
c        call prinl("near_nums = *", near_nums, 20)
c        call prin2("dists_near = *", dists_near(1,jjj), k) 

        call store_one(jjj,ipoint,d,near_suspects,dists_near,
     1      near_nums,k)

c        call prinl("jjj = *", jjj, 1)
c        call prinl("ns(jjj) = *", near_suspects(1,jjj), k)
c        if ((ipoint .EQ. 5) .AND. (jjj .EQ. 2)) then
c           stop
c        end if
c
 2600 continue

c
c ANDREI: PRINT POINT AND ITS SUSPECTS
c       call prinl("ipt = *", ipt, 1)
c       call prinl("suspects = *", suspects, nnj)
 
 2800 continue
c
 3000 continue

        it2=time()
        itimes(10)=itimes(10)+it2-it1

c
c        call prinl("boxes = *", boxes, nbout*10)
c        call prinl("near_suspects = *", near_suspects, k*n)
c        stop
        return
        end
c
c
c       INT8
c
c
        subroutine store_one(ipoint,jj,d,near_suspects,dists_near,
     1      near_nums,k)
        implicit none
        save
        integer *4 k,ifheap
        integer *8 ipoint,jj,near_suspects(k,1),near_nums(1),np,kl
        real *8 dists_near(k,1),d


c        call prinl("in so, ipoint = *", ipoint, 1)
c        call prinl("in so, jj = *", jj, 1)
c        call prin2("in so, d = *", d, 1)
c        call prinl("in so, ns = *", near_suspects(1,ipoint), k)
c        call prin2("in so, dn = *", dists_near(1,ipoint), k)
c        call prinl("in so, near_nums = *", near_nums, 24)
c        call prinf("in so, k = *", k, 1)        

c       
c       if this is the first neighbor for the point
c       ipoint - store it and exit
c
        np=near_nums(ipoint)
c        call prinl("in so, np = *", np, 1)
c
c       if the number of neighbors stored for the point
c       ipoint at this time is less than k, store it and 
c       exit
c
        ifheap=0
        if(np .eq. k-1) ifheap=1
c
        if(np .ge. k) goto 2000
c
        np=np+1
        near_suspects(np,ipoint)=jj
        dists_near(np,ipoint)=d
        near_nums(ipoint)=np
c
c       if this is the element number k, convert the array 
c       dists_near(1,ipoint) into a heap
c
        if(ifheap .eq. 0) return
c
        call peter_heap_build_down(dists_near(1,ipoint),
     1      near_suspects(1,ipoint),np)
c
        return
c
 2000 continue
c
c       there are k neighbors stored; store one more and 
c       reheapify
c
cccc        if(d .ge. dists_near(1,ipoint) ) return
c
c        call prin2("before prod, dn = *", dists_near(1,ipoint),k)
c        call prinl("before prod, ns = *", near_suspects(1,ipoint),k)
c        call prinf("before prod, k = *", k, 1)
c        call prinl("before prod, jj = *", jj, 1)
c        call prin2("before prod, d = *", d, 1)
        kl = k*1
        call peter_replace_one_down(dists_near(1,ipoint),
     1      near_suspects(1,ipoint),kl,d,jj)
c
c        call prin2("after prod, dn = *", dists_near(1,ipoint),k)
c        call prinl("after prod, ns = *", near_suspects(1,ipoint),k)
c        call prinf("after prod, k = *", k, 1)
c        call prinl("after prod, jj = *", jj, 1)
c        call prin2("after prod, d = *", d, 1)

        return
        end






c
c
c       INT8
c
c
        subroutine one_box_lists(ibox,
     1      levels,levmax,iconts)
        implicit none
        save
        integer *8 ia(100),iconts(1),ibox,i,ib,levels(1),levmax
c       
c       find all  boxes on the level levmax that have 
c       contact of type 1 with the box ibox
c
        iconts(1)=ibox
        call to_binary(ibox,levmax,ia)
c
        do 2000 i=1,levmax
c
        ia(i)=-ia(i)
        call find_box(ia,levels,levmax,ib)
c
        iconts(i+1)=ib
        ia(i)=-ia(i)
 2000 continue
c
        return
        end
c
c
c       INT8
c
c
        subroutine to_binary(ii,levmax,iarr)
        implicit none
        save
        integer *8 iarr(1),ii,i,iii,jj,j,levmax
c
        iii=ii
        do 1200 i=1,levmax
c
        j=iii/2
        jj=iii-j*2
        iarr(levmax-i+1)=jj
        iii=j
 1200 continue
c
        do 1400 i=1,levmax
c
        if(iarr(i) .eq. 0) iarr(i)=-1
 1400 continue
c
        return
        end
c
c
c       INT8
c
c
        subroutine find_box(indbox,levels,levmax,ib)
        implicit none
        save
        integer *8 indbox(1),ib,nlast,i,nn,levels(1),levmax,lastlev
c
        lastlev=levels(levmax+1)
        nlast=2**levmax
c
c       find the location in array boxes of the box
c       specified by its array indbox
c
        ib=0
        nn=1
        do 1200 i=levmax,1,-1
c
        ib=ib+(indbox(i)+1)*nn
        nn=nn*2
 1200 continue
c
        ib=ib/2+lastlev       
        return
        end
c
c
c       INT8
c
c
        subroutine struct_bld(n,m,a,
     1      boxes,levels,levmax,nbout,ipts,work,centers)
        implicit none
        save
        integer *4 m
        integer *8 boxes(10,1),ipts(1),n,nbout,i,levmax,lev,
     1             nb,ibox,ijk,nb2,iib1,iib2,ixs,lxs,iys,lys,
     2             iww,idad,icoord,lib1,lib2,lww,levels(1)
        real *8 a(m,n),centers(1),work(1),thresh
c
c       initialize the array ipts
c
        do 1400 i=1,n
c
        ipts(i)=i
 1400 continue
c
c        do 1600 i=1,ndmax
c ANDREI: BUG??? SHOULD BE levmax
c        do 1600 i=1,levmax
c
       do 1600 i = 1, 100
        levels(i)=0
 1600 continue
c
c       initialize the subdivision process
c
        boxes(1,1)=1
        boxes(2,1)=0
        boxes(3,1)=-7
        boxes(4,1)=-1
        boxes(5,1)=-1
        boxes(6,1)=1
        boxes(7,1)=n
        boxes(8,1)=-7
c
        boxes(9,1)=0
        boxes(10,1)=0
c
        nb=1
c
c       separate the array a into two groups with respect to 
c       the coordinate icoord
c
        ibox=1
        nbout=1
        do 2000 ijk=1,levmax

cccc        call prinf('ijk=*',ijk,1)
c
        nb2=nbout
        do 1800 i=1,nb2
c
        if(boxes(10,i) .eq. 7) goto 1800
c
        ibox=i
        iib1=1
        lib1=n+10
c
        iib2=iib1+lib1
        lib2=n+10
c
        ixs=iib2+lib2
        lxs=n+10
c
        iys=ixs+lxs
        lys=n+10
c
        iww=iys+lys
        lww=2*n+20
c
c        if (ibox .EQ. 249) then
c           call prinl("before os, ibox = *", ibox, 1)
c           call prinl("before os, boxes = *", boxes(1,495), 30)
c           stop
c        end if
        call one_subdiv(a,ipts,m,n,ibox,boxes,nb,nbout,
     1      work(iib1),work(iib2),
     2      icoord,thresh,work(ixs),work(iys),work(iww) )
c        if (ibox .EQ. 249) then
c           call prinl("after os, ibox = *", ibox, 1)
c           call prinl("after os, boxes = *", boxes(1,495), 30)
c           stop
c        end if
c

        centers(ibox)=thresh
        idad=boxes(3,ibox)
c
        nb=nbout
 1800 continue
c
 2000 continue
c
c       construct the array lev
c
        lev=-1
        do 2200 i=1,nbout
c
        if(boxes(2,i) .le. lev) goto 2200
c
        lev=boxes(2,i)
        levels(lev+1)=i
 2200 continue
c
        levmax=lev
c
        return
        end
c
c
c       INT8
c
c
        subroutine one_subdiv(a,ipts,m,n,ibox,boxes,nb,nbout,
     1      ib1,ib2,icoord,thresh,xs,ys,ww)
        implicit none
        save
        integer *4 m
        integer *8 n,boxes(10,1),ipts(1),ib1(1),ib2(1),
     1             ibox,nb,nbout,icoord,i1,n_ibox,
     2             ison1,ison2,n2,n1
        real *8 a(m,n),xs(1),ys(1),thresh,ww
c
c       This subroutine subdivides one box in the structure, 
c       obtaining the two sonnie boxes. The information about
c       this felicitous event is stored in the arrays boxes, 
c       ipts
c


c
c       . . . extract from array boxes the information about 
c             the box number ibox
c       
        i1=boxes(6,ibox)
        n_ibox=boxes(7,ibox)
c
c       determine the coordinate to be subdivided
c
        icoord=boxes(2,ibox)+1
c
c       perform the subdivision
c
c        call prinl("before one_sepa, ib1 = *", ib1, n)
c        call prinl("before one_sepa, ib2 = *", ib2, n)
        call one_sepa(a,ipts(i1),m,n_ibox,icoord,n1,n2,ib1,ib2,xs,
     1      thresh,ys,ww)
c
c       enter in array boxes, the information 
c       about the two sonnies                       
c
        nbout=nb+1
        ison1=nbout
        boxes(1,nbout)=nbout
        boxes(2,nbout)=boxes(2,ibox)+1
        boxes(3,nbout)=ibox
        boxes(4,nbout)=-1
        boxes(5,nbout)=-1
        boxes(6,nbout)=boxes(6,ibox)
        boxes(7,nbout)=n1
        boxes(8,nbout)=icoord
        boxes(9,nbout)=0
        boxes(10,nbout)=0
c
        boxes(4,ibox)=nbout
c
        nbout=nbout+1
        ison2=nbout
        boxes(1,nbout)=nbout
        boxes(2,nbout)=boxes(2,ibox)+1
        boxes(3,nbout)=ibox
        boxes(4,nbout)=-1
        boxes(5,nbout)=-1
        boxes(6,nbout)=boxes(6,ibox)+n1
        boxes(7,nbout)=n2
        boxes(8,nbout)=icoord
        boxes(9,nbout)=0
        boxes(10,nbout)=0
c
        boxes(5,ibox)=nbout
        boxes(10,ibox)=7
c
        return
        end
c
c
c       INT8
c
c
        subroutine one_sepa(a,ipts,m,n,icoord,n1,n2,ia,jpts,xs,
     1      thresh,ys,ww)
        implicit none
        save
        integer *4 m
        integer *8 n,ipts(1),jpts(1),ia(1),i,icoord,n1,n2,j,i1
        real *8 a(m,n),xs(1),ys(1),ww(1),thresh,thresh_new,
     1          xmax_lft,xmax_rht,x3,xmin_rht
c
c       sort the coordinates icoord of all points in this box
c
        do 1400 i=1,n
c
        j=ipts(i)
        xs(i)=a(icoord,j)
 1400 continue
c
cccc        call peter_heapsort_rea_part(xs,ia,n,n/2+1)
c
        call array_subdiv(xs,n,n/2,x3,i1,ia,ww)
c
        do 1600 i=1,n
c
        jpts(i)=ipts(ia(i))
        ys(i)=xs(ia(i))
 1600 continue
c
        n1=n/2
        n2=n-n1
c
        thresh=(xs(ia(n/2))+xs(ia(n/2+1)) )/2
c        call prin2("THRESHOLD = *", thresh, 1)
c        call prin2("COORDINATES = *", ys, n)
     
c ANDREI: MAKE THRESHOLD BE EQUAL TO THE MEDIAN
        xmax_lft = -1.0d50
        xmin_rht =  1.0d50
        do 1602 i = 1, n/2
           if (xs(ia(i)) .GE. xmax_lft) xmax_lft = xs(ia(i))
 1602   continue
        do 1603 i = n/2+1, n
           if (xs(ia(i)) .LE. xmin_rht) xmin_rht = xs(ia(i))
 1603   continue 
        thresh_new = (xmax_lft + xmin_rht) / 2
c        call prin2("THRESHOLD_NEW = *", thresh_new, 1)
        if (abs(thresh_new - thresh) .GT. 1D-10) then
c           call prinf("They were different*", i, 0)
        end if
        thresh = thresh_new
c ANDREI: FINISH THRESHOLD
c
        call peter_intcopy(jpts,ipts,n)
c
        return
        end
c
c
c       INT8
c
c

        subroutine array_subdiv(xs,n,nthresh,x3,i1,ia,w)
        implicit none
        save
        integer *8 n,w(1),ia(1),nthresh,i1,i,iiaa,liaa,iibb,libb,
     1             i2,ifout,nn1,nn,nthr
        real *8 xs(1),x1,x2,x3,xmin,xmax
c
c       initialize the array ia
c
        do 1100 i=1,n
c
        ia(i)=i
 1100 continue
c
c       find the minimum and maximum elements in array xs
c
        nthr=nthresh+1
        xmin=1.0d50
        xmax=-1.0d50
c
        iiaa=1
        liaa=n+2
c
        iibb=iiaa+liaa
        libb=n+2
c
        do 1200 i=1,n
c
        if(xs(i) .gt. xmax) xmax=xs(i)
        if(xs(i) .lt. xmin) xmin=xs(i)
 1200 continue
c
        x1=xmin
        x2=xmax
        i1=1
        i2=n
c
        ifout=0
c
        do 2000 i=1,100
c
        if( (i2-i1 .eq. 1) .and. (ifout .eq. 1) ) return
        if(i2-i1 .eq. 1) ifout=1
c
        x3=(x1+x2)/2
c
c       separate the subarray starting at i1 and ending at i2
c       into elements that are smaller than x3, and those otherwise
c
        nn=i2-i1+1
        call array_subdiv_onestep(ia(i1),xs,x1,x2,x3,nn,nn1,
     1      w(iiaa),w(iibb) )
c
        if(i1+nn1 .lt. nthr) then
c    
            i1=i1+nn1
            x1=x3
            goto 2000
        endif
c
        i2=i1+nn1
        x2=x3
c
 2000 continue
c
        return
        end      
c
c
c       INT8
c
c

        subroutine array_subdiv_onestep(ia,xs,xmin,xmax,x3,n,n1,
     1      iaa,ibb)
        implicit none
        save
        integer *8 iaa(1),ibb(1),ia(1),n,n1,ias,ibs,i
        real *8 xs(1),xmin,xmax,x3,x2,x1
c
c       conduct the bisection
c
        x1=xmin
        x2=xmax
        ias=0
        ibs=0
c
        do 1400 i=1,n
c
        if(xs(ia(i)) .le. x3) then
            ias=ias+1
            iaa(ias)=ia(i)
            goto 1400
        endif
c
        ibs=ibs+1
        ibb(ibs)=ia(i)
c
 1400 continue
c
        n1=ias
        do 1600 i=1,n1
c
        ia(i)=iaa(i)
 1600 continue
c
        do 1800 i=1,n-n1
c
        ia(n1+i)=ibb(i)
 1800 continue
c
        return
        end
c
C 
C       NO NEED INT8
C 
C 
        SUBROUTINE dist_compute(a1,a2,m,d)
        IMPLICIT REAL *8 (A-H,O-Z)
        save 
        dimension a1(1),a2(1)
c
        d=0
        do 1200 i=1,m
c
        d=d+(a1(i)-a2(i))**2
 1200 continue
c
cccc        d=sqrt(d)
        return
        end
c
c 
c       INT8
c 
c 
        subroutine peter_copy(a,b,n)
        implicit none
        save
        integer *8 n,i
        real *8 a(1),b(1)
c 
        do 1200 i=1,n
c 
        b(i)=a(i)
 1200 continue
        return
        end
c 
c 
c       INT8
c 
c 
        subroutine peter_intcopy(ia,ib,n)
        implicit none
        save
        integer *8 ia(1),ib(1),n,i
c 
        do 1200 i=1,n
c 
        ib(i)=ia(i)
 1200 continue
        return
        end

c
c
c
c                   Description of array boxes
c
c       boxes(1,i) =i
c       boxes(2,i) - box's level of subdivision
c       boxes(3,i) - box's daddy
c       boxes(4,i) - box's first sonny
c       boxes(5,i) - box's second sonny
c       boxes(6,i) - location in array ipts of box's first point
c       boxes(7,i) - number of points in the box
c       boxes(8,i) - sequence number of coordinate subdivided
c                    to obtain this box from his daddy
c
c       boxes(10,i) - indicator of the state of completion;
c                     boxes(10,i)=0 means that the i-th box
c                     has not been completely processed;
c                     boxes(10,i)=7 means that the i-th box
c                     has been completely processed

c
c
c       INT8
c
c
        subroutine peter_heapsort_rea_down(ra,ia,n)
        implicit none
        save
        integer *8 n,ia(n),i,ii,ii0,m,jj,lone
        real *8 ra(1),dd
c
        do 1600 ii=n,2,-2
c
        ii0=ii/2
        call peter_heapit_rea_down(ra,ia,n,ii0)
 1600 continue
c
c       . . . sort
c
        do 2000 i=1,n-1
c
        m=n-i+1
c
        jj=ia(m)
        ia(m)=ia(1)
        ia(1)=jj
c
        dd=ra(m)
        ra(m)=ra(1)
         ra(1)=dd 
c
        m=m-1
        lone = 1
        call peter_heapit_rea_down(ra,ia,m,lone)
 2000 continue
c
        return
        end
c
c
c       INT8
c
c
        subroutine peter_replace_one_down(ra,ia,n,x,ix)
        implicit none
        save
        integer *8 n,ia(n),ix,lone
        real *8 ra(1),x
c
c       if the new element is greater than ra(ia(1)), exit
c
c        call prin2("in prod, x = *", x, 1)
c        call prin2("in prod, ra(1) = *", ra(1), 1)
        if(x .ge. ra(1) ) return
c
c       replace and re-heapify
c
        ra(1)=x
        ia(1)=ix
        lone=1
        call peter_heapit_rea_down(ra,ia,n,lone)
c
        return
        end        
c
c
c       INT8
c
c
        subroutine peter_heap_build_down(ra,ia,n)
        implicit none
        save
        integer *8 n,ia(n),ii,ii0
        real *8 ra(1)
c
        do 1600 ii=n,2,-2
c
        ii0=ii/2
        call peter_heapit_rea_down(ra,ia,n,ii0)
 1600 continue
c
        return
        end
c
c
c       INT8
c
c
        subroutine peter_heapit_rea_down(ra,ia,n,ii)
        implicit none
        save
        integer *8 n,ia(n),ii,i,ijk,ison1,ison2,id0,id1,id2
        real *8 ra(1),d0,d1,d2

c        call prin2("in phrd, ra = *", ra, n)
c        call prinl("in phrd, ia = *", ia, n)
c        call prinl("in prhd, n = *", n, 1)
c        call prinl("in prhd, ii = *", ii, 1)

c
c       heapify
c
        i=ii
c        call prinl("in phrd, i = *", i, 1)
c        
        do 1400 ijk=1,100
c
        ison1=i*2
c        call prinl("in phrd, ison1 = *", ison1, 1)
        if(ison1 .gt. n) return
c
        ison2=i*2+1
        if(ison2 .gt. n) then
c
            d0=ra(i)
            d1=ra(ison1)
c
            id0=ia(i)
            id1=ia(ison1)
c
            if(d0 .ge. d1) return
c
            ra(i)=d1
            ra(ison1)=d0
            ia(i)=id1
            ia(ison1)=id0
c
            return
         endif
c
c
        d0=ra(i)
        d1=ra(ison1)
        d2=ra(ison2)
c
        id0=ia(i)
        id1=ia(ison1)
        id2=ia(ison2)

c        call prinl("in phrd, i = *", i, 1)
c        call prinl("in phrd, ison1 = *", ison1, 1)
c        call prinl("in phrd, ison2 = *", ison2, 1)
c        call prinl("in phrd, id0 = *", id0, 1)
c        call prinl("in phrd, id1 = *", id1, 1)
c        call prinl("in phrd, id2 = *", id2, 1)
c        call prin2("in phrd, d0 = *", d0, 1)
c        call prin2("in phrd, d1 = *", d1, 1)
c        call prin2("in phrd, d2 = *", d2, 1)

c
        if( (d0 .ge. d1) .and. (d0 .ge. d2) )  return
        if(d1 .ge. d2 ) then
            ra(i)=d1
            ra(ison1)=d0
c
            ia(i)=id1
            ia(ison1)=id0
            i=ison1
c
            goto 1400
        endif
c
            ra(i)=d2
            ra(ison2)=d0
c
            ia(i)=id2
            ia(ison2)=id0
            i=ison2
c
 1400 continue
            return
            end

c
c
c       INT8
c
c
        subroutine peter_heapsort_rea_up(ra,ia,n,k,kk)
        implicit none
        save
        integer *8 n,ia(n),k,kk,ii,ii0,jj,i,m,lone
        real *8 ra(1),dd
c
c        This subroutine sorts the array ra of n real numbers
c        in decreasing order starting from the bottom, and 
c        terminates the process when k numbers with different 
c        keys ia have been encountered. Thus, on exit the last
c        kk elements of ra are the smallest kk elements, with
c        the smallest one occupying the location ra(n). The 
c        keys ia are subjected to the same permutation as the
c        elements of ra.
c
c        . . . heapify
c
        do 1600 ii=n,2,-2
c
        ii0=ii/2
        call peter_heapit_rea_up(ra,ia,n,ii0)
 1600 continue
c
c       . . . sort
c
        kk=1
        do 2000 i=1,n-1
c
        m=n-i+1
c
        jj=ia(m)
        ia(m)=ia(1)
        ia(1)=jj
c
        dd=ra(m)
        ra(m)=ra(1)
        ra(1)=dd 
c
        if( (i .ne. 1) .and. (ia(m) .ne. ia(m+1)) ) kk=kk+1
        if(kk .eq. k) return
c
        m=m-1
        lone = 1
        call peter_heapit_rea_up(ra,ia,m,lone)
c
 2000 continue
c
        return
        end
c
c
c       INT8
c
c
        subroutine peter_heapit_rea_up(ra,ia,n,ii)
        implicit none
        save
        integer *8 n,ia(n),ii,i,ijk,ison1,ison2,id0,id1,id2
        real *8 ra(1),d0,d1,d2
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
            d0=ra(i)
            d1=ra(ison1)
c
            id0=ia(i)
            id1=ia(ison1)
c
            if(d0 .le. d1) return
c
            ra(i)=d1
            ra(ison1)=d0
            ia(i)=id1
            ia(ison1)=id0
c
            return
         endif
c
c
        d0=ra(i)
        d1=ra(ison1)
        d2=ra(ison2)
c
        id0=ia(i)
        id1=ia(ison1)
        id2=ia(ison2)

c
        if( (d0 .le. d1) .and. (d0 .le. d2) )  return
        if(d1 .le. d2 ) then
            ra(i)=d1
            ra(ison1)=d0
c
            ia(i)=id1
            ia(ison1)=id0
            i=ison1
c
            goto 1400
        endif
c
            ra(i)=d2
            ra(ison2)=d0
c
            ia(i)=id2
            ia(ison2)=id0
            i=ison2
c
 1400 continue
            return
            end

c
c
c       INT8
c       OLD (liberal) VERSION
c
c
        subroutine second_search_old(a,m,n,k,near_convicts,
     1      dists_convicts,
     2      iarr,rarr,iarr2,ia,rarr2)
        implicit none
        save
        integer *4 m,k
        integer *8 n,near_convicts(k,1),iarr(1),iarr2(1),ia(1),
     1             icc,i,j,ii,iii,jjj,jj,j1,kk,lk
        real *8 a(m,n),dists_convicts(k,1),
     1      rarr(1),rarr2(1),d
c
c       for each of the n points, collect together its 
c       putative neighbors and the putative neighbors of its
c       neighbors, find the nearest k neighbors among this
c       elite group of citizens, and put the results back 
c       in the arrays near_convicts, dists_convicts
c
        icc=0
        do 3000 i=1,n
c
        icc=icc+1
        if(icc .eq. 1000) then
            call prinl('in second_search2, i=*',i,1)
            icc=0
        endif
c
c        construct the derivative group for the i-th point
c
        do 1200 j=1,k
c
        iarr(j)=near_convicts(j,i)
        rarr(j)=dists_convicts(j,i)
 1200 continue
c  
        ii=k      
        do 1600 j=1,k
c
        j1=iarr(j)
        do 1400 jj=1,k
c
        if(near_convicts(jj,j1) .eq. i) goto 1400
        ii=ii+1
        iarr(ii)=near_convicts(jj,j1)
 1400 continue
 1600 continue
c
c       calculate the distances between the i-th point
c       and the points in the array iarr

        iii=0
        do 2200 j=k+1,ii
c
        
        jjj=iarr(j)

        call dist_compute(a(1,i),a(1,jjj),m,d)

        if(d .gt. rarr(k)) goto 2200
c
        iii=iii+1
        rarr(iii+k)=d
        iarr(iii+k)=iarr(j)



 2200 continue

        ii=iii+k


c
c       . . . sort
c

        lk = k
        call peter_heapsort_rea_up(rarr,iarr,ii,lk,kk)
c
c       extract from arrays iarr, rarr the non-repeating 
c       elements 
c
        jjj=1
        iarr2(1)=iarr(ii)
        rarr2(1)=rarr(ii)
c
        do 2300 j=ii-1,1,-1
c
        if(iarr(j) .eq. iarr(j+1)) goto 2300
c
        jjj=jjj+1
        iarr2(jjj)=iarr(j)
        rarr2(jjj)=rarr(j)
c
        if(jjj .eq. k) goto 2350
 2300 continue
 2350 continue
c
c
c
c       copy the first k elements of the sorted array into 
c       the array near_convicts
c
        do 2400 j=1,k
c
        jjj=ia(j)
        jjj=j
        near_convicts(j,i)=iarr2(jjj)
        dists_convicts(j,i)=rarr2(jjj)
 2400 continue
c
 3000 continue
c
c        stop

        return
        end


c
c
c       INT8
c
c
        subroutine second_search(a,m,n,k,near_convicts,
     1      dists_convicts,
     2      iarr,rarr,iarr2,ia,rarr2,near2,dists2,itimes)
        implicit none
        save
        integer *4 m,k,itimes(15),it1,it2
        integer *8 n,near_convicts(k,1),iarr(1),iarr2(1),ia(1),
     1             icc,i,j,ii,iii,jjj,jj,j1,kk,lk,near2(k,1)
        real *8 a(m,n),dists_convicts(k,1),
     1      rarr(1),rarr2(1),d,dists2(k,1)
c
c       for each of the n points, collect together its 
c       putative neighbors and the putative neighbors of its
c       neighbors, find the nearest k neighbors among this
c       elite group of citizens, and put the results back 
c       in the arrays near_convicts, dists_convicts
c
        icc=0
        do 3000 i=1,n
c
        icc=icc+1
        if(icc .eq. 1000) then
            call prinl('in second_search2, i=*',i,1)
            icc=0
        endif
c
c        construct the derivative group for the i-th point
c
        do 1200 j=1,k
c
        iarr(j)=near_convicts(j,i)
        rarr(j)=dists_convicts(j,i)
 1200 continue
c  
        ii=k      
        do 1600 j=1,k
c
        j1=iarr(j)
        do 1400 jj=1,k
c
        if(near_convicts(jj,j1) .eq. i) goto 1400
        ii=ii+1
        iarr(ii)=near_convicts(jj,j1)
 1400 continue
 1600 continue
c
c       calculate the distances between the i-th point
c       and the points in the array iarr

c APAR
c         if (i .LE. 320) then
c            call prinl("i = *", i, 1)
c            call prinl("iarr = *", iarr, ii)
c            if (i .GT. 10) stop
c         endif

c
cccc        call prin2('rarr=*',rarr,k)

        iii=0
        do 2200 j=k+1,ii
c
        
        jjj=iarr(j)
cccc        call dist_compute(a(1,i),a(1,jjj),m,rarr(j))

cccc        goto 2200

        call dist_compute(a(1,i),a(1,jjj),m,d)

        if(d .gt. rarr(k)) goto 2200
c
        iii=iii+1
        rarr(iii+k)=d
        iarr(iii+k)=iarr(j)



 2200 continue

cccc        call prinf('ii=*',ii,1)
cccc        call prinf('and iii+k=*',iii+k,1)


        ii=iii+k


c
c       . . . sort
c
cccc        call peter_copy(iarr,iarr7,ii)
cccc        call peter_copy(rarr,rarr7,ii)

        lk = k
        call peter_heapsort_rea_up(rarr,iarr,ii,lk,kk)
cccc        call peter_heapsort_rea_up(rarr7,iarr7,ii,k,kk)
c
c       extract from arrays iarr, rarr the non-repeating 
c       elements 
c
        jjj=1
        iarr2(1)=iarr(ii)
        rarr2(1)=rarr(ii)
c
        do 2300 j=ii-1,1,-1
c
        if(iarr(j) .eq. iarr(j+1)) goto 2300
c
        jjj=jjj+1
        iarr2(jjj)=iarr(j)
        rarr2(jjj)=rarr(j)
c
        if(jjj .eq. k) goto 2350
 2300 continue
 2350 continue
c
cccc        call prinf('and iarr2 after sorting*',iarr2,k)
cccc        call prin2('and rarr2 after sorting*',rarr2,k)
c
c
c       copy the first k elements of the sorted array into 
c       the array near_convicts
c
        do 2400 j=1,k
c
        jjj=ia(j)
        jjj=j
c        near_convicts(j,i)=iarr2(jjj)
c        dists_convicts(j,i)=rarr2(jjj)
        near2(j,i)=iarr2(jjj)
        dists2(j,i)=rarr2(jjj)
 2400 continue
c
 3000 continue
c
c        stop

        it1 = time()
        do 3100 i = 1, n
           do 3110 j = 1, k
              near_convicts(j,i) = near2(j,i)
              dists_convicts(j,i) = dists2(j,i)
 3110      continue
 3100   continue
        it2 = time()
        itimes(2) = it2-it1

        return
        end

c
c
c       INT8
c
c
        subroutine get_stats(n,m,a,numit,
     1      near_convicts,dists_convicts,k,w,lused,ntest,
     2      iw,istat)
        implicit none
        save
        integer *4 k,m,numit,istat
        integer *8 n,near_convicts(k,1),lused,ntest,iw,
     1             iia,lia,idists,ldists,idists2,ldists2,
     2             iiaa,liaa,iib,lib,iiperm,liperm
        real *8 a(m,n),w(1),
     1      dists_convicts(k,1)
c
c       . . . allocate memory
c
        iia=1
        lia=n+2
c
        idists=iia+lia
        ldists=n+2
c
        idists2=idists+ldists
        ldists2=n+2
c
        iiaa=idists2+ldists2
        liaa=n+2
c
        iib=iiaa+liaa
        lib=2*k+10

        iiperm=iib+lib

        call prinl("iia = *", iia, 1)
        call prinl("lia = *", lia, 1)
        call prinl("idists = *", idists, 1)
        call prinl("ldists = *", ldists, 1)
        call prinl("idists2 = *", idists2, 1)
        call prinl("ldists2 = *", ldists2, 1)
        call prinl("iiaa = *", iiaa, 1)
        call prinl("liaa = *", liaa, 1)
        call prinl("iib = *", iib, 1)
        call prinl("lib = *", lib, 1)
        call prinl("iiperm = *", iiperm, 1)
        call prinf("istat = *", istat, 1)

c
        call get_stats0(n,m,a,numit,
     1      near_convicts,dists_convicts,k,w,lused,ntest,iw,
     2      w(iia),w(idists),w(idists2),w(iiaa),w(iib),
     3      w(iiperm),istat)


        return
        end         
c
c
c       INT8
c
c
        subroutine get_stats0(n,m,a,numit,
     1      near_convicts,dists_convicts,k,w,lused,ntest,iw,
     2      ia,dists,dists2,iaa,ib,iperm,istat)
        implicit none
        save
        integer *4 k,m,numit,istat,ismall
        integer *8 n,near_convicts(k,1),ia(1),iaa(1),ib(1),iw,
     1             lused,ntest,nsmall,num_belong_tot,nwrong,ic,i,
     2             j,num_belong,nleft,kl,k2l,iperm(1),ii
        real *8 a(m,n),w(1),
     1      dists_convicts(k,1),
     2      dists(1),dists2(1),rat_nearest_aver,sum_true,sum_appr,
     3      sum_true2,sum_appr2,sum_true_sq,sum_appr_sq,
     4      d,dist_true_kth,dist_appr_kth,done,dnum_belong_average,
     5      prop_belong,var_true_k,rat_nearest,small_prop,var_appr_k
c
c        analyze the results and print a summary
c
c
c        . . . for the first ntest points, construct the
c              nearest k neighbors via brute force
c
        call prinf("istat = *", istat, 1)

        ntest=min(1000,n)
        call prinl("ntest = *", ntest, 1)

        if (istat .EQ. 1) then
           do 1011 i = 1, ntest
              iperm(i) = i
 1011      continue
           ismall = ntest
c           call prinl("iperm = *", iperm, ismall)
        end if

        if (istat .EQ. 2) then
           call corrand4_longs_no_reps(iperm, ntest, n)
           do 1012 i = 1, ntest
              iperm(i) = iperm(i) + 1
 1012      continue
           ismall = ntest
c           call prinl("iperm = *", iperm, ismall)
c           stop
        end if

c        stop

        nsmall=0
        num_belong_tot=0
        rat_nearest_aver=0
        nwrong=0
c
        sum_true=0
        sum_appr=0
c ANDREI
        sum_true2 = 0
        sum_appr2 = 0

        sum_true_sq = 0
        sum_appr_sq = 0
        dnum_belong_average = 0

        ic=0
        do 3000 ii=1,ntest
        
        i = iperm(ii)
c
        ic=ic+1
            if(ic .eq. 100) then
            call prinf("in get_stats, ii=*",ii,1)
            ic=0
        endif    
c
c       . . . construct the distances for the i-th point
c
        do 1600 j=1,n
c
        call dist_compute(a(1,i),a(1,j),m,d)
c
        dists(j)=d
        if(i .eq. j) dists(j)=1.0d30
        if(dists(j) .lt. dists_convicts(1,i)) nsmall=nsmall+1
 1600 continue
c        call prin2("dists = *", dists, n)
c
c       sort the distances from the i-th point to all other
c       points, and take the first k elements
c
        kl = k
        call peter_heapsort_rea_part(dists,ia,n,kl)
c        call prin2("dists = *", dists, n)
c        call prinl("ia = *", ia, k)
c        call prinl("i = *", i, 1)

        dist_true_kth = 0
        dist_appr_kth = 0
c
        do 1800 j=1,k
c
        dists2(j)=dists(ia(j))

c  ANDREI
        sum_true=sum_true+sqrt(dists2(j))
        sum_true2 = sum_true2 + dists2(j)
        dist_true_kth = dist_true_kth + dists2(j)

        sum_appr=sum_appr+sqrt(dists_convicts(j,i))
        sum_appr2 = sum_appr2 + dists_convicts(j,i)
        dist_appr_kth = dist_appr_kth + dists_convicts(j,i)
 1800 continue

        dist_true_kth = dist_true_kth / k
        dist_appr_kth = dist_appr_kth / k
        sum_true_sq = sum_true_sq + dist_true_kth**2
        sum_appr_sq = sum_appr_sq + dist_appr_kth**2

c
c ANDREI: PRINT THE DISTANCES
c        call prinf("i = *", i, 1)
c        call prin2("TRUE k distances = *", dists2, k)
c        call prin2("APPR k distances = *", dists_convicts(1, i), k)
c        call prinf("APPR k indices = *", near_convicts(1, i), k)
        

        rat_nearest=dists_convicts(1,i)/dists2(1)
        rat_nearest_aver=rat_nearest_aver+rat_nearest
c
        if(near_convicts(1,i) .ne. ia(1)) nwrong=nwrong+1
c
c       count the number of discrepancies between the "brute force"
c       and approximate nearest k neighbors for this point
c
        do 2000 j=1,k
c
        dists(j)=dists2(j)
        dists(j+k)=dists_convicts(j,i)
        ia(k+j)=near_convicts(j,i)
 2000 continue
c
c       . . . sort them things
c
        k2l = 2*k
        call peter_heapsort_inte(ia,iaa,k2l)
c
        do 2200 j=1,2*k
c
        ib(j)=ia(iaa(j))
 2200 continue
c
c       squeeze out the doubles
c
        nleft=1
        k2l = 2*k
        do 2400 j=2,k2l
        if(ib(j-1) .ne. ib(j)) nleft=nleft+1
 2400 continue
c
        done=1
        prop_belong=(nleft-k)*done/k        
        k2l = 2*k
        num_belong=k2l-nleft
        num_belong_tot=num_belong_tot+num_belong

c ANDREI: PRINT NUM_BELONG
c        call prin2("point coordinates = *", a(1, i), m)
c        call prinf("num_belong = *", num_belong, 1)
c
 3000 continue
c
        done=1
        small_prop=nsmall*done/n/ntest
c
        iw=13
c

c       ANDREI: COMMENT WRITE
c        write(6,3400)
c        write(iw,3400)
c        write(6,3400)
c        write(iw,3400)

 
 3050 format(15x,'Report on The Original Nearest Neighbor Search')
c       ANDREI: COMMENT WRITE
c        write(6,3050)       
c        write(iw,3050)       
c
c       ANDREI: COMMENT WRITE
c        write(6,3400)
c        write(iw,3400)
c        write(6,3400)
c        write(iw,3400)
c
 3100 format(6x,'In this run, the number of points n is',18x,i7,/,
     1       6x,'The dimension m is',36x,i7,/,
     1       6x,'The number of neighbors for each point k is',11x,i7,/,
     2       6x,'The number of iterations numit is',25x,i3,/,
     2       6x,'The number of test points ntest is',24x,i5)

c       ANDREI: COMMENT WRITE
c        write(6,3100) n,m,k,numit,ntest
c        write(iw,3100) n,m,k,numit,ntest

c       ANDREI: COMMENT WRITE
c        write(6,3400)
c        write(iw,3400)

 3200 format(6x,'Average proportion of points that are closer'
     1      ' to one',/,6x,'of the first NTEST points than the',
     2      ' nearest neighbor',/,6x,
     3      'as determined by the subroutine: ', 25x,E11.5)
c
c       ANDREI: COMMENT WRITE
c        write(6,3200) small_prop
c        write(iw,3200) small_prop

 3400 format(10x)

c       ANDREI: COMMENT WRITE
c        write(6,3400)
c        write(iw,3400)
c
 3600 format(6x,'The average ratio RAT between the TOTAL',
     1      ' sum of',/,6x, 'distances',
     2      ' of allegedly nearest k points to the same',/,6x,
     3      'sum of distances ','to true nearest k points,',
     4      ' taken',/,6x,'over the first NTEST',
     5      ' points:',30x,E11.5) 
c
c       ANDREI: COMMENT WRITE
c         write(6,3600) sum_appr/sum_true
c         write(iw,3600) sum_appr/sum_true
c
 3800 format(6x,'And RAT-1 is ',45x,E11.5)
c
c       ANDREI: COMMENT WRITE
c         write(6,3800) sum_appr/sum_true-1
c         write(iw,3800) sum_appr/sum_true-1
            
c       ANDREI: COMMENT WRITE
c         write(6,3400)
c         write(iw,3400)
c
        rat_nearest_aver=rat_nearest_aver/ntest
c
 4200 format(6x,'The average RAT_AVER of the individual',
     1      ' ratios between',/,6x,'the distance to the putative ',
     2      'nearest points',
     1      ' to the ',/,6x,
     2      'distances ','to their respective true nearest points,',
     3      /,6x,'taken over',' the first NTEST',
     3      ' points:',24x,E11.5) 
c
c       ANDREI: COMMENT WRITE
c         write(6,4200) rat_nearest_aver
c         write(iw,4200) rat_nearest_aver
c
 4400 format(6x,'And RAT_AVER-1 is ',40x,E11.5)
c
c       ANDREI: COMMENT WRITE
c         write(6,4400) rat_nearest_aver-1
c         write(iw,4400) rat_nearest_aver-1
            
c       ANDREI: COMMENT WRITE
c         write(6,3400)
c         write(iw,3400)

c
 4600 format(6x,'The number of nearest points determined ',
     2      'incorrectly',/,6x,'for the first NTEST points:',36x,i6)
c
c       ANDREI: COMMENT WRITE
c         write(6,4600) nwrong
c         write(iw,4600) nwrong
c
 4800 format(6x,'And the proportion of nearest points determined ',
     2      /,6x,'incorrectly:',46x,e11.5)
c
c       ANDREI: COMMENT WRITE
c         write(6,4800) nwrong*done/ntest
c         write(iw,4800) nwrong*done/ntest
            
c       ANDREI: COMMENT WRITE
c         write(6,3400)
c         write(iw,3400)
c
c
         dnum_belong_average=num_belong_tot*done/ntest
c
 5000 format(6x,'On the average, among the k allegedly nearest points ',
     1     /,6x,'found by the algorithm, the number that belong there:',
     2      5     x,E11.5)
c
c       ANDREI: COMMENT WRITE
c         write(6,5000) dnum_belong_average
c         write(iw,5000) dnum_belong_average
c
c       ANDREI: COMMENT WRITE
c         write(6,3400)
c         write(iw,3400)
c
 5400 format(6x,'The memory used in real *8 words: ',25x,i10.0)
c
c       ANDREI: COMMENT WRITE
c         write(6,5400) lused
c         write(iw,5400) lused
c
c       ANDREI: COMMENT WRITE
c         write(6,3400)
c         write(iw,3400)

         call prin2("SQ ratio(AVG) = *", sum_appr2/sum_true2, 1)

         call prin2("avg(true_k) = *", sum_true2/(ntest*k), 1)
         call prin2("avg(appr_k) = *", sum_appr2/(ntest*k), 1)

         call prinf("ntest = *", ntest, 1)
         call prin2("sum_true_sq = *", sum_true_sq, 1)

         var_true_k = ((sum_true_sq)/ntest - 
     1                 (sum_true2/(ntest*k))**2)*ntest/(ntest-1)

         var_appr_k = ((sum_appr_sq)/ntest - 
     1                 (sum_appr2/(ntest*k))**2)*ntest/(ntest-1)

         call prin2("var(true_k) = *", var_true_k, 1)
         call prin2("var(appr_k) = *", var_appr_k, 1)
         call prin2("err(true) = *", sqrt(var_true_k/ntest)*2, 1)
         call prin2("err(appr) = *", sqrt(var_appr_k/ntest)*2, 1)

        w(1) = sum_true2 / (ntest*k)
        w(2) = sum_appr2 / (ntest*k)
        w(3) = dnum_belong_average / k

        return
        end

 
c
c       INT8
c       all but m,k,numit
c
c

        subroutine get_memory_size(n,m,k,numit,isuper,nwords)
        implicit none

        integer *4 m,k,numit,isuper
        integer *8 n,nwords,levmax,nsteps,nboxes_max,keepw_max,
     1             lenstore,iindsv,lindsv,iv,lv,keep,iipts,lipts,
     2             iib,lib,irs,lrs,iib2,lib2,irs2,lrs2,iia,lia,
     3             inear_suspects,lnear_suspects,idists_near,
     4             ldists_near,iw,ibb,lbb,keep_appr,iboxes,lboxes,
     5             iinboxes,linboxes,iiconts,liconts,iiwork,liwork,
     6             isuspects,lsuspects,iinbox,linbox,
     7             inear_nums,lnear_nums,icenters,lcenters,ndmax,
     8             lused,iiarr,liarr,irarr,lrarr,iiarr2,liarr2,
     9             irarr2,lrarr2,nwords1,nwords2
        real *8 done,d

        done=1
        d=log(n*done/k)/log(2*done)
        levmax=d
        call prinl('levmax=*',levmax,1)
        if (levmax .GT. m) levmax = m

        nsteps=6
        nboxes_max=2**levmax*2 
        keepw_max=3*m*nsteps+10*m

        call prinl('keepw_max=*',keepw_max,1)

c
        lenstore=nboxes_max*10+levmax+nboxes_max+keepw_max
     1      +nboxes_max*(levmax+1)+n+m+1000
c
c        lenstore=lenstore*numit
c
        iindsv=100
        lindsv=numit+100
c
        iv=iindsv+lindsv
        lv=lenstore+1
c
        keep=iv+lv 
c
c       ANDREI: MEMORY
c
        call prinl("iindsv = *", iindsv, 1)
        call prinl("iv = *", iv, 1)
        call prinl("iv+lv = keep = *", keep, 1)

c
c        . . . allocate memory
c
        iipts=1 + lenstore+1100
c
c       ANDREI: MEMORY FIX
c       iipts might get less than iv+lv
c       (if numit > 900)
c       in this case w(iipts) overrides w(iv) !!!!
c       solution: make it later
c
c
        iipts = keep + 1000
        lipts=n+2
c
c
c
        call prinl("iipts = *", iipts, 1)

c
        iib=iipts+lipts
        lib=2*k+2
c
        irs=iib+lib
        lrs=2*k+2
c
        iib2=irs+lrs
        lib2=2*k+2
c
        irs2=iib2+lib2
        lrs2=2*k+2
c
        iia=irs2+lrs2
        lia=2*k+2
c
        inear_suspects=iia+lia  
        lnear_suspects=n*k+100
c
        idists_near=inear_suspects+lnear_suspects
        ldists_near=n*k+100
c
        iw=idists_near+ldists_near
c        lenw2=lenw-iw-1

c
c       ANDREI: SIZE
c       TO PREVENT FALLING FROM INSIDE ALL_ITER
c       EVEN IF REDUCE IT BY ONE: STOPPED
c
        nwords = iw+1+2000

c
c       NOW WITH ALL_ITER0
c
        ibb=1
        lbb=m*n+100
c
c       EXACTLY TO PREVENT FALLING FROM INSIDE ALL_ITER0
c
        nwords = nwords + ibb+lbb

c
c       CONSTRAINT 1
c
        keep_appr=3*m*nsteps+10*m
        nwords1 = nwords + keep_appr+100+m

c
c       CONSTRAINT 2
c

        iboxes=keep_appr+100
c       ANDREI: MEMORY MULTIPLY BY 2, NOT ENOUGH FOR BOXES
        lboxes=2*(2**levmax)*10+100
c
        iinboxes=iboxes+lboxes
        linboxes=n+10
c
        iiconts=iinboxes+linboxes
        liconts=1000
c
        iiwork=iiconts+liconts
        liwork=n*6+200
c
        isuspects=iiwork+liwork
        lsuspects=n+10
c
        irs=isuspects+lsuspects
        lrs=n+10
c
        iinbox=irs+lrs
        linbox=k*10+10
c
        iia=iinbox+linbox
        lia=n+10
c
        inear_nums=iia+lia
        lnear_nums=n+10
c
        icenters=inear_nums+lnear_nums
        lcenters=4*2**levmax+100
c
        ndmax=levmax+1
        lused=2**levmax*(10+ndmax)+7*n+10*k+2000+4*2**levmax+100

        nwords = nwords+max(lused,icenters+lcenters)
        nwords = max(nwords, nwords1)
        
        nwords = nwords + 5*k*(k+1)

c       
c       NEW SUPERCHARGING, REQUIRES MORE MEMORY
c
        if (isuper .LT. 2) goto 999

        iiarr = keep + 1000
        liarr=n+2
        if ( (k*(k+1)) .GT. n ) liarr = k*(k+1)+2
c
        irarr=iiarr+liarr
        lrarr=n+2             
        if ( (k*(k+1)) .GT. n ) lrarr = k*(k+1)+2
c
        iiarr2=irarr+lrarr
        liarr2=n+2
        if ( (k*(k+1)) .GT. n ) liarr2 = k*(k+1)+2
c
        iia=iiarr2+liarr2
        lia=n+2
        if ( (k*(k+1)) .GT. n ) lia = k*(k+1)+2
c
        irarr2=iia+lia
        lrarr2=n+2
        if ( (k*(k+1)) .GT. n ) lrarr2 = k*(k+1)+2
c
        inear_suspects=irarr2+lrarr2
        lnear_suspects=n*k+2
c
        idists_near=inear_suspects+lnear_suspects
        ldists_near=n*k+2

        nwords2 = idists_near + ldists_near
        call prinl("nwords = *", nwords, 1)
        call prinl("nwords2 = *", nwords2, 1)

        if (nwords2 .GT. nwords) nwords = nwords2
        call prinl("nwords = *", nwords, 1)

 999    continue

        return
        end

c
c       RANDOMIZATION ROUTINES
c       copied from test20.f (in /test_par/)
c

c
c       ixs, lixs, n_max: all integer *8
c       answer: ixs(i) = 0, ..., n_max-1
c       no repetitions
c       a rather slow algorithm
c
        subroutine corrand4_longs_no_reps(ixs, lixs, n_max)
        implicit none

        integer *8 lixs, ixs(1), n_max, i, il, jj
        integer *4 m, n1, n2, j, i1, i2, k

        if (lixs .GT. n_max) then
           call prinl("in no_reps, lixs = *", lixs, 1)
           call prinl("in no_reps, TOO SMALL n_max = *", n_max, 1)
           return
        end if

        m = 2**30
        n1 = n_max / m + 1
c        call prinf("n1 = *", n1, 1)
        n2 = m
        if (n1 .EQ. 1) n2 = n_max
c        call prinf("n2 = *", n2, 1)

        do 200 i = 1, lixs

           do 215 k = 1, 10000
c ONE NUMBER        
           do 210 j = 1, 1000
              call corrand4_integers(i1, 1, n2)
              call corrand4_integers(i2, 1, n1)
c              call prinf("i1 = *", i1, 1)
c              call prinf("i2 = *", i2, 1)
              ixs(i) = 1 * m
              ixs(i) = ixs(i) * i2
              ixs(i) = ixs(i) + i1
c              call prinl("ixs(i) = *", ixs(i), 1)
              if (ixs(i) .LT. n_max) goto 220
 210       continue
 220       continue

           do 230 jj = 1, i-1
              if (ixs(i) .EQ. ixs(jj)) goto 215
 230       continue
c NO REPETITION!
           goto 250

c TRY TO CHOOSE ONCE AGAIN
 215       continue

 250       continue

 200    continue
        
        return
        end



c
c       both n_max and ixs are integers *4
c       so is lixs
c       ixs(i) is a random number between 0 and n_max-1
c
        subroutine corrand4_integers(ixs, lixs, n_max)
        implicit none

        integer *4 lixs, ixs(1), n_max, i
        real *8 r

        do 100 i = 1, lixs
           call corrand4(1,r)
           ixs(i) = n_max*r
 100    continue

        return
        end

