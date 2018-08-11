      program Simulate
      
      implicit none
c 
c This program is to test the mechanics of adding a tubulin to the lattice
c
c variables:
c	state(i,j) is the lattice. 1 = GTP-tub, 0=emtpy, -1=GDP=tub
c	ev(i,j,k) returns the kth event number at site i,j
c 	nev(i,j) keeps the number of possible events at site i,j
c	evtype(i,j,k) returns the type of event at i,j: k=1 -> addition, k=2 -> dissoc, k=3 -> gtp hydrolysis, k=4 -> nucleotide exchange
c   xloc(k), yloc(k) return the lattice positions of a particular event k 
c 	addstate(i,j) = 1 if there is a possible addition event at i,j; it is 0 otherwise
c   dissocstate(i,j): 1-add, 2-1lat, 3-1long (gtp), 4-1long (gdp), 5-1each (GTP), 6-1 each gdp, 7-gtp hydro, 8-slow

      integer i,j,state(13,2000),totev
      integer evtype(13,2000,4),ev(13,2000,4),addstate(13,2000)
      integer hydrostate(13,2000),substate(13,2000)
      integer xloc(20000),yloc(20000),evno,nrecyc,recyc(20000)
      integer ntub,ontub,nev(13,2000),maxlen,k,l,xo,yo,xom,xop,yom,yop
      logical qinit,done,qdebug
   	  integer dissocstate(13,2000),dist(13,2000),dep(13,2000)   
	  integer temp,kk,x,y,jj,ll,n
	  integer test,event(20000),revmap(20000)
	  integer mono,di,tri,tet,pent,hex,hept,oct,ennea,dec
	  integer nhydro,lmax,corner_y,corner_n
	  double precision k_on,conc,GDP,K_lat,K_long,K_both,k_hyd,time
	  double precision K_thr_long,K_thr_lat,K_all,k_ne
	  double precision evtime(20000)
	  double precision get_evtime
c For seam multiplier
      double precision seammult(13,2000)
      
      double precision dummy
      logical depoly,poly		!to see whether depolymerization has started
      
c For depoly slope calculations...
      double precision length_d(300),time_d(300),tmlngth_d(300),tm_dsq(300)
      integer count_d,count_eqd
      double precision math_dt,pretm,tmlngth_dsm,tm_dsm,lngth_dsm,tm_dsqsm,rate_d
      
c For growth slope calculations...
      integer count_p,count_eqp
      double precision time_p(800),length_p(800),tmlngth_p(800),tm_psq(800)
      double precision tmlngth_psm,tm_psqsm,rate_p
      
c Initialization for depolymerization...      
      depoly=.false.
      count_d=0
      tmlngth_dsm=0
      tm_dsm=0
      lngth_dsm=0
      tm_dsqsm=0
      
c Initialization for growth...
      poly=.true.
      count_p=0
      tmlngth_psm=0
      tm_psqsm=0
	  
c initialize: no events to recycle and number of tubulin starts at 1
      time=0.
      maxlen=2000
      nrecyc = 0
      ntub = 1
      ontub = 0   !Old number of tubulin
      qinit=.TRUE.
      xo=7
      yo=100
c initialize: put a MTstub in the lattice, 5 tall on each protofilament
      do i=1,13
         do j=1,2000
            ev(i,j,1) = 0
            ev(i,j,2) = 0
            ev(i,j,3) = 0
            ev(i,j,4) = 0 !new for nuc ex
            nev(i,j)=0
            state(i,j)=0
            addstate(i,j) = 0
            dissocstate(i,j)=0
            dep(i,j)=9
            seammult(i,j)=1.0
         end do
      end do
      
      do i=1,13
         do j=98,102
            dep(i,j)=0
            state(i,j) = 1
            dissocstate(i,j)=8
         end do
      end do

      conc=10.0E-6
c      GDP=8.2
c      k_on=2.E6
c allowing K_lat for non-crashing purposes (hopefully)...
      K_lat=1.0
c      K_long=82.7E-6
c      K_both=3.4E-6
c      K_thr_long=2.7E-10
c      K_thr_lat=1.3E-7
c      K_all=1.1E-18
      GDP=12.2
      k_on=4.0E6
c      K_long=1.5E-7
      K_long=1.11E-3
c      K_both=4.0E-6
      K_both=3.72E-6
c      K_thr_long=1.1E-18	! guess
      K_thr_long=2.7E-11	! guess
c      K_thr_lat=1.1E-18
      K_thr_lat=1.24E-8
      K_all=1.1E-18
c      k_hyd=1.
      k_hyd=1.1E-18
      k_ne=0.2				!new for nuc ex (estimated, but still dummy)
c      k_ne=0.5

      write(6,*) 'rates',k_on*conc,k_on*K_lat,k_on*K_long,k_on*K_both
      call random_seed

c initializing counters for different n-mers      
      mono=0
      di=0
      tri=0
      tet=0
      pent=0
      hex=0
      hept=0
      oct=0
      ennea=0
      dec=0
      
c initializing counter for hydrolysis events
      nhydro=0
C
C create new events, adding on top
      do i=1,13
         if (i.lt.13) then
            addstate(i,103) = 1
            dissocstate(i,103) = 3
            nev(i,103) = 1
         else
            addstate(i,103) = 2
            dissocstate(i,103) = 5
            nev(i,103)=1
         endif
      end do
C
C create new events, adding on bottom
c      do i=1,13
c         if (i.eq.1) then
c            addstate(i,97) = 1
c            dissocstate(i,97) = 5
c            nev(i,97) = 1
c         else
c            addstate(i,97) = 1
c            dissocstate(i,97) = 3
c            nev(i,97)=1
c         endif
c      end do
c manually put list of possible events
      totev = 0
      do i=1,13
         j=103
c         do j=97,103,6
            totev = totev+1
            ev(i,j,1) = totev
            evtime(totev) = get_evtime(1,k_on,conc,
     @                       GDP,K_lat,K_long,K_both,k_hyd,
     @                       K_thr_long,K_thr_lat,K_all,k_ne,seammult(i,j))
            event(totev) = totev
            revmap(totev) = totev
            xloc(totev) = i
            yloc(totev) = j
c         end do
      end do

      
      do i=1,totev
         write(6,*) 'evtimes ',i,evtime(event(i)),event(i),revmap(i)
      end do
      
      call build_min_heap(event,revmap,evtime,totev)

      do i=1,totev
         write(6,*) 'evtimes ',i,evtime(event(i)),event(i),revmap(i)
      end do
            
c do fastest event
      done = .false.
      qdebug = .true.
      kk=0
      jj=0
      do while (.not.done)
      kk=kk+1
      
c For getting depolymerization rates...
      if (ntub.eq.350) then		!starting at 350 to allow for hydrolysis to occur & occasional dissociation
         !do once, shouldnt revert back
         write(6,*) 'start shrinking'
         depoly=.true.
         poly=.false.
         conc=1.1E-18   	!conc cannot be zero, error will occur
c         K_both=9.0E-6		!for testing...
         k_hyd=200   		!make everything GDP before dissociating, thus having pure GDP dissociations
         k_ne=1.1E-18		!nucleotide exchange rate cannot be zero, error will occur
      endif
            
      if (ntub.lt.1) then
      	write(6,*) 'zero or negative ntub'
      	stop
      endif

c      do i=1,totev
c         write(6,*) 'evtimes before recalc',i,evtime(event(i)),event(i),revmap(i)
c      end do
            
c      write(6,*) 'crude test, recalc all evtimes'
c Going to count all corner events
      corner_y=0
      corner_n=0
      
      do i=1, totev
c         write(6,*) 'xloc,yloc for event',xloc(event(i)),yloc(event(i)),i
         if (event(i).eq.ev(xloc(event(i)),yloc(event(i)),1)) then
c            write(6,*) 'is addition event'
            evtime(event(i)) = get_evtime(1,k_on,conc,
     @                       GDP,K_lat,K_long,K_both,k_hyd,
     @                       K_thr_long,K_thr_lat,K_all,k_ne,
     @                       seammult(xloc(event(i)),yloc(event(i))))
c Counting all empty subunits     
            if (dissocstate(xloc(event(i)),yloc(event(i))).eq.5) then
               if (xloc(event(i)).eq.1) then
c                  write(6,*) 'looking at row 1'
                  if (state(13,yloc(event(i))+2).ne.0) then
                     if ((dissocstate(13,yloc(event(i))+2).eq.3).or.
     @                  (dissocstate(13,yloc(event(i))+2).eq.4)) then
c                        write(6,*) 'dont add corner for: ',xloc(event(i)),yloc(event(i))
                     else
c                        write(6,*) 'add corner for: ',xloc(event(i)),yloc(event(i))
                        corner_n=corner_n+1
                     endif
                  elseif (state(xloc(event(i))+1,yloc(event(i))).ne.0) then
                     if ((dissocstate(xloc(event(i))+1,yloc(event(i))).eq.3).or.
     @                  (dissocstate(xloc(event(i))+1,yloc(event(i))).eq.4)) then
c                        write(6,*) 'dont add corner for: ',xloc(event(i)),yloc(event(i))
                     else
c                        write(6,*) 'add corner for: ',xloc(event(i)),yloc(event(i))
                        corner_n=corner_n+1
                     endif
                  endif
               elseif (xloc(event(i)).eq.13) then
c                  write(6,*) 'looking at row 13'
                  if (state(xloc(event(i))-1,yloc(event(i))).ne.0) then
                     if ((dissocstate(xloc(event(i))-1,yloc(event(i))).eq.3).or.
     @                  (dissocstate(xloc(event(i))-1,yloc(event(i))).eq.4)) then
c                        write(6,*) 'dont add corner for: ',xloc(event(i)),yloc(event(i))
                     else
c                        write(6,*) 'add corner for: ',xloc(event(i)),yloc(event(i))
                        corner_n=corner_n+1
                     endif
                  elseif (state(1,yloc(event(i))-1).ne.0) then
                     if ((dissocstate(1,yloc(event(i))-1).eq.3).or.
     @                  (dissocstate(1,yloc(event(i))-1).eq.4)) then
c                        write(6,*) 'dont add corner for: ',xloc(event(i)),yloc(event(i))
                     else
c                        write(6,*) 'add corner for: ',xloc(event(i)),yloc(event(i))
                        corner_n=corner_n+1
                     endif
                  endif
               else
c                  write(6,*) 'not looking at seam'
                  if (state(xloc(event(i))-1,yloc(event(i))).ne.0) then
                     if ((dissocstate(xloc(event(i))-1,yloc(event(i))).eq.3).or.
     @                  (dissocstate(xloc(event(i))-1,yloc(event(i))).eq.4)) then
c                        write(6,*) 'dont add corner for: ',xloc(event(i)),yloc(event(i))
                     else
c                        write(6,*) 'add corner for: ',xloc(event(i)),yloc(event(i))
                        corner_n=corner_n+1
                     endif
                  elseif (state(xloc(event(i))+1,yloc(event(i))).ne.0) then
                     if ((dissocstate(xloc(event(i))+1,yloc(event(i))).eq.3).or.
     @                  (dissocstate(xloc(event(i))+1,yloc(event(i))).eq.4)) then
c                        write(6,*) 'dont add corner for: ',xloc(event(i)),yloc(event(i))
                     else
c                        write(6,*) 'add corner for: ',xloc(event(i)),yloc(event(i))
                        corner_n=corner_n+1
                     endif
                  endif
               endif
            endif
         elseif (event(i).eq.ev(xloc(event(i)),yloc(event(i)),2)) then
c            write(6,*) 'is dissoc event'
            evtime(event(i)) = get_evtime(dissocstate(xloc(event(i)),yloc(event(i))),k_on,conc,
     @                       GDP,K_lat,K_long,K_both,k_hyd,
     @                       K_thr_long,K_thr_lat,K_all,k_ne,
     @                       seammult(xloc(event(i)),yloc(event(i))))
c Counting all present subunits     
            if ((dissocstate(xloc(event(i)),yloc(event(i))).eq.5).or.
     @         (dissocstate(xloc(event(i)),yloc(event(i))).eq.6)) then
               corner_y=corner_y+1
            endif
         elseif (event(i).eq.ev(xloc(event(i)),yloc(event(i)),3)) then
c            write(6,*) 'is hydro event'
            evtime(event(i)) = get_evtime(7,k_on,conc,
     @                       GDP,K_lat,K_long,K_both,k_hyd,
     @                       K_thr_long,K_thr_lat,K_all,k_ne,
     @                       seammult(xloc(event(i)),yloc(event(i))))
         elseif (event(i).eq.ev(xloc(event(i)),yloc(event(i)),4)) then								!new for nuc ex
c            write(6,*) 'is nucleotide exchange event'
            evtime(event(i)) = get_evtime(9,k_on,conc,
     @                       GDP,K_lat,K_long,K_both,k_hyd,
     @                       K_thr_long,K_thr_lat,K_all,k_ne,
     @                       seammult(xloc(event(i)),yloc(event(i))))
         else
            write(6,*) 'bad event retrieval',i,event(i),xloc(event(i)),yloc(event(i)),
     @                    ev(xloc(event(i)),yloc(event(i)),1),ev(xloc(event(i)),yloc(event(i)),2),
     @                    ev(xloc(event(i)),yloc(event(i)),3),ev(xloc(event(i)),yloc(event(i)),4)	!new for nuc ex
            stop
         endif
      end do
c Write out totals for corner interactions      
      if (poly) then
         do n=1,4000
            if (kk.eq.n*12.5) then   
               write(6,*) 'corner_y = ',corner_y,'corner_n = ',corner_n
            endif
         enddo
      endif
      
      if (depoly) then
         if ((ntub.le.300).and.(ntub.ge.14)) then
            write(6,*) 'corner_yd = ',corner_y,'corner_nd = ',corner_n
         endif
      endif

      call build_min_heap(event,revmap,evtime,totev)

c      do i=1,totev
c         write(6,*) 'evtimes after recalc',i,evtime(event(i)),event(i),revmap(i)
c      end do

      time = time+evtime(event(1))
      x = xloc(event(1))
      y = yloc(event(1))
      if ( (x.lt.1).or.(x.gt.13).or.(y.lt.1).or.(y.gt.2000) ) then
         write(6,*) 'Grew too wide/tall: x,y= ',x,y
         write(6,*) 'event info: ',event(1)
        write(6,*) 'position of event',xloc(event(1)),yloc(event(1))
        write(6,*) 'event,x,y,ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4)',
     @             event(1),x,y,ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4)				!new for nuc ex
     			do jj=1,totev
     			temp=event(jj)
         		write(6,'(A,I4,F10.6,4I4)') 'times ',
     @ jj,evtime(temp),event(jj),revmap(jj),xloc(temp),yloc(temp)
      			end do
         
         write(6,'(13A)') '  1 ',' 2 ',' 3 ',' 4 ',' 5 ',' 6 ',' 7 ',
     @                     ' 8 ',' 9 ',' 10',' 11',' 12',' 13'
         do ll=120,80,-1
            write(6,'(13I3)') state(1,ll),state(2,ll),state(3,ll),
     @ state(4,ll),state(5,ll),state(6,ll),state(7,ll),state(8,ll),
     @ state(9,ll),state(10,ll),state(11,ll),state(12,ll),state(13,ll)
         end do
         stop
      endif
      
      if (event(1).eq.ev(x,y,1)) then
         if (dissocstate(x,y).ne.2) then
            write(6,'(A,8I4)') 'Adding subunit at ',x,y,ntub,event(1),
     @ ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4),dissocstate(x,y)				!new for nuc ex
            evno = event(1)
            call AddTub(ntub,ontub,evno,xloc,yloc,state,nrecyc,recyc,
     @                  totev,addstate,evtype,nev,ev,maxlen,xo,yo,
     @                  qinit,dissocstate,event,revmap,
     @					evtime,k_on,conc,
     @                        GDP,K_lat,K_long,K_both,k_hyd,
     @                        K_thr_long,K_thr_lat,K_all,dist,qdebug,dep,seammult,k_ne)
c     		    if (qdebug) then
             do n=1,4000
                if (kk.eq.n*500) then
     			do jj=1,totev
     			temp=event(jj)
         		write(6,'(A,I4,F10.6,4I4)') 'times ',
     @ jj,evtime(temp),event(jj),revmap(jj),xloc(temp),yloc(temp)
      			end do
      			endif
      	     enddo
         else
            write(6,*) 'dstate=2, recycle to avoid lateral additions'
            time = time-evtime(event(1))
            evtime(ev(x,y,1))=10.
            call update(event,revmap,ev(x,y,1),evtime(ev(x,y,1)),evtime,totev)
            write(6,*) 'times(ev(x,y,1))',x,y,evtime(ev(x,y,1)),dissocstate(x,y)
            
         endif
      elseif (event(1).eq.ev(x,y,2)) then
         evno = event(1)
         write(6,'(A,7I4)') 'Losing subunit at ',x,y,ntub,event(1),
     @ ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4)							!new for nuc ex
         call Subtub(ntub,ontub,evno,xloc,yloc,state,nrecyc,recyc,
     @      			totev,addstate,nev,ev,maxlen,qinit,xo,yo,
     @ 			dissocstate,test,temp,event,evtime,revmap,k_on,conc,
     @                         GDP,K_lat,K_long,K_both,k_hyd,K_thr_long,
     @          K_thr_lat,K_all,dist,qdebug,dep,seammult,k_ne)
c     		    if (qdebug) then
             do n=1,4000
                if (kk.eq.n*500) then
     			do jj=1,totev
     			temp=event(jj)
         		write(6,'(A,I4,F10.6,4I4)') 'times ',
     @ jj,evtime(temp),event(jj),revmap(jj),xloc(temp),yloc(temp)
      			end do
      			endif
      		 enddo
      elseif (event(1).eq.ev(x,y,3)) then
         write(6,'(A,7I4)') 'Hydrolyzing subunit at ',x,y,ntub,event(1),
     @ ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4)								!new for nuc ex
         evno = event(1)
      call HydroTub(ntub,evno,xloc,yloc,state,nrecyc,recyc,
     @					  totev,addstate,nev,ev,maxlen,
     @ 					dissocstate,dep,k_on,conc,event,revmap,evtime,
     @                             GDP,K_lat,K_long,K_both,k_hyd,K_thr_long,
     @          K_thr_lat,K_all,dist,xo,yo,qdebug,seammult,k_ne)
c hydrolysis event counter
                nhydro=nhydro+1
c     		    if (qdebug) then
             do n=1,4000
                if (kk.eq.n*500) then
     			do jj=1,totev
     			temp=event(jj)
         		write(6,'(A,I4,F10.6,4I4)') 'times ',
     @ jj,evtime(temp),event(jj),revmap(jj),xloc(temp),yloc(temp)
      			end do
      			endif
      		 enddo
      elseif (event(1).eq.ev(x,y,4)) then										!new for nuc ex
         write(6,'(A,7I4)') 'Nucleotide exchange at ',x,y,ntub,event(1),
     @   ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4)
         evno = event(1)
         call NucEx(evno,xloc,yloc,state,nrecyc,recyc,
     @				   totev,nev,ev,maxlen,dissocstate,
     @ 				   event,revmap,evtime)
c         if (qdebug) then
         do n=1,4000
            if (kk.eq.n*500) then
     	       do jj=1,totev
     	          temp=event(jj)
                  write(6,'(A,I4,F10.6,4I4)') 'times ',
     @            jj,evtime(temp),event(jj),revmap(jj),xloc(temp),yloc(temp)
      		   end do
         	endif
         enddo
      	 
c         stop
      else 
        write(6,*) 'you messed up'
        write(6,*) 'position of event',xloc(event(1)),yloc(event(1))
        write(6,*) 'event,x,y,ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4)',
     @             event(1),x,y,ev(x,y,1),ev(x,y,2),ev(x,y,3),ev(x,y,4)				!new for nuc ex
     			do jj=1,totev
     			temp=event(jj)
         		write(6,'(A,I4,F10.6,4I4)') 'times ',
     @ jj,evtime(temp),event(jj),revmap(jj),xloc(temp),yloc(temp)
      			end do
       stop
      endif
      
      write(6,*) 'step number ',kk,'time = ',time,'ntub = ',ntub
      
      dummy=ntub/1625.0
      
      if (depoly) then
         if ((ntub.le.300).and.(ntub.ge.14)) then
            write(6,*) 'wanted2 = ',time,ntub,dummy
            count_d=count_d+1
            if (ntub.eq.300) then
               math_dt=0
               pretm=time
            elseif (ntub.eq.299) then
               math_dt=time-pretm
               pretm=time
            else
               math_dt=math_dt+(time-pretm)
               pretm=time
            endif
            time_d(count_d)=math_dt
            length_d(count_d)=dummy
         endif
      endif      
      
c Microsoft likes no more than 800 pts.
      if (poly) then
         do n=1,4000
            if (kk.eq.n*12.5) then
               write(6,*) 'wanted = ',time,ntub,dummy
               count_p=count_p+1
               time_p(count_p)=time
               length_p(count_p)=dummy
            endif
         enddo
      endif
      
      if (kk.eq.20000) then
         done=.true.
         write(6,*) 'SuperDone'

c For depolymerization rates...
         if (depoly) then
            do count_eqd=1,count_d
               !to get xy...
               tmlngth_d(count_eqd)=time_d(count_eqd)*length_d(count_eqd)
               !to get x^2...
               tm_dsq(count_eqd)=time_d(count_eqd)*time_d(count_eqd)
            enddo
            do count_eqd=1,count_d
               !to get sum(xy)...
               tmlngth_dsm=tmlngth_dsm+tmlngth_d(count_eqd)
               !to get sum(x)...
               tm_dsm=tm_dsm+time_d(count_eqd)
               !to get sum(y)...
               lngth_dsm=lngth_dsm+length_d(count_eqd)
               !to get sum(x^2)...
               tm_dsqsm=tm_dsqsm+tm_dsq(count_eqd)
            enddo
            
            !Rate is given in um/min (since time is in seconds, have to multiply by 60)
            rate_d=(((count_d*tmlngth_dsm)-(tm_dsm*lngth_dsm))/
     @             ((count_d*tm_dsqsm)-(tm_dsm*tm_dsm)))*60.0
            write(6,*) 'rate_d =',rate_d
c         else
c            write(6,*) 'rate_d = null'
         endif
         
c For growth rates...
         if (poly) then
            do count_eqp=1,count_p
               !to get xy
               tmlngth_p(count_eqp)=time_p(count_eqp)*length_p(count_eqp)
               !to get x^2
               tm_psq(count_eqp)=time_p(count_eqp)*time_p(count_eqp)
            enddo
            do count_eqp=1,count_p
               !to get sum(xy)
               tmlngth_psm=tmlngth_psm+tmlngth_p(count_eqp)
               !to get sum(x^2)
               tm_psqsm=tm_psqsm+tm_psq(count_eqp)
            enddo
            
            !Rate is given in um/min (since time is in seconds, have to multiply by 60)
            rate_p=(tmlngth_psm/tm_psqsm)*60.0
            write(6,*) 'rate_p =',rate_p
         endif

c         write(6,*) 'The number of different n-mers are:'
c         write(6,*) 'mono = ',mono
c         write(6,*) 'di = ',di
c         write(6,*) 'tri = ',tri
c         write(6,*) 'tet = ',tet
c         write(6,*) 'pent = ',pent
c         write(6,*) 'hex = ',hex
c         write(6,*) 'hept = ',hept
c         write(6,*) 'oct = ',oct
c         write(6,*) 'ennea = ',ennea
c         write(6,*) 'dec = ',dec
         
c         write(6,*) 'The number of hydrolysis events, nhydro = ',nhydro
      endif
      
      if (poly) then
      do n=1,4000
      if (kk.eq.n*500) then
                write(6,*) 'Printing configuration'
      			write(6,'(13A)') '  1 ',' 2 ',' 3 ',' 4 ',' 5 ',' 6 ',' 7 ',
     @                     ' 8 ',' 9 ',' 10',' 11',' 12',' 13'
            lmax = int(ntub/13)+112
         	do ll=lmax,102,-1
            	write(6,'(13I3)') state(1,ll),state(2,ll),state(3,ll),
     @ 	state(4,ll),state(5,ll),state(6,ll),state(7,ll),state(8,ll),
     @ 	state(9,ll),state(10,ll),state(11,ll),state(12,ll),state(13,ll)
       	 end do
       	 write(6,*) ''
       	  do ll=lmax,102,-1
     		write(6,'(13I3)') dissocstate(1,ll),dissocstate(2,ll),
     @		dissocstate(3,ll),
     @ 		dissocstate(4,ll),dissocstate(5,ll),dissocstate(6,ll),
     @		dissocstate(7,ll),
     @ 		dissocstate(8,ll),dissocstate(9,ll),dissocstate(10,ll),
     @		dissocstate(11,ll),
     @ 		dissocstate(12,ll),dissocstate(13,ll)
          end do
         write(6,*) ''
         	do ll=lmax,102,-1
     		write(6,'(13I3)') addstate(1,ll),addstate(2,ll),
     @		addstate(3,ll),
     @ 	addstate(4,ll),addstate(5,ll),addstate(6,ll),
     @		addstate(7,ll),
     @ 	addstate(8,ll),addstate(9,ll),addstate(10,ll),
     @	addstate(11,ll),
     @ 	addstate(12,ll),addstate(13,ll)
         end do
        write(6,*) ''
         	do ll=lmax,102,-1
     		write(6,'(13I3)') dep(1,ll),dep(2,ll),
     @		dep(3,ll),
     @ 	dep(4,ll),dep(5,ll),dep(6,ll),
     @		dep(7,ll),
     @ 	dep(8,ll),dep(9,ll),dep(10,ll),
     @	dep(11,ll),
     @ 	dep(12,ll),dep(13,ll)
         end do
        write(6,*) ''
        write(6,*) ''
        do ll=lmax,102,-1
     		write(6,'(13F5.1)') seammult(1,ll),seammult(2,ll),
     @		seammult(3,ll),
     @ 		seammult(4,ll),seammult(5,ll),seammult(6,ll),
     @		seammult(7,ll),
     @ 		seammult(8,ll),seammult(9,ll),seammult(10,ll),
     @		seammult(11,ll),
     @ 		seammult(12,ll),seammult(13,ll)
          end do
      endif
      enddo
      endif
      
      if (depoly) then
      if ((ntub.le.300).and.(ntub.ge.14)) then
                write(6,*) 'Printing configuration'
      			write(6,'(13A)') '  1 ',' 2 ',' 3 ',' 4 ',' 5 ',' 6 ',' 7 ',
     @                     ' 8 ',' 9 ',' 10',' 11',' 12',' 13'
            lmax = int(ntub/13)+112
         	do ll=lmax,102,-1
            	write(6,'(13I3)') state(1,ll),state(2,ll),state(3,ll),
     @ 	state(4,ll),state(5,ll),state(6,ll),state(7,ll),state(8,ll),
     @ 	state(9,ll),state(10,ll),state(11,ll),state(12,ll),state(13,ll)
       	 end do
       	 write(6,*) ''
       	  do ll=lmax,102,-1
     		write(6,'(13I3)') dissocstate(1,ll),dissocstate(2,ll),
     @		dissocstate(3,ll),
     @ 		dissocstate(4,ll),dissocstate(5,ll),dissocstate(6,ll),
     @		dissocstate(7,ll),
     @ 		dissocstate(8,ll),dissocstate(9,ll),dissocstate(10,ll),
     @		dissocstate(11,ll),
     @ 		dissocstate(12,ll),dissocstate(13,ll)
          end do
         write(6,*) ''
         	do ll=lmax,102,-1
     		write(6,'(13I3)') addstate(1,ll),addstate(2,ll),
     @		addstate(3,ll),
     @ 	addstate(4,ll),addstate(5,ll),addstate(6,ll),
     @		addstate(7,ll),
     @ 	addstate(8,ll),addstate(9,ll),addstate(10,ll),
     @	addstate(11,ll),
     @ 	addstate(12,ll),addstate(13,ll)
         end do
        write(6,*) ''
         	do ll=lmax,102,-1
     		write(6,'(13I3)') dep(1,ll),dep(2,ll),
     @		dep(3,ll),
     @ 	dep(4,ll),dep(5,ll),dep(6,ll),
     @		dep(7,ll),
     @ 	dep(8,ll),dep(9,ll),dep(10,ll),
     @	dep(11,ll),
     @ 	dep(12,ll),dep(13,ll)
         end do
        write(6,*) ''
        write(6,*) ''
        do ll=lmax,102,-1
     		write(6,'(13F5.1)') seammult(1,ll),seammult(2,ll),
     @		seammult(3,ll),
     @ 		seammult(4,ll),seammult(5,ll),seammult(6,ll),
     @		seammult(7,ll),
     @ 		seammult(8,ll),seammult(9,ll),seammult(10,ll),
     @		seammult(11,ll),
     @ 		seammult(12,ll),seammult(13,ll)
          end do
       endif
       endif      
      
c      write(6,*) 'total number of tubulins',ntub
      end do
      
      stop
      end
