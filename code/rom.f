c-----------------------------------------------------------------------
      subroutine rom_update

      ! set ROM operators and update rom solutions

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld, icalldmhd, z
      save    icalld
      data    icalld /0/

      logical ifmult,iftmp

      parameter (lt=lx1*ly1*lz1*lelt)

      common /romup/ rom_time
      common /mhdflag/ icalldmhd

      stime=dnekclock()

      ! TODO: work for lbb != lub

      if (icalld.eq.0) then
         rom_time=0.
         icalld=1
         icalldmhd=0
         call rom_setup
         if (ifcp) call cp_setup
      endif

      ad_step = istep
      jfield=ifield
      ifield=1

      ifmult=.not.ifrom(2).and.ifheat


      if (rmode.ne.'OFF'.or.rmode.eq.'AEQ') then
      if (ifmult) then
         if (ifflow) call exitti(
     $   'error: running rom_update with ifflow = .true.$',nelv)
         if (istep.gt.0) then
            if (ifrom(ifldb)) then
               call bdfext_step_mhd
            else
               call bdfext_step
            end if
            call post
            call reconv(vx,vy,vz,u) ! reconstruct velocity to be used in h-t
         endif
      else
         if (nio.eq.0) write (6,*) 'starting rom_step loop',ad_nsteps
         ad_step = 1

         cts='bdfext'
c        cts='copt  '
c        cts='rk1   '
c        cts='rkmp  '
c        cts='rk4   '
c        cts='rkck  '

         if (cts.ne.'bdfext'.and.cts.ne.'copt  ') then
            tinit=time
            tfinal=ad_nsteps*ad_dt+time
            ndump=nint(1.*ad_nsteps/ad_iostep)
            idump=1
            tnext=(tfinal-time)/ndump+time
         endif

         time=0
         ad_step=0
         call rom_userchk
         ad_step=1
         do i=1,ad_nsteps

            if (mod(ad_step,iostep).eq.0) then 
               do z=0,nb
                  if (nio.eq.0) write (6,*) 'b',b(z)
               enddo

               do z=0,nb
                  if (nio.eq.0) write (6,*) 'u',u(z)
               enddo
            endif

            if (cts.eq.'bdfext') then
               if (ifrom(ifldb)) then
                  call bdfext_step_mhd
               else
                  call bdfext_step
               end if
               time=time+ad_dt
               call post
            else if (cts.eq.'copt  ') then
               if (ifrom(2)) call rom_step_t_legacy
               if (ifrom(1)) call rom_step_legacy
               call postu_legacy
               call postt_legacy
               time=time+ad_dt
            else
               call rk_setup(cts)
               call copy(urki(1),u(1),nb)
               if (ifrom(2)) call copy(urki(nb+1),ut(1),nb)

               nrk=nb
               if (ifrom(2)) nrk=nb*2

               call rk_step(urko,rtmp1,urki,
     $            time,ad_dt,grk,rtmp2,nrk)

               call mxm(bu,nb,rtmp1,nb,rtmp2,1)
               call mxm(rtmp2,1,rtmp1,nb,err,1)
               err=sqrt(err)

               time=time+ad_dt

               if (nio.eq.0) write (6,1)
     $            ad_step,time,ad_dt,err/ad_dt,err

               if (rktol.ne.0.) then
                  ttmp=ad_dt*.95*((rktol*ad_dt)/err)**.2
                  ad_dt=min(2.*ad_dt,max(.5*ad_dt,ttmp))
               endif

               write (6,*) ad_step,time,tnext,ad_dt,'time'

               if (time*(1.+1.e-12).gt.tnext) then
                  idump=idump+1
                  tnext=(tfinal-tinit)*idump/ndump+tinit

                  call reconv(vx,vy,vz,u)
                  call recont(t,ut)
                  if (ifrom(ifldb)) call reconb(bx,by,bz,b)

                  ifto = .true. ! turn on temp in fld file
                  if (rmode.ne.'ON ') then
                     call outpost(vx,vy,vz,pr,t,'rom')
                     if (ifrom(ifldb)) then
                        write(*,*) 'ROM HERE', ''
                        call outpost(bx,by,bz,pr,t,'romb')
                     endif
                  endif
               endif

               if (cts.eq.'rkck'.and.rktol.ne.0.) then
               if (time+ad_dt.gt.(1.+1.e-12)*tnext) then
                  ttmp=ad_dt
                  ad_dt=tnext-time
                  write (6,3) ad_step,time,ad_dt,ttmp,tnext
               endif
               endif

               call copy(u(1),urko,nb)
               call copy(ut(1),urko(nb+1),nb)

               if (time*(1.+1.e-14).gt.tfinal) goto 2
            endif

            ad_step=ad_step+1
         enddo
         icalld=0
      endif
      endif

    2 continue

      if (ifei.and.(rmode.ne.'OFF'.and.rmode.ne.'AEQ')) then
         call cdres
c        call cres
      endif

      ifield=jfield

      dtime=dnekclock()-stime
      rom_time=rom_time+dtime

      if (ifmult.and.nio.eq.0) write (6,*) 'romd_time: ',dtime
      if (.not.ifmult.or.nsteps.eq.istep) call final

    1 format(i8,1p4e13.5,' err')
    3 format(i8,1p4e15.7,' modt')

      do i=0,nb
         if (nio.eq.0) write (6,*) 'b',b(i)
      enddo

      do i=0,nb
         if (nio.eq.0) write (6,*) 'u',u(i)
      enddo



      return
      end
c-----------------------------------------------------------------------
      subroutine offline_mode

      ! offline-wrapper for MOR

      include 'SIZE'
      include 'INPUT'

      param(173)=1.
      call rom_update

      return
      end
c-----------------------------------------------------------------------
      subroutine online_mode

      ! online-wrapper for MOR

      include 'SIZE'
      include 'INPUT'

      param(173)=2.
      call rom_update

      return
      end
c-----------------------------------------------------------------------
      subroutine recon_mode

      ! reconstruction online-wrapper for MOR

      include 'SIZE'
      include 'INPUT'

      param(173)=3.
      call rom_update

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup_params

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'AVG'
      include 'INPUT'

      logical iftmp,ifexist
      integer num_ts

      num_ts=1

      

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup_mhd

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'AVG'
      include 'INPUT'

      common /mhdflag/ icalldmhd

      if (nio.eq.0) write (6,*) 'inside rom_setup_mhd'
      jfield = ifield
      ifield=ifldmhd

      icalldmhd=1

      !TODO check on this one
c      call mor_set_params_uni_pre
      call mor_init_fields
      call mor_set_params_uni_post

      if (ifsetbases) call setbases
      call rom_userbases


      ubdim = lx1*ly1*lz1*lelm

      ! TODO: fix copy

      call copy(bs0,us0,lx1*ly1*lz1*lelm*ldim*lsu)

      do i=1,ubdim
      do j=0,lub
         bxb(i,j) = ub(i,j)
         byb(i,j) = vb(i,j)
         bzb(i,j) = wb(i,j)
      do k=0,ldim
         bxyzb(i,k,j) = uvwb(i,k,j)
      enddo
      enddo
      enddo
      
      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         call dump_bas
      endif

      
      icalldmhd=2
      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup

      ! set rom ic, ops, qoi, etc.

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'AVG'
      include 'INPUT'

      logical iftmp,ifexist
      integer num_ts, ubdim

      num_ts=1

      if (nio.eq.0) write (6,*) 'inside rom_setup'

      call nekgsync
      setup_start=dnekclock()

      ttime=time

      n=lx1*ly1*lz1*nelt

      call set_bdf_coef(ad_alpha,ad_beta)
      call mor_init_params

      if (param(170).lt.0) then
         call mor_set_params_par
      else
         call mor_set_params_rea
      endif

      if (nid.eq.0) write (6,*) 'post par rea'
      call mor_show_params
      if (nid.eq.0) write (6,*) 'post show params'
      
      call mor_set_params_uni_pre

      if (ifrom(ifldb)) call rom_setup_mhd

      call mor_init_fields
      call mor_set_params_uni_post

      if (ifsetbases) call setbases
      call rom_userbases

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         call dump_bas
      endif

c     call average_in_xy
c     call average_in_y
      call setops
      call setu
      if (ifrom(ifldb)) call setb_ic

      call setf

c      call setqoi
      call setmisc

      if (ifei) then
         call set_sigma
      endif

      if (nio.eq.0) write (6,*) 'end rom setup'

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $   call dump_misc

      time=ttime

      call nekgsync
      setup_end=dnekclock()

      if (nio.eq.0) write (6,*) 'exiting rom_setup'
      if (nio.eq.0) write (6,*) 'setup_time:', setup_end-setup_start

      return
      end
c-----------------------------------------------------------------------
      subroutine asnap(uuk,ttk)

      ! averaging of coefficients obtained from snapshots

      ! uuk := velocity coefficient set
      ! ttk := velocity coefficient set

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      common /scrasnap/ t1(0:lub)

      real uuk(0:nb,ns),ttk(0:nb,ns)

      call nekgsync
      asnap_time=dnekclock()

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (ifrom(1)) then
            call read_serial(uas,nb+1,'ops/uas ',t1,nid)
            call read_serial(uvs,nb+1,'ops/uvs ',t1,nid)
         endif
         if (ifrom(2)) then
            call read_serial(tas,nb+1,'ops/tas ',t1,nid)
            call read_serial(tvs,nb+1,'ops/tvs ',t1,nid)
         endif
      else
         s=1./real(ns)

         if (ifrom(1)) then
            call pv2b(uas,uavg,vavg,wavg,ub,vb,wb)
            call rzero(uvs,nb+1)
            do j=1,ns
               do i=0,nb
                  uvs(i)=uvs(i)+(uuk(i,j)-uas(i))**2
               enddo
            enddo
            call cmult(uvs,s,nb+1)

            call dump_serial(uas,nb+1,'ops/uas ',nid)
            call dump_serial(uvs,nb+1,'ops/uvs ',nid)
         endif

         if (ifrom(2)) then
            call ps2b(tas,tavg,tb)
            call rzero(tvs,nb+1)
            do j=1,ns
               do i=0,nb
                  tvs(i)=tvs(i)+(ttk(i,j)-tas(i))**2
               enddo
            enddo
            call cmult(tvs,s,nb+1)

            call dump_serial(tas,nb+1,'ops/tas ',nid)
            call dump_serial(tvs,nb+1,'ops/tvs ',nid)
         endif
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'asnap_time:',dnekclock()-asnap_time

      return
      end
c-----------------------------------------------------------------------
      subroutine setops

      ! set rom operators

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      common /tmpop/ ftmp0(0:lb,0:lb)

      real visc(10)

      logical iftmp

      if (nio.eq.0) write (6,*) 'inside setops'

      call nekgsync
      ops_time=dnekclock()

      jfield=ifield
      if (ifrom(1)) then
         ifield=1
         call seta(au,au0,'ops/au ')
         if (ifrom(ifldb)) then
            call seta_mhd(abu,ab0,'ops/ab ')
         endif
         if (rmode.eq.'AEQ') call setae(aue,'ops/aue ')
         call setb(bu,bu0,'ops/bu ')
         if (ifrom(ifldb)) then
            call setb_mhd(bbu,bb0,'ops/bb ')
         endif
         call setc(cul,'ops/cu ')
         if (ifrom(ifldb)) then
            call setc_mhd(cul,cbl,cubl,cbul
     $                    ,'ops/cuu', 'ops/cbb', 'ops/cub', 'ops/cbu')
            if (ifmhdfixskew) call fix_skew_mhd
            call dump_serial(cul,ncloc,'ops/cuu ',nid)
            call dump_serial(cbl,ncloc,'ops/cbb ',nid)
            call dump_serial(cbul,ncloc,'ops/cbu ',nid)
            call dump_serial(cubl,ncloc,'ops/cub ',nid)
         endif
      endif
      if (ifrom(2)) then
         ifield=2
         call seta(at,at0,'ops/at ')
         if (rmode.eq.'AEQ') call setae(ate,'ops/ate ')
         call setb(bt,bt0,'ops/bt ')
         call setc(ctl,'ops/ct ')
         call sets(st0,tb,'ops/ct ')
      endif
      if (ifedvs) then
         call seteddy(cedd,'ops/cedd ')
      endif

      if (rmode.eq.'AEQ') call setfluc(fv_op,ft_op,'fluc')

      if (ifbuoy.and.ifrom(1).and.ifrom(2))
     $   call setbut(buxt0,buyt0,buzt0)

      ifield=jfield

      call set_binv2

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $   call dump_ops

      call nekgsync
      if (nio.eq.0) write (6,*) 'ops_time:',dnekclock()-ops_time

      if (nio.eq.0) write (6,*) 'exiting setops'

      return
      end
c-----------------------------------------------------------------------
      subroutine setqoi

      ! set qoi factors

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk2/ a1(lt),a2(lt),a3(lt),wk(lub+1)

      if (nio.eq.0) write (6,*) 'begin setup for qoi'

      call cvdrag_setup
      call cnuss_setup
      call cubar_setup

      if (ifcdrag) then
         if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
            call read_serial(fd1,ldim*(nb+1),'qoi/fd1 ',wk,nid)
            call read_serial(fd3,ldim*(nb+1),'qoi/fd3 ',wk,nid)
         else
            do i=0,nb
               call lap2d(a1,ub(1,i))
               call lap2d(a2,vb(1,i))
               if (ldim.eq.3) call lap2d(a3,wb(1,i))
               call outpost(a1,a2,a3,a3,tb,'lap')
               call cint(fd3(1+ldim*i),a1,a2,a3)
            enddo
         endif
      endif

      if (nintp.gt.0) then
         m=0
         if (nid.eq.0) m=nintp
         do i=0,nb
            call gfldi(tbintp(1+nintp*i),tb(1,i,1),
     $                 xintp,yintp,zintp,m,1)
         enddo
         call dump_serial(xintp,nintp*(nb+1),'qoi/xintp',nid)
         call dump_serial(yintp,nintp*(nb+1),'qoi/yintp',nid)
         call dump_serial(zintp,nintp*(nb+1),'qoi/zintp',nid)
         call dump_serial(tbintp,nintp*(nb+1),'qoi/tintp',nid)
      endif

      if (cftype.eq.'POLY') then
         call set_les_imp(aules,atles)
      else if (cftype.eq.'SMAG') then
         call exitti('invalid filter type SMAG$',1)
      endif

      if (nio.eq.0) write (6,*) 'end setup for qoi'

      return
      end
c-----------------------------------------------------------------------
      subroutine update_k(uuk,uukp,ttk,ttkp)

      ! update snapshot coefficients with transformation wt

      ! uuk  := original velocity coefficient set
      ! uukp := new velocity coefficient set
      ! ttk  := original temperature coefficient set
      ! ttkp := new temperature coefficient set

      include 'SIZE'
      include 'MOR'

      real uuk(0:nb,ns),uukp(0:nb,ns),ttk(0:nb,ns),ttkp(0:nb,ns)

      if (nio.eq.0) write (6,*) 'begin update setup'

      call nekgsync
      kupd_time=dnekclock()

      if (ifpod(1)) then
         do i=1,ns
            call mxm(wt,nb,uuk(1,i),nb,uukp(1,i),1)
         enddo
      endif

      if (ifpod(2)) then
         do i=1,ns
            call mxm(wt(1,2),nb,ttk(1,i),nb,ttkp(1,i),1)
         enddo
      endif

      call nekgsync
      if (nio.eq.0) write (6,*) 'kupd_time:',dnekclock()-kupd_time

c     call hyperpar
      call update_hyper(ukp,tkp)

      return
      end
c-----------------------------------------------------------------------
      subroutine update_hyper(uukp,ttkp)

      ! update hyper parameters

      ! uukp := transformed velocity coefficient set
      ! ttkp := transformed thermal coefficient set

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uukp(0:nb,ns),ttkp(0:nb,ns)

      real ep

      call nekgsync
      hpar_time=dnekclock()

      ! eps is the free parameter
      ! 1e-2 is used in the paper
      ep = 1.e-2

      n  = lx1*ly1*lz1*nelt
      if (ifpod(1)) then
         call cfill(upmin,1.e9,nb)
         call cfill(upmax,-1.e9,nb)
         do j=1,ns
         do i=1,nb
            if (uukp(i,j).lt.upmin(i)) upmin(i)=uukp(i,j)
            if (uukp(i,j).gt.upmax(i)) upmax(i)=uukp(i,j)
         enddo
         enddo
         do j=1,nb                    ! compute hyper-parameter
            d= upmax(j)-upmin(j)
            upmin(j) = upmin(j) - ep * d
            upmax(j) = upmax(j) + ep * d
            if (nio.eq.0) write (6,*) j,upmin(j),upmax(j),'upminmax'
         enddo

         ! compute distance between umax and umin
         call sub3(updis,upmax,upmin,nb)

         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,updis(i),'updis'
            enddo
         endif
      endif   

      if (ifpod(2)) then
         call cfill(tpmin,1.e9,nb)
         call cfill(tpmax,-1.e9,nb)
         do j=1,ns
         do i=1,nb
            if (ttkp(i,j).lt.tpmin(i)) tpmin(i)=ttkp(i,j)
            if (ttkp(i,j).gt.tpmax(i)) tpmax(i)=ttkp(i,j)
         enddo
         enddo
         do j=1,nb                    ! compute hyper-parameter
            d= tpmax(j)-tpmin(j)
            tpmin(j) = tpmin(j) - ep * d
            tpmax(j) = tpmax(j) + ep * d
            if (nio.eq.0) write (6,*) j,tpmin(j),tpmax(j),'tpminmax'
         enddo

         ! compute distance between tmax and tmin
         call sub3(tpdis,tpmax,tpmin,nb)
         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,tpdis(i),'tpdis'
            enddo
         endif

      endif   

      call nekgsync
      if (nio.eq.0) write (6,*) 'hpar_time',dnekclock()-hpar_time

      return
      end
c-----------------------------------------------------------------------
      subroutine update_hyper_mhd(bbkp)

      ! update hyper parameters

      ! bbkp := transformed velocity coefficient set

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real bbkp(0:nb,ns)

      real ep

      call nekgsync
      hpar_time=dnekclock()

      ! eps is the free parameter
      ! 1e-2 is used in the paper
      ep = 1.e-2

      n  = lx1*ly1*lz1*nelt
      if (ifpod(ifldb)) then
         call cfill(bpmin,1.e9,nb)
         call cfill(bpmax,-1.e9,nb)
         do j=1,ns
         do i=1,nb
            if (bbkp(i,j).lt.bpmin(i)) bpmin(i)=bbkp(i,j)
            if (bbkp(i,j).gt.bpmax(i)) bpmax(i)=bbkp(i,j)
         enddo
         enddo
         do j=1,nb                    ! compute hyper-parameter
            d= bpmax(j)-bpmin(j)
            bpmin(j) = bpmin(j) - ep * d
            bpmax(j) = bpmax(j) + ep * d
            if (nio.eq.0) write (6,*) j,bpmin(j),bpmax(j),'bpminmax'
         enddo

         ! compute distance between bmax and bmin
         call sub3(bpdis,bpmax,bpmin,nb)

         if (nio.eq.0) then
            do i=1,nb
               write (6,*) i,bpdis(i),'bpdis'
            enddo
         endif
      endif   

      call nekgsync
      if (nio.eq.0) write (6,*) 'hpar_time',dnekclock()-hpar_time

      return
      end
c-----------------------------------------------------------------------
      subroutine setmisc

      ! set miscellaneous quantities

      include 'SIZE'
      include 'MOR'

      logical ifexist

      if (nio.eq.0) write (6,*) 'begin range setup'

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         call nekgsync
         proj_time=dnekclock()

         if (ifpod(1)) call p2k(uk,us0,uvwb,ndim,ukp)
         if (ifpod(2)) call p2k(tk,ts0(1,1,1),tb(1,0,1),1,tkp)
         if (ifpod(ifldb))  call p2k(bk,bs0,bxyzb,ndim,bkp)
         if (ifedvs)   call p2k(edk,ts0(1,1,4),tb(1,0,4),1,ukp)

         call nekgsync
         if (nio.eq.0) write (6,*) 'proj_time:',dnekclock()-proj_time
      else if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP') then
         inquire (file='ops/uk',exist=ifexist)
         if (ifexist)
     $      call read_mat_serial(uk,nb+1,ns,'ops/uk ',mb+1,nns,stmp,nid)

         inquire (file='ops/tk',exist=ifexist)
         if (ifexist)
     $      call read_mat_serial(tk,nb+1,ns,'ops/tk ',mb+1,nns,stmp,nid)
         
         inquire (file='ops/bk',exist=ifexist)
         if (ifexist)
     $      call read_mat_serial(bk,nb+1,ns,'ops/bk ',mb+1,nns,stmp,nid)
     
      endif

      if (ifpod(1)) call copy(ukp,uk,(nb+1)*ns)
      if (ifpod(2)) call copy(tkp,tk,(nb+1)*ns)
      if (ifpod(ifldb)) call copy(bkp,bk,(nb+1)*ns)

      !TODO check what this is
      call asnap(uk,tk)
      call hyperpar(uk,tk)

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_init_params

      ! initialize rom parameters before .rea / .par read

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      character*3 chartmp

      if (nio.eq.0) write (6,*) 'inside mor_init_params'

      ad_nsteps=nsteps
      ad_iostep=iostep

      ubarr0=1e-1
      ubarrseq=5
      tbarr0=1e-1
      tbarrseq=5

      ustep_time = 0.
      solve_time=0.
      lu_time=0.
      ucopt_count=0
      copt_time=0.
      quasi_time=0.
      lnsrch_time=0.
      compgf_time=0.
      compf_time=0.
      invhm_time=0.

      anum_galu=0.
      anum_galt=0.

      ad_dt = dt
      ad_re = 1/param(2)
      ad_pe = 1/param(8)
      ad_mag = 1/param(29)

      ifavg0=.false.

      ifctke=.false.
      ifcdrag=.false.
      iffastc=.false.
      iffasth=.false.
      ifcintp=.false.
      ifavisc=.false.
      ifdecpl=.false.
      ifsub0=.true.
      ifcomb=.false.
      ifpb=.true.
      ifcp=.false.
      ifcore=.true.
      ifquad=.false.
      ifsetbases=.true.

      do i=0,ldimt3
         ifpod(i)=.false.
         ifrom(i)=.false.
      enddo

      ifvort=.false. ! default to false for now

      ifpart=.false.
      ifcintp=.false.
      ifcflow=.false.

      ips='L2 '
      isolve=0
      rmode='OFF'
      ifei=.false.
      navg_step=1
      nb=lb
      ns=ls
      ntr=ltr
      nskip=0
      rktol=0.
      ad_qstep=ad_iostep
      iftneu=.false.
      inus=0
      iaug=0
      ifsub0=.true.

      ifforce=.false.
      ifsource=.false.
      ifbuoy=.false.
      ifplay=.false.
      nplay=0

      icopt=0

      podrat=0.5

      cfloc='NONE'
      cftype='NONE'
      rbf=0.5
      rdft=0.5

      gx=0.
      gy=0.
      gz=0.

      nintp=0
      nbat=200

      if (nio.eq.0) write (6,*) 'exiting mor_init_params'

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_set_params_rea

      ! set rom parameters with .rea

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      character*3 chartmp

      isolve=nint(param(170))

      if (param(171).eq.0) then
         ips='L2 '
      else if (param(171).eq.1) then
         ips='H10'
      else
         ips='HLM'
      endif

      if (param(172).ne.0) ifavg0=.true.

      np173=nint(param(173))

      if (np173.eq.0) then
         rmode='ALL'
      else if (np173.eq.1) then
         rmode='OFF'
      else if (np173.eq.2) then
         rmode='ON '
      else if (np173.eq.3) then
         rmode='ONB'
         call exitti('unsupported param(173), exiting...$',np173)
      endif

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         open (unit=10,file='ops/ips')
         read (10,*) chartmp
         close (unit=10)
         if (chartmp.ne.ips)
     $      call exitti('online ips does not match offline ips$',nb)
      endif

      ifrecon=(rmode.ne.'ON '.and.rmode.ne.'CP ')

      ifei=nint(param(175)).ne.0
      if (nint(param(175)).eq.0) then
         continue
      else if (nint(param(175)).eq.1) then
         eqn='NS  '
      else if (nint(param(175)).eq.2) then
         eqn='SNS '
      else if (nint(param(175)).eq.3) then
         eqn='POIS'
      else if (nint(param(175)).eq.4) then
         eqn='HEAT'
      else if (nint(param(175)).eq.5) then
         eqn='ADVE'
      else if (nint(param(175)).eq.6) then
         eqn='VNS '
      else if (nint(param(175)).eq.7) then
         eqn='VSNS'
      else
         call exitti('invalid option for eqn$',eqn)
      endif

      navg_step=nint(max(1.,param(176)))
      do i=1,3
         its(i)=min(navg_step+i-1,3)
      enddo

      nb=min(nint(param(177)),lb)
      if (nb.eq.0) nb=lb

      if (ifavg0.and.(nb.eq.ls))
     $   call exitti('nb == ls results in linear dependent bases$',nb)

      if (nb.gt.ls)
     $   call exitti('nb > ls is undefined configuration$',nb)

      rktol=param(179)
      if (rktol.lt.0.) rktol=10.**rktol

      ad_qstep=nint(param(180))+ad_iostep*max(1-nint(param(180)),0)

      iftneu=(param(178).ne.0.and.param(174).ne.0)

      if (param(181).ne.0) ifctke=.true.

      if (param(182).ne.0) ifcdrag=.true.

      inus=min(max(nint(param(183)),0),5)

      if (param(191).ne.0) iffastc=.true.

      if (param(192).ne.0.and.ips.eq.'HLM') iffasth=.true.

      if (param(197).ne.0.) ifsub0=.false.

      ifpod(1)=param(174).ge.0.
      ifpod(2)=(ifheat.and.param(174).ne.0.)
c     ifrom(1)=(ifpod(1).and.eqn.ne.'ADE')
      ifrom(1)=ifpod(1)
      ifrom(2)=ifpod(2)
      ifrom(ifldb)=ifpod(ifldb)

      ifpod(1)=ifpod(1).or.ifrom(2)

      ifforce=param(193).ne.0.
      ifsource=param(194).ne.0.
      ifbuoy=param(195).ne.0.

      nplay=nint(param(196))

      if (nplay.ne.0) then
         nplay=max(nplay,0)
         ifplay=.true.
      else
         ifplay=.false.
      endif

      ! constrained optimization parameter
      icopt=nint(param(184))
      barr_func=param(185)
      box_tol=param(186)

      ! barrier initial parameter and number of loop
      ubarr0  =param(187)
      ubarrseq=nint(param(188))
      tbarr0  =param(189)
      tbarrseq=nint(param(190))

      ! filter technique
      np198=nint(param(198))
      if (np198.eq.0) then
         cfloc='NONE'
         cftype='NONE'
      else if (np198.eq.1) then
         cfloc='CONV'
         cftype='TFUN'
      else if (np198.eq.2) then
         cfloc='POST'
         cftype='DIFF'
      else
         call exitti('unsupported param(198), exiting...$',np198)
      endif

      ! POD spatial filter
      if (param(199).ge.0) then
         rbf=min(nint(param(199)),nb-1)
      else
         rbf=nint(nb*min(1.,-param(199)))
      endif

      ! POD radius of the filter for differential filter
      rdft=param(200)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         rtmp1(1,1)=nb*1.
         call dump_serial(rtmp1(1,1),1,'ops/nb ',nid)
      else
         call read_serial(rtmp1(1,1),1,'ops/nb ',b,nid)
         mb=rtmp1(1,1)
         if (mb.lt.nb) call exitti('mb less than nb, exiting...$',mb)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_set_params_uni_pre

      ! set parameters universal to .rea and .par before mor_init_fields

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      nfbc=0

      nv=lx1*ly1*lz1*nelv

      do ie=1,nelt
      do ifc=1,2*ldim
         if (cbc(ifc,ie,2).eq.'f  ') nfbc=nfbc+1
      enddo
      enddo

      nfbc=iglsum(nfbc,1)

      if (nfbc.gt.0) iftflux=.true.
      if (iftflux) iftneu=.true.

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_set_params_uni_post

      ! set parameters universal to .rea and .par after mor_init_fields

      include 'SIZE'
      include 'MASS'
      include 'MOR'

      nv=lx1*ly1*lz1*nelv

      ux=glsc2(ub,bm1,nv)
      uy=glsc2(vb,bm1,nv)
      if (ldim.eq.3) then
         uz=glsc2(wb,bm1,nv)
      else
         uz=0.
      endif

      if (abs(ux).gt.max(abs(uy),abs(uz))) idirf=1
      if (abs(uy).gt.max(abs(ux),abs(uz))) idirf=2
      if (abs(uz).gt.max(abs(ux),abs(uy))) idirf=3

      if (nio.eq.0) write (6,*) idirf,ux,uy,uz,' idirf'

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_show_params

      ! dump out parameter settings

      include 'SIZE'
      include 'MOR'

      if (nid.eq.0) write (6,*) 'inside mor show params'

      if (nio.eq.0) then
         write (6,*) 'mp_nb         ',nb
         write (6,*) 'mp_lub        ',lub
         write (6,*) 'mp_ltb        ',ltb
         write (6,*) ' '
         write (6,*) 'mp_ls         ',ls
         write (6,*) 'mp_lsu        ',lsu
         write (6,*) 'mp_lst        ',lst
         write (6,*) 'mp_ltr        ',ltr
         write (6,*) 'mp_ns         ',ns
         write (6,*) 'mp_nskip      ',nskip
         write (6,*) ' '
         write (6,*) 'mp_ad_re      ',ad_re
         write (6,*) 'mp_ad_pe      ',ad_pe
         write (6,*) ' '
         write (6,*) 'mp_isolve     ',isolve
         write (6,*) 'mp_ips        ',ips
         write (6,*) 'mp_ifavg0     ',ifavg0
         write (6,*) 'mp_ifsub0     ',ifsub0
         write (6,*) 'mp_rmode      ',rmode
         write (6,*) 'mp_ad_qstep   ',ad_qstep
         write (6,*) 'mp_ifctke     ',ifctke
         write (6,*) 'mp_ifcdrag    ',ifcdrag
         write (6,*) 'mp_inus       ',inus
         write (6,*) 'mp_iffastc    ',iffastc
         write (6,*) 'mp_iffasth    ',iffasth
         write (6,*) 'mp_nplay      ',nplay
         write (6,*) 'mp_ifplay     ',ifplay
         write (6,*) 'mp_ifei       ',ifei
         write (6,*) 'mp_iaug       ',iaug
         write (6,*) 'mp_ifedvs     ',ifedvs
         write (6,*) ' '
         write (6,*) 'mp_ifforce    ',ifforce
         write (6,*) 'mp_ifsource   ',ifsource
         write (6,*) 'mp_ifbuoy     ',ifbuoy
         write (6,*) 'mp_ifpart     ',ifpart
         write (6,*) 'mp_ifvort     ',ifvort
         write (6,*) 'mp_ifcintp    ',ifcintp
         write (6,*) 'mp_ifdecpl    ',ifdecpl
         write (6,*) 'mp_ifcp       ',ifcp
         write (6,*) 'mp_ifcore     ',ifcore
         write (6,*) 'mp_ifquad     ',ifquad
         write (6,*) ' '
         do i=0,ldimt3
            write (6,*) 'mp_ifpod(',i,')   ',ifpod(i)
         enddo
         do i=0,ldimt3
            write (6,*) 'mp_ifrom(',i,')   ',ifrom(i)
         enddo
         write (6,*) ' '
         write (6,*) 'mp_gx          ',gx
         write (6,*) 'mp_gy          ',gy
         write (6,*) 'mp_gz          ',gz
         write (6,*) ' '
         write (6,*) 'mp_icopt       ',icopt
         write (6,*) 'mp_barr_func   ',barr_func
         write (6,*) 'mp_box_tol     ',box_tol
         write (6,*) 'mp_ubarr0      ',ubarr0
         write (6,*) 'mp_ubarrseq    ',ubarrseq
         write (6,*) 'mp_tbarr0      ',tbarr0
         write (6,*) 'mp_tbarrseq    ',tbarrseq
         write (6,*) ' '
         write (6,*) 'mp_cfloc       ',cfloc
         write (6,*) 'mp_cftype      ',cftype
         write (6,*) 'mp_rbf         ',rbf
         write (6,*) 'mp_rdft        ',rdft
         write (6,*) ' '
         write (6,*) 'mp_navg_step   ',navg_step
         write (6,*) 'mp_rk_tol      ',rk_tol
         write (6,*) 'mp_iftneu      ',iftneu
         write (6,*) 'mhd do cross terms      ',ifmhdcrossterms
         write (6,*) 'mhd do fix skew      ',ifmhdfixskew
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mor_init_fields

      ! initialize FOM fields needed for ROM generation

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical alist,iftmp

      character*128 fname1
      
      common /mhdflag/ icalldmhd

      if (nio.eq.0) write (6,*) 'inside rom_init_fields'

      n=lx1*ly1*lz1*nelt
      jfield=ifield
      ifield=1

      call rone(wm1,n)
      call rone(ones,n)
      call rzero(zeros,n)

      if (.not.iftneu) call rzero(gn,n)

      call rzero(num_galu,nb)
      call rzero(num_galt,nb)

      if (ifrom(1)) call opcopy(uic,vic,wic,vx,vy,vz)
      if (ifrom(2)) call copy(tic,t,n)
      if(ifrom(ifldb)) call opcopy(bxic,byic,bzic,bx,by,bz)
      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         fname1='file.list'
          if (ifrom(ifldb).and.icalldmhd.eq.1) then
            fname1='fileb.list'
         endif

         do i=0,ldimt3
            ifreads(i)=ifrom(i)
            if (nio.eq.0) write(6,*)'ifreads',ifreads(i)
         enddo

         call read_fields(
     $      us0,prs,ts0,ns,ls,nskip,ifreads,timek,fname1,.true.)
         
         !TODO: can maybe read b fields in here too

         rtmp1(1,1)=1.0*ns
         call dump_serial(rtmp1(1,1),1,'ops/ns ',nid)

         fname1='avg.list'
         inquire (file=fname1,exist=alist)
         if (alist) then
            call push_sol(vx,vy,vz,pr,t)
            ttime=time
            call auto_averager(fname1)
            time=ttime
            call copy_sol(uavg,vavg,wavg,pavg,tavg,vx,vy,vz,pr,t)
            call pop_sol(vx,vy,vz,pr,t)
         endif

         if (ifavg0) then
            idc_u=0
            idc_t=0
            do ie=1,nelt
            do ifc=1,2*ldim
               if (cbc(ifc,ie,1).eq.'W  '.or.
     $             cbc(ifc,ie,1).eq.'V  '.or.
     $             cbc(ifc,ie,1).eq.'v  ') idc_u=idc_u+1
               if (cbc(ifc,ie,2).eq.'T  '.or.
     $             cbc(ifc,ie,2).eq.'t  ') idc_t=idc_t+1
            enddo
            enddo
            idc_u=iglsum(idc_u,1)
            idc_t=iglsum(idc_t,1)

            call opcopy(ub,vb,wb,uavg,vavg,wavg)

            n2=lx2*ly2*lz2*nelv
            call copy(pb,pavg,n2)
            do j=1,npscal+1
               call copy(tb(1,0,j),tavg(1,1,1,1,j),n)
            enddo

            call opcopy(uvwb(1,1,0),uvwb(1,2,0),uvwb(1,ldim,0),
     $         ub,vb,wb)
c           if (idc_u.gt.0) call opzero(ub,vb,wb)
c           if (idc_t.gt.0) call rzero(tb,n)
         endif

         iftmp=ifxyo
         ifxyo=.true.
         if (ifrecon) then
            call outpost(uavg,vavg,wavg,pavg,tavg,'avg')
            call outpost(urms,vrms,wrms,prms,trms,'rms')
            call outpost(vwms,wums,uvms,prms,trms,'tmn')
         endif
         ifxyo=iftmp

         if (nio.eq.0) write (6,*) 'pre 0 mean'

c        if (ifrom(2).and.iftflux) then
c           call set0mean(tb)
c           do i=1,ns
c              call set0mean(ts0(1,i))
c           enddo
c           call set1dudn(tb)
c           do i=1,ns
c              call set1dudn(ts0(1,i))
c           enddo
c        endif

         if (nio.eq.0) write (6,*) 'post 0 mean'

         if (ifsub0) then
            do i=1,ns
               if (ifrom(1)) then
                  call sub2(us0(1,1,i),ub,n)
                  call sub2(us0(1,2,i),vb,n)
                  if (ldim.eq.3) call sub2(us0(1,ldim,i),wb,n)
               endif
              if (ifrom(ifldb)) then
                 call sub2(bs0(1,1,i),bxb,n)
                 call sub2(bs0(1,2,i),byb,n)
                 if (ldim.eq.3) call sub2(bs0(1,ldim,i),bzb,n)
              endif 
               if (ifrom(2)) call sub2(ts0(1,i,1),tb(1,0,1),n)
               if (ifedvs) call sub2(ts0(1,i,4),tb(1,0,4),n)
            enddo
            call sub2(uavg,ub,n)
            call sub2(vavg,vb,n)
            if (ldim.eq.3) call sub2(wavg,wb,n)
            if (ifrom(ifldb)) then
               call sub2(bxavg,bxb,n) ! these are not yet a thing
               call sub2(byavg,byb,n)
               if (ldim.eq.3) call sub2(bzavg,bzb,n)
            endif
            if (ifpod(2)) call sub2(tavg,tb(1,0,1),n)
            if (ifedvs) call sub2(tavg(1,1,1,1,4),tb(1,0,4),n)
         endif
      endif

      if (rmode.eq.'AEQ') then
         iftmp=.false.
         if (ifpod(1)) then
            fname1='uavg.list '
            ifreads(1)=.true.
            ifreads(2)=.true.
            ifreads(3)=.false.
            call read_fields(uafld,pafld,fldtmp,
     $         nbavg,lbavg,0,ifreads,1,tk,fname1,iftmp)

            fname1='urms.list'
            ifreads(2)=.false.
            call read_fields(uufld,fldtmp,fldtmp,
     $         nbavg,lbavg,0,ifreads,tk,fname1,iftmp)

            fname1='urm2.list'
            call read_fields(uvfld,fldtmp,fldtmp,
     $         nbavg,lbavg,0,ifreads,tk,fname1,iftmp)

            ifxyo=.true.

            do i=1,nbavg
               call setupvp(upup(1,1,i),
     $            uufld(1,1,i),uafld(1,1,i),uafld(1,1,i))
               call setupvp(upup(1,2,i),
     $            uufld(1,2,i),uafld(1,2,i),uafld(1,2,i))
               if (ldim.eq.3) call setupvp(upup(1,3,i),
     $            uufld(1,3,i),uafld(1,3,i),uafld(1,3,i))

               call setupvp(upvp(1,1,i),
     $            uvfld(1,1,i),uafld(1,1,i),uafld(1,2,i))
               if (ldim.eq.3) call setupvp(upvp(1,2,i),
     $            uvfld(1,2,i),uafld(1,2,i),uafld(1,3,i))
               if (ldim.eq.3) call setupvp(upvp(1,3,i),
     $            uvfld(1,3,i),uafld(1,3,i),uafld(1,1,i))

               call divm1(flucv(1,1,i),
     $            upup(1,1,i),upvp(1,1,i),upvp(1,ldim,i))
               call divm1(flucv(1,2,i),
     $            upvp(1,1,i),upup(1,2,i),upvp(1,2,i))
               call divm1(flucv(1,3,i),
     $            upvp(1,3,i),upvp(1,2,i),upup(1,3,i))
            enddo
         endif

         if (ifpod(2)) then
            fname1='tavg.list'
            ifreads(1)=.false.
            ifreads(2)=.false.
            ifreads(3)=.true.
            call read_fields(fldtmp,fldtmp,tafld,
     $         nbavg,lbavg,0,ifreads,tk,fname1,iftmp)

            ifreads(1)=.true.
            ifreads(3)=.false.
            fname1='utms.list'
            call read_fields(utfld,fldtmp,fldtmp,
     $         nbavg,lbavg,0,ifreads,tk,fname1,iftmp)

            do i=1,nbavg
               call setupvp(uptp(1,1,i),
     $            utfld(1,1,i),uafld(1,1,i),tafld(1,i))
               call setupvp(uptp(1,2,i),
     $            utfld(1,2,i),uafld(1,2,i),tafld(1,i))
               if (ldim.eq.3) call setupvp(uptp(1,3,i),
     $            utfld(1,3,i),uafld(1,3,i),tafld(1,i))

               call divm1(fluct(1,i),
     $            uptp(1,1,i),uptp(1,2,i),uptp(1,ldim,i))
            enddo
         endif

         do i=1,nbavg
            call setdiff(udfld(1,1,i),tdfld(1,i),uafld(1,1,i),tafld(1,i)
     $         ,upup(1,1,i),upvp(1,1,i),uptp(1,1,i))
         enddo
      endif

      ifield=jfield

      if (nio.eq.0) write (6,*) 'exiting rom_init_fields'

      return
      end
c-----------------------------------------------------------------------
      subroutine setc_mhd(cuul, cbbl, cuubl, cbbul, fn1, fn2, fn3, fn4)

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      integer lxy
      integer nd

      real udu(lx1*ly1*lz1*lelt,ldim)
      real bdb(lx1*ly1*lz1*lelt,ldim)
      real bdu(lx1*ly1*lz1*lelt,ldim)
      real udb(lx1*ly1*lz1*lelt,ldim)
      
      common /convectb/ bc1v(ltd),bc2v(ltd),bc3v(ltd)
     $                 ,bu1v(ltd),bu2v(ltd),bu3v(ltd)
      common /convect/ c1v(ltd),c2v(ltd),c3v(ltd)
     $                 ,u1v(ltd),u2v(ltd),u3v(ltd)
      common /convect_all/ u1va(ltd,0:lb),u2va(ltd,0:lb),u3va(ltd,0:lb)
      common /convect_allb/ bu1va(ltd,0:lb),bu2va(ltd,0:lb)
     $                     ,bu3va(ltd,0:lb)
      

      character*128 fn1, fn2, fn3, fn4
      character*128 fnlint

      lxy = lx1*ly1*lz1*lelt
      nd=lxd*lyd*lzd*nelv

      !call lints(fnlint,fn1,128)
      !if (nid.eq.0) open (unit=100,file=fnlint)
      !call lints(fnlint,fn2, 128)
      !if (nid.eq.0) open (unit=101,file=fnlint)
      !call lints(fnlint,fn3,128)
      !if (nid.eq.0) open (unit=102,file=fnlint)
      !call lints(fnlint,fn4,128)
      !if (nid.eq.0) open (unit=103,file=fnlint)

      write(*,*) 'inside setc_mhd', ''

      call cpart(kc1,kc2,jc1,jc2,ic1,ic2,ncloc,nb,np,nid+1)

   !  setup bases for convection
      do i=0, nb
         call setcnv_u(ub(1,i),vb(1,i),wb(1,i))
         call copy(u1va(1,i),u1v, nd)
         call copy(u2va(1,i),u2v,nd)
         if (ldim.eq.3) call copy(u3va(1,i),u3v,nd)

         call setcnv_bu(bxb(1,i),byb(1,i),bzb(1,i))
         call copy(bu1va(1,i),bu1v,nd)
         call copy(bu2va(1,i),bu2v,nd)
         if (ldim.eq.3) call copy(bu3va(1,i),bu3v,nd)
      enddo


      !TODO: adjust for loop for lbb != lub
      
      do k=0,nb
         call setcnv_c(ub(1,k),vb(1,k),wb(1,k))
         call setcnv_bc(bxb(1,k),byb(1,k),bzb(1,k))
      do j=0,nb
      
      call copy(u1v,u1va(1,j),nd)
      call copy(u2v,u2va(1,j),nd)
      if (ldim.eq.3) call copy(u3v,u3va(1,j),nd)

      ! udu
      call convect_new(udu(1,1), u1v,.true.
     $                       ,c1v,c2v,c3v,.true.)
      call convect_new(udu(1,2), u2v,.true. 
     $                       ,c1v,c2v,c3v,.true.)
      call convect_new(udu(1,3), u3v,.true. 
     $                       ,c1v,c2v,c3v,.true.)

      call copy(bu1v,bu1va(1,j),ltd)
      call copy(bu2v,bu2va(1,j),ltd)

      ! bdb
      if (ldim.eq.3) call copy(bu3v,bu3va(1,j),ltd)

      call convect_new(bdb(1,1), bu1v,.true.
     $                       ,bc1v,bc2v,bc3v,.true.)
      call convect_new(bdb(1,2), bu2v,.true. 
     $                       ,bc1v,bc2v,bc3v,.true.)
      call convect_new(bdb(1,3), bu3v,.true. 
     $                       ,bc1v,bc2v,bc3v,.true.)

      ! udb
           call convect_new(udb(1,1), bu1v,.true.
     $                       ,c1v,c2v,c3v,.true.)
      call convect_new(udb(1,2), bu2v,.true. 
     $                       ,c1v,c2v,c3v,.true.)
      call convect_new(udb(1,3), bu3v,.true. 
     $                       ,c1v,c2v,c3v,.true.)

      ! bdu
           call convect_new(bdu(1,1), u1v,.true.
     $                       ,bc1v,bc2v,bc3v,.true.)
      call convect_new(bdu(1,2), u2v,.true. 
     $                       ,bc1v,bc2v,bc3v,.true.)
      call convect_new(bdu(1,3), u3v,.true. 
     $                       ,bc1v,bc2v,bc3v,.true.)


      ! TODO: check dimensions of rtmp5, rtmp6

      call uip(rtmp1,udu,uvwb(1,1,1),nb,0,ndim,nbat,fldtmp)
      call uip(rtmp4,bdb,uvwb(1,1,1),nb,0,ndim,nbat,fldtmp)
      call uip(rtmp5,udb,bxyzb(1,1,1),nb,0,ndim,nbat,fldtmp)
      call uip(rtmp6,bdu,bxyzb(1,1,1),nb,0,ndim,nbat,fldtmp) 

      do i=1,nb
         call setc_local(cuul,rtmp1(i,1),
     $                    ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
         call setc_local(cbbl,rtmp4(i,1),
     $                    ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
         call setc_local(cuubl,rtmp5(i,1),
     $                    ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
         call setc_local(cbbul,rtmp6(i,1),
     $                    ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)

      enddo
      enddo
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine fix_skew_mhd

      include 'SIZE'
      include 'MOR'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      real tmp(nb, nb)

      write(*,*) 'inside fix_skew_mhd', ''

      do i=1, nb+1
         call transpose(tmp,nb,cul((i-1)*nb*(nb+1) + nb + 1),nb)
         call cmult(tmp,-1.0,nb*nb)
         call add2(cul((i-1)*nb*(nb+1) + nb + 1),tmp,nb*nb)
         call cmult(cul((i-1)*nb*(nb+1) + nb + 1),0.5,nb*nb)

         call transpose(tmp,nb,cubl((i-1)*nb*(nb+1) + nb + 1),nb)
         call cmult(tmp,-1.0,nb*nb)
         call add2(cubl((i-1)*nb*(nb+1) + nb + 1),tmp,nb*nb)
         call cmult(cubl((i-1)*nb*(nb+1) + nb + 1),0.5,nb*nb)

         call transpose(tmp,nb,cbl((i-1)*nb*(nb+1) + nb + 1),nb)
         call cmult(tmp,-1.0,nb*nb)
         call add2(cbul((i-1)*nb*(nb+1) + nb + 1),tmp,nb*nb)
         call cmult(cbul((i-1)*nb*(nb+1) + nb + 1),0.5,nb*nb)
         call transpose(cbl((i-1)*nb*(nb+1) + nb + 1),nb
     $                  ,cbul((i-1)*nb*(nb+1) + nb + 1),nb)
         call cmult(cbl((i-1)*nb*(nb+1) + nb + 1),-1.0,nb*nb)
      enddo

c      call dump_serial(cul,ncloc,'ops/cuu ',nid)
c      call dump_serial(cbl,ncloc,'ops/cbb ',nid)
c      call dump_serial(cbul,ncloc,'ops/cbu ',nid)
c      call dump_serial(cubl,ncloc,'ops/cub ',nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine setc(cl,fname)

      ! set local convection tensor C

      ! cl    := local partition of C
      ! fname := name of C location for writing / reading

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)
      parameter (ltd=lxd*lyd*lzd*lelt)

      real cux(lt),cuy(lt),cuz(lt)

      common /scrcwk/ wk(lcloc),wk2(0:lub),cu(lt,ldim),wku(lt,ldim)
      common /convect_all/ u1va(ltd,0:lb),u2va(ltd,0:lb),u3va(ltd,0:lb)
      common /convect/ c1v(ltd),c2v(ltd),c3v(ltd),
     $                 u1v(ltd),u2v(ltd),u3v(ltd)

      real cl(lcloc)

      character*128 fname
      character*128 fnlint

      if (nio.eq.0) write (6,*) 'inside setc'

      call nekgsync
      conv_time=dnekclock()

      if (iffastc) call exitti('fastc not supported in setc_new$',nb)

      call cpart(kc1,kc2,jc1,jc2,ic1,ic2,ncloc,nb,np,nid+1) ! old indexing
c     call cpart(ic1,ic2,jc1,jc2,kc1,kc2,ncloc,nb,np,nid+1) ! new indexing

      n=lx1*ly1*lz1*nelv

      call lints(fnlint,fname,128)
      if (nid.eq.0) open (unit=100,file=fnlint)
      if (nio.eq.0) write (6,*) 'setc file:',fnlint

      nd=lxd*lyd*lzd*nelv

      if (rmode.eq.'ON '.or.rmode.eq.'ONB') then
         do k=0,nb
         do j=0,mb
         do i=1,mb
            cel=0.
            if (nid.eq.0) read(100,*) cel
            cel=glsum(cel,1)
            call setc_local(cl,cel,ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
         enddo
         enddo
         enddo
      else
         do i=0,nb
            if (ifield.eq.1) then
               call setcnv_u(ub(1,i),vb(1,i),wb(1,i))
               call copy(u1va(1,i),u1v,nd)
               call copy(u2va(1,i),u2v,nd)
               if (ldim.eq.3) call copy(u3va(1,i),u3v,nd)
            else
               call setcnv_u1(tb(1,i,1))
               call copy(u1va(1,i),u1v,nd)
            endif
         enddo
         do k=0,nb
            if (nio.eq.0) write (6,*) 'setc: ',k,'/',nb
            call setcnv_c(ub(1,k),vb(1,k),wb(1,k))
            do j=0,nb
               if (ifaxis) then
                  if (ifield.eq.1) then
                     call opcopy(wku(1,1),wku(1,2),wku(1,ldim),
     $                  ub(1,j),vb(1,j),wb(1,j))
                     call convect_axis(cu,ldim,
     $                  ub(1,k),vb(1,k),wb(1,k),wku)
                  else
                     call convect_axis(cu,1,
     $                  ub(1,k),vb(1,k),wb(1,k),tb(1,j,1))
                  endif
               else
                  if (ifield.eq.1) then
c                    call setcnv_u(ub(1,j),vb(1,j),wb(1,j))
                     call copy(u1v,u1va(1,j),nd)
                     call copy(u2v,u2va(1,j),nd)
                     if (ldim.eq.3) call copy(u3v,u3va(1,j),nd)
                     call cc(cu,ldim)
                  else
c                    call setcnv_u1(tb(1,j))
                     call copy(u1v,u1va(1,j),nd)
                     call cc(cu,1)
                  endif
               endif
               if (ifield.eq.1) then
                  call uip(rtmp1,cu,uvwb(1,1,1),nb,0,ndim,nbat,fldtmp)
               else
                  call uip(rtmp1,cu,tb(1,1,1),nb,0,1,nbat,fldtmp)
               endif
               do i=1,nb
                  call setc_local(cl,rtmp1(i,1),
     $               ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
                  if (nid.eq.0) write (100,*) rtmp1(i,1)
               enddo
            enddo
         enddo
      endif

      if (nid.eq.0) close (unit=100)

      call nekgsync
      if (nio.eq.0) write (6,*) 'conv_time: ',dnekclock()-conv_time
      if (nio.eq.0) write (6,*) 'ncloc=',ncloc

      if (nio.eq.0) write (6,*) 'exiting setc'

      return
      end
c-----------------------------------------------------------------------
      subroutine seta(a,a0,fname)

      ! set diffusion operator A

      ! a     := rom operator A w/o 0th mode interactions
      ! a0    := rom operator A w/  0th mode interactions
      ! fname := read target

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrseta/ wk1(lt)

      real a0(0:nb,0:nb),a(nb,nb)

      character*128 fname

      n=lx1*ly1*lz1*nelt

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading a...'
         call read_mat_serial(a0,nb+1,nb+1,fname,mb+1,nb+1,wk1,nid)
      else
         if (nio.eq.0) write (6,*) 'forming a...'
         if (ifield.eq.2.and.(nelgt.ne.nelvt)) then
            call copy(vdm1,vdiff(1,1,1,1,2),n)
            sc=1./param(8)
            call cmult(vdm1,sc,n)
         endif
         do j=0,nb
            if (nio.eq.0) write (6,*) 'seta: ',j,'/',nb
            nio=-1
            do i=0,nb
               if (ifield.eq.1) then
                  a0(i,j)=h10vip(ub(1,i),vb(1,i),wb(1,i),
     $                           ub(1,j),vb(1,j),wb(1,j))
               else
                  if (nelgt.ne.nelvt) then
                     a0(i,j)=h10sip_vd(tb(1,i,1),tb(1,j,1),vdm1)
                  else
                     a0(i,j)=h10sip(tb(1,i,1),tb(1,j,1))
                  endif
               endif
            enddo
            nio=nid
         enddo
      endif

      do j=1,nb
      do i=1,nb
         a(i,j)=a0(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine seta_mhd(a,a0,fname)

      ! set diffusion operator A

      ! a     := rom operator A w/o 0th mode interactions
      ! a0    := rom operator A w/  0th mode interactions
      ! fname := read target

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrseta/ wk1(lt)

      real a0(0:nb,0:nb),a(nb,nb)

      character*128 fname

      n=lx1*ly1*lz1*nelt

      if (nio.eq.0) write (6,*) 'inside seta_mhd'

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading a...'
         call read_mat_serial(a0,nb+1,nb+1,fname,mb+1,nb+1,wk1,nid)
      else
         if (nio.eq.0) write (6,*) 'forming a...'
         if (ifield.eq.2.and.(nelgt.ne.nelvt)) then
            call copy(vdm1,vdiff(1,1,1,1,2),n)
            sc=1./param(8)
            call cmult(vdm1,sc,n)
         endif
         do j=0,nb
            if (nio.eq.0) write (6,*) 'seta: ',j,'/',nb
            nio=-1
            do i=0,nb
               if (ifield.eq.1) then
                  a0(i,j)=h10vip(bxb(1,i),byb(1,i),bzb(1,i),
     $                           bxb(1,j),byb(1,j),bzb(1,j))
               else
                  if (nelgt.ne.nelvt) then
                     a0(i,j)=h10sip_vd(tb(1,i,1),tb(1,j,1),vdm1)
                  else
                     a0(i,j)=h10sip(tb(1,i,1),tb(1,j,1))
                  endif
               endif
            enddo
            nio=nid
         enddo
      endif

      do j=1,nb
      do i=1,nb
         a(i,j)=a0(i,j)
      enddo
      enddo

      call dump_serial(ab0,(nb+1)**2,'ops/ab ',nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine sets(s0,tt,fname)

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'GEOM'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lub+1)**2)
      common /scrsets/ tx(lx1,ly1,lx1,lelt),
     $                 ty(lx1,ly1,lx1,lelt),
     $                 tz(lx1,ly1,lx1,lelt)

      real s0(0:nb),tt(lx1,ly1,lz1,lelt,0:nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside sets'

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading s...'
         call read_mat_serial(s0,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
      else
         if (nio.eq.0) write (6,*) 'forming s...'
         nio=-1
         do i=0,nb
            if (nio.eq.0) write (6,*) 'sets: ',i,'/',nb
            s=0.
            do ie=1,nelt
            do ifc=1,2*ldim
               if (cbc(ifc,ie,2).ne.'E  '.and.
     $             cbc(ifc,ie,2).ne.'P  ') then
                  call facind(kx1,kx2,ky1,ky2,kz1,kz2,
     $                        lx1,ly1,lz1,ifc)
                  l=0
                  do iz=kz1,kz2
                  do iy=ky1,ky2
                  do ix=kx1,kx2
                     l=l+1
                     s=s+area(l,1,ifc,ie)*tt(ix,iy,iz,ie,i)
     $                                   *gn(ix,iy,iz,ie)
                  enddo
                  enddo
                  enddo
               endif
            enddo
            enddo
            s0(i)=glsum(s,1)
         enddo
         nio=nid
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setb(bop,b0,fname)

      ! set mass operator B

      ! bop     := rom operator B w/o 0th mode interactions
      ! b0    := rom operator B w/  0th mode interactions
      ! fname := read target

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'
      include 'SOLN'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lub+1)**2)

      real bop(nb,nb),b0(0:nb,0:nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside setb'

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading b...'
         call read_mat_serial(b0,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
      else
         if (ifield.eq.2.and.(nelgt.ne.nelvt)) then
            call copy(brhom1,vtrans(1,1,1,1,2),n)
c           sc=1./param(8)
c           call cmult(vdm1,sc,n)
         endif
         if (nio.eq.0) write (6,*) 'forming b...'
         do j=0,nb
            if (nio.eq.0) write (6,*) 'setb: ',j,'/',nb
            nio=-1
            do i=0,nb
               if (ifield.eq.1) then
                  b0(i,j)=wl2vip(ub(1,i),vb(1,i),wb(1,i),
     $                           ub(1,j),vb(1,j),wb(1,j))
               else
                  if (nelgt.eq.nelvt) then
                     b0(i,j)=wl2sip(tb(1,i,1),tb(1,j,1))
                  else
                     b0(i,j)=wl2sip_vd(tb(1,i,1),tb(1,j,1),brhom1)
                  endif
               endif
            enddo
            nio=nid
         enddo
      endif

      do j=1,nb
      do i=1,nb
         bop(i,j)=b0(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setb_mhd(bop,b0,fname)

      ! set mass operator B

      ! bop     := rom operator B w/o 0th mode interactions
      ! b0    := rom operator B w/  0th mode interactions
      ! fname := read target

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'
      include 'SOLN'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lub+1)**2)

      real bop(nb,nb),b0(0:nb,0:nb)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside setb_mhd'

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading b...'
         call read_mat_serial(b0,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
      else
         if (ifield.eq.2.and.(nelgt.ne.nelvt)) then
            call copy(brhom1,vtrans(1,1,1,1,2),n)
c           sc=1./param(8)
c           call cmult(vdm1,sc,n)
         endif
         if (nio.eq.0) write (6,*) 'forming b...'
         do j=0,nb
            if (nio.eq.0) write (6,*) 'setb: ',j,'/',nb
            nio=-1
            do i=0,nb
               if (ifield.eq.1) then      
                  b0(i,j)=wl2vip(bxb(1,i),byb(1,i),bzb(1,i),
     $                           bxb(1,j),byb(1,j),bzb(1,j))
               else        ! shouldn't run
                  if (nelgt.eq.nelvt) then   
                     b0(i,j)=wl2sip(tb(1,i,1),tb(1,j,1))
                  else
                     b0(i,j)=wl2sip_vd(tb(1,i,1),tb(1,j,1),brhom1)
                  endif
               endif
            enddo
            nio=nid
         enddo
      endif

      do j=1,nb
      do i=1,nb
         bop(i,j)=b0(i,j)
      enddo
      enddo

      call dump_serial(bb0,(nb+1)**2,'ops/bb ',nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine setae(a,fname)

      ! set eddy-viscosity/diffusivity operator

      ! a     := eddy-diffusion matrix
      ! fname := read target

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrseta/ wk1(lt)

      real a(0:nb,0:nb,nb)

      character*128 fname

      n=lx1*ly1*lz1*nelt

      if (rmode.eq.'ON '.or.rmode.eq.'ONB'.or.rmode.eq.'CP ') then
         if (nio.eq.0) write (6,*) 'reading ae...'
         call read_mat_serial(aue,nb+1,nb+1,fname,mb+1,nb+1,wk1,nid)
      else
         if (nio.eq.0) write (6,*) 'forming ae...'
         do k=1,nb
         do j=0,nb
            if (nio.eq.0) write (6,*) 'setae: ',k,j,'/',nb
            nio=-1
            do i=0,nb
               if (ifield.eq.1) then
                  a(i,j,k)=h10vip_vd(ub(1,i),vb(1,i),wb(1,i),
     $                             ub(1,j),vb(1,j),wb(1,j),udfld(1,1,k))
               else
                  a(i,j,k)=h10sip_vd(tb(1,i,1),tb(1,j,1),tdfld(1,k))
               endif
            enddo
            nio=nid
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setu

      ! set initial condition for ROM coefficients

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetu/ uu(lt),vv(lt),ww(lt),tt(lt),wk(lb+1)

      logical iftmp,ifexist

      if (nio.eq.0) write (6,*) 'inside setu'

      n=lx1*ly1*lz1*nelv

      if (rmode.eq.'ON '.or.rmode.eq.'CP ') then
         inquire (file='ops/u0',exist=ifexist)
         if (ifexist) call read_serial(u,nb+1,'ops/u0 ',wk,nid)

         inquire (file='ops/t0',exist=ifexist)
         if (ifexist) call read_serial(ut,nb+1,'ops/t0 ',wk,nid)
      else
         jfield=ifield

         if (ifrom(1)) then
            ifield=1
            call opsub2(uic,vic,wic,ub,vb,wb)
            if (ips.eq.'H10') then
               call h10pv2b(u,uic,vic,wic,ub,vb,wb)
            else if (ips.eq.'HLM') then
               call hlmpv2b(u,uic,vic,wic,ub,vb,wb)
            else
               call pv2b(u,uic,vic,wic,ub,vb,wb)
            endif
            do i=0,nb
               if (nio.eq.0) write (6,*) 'u',u(i)
            enddo
            call opadd2(uic,vic,wic,ub,vb,wb)
         else
            call rzero(u,(nb+1)*3)
            u((nb+1)*0)=1.
            u((nb+1)*1)=1.
            u((nb+1)*2)=1.
         endif

         if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $      call dump_serial(u,nb+1,'ops/u0 ',nid)

         if (ifrom(2)) then
            ifield=2
            call sub2(tic,tb(1,0,1),n)
            call ps2b(ut,tic,tb(1,0,1))
            do i=0,nb
               if (nio.eq.0) write (6,*) 'ut',ut(i)
            enddo
            call add2(tic,tb,n)
            if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $         call dump_serial(ut,nb+1,'ops/t0 ',nid)
         endif
         ifield=jfield

         if (ifcomb) then
            call opsub2(uic,vic,wic,ub,vb,wb)
            call sub2(tic,tb,n)
            call pc2b(u,ut,uic,vic,wic,tic,ub,vb,wb,tb)
            do i=0,nb
               if (nio.eq.0) write (6,*) 'u&ut',u(i)
            enddo
            call opadd2(uic,vic,wic,ub,vb,wb)
            call add2(tic,tb,n)
         endif

         call reconv(uu,vv,ww,u)
         call recont(tt,ut)

         iftmp=ifxyo
         ifxyo=.true.
         if (ifrecon) call outpost(uu,vv,ww,pr,tt,'rom')

         ttime=time
         jstep=istep
         time=1.
         istep=1
         call outpost(uic,vic,wic,pr,tic,'uic')
         ifxyo=.false.
         time=2.
         istep=2
         call outpost(uu,vv,ww,pr,tt,'uic')
         call opsub2(uu,vv,ww,uic,vic,wic)
         call sub2(tt,tic,n)
         time=3.
         istep=3
         call outpost(uu,vv,ww,pr,tt,'uic')
         time=ttime
         istep=jstep

         ifxyo=iftmp
      endif

      if (nio.eq.0) write (6,*) 'exiting setu'

      return
      end
c-----------------------------------------------------------------------
      subroutine setb_ic

      ! set initial condition for B field ROM coefficients

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetb/ bxx(lt),byy(lt),bzz(lt),bwk(lb+1)

      logical iftmp,ifexist

      if (nio.eq.0) write (6,*) 'inside setb_ic'

      n=lx1*ly1*lz1*nelv

      if (rmode.eq.'ON '.or.rmode.eq.'CP ') then
         inquire (file='ops/b0',exist=ifexist)
         if (ifexist) call read_serial(b,nb+1,'ops/b0 ',bwk,nid)

      else
         jfield=ifield

         if (ifrom(ifldb)) then
            ifield=1
            call opsub2(bxic,byic,bzic,bxb,byb,bzb)
            if (ips.eq.'H10') then
               call h10pv2b(b,bxic,byic,bzic,bxb,byb,bzb)
            else if (ips.eq.'HLM') then
               call hlmpv2b(b,bxic,byic,bzic,bxb,byb,bzb)
            else
               call pv2b(b,bxic,byic,bzic,bxb,byb,bzb)
            endif
            do i=0,nb
               if (nio.eq.0) write (6,*) 'b',b(i)
            enddo
            if (ifrom(ifldb)) then
               do i=0,nb
                  if (nio.eq.0) write (6,*) 'b',b(i)
               enddo
            endif
            call opadd2(bxic,byic,bzic,bxb,byb,bzb)
         else
            call rzero(b,(nb+1)*3)
            b((nb+1)*0)=1.
            b((nb+1)*1)=1.
            b((nb+1)*2)=1.
         endif

         if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ')
     $      call dump_serial(b,nb+1,'ops/b0 ',nid)

        

         ifield=jfield

         if (ifcomb) then
         !TODO: if combined
            call opsub2(uic,vic,wic,ub,vb,wb)
            call sub2(tic,tb,n)
            call pc2b(u,ut,uic,vic,wic,tic,ub,vb,wb,tb)
            do i=0,nb
               if (nio.eq.0) write (6,*) 'u&ut',u(i)
            enddo
            call opadd2(uic,vic,wic,ub,vb,wb)
            call add2(tic,tb,n)
         endif

         call reconb(bxx,byy,bzz,b)

         iftmp=ifxyo
         ifxyo=.true.
         if (ifrecon) call outpost(bxx,byy,bzz,pr,tt,'romb')

         ttime=time
         jstep=istep
         time=1.
         istep=1
         call outpost(bxic,byic,bzic,pr,tic,'bic')
         ifxyo=.false.
         time=2.
         istep=2
         call outpost(bxx,byy,bzz,pr,tt,'bic')
         call opsub2(bxx,byy,bzz,bxic,byic,bzic)
         time=3.
         istep=3
         call outpost(bxx,byy,bzz,pr,tt,'bic')
         time=ttime
         istep=jstep

         ifxyo=iftmp
      endif

      if (nio.eq.0) write (6,*) 'exiting setb_ic'

      return
      end
c-----------------------------------------------------------------------
      subroutine setf

      ! set forcing terms

      include 'SIZE'
      include 'SOLN'
      include 'MOR'
      include 'MASS'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetf/ wk1(lt),wk2(lt),wk3(lt),wk(lb+1)

      logical ifexist

      if (nio.eq.0) write (6,*) 'inside setf'

      call rzero(rf,nb)
      call rzero(rq,nb)
      call rzero(rg,nb)

      n=lx1*ly1*lz1*nelv

      if (ifforce.and.ifrom(1)) then ! assume fx,fy,fz has mass
         if (rmode.eq.'ON '.or.rmode.eq.'CP ') then
            inquire (file='ops/rf',exist=ifexist)
            if (ifexist) call read_serial(rf,nb,'ops/rf ',wk,nid)
         else
            do i=1,nb
               rf(i)=glsc2(ub(1,i),fx,n)+glsc2(vb(1,i),fy,n)
               if (ldim.eq.3) rf(i)=rf(i)+glsc2(wb(1,i),fz,n)
               if (nio.eq.0) write (6,*) rf(i),i,'rf'
            enddo
            call opcopy(wk1,wk2,wk3,fx,fy,fz)
            call opbinv1(wk1,wk2,wk3,wk1,wk2,wk3,1.)
            call outpost(wk1,wk2,wk3,pavg,tavg,'fff')
         endif
      endif

      if (ifsource.and.ifrom(2)) then ! assume qq has mass
         if (rmode.eq.'ON '.or.rmode.eq.'CP ') then
            inquire (file='ops/rq',exist=ifexist)
            if (ifexist) call read_serial(rq,nb,'ops/rq ',wk,nid)
         else
            do i=1,nb
               rq(i)=glsc2(qq,tb(1,i,1),n)
               if (nio.eq.0) write (6,*) rq(i),i,'rq'
            enddo
            call copy(wk1,qq,n)
            call binv1(wk1)
            call outpost(vx,vy,vz,pavg,wk1,'qqq')
            if (ifsrct) then
               do i=1,nb
                  rqt(i)=glsc2(qqxyz,tb(1,i,1),n)
                  if (nio.eq.0) write (6,*) rqt(i),i,'rqt'
               enddo
               call copy(wk1,qqxyz,n)
               call binv1(wk1)
               call outpost(vx,vy,vz,pavg,wk1,'qqq')
            endif
         endif
      endif

      if (nio.eq.0) write (6,*) 'exiting setf'

      return
      end
c-----------------------------------------------------------------------
      subroutine setupvp(upvp,uvfld,ufld,vfld)

      ! set fluctuation correlation field

      ! ufld  := <u>, temporal mean of first field
      ! vfld  := <v>, temporal mean of second field
      ! uvfld := <uv>, temporal mean of product of both fields
      ! upvp  := <(u-<u>)(v-<v>)>, temporal mean of product of fluctuations

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real upvp(lt),uvfld(lt),ufld(lt),vfld(lt)

      n=lx1*ly1*lz1*nelv

      call copy(upvp,uvfld,n)
      call admcol3(upvp,ufld,vfld,-1.,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine setfluc(fv,ftt,fname)

      ! set fluctuation ROM operators

      ! fv    := inner product of basis and fluctuation field
      ! ftt   := inner product of basis and fluctuation field
      ! fname := name

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lub+1)**2)

      real fv(nbavg,nbavg),ftt(nbavg,nbavg)

      character*128 fname

      if (nio.eq.0) write (6,*) 'inside setfluc',nbavg

      do j=1,nbavg
      do i=1,nbavg
         fv(i,j)=wl2vip(ub(1,i),vb(1,i),wb(1,i),
     $                  flucv(1,1,j),flucv(1,2,j),flucv(1,ldim,j))
         ftt(i,j)=wl2sip(tb(1,i,1),fluct(1,j))
      enddo
      enddo

      call dump_serial(fv,nbavg**2,'ops/fv ',nid)
      call dump_serial(ftt,nbavg**2,'ops/ft ',nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine final

      ! dump out coefficients, time, and misc. stats

      include 'SIZE'
      include 'MOR'

      common /romup/ rom_time,t1(0:lb) ,t2(0:lb)

      if (nio.eq.0) write (6,*) 'final...'
      if (nid.eq.0) then
         write (6,*) 'evalc_time:  ',evalc_time
         write (6,*) 'lu_time:     ',lu_time
         write (6,*) 'solve_time:  ',solve_time
         write (6,*) 'ustep_time:  ',ustep_time
         write (6,*) 'copt_time:   ',copt_time
         write (6,*) 'quasi_time:  ',quasi_time
         write (6,*) 'lnsrch_time: ',lnsrch_time
         write (6,*) 'compf_time:  ',compf_time
         write (6,*) 'compgf_time: ',compgf_time
         write (6,*) 'invhm_time:  ',invhm_time
         write (6,*) 'ucopt_active:',ucopt_count,
     $         '/',ad_step-1
         if (ifrom(2)) then
            write (6,*) 'tsolve_time: ',tsolve_time
            write (6,*) 'tstep_time:  ',tstep_time
            write (6,*) 'tcopt_active:',tcopt_count,
     $         '/',ad_step-1
         endif
         write (6,*) 'rom_time:    ',rom_time
         write (6,*) 'postu_time:  ',postu_time
         write (6,*) 'postt_time:  ',postt_time
         write (6,*) 'cp_time:     ',cp_time
      endif

      call dump_serial(u,nb+1,'ops/uf ',nid)
      if (ifrom(2)) then
         call dump_serial(ut,nb+1,'ops/tf ',nid)
      endif
      do i=0,nb
         t1(i)=u2a(1+i+i*(nb+1))-ua(i)*ua(i)
      enddo

      call dump_serial(t1,nb+1,'ops/uv ',nid)

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         call dump_serial(uas,nb+1,'ops/uas ',nid)
         call dump_serial(uvs,nb+1,'ops/uvs ',nid)
      endif

      call dump_serial(ua,nb+1,'ops/ua ',nid)
      call dump_serial(u2a,(nb+1)**2,'ops/u2a ',nid)
      if (ifrom(2)) then
         call dump_serial(uta,nb+1,'ops/uta ',nid)
         call dump_serial(ut2a,(nb+1)**2,'ops/ut2a ',nid)
      endif

      if (ifrecon) call dump_sfld

      return
      end
c-----------------------------------------------------------------------
      subroutine setbut(bx,by,bz)

      ! set buoyancy operators

      ! <bx,by,bz> := buoyancy operator in x,y,z

      include 'SIZE'
      include 'MOR'
      include 'AVG'
      include 'INPUT'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((lb+1)**2)
      common /scrk/ wk1(lt),wk2(lt),wk3(lt)

      logical iftmp

      real bx(0:nb,0:nb),by(0:nb,0:nb),bz(0:nb,0:nb)

      character*128 fname

      n=lx1*ly1*lz1*nelv

      if (rmode.eq.'ALL'.or.rmode.eq.'OFF'.or.rmode.eq.'AEQ') then
         do j=0,nb
         do i=0,nb
            bx(i,j)=glsc3(ub(1,i),tb(1,j,1),bm1,n)
            by(i,j)=glsc3(vb(1,i),tb(1,j,1),bm1,n)
            if (ldim.eq.3) bz(i,j)=glsc3(wb(1,i),tb(1,j,1),bm1,n)
         enddo
         enddo
         call dump_serial(bx,(nb+1)**2,'ops/buxt ',nid)
         call dump_serial(by,(nb+1)**2,'ops/buyt ',nid)
         if (ldim.eq.3) call dump_serial(bz,(nb+1)**2,'ops/buzt ',nid)
      else
         if (nio.eq.0) write (6,*) 'reading but...'
         fname='ops/buxt '
         call read_mat_serial(bx,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
         fname='ops/buyt '
         call read_mat_serial(by,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
         if (ldim.eq.3) then
            fname='ops/buzt '
            call read_mat_serial(bz,nb+1,nb+1,fname,mb+1,nb+1,tab,nid)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_tmean(ttmean,tt)

      ! set mean weighted scalar field

      ! ttmean := mean value
      ! tt     := target scalar field

      include 'SIZE'
      include 'TOTAL'

      real tt(lx1,ly1,lz1,lelt)

      nt=lx1*ly1*lz1*nelt

      ttmean=glsc3(tt,bm1,vtrans(1,1,1,1,2),nt)/
     $       glsc2(bm1,vtrans(1,1,1,1,2),nt)

      return
      end
c-----------------------------------------------------------------------
      subroutine set0mean(tt)

      ! set remove weighted-mean from scalar field

      ! tt := target scalar field

      include 'SIZE'
      include 'TOTAL'

      real tt(lx1,ly1,lz1,lelt)

      nt=lx1*ly1*lz1*nelt
      call calc_tmean(ttmean,tt)
      call cadd(tt,-ttmean,nt)

      if (nid.eq.0) write (6,*) ttmean,'tmean'

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_dudn(dudn,tt)

      ! calculate mean boundary gradient of scalar field

      ! dudn := mean gradient value
      ! tt   := target scalar field

      include 'SIZE'
      include 'TOTAL'

      common /mydudn/ tx(lx1,ly1,lz1,lelt),
     $                ty(lx1,ly1,lz1,lelt),
     $                tz(lx1,ly1,lz1,lelt)

      real tt(lx1,ly1,lz1,lelt)

      call gradm1(tx,ty,tz,tt)

      nt=lx1*ly1*lz1*nelt

      s=0.
      a=0.

      do ie=1,nelt
      do ifc=1,2*ldim
         if (cbc(ifc,ie,2).eq.'f  ') then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,ifc)
            l=0
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               l=l+1
               a=a+area(l,1,ifc,ie)
               s=s+area(l,1,ifc,ie)*(tx(ix,iy,iz,ie)*unx(l,1,ifc,ie)+
     $                               ty(ix,iy,iz,ie)*uny(l,1,ifc,ie)+
     $                               tz(ix,iy,iz,ie)*unz(l,1,ifc,ie))
            enddo
            enddo
            enddo
         endif
      enddo
      enddo

      s=glsum(s,1)
      a=glsum(a,1)

      if (nid.eq.0) write (6,*) s,a,s/a,'dudn'
      dudn=s/a

      return
      end
c-----------------------------------------------------------------------
      subroutine set1dudn(tt)

      ! set unit-mean boundary gradient for a scalar field

      ! tt   := target scalar field

      include 'SIZE'
      include 'TOTAL'

      real tt(lx1,ly1,lz1,lelt)

      nt=lx1*ly1*lz1*nelt
      call calc_dudn(dudn,tt)
      sf=1./dudn
      call cmult(tt,sf,nt)

      return
      end
c-----------------------------------------------------------------------
      subroutine seteddy(cl,fname)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt),cuy(lt),cuz(lt)

      common /scrcwk/ wk(lcloc),wk2(0:lub)

      real cl(lcloc)

      character*128 fname
      character*128 fnlint

      if (nio.eq.0) write (6,*) 'inside seteddy'

      call nekgsync
      conv_time=dnekclock()

      call cpart(kc1,kc2,jc1,jc2,ic1,ic2,ncloc,nb,np,nid+1) ! old indexing
c     call cpart(ic1,ic2,jc1,jc2,kc1,kc2,ncloc,nb,np,nid+1) ! new indexing

      n=lx1*ly1*lz1*nelv

      call lints(fnlint,fname,128)
      if (nid.eq.0) open (unit=100,file=fnlint)
      if (nio.eq.0) write (6,*) 'seteddy file:',fnlint

      if (rmode.eq.'ON '.or.rmode.eq.'ONB') then
         do k=0,nb
         do j=0,mb
         do i=1,mb
            cel=0.
            if (nid.eq.0) read(100,*) cel
            cel=glsum(cel,1)
            call setc_local(cl,cel,ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
         enddo
         enddo
         enddo
      else
         ifstrs=.true.
         iftmp=ifxyo
         ifxyo=.true.
         do k=0,nb
            if (nio.eq.0) write (6,*) 'seteddy: ',k,'/',nb
            do j=0,nb
               call ophx(cux,cuy,cuz,ub(1,j),vb(1,j),
     $                  wb(1,j),tb(1,k,4),zeros)
c              call outpost(cux,cuy,cuz,pr,tb(1,k,4),'stt')
               do i=1,nb
                  cel=op_glsc2_wt(
     $                  ub(1,i),vb(1,i),wb(1,i),cux,cuy,cuz,ones)
                  call setc_local(cl,cel,ic1,ic2,jc1,jc2,kc1,kc2,i,j,k)
                  if (nid.eq.0) write (100,*) cel
               enddo
            enddo
         enddo
      endif
         ifxyo=iftmp
         ifstrs=.false.

      if (nid.eq.0) close (unit=100)

      call nekgsync
      if (nio.eq.0) write (6,*) 'conv_time: ',dnekclock()-conv_time
      if (nio.eq.0) write (6,*) 'ncloc=',ncloc

      if (nio.eq.0) write (6,*) 'exiting seteddy'

      return
      end
c-----------------------------------------------------------------------
      subroutine copy2(a,b,nx,ny)
      real a(1,1),b(1,1)

      do i=1,nx
      do j=0,ny-1
         a(i,j)=b(i,j)
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copy3(a,b,nx,ny,nz)
      real a(1,1,1),b(1,1,1)

      do i=1,nx
      do j=1,ny
      do k=0,nz-1
         a(i,j,k)=b(i,j,k)
      enddo
      enddo
      enddo
      return
      end
