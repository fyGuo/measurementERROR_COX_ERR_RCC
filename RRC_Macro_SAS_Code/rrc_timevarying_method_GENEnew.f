cccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssc
c     This program is written by Xiaomei Liao.                         c
c     Last modified: Oct,16,2009                                       c
c                                                                      c
c     It applies on HPFS data the new method developed in the paper of c
c     "Survival analysis with error-prone time-varying covariates: a   c
c     risk set calibration approach" by Xiaomei Liao, David Zucker,    c
c     Yi Li, Donna Spiegelman.                                         c
c                                                                      c
c     The user can change the value of 'nex' and 'ncov' to define the  c
c     dimension of exposures and covariates. nex>0 and ncov>=0.        c 
c                                                                      c
ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssc
c                                                                      c 
c     This program is now applied on the HPFS lawrence paper data.     c
c     Last modified: March, 12, 2012                                   c                                                                     
c                                                                      c
c     It takes into account the correlation structures of the repeated c
c     measurements.                                                    c
c                                                                      c
c     It's for the general measurement error model.                    c
c     c=\alpha_0+\alpha_1*C                                            c
c                                                                      c
c     It only works for one exposure so far.                           c
c                                                                      c
ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssc
c
c     Defintion of the variables:
c     nmain: the number of subjects in the main study. (exact)
c     nval: the number of subjects in the validation study. (exact)
c     nageobs: the number of unique failure time in the main study.
c     nex: the number of exposures. (exact)
c     ncov: the number of covariates. (exact)
c     nage_val: the maximum number of observation for each subject in
c               the validation study.
c     nage_main: the maximum number of observation for each subject in
c               the main study.
c     mainobs: the total number of observation in the main study.(exact)
c     valobs: the total number of observation in the valdiation study.
c             (exact)
c
c     tf: the vector contains all the unique failure time in years.
c     groupnum: the number of subjects in each age group.
c     ngrp: the number of the age group.
c     age_low, age_up: the cutoff points of the age group.
c
c     exunit: the increment of the each exposure.
c 
c     icc : the correlation of the repeated true exposure.
c     tau : the damped exponential factor.
ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssc


      Module variable

      implicit none

      real*8 pi
      parameter (pi=3.14159265358979323846)

      integer idbgstop
      integer nmain, nval, nage_val, nage_main, nageobs, ncov,nex
c     parameter for residual analysis
      integer isres
c     parameter for the indicator of measurement error model
      integer ismodel
c
      integer, dimension(:), allocatable :: j_val, j_main

      integer mainobs, valobs
      integer outfile1
c     variables for the age group
      integer   ngrp, groupnum

      double precision, dimension(:), allocatable :: 
     +  age_low, age_up 

c     variable for the increment
      double precision, dimension(:), allocatable :: 
     +  exunit

c     variable for exposure and confounder names
      character(len=10), dimension(:), allocatable ::
     +  expname, covname

c     Covariance matrix parameters
      double precision icc, tau, slp

      integer, dimension(:), allocatable :: sub_id_val
      double precision, dimension(:,:,:), allocatable :: 
     +  dtcalc_val, ffqcalc_val 
      double precision, dimension(:,:), allocatable :: 
     +  age_val, bc_val
      double precision, dimension(:), allocatable ::
     +  t_val, del_val
      double precision, dimension(:,:,:), allocatable :: 
     +  zt_depval, zs_depval

cc    add in covariates
      double precision, dimension(:,:,:), allocatable :: 
     +  covar_val
      double precision, dimension(:,:), allocatable ::
     +  dtcalc_valtmp, ffqcalc_valtmp
      integer, dimension(:), allocatable ::
     +  sub_id_main
      double precision, dimension(:,:), allocatable ::
     +  age_main, bc_main
      double precision, dimension(:), allocatable ::
     +  t_main, del_main
      double precision, dimension(:,:,:), allocatable :: 
     +  zs_depmain, zs_depmain_po

cc    add in covariates
      double precision, dimension(:,:,:), allocatable :: 
     +  covar_main, zs_common, covar_common

      integer ndist,nfa,nal
      integer, dimension(:), allocatable ::
     +  ilow, iup, fnum
      integer, dimension(:,:), allocatable ::
     +  fsub
      integer nfa2, vaintl

      double precision, dimension(:), allocatable ::
     +  evcnt, tf, ageobs, tf_ageyr, tf_uyr
      integer, dimension(:), allocatable ::
     +  ry_val, ry_main, cntcom, count1, count2,
     +  count3, count4
      integer case_cnt
 
      double precision, dimension(:,:,:), allocatable ::
     +  zshat_depmain, covarhat_main, xxcom, yycom
      double precision, dimension(:,:), allocatable ::
     +  xcom, ycom
      double precision, dimension(:,:,:), allocatable ::
     +  acom
      double precision, dimension(:), allocatable ::
     +  r22

      integer, dimension(:,:), allocatable :: id_val
      double precision, dimension(:), allocatable ::
     +  betaes_naiv, serror_naiv, betaes_rrc, serror_rrc
      double precision, dimension(:), allocatable ::
     + meanx, meany, varx, vary, covxy

      integer maxi_case
      double precision, dimension(:,:), allocatable ::
     + rsquare

c These are for filenames

      character*255 pathname
      character*20 procid


      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine allocvar
      use variable
      implicit none

      character*255 filename
      filename=trim(pathname)//'/runparams.'//trim(procid)//'.dat'
      open(1013,file=filename,action='read')
      read(1013,*) nmain
      read(1013,*) nval
      read(1013,*) nageobs
      read(1013,*) nage_val
      read(1013,*) nage_main
      read(1013,*) nex
      read(1013,*) ncov
      read(1013,*) isres
      read(1013,*) ismodel
      read(1013,*) mainobs
      read(1013,*) valobs
      read(1013,*) groupnum
      close(1013)

      allocate(j_val(nval))
      allocate(j_main(nmain))
      allocate(age_low(nageobs))
      allocate(age_up(nageobs))
      allocate(exunit(nex))
      allocate(expname(nex))
      allocate(covname(ncov))
      allocate(sub_id_val(nval))
      allocate(dtcalc_val(nval,nage_val,nex))
      allocate(ffqcalc_val(nval,nage_val,nex))
      allocate(age_val(nval,nage_val))
      allocate(bc_val(nval,nage_val))
      allocate(t_val(nval))
      allocate(del_val(nval))
      allocate(zt_depval(nval,nage_val,nex))
      allocate(zs_depval(nval,nage_val,nex+ncov))
      allocate(covar_val(nval,nage_val,0:ncov))
      allocate(dtcalc_valtmp(valobs,nex))
      allocate(ffqcalc_valtmp(valobs,nex))
      allocate(sub_id_main(nmain))
      allocate(age_main(nmain,nage_main))
      allocate(bc_main(nmain,nage_main))
      allocate(t_main(nmain))
      allocate(del_main(nmain))
      allocate(zs_depmain(nmain,nage_main,nex))
      allocate(zs_depmain_po(nmain,nage_main,nex))
      allocate(covar_main(nmain,nage_main,0:ncov))
      allocate(zs_common(nmain,nage_main,nex))
      allocate(covar_common(nmain,nage_main,0:ncov))
      allocate(ilow(nmain))
      allocate(iup(nmain))
      allocate(fnum(nmain))
      allocate(fsub(nageobs,2))
      allocate(evcnt(nmain))
      allocate(tf(nageobs))
      allocate(ageobs(nmain))
      allocate(tf_ageyr(nageobs))
      allocate(tf_uyr(nageobs))
      allocate(ry_val(nageobs))
      allocate(ry_main(nageobs))
      allocate(cntcom(nageobs))
      allocate(count1(nageobs))
      allocate(count2(nageobs))
      allocate(count3(nageobs))
      allocate(count4(nageobs))
      allocate(zshat_depmain(nmain,nage_main,nex))
      allocate(covarhat_main(nmain,nage_main,0:ncov))
      allocate(xxcom(valobs,nageobs,nex+ncov))
      allocate(yycom(valobs,nageobs,nex))
      allocate(xcom(valobs,nageobs))
      allocate(ycom(valobs,nageobs))
      allocate(acom(nageobs,1+nex+ncov,nex))
      allocate(r22(nageobs))
      allocate(id_val(valobs,nageobs))
      allocate(betaes_naiv(nex+ncov))
      allocate(serror_naiv(nex+ncov))
      allocate(betaes_rrc(nex+ncov))
      allocate(serror_rrc(nex+ncov))
      allocate(meanx(nageobs))
      allocate(meany(nageobs))
      allocate(varx(nageobs))
      allocate(vary(nageobs))
      allocate(covxy(nageobs))
      maxi_case=nageobs
      allocate(rsquare(nageobs,nex))

      return
      end
ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine deallocvar
      use variable
      implicit none

      deallocate(j_val)
      deallocate(j_main)
      deallocate(age_low)
      deallocate(age_up)
      deallocate(exunit)
      deallocate(expname)
      deallocate(covname)
      deallocate(sub_id_val)
      deallocate(dtcalc_val)
      deallocate(ffqcalc_val)
      deallocate(age_val)
      deallocate(bc_val)
      deallocate(t_val)
      deallocate(del_val)
      deallocate(zt_depval)
      deallocate(zs_depval)
      deallocate(covar_val)
      deallocate(dtcalc_valtmp)
      deallocate(ffqcalc_valtmp)
      deallocate(sub_id_main)
      deallocate(age_main)
      deallocate(bc_main)
      deallocate(t_main)
      deallocate(del_main)
      deallocate(zs_depmain)
      deallocate(zs_depmain_po)
      deallocate(covar_main)
      deallocate(zs_common)
      deallocate(covar_common)
      deallocate(ilow)
      deallocate(iup)
      deallocate(fnum)
      deallocate(fsub)
      deallocate(evcnt)
      deallocate(tf)
      deallocate(ageobs)
      deallocate(tf_ageyr)
      deallocate(tf_uyr)
      deallocate(ry_val)
      deallocate(ry_main)
      deallocate(cntcom)
      deallocate(count1)
      deallocate(count2)
      deallocate(count3)
      deallocate(count4)
      deallocate(zshat_depmain)
      deallocate(covarhat_main)
      deallocate(xxcom)
      deallocate(yycom)
      deallocate(xcom)
      deallocate(ycom)
      deallocate(acom)
      deallocate(r22)
      deallocate(id_val)
      deallocate(betaes_naiv)
      deallocate(serror_naiv)
      deallocate(betaes_rrc)
      deallocate(serror_rrc)
      deallocate(meanx)
      deallocate(meany)
      deallocate(varx)
      deallocate(vary)
      deallocate(covxy)
      deallocate(rsquare)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      program hpfs_data

      use variable
      implicit none

      integer i,j,k,k1,b, jtemp

      integer, dimension(:), allocatable ::
     +  sub_id_valtmp
      double precision, dimension(:), allocatable ::
     +  age_valtmp, bc_valtmp
      integer, dimension(:), allocatable ::
     +  sub_id_maintmp
      double precision, dimension(:,:), allocatable ::
     +  ffqcalc_maintmp, ind_main

cc    add in covariates
      double precision, dimension(:,:), allocatable ::
     +  covar_maintmp, covar_valtmp
      double precision, dimension(:), allocatable ::
     +  age_maintmp, bc_maintmp
      double precision a, age_temp(43), age_start, age_end
      double precision, dimension(:), allocatable ::
     +  xtmp, xhattmp, xave, bxave, c_corr, x_corr
      integer, dimension(:), allocatable ::
     +  xcnt, xhatcnt, cnt2, cnmiss, xnmiss, casecnt_age
      double precision, dimension(:), allocatable ::
     +  a0temp, a1temp 

c  Get path and proc id for reconstructing input filenames
      character*255 arg
      character*255 filename

      integer narg
      integer iargc

      narg = iargc()
      if ( narg .ne. 2 ) then
        write(*,*) 'usage: testprog < path > < pid >'
        stop
      endif

      call getarg(1,arg)
      pathname = arg
      call getarg(2,arg)
      procid = arg

c     allocate memory for module variable
      call allocvar

      allocate(sub_id_valtmp(valobs))
      allocate(age_valtmp(valobs))
      allocate(bc_valtmp(valobs))
      allocate(sub_id_maintmp(mainobs))
      allocate(ffqcalc_maintmp(mainobs,nex))
      allocate(ind_main(nmain,nage_main))
      allocate(covar_maintmp(mainobs,0:ncov))
      allocate(covar_valtmp(valobs,0:ncov))
      allocate(age_maintmp(mainobs))
      allocate(bc_maintmp(mainobs))
      allocate(xtmp(nageobs))
      allocate(xhattmp(nageobs))
      allocate(xave(nageobs))
      allocate(bxave(nageobs))
      allocate(c_corr(nageobs))
      allocate(x_corr(nageobs))
      allocate(xcnt(nageobs))
      allocate(xhatcnt(nageobs))
      allocate(cnmiss(nageobs))
      allocate(xnmiss(nageobs))
      allocate(cnt2(nageobs))
      allocate(casecnt_age(nageobs))
      allocate(a0temp(nageobs))
      allocate(a1temp(nageobs))

c     outfile1 is the unit number for the output file
      outfile1=1011

      count2=-1
      count3=-1
      count4=-1



c     Read in the icc and tau for the covariance matrix
      filename=trim(pathname)//'/covpara.'//trim(procid)//'.dat'
      open(41,file=filename,action='read')
      read(41,*) icc, tau
      close(41)

c     Read in the units for each exposure.
      filename=trim(pathname)//'/exunit.'//trim(procid)//'.dat'
      open(51,file=filename,action='read')
      do i=1,nex
         read(51,*) exunit(i)
      enddo
      close(51)

c     Read in the variable name for each exposure.
      filename=trim(pathname)//'/expnames.'//trim(procid)//'.dat'
      open(61,file=filename,action='read')
      do i=1,nex
         read(61,*) expname(i)
      enddo
      close(61)

c     Read in the variable name for each confounder.
      filename=trim(pathname)//'/covnames.'//trim(procid)//'.dat'
      open(71,file=filename,action='read')
      do i=1,ncov
         read(71,*) covname(i)
      enddo
      close(71)


c     Read in the main study data set.

      filename=trim(pathname)//'/main.'//trim(procid)//'.dat'
      open(101,file=filename, action='read')
      do i=1,mainobs
         if (ncov.ge.1) then
            read(101,*) 
     *           sub_id_maintmp(i),age_maintmp(i),bc_maintmp(i),
     *           ffqcalc_maintmp(i,1:nex),covar_maintmp(i,1:ncov)  
         else
            read(101,*) 
     *           sub_id_maintmp(i),age_maintmp(i),bc_maintmp(i),
     *           ffqcalc_maintmp(i,1:nex)
 
         endif  
      enddo      

      close(101)

c     scaling by the user input increment
      do j=1,nex
         ffqcalc_maintmp(:,j)=ffqcalc_maintmp(:,j)/exunit(j)
      enddo
 

c      way to check the correctness of reading of the main study data
       print *,'i=',i
       print *,sub_id_maintmp(i-1),age_maintmp(i-1),bc_maintmp(i-1),
     *     ffqcalc_maintmp(i-1,1:nex),covar_maintmp(i-1,0:ncov)



c     transformation of main study
      k=1
      j=1
      sub_id_main(k)=sub_id_maintmp(k)
      zs_depmain(k,j,1:nex)=ffqcalc_maintmp(k,1:nex)
      age_main(k,j)=age_maintmp(k)
      bc_main(k,j)=bc_maintmp(k)
      if (ncov.ge.1) then
         covar_main(k,j,1:ncov)=covar_maintmp(k,1:ncov)
      endif
      j_main(k)=j

      do i=2,mainobs
         if(sub_id_maintmp(i).eq.sub_id_maintmp(i-1)) then
            j=j+1
            sub_id_main(k)=sub_id_maintmp(i)
            zs_depmain(k,j,1:nex)=ffqcalc_maintmp(i,1:nex)
            age_main(k,j)=age_maintmp(i)
            bc_main(k,j)=bc_maintmp(i)
            if (ncov.ge.1) then
               covar_main(k,j,1:ncov)=covar_maintmp(i,1:ncov)
            endif
            j_main(k)=j
         else         
            k=k+1
            j=1
            j_main(k)=j
            sub_id_main(k)=sub_id_maintmp(i)
            zs_depmain(k,j,1:nex)=ffqcalc_maintmp(i,1:nex)
            age_main(k,j)=age_maintmp(i)
            bc_main(k,j)=bc_maintmp(i)
            if (ncov.ge.1) then
               covar_main(k,j,1:ncov)=covar_maintmp(i,1:ncov)
            endif
         endif
      enddo
      
c      print *, '# of subjects in the main study: nmain=', k
c      print *, 'The last subject in the main study:', sub_id_main(k)
c      print *, 'The last subject data:', zs_depmain(k,1:j,1:nex), 
c     *age_main(k,1:j)


c     To determine follow-up time and censoring indicator 
      do i=1,nmain
         t_main(i)=age_main(i,j_main(i))
         del_main(i)=bc_main(i,j_main(i))
      enddo
       

      call sort_main(nmain,nage_main)
      
      print *, 'The following results are for the naive cox method: '
      
      zs_depmain_po=zs_depmain

      call cal_cum_ave(zs_depmain_po,zs_depmain)

      call naivecoxmethod

      zs_depmain=zs_depmain_po


c     Read in the validation study data set 

      filename=trim(pathname)//'/validation.'//trim(procid)//'.dat'
      open(201,file=filename, action='read')
      do i=1,valobs
         if (ncov.ge.1) then
            read(201,*) sub_id_valtmp(i),age_valtmp(i),
     *           (dtcalc_valtmp(i,k),ffqcalc_valtmp(i,k), k=1,nex),
     *           (covar_valtmp(i,j),j=1,ncov)
         else
            read(201,*) sub_id_valtmp(i),age_valtmp(i),
     *           (dtcalc_valtmp(i,k),ffqcalc_valtmp(i,k),k=1,nex)
         endif
      enddo

      close(201)

c     scaling by the user input increment
      do j=1,nex
         ffqcalc_valtmp(:,j)=ffqcalc_valtmp(:,j)/exunit(j)
         dtcalc_valtmp(:,j)=dtcalc_valtmp(:,j)/exunit(j)
      enddo


c     tranformation of validation study
      k=1
      j=1
      sub_id_val(k)=sub_id_valtmp(k)
      dtcalc_val(k,j,1:nex)=dtcalc_valtmp(k,1:nex)
      ffqcalc_val(k,j,1:nex)=ffqcalc_valtmp(k,1:nex)
      age_val(k,j)=age_valtmp(k)
      if (ncov.ge.1) then
         covar_val(k,j,1:ncov)=covar_valtmp(k,1:ncov)
      endif
      j_val(k)=j

      do i=2,valobs
         if (sub_id_valtmp(i).eq.sub_id_valtmp(i-1)) then
            j=j+1
            sub_id_val(k)=sub_id_valtmp(i)
            dtcalc_val(k,j,1:nex)=dtcalc_valtmp(i,1:nex)
            ffqcalc_val(k,j,1:nex)=ffqcalc_valtmp(i,1:nex)
            age_val(k,j)=age_valtmp(i)
            if (ncov.ge.1) then
               covar_val(k,j,1:ncov)=covar_valtmp(i,1:ncov)
            endif
            j_val(k)=j
         else
            k=k+1
            j=1
            j_val(k)=j
            sub_id_val(k)=sub_id_valtmp(i)
            dtcalc_val(k,j,1:nex)=dtcalc_valtmp(i,1:nex)
            ffqcalc_val(k,j,1:nex)=ffqcalc_valtmp(i,1:nex)
            age_val(k,j)=age_valtmp(i)
            if (ncov.ge.1) then
               covar_val(k,j,1:ncov)=covar_valtmp(i,1:ncov)
            endif
         endif
      enddo


c     To determine follow-up time and censoring indicator 
      do i=1,nval
         t_val(i)=age_val(i,j_val(i))
      enddo

      call sort_val(nval,nage_val)

      zt_depval(:,:,1:nex)=dtcalc_val
      zs_depval(:,:,1:nex)=ffqcalc_val

      if (ncov.ge.1) then
         zs_depval(:,:,(1+nex):(nex+ncov))=covar_val(:,:,1:ncov)
      endif

      print *, 'The following results are for the RRC method: '

      call rrc
      call rrccoxmethod

c     summarizing the result

      call data_summary

      deallocate(sub_id_valtmp)
      deallocate(age_valtmp)
      deallocate(bc_valtmp)
      deallocate(sub_id_maintmp)
      deallocate(ffqcalc_maintmp)
      deallocate(ind_main)
      deallocate(covar_maintmp)
      deallocate(covar_valtmp)
      deallocate(age_maintmp)
      deallocate(bc_maintmp)
      deallocate(xtmp)
      deallocate(xhattmp)
      deallocate(xave)
      deallocate(bxave)
      deallocate(c_corr)
      deallocate(x_corr)
      deallocate(xcnt)
      deallocate(xhatcnt)
      deallocate(cnmiss)
      deallocate(xnmiss)
      deallocate(cnt2)
      deallocate(casecnt_age)
      deallocate(a0temp)
      deallocate(a1temp)

      call deallocvar

c     corresponding to the program 
      end


ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c     this subroutine is to calculate the 12-month average exposure.

      subroutine cal_12mon_ave(zp,zs)

      use variable
      implicit none

      integer i,j,k,j1
      double precision zp(nmain,nage_main,nex),zs(nmain,nage_main,nex),
     +  s

      do k=1,nex
         do i=1,nmain
            do j=1,j_main(i)
               s=0.0d0
               if (j.lt.12) then
                  do j1=1,j
                     s=s+zp(i,j1,k)
                  enddo
                  zs(i,j,k)=s/dble(j)
               else
                  do j1=j-11,j
                     s=s+zp(i,j1,k)
                  enddo
                  zs(i,j,k)=s/12.0d0
               endif
            enddo
         enddo
      enddo
      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c     this subroutine is to calculate the cumulative average exposure.

      subroutine cal_cum_ave(zp,zs)

      use variable
      implicit none

      integer i,j,k,j1
      double precision zp(nmain,nage_main,nex),zs(nmain,nage_main,nex),
     +  s

      do k=1,nex
         do i=1,nmain
            do j=1,j_main(i)
               s=0.0d0
               do j1=1,j
                  s=s+zp(i,j1,k)
               enddo
               zs(i,j,k)=s/dble(j)
            enddo
         enddo
      enddo
      
      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c     this subroutine is to calculate the cumulative total exposure.

      subroutine cal_cum_tot(zp,zs)

      use variable
      implicit none

      integer i,j,k,j1
      double precision zp(nmain,nage_main,nex),zs(nmain,nage_main,nex),
     +  s

      do k=1,nex
         do i=1,nmain
            do j=1,j_main(i)
               s=0.0d0
               do j1=1,j
                  s=s+zp(i,j1,k)
               enddo
               zs(i,j,k)=s
            enddo
         enddo
      enddo
      
      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine sort_main(nmax,ntm)

      use variable
      implicit none

c     on the entry variable
      integer nmax, ntm

c     dummy variable
      integer wksp(nmax)

      integer, dimension(:), allocatable ::
     +  indx, tt_main
      double precision, dimension(:), allocatable ::
     +  wksp1, rsk
      double precision, dimension(:,:), allocatable ::
     +  wksp2
      double precision, dimension(:,:,:), allocatable ::
     +  wksp3, wksp4

      real*8 told,rsknum,rm
      integer ix,i,j
 
c     external dsvrgp
      external quicksort

      allocate(wksp1(nmax))
      allocate(wksp2(nmax,ntm))
      allocate(wksp3(nmax,ntm,0:ncov))
      allocate(wksp4(nmax,ntm,nex))
      allocate(rsk(nmax))
      allocate(indx(nmax))
      allocate(tt_main(nmax))
  
c     SORT DATA
      do i = 1, nmax
         indx(i) = i
      enddo

      do i = 1, nmax
        tt_main(i) = t_main(i)
      enddo
      call quicksort(nmax,t_main, indx)
      do i = 1, nmax
        t_main(i) = tt_main(indx(i))
      enddo
c     call dsvrgp(nmax,t_main,t_main,indx)

      do i = 1, nmax
         wksp(i) = j_main(i)
      enddo
      do i = 1, nmax
         j_main(i) = wksp(indx(i))
      enddo
      do i = 1, nmax 
         wksp(i) = sub_id_main(i) 
      enddo 
      do i = 1, nmax 
         sub_id_main(i) = wksp(indx(i)) 
      enddo 

      do i = 1, nmax
         wksp1(i) = del_main(i)
      enddo
      do i = 1, nmax
         del_main(i) = wksp1(indx(i))
      enddo
      do i = 1, nmax
         wksp4(i,:,:) = zs_depmain(i,:,:)
      enddo
      do i = 1, nmax
         zs_depmain(i,:,:) = wksp4(indx(i),:,:)
      enddo
      do i = 1, nmax
         wksp2(i,:) = age_main(i,:)
      enddo
      do i = 1, nmax
         age_main(i,:) = wksp2(indx(i),:)
      enddo
      do i = 1, nmax
         wksp2(i,:) = bc_main(i,:)
      enddo
      do i = 1, nmax
         bc_main(i,:) = wksp2(indx(i),:)
      enddo
      if (ncov.ge.1) then
         do i = 1, nmax
            wksp3(i,:,1:ncov)= covar_main(i,:,1:ncov)
         enddo
         do i = 1, nmax
            covar_main(i,:,1:ncov) = wksp3(indx(i),:,1:ncov)
         enddo
      endif


c     IDENTIFY UNIQUE FOLLOW-UP TIME
      ix = 1
      told = t_main(1)
      ilow(1) = 1
      rsknum = nmax
      rsk(1) = rsknum
      rm = 1.0d0
      evcnt(1) = del_main(1)
      ageobs(1)=told
      do i  = 2, nmax
         if (t_main(i) .ne. told) then
            iup(ix) = i-1
            rsk(ix) = rsknum
            rsknum = rsknum - rm
            ix = ix+1
            ilow(ix) = i
            rm = 0.0d0
            evcnt(ix) = 0.0d0
            told = t_main(i)
            ageobs(ix)=told
         endif
         rm = rm + 1.0d0
         evcnt(ix) = evcnt(ix) + del_main(i)
         if (i .eq. nmax)  iup(ix) = nmax
      enddo
      ndist = ix
c      print *, 'ndist=', ndist

c     IDENTIFY UNIQUE FAILURE TIME
      nfa=0
      tf=0.0d0
      fnum=0
      do i=1,ndist
         if (evcnt(i).ne.0.0d0) then
            nfa=nfa+1
            fsub(nfa,1)=ilow(i)
            fsub(nfa,2)=iup(i)
            j=ilow(i)
            tf(nfa)=t_main(j)
c            fnum(nfa)=nint(t_main(j)/dt1)
         endif
      enddo
c      nal=sum(fnum)

      print *, '# of unique failure time=', nfa

      deallocate(wksp1)
      deallocate(wksp2)
      deallocate(wksp3)
      deallocate(wksp4)
      deallocate(rsk)
      deallocate(indx)
      deallocate(tt_main)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine sort_val(nmax,ntm)
c     this subroutine is used to sort the data in the validation study

      use variable

      implicit none
  
c     variables on the entry
      integer nmax, ntm

c     dummy variables
      integer i

      integer, dimension(:), allocatable ::
     +  indx, wksp
      double precision, dimension(:), allocatable ::
     +  wksp1, tt_val
      double precision, dimension(:,:), allocatable ::
     +  wksp2
      double precision, dimension(:,:,:), allocatable ::
     +  wksp3, wksp4

      external quicksort

      allocate(indx(nmax))
      allocate(wksp(nmax))
      allocate(wksp1(nmax))
      allocate(wksp2(nmax,ntm))
      allocate(wksp3(nmax,ntm,0:ncov))
      allocate(wksp4(nmax,ntm,nex))
      allocate(tt_val(nmax))
  

c     SORT DATA
      do i = 1,nmax
         indx(i) = i
      enddo
      do i = 1, nmax
        tt_val(i) = t_val(i)
      enddo
      call quicksort(nmax,t_val, indx)
      do i = 1, nmax
        t_val(i) = tt_val(indx(i))
      enddo
c     call dsvrgp(nmax,t_val,t_val,indx)

      do i = 1, nmax
         wksp(i) = j_val(i)
      enddo
      do i = 1, nmax
         j_val(i) = wksp(indx(i))
      enddo
      do i = 1, nmax 
         wksp(i) = sub_id_val(i) 
      enddo 
      do i = 1, nmax 
         sub_id_val(i) = wksp(indx(i)) 
      enddo 

      do i = 1, nmax
         wksp4(i,:,:) = dtcalc_val(i,:,:)
      enddo
      do i = 1, nmax
         dtcalc_val(i,:,:) = wksp4(indx(i),:,:)
      enddo
      do i = 1, nmax
         wksp4(i,:,:) = ffqcalc_val(i,:,:)
      enddo
      do i = 1, nmax
         ffqcalc_val(i,:,:) = wksp4(indx(i),:,:)
      enddo
      
      do i = 1, nmax
         wksp2(i,:) = age_val(i,:)
      enddo
      do i = 1, nmax
         age_val(i,:) = wksp2(indx(i),:)
      enddo
    
      if (ncov.ge.1) then
         do i = 1, nmax 
            wksp3(i,:,1:ncov)= covar_val(i,:,1:ncov) 
         enddo 
         do i = 1, nmax 
            covar_val(i,:,1:ncov) = wksp3(indx(i),:,1:ncov) 
         enddo 
      endif
    
      deallocate(indx)
      deallocate(wksp)
      deallocate(wksp1)
      deallocate(wksp2)
      deallocate(wksp3)
      deallocate(wksp4)
      deallocate(tt_val)

      return
      end

ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine naivecoxmethod

      use variable
      implicit none

      integer nd, ifail

      double precision, dimension(:), allocatable ::
     +  betaes, st_error

      nd=nex+ncov
      allocate(betaes(nex+ncov))
      allocate(st_error(nex+ncov))

      zs_common=zs_depmain
      covar_common=covar_main

      if (nd.ge.2) then
         call cox_bcoah(nd,betaes,st_error,ifail)
      else if (nd.eq.1) then
         idbgstop = 0
         call cox_zbren(nd,betaes,st_error,ifail)
      else
         print *, 'The dimension of data is zero!'
      endif
         
      betaes_naiv=betaes
      serror_naiv=st_error

      if (ifail.eq.0) then
         print *,'naive','beta=', betaes
         print *,'se=', st_error
      else
         print *,'The method does not converge.'
      endif

      deallocate(betaes)
      deallocate(st_error)
      
      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
 
      subroutine cox_bcoah(nd,betaes,st_error,ifail)

      use variable
      implicit none

c     variables on the entry
      integer nd, ifail
      double precision betaes(nd), st_error(nd)
      
c     dummy variables
      real*8 fvalue,fscale 

      integer, dimension(:), allocatable ::
     +  ipiv
      double precision, dimension(:), allocatable ::
     +  xguess, xlb, xub, xscale, x, eval, work
      double precision, dimension(:,:), allocatable ::
     +  hes_fval, hesinv, r

      integer iparam(7), rparam(7),ibtype
      integer isev, iwk(nd),i, ier 
      integer info, ir, ic, nw

c ----------------------------------------------------------------------------------
      double precision wrkf(5000)
      double precision fval(1), gval(1)
      integer iwrkf(5000), iprint
      external obj,constr,gradobj,grcnfd
c ----------------------------------------------------------------------------------

      external rrccox_likelihood, rrcgradient_cox
      external rrchessian_cox, dbcoah

      allocate(xguess(nd))
      allocate(xlb(nd))
      allocate(xub(nd))
      allocate(xscale(nd))
      allocate(x(nd))
      allocate(hes_fval(nd,nd))
      allocate(hesinv(nd,nd))
      allocate(eval(nd))
      allocate(r(nd,nd))
      allocate(work(1000))
      allocate(ipiv(nd))

      xguess= 0.02d0 
c      print *, 'initial guess=', xguess
  
      ibtype=0
      xlb(1:2)=-1.0d0 
      xub(1:2)=1.0d0
      xlb(3:nd)=-2.0d0
      xub(3:nd)=2.0d0


      xscale=1.0d0
      fscale=1.0d0
      iparam(1)=0

c ----------------------------------------------------------------------------------
      iprint = 0
      call ffsqp(nd,1,0,0,0,0,100,iprint,500,info,1d+10,1d-4,0,0,
     *  xlb, xub, x, fval, gval, iwrkf, 5000, wrkf, 5000,
     *  obj,constr,gradobj,grcnfd)
c ----------------------------------------------------------------------------------

c     call dbcoah(rrccox_likelihood,rrcgradient_cox, rrchessian_cox
c    *     ,nd,xguess,ibtype,xlb,xub,xscale,fscale,iparam,rparam,x,
c    *     fvalue)
 
c     ier = iercd()
c     isev = n1rty(1)
c     ifail = 1
c      print *, 'isev=', isev

      isev = 0
      if (isev .eq. 0) then
         
         betaes=x
c        print *,'bcoah estimates=', betaes

         call rrchessian_cox(nd,betaes,hes_fval,nd)
c        print *, 'hessian_ana=', hes_fval

         do ir=1,nd
         do ic=ir,nd
           r(ir,ic) = hes_fval(ir,ic)
         enddo
         enddo

         nw = 1000
         call dsyev('N', 'U', nd, r, nd, eval, work, nw, info)

c        call devlsf(nd,hes_fval,nd,eval)
c         print *, 'eigenvalue=',eval
         if (minval(eval).gt.0.0d0) then
            ifail=0

            do ir=1,nd
            do ic=1,nd
              r(ir,ic) = hes_fval(ir,ic)
            enddo
            enddo

            call dpotrf('U', nd, r, nd, info)
            call dpotri('U', nd, r, nd, info)

            do ir=1,nd
            do ic=ir,nd
              hesinv(ir,ic) = r(ir,ic)
              hesinv(ic,ir) = hesinv(ir,ic)
            enddo
            enddo

c           call dlinds(nd,hes_fval,nd,hesinv,nd)
            Do i=1,nd
               st_error(i)=dsqrt(hesinv(i,i))
            End Do    
         endif

c         print *, 's.e.=', st_error
         
      endif

      deallocate(xguess)
      deallocate(xlb)
      deallocate(xub)
      deallocate(xscale)
      deallocate(x)
      deallocate(hes_fval)
      deallocate(hesinv)
      deallocate(eval)
      deallocate(r)
      deallocate(work)
      deallocate(ipiv)
 
      return
      end 

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c ----------------------------------------------------------------------------------
      subroutine obj(nparam,j,x,fj)
      implicit none
      integer nparam, j
      double precision x(nparam), fj
      call rrccox_likelihood(nparam,x,fj)
      return
      end

c ----------------------------------------------------------------------------------
      subroutine gradobj(nparam,j,x,gradfj,dummy)
      implicit none
      integer nparam, j
      double precision x(nparam), gradfj(nparam), dummy
      call rrcgradient_cox(nparam,x,gradfj)
      return
      end

c ----------------------------------------------------------------------------------
      subroutine constr(nparam,j,x,gj)
      implicit none
      integer nparam,j
      double precision x(nparam), gj
c     null function - required for compilation
      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine cox_zbren(nd,betaes,st_error,ifail)

      use variable
      implicit none

c     variables on the entry
      integer nd, ifail
      double precision betaes(nd), st_error(nd)
      
c     dummy variables
      integer maxfn, isev 
      real*8 dx,abmx,za(1),zb(1),errabs,errrel
      real*8 derv
      real*8 coxscore

      double precision, dimension(:,:), allocatable ::
     +  hess
      
      external coxscore, fzero
      real*8 rrr
      integer istat

      allocate(hess(nd,nd))

      dx=0.10d0
      abmx=0.1d0
      za = abmx-dx
      zb = abmx+dx
      do while (coxscore(za)*coxscore(zb).gt.0.0)
         za = za-dx
         zb = zb+dx
      end do
      
      errabs = 0.0d0
      errrel = 0.00001d0
      maxfn = 300
      call fzero(coxscore,za(1),zb(1),rrr,errrel,errabs,istat)

c     call dzbren(coxscore,errabs,errrel,za,zb,maxfn)

c     isev = n1rty(1)
c     ifail = 1
c     st_error=0.0
      isev = 0
      if (isev .eq. 0) then 
         betaes = zb
         call rrchessian_cox(nd,betaes,hess,nd)
         derv=hess(1,1)
         if (derv.gt.0.0) then
            ifail=0
            st_error=dsqrt(1.0d0/derv)
         endif
      endif

      deallocate(hess)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      double precision function coxscore(bval)
      use variable
      implicit none

c     variables on the entry
      real*8 bval(1), g_fval(1)

c     dummy variables
      integer nd
      nd=1

      call rrcgradient_cox(nd,bval,g_fval)
      coxscore=g_fval(1)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine rrc

      use variable
      implicit none

      integer pp,qq
      integer i,j,k,k1,i1,j1,j2,j3,ctnum,ir,ic,il,info
      integer ifail,isev, groupsize,newgroup, brn

      integer idbg, jdbg
      
      real*8 beta_rrc, st_rrc, st_rrcs, st_rrcs1

      double precision, dimension(:), allocatable ::
     +  aprime, rs
      double precision, dimension(:,:), allocatable ::
     +  g, vargam, s, gamma, temp1, temp2,
     +  stemp, matx, maty, matxy, matxinv, tempc, 
     +  tempx, tempy, tempsh, tempch, xtemp, ytemp, r
      double precision, dimension(:,:,:), allocatable ::
     +  acom1

      real*8 rhoi, tha
c      parameter (rhoi=0.3d0,tha=1.0d0)

      pp=nex
      qq=nex+ncov

      allocate(g(qq,pp))
      allocate(aprime(pp))
      allocate(vargam(pp*qq,pp*qq))
      allocate(s(pp,pp))
      allocate(rs(pp))
      allocate(gamma(qq+1,pp))
      allocate(temp1(pp,1))
      allocate(temp2(qq+1,1))
      allocate(stemp(pp,1))
      allocate(acom1(nageobs,1+nex+ncov,nex))
      allocate(matx(nageobs,nageobs))
      allocate(maty(nageobs,nageobs))
      allocate(matxy(nageobs,nageobs))
      allocate(matxinv(nageobs,nageobs))
      allocate(tempc(nageobs,1))
      allocate(tempx(nageobs,1))
      allocate(tempy(nageobs,1))
      allocate(tempsh(nageobs,1))
      allocate(tempch(nageobs,1))
      allocate(xtemp(valobs,qq))
      allocate(ytemp(valobs,pp))
      allocate(r(nageobs,nageobs))

      rhoi=icc
      tha=tau

c     Define the risk process ry_main and ry_val for each observation time
      ry_main=0
      i=1
      do j=1,nfa
         do while (t_main(i).lt.tf(j)) 
            if (i.lt.nmain) then
               i=i+1
            else
               ry_main(j)=nmain+1
               go to 10
            endif              
         enddo
         ry_main(j)=i
 10   enddo
      
      ry_val=0
      i=1
      do j=1,nfa
         do while (t_val(i).lt.tf(j)) 
            if (i.lt.nval) then
               i=i+1
            else 
               ry_val(j)=nval+1
               go to 20
            endif
         enddo
         ry_val(j)=i
 20   enddo



c     Calculate the regression coefficients in the risk sets of 
c     the validation study.
c     Grouping idea is using in this version . 

c      age_low(1)=tf(1)
c      age_up(1)=tf(vaintl)
c 
c      do i=2,ngrp-1
c         age_low(i)=tf(1+(i-1)*vaintl)
c         age_up(i)=tf(i*vaintl)
c      enddo
c
c      age_low(ngrp)=tf(1+(ngrp-1)*vaintl)
c      age_up(ngrp)=tf(nfa)


c
cc     age grouping
c
      groupsize=0
      ngrp=0
      newgroup=1
      do i=1,nageobs
         if (newgroup==1) then  
            age_low(ngrp+1)=tf(i)
            groupsize=0
c            print *,'age_low',ngrp+1,age_low(ngrp+1)
         endif
         do j1=1,nval
            do j2=1,j_val(j1)
               if (tf(i)==age_val(j1,j2)) then  
                  groupsize=groupsize+1
               endif
            enddo
         enddo
         if (groupsize>groupnum-1) then 
            ngrp=ngrp+1
            age_up(ngrp)=tf(i)
            newgroup=1
c            print *,'age_up',ngrp,age_up(ngrp)
         else if ((groupsize<=groupnum-1) .and. (i/=nageobs)) then
            newgroup=0
         else 
            newgroup=0
            age_up(ngrp)=tf(nageobs)
c            print *,'age_up',ngrp,age_up(ngrp),tf(i)
        
         endif
      enddo



c     Calculate the regression coefficients in the baseline set of
c     validation study

cc    For general linear measurement error model, we have derived another way to work it out.
c      if (ismodel.eq.1) then
cc     For general error model
c         i1=0
cc     the number of the baseline risk set
c         brn=1
c         do i=1,nval
c            if ((t_val(i)).ge.age_low(1)) then
c               do j1=1,j_val(i)
c                  if (((age_val(i,j1)).ge.age_low(1))
c     *              .and.((age_val(i,j1)).le.age_up(brn))
c     *              ) then
c                     i1=i1+1
c                     xtemp(i1,1:(nex+ncov))=
c     *                    zs_depval(i,j1,1:(nex+ncov))
c                     ytemp(i1,1:nex)=zt_depval(i,j1,1:nex)
c 
c                  endif
c               enddo
c            endif
c         enddo
c
c         call mvlm2(pp,qq,i1,ytemp(1:i1,:),xtemp(1:i1,:),g,rs)
c
c         print *, ''
c         print *, 'i1=', i1
c         print *, 'a0=', g(1,1)
c         print *, 'a1=', g(2,1)
c         print *, 'rs=', rs(1)
c
c         slp=g(2,1)
c      endif

      if (ismodel.eq.2) then
c     For additive error model
         slp=1.0d0
      endif


c     if the censoring age of a subject is not younger than age_low()
c     then this subject is in the risk set (age_low, age_up).

      id_val=0
c     special dealing with the first entry
      k=1
      i1=0
      do i=1,nval
         if ((t_val(i)).ge.age_low(k)) then
            do j1=1,j_val(i)
               if (((age_val(i,j1)).ge.age_low(k))
     *              .and.((age_val(i,j1)).le.age_up(k))
     *              .and.(zt_depval(i,j1,1).ne.-100.0d0)
     *              .and.(zs_depval(i,j1,1).ne.-100.0d0)
     *              ) then
                  i1=i1+1
                  xxcom(i1,k,1:(nex+ncov))=
     *                 zs_depval(i,j1,1:(nex+ncov))
                  yycom(i1,k,1:nex)=zt_depval(i,j1,1:nex)
c     add for variance derivation
                  xcom(i1,k)=xxcom(i1,k,1)
                  ycom(i1,k)=yycom(i1,k,1)

                  id_val(i1,k)=sub_id_val(i)
               endif
            enddo
         endif
      enddo
      count4(k)=i1
      

      do k=2, ngrp-1
         i1=0
         do i=1,nval
            if ((t_val(i)).ge.age_low(k)) then
               do j1=1,j_val(i)
                  if (((age_val(i,j1)).ge.age_low(k))
     *                 .and.((age_val(i,j1)).le.age_up(k))
     *                 .and.(zt_depval(i,j1,1).ne.-100.0d0)
     *                 .and.(zs_depval(i,j1,1).ne.-100.0d0)
     *                 ) then
                     i1=i1+1
                     xxcom(i1,k,1:(nex+ncov))=
     *                    zs_depval(i,j1,1:(nex+ncov))
                     yycom(i1,k,1:nex)=zt_depval(i,j1,1:nex)
c     add for variance derivation
                  xcom(i1,k)=xxcom(i1,k,1)
                  ycom(i1,k)=yycom(i1,k,1)

                     id_val(i1,k)=sub_id_val(i)
                  endif
               enddo
            endif
         enddo
         count4(k)=i1

      enddo


c     special dealing with the last entry
      k=ngrp
      i1=0
      do i=1,nval
         if ((t_val(i)).ge.age_low(k)) then
            do j1=1,j_val(i)
               if (((age_val(i,j1)).ge.age_low(k))
     *              .and.((age_val(i,j1)).le.age_up(k))
     *              .and.(zt_depval(i,j1,1).ne.-100.0d0)
     *              .and.(zs_depval(i,j1,1).ne.-100.0d0)
     *              ) then
                  i1=i1+1
                  xxcom(i1,k,1:(nex+ncov))=
     *                 zs_depval(i,j1,1:(nex+ncov))
                  yycom(i1,k,1:nex)=zt_depval(i,j1,1:nex)
c     add for variance derivation
                  xcom(i1,k)=xxcom(i1,k,1)
                  ycom(i1,k)=yycom(i1,k,1)

                  id_val(i1,k)=sub_id_val(i)
               endif
            enddo
         endif
      enddo
      count4(k)=i1

      cntcom=count4

      print *, 'count4=', count4(1:(ngrp+1))

c      do k=1,ngrp
c         if (count4(k).ge.2) then
c            print *,''
c            print *,'k=', k
c            print *,'age_low=',age_low(k)
c            print *,'age_up=',age_up(k)
c	    print *,'count4=',count4(k)
c            print *,'xxcom=',xxcom(1:count4(k),k,1:(nex+ncov))
c            print *,'yycom=',yycom(1:count4(k),k,1:(nex))
c	    print *, ''
c         endif
c      enddo

      do j=1,ngrp
         if (count4(j).ge.2) then
c     Estimate the first moment and the second moment instead of 
c     \alpha_0 and \alpha_1

            meanx(j)=0.0d0
            meany(j)=0.0d0
            varx(j)=0.0d0
            vary(j)=0.0d0
            covxy(j)=0.0d0
            do i=1,count4(j)
               meanx(j)=meanx(j)+xxcom(i,j,1)/dble(count4(j))
               meany(j)=meany(j)+yycom(i,j,1)/dble(count4(j))
            enddo
            
c            print *, 'meanx=', meanx(j)
c            print *, 'meany=', meany(j)

            do i=1,count4(j)
              varx(j)=varx(j)+(xxcom(i,j,1)-meanx(j))**2/dble(count4(j))
              vary(j)=vary(j)+(yycom(i,j,1)-meany(j))**2/dble(count4(j))
              covxy(j)=covxy(j)+(xxcom(i,j,1)-meanx(j))*
     *(yycom(i,j,1)-meany(j))/dble(count4(j))
            enddo

c            print *, 'varx=', varx(j)
c            print *, 'vary=', vary(j)
c            print *, 'covxy=',covxy(j)

            acom1(j,2,1)=covxy(j)/varx(j)
            acom1(j,1,1)=meany(j)-covxy(j)/varx(j)*meanx(j)
            
c            print *, 'acom1(j,k1)=', j, acom1(j,:,:)
            
            call mvlm2(pp,qq,count4(j),yycom(1:count4(j),j,:),
     *           xxcom(1:count4(j),j,:),gamma,rsquare(j,:))

            acom(j,:,:)=gamma 
c            print *, 'j,acom(j,k1)=', j, acom(j,:,:)
         endif
      enddo
     
c      stop

c     construct the x_hat value from the regression coefficients
c     acom(60,1+nex+ncov,nex)
      zshat_depmain=-100.0d0

      do j=1,nfa

         do k=1,ngrp
            if (((tf(j)).ge.age_low(k)).and.
     *           ((tf(j)).le.age_up(k))) then
               do i=ry_main(j),nmain
                  do j1=1,j_main(i)
                     stemp=0.0d0
                     if (age_main(i,j1).eq.tf(j)) then                        
                        do j2=1,j1
c     so far for one exposure only, no covariates
                           matx(j2,j2)=varx(k)
                           maty(j2,j2)=covxy(k)
                           do j3=1,j1
                              if (j3.lt.j2) then
                                 matx(j2,j3)=
     &   varx(k)*(rhoi**(dble(abs(j2-j3))**tha))
                                 matx(j3,j2)=matx(j2,j3)
                                 maty(j2,j3)=
     &   slp*varx(k)*(rhoi**(dble(abs(j2-j3))**tha))
                                 maty(j3,j2)=maty(j2,j3)
                              endif
                           enddo
                           tempc(j2,1)=zs_depmain(i,j2,1)
                           tempx(j2,1)=meanx(k)
                           tempy(j2,1)=meany(k)
                        enddo

c                        print *, 'matx=', matx(1:j1,1:j1)

                         do ir=1,j1
                         do ic=ir,j1
                           r(ir,ic) = matx(ir,ic)
                         enddo
                         enddo

                         call dpotrf('U', j1, r, nageobs, info)
                         call dpotri('U', j1, r, nageobs, info)

                         do ir=1,j1
                         do ic=ir,j1
                           matxinv(ir,ic) = r(ir,ic)
                           matxinv(ic,ir) = matxinv(ir,ic)
                         enddo
                         enddo

c                      call dlinds(j1,matx(1:j1,1:j1),j1,
c    *matxinv(1:j1,1:j1),j1) 
                       
                       isev = 0
                       if (isev.ne.0) then 
                          print *, 'ERROR',isev
                          print *, 'j1=', j1
                          print *, 'k=', k
                          print *, 'matx=', matx(1:j1,1:j1)
                          stop
                       endif

                       do ir=1,j1
                       do ic=1,j1
                         matxy(ir,ic) = 0.0
                         do il=1,j1
                           matxy(ir,ic) = matxy(ir,ic) +
     *                       maty(ir,il)*matxinv(il,ic)
                         enddo
                       enddo
                       enddo

c                      call dmrrrr(j1,j1,maty(1:j1,1:j1),j1,
c    *j1,j1,matxinv(1:j1,1:j1),j1,j1,j1,matxy(1:j1,1:j1),j1)

                       do ir=1,j1
                         tempch(ir,1) = 0.0
                         do il=1,j1
                           tempch(ir,1) = tempch(ir,1) +
     *                       matxy(ir,il)*tempc(il,1)
                         enddo
                       enddo

c                      call dmrrrr(j1,j1,matxy(1:j1,1:j1),j1,j1,1,
c    *tempc(1:j1,1),j1,j1,1,tempch(1:j1,1),j1)

                       do ir=1,j1
                         tempsh(ir,1) = 0.0
                         do il=1,j1
                           tempsh(ir,1) = tempsh(ir,1) +
     *                       matxy(ir,il)*tempx(il,1)
                         enddo
                       enddo


c                      call dmrrrr(j1,j1,matxy(1:j1,1:j1),j1,j1,1,
c    *tempx(1:j1,1),j1,j1,1,tempsh(1:j1,1),j1)


                       tempsh(1:j1,1)=tempy(1:j1,1)-tempsh(1:j1,1)
                       tempch(1:j1,1)=tempsh(1:j1,1)+tempch(1:j1,1)


                       do j2=1,j1
                          stemp(1,1)=stemp(1,1)+tempch(j2,1)
                       enddo

                       zshat_depmain(i,j1,:)=stemp(:,1)/dble(j1)
                     endif
                  enddo
               enddo
            endif
         enddo
 
      enddo


c      print *, 'rrc1'
c      do i=1,10
c         print *, 'i=',i
c         print *, 'age_main=', age_main(i,:)
c         print *, 'zshat_depmain=', zshat_depmain(i,:,:)
c      enddo
c      print *, 'rrc1'

      deallocate(g)
      deallocate(aprime)
      deallocate(vargam)
      deallocate(s)
      deallocate(rs)
      deallocate(gamma)
      deallocate(temp1)
      deallocate(temp2)
      deallocate(stemp)
      deallocate(acom1)
      deallocate(matx)
      deallocate(maty)
      deallocate(matxy)
      deallocate(matxinv)
      deallocate(tempc)
      deallocate(tempx)
      deallocate(tempy)
      deallocate(tempsh)
      deallocate(tempch)
      deallocate(xtemp)
      deallocate(ytemp)
      deallocate(r)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      
      subroutine rrccoxmethod
      
      use variable
      implicit none
         
      integer nd, ifail
      double precision, dimension(:), allocatable ::
     +  betaes, st_error, st_errors, st_errors1

      nd=nex+ncov
      allocate(betaes(nd)) 
      allocate(st_error(nd))
      allocate(st_errors(nd))
      allocate(st_errors1(nd))

      zs_common=zshat_depmain
      covar_common=covar_main

      if (nd.ge.2) then
         call cox_bcoah(nd,betaes,st_error,ifail)
      else if (nd.eq.1) then
         call cox_zbren(nd,betaes,st_error,ifail)
      else
         print *, 'The dimension of data is zero!'
      endif
      betaes_rrc=betaes

      if (ifail.eq.0) then
         call rrc_sderror(betaes, st_error, st_errors, st_errors1)
         serror_rrc=st_errors
c         print *, 'All variance=', serror_rrc
         print *, 'rrc','beta=', betaes_rrc
         print *, 'se=', serror_rrc
      else     
         print *,'The method does not converge.'
      endif

      deallocate(betaes) 
      deallocate(st_error)
      deallocate(st_errors)
      deallocate(st_errors1)
      
      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrc_sderror(bval,st_error, st_errors, st_errors1)
c     this subroutine is used to calculate the variance estimator of
c     rrc estimates. 

      use variable
      implicit none

c     variables on the entry
      integer nd
      double precision bval(nex+ncov),st_error(nex+ncov),
     +  st_errors(nex+ncov),st_errors1(nex+ncov)

c     dummy variables
      integer i,j, npa, ir, ic, il, info, ix, iy

      double precision, dimension(:), allocatable ::
     +  st_erbe1,st_erbe1s
      double precision, dimension(:,:), allocatable ::
     +  fchat,fbhat,fbinv,fbsand,ftmp,paramatrix,
     +  paratemp1,g_score,u2matrix,tempmatrix,
     +  st_erpsi,g_fval,tempmatrix2,r

      nd=nex+ncov
      npa=5

      allocate(fchat(nd,nd))
      allocate(fbhat(nd,nd))
      allocate(fbinv(nd,nd))
      allocate(fbsand(nd,nd))
      allocate(ftmp(nd,nd))
      allocate(paramatrix(nd,nd))
      allocate(paratemp1(nd,nd))
      allocate(st_erbe1(nd))
      allocate(st_erbe1s(nd))
      allocate(g_score(npa*ngrp,npa*ngrp))
      allocate(u2matrix(npa*ngrp,npa*ngrp))
      allocate(tempmatrix(npa*ngrp,npa*ngrp))
      allocate(st_erpsi(npa*ngrp,npa*ngrp))
      allocate(g_fval(npa*ngrp,nd))
      allocate(tempmatrix2(nd,npa*ngrp))
      allocate(r(nd,nd))



      call score2rrc(bval,fchat)
      call rrchessian_cox(nd,bval,fbhat,nd)

c     inverse hessian matrix part 

      do ir=1,nd
      do ic=ir,nd
        r(ir,ic) = fbhat(ir,ic)
      enddo
      enddo

      call dpotrf('U', nd, r, nd, info)
      call dpotri('U', nd, r, nd, info)

      do ir=1,nd
      do ic=ir,nd
        fbinv(ir,ic) = r(ir,ic)
        fbinv(ic,ir) = fbinv(ir,ic)
      enddo
      enddo
c     call dlinds(nd,fbhat,nd,fbinv,nd) 
      Do i=1,nd
         st_erbe1(i)=dsqrt(fbinv(i,i))
      End Do

c     pure sandwich variance estimator part
      do ir=1,nd
      do ic=1,nd
        ftmp(ir,ic) = 0
        do il = 1,nd
          ftmp(ir,ic) = ftmp(ir,ic) + fbinv(ir,il)*fchat(il,ic)
        enddo
      enddo
      enddo
c     call dmrrrr(nd,nd,fbinv,nd,nd,nd,fchat,nd,nd,nd,ftmp,nd)

      do ir=1,nd
      do ic=1,nd
        fbsand(ir,ic) = 0
        do il = 1,nd
          fbsand(ir,ic) = fbsand(ir,ic) + ftmp(ir,il)*fbinv(il,ic)
        enddo
      enddo
      enddo
c     call dmrrrr(nd,nd,ftmp,nd,nd,nd,fbinv,nd,nd,nd,fbsand,nd)
      Do i=1,nd 
         st_erbe1s(i)=dsqrt(fbsand(i,i)) 
      End Do 


c     test the derivative
c      call test_rrcscore_para(1,bval)
c      stop


c     parameters' variability part
      call parascore_grad(ngrp, st_erpsi)
      
      idbgstop = 1
      call grrcscore_para(ngrp,bval,g_fval)
c     print *, '2 g_fval=', g_fval

      do ir = 1,nd
      do ic = 1,npa*ngrp
        tempmatrix2(ir,ic) = 0.0
        do il = 1,npa*ngrp
          tempmatrix2(ir,ic) = tempmatrix2(ir,ic)  +
     *      g_fval(il,ir)*st_erpsi(il,ic)
        enddo
      enddo
      enddo
      
c     call dmxtyf(npa*ngrp,nd,g_fval,npa*ngrp,npa*ngrp,npa*ngrp,
c    *     st_erpsi,npa*ngrp,nd,npa*ngrp,tempmatrix2,nd)

      do ir=1,nd
      do ic=1,nd
        paramatrix(ir,ic) = 0
        do il = 1,npa*ngrp
          paramatrix(ir,ic) = paramatrix(ir,ic) +
     *        tempmatrix2(ir,il)*g_fval(il,ic)
        enddo
      enddo
      enddo

c     call dmrrrr(nd,npa*ngrp,tempmatrix2,nd,npa*ngrp,nd,g_fval,
c    *     npa*ngrp,nd,nd,paramatrix,nd)
      

      do ir=1,nd
      do ic=1,nd
        paratemp1(ir,ic) = 0
        do il = 1,nd
          paratemp1(ir,ic) = paratemp1(ir,ic) +
     *        fbinv(ir,il)*paramatrix(il,ic)
        enddo
      enddo
      enddo

c     call dmrrrr(nd,nd,fbinv,nd,nd,nd,paramatrix,nd,nd,nd,
c    *paratemp1,nd)

      do ir=1,nd
      do ic=1,nd
        paramatrix(ir,ic) = 0
        do il = 1,nd
          paramatrix(ir,ic) = paramatrix(ir,ic) +
     *        paratemp1(ir,il)*fbinv(il,ic)
        enddo
      enddo
      enddo

c     call dmrrrr(nd,nd,paratemp1,nd,nd,nd,fbinv,nd,nd,nd,
c    *paramatrix,nd)


c     three kinds of variance estimators
c     st_errors1: pure sandwich variance estimator
c     st_errors:  sandwich variance estimator + variability of parameters
c              :  the most accurate one.
c     st_error:   non-sandwich variance estimator + variability of 
c             :   parameters. not so mathematically correctly but acceptable
c             :   in this situation.

      do i=1,nd
         st_error(i)=dsqrt(fbinv(i,i)+paramatrix(i,i))
         st_errors(i)=dsqrt(fbsand(i,i)+paramatrix(i,i))
      enddo

      st_errors1=st_erbe1s

      print *, 'Classical variance=', st_erbe1
      print *, 'Pure sandwich variance=', st_errors1
      print *, 'non_sandwich variance + variability of estimator=', 
     *     st_error

      deallocate(fchat)
      deallocate(fbhat)
      deallocate(fbinv)
      deallocate(fbsand)
      deallocate(ftmp)
      deallocate(paramatrix)
      deallocate(paratemp1)
      deallocate(st_erbe1)
      deallocate(st_erbe1s)
      deallocate(g_score)
      deallocate(u2matrix)
      deallocate(tempmatrix)
      deallocate(st_erpsi)
      deallocate(g_fval)
      deallocate(tempmatrix2)
      deallocate(r)
      
      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine score2rrc(bval,fval)
c     This is the updated one with the second term suggested by David.
c     02/21/09, finished the multi-dimensional version.

      use variable  
      implicit none

c     variables on the entry
      integer nd
      real*8 bval(nex+ncov),fval(nex+ncov,nex+ncov)

c     dummy variables
      real*8 didl
      parameter (didl=1.0d-8)
      integer j, jtemp, i, k, j1, indk, il, ir, ic

      double precision, dimension(:), allocatable ::
     +  s0, s0a
      double precision, dimension(:,:), allocatable ::
     +  rrf1, s1, ftmp, fvaltmp, temp1, temp2, ss2
      double precision, dimension(:,:,:), allocatable ::
     +  rrf2
      double precision s0b, temp

      nd=nex+ncov
      allocate(rrf1(nmain,nageobs))
      allocate(rrf2(nmain,nageobs,nd))
      allocate(s0(nageobs))
      allocate(s1(nageobs,nex+ncov))
      allocate(s0a(nageobs))
      allocate(ss2(nmain,nex+ncov))
      allocate(ftmp(nd,1))
      allocate(fvaltmp(nd,nd))
      allocate(temp1(nd,1))
      allocate(temp2(nd,1))

      rrf1=0.0d0; rrf2=0.0d0
      temp1(:,1)=bval
      do j=1,nfa
         do i=ry_main(j),nmain
            do j1=1,j_main(i)
               if ((age_main(i,j1).eq.tf(j)).and.
     *              (zshat_depmain(i,j1,1).ne.-100.0d0)) then
                  if (ncov.ge.1) then
                     temp2(:,1)=(/zshat_depmain(i,j1,1:nex),
     *                    covar_main(i,j1,1:ncov)/)
                  else
                     temp2(1:nex,1)=zshat_depmain(i,j1,1:nex)
                  endif

                  temp = 0.0
                  do il=1,nd
                    temp = temp + temp1(il,1)*temp2(il,1)
                  enddo

c                 call dmxtyf(nd,1,temp1,nd,nd,1,temp2,nd,1,1,temp,1)

                  rrf1(i,j)=dexp(temp)
                  rrf2(i,j,1:nex)=rrf1(i,j)*zshat_depmain(i,j1,1:nex)
                  if (ncov.ge.1) then
                     rrf2(i,j,(nex+1):nd)=rrf1(i,j)*
     *                    covar_main(i,j1,1:ncov)
                  endif
               endif
            enddo
         enddo
      enddo


      s0=0.0d0
      s1=0.0d0
      do j=1,nfa
         do i=ry_main(j), nmain
            if (rrf1(i,j).ne.0.0d0) then
               s0(j)=s0(j)+rrf1(i,j)
               do k=1,nd
                  s1(j,k)=s1(j,k)+rrf2(i,j,k)
               enddo
            endif
         enddo
      enddo

 
      ss2=0.0d0
      Do i = 1,nmain
         Do j= 1,nfa
            if (rrf1(i,j).ne.0.0d0) then
               s0a(j) = dmax1(s0(j),didl)
               if (t_main(i).ge.tf(j)) then
                  do k=1,nd
                     ss2(i,k)=ss2(i,k)+rrf1(i,j)* 
     *   (rrf2(i,j,k)/rrf1(i,j)-s1(j,k)/s0a(j))/(dble(nmain)*s0a(j))
                  enddo
               endif
            endif
         End Do
      End Do
  

      fval=0.0d0
      do j=1,ndist         
         if (evcnt(j).ne.0.0d0) then
           
            do j1=1,nfa
               if (t_main(ilow(j)).eq.tf(j1)) indk=j1
            enddo
           
            do i=ilow(j),iup(j)
               if (rrf1(i,indk).ne.0.0d0) then
                  if (del_main(i).eq.1.0d0) then
                     do k=1,nd
                        ftmp(k,1)=
     * (rrf2(i,indk,k)/rrf1(i,indk)-ss2(i,k)-s1(indk,k)/s0a(indk))
                     enddo

                     do ir = 1,nd
                     do ic = 1,nd
                       fvaltmp(ir,ic) = ftmp(ir,1)*ftmp(ic,1)
                     enddo
                     enddo
c                    call dmxytf(nd,1,ftmp,nd,nd,1,ftmp,nd,nd,nd,
c    *                    fvaltmp,nd) 
                     fval=fval+fvaltmp
                  endif
               endif
            enddo

         endif        
      enddo


      Do i = 1,nmain
         if (del_main(i).eq.0.0) then
            do k=1,nd
               ftmp(k,1)=-ss2(i,k)
            enddo

            do ir = 1,nd
            do ic = 1,nd
              fvaltmp(ir,ic) = ftmp(ir,1)*ftmp(ic,1)
            enddo
            enddo
c           call dmxytf(nd,1,ftmp,nd,nd,1,ftmp,nd,nd,nd,fvaltmp,nd)
            fval=fval+fvaltmp
         endif
      End Do
     
      deallocate(rrf1)
      deallocate(rrf2)
      deallocate(s0)
      deallocate(s1)
      deallocate(s0a)
      deallocate(ss2)
      deallocate(ftmp)
      deallocate(fvaltmp)
      deallocate(temp1)
      deallocate(temp2)

      return
      end



ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine grrcscore_para(ng,bval,g_fval)

      use variable
      implicit none

      integer i,j,k,ix,iy
      integer nd,ng
      real*8 g_fval(ng*5,nex+ncov),bval(nex+ncov), da

      double precision, dimension(:), allocatable ::
     +  fval, fval1
      double precision, dimension(:,:), allocatable ::
     +  par

      nd=nex+ncov
      allocate(par(ng,5))
      allocate(fval(nd))
      allocate(fval1(nd))
      
      da=0.0001d0
      do j=1,ng
         par(j,1)=meanx(j)
         par(j,2)=meany(j)
         par(j,3)=varx(j)
         par(j,4)=vary(j)
         par(j,5)=covxy(j)
      enddo

      call scorerrc_para(ng,bval,par,fval)
      do k=1,ng
         do i=1,5
            par(k,i)=par(k,i)+da
            call scorerrc_para(ng,bval,par,fval1)
            g_fval((k-1)*5+i,:)=(fval1-fval)/da
            par(k,i)=par(k,i)-da
         enddo
      enddo
 
      deallocate(par)
      deallocate(fval)
      deallocate(fval1)

      return
      end
ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine scorerrc_para(ng,bval,par,fval)
   
      use variable  
      implicit none
   
      integer nd,ng, pp, qq
      Integer i, j, k,k1,i1,j1,j2,j3, isev, ir, ic, il, info
      double precision par(ng,5)
      double precision fval(nex+ncov),bval(nex+ncov)

      double precision, dimension(:,:), allocatable ::
     +  temp1,temp2,stemp,matx,maty,matxy,matxinv,tempc,
     +  tempx,tempy,tempsh,tempch,xtemp,ytemp,r
      double precision, dimension(:,:,:), allocatable ::
     +  zshat_tmp_po,zshat_tmp
      real*8 rhoi, tha
      integer ix, iy, iz

      nd=nex+ncov
      pp=nex
      qq=nex+ncov
      rhoi=icc
      tha=tau

      allocate(zshat_tmp_po(nmain,nage_main,nex))
      allocate(zshat_tmp(nmain,nage_main,nex))
      allocate(temp1(pp,1))
      allocate(temp2(qq+1,1))
      allocate(stemp(pp,1))
      allocate(matx(nageobs,nageobs))
      allocate(maty(nageobs,nageobs))
      allocate(matxy(nageobs,nageobs))
      allocate(matxinv(nageobs,nageobs))
      allocate(tempc(nageobs,1))
      allocate(tempx(nageobs,1))
      allocate(tempy(nageobs,1))
      allocate(tempsh(nageobs,1))
      allocate(tempch(nageobs,1))
      allocate(xtemp(valobs,qq))
      allocate(ytemp(valobs,pp))
      allocate(r(nageobs,nageobs))

c     construct the x_hat value from the regression coefficients
c     acom(60,1+nex+ncov,nex)
      zshat_tmp=-100.0d0

      do j=1,nfa

         do k=1,ngrp
            if (((tf(j)).ge.age_low(k)).and.
     *           ((tf(j)).le.age_up(k))) then
               do i=ry_main(j),nmain
                  do j1=1,j_main(i)
                      stemp=0.0d0
                     if (age_main(i,j1).eq.tf(j)) then                        
                        do j2=1,j1
c     so far for one exposure only, no covariates
                           matx(j2,j2)=par(k,3)
                           maty(j2,j2)=par(k,5)
                           do j3=1,j1
                              if (j3.lt.j2) then
                                 matx(j2,j3)=
     &  par(k,3)*(rhoi**(dble(abs(j2-j3))**tha))
                                 matx(j3,j2)=matx(j2,j3)
                                 maty(j2,j3)=
     &  slp*par(k,3)*(rhoi**(dble(abs(j2-j3))**tha))
                                 maty(j3,j2)=maty(j2,j3)
                              endif
                           enddo
                           tempc(j2,1)=zs_depmain(i,j2,1)
                           tempx(j2,1)=par(k,1)
                           tempy(j2,1)=par(k,2)
                        enddo

c                        print *, 'matx=', matx(1:j1,1:j1)

                         do ir=1,j1
                         do ic=ir,j1
                           r(ir,ic) = matx(ir,ic)
                         enddo
                         enddo

                         call dpotrf('U', j1, r, nageobs, info)
                         call dpotri('U', j1, r, nageobs, info)

                         do ir=1,j1
                         do ic=ir,j1
                           matxinv(ir,ic) = r(ir,ic)
                           matxinv(ic,ir) = matxinv(ir,ic)
                         enddo
                         enddo

c                      call dlinds(j1,matx(1:j1,1:j1),j1,
c    *matxinv(1:j1,1:j1),j1)

                       isev = 0
                       if (isev.ne.0) then 
                          print *, 'ERROR',isev
                          print *, 'j1=', j1
                          print *, 'k=', k
                          print *, 'matx=', matx(1:j1,1:j1)
                          stop
                       endif

                       do ir=1,j1
                       do ic=1,j1
                         matxy(ir,ic) = 0.0
                         do il=1,j1
                           matxy(ir,ic) = matxy(ir,ic) + 
     *                       maty(ir,il)*matxinv(il,ic)
                         enddo
                       enddo
                       enddo

c                      call dmrrrr(j1,j1,maty(1:j1,1:j1),j1,
c    *j1,j1,matxinv(1:j1,1:j1),j1,j1,j1,matxy(1:j1,1:j1),j1)

                       do ir=1,j1
                         tempch(ir,1) = 0.0
                         do il=1,j1
                           tempch(ir,1) = tempch(ir,1) + 
     *                       matxy(ir,il)*tempc(il,1)
                         enddo
                       enddo

c                      call dmrrrr(j1,j1,matxy(1:j1,1:j1),j1,j1,1,
c    *tempc(1:j1,1),j1,j1,1,tempch(1:j1,1),j1)

                       do ir=1,j1
                         tempsh(ir,1) = 0.0
                         do il=1,j1
                           tempsh(ir,1) = tempsh(ir,1) + 
     *                       matxy(ir,il)*tempx(il,1)
                         enddo
                       enddo

c                      call dmrrrr(j1,j1,matxy(1:j1,1:j1),j1,j1,1,
c    *tempx(1:j1,1),j1,j1,1,tempsh(1:j1,1),j1)
                       tempsh(1:j1,1)=tempy(1:j1,1)-tempsh(1:j1,1)
                       tempch(1:j1,1)=tempsh(1:j1,1)+tempch(1:j1,1)

                       do j2=1,j1
                          stemp(1,1)=stemp(1,1)+tempch(j2,1)
                       enddo

                       zshat_tmp(i,j1,:)=stemp(:,1)/dble(j1)
                     endif
                  enddo
               enddo
            endif
         enddo
 
      enddo

      zs_common=zshat_tmp
      covar_common=covar_main
      call rrcgradient_cox(nd,bval,fval)
    
      deallocate(zshat_tmp_po)
      deallocate(zshat_tmp)
      deallocate(temp1)
      deallocate(temp2)
      deallocate(stemp)
      deallocate(matx)
      deallocate(maty)
      deallocate(matxy)
      deallocate(matxinv)
      deallocate(tempc)
      deallocate(tempx)
      deallocate(tempy)
      deallocate(tempsh)
      deallocate(tempch)
      deallocate(xtemp)
      deallocate(ytemp)
      deallocate(r)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine parascore_grad(nd,st_erpsi)
c     this subroutine is used to calculate the variance of 
c     \psi=(\mu_C, \mu_c,\sigma^2_C,\sigma^2_c,\sigma_cC)
c     this way consists of stacking of matrices.
c     It's probably time-consuming and I want to compare with the 
c     results from the subroutine st_para.
      
      use variable
      implicit none
     
c     variables on the entry
      integer nd,npa
      real*8 st_erpsi(5*nd,5*nd)
     
c     dummy variables
      integer i,j,k,k1,k2,ir,ic,il
      integer g_p_, ldg_aa, ldg_ua
      integer info

      integer, dimension(:), allocatable :: ipiv
      double precision, dimension(:), allocatable ::
     +  aa, ua 
      double precision, dimension(:,:), allocatable ::
     +  g_aa,g_ua,g_uas,g_temp,uwork1,uwork2,temp1,temp2,eye
      double precision, dimension(:,:,:), allocatable ::
     +  g_tempj,uall
      double precision, dimension(:,:,:,:), allocatable ::
     +  u2matrix, uwork

      npa=5

      allocate(aa(npa))
      allocate(ua(npa))
      allocate(g_aa(npa,npa))
      allocate(g_ua(npa,npa))
      allocate(g_uas(npa,npa))
      allocate(temp1(1,1))
      allocate(temp2(1,1))
      allocate(g_temp(npa,npa))
      allocate(g_tempj(nd,npa,npa))
      allocate(uall(nd,nval,npa))
      allocate(u2matrix(nd,nd,npa,npa))
      allocate(uwork1(1,npa))
      allocate(uwork2(1,npa))
      allocate(uwork(nd,nd,npa,npa))
      allocate(eye(npa,npa))
      allocate(ipiv(npa))
      uall=0.0d0
      
      Do j = 1,nd
         g_ua=0.0d0

         aa(1)=meanx(j)
         aa(2)=meany(j)
         aa(3)=varx(j)
         aa(4)=vary(j)
         aa(5)=covxy(j)

         Do i= 1,cntcom(j)
            uall(j,i,1)=xcom(i,j)-meanx(j)
            uall(j,i,2)=ycom(i,j)-meany(j)
            uall(j,i,3)=xcom(i,j)**2-meanx(j)**2-varx(j)
            uall(j,i,4)=ycom(i,j)**2-meany(j)**2-vary(j)
            uall(j,i,5)=ycom(i,j)*xcom(i,j)-meany(j)*meanx(j)-covxy(j)
         Enddo

c     Calculate the derivative by hand, not using adifor

         do k=1,npa
            g_ua(k,k)=dble(cntcom(j))
         enddo
         
         g_ua(3,1)=2.0d0*dble(cntcom(j))*meanx(j)
         g_ua(4,2)=2.0d0*dble(cntcom(j))*meany(j)
         g_ua(5,1)=dble(cntcom(j))*meany(j)
         g_ua(5,2)=dble(cntcom(j))*meanx(j)
 
         g_temp=-g_ua

c     Compute the inverse of a real matrix
 
         do ir=1,npa
         do ic=1,npa
           if ( ir .eq. ic ) then
             eye(ir,ic) = 1.0
           else
             eye(ir,ic) = 0.0
           endif
         enddo
         enddo

         call dgesv( npa, npa, g_temp, npa, ipiv, eye, npa, info )

         do ir=1,npa
         do ic=1,npa
           g_tempj(j,ir,ic) = eye(ir,ic)
         enddo
         enddo
c        call dlinrg(npa,g_temp,npa,g_tempj(j,:,:),npa)
                     
      Enddo

c     This part calculates \sum U*U^t.
      Do j=1,nd
         Do i=j,nd
            u2matrix(j,i,:,:)=0.0d0
            Do k1=1,max0(cntcom(i),cntcom(j))
               Do k2=1,max0(cntcom(i),cntcom(j))
                  if ((id_val(k1,i).ne.0).and.(id_val(k2,j).ne.0).and.
     *                 (id_val(k1,i).eq.id_val(k2,j))) then
                     uwork1(1,:)=uall(j,k2,:)
                     uwork2(1,:)=uall(i,k1,:)

                     do ir=1,npa
                     do ic=1,npa
                        uwork(j,i,ir,ic) = uwork1(1,ir)*uwork2(1,ic)
                     enddo
                     enddo
c                    call dmxtyf(1,npa,uwork1,1,1,npa,uwork2,1,npa,npa,
c    *                    uwork(j,i,:,:) ,npa)

                     u2matrix(j,i,:,:)=u2matrix(j,i,:,:)+
     * uwork(j,i,:,:)
                  endif
               Enddo
            Enddo
         Enddo

         Do i=1,j-1
c           print *, 'i=', i, 'j=', j
c           read(*,*)

            do ir=1,npa
            do ic=1,npa
              u2matrix(j,i,ir,ic) = u2matrix(i,j,ic,ir)
            enddo
            enddo

c           call dtrnrr(npa,npa,u2matrix(i,j,:,:),npa,npa,npa,
c    *           u2matrix(j,i,:,:),npa)
         Enddo
      Enddo

c      do i=1,nd
c         do j =1,nd
c            print *, 'u2=',i,j,u2matrix(i,j,:,:)
c            read(*,*)
c         enddo
c      enddo


c     Multiplying the matrices together.
      Do i=1,nd
         Do j=i,nd

            do ir=1,npa
            do ic=1,npa
              uwork(i,j,ir,ic) = 0
              do il = 1,npa
                uwork(i,j,ir,ic) = uwork(i,j,ir,ic) +
     *              g_tempj(i,ir,il)*u2matrix(i,j,il,ic)
              enddo
            enddo
            enddo

c           call dmrrrr(npa,npa,g_tempj(i,:,:),npa,npa,npa,
c    *u2matrix(i,j,:,:),npa,npa,npa,uwork(i,j,:,:),npa)


            do ir=1,npa
            do ic=1,npa
              u2matrix(i,j,ir,ic) = 0.0
              do il=1,npa
                u2matrix(i,j,ir,ic) = u2matrix(i,j,ir,ic) +  
     *            uwork(i,j,ir,il)*g_tempj(j,ic,il)
              enddo
            enddo
            enddo

c           call dmxytf(npa,npa,uwork(i,j,:,:),npa,npa,npa,
c    *g_tempj(j,:,:),npa,npa,npa,u2matrix(i,j,:,:),npa)
         Enddo

         Do j=1,i-1
            do ir=1,npa
            do ic=1,npa
              u2matrix(i,j,ir,ic) = u2matrix(j,i,ic,ir)
            enddo
            enddo
c           call dtrnrr(npa,npa,u2matrix(j,i,:,:),npa,npa,npa,
c    *u2matrix(i,j,:,:),npa)
         Enddo
      Enddo

      Do i=1,nd
         Do j=1,nd
            st_erpsi(1+npa*(i-1):npa*i,1+npa*(j-1):npa*j)=
     *u2matrix(i,j,:,:)
         Enddo
      Enddo
      
      deallocate(aa)
      deallocate(ua)
      deallocate(g_aa)
      deallocate(g_ua)
      deallocate(g_uas)
      deallocate(temp1)
      deallocate(temp2)
      deallocate(g_temp)
      deallocate(g_tempj)
      deallocate(uall)
      deallocate(u2matrix)
      deallocate(uwork1)
      deallocate(uwork2)
      deallocate(uwork)
      deallocate(eye)

      return
      end

Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrccox_likelihood(nd,betval,fcn_ran)

      use variable
      implicit none

c     variables on the entry
      integer nd
      real*8 betval(nd),fcn_ran

      call rrc_cum_cox(nmain,nage_main,nd,betval,fcn_ran)
      fcn_ran=-fcn_ran
c      print *, '-2Log L=', 2.0d0*fcn_ran

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrc_cum_cox(nmax,ntm,nd,betval,fcn_ran)

      use variable
      implicit none

c     variables on the entry
      integer nmax, ntm, nd
      double precision betval(nd), fcn_ran

c     dummy variables
      integer i,j,k,j1,k1,indk,kch
      double precision eps
      parameter(eps=1.0d-8)

      double precision temp
      double precision, dimension(:), allocatable :: s0
      double precision, dimension(:,:), allocatable :: rrf

      allocate(rrf(nmax,ntm))
      allocate(s0(nageobs))

      do i=1,nmax
         do k=1,j_main(i)
            if (zs_common(i,k,1).ne.-100.0d0) then
               temp=0.0d0
               do j1=1,nex
                  temp=temp+betval(j1)*zs_common(i,k,j1)
               enddo
               if (ncov.ge.1) then
                  do j1=1,ncov
                     temp=temp+betval(nex+j1)*dble(covar_common(i,k,j1))
                  enddo
               endif
               rrf(i,k)=dexp(temp)
            endif
         enddo
      enddo
      
      do k=1,nfa
         s0(k)=0.0d0
         do i=1,nmax
            do j=1,j_main(i)
               if ((age_main(i,j).eq.tf(k))
     *              .and.(zs_common(i,j,1).ne.-100.0d0)
     *              ) then
                  s0(k)=s0(k)+rrf(i,j)
               endif
            enddo
         enddo
      enddo


      fcn_ran=0.0d0
      do i=1,ndist

         if (evcnt(i).ne.0.0d0) then
            do j1=1,nfa
               if (t_main(ilow(i)).eq.tf(j1)) indk=j1
            enddo

            do j=ilow(i),iup(i)
               if ((del_main(j).eq.1.0d0)
     *              .and.(zs_common(j,j_main(j),1).ne.-100.0d0)
     *              )then
                  do k1=1,nex
                     fcn_ran=fcn_ran+betval(k1)*
     *                    zs_common(j,j_main(j),k1)
                  enddo

                  if (ncov.ge.1.0d0) then
                     do k1=1,ncov
                        fcn_ran=fcn_ran+betval(nex+k1)*
     *                       dble(covar_common(j,j_main(j),k1))
                     enddo
                  endif
               endif
            enddo
            fcn_ran=fcn_ran-evcnt(i)*dlog(dmax1(s0(indk),eps))
            
         endif

      enddo

      deallocate(rrf)
      deallocate(s0)
      
      return
      end

ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrcgradient_cox(nd,betval,g_fval)

      use variable
      implicit none

c     variables on the entry
      integer nd
      real*8 betval(nd),g_fval(nd)

c     dummy variables
      integer g_p_, ldg_betval,ldg_fcn_ran, i

      double precision, dimension(:), allocatable :: g_fcn_ran
      double precision, dimension(:,:), allocatable :: g_betval
      double precision fcn_ran

      g_p_=nex+ncov
      ldg_betval=nex+ncov
      ldg_fcn_ran=nex+ncov

      allocate(g_betval(ldg_betval,nd))
      allocate(g_fcn_ran(ldg_fcn_ran))

      g_betval=0.0d0
      do i=1,g_p_
         g_betval(i,i)=1.0d0
      enddo

      call  g_rrccox_likelihood(g_p_,nd,betval, g_betval, ldg_betval, fc
     *n_ran, g_fcn_ran, ldg_fcn_ran)
      g_fval=g_fcn_ran

      if ( idbgstop .eq. 1 ) then
c      print *, 'Gradient=', g_fval
      endif
c      print *, '-2Log L=', 2.0d0*fcn_ran

      deallocate(g_betval)
      deallocate(g_fcn_ran)

      return
      end

C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_rrccox_likelihood(g_p_, nd, betval, g_betval, 
     +  ldg_betval, fcn_ran, g_fcn_ran, ldg_fcn_ran)
C
      use variable
      implicit none

C     variables on the entry
        integer g_p_, nd, ldg_betval, ldg_fcn_ran 
        double precision betval(nd), g_betval(ldg_betval, nd), fcn_ran
        double precision g_fcn_ran(ldg_fcn_ran) 
C
        integer g_pmax_, g_i_
        external g_rrc_cum_cox

        integer g_ehfid
        data g_ehfid /0/
C
        g_pmax_ = nex+ncov
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        call g_rrc_cum_cox(g_p_, nmain, nage_main, betval, g_betval, 
     *ldg_betval, fcn_ran, g_fcn_ran, ldg_fcn_ran)
        do g_i_ = 1, g_p_
          g_fcn_ran(g_i_) = -g_fcn_ran(g_i_)
        enddo
        fcn_ran = -fcn_ran
C--------
C
        return
      end
C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_rrc_cum_cox(g_p_, nmax, ntm, betval, g_betval, 
     * ldg_betval, fcn_ran, g_fcn_ran, ldg_fcn_ran)
C
      use variable
      implicit none

C     variables on the entry
      integer g_p_, nmax, ntm, ldg_betval, ldg_fcn_ran
      double precision betval(nex + ncov), fcn_ran, 
     +  g_betval(ldg_betval,nex+ncov),g_fcn_ran(ldg_fcn_ran)
      
C
C     dummy variables

      integer i, j, k, j1, k1, indk, kch, ix, iy
      double precision eps, temp
      parameter (eps = 1.0d-8)

      double precision, dimension(:), allocatable :: s0
      double precision, dimension(:,:), allocatable :: rrf

      double precision, dimension(:), allocatable, save :: g_temp 
      double precision, dimension(:), allocatable, save :: g_d1_w
      double precision, dimension(:,:), allocatable, save :: g_s0 
      double precision, dimension(:,:,:), allocatable, save :: g_rrf 
C
      integer g_i_, g_pmax_
      double precision d5_b, d2_v, d2_p, d1_w, d1_p, d3_v


      integer g_ehfid
      data g_ehfid /0/
C

       g_pmax_ = nex+ncov
       allocate(rrf(nmain, nage_main))
       allocate(s0(nageobs))
       if ( .not. allocated(g_temp) )  then
         allocate(g_temp(g_pmax_))
         allocate(g_rrf(g_pmax_, nmain, nage_main))
         allocate(g_s0(g_pmax_, nageobs))
         allocate(g_d1_w(g_pmax_))
       endif 

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do i = 1, nmax
          do k = 1, j_main(i)
            if (zs_common(i, k, 1) .ne. (-100.0d0)) then
              do g_i_ = 1, g_p_
                g_temp(g_i_) = 0.0d0
              enddo
              temp = 0.0d0
C--------

              do j1 = 1, nex
                do g_i_ = 1, g_p_
                  g_temp(g_i_) = zs_common(i, k, j1) * 
     *g_betval(g_i_, j1) + g_temp(g_i_)
                enddo
                temp = temp + betval(j1) * zs_common(i, k, j1)
C--------
              enddo

              if (ncov .ge. 1) then
                do j1 = 1, ncov
                  do g_i_ = 1, g_p_
                    g_temp(g_i_) = covar_common(i, k, j1) * 
     *g_betval(g_i_, nex + j1) + g_temp(g_i_)
                  enddo
                  temp = temp + betval(nex + j1)*covar_common(i, k, j1)
C--------
                enddo
              endif
              d2_v = exp(temp)
              d1_p =  d2_v
              do g_i_ = 1, g_p_
                g_rrf(g_i_, i, k) = d1_p * g_temp(g_i_)
              enddo
              rrf(i, k) = d2_v

C--------
            endif
          enddo
        enddo
C
        do k = 1, nfa
          do g_i_ = 1, g_p_
            g_s0(g_i_, k) = 0.0d0
          enddo
          s0(k) = 0.0d0
C--------
          do i = 1, nmax
            do j = 1, j_main(i)
              if ((age_main(i, j) .eq. tf(k)) .and. 
     * (zs_common(i, j, 1) .ne. (-100.0d0))) then
                do g_i_ = 1, g_p_
                  g_s0(g_i_, k) = g_rrf(g_i_, i, j) + g_s0(g_i_, k)
                enddo
                s0(k) = s0(k) + rrf(i, j)
C--------
              endif
            enddo
          enddo
        enddo
C
C
        do g_i_ = 1, g_p_
          g_fcn_ran(g_i_) = 0.0d0
        enddo
        fcn_ran = 0.0d0
C--------
        do i = 1, ndist
C
          if (evcnt(i) .ne. 0.0d0) then
            do j1 = 1, nfa
              if (t_main(ilow(i)) .eq. tf(j1)) then
                indk = j1
              endif
            enddo
C
            do j = ilow(i), iup(i)
              if ((del_main(j) .eq. 1.0d0) .and. 
     *(zs_common(j, j_main(j), 1) .ne. (-100.0d0))) then
                do k1 = 1, nex
                  do g_i_ = 1, g_p_
                    g_fcn_ran(g_i_) = zs_common(j, j_main(j), k1) * 
     *g_betval(g_i_, k1) + g_fcn_ran(g_i_)
                  enddo
                  fcn_ran = fcn_ran + betval(k1) * 
     *zs_common(j, j_main(j), k1)
C--------
                enddo
C
                if (ncov .ge. 1.0d0) then
                  do k1 = 1, ncov
                    do g_i_ = 1, g_p_
                      g_fcn_ran(g_i_) = covar_common(j, j_main(j), k1) * 
     *g_betval(g_i_, nex + k1) + g_fcn_ran(g_i_)
                    enddo
                    fcn_ran = fcn_ran + betval(nex + k1) * 
     *covar_common(j, j_main(j), k1)
C--------
                  enddo
                endif
              endif
            enddo
            d2_v = max (s0(indk), eps)
            if (s0(indk) .gt.  eps) then
               d1_p = 1.0d0
               d2_p = 0.0d0
            else if (s0(indk) .lt.  eps) then
               d1_p = 0.0d0
               d2_p = 1.0d0
            else
               d1_p = 0.5d0
               d2_p = 0.5d0
            endif
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d1_p * g_s0(g_i_, indk)
            enddo
            d1_w = d2_v
            d3_v = log(d1_w)
            d1_p = 1.0d0 / d1_w
            d5_b = -evcnt(i) * d1_p
            do g_i_ = 1, g_p_
              g_fcn_ran(g_i_) = d5_b * g_d1_w(g_i_) + g_fcn_ran(g_i_)
            enddo
            fcn_ran = fcn_ran - evcnt(i) * d3_v
C--------
C
          endif
C
        enddo
C
      deallocate(rrf)
      deallocate(s0)
      return
      end
C
C
Ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrchessian_cox(nd,betval,hess_fval,ldh)

      use variable
      implicit none

c     variables on the entry
      integer nd,ldh
      real*8 betval(nd), hess_fval(ldh,nd)
      
c     dummy variables
      integer g__p_, ldg__betval, ldg_g_fval,i

      double precision, dimension(:), allocatable ::
     +  g_fval
      double precision, dimension(:,:), allocatable ::
     +  g__betval, g_g_fval

      g__p_=nex+ncov
      ldg__betval=nex+ncov
      ldg_g_fval=nex+ncov

      allocate(g_fval(nd))
      allocate(g_g_fval(ldg_g_fval, nd))
      allocate(g__betval(ldg__betval,nd))

      g__betval=0.0d0
      do i=1,g__p_
         g__betval(i,i)=1.0d0
      enddo
      call g_rrcgradient_cox(g__p_,nd,betval, g__betval, ldg__betval
     *, g_fval, g_g_fval, ldg_g_fval)
      hess_fval=g_g_fval

      deallocate(g_fval)
      deallocate(g_g_fval)
      deallocate(g__betval)

      return
      end

C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_rrcgradient_cox(g__p_, nd, betval, g__betval, 
     +  ldg__betval, g_fval, g_g_fval, ldg_g_fval)
C
      use variable
      implicit none
C
C     variables on the entry

      integer nd, g__p_, ldg__betval, ldg_g_fval
      double precision betval(nd),g__betval(ldg__betval, nd),g_fval(nd),
     +  g_g_fval(ldg_g_fval, nd)

C
C     dummy variables

      integer g_p_, ldg_betval, ldg_fcn_ran, i, j
      real*8 fcn_ran
C
      integer g_pmax_
      integer g_i_

      double precision, dimension(:), allocatable, save ::
     + g_fcn_ran
      double precision, dimension(:,:), allocatable, save ::
     +  g_g_fcn_ran, g_betval

      external g_g_rrccox_likelihood
      integer g_ehfid
      data g_ehfid /0/
C
      g_p_=nex+ncov
      ldg_betval=nex+ncov
      ldg_fcn_ran=nex+ncov
      g_pmax_ = nex+ncov

      allocate(g_betval(ldg_betval, ldg_betval))
      allocate(g_fcn_ran(ldg_fcn_ran))
      if ( .not. allocated(g_g_fcn_ran) ) then
        allocate(g_g_fcn_ran(g_pmax_, ldg_fcn_ran))
      endif

C
        if (g__p_ .gt. g_pmax_) then
          print *, 'Parameter g__p_ is greater than g_pmax_'
          stop
        endif
        do i = 1, g_p_
          do j = 1, g_p_
            if (i .ne. j) then
              g_betval(i, j) = 0.0d0
            else
              g_betval(i, j) = 1.0d0
            endif
          enddo
        enddo
C
        call g_g_rrccox_likelihood(g__p_, g_p_, nd, betval, g__betval, l
     *dg__betval, g_betval, ldg_betval, fcn_ran, g_fcn_ran, g_g_fcn_ran,
     * g_pmax_, ldg_fcn_ran)
C
        do i = 1, g_p_
          do g_i_ = 1, g__p_
            g_g_fval(g_i_, i) = g_g_fcn_ran(g_i_, i)
          enddo
          g_fval(i) = g_fcn_ran(i)
C--------
        enddo
C
      deallocate(g_betval)
      deallocate(g_fcn_ran)
        return
      end
C
C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_g_rrccox_likelihood(g__p_, g_p_, nd, betval, 
     +  g__betval,ldg__betval,g_betval,ldg_betval,fcn_ran,g_fcn_ran, 
     +  g_g_fcn_ran, ldg_g_fcn_ran, ldg_fcn_ran)
C
      use variable
      implicit none

C     variables on the entry
      integer g__p_, g_p_, nd, ldg__betval, ldg_g_fcn_ran, ldg_fcn_ran
      double precision betval(nd), fcn_ran, g__betval, 
     +  g_fcn_ran(ldg_fcn_ran), g_g_fcn_ran(ldg_g_fcn_ran,ldg_fcn_ran),
     +  g_betval(ldg_betval, nd)
C
      integer g_pmax_,g_i_, ldg_betval,g__pmax_,g__i_

      external g_g_rrc_cum_cox

      integer g_ehfid
      data g_ehfid /0/

      g_pmax_ = nex+ncov
      g__pmax_ = nex+ncov
C
        if (g__p_ .gt. g__pmax_) then
          print *, 'Parameter g__p_ is greater than g__pmax_'
          stop
        endif
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        call g_g_rrc_cum_cox(g__p_, g_p_, nmain, nage_main,
     * betval, g__betval, ldg__betval, g_betval, ldg_betval, fcn_ran, 
     *g_fcn_ran, g_g_fcn_ran, ldg_g_fcn_ran, ldg_fcn_ran)
        do g_i_ = 1, g_p_
          do g__i_ = 1, g__p_
            g_g_fcn_ran(g__i_, g_i_) = -g_g_fcn_ran(g__i_, g_i_)
          enddo
          g_fcn_ran(g_i_) = -g_fcn_ran(g_i_)
C--------
        enddo
        fcn_ran = -fcn_ran
C--------

      return
      end
C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_g_rrc_cum_cox(g__p_, g_p_, nmax, ntm, 
     *  betval, g__betval, ldg__betval, g_betval, ldg_betval, 
     *  fcn_ran, g_fcn_ran, g_g_fcn_ran, ldg_g_fcn_ran, ldg_fcn_ran)
C
      use variable
      implicit none

C     variables on the entry
      integer g__p_, g_p_, nmax, ntm, ldg__betval, ldg_betval, 
     *  ldg_g_fcn_ran, ldg_fcn_ran

      double precision betval(nex + ncov), fcn_ran,
     *  g__betval(ldg__betval, nex + ncov), 
     *  g_betval(ldg_betval, nex + ncov), g_fcn_ran(ldg_fcn_ran), 
     *  g_g_fcn_ran(ldg_g_fcn_ran, ldg_fcn_ran)

C     dummy variables

      integer i, j, k, j1, k1, indk, kch
      integer g_pmax_, g_i_,g__pmax_,g__i_,g_ehfid
      double precision eps, temp
      double precision d5_b, d2_v, d2_p, d1_w, d1_p, d3_v

      double precision, dimension(:), allocatable ::
     +  s0

      double precision, dimension(:), allocatable, save ::
     +  g_temp, g_d1_w
      double precision, dimension(:,:), allocatable, save ::
     +  g_s0, rrf 
      double precision, dimension(:,:,:), allocatable, save ::
     +  g_rrf
C
      double precision d2_b, d4_v, d4_p, d3_p
      double precision, dimension(:), allocatable, save ::
     *  g__temp, g_d2_v, g_d1_p, g__d1_w, g_d5_b
      double precision, dimension(:,:), allocatable, save ::
     *  g__s0, g_g_d1_w
      double precision, dimension(:,:,:), allocatable, save ::
     *  g__rrf,g_g_s0
      double precision, dimension(:,:,:,:), allocatable, save ::
     *  g_g_rrf

C
      data g_ehfid /0/
      eps = 1.0d-8
      g_pmax_ = nex+ncov
      g__pmax_ = nex+ncov

      allocate(rrf(nmain, nage_main)) 
      allocate(s0(nageobs))
      if ( .not. allocated(g_temp) ) then
        allocate(g_temp(g_pmax_))
        allocate(g_rrf(g_pmax_, nmain, nage_main))
        allocate(g_s0(g_pmax_, nageobs))
        allocate(g_d1_w(g_pmax_))
        allocate(g__temp(g__pmax_))
        allocate(g_d2_v(g__pmax_))
        allocate(g_d1_p(g__pmax_))
        allocate(g_g_rrf(g__pmax_, g_pmax_, nmain, nage_main))
        allocate(g__rrf(g__pmax_, nmain, nage_main))
        allocate(g_g_s0(g__pmax_, g_pmax_, nageobs))
        allocate(g__s0(g__pmax_, nageobs))
        allocate(g_g_d1_w(g__pmax_, g_pmax_))
        allocate(g__d1_w(g__pmax_))
        allocate(g_d5_b(g__pmax_))
      endif
C
        if (g__p_ .gt. g__pmax_) then
          print *, 'Parameter g__p_ is greater than g__pmax_'
          stop
        endif
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do i = 1, nmax
          do k = 1, j_main(i)
            if (zs_common(i, k, 1) .ne. (-100.0d0)) then
              do g_i_ = 1, g_p_
                g_temp(g_i_) = 0.0d0
              enddo
              do g__i_ = 1, g__p_
                g__temp(g__i_) = 0.0d0
              enddo
              temp = 0.0d0
C--------
C--------
              do j1 = 1, nex
                do g_i_ = 1, g_p_
                  g_temp(g_i_) = zs_common(i, k, j1) * 
     *g_betval(g_i_, j1) + g_temp(g_i_)
                enddo
                do g__i_ = 1, g__p_
                  g__temp(g__i_) = zs_common(i, k, j1) * 
     *g__betval(g__i_, j1) + g__temp(g__i_)
                enddo
                temp = temp + betval(j1) * zs_common(i, k, j1)
C--------
C--------
              enddo
              if (ncov .ge. 1) then
                do j1 = 1, ncov
                  do g_i_ = 1, g_p_
                    g_temp(g_i_) = covar_common(i, k, j1) * 
     *g_betval(g_i_, nex + j1) + g_temp(g_i_)
                  enddo
                  do g__i_ = 1, g__p_
                    g__temp(g__i_) = covar_common(i, k, j1) * 
     *g__betval(g__i_, nex + j1) + g__temp(g__i_)
                  enddo
                  temp = temp + betval(nex + j1)*covar_common(i, k, j1)
C--------
C--------
                enddo
              endif
              d4_v = exp(temp)
              d3_p =  d4_v
              do g__i_ = 1, g__p_
                g_d2_v(g__i_) = d3_p * g__temp(g__i_)
              enddo
              d2_v = d4_v
C--------
              do g__i_ = 1, g__p_
                g_d1_p(g__i_) = g_d2_v(g__i_)
              enddo
              d1_p = d2_v
C--------
              do g_i_ = 1, g_p_
                do g__i_ = 1, g__p_
                  g_g_rrf(g__i_, g_i_, i, k) = g_temp(g_i_) * g_d1_p(g__
     *i_)
                enddo
                g_rrf(g_i_, i, k) = d1_p * g_temp(g_i_)
C--------
              enddo
              do g__i_ = 1, g__p_
                g__rrf(g__i_, i, k) = g_d2_v(g__i_)
              enddo
              rrf(i, k) = d2_v
C--------
C--------
            endif
          enddo
        enddo
C
        do k = 1, nfa
          do g_i_ = 1, g_p_
            do g__i_ = 1, g__p_
              g_g_s0(g__i_, g_i_, k) = 0.0d0
            enddo
            g_s0(g_i_, k) = 0.0d0
C--------
          enddo
          do g__i_ = 1, g__p_
            g__s0(g__i_, k) = 0.0d0
          enddo
          s0(k) = 0.0d0
C--------
C--------
          do i = 1, nmax
            do j = 1, j_main(i)
              if ((age_main(i, j) .eq. tf(k)) .and. 
     *(zs_common(i, j, 1) .ne. (-100.0d0))) then
                do g_i_ = 1, g_p_
                  do g__i_ = 1, g__p_
                    g_g_s0(g__i_, g_i_, k) = g_g_s0(g__i_, g_i_, k) + g_
     *g_rrf(g__i_, g_i_, i, j)
                  enddo
                  g_s0(g_i_, k) = g_rrf(g_i_, i, j) + g_s0(g_i_, k)
C--------
                enddo
                do g__i_ = 1, g__p_
                  g__s0(g__i_, k) = g__rrf(g__i_, i, j) + g__s0(g__i_, k
     *)
                enddo
                s0(k) = s0(k) + rrf(i, j)
C--------
C--------
              endif
            enddo
          enddo
        enddo
C
C
        do g_i_ = 1, g_p_
          do g__i_ = 1, g__p_
            g_g_fcn_ran(g__i_, g_i_) = 0.0d0
          enddo
          g_fcn_ran(g_i_) = 0.0d0
C--------
        enddo
        fcn_ran = 0.0d0
C--------
        do i = 1, ndist
C
          if (evcnt(i) .ne. 0.0d0) then
            do j1 = 1, nfa
              if (t_main(ilow(i)) .eq. tf(j1)) then
                indk = j1
              endif
            enddo
C
            do j = ilow(i), iup(i)
              if ((del_main(j) .eq. 1.0d0) .and. 
     *(zs_common(j, j_main(j), 1) .ne. (-100.0d0))) then
                do k1 = 1, nex
                  do g_i_ = 1, g_p_
                    do g__i_ = 1, g__p_
                      g_g_fcn_ran(g__i_, g_i_) = g_g_fcn_ran(g__i_, g_i_
     *)
                    enddo
                    g_fcn_ran(g_i_) = zs_common(j, j_main(j), k1) * 
     *g_betval(g_i_, k1) + g_fcn_ran(g_i_)
C--------
                  enddo
                  fcn_ran = fcn_ran + betval(k1) * 
     *zs_common(j, j_main(j), k1)
C--------
                enddo
C
                if (ncov .ge. 1.0d0) then
                  do k1 = 1, ncov
                    do g_i_ = 1, g_p_
                      do g__i_ = 1, g__p_
                        g_g_fcn_ran(g__i_, g_i_) = g_g_fcn_ran(g__i_, g_
     *i_)
                      enddo
                      g_fcn_ran(g_i_) = covar_common(j, j_main(j), k1) * 
     *g_betval(g_i_, nex + k1) + g_fcn_ran(g_i_)
C--------
                    enddo
                    fcn_ran = fcn_ran + betval(nex + k1) * 
     *covar_common(j, j_main(j), k1)
C--------
                  enddo
                endif
              endif
            enddo
            d4_v = max (s0(indk), eps)
            if (s0(indk) .gt.  eps) then
               d3_p = 1.0d0
               d4_p = 0.0d0
            else if (s0(indk) .lt.  eps) then
               d3_p = 0.0d0
               d4_p = 1.0d0
            else
               d3_p = 0.5d0
               d4_p = 0.5d0
            endif
            do g__i_ = 1, g__p_
              g_d2_v(g__i_) = d3_p * g__s0(g__i_, indk)
            enddo
            d2_v = d4_v
C--------
            if (s0(indk) .gt. eps) then
              do g__i_ = 1, g__p_
                g_d1_p(g__i_) = 0.0d0
              enddo
              d1_p = 1.0d0
C--------
              d2_p = 0.0d0
            else
              if (s0(indk) .lt. eps) then
                do g__i_ = 1, g__p_
                  g_d1_p(g__i_) = 0.0d0
                enddo
                d1_p = 0.0d0
C--------
                d2_p = 1.0d0
              else
                do g__i_ = 1, g__p_
                  g_d1_p(g__i_) = 0.0d0
                enddo
                d1_p = 0.5d0
C--------
                d2_p = 0.5d0
              endif
            endif
            do g_i_ = 1, g_p_
              do g__i_ = 1, g__p_
                g_g_d1_w(g__i_, g_i_) = d1_p * g_g_s0(g__i_, g_i_, indk)
     * + g_s0(g_i_, indk) * g_d1_p(g__i_)
              enddo
              g_d1_w(g_i_) = d1_p * g_s0(g_i_, indk)
C--------
            enddo
            do g__i_ = 1, g__p_
              g__d1_w(g__i_) = g_d2_v(g__i_)
            enddo
            d1_w = d2_v
C--------
            d3_v = log(d1_w)
            d4_v = 1.0d0 / d1_w
            d2_b = -d4_v / d1_w
            do g__i_ = 1, g__p_
              g_d1_p(g__i_) = d2_b * g__d1_w(g__i_)
            enddo
            d1_p = d4_v
C--------
            do g__i_ = 1, g__p_
              g_d5_b(g__i_) = -evcnt(i) * g_d1_p(g__i_)
            enddo
            d5_b = -evcnt(i) * d1_p
C--------
            do g_i_ = 1, g_p_
              do g__i_ = 1, g__p_
                g_g_fcn_ran(g__i_, g_i_) = g_g_fcn_ran(g__i_, g_i_) + d5
     *_b * g_g_d1_w(g__i_, g_i_) + g_d1_w(g_i_) * g_d5_b(g__i_)
              enddo
              g_fcn_ran(g_i_) = d5_b * g_d1_w(g_i_) + g_fcn_ran(g_i_)
C--------
            enddo
            fcn_ran = fcn_ran - evcnt(i) * d3_v
C--------
C
          endif
C
        enddo
C
      deallocate(rrf)
      deallocate(s0)
      return
      end
C
C
Ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


	subroutine mvlm2(p,q,n,y,xin,gamma,r)
        implicit none
c
c this routine computes the parameter estimates and their variance-covariance
c matrix for the multivariate linear model of p-vector y on (p+q)-vector xin
c the slopes are in (p+q)xp matrix g; aprim has the intercepts
c
c	n = number of observations in the dataset
c
c	p = dim(Y) ! which are the true covariates 
c
c	q = dim(X) ! which are the error measured covariates and other explanatory variable
c
      integer p,q,n
      double precision y(n,p),xin(n,q),gamma(q+1,p),r(p)

      integer, dimension(:), allocatable :: ipiv
      double precision, dimension(:), allocatable ::
     +  aprime, xbar, ybar
      double precision, dimension(:,:), allocatable ::
     +  g, x, xtx, xix, vargam, e, proj, yp, s, xtx2, temp, isubn, eye

      double precision xsum,ysum,part1,part2,part3
      integer q1,ir,ic,il,i,j,info

      allocate(g(q,p))
      allocate(aprime(p))
      allocate(x(n,q+1))
      allocate(xbar(q))
      allocate(ybar(p))
      allocate(xtx(q+1,q+1))
      allocate(xix(q+1,n))
      allocate(vargam(p*q,p*q))
      allocate(e(n,n))
      allocate(proj(n,n))
      allocate(yp(p,n))
      allocate(s(p,p))
      allocate(xtx2(q,q))
      allocate(temp(q,q))
      allocate(isubn(n,n))
      allocate(ipiv(q+1))
      allocate(eye(q+1,q+1))

c
c add 1 to the dimension of the problem to account for intercepts
c
c	augment x matrix with a vector of 1's
c
c       Compute r-square        
        do i=1,p
           xsum=0
           ysum=0
           do j=1,n
              xsum=xsum+xin(j,i)
              ysum=ysum+y(j,i)
           enddo
           xbar(i)=xsum/n
           ybar(i)=ysum/n
        enddo
        do i=1,p
           part1=0
           part2=0
           part3=0
           do j=1,n
              part1=part1+(xin(j,i)-xbar(i))*(y(j,i)-ybar(i))
              part2=part2+(xin(j,i)-xbar(i))**2
              part3=part3+(y(j,i)-ybar(i))**2
           enddo
           r(i)=part1**2/(part2*part3)
        enddo

	q1=q+1
	do 10 i=1,n
	do 11 j=1,q
		x(i,j+1)=xin(i,j)
11	continue
		x(i,1)=1.d0
10	continue
c
c get means of x,y
c
c	do 83 j=1,q+1
c	  xb(j)=0.d0
c83	continue
c	do 82 i=1,n
c	do 82 j=1,q+1
c	  xb(j)=xb(j)+x(i,j)
c82	continue
c	do 87 j=1,q+1
c		xb(j)=xb(j)/float(n)
c87	continue
c	write(*,*) ' means of x ',xb
c
c compute gammas
c
c first, get x'x
c
        do ir=1,q1
        do ic=1,ir
          xtx(ir,ic) = 0.0
          do il = 1,n
            xtx(ir,ic) = xtx(ir,ic) + x(il,ir)*x(il,ic)
          enddo
          xtx(ic,ir) = xtx(ir,ic)
        enddo
        enddo
c	call dmxtxf(n,q1,x,n,q1,xtx,q1)
c
c get inverse (x'x)-1
c
c        print *, 'x=', x
c        print *, 'xtx=', xtx

        do ir=1,q1
        do ic=1,q1
          if ( ir .eq. ic ) then
            eye(ir,ic) = 1.0
          else
            eye(ir,ic) = 0.0
          endif
        enddo
        enddo
        call dgesv( q1, q1, xtx, q1, ipiv, eye, q1, info )
        do ir=1,q1
        do ic=1,q1
          xtx(ir,ic) = eye(ir,ic)
        enddo
        enddo

c       call dlinrg(q1,xtx,q1,xtx,q1)
c
c (x'x)-1 x'
c
        do ir=1,q1
        do ic=1,n
          xix(ir,ic) = 0.0
          do il = 1,q1
            xix(ir,ic) = xix(ir,ic) + xtx(ir,il)*x(ic,il)
          enddo
        enddo
        enddo
c       call dmxytf(q1,q1,xtx,q1,n,q1,x,n,q1,n,xix,q1)
c
c (x'x)-1 x' y
c
        do ir=1,q1
        do ic=1,p
          gamma(ir,ic) = 0
          do il = 1,n
            gamma(ir,ic) = gamma(ir,ic) + xix(ir,il)*y(il,ic)
          enddo
        enddo
        enddo

c       call dmrrrr(q1,n,xix,q1,n,p,y,n,q1,p,gamma,q1)

      deallocate(g)
      deallocate(aprime)
      deallocate(x)
      deallocate(xbar)
      deallocate(ybar)
      deallocate(xtx)
      deallocate(xix)
      deallocate(vargam)
      deallocate(e)
      deallocate(proj)
      deallocate(yp)
      deallocate(s)
      deallocate(xtx2)
      deallocate(temp)
      deallocate(isubn)
      deallocate(ipiv)
      deallocate(eye)

      return
      end

Ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine data_summary
c     this subroutine is used to summarize the data statistics.

      use variable
      implicit none

      integer i,nd,k

      double precision, dimension(:), allocatable ::
     +  naivci_left,naivci_right,or_naiv,rrcci_left,rrcci_right,
     +  or_rrc,pval
      character(len=10), dimension(:), allocatable :: names
      character*255 filename

      integer outfile2
      nd=nex+ncov

      allocate(naivci_left(nd))
      allocate(naivci_right(nd))
      allocate(or_naiv(nd))
      allocate(rrcci_left(nd))
      allocate(rrcci_right(nd))
      allocate(or_rrc(nd))
      allocate(pval(nd))
      allocate(names(nd))

      filename=trim(pathname)//'/Output.'//trim(procid)//'.txt'
      open(unit=outfile1,file=filename,status='replace')
      filename=trim(pathname)//'/Output.'//trim(procid)//'.csv'
      open(unit=1016,file=filename,status='replace')

      write(outfile1,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'  
      write(outfile1,*) 
     * 'This analysis is for cumulatively average updated exposure'
      write(outfile1,*)
     * 'with classical additive measurement error model.'
      write(outfile1,*) ''
      write(outfile1,*) 'The increment for each exposure is:',exunit
      write(outfile1,*) ''  

      if ((nex.eq.1).and.(ncov.eq.0).and.(isres.eq.0)) then
         write(outfile1,*)
     *'There is only one exposure in the analysis and no covariate.'
         write(outfile1,*) 'The main exposure is measured with error.'
      else if ((nex.eq.1).and.(ncov.ge.1)) then
          write(outfile1,*)      'There is one exposure and one or more
     *covariates in the analysis.'
          write(outfile1,*)   'The exposure is measured with error, but
     *the covariates are error-free.'
      endif

      if(dabs(tau)<1.0d-9) then
        write(outfile1,*) 'It used compound symmetry covariance matrix'
        write(outfile1,*) 'with ICC=',icc,'.'
      else if (dabs(tau-1.0d0)<1.0d-9) then
        write(outfile1,*) 'It used AR(1) covariance matrix with'
        write(outfile1,*) 'ICC=',icc,'.'
      else
        write(outfile1,*) 'It used damped exponential covariance matrix'
        write(outfile1,*) 'with ICC=',icc, ',',
     *'damp factor THETA=', tau,'.'
      endif

      if ((isres.eq.1)) then
          write(outfile1,*)
     *"The residual analysis has been used in this analysis"
          write(outfile1,*) "for the dimension reduction."
      endif

      write(outfile1,*) ''

      write(outfile1,*) '# of subjects in the main study:', nmain
      write(outfile1,*) '# of person-year observations in the main',
     *     ' study:', mainobs   
      write(outfile1,*) '# of cases in the main study:',
     *     nint(sum(del_main))
      
      write(outfile1,*) '# of unique failure time:', nfa
      write(outfile1,*) '' 
      write(outfile1,*) 'The unique failure times are:'
      do k=1,int(dble(nfa)/10.0d0)
         write(outfile1,*) tf((1+(k-1)*10):(10*k))
      enddo
      write(outfile1,*) tf((1+(k-1)*10):nfa)
      write(outfile1,*) '' 

      write(outfile1,*) '# of subjects in the validation study:',nval
      write(outfile1,*) '# of person-year observations in the ',
     *     'validation study:', valobs
c      write(outfile1,*) '# of cases in the validation',
c     *     ' study:',nint(sum(del_val))
      write(outfile1,*) '# of risk sets in the validation study:', ngrp
      write(outfile1,*) ''

      write(outfile1,*) '----------------------------------------------'
      write(outfile1,*) 'Information about the measurement error model' 
      write(outfile1,*) 'fitting in the validation study:'
      write(outfile1,*) '----------------------------------------------'

      write(outfile1,*) ''

      write(outfile1,199) 'Risk_set', 'Age_low', 'Age_up',
     *'Person_year','Intercept', 'slope','   R^2' 

 199  format(7(A,2X))

      do k=1,ngrp
         write(outfile1,200) k, age_low(k), age_up(k), count4(k), 
     *acom(k,1,1),acom(k,2:(1+nex+ncov),1),rsquare(k,1)
      enddo

 200  format(I4,6X,2(F6.1,2X),I8,6X,3(F7.3,2X))

      write(outfile1,*) ''
      write(outfile1,*) '----------------------------------------------'
      write(outfile1,*) 'Results for uncorrected estimates (naive cox)' 
      write(outfile1,*) 'and corrected estimates (RRC):'
      write(outfile1,*) '----------------------------------------------'
      write(outfile1,*) ''
                  
      names(1:nex)=expname(1:nex)
      names((nex+1):nd)=covname(1:ncov)

      write(1016,*) 'variable,modeltype,para_es,stderr,pval,hr,',
     *'lowerci,upperci,groupnum,numrisksets'

      do i=1,nd
         call CI_95_alco(betaes_naiv(i), serror_naiv(i),or_naiv(i),
     *        naivci_left(i), naivci_right(i), pval(i))
      enddo

      write(outfile1,*) 'Uncorrected result:'
      write(outfile1,90) 'Variable','Para_es','S.E.','p_val','H.R.',
     *     '95% C.I.'
      do i=1,nd
         write(outfile1,100) names(i), betaes_naiv(i),serror_naiv(i), 
     *  pval(i),or_naiv(i),'[',naivci_left(i),',',naivci_right(i),']'
         write(1016,101) trim(names(i)), betaes_naiv(i),serror_naiv(i),
     *  pval(i),or_naiv(i),naivci_left(i),naivci_right(i),groupnum,ngrp
      enddo


      write(outfile1,*) '----------------------------------------------'
      do i=1,nd
         call CI_95_alco(betaes_rrc(i), serror_rrc(i), or_rrc(i), 
     *        rrcci_left(i),rrcci_right(i), pval(i))
      enddo

      write(outfile1,*) 'RRC result:'
      write(outfile1,90) 'Variable','Para_es','S.E.','p_val','H.R.',
     *     '95% C.I.'
      do i=1,nd
         write(outfile1,100) names(i), betaes_rrc(i),serror_rrc(i), 
     *  pval(i),or_rrc(i),'[',rrcci_left(i),',',rrcci_right(i),']'
         write(1016,101) trim(names(i)), betaes_rrc(i),serror_rrc(i),
     *  pval(i),or_rrc(i),rrcci_left(i),rrcci_right(i),groupnum,ngrp
      enddo 


      write(outfile1,*) ''
      write(outfile1,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

 90   format(1X,A,6X,A7,1X,A6,2X,A5,2X,A5,2X,A14)
 100  format(1X,A,4X,F6.3,2X,F6.3,2X,F5.3,2X,F5.3,2X, A, F6.2, A, F6.2,
     *     A)
 101  format(A,',uncorrected',6(',',F6.3),2(',',I4))
 102  format(A,',rrc',6(',',F6.3),2(',',I4))

      close(1016)

      close(outfile1)


cccc temporarily added for the bootstrap variance
c
c      outfile2=2011
c      open(unit=outfile2,file='Outboot.txt',access='append',
c     *     status='unknown')
c      write(outfile2,*) betaes_naiv(1),serror_naiv(1)
c      write(outfile2,*) betaes_rrc(1),serror_rrc(1)      
c      close(outfile2)

      deallocate(naivci_left)
      deallocate(naivci_right)
      deallocate(or_naiv)
      deallocate(rrcci_left)
      deallocate(rrcci_right)
      deallocate(or_rrc)
      deallocate(pval)
      deallocate(names)
      
      return
      end


ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine CI_95_alco(esti, serror, oddr, in_left, in_right, pval)
c     this subroutine is used to calculate the coverage probability
c     and the relative percentage bias.

      implicit none

c     variables on the entry
      real*8 esti, serror,oddr,in_left, in_right, pval 

c     dummy variables
      real*8 cr, scale, dchidf
      real*8 xp, xq, bound
      integer istat

      cr=1.960d0
      scale=1.0d0
      oddr=dexp(esti*scale)
      in_left=dexp((esti-cr*serror)*scale)
      in_right=dexp((esti+cr*serror)*scale)
c     pval=1.0d0-dchidf(esti**2/serror**2,1.0d0)
      call cdfchi(1,xp,pval,esti**2/serror**2,1.0d0,istat,bound)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

c     subroutine cal_cov(ldx,x1,x2,x1mean,x2mean,corr, nmiss)

c     use variable
c     implicit none

c     variables for IMSL subroutine corvc
c-------------------------------------------------------------
c     integer i
c     integer ido,nrow,nvar,ldx,ifrq,iwt,mopt,icopt,ldcov
c     integer ldincd,nobs,nmiss
c     parameter (nvar=2,ldcov=2,ldincd=1)
c     real*8 x(ldx,nvar),cov(ldcov,nvar),sumwt,xmean(nvar)
c     integer incd(nvar,nvar)

c     variables on the entry
c     real*8 x1(ldx), x2(ldx), x1mean, x2mean, corr

c     external dcorvc
c-------------------------------------------------------------

c     ido=0
c     nrow=ldx
c     ifrq=0
c     iwt=0
c     mopt=3

c     COV contains the correlation matrix, except
c     for the diagonal elements, which are the standard deviations.
c     icopt=3


c     do i=1,ldx
c        x(i,1)=x1(i)
c        x(i,2)=x2(i)
c     enddo

c     call dcorvc(ido,nrow,nvar,x,ldx,ifrq,iwt,mopt,icopt,xmean,cov,
c    *     ldcov,incd,ldincd,nobs,nmiss,sumwt)

c     x1mean=xmean(1)
c     x2mean=xmean(2)
c     corr=cov(1,2)

c     return
c     end
ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

