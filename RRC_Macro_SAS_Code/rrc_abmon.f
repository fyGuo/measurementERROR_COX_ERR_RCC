
ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssc
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
ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssc
c                                                                      c
c     This program is good for the measurement error correction with   c
c     12-month average exposure.                                       c
c                                                                      c
c
c     It works for multi-dimensional case.                             c
c                                                                      c
c     Last modified: March 5, 2013
c
c     It used dbcoah subroutine instead of ffsqp
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssc


      module variable

      implicit none

      real*8 pi
      parameter (pi=3.14159265358979323846)

      integer nmain,nval,nage_val,nage_main,nageobs,ncov,ncovm,nex
c     07/19/2013
c     Xiaomei added new variables for the new number of covariates.
      integer ncovvm, ncove, ncovem

      integer isres, mainobs, valobs, groupnum, istype
      integer ndist,nfa,nal,nfa2, case_cnt
      integer outfile1
      integer ngrp, month1, month
c     month1<month
      integer nmain_new,nmain_newcp, mainobscp

      integer, dimension(:), allocatable ::
     +  j_val, j_main, sub_id_val, sub_id_main, ilow, iup, fnum,
     +  ry_val, ry_main, cntcom, count1, count2, count3, count4,
     +  pd_main, mis_main, idmain_new, jmain_new

      integer, dimension(:,:), allocatable ::
     +  fsub, id_val, period, missing

      double precision, dimension(:), allocatable ::
     +  age_low, age_up, t_val, del_val, t_main, del_main,
     +  evcnt, tf, ageobs, r22, exunit,exunit2,
     +  betaes_naiv, serror_naiv,
     +  betaes_rrc, serror_rrc, tmain_new, delmain_new,
     +  xguess, xupbd, xlowbd

      double precision, dimension(:,:), allocatable ::
     +  age_main, bc_main, age_val, bc_val, rsquare,
     +  agemain_new, bcmain_new

      double precision, dimension(:,:,:), allocatable ::
     +  dtcalc_val, ffqcalc_val, zt_depval, zs_depval,
     +  covar_val, zs_depmain, zs_depmain_po, covar_main,
     +  zs_common, covar_common, zshat_depmain,
     +  covarhat_main, xxcom, yycom, acom, zshat_depmain_po,
     +  zsdepmain_new, covmain_new, zshat_depmainnew,
     +  zs_depmain_res
   

c These are for filenames

      character*255 pathname
      character*20 procid

c     variable for exposure and confounder names
      character(len=10), dimension(:), allocatable ::
     +  expname, covname

      end

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
c      read(1013,*) ncov
      read(1013,*) isres
      read(1013,*) istype
      read(1013,*) mainobs
      read(1013,*) valobs
      read(1013,*) groupnum
      read(1013,*) month
      read(1013,*) month1
      read(1013,*) ncovm
      read(1013,*) ncove
      read(1013,*) ncovem
      close(1013)

c     the new number of covariates=ncov(the total number of covariates
c     in the validation study)+ncovm(the total number of covariates in
c     the mains tudy)
      if (isres==1) then
         ncov=ncove
         ncovm=ncovem+ncovm
         ncovvm=ncove+ncovm
      else
         ncov=ncove+ncovm
         ncovm=ncovem+ncovm
         ncovvm=ncove+ncovm
      endif

      allocate(j_val(nval))
      allocate(j_main(nmain))
      allocate(sub_id_val(nval))
      allocate(sub_id_main(nmain))
      allocate(ilow(nmain))
      allocate(iup(nmain))
      allocate(fnum(nmain))
      allocate(ry_val(nageobs))
      allocate(ry_main(nageobs))
      allocate(cntcom(nageobs))
      allocate(count1(nageobs))
      allocate(count2(nageobs))
      allocate(count3(nageobs))
      allocate(count4(nageobs))

      allocate(fsub(nageobs,2))
      allocate(id_val(valobs,nageobs))

      allocate(age_low(nageobs))
      allocate(age_up(nageobs))
      allocate(t_val(nval))

      allocate(del_val(nval))
      allocate(t_main(nmain))
      allocate(del_main(nmain))
      allocate(evcnt(nmain))
      allocate(tf(nageobs))
      allocate(ageobs(nmain))
      allocate(r22(nageobs))
      allocate(exunit(nex))
      allocate(exunit2(nex))
      allocate(betaes_naiv(nex+ncovm))
      allocate(serror_naiv(nex+ncovm))
      allocate(betaes_rrc(nex+ncovm))
      allocate(serror_rrc(nex+ncovm))

      allocate(age_main(nmain,nage_main))
      allocate(bc_main(nmain,nage_main))
      allocate(age_val(nval,nage_val))
      allocate(bc_val(nval,nage_val))
      allocate(rsquare(nageobs,nex))

      allocate(dtcalc_val(nval,nage_val,nex))
      allocate(ffqcalc_val(nval,nage_val,nex))
      allocate(zt_depval(nval,nage_val,nex))
      allocate(zs_depval(nval,nage_val,nex+ncov))
      allocate(covar_val(nval,nage_val,0:ncov))
      allocate(zs_depmain(nmain,nage_main,nex))
      allocate(zs_depmain_res(nmain,nage_main,nex))
      allocate(zs_depmain_po(nmain,nage_main,nex))
      allocate(covar_main(nmain,nage_main,0:ncovvm))
      allocate(zs_common(nmain,nage_main,nex))
      allocate(covar_common(nmain,nage_main,0:ncovvm))
      allocate(zshat_depmain(nmain,nage_main,nex))
      allocate(zshat_depmain_po(nmain,nage_main,nex))
      allocate(covarhat_main(nmain,nage_main,0:ncovvm))
      allocate(xxcom(valobs,nageobs,nex+ncov))
      allocate(yycom(valobs,nageobs,nex))
      allocate(acom(nageobs,1+nex+ncov,nex))

      allocate(period(nmain,nage_main))
      allocate(missing(nmain,nage_main))
      allocate(pd_main(nmain))
      allocate(mis_main(nmain))
      allocate(idmain_new(nmain))
      allocate(jmain_new(nmain))
      allocate(agemain_new(nmain,nage_main))
      allocate(bcmain_new(nmain,nage_main))
      allocate(tmain_new(nmain))
      allocate(delmain_new(nmain))
      allocate(zsdepmain_new(nmain,nage_main,nex))
      allocate(covmain_new(nmain,nage_main,0:ncovvm))
      allocate(zshat_depmainnew(nmain,nage_main,nex))

      allocate(expname(nex))
      allocate(covname(ncovvm))

      allocate(xguess(nex+ncovm))
      allocate(xupbd(nex+ncovm))
      allocate(xlowbd(nex+ncovm))

      return
      end

      subroutine deallocvar
      use variable
      implicit none

      deallocate(j_val)
      deallocate(j_main)
      deallocate(sub_id_val)
      deallocate(sub_id_main)
      deallocate(ilow)
      deallocate(iup)
      deallocate(fnum)
      deallocate(ry_val)
      deallocate(ry_main)
      deallocate(cntcom)
      deallocate(count1)
      deallocate(count2)
      deallocate(count3)
      deallocate(count4)

      deallocate(fsub)
      deallocate(id_val)

      deallocate(age_low)
      deallocate(age_up)
      deallocate(t_val)
      deallocate(del_val)
      deallocate(t_main)
      deallocate(del_main)
      deallocate(evcnt)
      deallocate(tf)
      deallocate(ageobs)
      deallocate(r22)
      deallocate(exunit)
      deallocate(exunit2)
      deallocate(betaes_naiv)
      deallocate(serror_naiv)
      deallocate(betaes_rrc)
      deallocate(serror_rrc)

      deallocate(age_main)
      deallocate(bc_main)
      deallocate(age_val)
      deallocate(bc_val)
      deallocate(rsquare)

      deallocate(dtcalc_val)
      deallocate(ffqcalc_val)
      deallocate(zt_depval)
      deallocate(zs_depval)
      deallocate(covar_val)
      deallocate(zs_depmain)
      deallocate(zs_depmain_res)
      deallocate(zs_depmain_po)
      deallocate(covar_main)
      deallocate(zs_common)
      deallocate(covar_common)
      deallocate(zshat_depmain)
      deallocate(zshat_depmain_po)
      deallocate(covarhat_main)
      deallocate(xxcom)
      deallocate(yycom)
      deallocate(acom)

      deallocate(period)
      deallocate(missing)
      deallocate(pd_main)
      deallocate(mis_main)
      deallocate(idmain_new)
      deallocate(jmain_new)
      deallocate(agemain_new)
      deallocate(bcmain_new)
      deallocate(tmain_new)
      deallocate(delmain_new)
      deallocate(zsdepmain_new)
      deallocate(covmain_new)
      deallocate(zshat_depmainnew)

      deallocate(expname)
      deallocate(covname)

      deallocate(xguess)
      deallocate(xupbd)
      deallocate(xlowbd)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      program hpfs_data

      use variable
      implicit none

      integer i,j,k,k1,b, jtemp

      integer, dimension(:), allocatable ::
     +  sub_id_valtmp, sub_id_maintmp, xcnt,
     +  xhatcnt, cnt2, cnmiss, xnmiss, casecnt_age,
     +  periodtmp, missingtmp

      double precision, dimension(:), allocatable ::
     +  age_valtmp, bc_valtmp, age_maintmp,
     +  bc_maintmp, xtmp, xhattmp, xave, bxave,
     +  c_corr, x_corr, a0temp, a1temp

      double precision, dimension(:,:), allocatable ::
     +  dtcalc_valtmp, ffqcalc_valtmp, ffqcalc_maintmp,
     +  covar_maintmp, covar_valtmp

      real*8 a, age_temp(43), age_start, age_end

c  get path and proc ID to reconstruct filenames
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

c  allocate arrays in module variable
      call allocvar

      allocate(sub_id_valtmp(valobs))
      allocate(dtcalc_valtmp(valobs,nex))
      allocate(ffqcalc_valtmp(valobs,nex))
      allocate(age_valtmp(valobs))
      allocate(bc_valtmp(valobs))
      allocate(sub_id_maintmp(mainobs))
      allocate(ffqcalc_maintmp(mainobs,2*nex))


cc    add in covariates
      allocate(covar_maintmp(mainobs,0:ncovvm))
      allocate(covar_valtmp(valobs,0:ncov))
      allocate(age_maintmp(mainobs))
      allocate(bc_maintmp(mainobs))
      allocate(xtmp(nageobs))
      allocate(xhattmp(nageobs))
      allocate(xave(nageobs))
      allocate(bxave(nageobs))
      allocate(c_corr(nageobs))
      allocate(x_corr(nageobs))
      allocate(a0temp(nageobs))
      allocate(a1temp(nageobs))

      allocate(xcnt(nageobs))
      allocate(xhatcnt(nageobs))
      allocate(cnt2(nageobs))
      allocate(cnmiss(nageobs))
      allocate(xnmiss(nageobs))
      allocate(casecnt_age(nageobs))

      allocate(periodtmp(mainobs))
      allocate(missingtmp(mainobs))


c     outfile1 is the unit number for the output file
      outfile1=1011

      count2=-1
      count3=-1
      count4=-1


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
      do i=1,ncovm
         read(71,*) covname(i)
      enddo
      close(71)


c     Read in the initial value for each variable.
      filename=trim(pathname)//'/xinitial.'//trim(procid)//'.dat'
      open(81,file=filename,action='read')
      do i=1,nex+ncovm
         read(81,*) xguess(i)
      enddo
      close(81)

c     Read in the upper bound for each variable.
      filename=trim(pathname)//'/xupbound.'//trim(procid)//'.dat'
      open(81,file=filename,action='read')
      do i=1,nex+ncovm
         read(81,*) xupbd(i)
      enddo
      close(81)

c     Read in the lower bound for each variable.
      filename=trim(pathname)//'/xlowbound.'//trim(procid)//'.dat'
      open(81,file=filename,action='read')
      do i=1,nex+ncovm
         read(81,*) xlowbd(i)
      enddo
      close(81)


c     Read in the main study data set.

      filename=trim(pathname)//'/main.'//trim(procid)//'.dat'
      open(101,file=filename, action='read')
      do i=1,mainobs
       if (isres==0) then
         if (ncovvm.ge.1) then
            read(101,*)
     *      sub_id_maintmp(i),age_maintmp(i),bc_maintmp(i),
     *      ffqcalc_maintmp(i,1:nex),covar_maintmp(i,1:ncovvm),
     *      periodtmp(i), missingtmp(i)
         else
            read(101,*)
     *      sub_id_maintmp(i),age_maintmp(i),bc_maintmp(i),
     *      ffqcalc_maintmp(i,1:nex),periodtmp(i), missingtmp(i)
         endif
       else
         if (ncovvm.ge.1) then
            read(101,*)
     *      sub_id_maintmp(i),age_maintmp(i),bc_maintmp(i),
     *      ffqcalc_maintmp(i,1:(2*nex)),covar_maintmp(i,1:ncovvm),
     *      periodtmp(i), missingtmp(i)
         else
            read(101,*)
     *      sub_id_maintmp(i),age_maintmp(i),bc_maintmp(i),
     *      ffqcalc_maintmp(i,1:(2*nex)),periodtmp(i), missingtmp(i)
         endif
       endif
      enddo

      close(101)

c     scaling by the user input increment
      
      do j=1,nex
         ffqcalc_maintmp(:,j)=ffqcalc_maintmp(:,j)/exunit(j)
      enddo
      if (isres==1) then
         exunit2=exunit
         do j=1,nex
            ffqcalc_maintmp(:,(nex+j))=ffqcalc_maintmp(:,(nex+j))
     *           /exunit2(j)
         enddo
      endif	


c      way to check the correctness of reading of the main study data
c       print *,'i=',i
c       print *,sub_id_maintmp(i-1),age_maintmp(i-1),bc_maintmp(i-1),
c     *     ffqcalc_maintmp(i-1,1:nex),covar_maintmp(i-1,0:ncov)



c     transformation of main study
      k=1
      j=1
      sub_id_main(k)=sub_id_maintmp(k)
      zs_depmain(k,j,1:nex)=ffqcalc_maintmp(k,1:nex)
      if (isres==1) then
         zs_depmain_res(k,j,1:nex)=ffqcalc_maintmp(k,(nex+1):2*nex)
      endif      
      age_main(k,j)=age_maintmp(k)
      bc_main(k,j)=bc_maintmp(k)
      period(k,j)=periodtmp(k)
      missing(k,j)=missingtmp(k)


      if (ncovvm.ge.1) then
         covar_main(k,j,1:ncovvm)=covar_maintmp(k,1:ncovvm)
      endif
      j_main(k)=j

      do i=2,mainobs
         if(sub_id_maintmp(i).eq.sub_id_maintmp(i-1)) then
            j=j+1
            sub_id_main(k)=sub_id_maintmp(i)
            zs_depmain(k,j,1:nex)=ffqcalc_maintmp(i,1:nex)
            if (isres==1) then
            zs_depmain_res(k,j,1:nex)=ffqcalc_maintmp(i,(nex+1):2*nex)
            endif
            age_main(k,j)=age_maintmp(i)
            bc_main(k,j)=bc_maintmp(i)
            period(k,j)=periodtmp(i)
            missing(k,j)=missingtmp(i)
            if (ncovvm.ge.1) then
             covar_main(k,j,1:ncovvm)=covar_maintmp(i,1:ncovvm)
            endif
            j_main(k)=j
         else
            k=k+1
            j=1
            j_main(k)=j
            sub_id_main(k)=sub_id_maintmp(i)
            zs_depmain(k,j,1:nex)=ffqcalc_maintmp(i,1:nex)
            if (isres==1) then
            zs_depmain_res(k,j,1:nex)=ffqcalc_maintmp(i,(nex+1):2*nex)
            endif
            age_main(k,j)=age_maintmp(i)
            bc_main(k,j)=bc_maintmp(i)
            period(k,j)=periodtmp(i)
            missing(k,j)=missingtmp(i)
            if (ncovvm.ge.1) then
             covar_main(k,j,1:ncovvm)=covar_maintmp(i,1:ncovvm)
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
         pd_main(i)=period(i,j_main(i))
         mis_main(i)=missing(i,j_main(i))
      enddo

      call sort_main_simple(nmain,nage_main)

      print *, 'The following results are for the naive cox method: '


      zs_depmain_po=zs_depmain

      call cal_abmon_ave(month1,month,zs_depmain_po,zs_depmain,j_main)


c     create a new main study with only the subjects and exposure after
c     the first 12-month.
c     new variable name in the new main data set:
c     idmain_new, zsdepmain_new, agemain_new, bcmain_new, period_new
c     missing_new
c     use the new variable to store the old copy
c     use the old variable to store the new data for the convenience

      idmain_new=sub_id_main
      jmain_new=j_main
      zsdepmain_new=zs_depmain
      agemain_new=age_main
      bcmain_new=bc_main
c      period_new=period
c      missing_new=missing
      covmain_new=covar_main
      tmain_new=t_main
      delmain_new=del_main

      nmain_new=0
      do i=1,nmain
         if ((pd_main(i)>0).and.(mis_main(i).eq.0)) then
            nmain_new=nmain_new+1
            sub_id_main(nmain_new)=idmain_new(i)
            j_main(nmain_new)=0
            do j=1,jmain_new(i)
               if ((period(i,j)>0).and.(missing(i,j).eq.0)) then
                  j_main(nmain_new)=j_main(nmain_new)+1
                  if (j_main(nmain_new).eq.0) then
                     print *, 'nmain_new=', nmain_new
                     print *, 'j_main=', j_main(nmain_new)
                  endif
                  zs_depmain(nmain_new,j_main(nmain_new),:)
     *                 =zsdepmain_new(i,j,:)
                  age_main(nmain_new,j_main(nmain_new))
     *                 =agemain_new(i,j)
                  bc_main(nmain_new,j_main(nmain_new))
     *                 =bcmain_new(i,j)
c                  period(nmain_new,j_main(nmain_new))
c     *                 =period_new(i,j)
c                  missing(nmain_new,j_main(nmain_new))
c     *                 =missing_new(i,j)

                  if (ncovvm.ge.1) then
                   covar_main(nmain_new,j_main(nmain_new),1:ncovvm)
     *                    =covmain_new(i,j,1:ncovvm)
                  endif
               endif
            enddo
         endif
      enddo

      print *, "nmain_new=", nmain_new
      print *, "mainobs=", sum(j_main(1:nmain_new))

      do i=1,nmain_new
         if (j_main(i).eq.0) then
            print *, 'i=', i
            print *, 'id=',sub_id_main(i)
            print *, 'j_main=', j_main(i)
         endif
         t_main(i)=age_main(i,j_main(i))
         del_main(i)=bc_main(i,j_main(i))
      enddo

c     some statistics for data summary
      nmain_newcp=nmain_new
      mainobscp=sum(j_main(1:nmain_new))

      print *, "nmain_new=", nmain_new, nmain_newcp
      print *, "mainobscp=", mainobscp

      print *, "# of cases=",nint(sum(del_main(1:nmain_newcp)))
c      print *, "2 # of cases=",nint(sum(del_main))


      filename=trim(pathname)//'/xsurrogate.'//trim(procid)//'.dat'
      open(999,file=filename, action='write')
      do i=1,nmain_new
         do j=1,j_main(i)
            write(999,*) sub_id_main(i), age_main(i,j), 
     *           zs_depmain(i,j,1:nex)
         enddo
      enddo
      close(999)


      call sort_main(nmain_new,nage_main)

      call naivecoxmethod

c      stop

      zsdepmain_new=zs_depmain_po


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

c      stop

      call rrccoxmethod

c      stop

c     summarizing the result

      call data_summary

      deallocate(sub_id_valtmp)
      deallocate(dtcalc_valtmp)
      deallocate(ffqcalc_valtmp)
      deallocate(age_valtmp)
      deallocate(bc_valtmp)
      deallocate(sub_id_maintmp)
      deallocate(ffqcalc_maintmp)

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
      deallocate(a0temp)
      deallocate(a1temp)

      deallocate(xcnt)
      deallocate(xhatcnt)
      deallocate(cnt2)
      deallocate(cnmiss)
      deallocate(xnmiss)
      deallocate(casecnt_age)
      deallocate(periodtmp)
      deallocate(missingtmp)


      call deallocvar

c     corresponding to the program
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c     this subroutine is to calculate the b~a month average exposure,
c     a>b
c     for example, a=12, b=1 ---> 12-month moving average exposure
c     for example, a=12, b=6 ---> 6~12 month average exposure

      subroutine cal_abmon_ave(b,a,zp,zs,jlim)

      use variable
      implicit none

      integer i,j,k,j1,a,b,ctp,jlim(nmain)
      real*8 zp(nmain,nage_main,nex), zs(nmain,nage_main,nex), s

      do k=1,nex
         do i=1,nmain
            do j=1,jlim(i)
               s=0.0d0
               if (j.lt.b) then
                  zs(i,j,k)=zp(i,j,k)
               else if (j.lt.a) then
                  ctp=0
                  do j1=1,j-(b-1)
                     s=s+zp(i,j1,k)
                     if (missing(i,j1).eq.0) then
                        ctp=ctp+1
                     endif
                  enddo
                  if (ctp.ne.0) then
                     zs(i,j,k)=s/dble(ctp)
                  else
                     zs(i,j,k)=0.0d0
                  endif
               else
                  ctp=0
                  do j1=(j-(a-1)),(j-(b-1))
                     s=s+zp(i,j1,k)
                     if (missing(i,j1).eq.0) then
                        ctp=ctp+1
                     endif
                  enddo
                  if (ctp.ne.0) then
                     zs(i,j,k)=s/dble(ctp)
                  else
                     zs(i,j,k)=0.0d0
                  endif
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
      real*8 zp(nmain,nage_main,nex), zs(nmain,nage_main,nex), s

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
      real*8 zp(nmain,nage_main,nex), zs(nmain,nage_main,nex),s

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

      subroutine sort_main_simple(nmax,ntm)

      use variable
      implicit none

c     on the entry variable
      integer nmax, ntm

c     dummy variable

      integer, dimension(:), allocatable :: wksp, indx
      integer, dimension(:,:), allocatable :: wksp5
      double precision, dimension(:), allocatable :: wksp1, rsk, tt_main
      double precision, dimension(:,:), allocatable :: wksp2
      double precision, dimension(:,:,:), allocatable :: wksp3,
     +  wksp4

      real*8 told,rsknum,rm
      integer ix,i,j

      external quicksort

      allocate(wksp(nmax))
      allocate(wksp5(nmax,ntm))
      allocate(wksp1(nmax))
      allocate(wksp2(nmax,ntm))
      allocate(wksp3(nmax,ntm,0:ncovvm))
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
      call quicksort(nmax,t_main,indx)
      do i = 1, nmax
        t_main(i) = tt_main(indx(i))
      enddo

c  ... (debugging) read sort order

c       open(1101,file='sort.simple.orig.out',action='read')
c       do idebug = 1,nmax
c         read(1101,*) indx(idebug), t_main(idebug)
c       enddo
c 404   format(i8,e20.10)
c       close(1101)

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
         wksp(i) = pd_main(i)
      enddo
      do i = 1, nmax
         pd_main(i) = wksp(indx(i))
      enddo
      do i = 1, nmax
         wksp(i) = mis_main(i)
      enddo
      do i = 1, nmax
         mis_main(i) = wksp(indx(i))
      enddo

      do i = 1, nmax
         wksp5(i,:) = period(i,:)
      enddo
      do i = 1, nmax
         period(i,:) = wksp5(indx(i),:)
      enddo
      do i = 1, nmax
         wksp5(i,:) = missing(i,:)
      enddo
      do i = 1, nmax
         missing(i,:) = wksp5(indx(i),:)
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
         wksp4(i,:,:) = zs_depmain_res(i,:,:)
      enddo
      do i = 1, nmax      
         zs_depmain_res(i,:,:) = wksp4(indx(i),:,:)    
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
      if (ncovvm.ge.1) then
         do i = 1, nmax
            wksp3(i,:,1:ncovvm)= covar_main(i,:,1:ncovvm)
         enddo
         do i = 1, nmax
            covar_main(i,:,1:ncovvm) = wksp3(indx(i),:,1:ncovvm)
         enddo
      endif

      deallocate(wksp)
      deallocate(wksp5)
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


      subroutine sort_main(nmax,ntm)

      use variable
      implicit none

c     on the entry variable
      integer nmax, ntm

c     dummy variable

      real*8 told,rsknum,rm
      integer ix,i,j

      double precision, dimension(:), allocatable :: rsk

      allocate(rsk(nmax))

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

      print *, '# of unique failure time =', nfa
      print *, 'The unique failure times are',tf(1:nfa)

      deallocate(rsk)

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
c  ... (debugging) read sort order

c       open(1101,file='sort.val.orig.out',action='read')
c       do idebug = 1,nmax
c         read(1101,*) indx(idebug), t_val(idebug)
c       enddo
c 404   format(i8,e20.10)
c       close(1101)

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

      integer nd, ifail, i
      double precision, dimension(:), allocatable ::
     +  betaes, st_error

      double precision, dimension(:,:), allocatable :: hesinv

      nd=nex+ncovm
      allocate(betaes(nd))
      allocate(st_error(nd))
      allocate(hesinv(nd,nd))

      zs_common=zs_depmain
      covar_common=covar_main

      if (nd.ge.2) then
         call cox_bcoah(nd,betaes,st_error,ifail,hesinv)
      else if (nd.eq.1) then
         call cox_zbren(nd,betaes,st_error,ifail)
      else
         print *, 'The dimension of data is zero!'
      endif

      betaes_naiv=betaes
      serror_naiv=st_error

      if (ifail.eq.0) then
         print *,'naive','beta=', betaes
         print *,'se=', st_error

         print *, "This is for the cov(1,2)"
         
         print *, "Cov(1,2)=", hesinv(1,2)

         print *, "*************************"

         print *, "The covariance matrix"

         do i=1,nd
            print *,hesinv(i,1:nd)
            print *," "
         enddo

         print *, "*************************"

      else
         print *,'The method does not converge.'
      endif

      deallocate(betaes)
      deallocate(st_error)
      deallocate(hesinv)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine cox_bcoah(nd,betaes,st_error,ifail,hesinv)

      use variable
      implicit none

c     variables on the entry
      integer nd, ifail
      real*8 betaes(nd), st_error(nd), hesinv(nd,nd)

c     dummy variables
      real*8 fvalue,fscale

      integer, dimension(:), allocatable ::
     +  ipiv, iwrkf
      double precision, dimension(:), allocatable ::
     +  xlb, xub, xscale, x, eval, wrkf
      double precision, dimension(:,:), allocatable ::
     +  hes_fval, r


      integer iparam(7), rparam(7),ibtype
      integer isev, iwk(nd),i, ier
      integer iercd, n1rty

      integer info, ir, ic, nw, nwsize, iwsize
      real*8 work(1000)

      external rrccox_likelihood, rrcgradient_cox
      external rrchessian_cox, ffsqp

c ----------------------------------------------------------------------------------
      double precision fval(1), gval(1)
      integer iprint
      external obj,constr,gradobj,grcnfd
c ----------------------------------------------------------------------------------

      iwsize=6*nd+8+7+30
      nwsize=4*(nd**2)+5*nd+3*nd+26*(nd+1)+45+100

c      allocate(xguess(nd))
      allocate(xlb(nd))
      allocate(xub(nd))
      allocate(xscale(nd))
      allocate(x(nd))
      allocate(hes_fval(nd,nd))
c      allocate(hesinv(nd,nd))
      allocate(eval(nd))
      allocate(r(nd,nd))
      allocate(ipiv(nd))
      allocate(iwrkf(iwsize))
      allocate(wrkf(nwsize))



c      xguess= 0.02d0
      x=xguess
      print *, 'initial guess=', xguess

      ibtype=0
c      xlb(1)=-10.0d0
c      xub(1)=10.0d0
c      xlb(2)=-10.0d0  
c      xub(2)=10.0d0
c      xlb(3:nd)=-10.0d0
c      xub(3:nd)=10.0d0

c      xlb(1)=-5.0d0
c      xub(1)=5.0d0                 
c      xlb(2)=-5.0d0            
c      xub(2)=5.0d0
c      xlb(3:nd)=-5.0d0
c      xub(3:nd)=5.0d0

      xlb=xlowbd
      xub=xupbd

      xscale=1.0d0
      fscale=1.0d0
      iparam(1)=0

c ----------------------------------------------------------------------------------
      iprint = 3
      write(*,*) 'Calling ffsqp'
      call ffsqp(nd,1,0,0,0,0,100,iprint,500,info,1d+10,1d-4,0,0,
     *  xlb, xub, x, fval, gval, iwrkf,iwsize, wrkf,nwsize,
     *  obj,constr,gradobj,grcnfd)
c ----------------------------------------------------------------------------------

c      call dbcoah(rrccox_likelihood,rrcgradient_cox, rrchessian_cox
c     *     ,nd,xguess,ibtype,xlb,xub,xscale,fscale,iparam,rparam,x,
c     *     fvalue)
c
c      ier = iercd()
c      isev = n1rty(1)
c      ifail = 1
c      print *, 'isev=', isev

      	
c      if (isev .eq. 0) then
      if (info .eq. 0 ) then

         betaes=x
c         print *,'bcoah estimates=', betaes

         call rrchessian_cox(nd,betaes,hes_fval,nd)
c         print *, 'hessian_ana=', hes_fval


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

c      deallocate(xguess)
      deallocate(xlb)
      deallocate(xub)
      deallocate(xscale)
      deallocate(x)
      deallocate(hes_fval)
c      deallocate(hesinv)
      deallocate(eval)
      deallocate(r)
      deallocate(ipiv)
      deallocate(iwrkf)
      deallocate(wrkf)

      return
      end

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
      real*8 betaes(nd), st_error(nd)

c     dummy variables
      integer maxfn, isev, n1rty
      real*8 dx,abmx,za(1),zb(1),errabs,errrel, derv
      real*8 coxscore

      double precision, dimension(:,:), allocatable ::
     +  hess

      real*8 rrr
      integer istat
      external coxscore, fzero

      allocate(hess(nd,nd))

      dx=0.10d0
      abmx=xguess(1)
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
      isev = 0
      ifail = 1
      st_error=0.0
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

c     dummy variables
      integer nd
      real*8 bval(1), g_fval(1)
      nd=1

      call rrcgradient_cox(nd,bval,g_fval)
      coxscore=g_fval(1)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine rrc

      use variable
      implicit none

      integer pp,qq,ir,il
      integer i,j,k,k1,i1,j1,j2
      integer ifail,groupsize,newgroup

      real*8 beta_rrc, st_rrc, st_rrcs, st_rrcs1, scount
      double precision, dimension(:,:), allocatable ::
     +  gamma, temp1, temp2
      character*255 filename

      pp=nex
      qq=nex+ncov

      allocate(gamma(qq+1,pp))
      allocate(temp1(pp,1))
      allocate(temp2(qq+1,1))

c     Define the risk process ry_main and ry_val for each observation time
      ry_main=0
      i=1
      do j=1,nfa
         do while (t_main(i).lt.tf(j))
            if (i.lt.nmain_new) then
               i=i+1
            else
               ry_main(j)=nmain_new+1
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
cc     Still use this version
c
      groupsize=0
      ngrp=0
      newgroup=1
      do i=1,nfa
         if (newgroup==1) then
            age_low(ngrp+1)=tf(i)
            groupsize=0
c             print *,'age_low',ngrp+1,age_low(ngrp+1)
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
c             print *,'age_up',ngrp,age_up(ngrp)
         else if ((groupsize<=groupnum-1) .and. (i/=nfa)) then
            newgroup=0
         else if (ngrp==0) then
            newgroup=0
            ngrp=ngrp+1
            age_up(ngrp)=tf(nfa)
c            print *,'age_up',ngrp,age_up(ngrp),tf(i)
         else
            newgroup=0
            age_up(ngrp)=tf(nfa)
         endif
      enddo
c
c      print *, 'test1'
c      print *, "ngrp=", ngrp

c     if the censoring age of a subject is not younger than age_up()
c     then this subject is in the risk set (age_low, age_up).

      id_val=0
c     special dealing with the first entry
      k=1
      i1=0
      do i=1,nval
         if (t_val(i).ge.age_low(k)) then
            do j1=1,j_val(i)
               if ((age_val(i,j1).ge.age_low(k))
     *              .and.(age_val(i,j1).le.age_up(k))
     *              .and.(zt_depval(i,j1,1).ne.-1.0d0)
     *              .and.(zs_depval(i,j1,1).ne.-1.0d0)
     *              ) then
                  i1=i1+1
                  xxcom(i1,k,1:(nex+ncov))=
     *                 zs_depval(i,j1,1:(nex+ncov))
                  yycom(i1,k,1:nex)=zt_depval(i,j1,1:nex)
                  id_val(i1,k)=sub_id_val(i)
               endif
            enddo
         endif
      enddo
      count4(k)=i1


      do k=2, ngrp-1
         i1=0
         do i=1,nval
            if (t_val(i).ge.age_low(k)) then
               do j1=1,j_val(i)
                  if ((age_val(i,j1).ge.age_low(k))
     *                 .and.(age_val(i,j1).le.age_up(k))
     *                 .and.(zt_depval(i,j1,1).ne.-1.0d0)
     *                 .and.(zs_depval(i,j1,1).ne.-1.0d0)
     *                 ) then
                     i1=i1+1
                     xxcom(i1,k,1:(nex+ncov))=
     *                    zs_depval(i,j1,1:(nex+ncov))
                     yycom(i1,k,1:nex)=zt_depval(i,j1,1:nex)
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
         if (t_val(i).ge.age_low(k)) then
            do j1=1,j_val(i)
               if ((age_val(i,j1).ge.age_low(k))
     *              .and.(age_val(i,j1).le.age_up(k))
     *              .and.(zt_depval(i,j1,1).ne.-1.0d0)
     *              .and.(zs_depval(i,j1,1).ne.-1.0d0)
     *              ) then
                  i1=i1+1
                  xxcom(i1,k,1:(nex+ncov))=
     *                 zs_depval(i,j1,1:(nex+ncov))
                  yycom(i1,k,1:nex)=zt_depval(i,j1,1:nex)
                  id_val(i1,k)=sub_id_val(i)
               endif
            enddo
         endif
      enddo
      count4(k)=i1


      do j=1,ngrp
         if (count4(j).ge.2) then
            call mvlm2(pp,qq,count4(j),yycom(1:count4(j),j,:),
     *           xxcom(1:count4(j),j,:),gamma,rsquare(j,:))

            acom(j,:,:)=gamma
c            print *, 'j,acom(j,k1)=', j, acom(j,:,:)
            print *, 'j,age_low_up', j, age_low(j), age_up(j)
            print *, 'j,acom(j,k1)=', j,count4(j),acom(j,:,:)
         endif
      enddo

      print *, 'test2'


c     count how many subjects used in the risk set regression
      scount=0.0d0
      do k=1,ngrp
         scount=scount+count4(k)
      enddo

      print *, 'sum of count=', scount

c     construct the c_hat value from the regression coefficients
c     then calculate x_hat from c_hat
      zshat_depmain_po=0.0d0
      do i=1,nmain
         do j1=1,jmain_new(i)
           if (missing(i,j1).eq.0) then
            if (agemain_new(i,j1).lt.tf(1)) then
               k=1
               if (isres==0) then
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zsdepmain_new(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)
                  else
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zsdepmain_new(i,j1,1:nex)/)
                  endif
               else
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zs_depmain_res(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)      
                  else
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zs_depmain_res(i,j1,1:nex)/)
                  endif
               endif

               do ir = 1,pp
                 temp1(ir,1) = 0.0
                 do il = 1,qq+1
                   temp1(ir,1) = temp1(ir,1)+acom(k,il,ir)*temp2(il,1)
                 enddo
               enddo
c              call dmxtyf(qq+1,pp,acom(k,:,:),qq+1,qq+1,1,
c    *              temp2,qq+1,pp,1,temp1,pp)
               zshat_depmain_po(i,j1,:)=temp1(:,1)

            else if (agemain_new(i,j1).gt.tf(nfa)) then
               k=ngrp
	       if (isres==0) then
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zsdepmain_new(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)
                  else
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zsdepmain_new(i,j1,1:nex)/)
                  endif
               else
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zs_depmain_res(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)
                  else               
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zs_depmain_res(i,j1,1:nex)/)              
                  endif
               endif
               do ir = 1,pp
                 temp1(ir,1) = 0.0
                 do il = 1,qq+1
                   temp1(ir,1) = temp1(ir,1) + acom(k,il,ir)*temp2(il,1)
                 enddo
               enddo
c              call dmxtyf(qq+1,pp,acom(k,:,:),qq+1,qq+1,1,
c    *              temp2,qq+1,pp,1,temp1,pp)
               zshat_depmain_po(i,j1,:)=temp1(:,1)

            else
               j=1
               do while(agemain_new(i,j1)>tf(j))
                  j=j+1
               enddo
               if(agemain_new(i,j1)<tf(j)) j=j-1

               do k=1,ngrp
                  if ((tf(j).ge.age_low(k)).and.(tf(j).le.age_up(k)))
     *              then
                     if (isres==0) then
                        if (ncov.ge.1) then
                           temp2(:,1)=(/1.0d0,zsdepmain_new(i,j1,1:nex),
     *                          covmain_new(i,j1,1:ncov)/)
                        else
                           temp2(1:(1+nex),1)=
     *                          (/1.0d0,zsdepmain_new(i,j1,1:nex)/)
                        endif
		     else
                       if (ncov.ge.1) then
                          temp2(:,1)=(/1.0d0,zs_depmain_res(i,j1,1:nex),
     *                          covmain_new(i,j1,1:ncov)/)
                       else
                          temp2(1:(1+nex),1)=                  
     *                         (/1.0d0,zs_depmain_res(i,j1,1:nex)/)
                       endif
                     endif
                     do ir = 1,pp
                       temp1(ir,1) = 0.0
                       do il = 1,qq+1
                         temp1(ir,1) = temp1(ir,1)
     *                       + acom(k,il,ir)*temp2(il,1)
                       enddo
                     enddo
c                    call dmxtyf(qq+1,pp,acom(k,:,:),qq+1,qq+1,1,
c    *                    temp2,qq+1,pp,1,temp1,pp)
                     zshat_depmain_po(i,j1,:)=temp1(:,1)
                  endif
               enddo
            endif
           endif
         enddo
      enddo

c      print *, 'test3'

      call cal_abmon_ave(month1,month,zshat_depmain_po,zshat_depmain,
     *     jmain_new)

c      print *, 'test4'


c  ..... debugging ..... name=100017 index=64540 ndata=72

c     keep the wanted data
      zshat_depmainnew=zshat_depmain

      nmain_new=0
      do i=1,nmain
         if ((pd_main(i)>0).and.(mis_main(i).eq.0)) then
            nmain_new=nmain_new+1
            j_main(nmain_new)=0
            do j=1,jmain_new(i)
               if ((period(i,j)>0).and.(missing(i,j).eq.0)) then
                  j_main(nmain_new)=j_main(nmain_new)+1
                  zshat_depmain(nmain_new,j_main(nmain_new),:)
     *                 =zshat_depmainnew(i,j,:)

               endif
            enddo

c            if ( idmain_new(i) .eq. 100017 ) then
c              write(*,*) 'squeezing in rrc '
c              write(*,*) 'i =',i,' nmain_new=',nmain_new
c              write(*,*) 'j_main = ', j_main(nmain_new)
c              do j = 1,jmain_new(i)
c                write(*,'(e16.6)',advance='no')
c     *             zshat_depmain(nmain_new,j,1)
c              enddo
c              write(*,*)
c              do j = 1,jmain_new(i)
c                write(*,'(e16.6)',advance='no')
c     *             zshat_depmain_po(nmain_new,j,1)
c              enddo
c            endif

         endif
      enddo

      print *, "nmain_new=", nmain_new
      print *, "nmainobs=", sum(j_main(1:nmain_new))


      filename=trim(pathname)//'/xhat.'//trim(procid)//'.dat'
      open(1000,file=filename, action='write')
      do i=1,nmain_new
         do j=1,j_main(i)
            write(1000,*) sub_id_main(i), age_main(i,j), 
     *           zshat_depmain(i,j,1:nex)
         enddo
      enddo
      close(1000)
      


      deallocate(gamma)
      deallocate(temp1)
      deallocate(temp2)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      
      subroutine rrccoxmethod
      
      use variable
      implicit none
         
      integer nd, ifail, i

      double precision, dimension(:), allocatable ::
     +  betaes, st_error, st_errors, st_errors1

      double precision, dimension(:,:), allocatable :: hesinv

      nd=nex+ncovm
      allocate(betaes(nd)) 
      allocate(st_error(nd))
      allocate(st_errors(nd))
      allocate(st_errors1(nd))
      allocate(hesinv(nd,nd))

      zs_common=zshat_depmain
      covar_common=covar_main

      if (nd.ge.2) then
         call cox_bcoah(nd,betaes,st_error,ifail,hesinv)
      else if (nd.eq.1) then
         call cox_zbren(nd,betaes,st_error,ifail)
      else
         print *, 'The dimension of data is zero!'
      endif
      betaes_rrc=betaes

      print *, 'test5'
      print *, 'betaes_rrc=', betaes_rrc


      if (ifail.eq.0) then
         call rrc_sderror(betaes, st_error, st_errors, st_errors1,
     +        hesinv)

         serror_rrc=st_errors
c         print *, 'All variance=', serror_rrc
         print *, 'rrc','beta=', betaes_rrc
         print *, 'se=', serror_rrc

         print *, "This is for the cov(1,2)"
         
         print *, "Cov(1,2)=", hesinv(1,2)

         print *, "*************************"

         print *, "The covariance matrix"

         do i=1,nd
            print *,hesinv(i,1:nd)
            print *," "
         enddo

         print *, "*************************"

      else     
         print *,'The method does not converge.'
      endif

      deallocate(betaes) 
      deallocate(st_error)
      deallocate(st_errors)
      deallocate(st_errors1)
      deallocate(hesinv)
      
      return
      end


ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrc_sderror(bval,st_error, st_errors, st_errors1,
     +     covmatrix)

c     this subroutine is used to calculate the variance estimator of
c     rrc estimates. 

      use variable
      implicit none

c     variables on the entry
      integer nd, ndm
      real*8 bval(nex+ncovm),st_error(nex+ncovm),st_errors(nex+ncovm),
     +  st_errors1(nex+ncovm), covmatrix(nex+ncovm,nex+ncovm)

c     dummy variables
      integer i,j, npa, ir, ic, il

      double precision, dimension(:,:), allocatable ::
     +  fchat, fbhat, fbinv, fbsand, ftmp, paramatrix,
     +  paratemp1, g_score, st_erpsi,
     +  g_fval, tempmatrix2, r

      double precision, dimension(:), allocatable ::
     +  st_erbe1, st_erbe1s
      integer info

      nd=nex+ncov
      ndm=nex+ncovm
      npa=(1+nex+ncov)*nex

      allocate(fchat(ndm,ndm))
      allocate(fbhat(ndm,ndm))
      allocate(fbinv(ndm,ndm))
      allocate(fbsand(ndm,ndm))
      allocate(ftmp(ndm,ndm))
      allocate(paramatrix(ndm,ndm)) 
      allocate(paratemp1(ndm,ndm))
      allocate(st_erbe1(ndm))
      allocate(st_erbe1s(ndm))
      allocate(g_score(npa*ngrp,npa*ngrp))
*      allocate(u2matrix(npa*ngrp,npa*ngrp))
*      allocate(tempmatrix(npa*ngrp,npa*ngrp))
      allocate(st_erpsi(npa*ngrp,npa*ngrp))
      allocate(g_fval(npa*ngrp,ndm))
      allocate(tempmatrix2(ndm,npa*ngrp))
      allocate(r(ndm,ndm))


      call score2rrc(bval,fchat)
      call rrchessian_cox(ndm,bval,fbhat,ndm)

c     inverse hessian matrix part 
      print *, "test6"
      do ir=1,ndm
      do ic=ir,ndm
        r(ir,ic) = fbhat(ir,ic)
      enddo
      enddo

      call dpotrf('U', ndm, r, ndm, info)
      call dpotri('U', ndm, r, ndm, info)

      do ir=1,ndm
      do ic=ir,ndm
        fbinv(ir,ic) = r(ir,ic)
        fbinv(ic,ir) = fbinv(ir,ic)
      enddo
      enddo
c     call dlinds(nd,fbhat,nd,fbinv,nd) 

      print *, 'fbhat=', fbhat
      Do i=1,ndm
         st_erbe1(i)=dsqrt(fbinv(i,i))
      End Do

c     pure sandwich variance estimator part

      do ir=1,ndm
      do ic=1,ndm
        ftmp(ir,ic) = 0.0;
        do il=1,ndm
          ftmp(ir,ic) = ftmp(ir,ic) + fbinv(ir,il)*fchat(il,ic)
        enddo
      enddo
      enddo
c     call dmrrrr(nd,nd,fbinv,nd,nd,nd,fchat,nd,nd,nd,ftmp,nd)

      do ir=1,ndm
      do ic=1,ndm
        fbsand(ir,ic) = 0.0;
        do il=1,ndm
          fbsand(ir,ic) = fbsand(ir,ic) + ftmp(ir,il)*fbinv(il,ic)
        enddo
      enddo
      enddo
c     call dmrrrr(nd,nd,ftmp,nd,nd,nd,fbinv,nd,nd,nd,fbsand,nd)
      Do i=1,ndm 
         st_erbe1s(i)=dsqrt(fbsand(i,i)) 
      End Do 

c     parameters' variability part
      call parascore_grad(ngrp, st_erpsi)

      print *, 'st_erpsi=',  st_erpsi      

      call grrcscore_para(ngrp,bval,g_fval)
      print *, 'g_fval=', g_fval

      
      do ir = 1,ndm
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

      do ir=1,ndm
      do ic=1,ndm
        paramatrix(ir,ic) = 0.0;
        do il=1,npa*ngrp
          paramatrix(ir,ic) =
     *        paramatrix(ir,ic) + tempmatrix2(ir,il)*g_fval(il,ic)
        enddo
      enddo
      enddo
c     call dmrrrr(nd,npa*ngrp,tempmatrix2,nd,npa*ngrp,nd,g_fval,
c    *     npa*ngrp,nd,nd,paramatrix,nd)
      do ir=1,ndm
      do ic=1,ndm
        paratemp1(ir,ic) = 0.0;
        do il=1,ndm
          paratemp1(ir,ic) =
     *        paratemp1(ir,ic) + fbinv(ir,il)*paramatrix(il,ic)
        enddo
      enddo
      enddo
c     call dmrrrr(nd,nd,fbinv,nd,nd,nd,paramatrix,nd,nd,nd,
c    *paratemp1,nd)
      do ir=1,ndm
      do ic=1,ndm
        paramatrix(ir,ic) = 0.0;
        do il=1,ndm
          paramatrix(ir,ic) =
     *        paramatrix(ir,ic) + paratemp1(ir,il)*fbinv(il,ic)
        enddo
      enddo
      enddo
c     call dmrrrr(nd,nd,paratemp1,nd,nd,nd,fbinv,nd,nd,nd,
c    *paramatrix,nd)

c     print *, 'paramatrix=', paramatrix


c     Since the "paramatrix" contributes very big variation, I just 
c     output the variation from the main study "fbsand".

c      covmatrix=fbsand+paramatrix
      covmatrix=fbsand
      
c     three kinds of variance estimators
c     st_errors1: pure sandwich variance estimator
c     st_errors:  sandwich variance estimator + variability of parameters
c              :  the most accurate one.
c     st_error:   non-sandwich variance estimator + variability of 
c             :   parameters. not so mathematically correctly but acceptable
c             :   in this situation.

      do i=1,ndm
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
*      deallocate(u2matrix)
*      deallocate(tempmatrix)
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
      real*8 bval(nex+ncovm),fval(nex+ncovm,nex+ncovm)

c     dummy variables
      real*8 didl
      parameter (didl=1.0d-8)
      integer j, jtemp, i, k, j1, indk, ir, ic, il

      double precision, dimension(:), allocatable ::
     +  s0, s0a
      double precision, dimension(:,:), allocatable ::
     +  rrf1, s1, ss2,ftmp, fvaltmp, temp1, temp2
      double precision, dimension(:,:,:), allocatable ::
     +  rrf2

      real*8 s0b
      real*8 temp

      nd=nex+ncovm

      print *, "IN score2rrc"
      print *, "nd=", nd
      print *, "nex=", nex
      print *, "ncovm=", ncovm

      allocate(rrf1(nmain_new,nageobs))
      allocate(rrf2(nmain_new,nageobs,nd))
      allocate(s0(nageobs))
      allocate(s1(nageobs,nex+ncovm))
      allocate(s0a(nageobs))
      allocate(ss2(nmain_new,nex+ncovm))
      allocate(ftmp(nd,1))
      allocate(fvaltmp(nd,nd))
      allocate(temp1(nd,1))
      allocate(temp2(nd,1))

      rrf1=0.0d0; rrf2=0.0d0
      temp1(:,1)=bval
      do j=1,nfa
         do i=ry_main(j),nmain_new
            do j1=1,j_main(i)
               if ((age_main(i,j1).eq.tf(j)).and.
     *              (zshat_depmain(i,j1,1).ne.-100.0d0)) then
                  if (ncovm.ge.1) then
                     temp2(:,1)=(/zshat_depmain(i,j1,1:nex),
     *                    covar_main(i,j1,(ncove+1):(ncove+ncovm))/)
                  else
                     temp2(1:nex,1)=zshat_depmain(i,j1,1:nex)
                  endif
                  temp = 0.0
                  do il = 1,nd
                    temp = temp + temp1(il,1)*temp2(il,1)
                  enddo
c                 call dmxtyf(nd,1,temp1,nd,nd,1,temp2,nd,1,1,temp,1)

                  rrf1(i,j)=dexp(temp)
                  rrf2(i,j,1:nex)=rrf1(i,j)*zshat_depmain(i,j1,1:nex)
                  if (ncovm.ge.1) then
                     rrf2(i,j,(nex+1):nd)=rrf1(i,j)*
     *                    covar_main(i,j1,(ncove+1):(ncove+ncovm))
                  endif
               endif
            enddo
         enddo
      enddo


      s0=0.0d0
      s1=0.0d0
      do j=1,nfa
         do i=ry_main(j), nmain_new
            if (rrf1(i,j).ne.0.0d0) then
               s0(j)=s0(j)+rrf1(i,j)
               do k=1,nd
                  s1(j,k)=s1(j,k)+rrf2(i,j,k)
               enddo
            endif
         enddo
      enddo

 
      ss2=0.0d0
      Do i = 1,nmain_new
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


      Do i = 1,nmain_new
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
      end Do

      print *, "ftmp=", ftmp
     
      deallocate(rrf1)
      deallocate(rrf2)
      deallocate(s0)
      deallocate(s1)
      deallocate(s0a)
      deallocate(ss2)

      if (allocated(ftmp)) deallocate(ftmp)

      deallocate(fvaltmp)
      deallocate(temp1)
      deallocate(temp2)

      print *, "End score2rrc"

      return
      end


ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c     This is the subroutine of the numerical derivative for calculating
c     the gradient of rrc score w.r.t the regression parameters
c     (\alpha_0,\alpha_1,\alpha_2).
ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine grrcscore_para(ng,bval,g_fval)

      use variable
      implicit none

      integer i,j,k
      integer nd,ng
      real*8 g_fval(ng*(1+nex+ncov)*nex,nex+ncovm),bval(nex+ncovm)
      real*8 da

      double precision, dimension(:), allocatable ::
     +  fval, fval1
      double precision, dimension(:,:,:), allocatable ::
     +  a

      nd=nex+ncov
      allocate(a(ng,1+nd,nex))
      allocate(fval(nex+ncovm))
      allocate(fval1(nex+ncovm))

      da=0.0001d0
      a=acom(1:ng,:,:)
      call scorerrc_para(ng,bval,a,fval)

      do k=1,ng
         do i=1,(1+nd)
            do j=1,nex
               a(k,i,j)=a(k,i,j)+da
               call scorerrc_para(ng,bval,a,fval1)
               g_fval((k-1)*(1+nd)+(i-1)*nex+j,:)=(fval1-fval)/da
               a(k,i,j)=a(k,i,j)-da
            enddo
         enddo
      enddo

  601 format(4e20.8)

      deallocate(a)
      deallocate(fval)
      deallocate(fval1)

      return
      end
ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


      subroutine scorerrc_para(ng,bval,a,fval)

      use variable
      implicit none

      integer nd,ng, pp, qq, ndm
      integer i, j, k, j1, ir, il
      real*8 fval(nex+ncovm),bval(nex+ncovm), a(ng,1+nex+ncov,nex)

      double precision, dimension(:,:), allocatable ::
     +  temp1, temp2
      double precision, dimension(:,:,:), allocatable ::
     +  zshat_tmp, zshat_tmp_po, zshat_tmpnew

      nd=nex+ncov
      ndm=nex+ncovm
      pp=nex
      qq=nex+ncov

      allocate(temp1(pp,1))
      allocate(temp2(qq+1,1))
      allocate(zshat_tmp_po(nmain,nage_main,nex))
      allocate(zshat_tmp(nmain,nage_main,nex))
      allocate(zshat_tmpnew(nmain,nage_main,nex))

      zshat_tmp_po = 0.0d0
      zshat_tmp = 0.0d0
      zshat_tmpnew = 0.0d0

c     construct the c_hat value from the regression coefficients
c     then calculate x_hat from c_hat

      zshat_tmp_po=0.0d0
      do i=1,nmain
         do j1=1,jmain_new(i)
           if (missing(i,j1).eq.0) then
            if (agemain_new(i,j1).lt.tf(1)) then
               k=1
               if (isres==0) then
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zsdepmain_new(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)
                  else
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zsdepmain_new(i,j1,1:nex)/)
                  endif
               else
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zs_depmain_res(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)      
                  else
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zs_depmain_res(i,j1,1:nex)/)
                  endif
               endif   
               do ir=1,pp
                 temp1(ir,1) = 0.0
                 do il=1,qq+1
                   temp1(ir,1) = temp1(ir,1) +
     *               a(k,il,ir)*temp2(il,1)
                 enddo
               enddo
               zshat_tmp_po(i,j1,:)=temp1(:,1)

            else if (agemain_new(i,j1).gt.tf(nfa)) then
               k=ngrp
               if (isres==0) then
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zsdepmain_new(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)
                  else
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zsdepmain_new(i,j1,1:nex)/)
                  endif
               else
                  if (ncov.ge.1) then
                     temp2(:,1)=(/1.0d0,zs_depmain_res(i,j1,1:nex),
     *                    covmain_new(i,j1,1:ncov)/)
                  else               
                     temp2(1:(1+nex),1)=
     *                    (/1.0d0,zs_depmain_res(i,j1,1:nex)/)              
                  endif
               endif

               do ir=1,pp
                 temp1(ir,1) = 0.0
                 do il=1,qq+1
                   temp1(ir,1) = temp1(ir,1) +
     *               a(k,il,ir)*temp2(il,1)
                 enddo
               enddo
               zshat_tmp_po(i,j1,:)=temp1(:,1)

            else
               j=1
               do while(agemain_new(i,j1)>tf(j))
                  j=j+1
               enddo
               if(agemain_new(i,j1)<tf(j)) j=j-1

               do k=1,ngrp
                  if ((tf(j).ge.age_low(k)).and.(tf(j).le.age_up(k)))
     *              then
                     if (isres==0) then
                        if (ncov.ge.1) then
                           temp2(:,1)=(/1.0d0,zsdepmain_new(i,j1,1:nex),
     *                          covmain_new(i,j1,1:ncov)/)
                        else
                           temp2(1:(1+nex),1)=
     *                          (/1.0d0,zsdepmain_new(i,j1,1:nex)/)
                        endif
                     else
                        if (ncov.ge.1) then
                          temp2(:,1)=(/1.0d0,zs_depmain_res(i,j1,1:nex),
     *                          covmain_new(i,j1,1:ncov)/)
                       else
                          temp2(1:(1+nex),1)=                  
     *                         (/1.0d0,zs_depmain_res(i,j1,1:nex)/)
                       endif
                     endif

                     do ir=1,pp
                       temp1(ir,1) = 0.0
                       do il=1,qq+1
                         temp1(ir,1) = temp1(ir,1) +
     *                     a(k,il,ir)*temp2(il,1)
                       enddo
                     enddo
                     zshat_tmp_po(i,j1,:)=temp1(:,1)
                  endif
               enddo
            endif
           endif
         enddo
      enddo

      call cal_abmon_ave(month1,month,zshat_tmp_po,zshat_tmp,jmain_new)

c     keep the wanted data

      nmain_new=0
      do i=1,nmain
         if ((pd_main(i)>0).and.(mis_main(i).eq.0)) then
            nmain_new=nmain_new+1
            j_main(nmain_new)=0
            do j=1,jmain_new(i)
               if ((period(i,j)>0).and.(missing(i,j).eq.0)) then
                  j_main(nmain_new)=j_main(nmain_new)+1
                  zshat_tmpnew(nmain_new,j_main(nmain_new),:)
     *                 =zshat_tmp(i,j,:)

               endif
            enddo
         endif
      enddo

      zs_common=zshat_tmpnew
      covar_common=covar_main

      call rrcgradient_cox(ndm,bval,fval)

      deallocate(temp1)
      deallocate(temp2)
      deallocate(zshat_tmp_po)
      deallocate(zshat_tmp)
      deallocate(zshat_tmpnew)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine parascore_grad(nd,st_erpsi)
c     this subroutine is used to calculate the variance of
c     \psi=(\alpha_0,\alpha_1, \alpha_2)
c     this way consists of stacking of matrices, which is time consuming.

      use variable
      implicit none

c     variables on the entry
      integer nd, npa, ir, ic, il

      real*8 st_erpsi((1+nex+ncov)*nex*nd,(1+nex+ncov)*nex*nd)

c     dummy variables
      integer i,j,k,k1,k2,i1
      integer pp,qq
      integer g_p_, ldg_aa, ldg_ua

      double precision, dimension(:), allocatable ::
     +  aa, ua
      double precision, dimension(:,:), allocatable ::
     +  g_aa, g_ua, g_uas, temp1, temp2, g_temp,
     +  uwork1, uwork2, gamma, rtmp, r
      double precision, dimension(:,:,:), allocatable ::
     +  g_tempj, uall
      double precision, dimension(:,:,:,:), allocatable ::
     + uwork, xtmp, ytmp, u2matrix


      integer info

      external dpotrf dpotri

      npa=(1+nex+ncov)*nex
      pp=nex
      qq=nex+ncov

      allocate(aa(npa))
      allocate(ua(npa))
      allocate(g_aa(npa,npa))
      allocate(g_ua(npa,npa))
      allocate(g_uas(npa,npa))
      allocate(temp1(nex+ncov,1))
      allocate(temp2(nex,1))
      allocate(g_temp(npa,npa))
      allocate(g_tempj(nd,npa,npa))
      allocate(uall(nd,valobs,npa))
      allocate(u2matrix(nd,nd,npa,npa))
      allocate(uwork1(1,npa))
      allocate(uwork2(1,npa))
      allocate(uwork(nd,nd,npa,npa))
      allocate(gamma(qq+1,pp))
      allocate(rtmp(nageobs,nex))
      allocate(xtmp(valobs,valobs,nageobs,nex+ncov))
      allocate(ytmp(valobs,valobs,nageobs,nex))
      allocate(r(npa,npa))

      g_p_=npa
      ldg_aa=g_p_
      ldg_ua=g_p_

      g_aa=0.0d0

      do j=1,npa
         g_aa(j,j)=1.0d0
      enddo

      uall=0.0d0

c     use cross-validation idea to get uall,
c     (Donna's, leave-one-out residual)

      do j=1,nd
         do i=1,count4(j)
            i1=0
            do k=1,count4(j)
               if (k.ne.i) then
                  i1=i1+1
                  xtmp(i1,i,j,:)=xxcom(k,j,:)
                  ytmp(i1,i,j,:)=yycom(k,j,:)
               endif
            enddo
            call mvlm2(pp,qq,i1,ytmp(1:i1,i,j,:),
     &           xtmp(1:i1,i,j,:),gamma,rtmp(j,:))

c           print *, 'gamma=', gamma

            do k1=1,1+nex+ncov
               do k2=1,nex
                  aa((k1-1)*nex+k2)=gamma(k1,k2)
               enddo
            enddo

            temp1(:,1)=xxcom(i,j,:)
            temp2(:,1)=yycom(i,j,:)

            call g_parascore(g_p_, temp1, temp2,
     *           aa, g_aa, ldg_aa, ua, g_uas,ldg_ua)

c            print *, 'ua=', ua
c            write(*,202) ((g_uas(ir,ic), ir=1,npa), ic=1,npa)
c 202        format(9f10.4)

            uall(j,i,:)=ua
         enddo
      enddo


      Do j = 1,nd
         g_ua=0.0d0
         do i=1,1+nex+ncov
            do k=1,nex
               aa((i-1)*nex+k)=acom(j,i,k)
            enddo
         enddo

         Do i= 1,count4(j)
            temp1(:,1)=xxcom(i,j,:)
            temp2(:,1)=yycom(i,j,:)
            call g_parascore(g_p_, temp1, temp2,
     *           aa, g_aa, ldg_aa, ua, g_uas,ldg_ua)

            g_ua=g_ua+g_uas

         Enddo


         g_temp=-g_ua

c     Compute the inverse of a real symmetric positive definite matrix

c         print *, "test7"

c         print *, "g_temp=", g_temp

        do ir=1,npa
        do ic=ir,npa
          r(ir,ic) = g_temp(ir,ic)
        enddo
        enddo

        call dpotrf('U', npa, r, npa, info)
        call dpotri('U', npa, r, npa, info)

        do ir=1,npa
        do ic=ir,npa
          g_tempj(j,ir,ic) = r(ir,ic)
          g_tempj(j,ic,ir) = g_tempj(j,ir,ic)
        enddo
        enddo
c       call dlinds(npa,g_temp,npa,g_tempj(j,:,:),npa)

      Enddo

c     This part calculates \sum U*U^t.
      Do j=1,nd
         Do i=j,nd
            u2matrix(j,i,:,:)=0.0d0
            Do k1=1,max0(count4(i),count4(j))
               Do k2=1,max0(count4(i),count4(j))
                  if ((id_val(k1,i).ne.0).and.(id_val(k2,j).ne.0).and.
     *                 (id_val(k1,i).eq.id_val(k2,j))) then

c                     print *, "k1=", k1
c                     print *, "k2=", k2
c                     print *, "i=",i
c                     print *, "j=",j
c                     print *, "id_val(k1,i)=",id_val(k1,i)
c                     print *, "id_val(k2,j)=",id_val(k2,j)

                     uwork1(1,:)=uall(j,k2,:)
                     uwork2(1,:)=uall(i,k1,:)
                     
c                     print *, "uwork1(1,:)=",uwork1(1,:)
c                     print *, "uwork2(1,:)=",uwork2(1,:)

                     do ir=1,npa
                     do ic=1,npa
                        uwork(j,i,ir,ic) = uwork1(1,ir)*uwork2(1,ic)

c                        print *,"ir=", ir
c                        print *, "ic=", ic
c                        print *,"uwork(j,i,ir,ic)=",uwork(j,i,ir,ic)
                     enddo
                     enddo
c                    call dmxtyf(1,npa,uwork1,1,1,npa,uwork2,1,npa,npa,
c    *                    uwork(j,i,:,:) ,npa)

c                     print *, "uwork(j,i,:,:)=",uwork(j,i,:,:)
c                     print *, "u2matrix(j,i,:,:)=",u2matrix(j,i,:,:)

                     u2matrix(j,i,:,:)=u2matrix(j,i,:,:)+uwork(j,i,:,:)
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
              uwork(i,j,ir,ic) = 0.0d0
              do il=1,npa
                uwork(i,j,ir,ic) = uwork(i,j,ir,ic) +
     *              g_tempj(i,ir,il)*u2matrix(i,j,il,ic)

c                print *, "uwork=",uwork(i,j,ir,ic)
c                print *, "g_tempj=", g_tempj(i,ir,il)
                
              enddo
            enddo
            enddo
c           call dmrrrr(npa,npa,g_tempj(i,:,:),npa,npa,npa,
c    *u2matrix(i,j,:,:),npa,npa,npa,uwork(i,j,:,:),npa)

            do ir=1,npa
            do ic=1,npa
              u2matrix(i,j,ir,ic) = 0.0d0
              do il=1,npa
                u2matrix(i,j,ir,ic) = u2matrix(i,j,ir,ic) +
     *              uwork(i,j,ir,il)*g_tempj(j,il,ic)

c                print *, "u2matrix(i,j,ir,ic)=",u2matrix(i,j,ir,ic)

              enddo
            enddo
            enddo
c           call dmrrrr(npa,npa,uwork(i,j,:,:),npa,npa,npa,
c    *g_tempj(j,:,:),npa,npa,npa,u2matrix(i,j,:,:),npa)
         enddo

         do j=1,i-1
            do ir=1,npa
            do ic=1,npa
              u2matrix(i,j,ir,ic) = u2matrix(j,i,ic,ir)
            enddo
            enddo
c           call dtrnrr(npa,npa,u2matrix(j,i,:,:),npa,npa,npa,
c    *u2matrix(i,j,:,:),npa)
         enddo
      enddo

      do i=1,nd
         do j=1,nd
            st_erpsi(1+npa*(i-1):npa*i,1+npa*(j-1):npa*j)=
     *u2matrix(i,j,:,:)
c            print *, "test_xm"
c            print *, "i=", i
c            print *, "j=", j
c            print *, "U2matrix=", u2matrix
c            print *, "st_erpsi=", st_erpsi
         enddo
      enddo

c      print *, "st_erpsi=", st_erpsi

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
      deallocate(gamma)
      deallocate(rtmp)
      deallocate(xtmp)
      deallocate(ytmp)
      deallocate(r)

c      print *, "END"
c      print *, "st_erpsi=", st_erpsi

      return
      end


ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine g_parascore(g_p_, x, y, aa, g_aa, ldg_aa, ua, g_ua,
     * ldg_ua)
C
      use variable
      implicit none

      integer g_p_, ldg_aa, ldg_ua
      real*8 x(nex+ncov, 1), y(nex, 1), aa((1 + nex+ncov) * nex),
     +  ua((1+nex+ncov)*nex), g_ua(ldg_ua,(1+nex+ncov)*nex),
     +  g_aa(ldg_aa, (1 + nex+ncov) * nex)

      integer i, j, nd, g_i_a
      double precision, dimension(:,:), allocatable ::
     +  a, atran, ua0, ua1, ua2, uatot, temp1, temp2, temp3,
     +  temp4, temp5

      double precision, dimension(:,:,:), allocatable, save ::
     +  g_a, g_ua2, g_temp4, g_uatot, g_ua0, g_ua1, g_atran,
     +  g_temp1, g_temp2, g_temp3

      integer g_i_, g_pmax_

      external g_matrixsubtract
      external g_matrixprod
      external g_matrixtranspose

      integer g_ehfid
      integer ir, ic, ip
      data g_ehfid /0/

      nd = nex + ncov
      allocate(a(1+nd,nex))
      allocate(atran(nex,1+nd))
      allocate(ua0(nex, 1))
      allocate(ua1(nex, nex))
      allocate(ua2(nex, 0:ncov))
      allocate(uatot(nex, 1 + nex+ncov))
      allocate(temp1(1+nex+ncov,1))
      allocate(temp2(nex, 1))
      allocate(temp3(nex, 1))
      allocate(temp4(1,nex))
      allocate(temp5(0:ncov, 1))

      g_pmax_ = (1+nd)*nex
      if ( .not. allocated(g_a) )  then
        allocate(g_a(g_pmax_, 1 + nd, nex))
        allocate(g_ua2(g_pmax_, nex, 0:ncov))
        allocate(g_temp4(g_pmax_, 1, nex))
        allocate(g_uatot(g_pmax_, nex, 1 + nd))
        allocate(g_ua0(g_pmax_, nex, 1))
        allocate(g_ua1(g_pmax_, nex, nex))
        allocate(g_atran(g_pmax_, nex,1 + nd))
        allocate(g_temp1(g_pmax_, 1 + nd, 1))
        allocate(g_temp2(g_pmax_, nex, 1))
        allocate(g_temp3(g_pmax_, nex,1))
      endif


        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do i = 1, 1 + nd
          do j = 1, nex
            do g_i_ = 1, g_p_
              g_a(g_i_, i, j) = g_aa(g_i_, (i - 1) * nex + j)
            enddo
            a(i, j) = aa((i - 1) * nex + j)
C--------
          enddo
        enddo
C
        temp1(1, 1) = 1.0d0
        do i = 1, nd
          temp1(1 + i, 1) = x(i, 1)
        enddo

c -------- DEBUG ------------------
        g_temp1 = 0d0
c -------- DEBUG ------------------
C
        call g_matrixtranspose(g_p_, 1 + nd, nex, a, g_a, g_pmax_, atran
     *, g_atran, g_pmax_)

        call g_matrixprod(g_p_, nex, 1 + nd, atran, g_atran, g_pmax_, 1
     *+ nd, 1, temp1, g_temp1, g_pmax_, nex, 1, temp2, g_temp2, g_pmax_)
        call g_matrixsubtract(g_p_, nex, 1, y, temp2, g_temp2, g_pmax_,
     *ua0, g_ua0, g_pmax_)

C
        do i = 1, nex
          temp3(i, 1) = x(i, 1)
        enddo
c -------- DEBUG ------------------
        g_temp3 = 0d0
c -------- DEBUG ------------------

        call g_matrixtranspose(g_p_, nex, 1, ua0, g_ua0, g_pmax_, temp4,
     * g_temp4, g_pmax_)
        call g_matrixprod(g_p_, nex, 1, temp3, g_temp3, g_pmax_, 1, nex,
     * temp4, g_temp4, g_pmax_, nex, nex, ua1, g_ua1, g_pmax_)
C
        if (ncov .ge. 1) then
          do i = 1, ncov
            temp5(i, 1) = x(nex + i, 1)
          enddo
C
          do i = 1, ncov
            do j = 1, nex
              do g_i_ = 1, g_p_
                g_ua2(g_i_, j, i) = temp5(i, 1) * g_temp4(g_i_, 1, j)
              enddo
              ua2(j, i) = temp5(i, 1) * temp4(1, j)
C--------
            enddo
          enddo
        endif
C
        do i = 1, nex
          do g_i_ = 1, g_p_
            g_uatot(g_i_, i, 1) = g_ua0(g_i_, i, 1)
          enddo
          uatot(i, 1) = ua0(i, 1)
C--------
          do j = 1, nex
            do g_i_ = 1, g_p_
              g_uatot(g_i_, i, 1 + j) = g_ua1(g_i_, j, i)
            enddo
            uatot(i, 1 + j) = ua1(j, i)
C--------
          enddo
          if (ncov .ge. 1) then
            do j = 1, ncov
              do g_i_ = 1, g_p_
                g_uatot(g_i_, i, 1 + nex + j) = g_ua2(g_i_, i, j)
              enddo
              uatot(i, 1 + nex + j) = ua2(i, j)
C--------
            enddo
          endif
        enddo
C
        do i = 1, 1 + nd
          do j = 1, nex
            do g_i_ = 1, g_p_
              g_ua(g_i_, (i - 1) * nex + j) = g_uatot(g_i_, j, i)
            enddo
            ua((i - 1) * nex + j) = uatot(j, i)
C--------
          enddo
        enddo
C
      deallocate(a)
      deallocate(atran)
      deallocate(ua0)
      deallocate(ua1)
      deallocate(ua2)
      deallocate(uatot)
      deallocate(temp1)
      deallocate(temp2)
      deallocate(temp3)
      deallocate(temp4)
      deallocate(temp5)

      return
      end
C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_matrixprod(g_p_, nar, nac, a, g_a, ldg_a, nbr, nbc, b
     *, g_b, ldg_b, ncr, ncc, c, g_c, ldg_c)
C
C     used to compute the real matrix(or vector) product c=a*b
C
      use variable
      implicit none
C
C     on the entry variables
      integer g_p_, nar, nac, ldg_a, ldg_b, ldg_c
      integer nbr, nbc, ncr, ncc

      real*8 a(nar, nac), b(nbr, nbc), c(ncr, ncc)
      double precision g_c(ldg_c, ncr, ncc), g_a(ldg_a, nar, nac),
     *  g_b(ldg_b, nbr, nbc)

      integer i, j, k, g_pmax_, g_i_

      integer g_ehfid
      data g_ehfid /0/

      g_pmax_ = (1+nex+ncov)*nex
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        if ((nac .ne. nbr) .or. (nar .ne. ncr) .or. (nbc .ne. ncc)) then
          print *, 'Dimension does not match'
          stop
        endif
C
        do i = 1, nar
          do j = 1, ncc
            do g_i_ = 1, g_p_
              g_c(g_i_, i, j) = 0.0d0
            enddo
            c(i, j) = 0.0d0
C--------
            do k = 1, nac
              do g_i_ = 1, g_p_
                g_c(g_i_, i, j) = a(i, k) * g_b(g_i_, k, j) + b(k, j) *
     *g_a(g_i_, i, k) + g_c(g_i_, i, j)
              enddo
              c(i, j) = c(i, j) + a(i, k) * b(k, j)
C--------
            enddo
          enddo
        enddo
C
        return
      end
C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_matrixsubtract(g_p_, nr, nc, a, b, g_b, ldg_b, c,
     *  g_c, ldg_c)

C     used to compute the real matrix(or vector) subtraction c=a-b
C
      use variable
      implicit none
C
C     on the entry variables
        integer nr, nc
        real*8 a(nr, nc), b(nr, nc), c(nr, nc)
        integer g_p_, ldg_c, ldg_b
        double precision g_c(ldg_c, nr, nc), g_b(ldg_b, nr, nc)

C     dummy variables
        integer i, j, k, g_pmax_, g_i_
        integer g_ehfid
        data g_ehfid /0/
C
        g_pmax_ = (1+nex+ncov)*nex
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do i = 1, nr
          do j = 1, nc
            do g_i_ = 1, g_p_
              g_c(g_i_, i, j) = -g_b(g_i_, i, j)
            enddo
            c(i, j) = a(i, j) - b(i, j)
C--------
          enddo
        enddo
C
        return
      end
C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
C
      subroutine g_matrixtranspose(g_p_, nar, nac, a, g_a, ldg_a, at,
     *  g_at, ldg_at)
C     used to transpose a real matrix(or vector)
C
      use variable
      implicit none
C
C     on the entry variables
        integer g_p_,nar, nac, ldg_at, ldg_a
        real*8 a(nar, nac), at(nac, nar)
        double precision g_at(ldg_at, nac, nar), g_a(ldg_a, nar, nac)
C
C     dummy variables
        integer i, j, g_pmax_, g_i_
        integer g_ehfid
        data g_ehfid /0/
C

C
        g_pmax_ = (1+nex+ncov)*nex
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do i = 1, nac
          do j = 1, nar
            do g_i_ = 1, g_p_
              g_at(g_i_, i, j) = g_a(g_i_, j, i)
            enddo
            at(i, j) = a(j, i)
C--------
          enddo
        enddo
C
        return
      end
C
C
Ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrccox_likelihood(nd,betval,fcn_ran)

      use variable
      implicit none

c     variables on the entry
      integer nd
      real*8 betval(nd),fcn_ran

      call rrc_cum_cox(nmain_new,nage_main,nd,betval,fcn_ran)
      fcn_ran=-fcn_ran
      print *, '-2Log L=', 2.0d0*fcn_ran
      print *, 'betval=', betval
      write(*,*)

      return
      end

ccccccssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine rrc_cum_cox(nmax,ntm,nd,betval,fcn_ran)

      use variable
      implicit none

c     variables on the entry
      integer nmax, ntm, nd
      real*8 betval(nd), fcn_ran


c     dummy variables
      integer i,j,k,j1,k1,indk,kch
      real*8 eps
      parameter(eps=1.0d-8)


      double precision, dimension(:), allocatable ::
     +  s0
      double precision, dimension(:,:), allocatable ::
     +  rrf

      real*8 temp
      allocate(rrf(nmax,ntm))
      allocate(s0(nageobs))

      do i=1,nmax
         do k=1,j_main(i)
            temp=0.0d0
            do j1=1,nex
               temp=temp+betval(j1)*zs_common(i,k,j1)
            enddo
            if (ncovm.ge.1) then
               do j1=1,ncovm
                  temp=temp+betval(nex+j1)*
     *                 dble(covar_common(i,k,ncove+j1))
               enddo
            endif
            rrf(i,k)=dexp(temp)
         enddo
      enddo
      
      do k=1,nfa
         s0(k)=0.0d0
         do i=1,nmax
            do j=1,j_main(i)
               if (age_main(i,j).eq.tf(k)) then
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
               if (del_main(j).eq.1.0d0) then
                  do k1=1,nex
                     fcn_ran=fcn_ran+betval(k1)*
     *                    zs_common(j,j_main(j),k1)
                  enddo

                  if (ncovm.ge.1) then
                     do k1=1,ncovm
                        fcn_ran=fcn_ran+betval(nex+k1)*
     *                       dble(covar_common(j,j_main(j),ncove+k1))
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

      call g_rrc_cum_cox(nmain_new,nage_main,nd,betval,g_fval)
      g_fval=-g_fval

      print *, 'Gradient=', g_fval

      return
      end

ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c     score of Cox partial likelihood

      subroutine g_rrc_cum_cox(nmax,ntm,nd,betval,fcn_ran)

      use variable
      implicit none

c     variables on the entry
      integer nmax, ntm, nd
      real*8 betval(nd), fcn_ran(nd)


c     dummy variables
      integer i,j,k,j1,k1,indk,kch
      real*8 eps
      parameter(eps=1.0d-8)


      double precision, dimension(:), allocatable ::
     +  s0
      double precision, dimension(:,:), allocatable ::
     +  rrf, s1

      real*8 temp
      allocate(rrf(nmax,ntm))
      allocate(s0(nageobs))
      allocate(s1(nageobs,nd))

      do i=1,nmax
         do k=1,j_main(i)
            temp=0.0d0
            do j1=1,nex
               temp=temp+betval(j1)*zs_common(i,k,j1)
            enddo
            if (ncovm.ge.1) then
               do j1=1,ncovm
                  temp=temp+betval(nex+j1)*
     *                 dble(covar_common(i,k,ncove+j1))
               enddo
            endif
            rrf(i,k)=dexp(temp)
         enddo
      enddo

      s0=0.0d0
      do k=1,nfa
         do i=1,nmax
            do j=1,j_main(i)
               if (age_main(i,j).eq.tf(k)) then
                  s0(k)=s0(k)+rrf(i,j)
               endif
            enddo
         enddo
      enddo

      s1=0.0d0
      do k=1,nfa
         do j1=1,nex
            do i=1,nmax
               do j=1,j_main(i)
                  if (age_main(i,j).eq.tf(k)) then
                     s1(k,j1)=s1(k,j1)+zs_common(i,j,j1)*rrf(i,j)
                  endif
               enddo
            enddo
         enddo
         if (ncovm.ge.1) then
            do j1=1,ncovm
               do i=1,nmax
                  do j=1,j_main(i)
                     if (age_main(i,j).eq.tf(k)) then
          s1(k,nex+j1)=s1(k,nex+j1)+covar_common(i,j,ncove+j1)*rrf(i,j)
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo


      fcn_ran=0.0d0
      do i=1,ndist

         if (evcnt(i).ne.0.0d0) then
            do j1=1,nfa
               if (t_main(ilow(i)).eq.tf(j1)) indk=j1
            enddo

            do j=ilow(i),iup(i)
               if (del_main(j).eq.1.0d0) then
                  do k1=1,nex
                     fcn_ran(k1)=fcn_ran(k1)+zs_common(j,j_main(j),k1)
                  enddo

                  if (ncovm.ge.1) then
                     do k1=1,ncovm
                        fcn_ran(nex+k1)=fcn_ran(nex+k1)+
     *                       dble(covar_common(j,j_main(j),ncove+k1))
                     enddo
                  endif
               endif
            enddo

            do k1=1,nd
               fcn_ran(k1)=fcn_ran(k1)-
     *              evcnt(i)*(s1(indk,k1)/dmax1(s0(indk),eps))
            enddo
            
         endif

      enddo
      deallocate(rrf)
      deallocate(s0)
      deallocate(s1)
      
      return
      end

c
Ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine rrchessian_cox(nd,betval,hess_fval,ldh)

      use variable
      implicit none

c     variables on the entry
      integer nd,ldh
      real*8 betval(nd), hess_fval(ldh,nd)
      
c     dummy variables

      call g_g_rrc_cum_cox(nmain_new,nage_main,nd,betval,hess_fval)

      print *, "hessian=",hess_fval

      return
      end


ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c     hessian of Cox partial likelihood
c     Xiaomei is working on it before.

      subroutine g_g_rrc_cum_cox(nmax,ntm,nd,betval,fcn_ran)

      use variable
      implicit none

c     variables on the entry
      integer nmax, ntm, nd
      real*8 betval(nd), fcn_ran(nd,nd)


c     dummy variables
      integer i,j,k,j1,k1,indk,kch
      real*8 eps, temp
      parameter(eps=1.0d-8)

      double precision, dimension(:), allocatable ::
     +  s0, temp2
      double precision, dimension(:,:), allocatable ::
     +  rrf, s1
      double precision, dimension(:,:,:), allocatable :: s2



      allocate(rrf(nmax,ntm))
      allocate(s0(nageobs))
      allocate(s1(nageobs,nd))
      allocate(s2(nageobs,nd,nd))
      allocate(temp2(nd))

      do i=1,nmax
         do k=1,j_main(i)
            temp=0.0d0
            do j1=1,nex
               temp=temp+betval(j1)*zs_common(i,k,j1)
            enddo
            if (ncovm.ge.1) then
               do j1=1,ncovm
                  temp=temp+betval(nex+j1)*
     *                 dble(covar_common(i,k,ncove+j1))
               enddo
            endif
            rrf(i,k)=dexp(temp)
         enddo
      enddo
      
      s0=0.0d0
      do k=1,nfa
         do i=1,nmax
            do j=1,j_main(i)
               if (age_main(i,j).eq.tf(k)) then
                  s0(k)=s0(k)+rrf(i,j)
               endif
            enddo
         enddo
      enddo

      s1=0.0d0
      do k=1,nfa
         do j1=1,nex
            do i=1,nmax
               do j=1,j_main(i)
                  if (age_main(i,j).eq.tf(k)) then
                     s1(k,j1)=s1(k,j1)+zs_common(i,j,j1)*rrf(i,j)
                  endif
               enddo
            enddo
         enddo
         if (ncovm.ge.1) then
            do j1=1,ncovm
               do i=1,nmax
                  do j=1,j_main(i)
                     if (age_main(i,j).eq.tf(k)) then
          s1(k,nex+j1)=s1(k,nex+j1)+covar_common(i,j,ncove+j1)*rrf(i,j)
                     endif
                  enddo
               enddo
            enddo
         endif
      enddo

      s2=0.0d0
      do k=1,nfa
         do i=1,nmax
            do j=1,j_main(i)
               if (ncovm .ge. 1) then
                  temp2(1:nd)=(/zs_common(i,j,1:nex),
     *                 covar_common(i,j,(ncove+1):(ncove+ncovm))/)
                  if (age_main(i,j).eq.tf(k)) then
                     do j1=1,nd
                        do k1=1,nd
                           s2(k,j1,k1)=s2(k,j1,k1)+temp2(j1)*temp2(k1)
     *                          *rrf(i,j)
                        enddo
                     enddo
                  endif
               else
                  temp2(1:nex)=zs_common(i,j,1:nex)
                  if (age_main(i,j).eq.tf(k)) then
                     do j1=1,nex
                        do k1=1,nex
                           s2(k,j1,k1)=s2(k,j1,k1)+temp2(j1)*temp2(k1)
     *                          *rrf(i,j)
                        enddo
                     enddo
                  endif
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

            do j1=1,nd
               do k1=1,nd
                  fcn_ran(j1,k1)=fcn_ran(j1,k1)
     *       +evcnt(i)*s2(indk,j1,k1)/dmax1(s0(indk),eps)  
     *       -evcnt(i)*(s1(indk,j1)*s1(indk,k1)/dmax1(s0(indk)**2,eps))
               enddo
            enddo
            
         endif

      enddo
      deallocate(rrf)
      deallocate(s0)
      deallocate(s1)
      deallocate(s2)
      deallocate(temp2)
      
      return
      end

Ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss


        subroutine mvlm2(p,q,n,y,xin,gamma,r)
        implicit none
c
c this routine computes the parameter estimates and their variance-covariance
c matrix for the multivariate linear model of p-vector y on (p+q)-vector xin
c the slopes are in (p+q)xp matrix g; aprim has the intercepts
c
c       n = number of observations in the dataset
c
c       p = dim(Y) ! which are the true covariates
c
c       q = dim(X) ! which are the error measured covariates and other explanatory variable
c
       integer n,p,q,q1
       real*8 y(n,p),xin(n,q),gamma(q+1,p), r(p)

       integer, dimension(:), allocatable ::
     +  ipiv

       double precision, dimension(:), allocatable ::
     +   aprime,xbar,ybar, yhat

       double precision, dimension(:,:), allocatable ::
     +   g,x,xtx,xix,vargam,e,proj,yp,s,xtx2,
     +   temp,isubn,eye

       double precision part1, part2, part3, xsum, ysum

       integer i, j, ir, ic, il, info
       external dgesv

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
       allocate(eye(q+1,q+1))
       allocate(ipiv(q+1))
       allocate(yhat(n))	

c
c add 1 to the dimension of the problem to account for intercepts
c
c       augment x matrix with a vector of 1's
c
c
        q1=q+1
        do 10 i=1,n
        do 11 j=1,q
                x(i,j+1)=xin(i,j)
11      continue
                x(i,1)=1.d0
10      continue

c get means of x,y
c
c       do 83 j=1,q+1
c         xb(j)=0.d0
c83     continue
c       do 82 i=1,n
c       do 82 j=1,q+1
c         xb(j)=xb(j)+x(i,j)
c82     continue
c       do 87 j=1,q+1
c               xb(j)=xb(j)/float(n)
c87     continue
c       write(*,*) ' means of x ',xb
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
c       call dmxtxf(n,q1,x,n,q1,xtx,q1)
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
          gamma(ir,ic) = 0.0
          do il = 1,n
            gamma(ir,ic) = gamma(ir,ic) + xix(ir,il)*y(il,ic)
          enddo
        enddo
        enddo
c       call dmrrrr(q1,n,xix,q1,n,p,y,n,q1,p,gamma,q1)
c
c       Compute r-square
	do j=1,n
	yhat(j)=0.0d0
        do ir=1,q1
        do ic=1,p
	   yhat(j)=yhat(j)+gamma(ir,ic)*x(j,ir)
	enddo
        enddo          
        enddo

       do i=1,p
           ysum=0
           do j=1,n
              ysum=ysum+y(j,i)
           enddo
           ybar(i)=ysum/n
        enddo 
	
	do i=1,p
           part1=0
           part3=0
           do j=1,n
              part1=part1+(y(j,i)-yhat(j))**2
              part3=part3+(y(j,i)-ybar(i))**2
           enddo
           r(i)=1.0d0-part1/part3
        enddo

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
       deallocate(eye)
       deallocate(yhat)	

       return
       end

Ccccccsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine data_summary
c     this subroutine is used to summarize the data statistics.

      use variable
      implicit none

      integer i,nd,k

      double precision, dimension(:), allocatable ::
     +  naivci_left, naivci_right, or_naiv, rrcci_left,
     +  rrcci_right, or_rrc, pval
      character(len=10), dimension(:), allocatable ::
     +  names
      character*255 filename

      nd=nex+ncovm
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
     * 'This analysis is for',month1,'-',month,'month average exposure'
      write(outfile1,*)
     * 'with general linear measurement error model.'
      write(outfile1,*) ''
      write(outfile1,*) 'The increment for each exposure is:',exunit
      write(outfile1,*) ''


      if ((nex.eq.1).and.(ncovm.eq.0).and.(isres.eq.0)) then
         write(outfile1,*)
     *'There is only one exposure in the analysis and no covariate.'
         write(outfile1,*) 'The main exposure is measured with error.'
      else if ((nex.eq.1).and.(ncovm.ge.1)) then
          write(outfile1,*)      'There is one exposure and one or more
     *covariates in the analysis.'
          write(outfile1,*)   'The exposure is measured with error, but
     *the covariates are error-free.'
      endif

      if ((isres.eq.1)) then
          write(outfile1,*)
     *"The residual analysis has been used in this analysis"
          write(outfile1,*) "to control for confounding."
      endif

      write(outfile1,*) ''

      write(outfile1,*) '# of subjects in the main study:', nmain_newcp
      write(outfile1,*) '# of person-time observations in the main',
     *     ' study:', mainobscp
      write(outfile1,*) '# of cases in the main study:',
     *     nint(sum(del_main(1:nmain_newcp)))

      write(outfile1,*) '# of unique failure time:', nfa
      write(outfile1,*) ''
      write(outfile1,*) 'The unique failure times are:'
      do k=1,int(dble(nfa)/10.0d0)
         write(outfile1,*) tf((1+(k-1)*10):(10*k))
      enddo
      write(outfile1,*) tf((1+(k-1)*10):nfa)
      write(outfile1,*) ''

      write(outfile1,*) '# of subjects in the validation study:',nval
      write(outfile1,*) '# of person-time observations in the ',
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
     *acom(k,1,1),acom(k,2:(1+nex),1),rsquare(k,1)
      enddo

 200  format(I4,6X,2(F6.1,2X),I8,6X,3(F7.3,2X))

      write(outfile1,*) ''
      write(outfile1,*) '----------------------------------------------'
      write(outfile1,*) 'Results for uncorrected estimates (naive cox)'
      write(outfile1,*) 'and corrected estimates (RRC):'
      write(outfile1,*) '----------------------------------------------'
      write(outfile1,*) ''


      names(1:nex)=expname(1:nex)
      if (nd >= nex+1) then
         names((nex+1):nd)=covname(1:ncovm)
      endif      


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
         write(1016,102) trim(names(i)), betaes_rrc(i),serror_rrc(i),
     *  pval(i),or_rrc(i),rrcci_left(i),rrcci_right(i),groupnum,ngrp
      enddo


      write(outfile1,*) ''
      write(outfile1,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      close(1016)

 90   format(1X,A,6X,A7,1X,A6,2X,A5,2X,A5,2X,A14)
 100  format(1X,A,4X,F6.3,2X,F6.3,2X,F5.3,2X,F5.3,2X, A, F6.2, A, F6.2,
     *     A)
 101  format(A,',uncorrected',6(',',F6.3),2(',',I4))
 102  format(A,',rrc',6(',',F6.3),2(',',I4))


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

      real*8 xp, bound
      integer istat
      external cdfchi

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










