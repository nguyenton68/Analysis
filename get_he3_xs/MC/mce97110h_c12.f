	program mce97110h
!________________________________________________________________________
! Check sieve IN
! Check collimator
! Check foil
! 11/14/2016
! Nguyen modified formula for scattering angle, without approximation

! Monte-Carlo of HRS hadron/right arm  spectrometer using uniform illumination.
! Nilanga Liyanage: took Alex's MC program and modified it for the broken setpum case
! please see the README file
! Since there is no magnetic model of the broken septum, the forward matrix elements
! used here from a fit to data.
! No transport through the spectrometer done, no multiple scattering effects taken into 
! account. 
!This  program has the structure the ancient SLAC monte-carlo program, uniform ( David Potterveld, Sept. 1985). 
!major changes:
!  - For UNIX and implicit none
!  - Not using the same coordinate systemes (now we are right handed)
!  - HRS setting and optics
!  - polarized 3He target setting
!  - get ride of the TopDraw setup (now PAW is used to plot output)
!  - new random numbers generator
!  - Raster option 
!  - Physics: C12 elastic cross sections
!             Radiative corrections and landau tail for elastic peak
!             3He quasi-elastic cross sections, asymmetries and external RC 
!A. Deur June 1999
!This version includes the septum magnet. March 2003
!Includes target collimators. April 2003

! N. Ton Jul 2017
! - This version includes the miswired septum magnet. 
! - Target collimators.
! - Modified formula for scattering angle, without approximation
! - Cross section for 3He
!________________________________________________________________________

	implicit none

	real*4 d_r,pi,r_d,res_x,res_y,res_th,res_ph,root,flagsv
C Type declarations
	character*21 qelname
	character*13 qelfile

	integer i,n_trials
	integer trial,foil

	real rndm,rand      !random number generator must link to CERN libs
	integer*4 today(3), now(3)
	integer*4 nseed

	real*4 Etest,Etest2,eqel
	real*4 Ei,Ep,Eps,th_spec !incident energy,scat energy,spectro setting
	real*4 inl,outl ! energy losses (incoming and outgoing) in rad length 
	real*4 radw			!temp. radiation length
	real*4 rc,dE1,dE2,dE3,dE4,detot,w ! if rc=0, no elastic rad cor,
        real*4 aav,avz,zav,ang,angstrag,ioni1,ioni2  ! if rc=1, rad cor turned on

	real*4 cos_HRS,sin_HRS,tan_HRS
	real*4 x,xtest,y,ytest,z,t1,t2,t3,t4,tt,u1,u2,etemp!temporaries
	real*4 t5,t6,t7,t8,et !more temporaries
	real*4 xtr,ytr,ztr			!HRS transport coords. (cm)
	real*4 xfoc,yfoc,phfoc,thfoc	!HRS coords. (cm) in focal plane
	real*4 deltat,thetat,phit,yt,zreact !HRS coords. (cm) at target (reconstructed)
	real*4 rec_th, rec_ph,rec_y !reconstructed, not use for w calculation
	real*8 dpor,thor,phor,yor,Epor   !HRS coords. (cm) at target (original)
	real*8 wor
	real*4 dpnstrg ! delta P before straggling      
	real*4 dpp_ac				!HRS deltap/p (percent)
	real*4 dth_ac				!HRS delta theta (mr)
	real*4 dph_ac                           !HRS delta phi   (mr)
	real*4 dpp_max,dpp_min			!cut deltap/p (percent)
	real*4 dth_max,dth_min			!cut delta theta (mr)
	real*4 dph_max,dph_min			!cut delta phi   (mr)
	real*4 tgt_max,tgt_min
	real*4 w_max,w_min
	real*4 dpp				!rnd deltap/p (percent)
	real*4 dth				!rnd delta theta (mr)
	real*4 dph				!rnd delta phi   (mr)

	real*4 bcol,ccol,ecol,fcol,tanphi !variables for target collimators
	real*4 offbc,offcc,offec,offfc,ycoll !offsets collimators
	integer sieve

	real spot_x,spot_y,tgt_l,xoff,yoff,zoff,xbeam,ybeam
	real asym_calc3,cross_section,aspin,choice,xs,asy,avasy ! variables for Xsection and asy computation
	real xom,sigma,qsq,atp,atlp,atpp,atpn,atlpp,atlpn,al,app ! variables for quasi-elastic Xsection and asy computation
	REAL XSPREV,asyprev,xomprev !temporaries for quasi-elastic interpolation
       
	real tacc,xsc,cacc,full_sa,xscount,xs_sum,cutc	! acceptance and cross section counters
	real r_phi,r_theta !angle in radian to calculate real scattering angle
	real QQ,mott
	real c_acc,c_tacc
! Function for rad corr
	real incl,intrad,ioni!,deltazero,deltazero1,deltazero2
	real xdi,xdo !average density*thickness for incoming and outgoing e-
	common/landau/ld(2,68)
	real ld ! contain the landau tail distribution
c
c NL arrays for the corners of the solid angle
c
	real y_corner(6)
	integer ii
	real ph_corner(6)
	real edge(6)
	real slope, intercept

C NT: array for corners of phase space cut
C Use to determine analysis cut in simulation
	real ph_pos(4)
	real th_pos(4)
	real ps_edge(4)
	real ps_slope,ps_intercept
	integer ll
C NT: array for corner of xs cut
	real xs_th(4)
	real xs_ph(4)
	real xs_edge(4)
	real xs_slope,xs_intercept
	integer kk

c
C Coefficient for forward matrix
	real coef_x(14),coef_y(15),coef_phi(12),coef_th(10)
	integer tmp_x,tmp_y,tmp_phi,tmp_th

c	
!!!declare stuff to be able to make Raw Wise Ntuple!!

	common/pawc/blanc(150000)
	real blanc 
        real ztup(26)
        call initpaw

!zero acceptance and cross section counters
        tacc=0.
	xsc=0.
	cacc=0.
	avasy=0.

C Math constants
	pi = 3.141592654
	d_r = pi/180.
	r_d = 180./pi
	root = 0.707106781 !square root of 1/2
C HRS Wire chamber resolutions (half width at one sigma.)

	res_x = 0.013	!cm
	res_y = 0.013	!cm
	res_th = 0.3	!mr
	res_ph = 0.3	!mr

C Beam position offsets 2.134 GeV, 6 degrees.
	
	xoff = -0.05043! cm
	yoff = 0.05985!cm
	zoff = -0.057 ! cm foil position offset

	open (20,file='c12.inp',status='old')
	open (25,file='cutc12.inp',status='old')
! Random seed
        call idate(today) ! select a seed for rndm
        call itime(now)
        nseed =today(1)+today(2)+now(1)+now(2)+now(3) ! 
        u1=rand(nseed)
        print*,'Print u1= ',nseed,' ',u1

! Read data lines

	read (20,*) n_trials,Ei,Ep,th_spec,dpp_ac,dth_ac,dph_ac,
     +  spot_x,spot_y,tgt_l,aspin,choice,inl,outl,xdi,xdo,rc
	read (25,*)dpp_max,dpp_min,dth_max,dth_min,dph_max,dph_min,
     +  tgt_max,tgt_min,w_max,w_min
        if (tgt_max.lt.tgt_min) then ! security
           print*, 'tgt_max smaller than tgt_min'
        endif
         if (dpp_max.lt.dpp_min) then
           print*, 'dpp_max smaller than dpp_min'
        endif
          if (dph_max.lt.dph_min) then
           print*, 'dph_max smaller than dph_min'
        endif

 	if (dth_max.lt.dth_min) then
           print*, 'dth_max smaller than dth_min'
        endif
	if (tgt_max.gt.tgt_l/2.) then
	   print*, 'target is not completely illuminated'
	endif
	if (abs(tgt_min).gt.tgt_l/2.) then
	   print*, 'target is not completely illuminated'
	endif
	if (dpp_max.gt.dpp_ac/2.) then
	   print*, 'momentum range is not completely illuminated'
	endif
	if (abs(dpp_min).gt.dpp_ac/2.) then
	   print*, 'momentum range is not completely illuminated'
	endif
	if (dph_max.gt.dph_ac/2.) then
	   print*, 'phi range is not completely illuminated'
	endif
	if (abs(dph_min).gt.dph_ac/2.) then
	   print*, 'phi range is not completely illuminated'
	endif
	if (dth_max.gt.dth_ac/2.) then
	   print*, 'theta range is not completely illuminated'
	endif
	if (abs(dth_min).gt.dth_ac/2.) then
	   print*, 'theta range is not completely illuminated'
	endif			! end of security checks

! Read coefficient for forward matrix
	open(30,file='infiles/x_coef.inp',status='old')
	open(31,file='infiles/y_coef.inp',status='old')
	open(32,file='infiles/phi_coef.inp',status='old')
	open(33,file='infiles/th_coef.inp',status='old')
	do tmp_x = 1,14
	   read(30,*) coef_x(tmp_x)
	enddo
	do tmp_y = 1,15
	   read(31,*) coef_y(tmp_y)
	enddo
	do tmp_phi = 1,12
	   read(32,*) coef_phi(tmp_phi)
	enddo
	do tmp_th = 1,10
	   read(33,*) coef_th(tmp_th)
	enddo

! End read coefficent	   
 	
	etemp=Ei

	
	th_spec = th_spec*d_r ! in rad
	cos_HRS = cos(th_spec)
	sin_HRS = sin(th_spec)
	tan_HRS = sin_HRS/cos_HRS
	spot_x=spot_x*100. ! in cm.
	spot_y=spot_y*100. ! in cm.
	
	!Keep cuts on target length (z) and not on y_target 12/20/2005
	!tgt_min=tgt_min*abs(sin_HRS)!tgt_min is now the cut on y_target (HRS
	!tgt_max=tgt_max*abs(sin_HRS)!frame) not a cut on target length (z) anymore 
	if (rc.eq.1.) then ! fill the common block for landau distribution
	   call stuffit
	endif
C look at physics database
****elastic is in a subroutine *************************
***quasi-elastic****************************************
**this is now not valid anymore***************
	if (choice.eq.2) then	!quasi-elastic
	   etest2=0.
	   open(23,file='qel.dir/qel_en.dat',status='old')
	do i=1,100 ! chose the right incident energy file
 	   read(23,*,end=998) Etest,qelfile
	   if (abs(etest-Ei).le.abs(etest2-Ei)) then
	   eqel=Etest
	   qelname= 'qel.dir/' //  qelfile
	   endif
	   etest2=etest
 998	continue
	enddo

	open(24,file=qelname,status='old')!contains the quasi-elastic data
	!print*,qelname 
	endif
***end of quasi-elastic****************************************
	if (choice.eq.0) then   !phase-space
c
c Calculate the full solid angle for the box we have defined
c dp values are given as percent, need to be divided by 100
c theta and phi are given as mrad, need to be divided by 1000.	   
c
	   full_sa = (Ep*dpp_max-Ep*dpp_min)/100.0
     +               *(dth_max-dth_min)/1000.0
     +               *(dph_max-dph_min)/1000.0     
c
	   write(6,*) "Full solid angle is ",full_sa," Gev. st_rad"
	endif
c
	if (choice.eq.1) then   !phase-space for elastic
c
c Calculate the full solid angle for the box we have defined
c	   
	   full_sa = (dth_max-dth_min)/1000.0
     +               *(dph_max-dph_min)/1000.0     
c
	   write(6,*) "Full solid angle is ",full_sa," st_rad"
	endif

c
C------------------------------------------------------------------------C
C                       Top of Monte-Carlo loop                          C
C------------------------------------------------------------------------C
	do trial = 1,n_trials

C Pick two indep., normal dist. numbers, used only if raster is off.
 56	u1=rand()!rndm(trial)
	u2=rand()!rndm(trial+1)
!	t1=sqrt(-2.*log(u1))*cos(2*pi*u2) !t1 norm dist (cf notes Fonvieille on  MC)
!	u1=rndm(trial+n_trials) 
!	u2=rndm(n_trials+trial+1)
!	t2=sqrt(-2*log(u1))*cos(2*pi*u2)

C Pick starting point within target. Z is picked uniformly, X and Y are
C chosen as truncated gaussians. (Truncated to 3 sigma) if raster off
C If CIRCULAR raster is on, uniform dist. 
C Units are cm.
     
!	x = t1*spot_x+xoff	! if raster off, unit in cm
!	y = t2*spot_y+yoff	! if raster off, in cm

!	if (((2*u1-1)**2.+(2*u2-1)**2.).gt.1.) then! if circular raster on (uniform dist), 
                   ! in cm. for a squared uniform dist., comment the 3 lines 
!	goto 56
!	endif

	x=(u1-0.5)*spot_x+xoff !*2.
 	y=(u2-0.5)*spot_y+yoff !*2.
	xbeam=x
	ybeam=y

	z=rand()!rndm(trial+2)
	z=(z-0.5)*tgt_l ! in cm 
	
C       Transform from target to HRS transport coordinate 

!	   xtr = -y
!	   ytr = x*cos_HRS - z * sin_HRS ! changed -x to +x 12/07/2005
	   !ztr = z * cos_HRS - x * sin_HRS

C Pick scattering angles and DP/P after to have computed energy losses
C picks a gaussian dp/p center around dp/p_scat
C pick a gaussian for Beam energy dispertion
	    Ei=Etemp !last Ei was shift by energy loss and dispertion. Come back at the set one 
C********************
C This is where problem with radiation correction come from
C If I use rand(), it will create negative energy
C If I use rndm, no problem 
	    u1=rndm(trial+10)
	    u2=rndm(trial+11)
  	    t1=sqrt(-2.*log(u1))*cos(2.*pi*u2) ! gaussian dist

	    Ei=t1*3*10.**(-5.)*Ei+Ei ! beam dispertion is taken as 3.e-5
	    Et=Ei		! used in W - M calculation
	    
	    
	    dth =rand()!rndm(n_trials+4)
	    dth =(dth-0.5)*dth_ac

	    dph =rand()!rndm(n_trials+5) 
            dph =(dph-0.5)*dph_ac



	    xtr = -y -dth/1000.0*(z*cos_HRS/cos(atan(dph/1000.0)))!transport coords
	    ytr=sin(th_spec+atan(dph/1000.0))/cos(atan(dph/1000.0))*
     >           (x/tan(th_spec+atan(dph/1000.0))-z)

! Without approximation
	    ang=acos((1+abs(tan_HRS)*tan(dph/1000.0))
     +  /sqrt((1+tan_HRS**2)*(1+tan(dph/1000.0)**2+tan(dth/1000.0)**2)))

! With approximation in small angle
C	    ang=acos((cos_HRS + abs(sin_HRS)*dph/1000.)
C     +          /sqrt(1+(dph/1000.)**2.+(dth/1000.)**2.))

	      if (rc.eq.1.) then ! Rad cor enable. pick up an
		 !!energy loss on incoming particle!!
	         ! pick up an dEi energy loss du to external 
		 aav=18.091 ! average A for incoming e-
		 zav=8.787  ! average Z for incoming e-
		 dE1=incl(Ei,zav,inl)  
		 Ei=Ei-dE1 
	         ! pick up an dEi energy loss du to internal 	 	 
 		 dE2=intrad(Ei,abs(ang))
		 ioni1=ioni(Ei,aav,zav,xdi)! ioni is loss by ionization 
                 !deltazero1=deltazero(Ei,aav,zav,xdi)!most probab loss by ioni
		 Ei=Ei-dE2-ioni1 !-0.000207

                 !now chose or compute scattered particle energy
		 if (choice.eq.1) then ! elastic
		 !have to compute Eps with new Ei
	         Eps=Ei/(1+2*Ei/2.808*(sin(ang/2.))**2.)!2.808 GeV mass He3
	         Epor=Etemp/(1+2*Etemp/2.808*(sin(ang/2.))**2.)!without RC and resolution
		  else ! quasi elastic or phase space 
	          dpp =rand()!rndm(n_trials+3)
	          Eps =Ep*(1+(dpp-0.5)*dpp_ac/100.)
	          Epor =Eps !without RC and resolution
		 endif   
		 !!!!!!!

		 ! now outgoing e- energy loss due to internal
 		 dE3=intrad(Eps,abs(ang)) 
		 Eps=Eps-dE3	
		! external loss of outgoing e-
		 aav=22.503 ! average A for outgoing e-
		 zav=11.257  ! average Z for outgoing e-

		 dE4=incl(Eps,zav,outl)
		 ioni2=ioni(Eps,aav,zav,xdo)
                 !deltazero2=deltazero(Eps,aav,zav,xdo)
		 Eps=Eps-dE4-ioni2 ! -0.00007
	         ! Now compute the final relative momentum dpp with resolution
!                  Will correct for HRS mom resolution when vdc resolution is 
!                  taken care of on lines 650-676
	         u1=rndm(trial+16)
	         u2=rndm(trial+17)
  	         t1=sqrt(-2.*log(u1))*cos(2.*pi*u2) ! gaussian dist
	         t1=t1*4.0*10.**(-4.)*Ep !relative HRS mom res taken as 3e-4
		 Eps=Eps+t1
		 dpp=(Eps-Ep)*100./Ep	!   dpp has to be in % 
	                   !   (converted in fractional value in 
	                   !    subroutine scavec for JJL function)
		  !print*, de1,de2,de3,dpp

	!write energy loss and equialent w to file. Effect of the incoming 
        !energy loss is corrected for recoil.
		 dEtot=(de1+de2)/(1+2*Ei/2.808*(sin(ang/2.))**2.)+de3+de4
     ++ioni1+ioni2

	    else ! rad Cor not on
	       u1=rand()!rndm(trial+16)
	       u2=rand()!rndm(trial+17)
	       t1=sqrt(-2.*log(u1))*cos(2.*pi*u2) ! gaussian dist
	       t1=t1*4.0*10.**(-4.)*Ep !relative HRS mom res taken as 3e-4
	       if (choice.eq.1) then ! elastic
	        Epor=Etemp/(1+2*Etemp/2.808*(sin(ang/2.))**2.)!without RC and resolution 
	        Eps=Ei/(1+2*Ei/2.808*(sin(ang/2.))**2.)!2.808 GeV mass He3
	        Eps=Eps+t1
	        dpp=(Eps-Ep)*100./Ep 
	        else ! quasi elastic or phase space 
	        dpp =rand()!rndm(n_trials+3)
	        Eps =Ep*(1+(dpp-0.5)*dpp_ac/100.)
	        Epor =Eps!without RC and resolution
	        dpp=(Eps-Ep)*100./Ep 
	       endif 

	    endif ! end of energy loss computation due to RC


C origin coordinates in the HRS frame

	   yor=ytr
	   phor=dph
	   thor=dth
	   dpor=(Epor-Ep)*100./Ep ! otherwise, shifted by radiative proces

	   dpnstrg=dpp
	   wor=sqrt(2.808**2+2*2.808*(Etemp-Epor)
     +	   -4*Etemp*Epor*(sin(-1*ang/2.))**2.)-2.808
C
C  Calculate reactz from reconstructed quantities at target
	  zreact = -(yor/100.0)*cos(atan(phor/1000.0))/(sin(th_spec + 
     +    atan(phor/1000.0)))+xbeam/(100.0*tan(th_spec 
     +    + atan(phor/1000.0)))

C Project trajectory to sieve slit
C sieve slit box to sieve slit location

c$$$	tt = 80.06		! cm 
c$$$	call project(xtr,ytr,tt,dth,dph)
c$$$
c$$$	   sieve=1
c$$$	   if (sieve.eq.1) then ! sieve slit is in
c$$$	      call sieve_slit(xtr/100.0,ytr/100.0,flagsv)
c$$$	      if (flagsv.lt.1.0) then
c$$$		 goto 500
c$$$	      endif
c$$$	   endif
C Cut on momentum range for phase space: -3.5 to 4.5
        if(dpor.le.-3.5.or.dpor.ge.4.5) then
	   goto 500
	endif


*------------------------------------------------------------------------*
***** This is for count tacc
* Here use mrad unit: phor, thor and xs_ph, xs_th: mrad
	  xs_ph(1) =  14.405!13.504
	  xs_ph(2) =   3.680!4.595
	  xs_ph(3) =  -9.669!-8.688
	  xs_ph(4) =   0.861!3.594

	  xs_th(1) = 44.359!44.876
	  xs_th(2) = 29.856!31.774
	  xs_th(3) = 37.212!36.295
	  xs_th(4) = 53.972!55.768
	do kk=1,4
	   if(kk==4) then
	      xs_slope=(xs_th(1)-xs_th(kk))/(xs_ph(1)-xs_ph(kk))
	   else
	      xs_slope=(xs_th(kk+1)-xs_th(kk))/(xs_ph(kk+1)-xs_ph(kk))
	   endif
	   xs_intercept=xs_th(kk)-xs_slope*xs_ph(kk)
	   xs_edge(kk)=xs_slope*phor+xs_intercept
	enddo

	if(thor.ge.xs_edge(1).and.thor.ge.xs_edge(2).and.thor.le.
     + xs_edge(3).and.thor.le.xs_edge(4).and.wor.le.0.01)then
	   tacc=tacc+1
	endif
C Cut on z_target: +/- 18cm
!        if(z.le.-18.0.or.z.ge.18.0) then
!	   goto 500
!	endif
c use this cut for foil target runs only, 3 12C foils at +/- 10 cm and 0
c
c	  if ((abs(zreact-0.1).ge.0.01).and.(abs(zreact+0.1).ge.0.01)
c     +         .and.(abs(zreact).ge.0.01))goto 500

c NL: Apply box cuts given in the cuts file for zreact,th and ph

c Ton: make sure use correct collimator for each case
* optic: sieve collimator: only cut for optic
* acceptance: target collimator
c Apply 6 degrees collimator block cuts
c
	if ((abs(th_spec).ge.(5.75*d_r)).and.
     @     (abs(th_spec).le.(6.25*d_r))) then    
!!!!!!!!ice cone collimator!!!!!!!!!
!C check if it going through the middle collimator
	   ecol=18.99		!position collimator along the beam line
	   fcol=2.326		!dist from collimator to beam axis
	   tanphi=(fcol+x)/abs(z-ecol)/cos(thor/1000.)
	   if((tan(abs(th_spec+phor/1000.))).ge.tanphi) goto 500

!C check if it going through the downstream window collimator
	   bcol=33.96		!z position of collimator
	   ccol=2.85		!dist top of collimator to beam axis
	   tanphi=(ccol+x)/abs(bcol-z)/cos(thor/1000.)
	   if((tan(abs(th_spec+phor/1000.))).le.tanphi) goto 500
c	    endif !end of 6 degree 

	endif			!end of 6 degree 



C
C Now calculate the real focal plane quantities for the broken septum
C case using NL forward matrix from target to focal plane
C
C	  convert the target variables  to correct units first

	  phit=phor/1000.0
	  thetat=thor/1000.0
	  yt = yor/100.0
	  deltat=dpp/100.0

C NT: complete transport matrix from target to focal plane
! x: 10 coefficients
! y: 12 coeff
! theta: 7
! phi: 10
! These coeffs are obtained from 97 data points fitting 
C Start of forward transport 
c============== x focal ========================
          xfoc =  coef_x(1)
     +           +coef_x(2)*deltat
     +           +coef_x(3)*deltat*deltat
     +           +coef_x(4)*phit
     +           +coef_x(5)*phit*phit
     +           +coef_x(6)*yt
     +           +coef_x(7)*thetat
     +           +coef_x(8)*thetat*thetat
     +           +coef_x(9)*thetat*yt
     +           +coef_x(10)*thetat*phit
     +           +coef_x(11)*Ep
     +           +coef_x(12)*Ep*thetat
     +           +coef_x(13)*Ep*phit
     +           +coef_x(14)*Ep*yt	  
c============== y focal ========================
          yfoc =  coef_y(1)
     +           +coef_y(2)*deltat
     +           +coef_y(3)*phit
     +           +coef_y(4)*phit*deltat
     +           +coef_y(5)*phit*phit
     +           +coef_y(6)*yt
     +           +coef_y(7)*yt*deltat
     +           +coef_y(8)*phit*yt
     +           +coef_y(9)*thetat
     +           +coef_y(10)*thetat*thetat
     +           +coef_y(11)*thetat*deltat
     +           +coef_y(12)*thetat*phit
     +           +coef_y(13)*Ep
     +           +coef_y(14)*Ep*thetat
     +           +coef_y(15)*Ep*phit   

c============== theta focal ========================
          thfoc = coef_th(1)
     +           +coef_th(2)*deltat
     +           +coef_th(3)*phit
     +           +coef_th(4)*yt
     +           +coef_th(5)*thetat
     +           +coef_th(6)*thetat*yt
     +           +coef_th(7)*thetat*phit  
     +           +coef_th(8)*Ep 
     +           +coef_th(9)*Ep*thetat 
     +           +coef_th(10)*Ep*phit                     
c==============phi focal ==========
          phfoc = coef_phi(1)
     +           +coef_phi(2)*deltat
     +           +coef_phi(3)*phit
     +           +coef_phi(4)*phit*deltat
     +           +coef_phi(5)*phit*phit
     +           +coef_phi(6)*yt
     +           +coef_phi(7)*yt*deltat
     +           +coef_phi(8)*thetat
     +           +coef_phi(9)*thetat*deltat
     +           +coef_phi(10)*thetat*phit
     +           +coef_phi(11)*Ep
     +           +coef_phi(12)*Ep*thetat

C End of forward transport matrix


c NT: Start the acceptance cut for first period
! This acceptance cut should work for dp = [-3,2]%
! and z=[-10,10]cm.
! How to determine this can be found below link:
!https://hallaweb.jlab.org/dvcslog/SaGDH/170809_154651/
!how_To_determine_acceptance_function.pdf
c
c This is a momentum dependant cut applied in y_rot (or yfoc) 
c and ph_rot space; a shape with 6 corners (6 sides)
c 
c Note that the same set of cuts should be applied to data: here is the 
c correlation between the varialbes here and the ones in the data root file
c
c deltat = R.gold.dp
c phfoc = R.tr.r_ph
c yfoc = R.tr.r_y
        y_corner(1) = -0.227484*deltat + 0.00965004
        y_corner(2) = -0.162115*deltat - 0.017178
        y_corner(3) =  0.203004*deltat - 0.0157547
        y_corner(4) =  0.524234*deltat + 0.0129056
        y_corner(5) = -0.260981*deltat + 0.024819
        y_corner(6) = -0.394433*deltat + 0.0205985
        
        ph_corner(1) = -0.120625*deltat + 0.0271181
        ph_corner(2) = -0.0885365*deltat + 0.012184
        ph_corner(3) = 0.209205*deltat - 0.0203365
        ph_corner(4) = 0.677652*deltat - 0.0290798
        ph_corner(5) = 0.105192*deltat + 0.0445469
        ph_corner(6) = -0.0579263*deltat + 0.0470789
c Get the values for the 6 edges of the acceptance 
c
c Edges 2 and 4 are vertical lines in ph vs. y plot; used for cuts in y
c The rests are horizontal lines, used for cuts in ph.
c

        do ii=1,6
           if((ii.eq.2).or.(ii.eq.4)) then
              slope = (y_corner(ii+1)-y_corner(ii))/
     +               (ph_corner(ii+1)-ph_corner(ii))
              intercept = y_corner(ii) - slope*ph_corner(ii)
              edge(ii) = slope*phfoc + intercept
           else
              if(ii.eq.6) then
                 slope = (ph_corner(1)-ph_corner(ii))/
     +                  (y_corner(1)-y_corner(ii))
              else
                 slope = (ph_corner(ii+1)-ph_corner(ii))/
     +                  (y_corner(ii+1)-y_corner(ii))
              endif
              intercept = ph_corner(ii) - slope*y_corner(ii)
              edge(ii) = slope*yfoc + intercept
           endif
        enddo
c Now apply the acceptance cut
        if((phfoc.gt.edge(1).and.phfoc.gt.edge(6)).or.
     + phfoc.lt.edge(3).or.phfoc.gt.edge(5).or.
     + yfoc.lt.edge(2).or.yfoc.gt.edge(4)) then
           goto 500
        endif
c

c NL,NT: End of First period cuts

!*********************** Reconstructed theta, phi, y ********************!
	rec_y = 0.00150805
     +	+0.38421*yfoc
     +	-0.39521*phfoc
     +	+0.00557369*xfoc
     +	+0.513705*thfoc
     +	+3.64968*yfoc*phfoc
     +	+7.22965*yfoc*yfoc
     +	-16.369*phfoc*thfoc

	rec_th = 0.0378725
     + -0.00376131*xfoc
     + +0.693092*yfoc
     + -0.269014*thfoc
     + +0.0412648*phfoc

	rec_ph = 0.00237242
     + -0.00698389*xfoc
     + +0.234642*yfoc
     + -0.0392215*thfoc
     + +0.40894*phfoc
     + -16.3558*yfoc*phfoc
     + +16.1314*yfoc*thfoc
     + -1.7407*yfoc*yfoc
     + +5.57496*phfoc*phfoc
*---------------------------------------------------------------------*
*** THis is for phase space cut to determine boundary for cacc, tacc
c$$$	ph_pos(1) = 0.012683!0.013065
c$$$	ph_pos(2) = 0.009747!0.008647
c$$$	ph_pos(3) = -0.004306!-0.005693
c$$$	ph_pos(4) = 0.001018!0.0
c$$$
c$$$	th_pos(1) = 0.045982!0.044567
c$$$	th_pos(2) = 0.031564!0.032660
c$$$	th_pos(3) = 0.035552!0.037486
c$$$	th_pos(4) = 0.047413!0.047139
c$$$
c$$$	do ll=1,4
c$$$	   if(ll.eq.4) then
c$$$	      ps_slope=(th_pos(1)-th_pos(ll))/(ph_pos(1)-ph_pos(ll))
c$$$	   else
c$$$	      ps_slope=(th_pos(ll+1)-th_pos(ll))/(ph_pos(ll+1)-ph_pos(ll))
c$$$	   endif
c$$$	   ps_intercept=th_pos(ll)-ps_slope*ph_pos(ll)
c$$$	   ps_edge(ll)=rec_ph*ps_slope + ps_intercept
c$$$	enddo
c$$$
c$$$	if(rec_th.lt.ps_edge(1).or.rec_th.lt.ps_edge(2).or.rec_th.gt.
c$$$     +   ps_edge(3).or.rec_th.gt.ps_edge(4)) then
c$$$	   goto 500
c$$$	endif
c convert the units back to the ones used in this program 
	  phit=phit*1000.0
	  thetat=thetat*1000.0
	  yt = yt*100.0
	  deltat=deltat*100.0
C Monte-Carlo trial is a 'GOOD' event

! Compute theoritical cross section and assymetries for each event
! kinematic Need to use dp/p before outgoing energy loss for Xsec computation
!	angstrag=th_spec+phor/1000.+(thor/1000.)**2. 

! Angle formula without approximation
	  angstrag=acos((1+abs(tan_HRS)*tan(phit/1000.0))
     +  /sqrt((1+tan_HRS**2)*(1+tan(phit/1000.0)**2
     +  +tan(thetat/1000.0)**2)))

! With approximation in calculate angle
C	  angstrag=acos((cos_HRS + abs(sin_HRS)*phit/1000.)
C     +           /sqrt(1+(phit/1000.)**2.+(thetat/1000.)**2.))

	  Eps=Ep*(1+dpnstrg/100.)
	  w=sqrt(2.808**2+2*2.808*(Etemp-Ep*(1+dpp/100.))
     +    -4*Etemp*Ep*(1+dpp/100.)*(sin(-1*angstrag/2.))**2.)-2.808

          QQ=4*Ei*Eps*(sin(-1*ang/2.))**2
! for c12
!	  mott=(6.*1./137.*cos(ang/2)/2./Ei/(sin(ang/2))**2.)**2.	

! for he3
	  mott=(1.*1./137.*cos(ang/2)/2./Ei/(sin(ang/2))**2.)**2.
!???          mottc=mottc*hbc2*1.e6 !convert to nanobarns per sr	
! 	  print*, 'Mott =', hbc2

c$$$	   if(rec_th.gt.0.048
c$$$     + .or.rec_ph.gt.0.016.or.rec_ph.lt.-0.015
c$$$     + )then
c$$$	     goto 500
c$$$	  endif




***********elastic**********************************
	if (choice.eq.1) then ! compute elastic
	   xs=cross_section(Ei,abs(ang))
	   if(rc.eq.1)then	!correct for vaccum polar and vertex 
	      xs=xs*(1+0.004647*(13./12.*log(QQ/(0.000511)**2.)-14./9.))
	   endif
	   asy=asym_calc3(Ei,ang,aspin)	

* test inf xs
	   if(xs.gt.HUGE(0.0)) then
	      print*,'************* INF cross section ********'
	      print*,'inf xs',ang*180./3.1415,
     + 'Ei=',Ei,'eps=',Eps,'Ep',Ep,'dp',dpp
	      goto 500
	   endif
  

	   if(thor.ge.xs_edge(1).and.thor.ge.xs_edge(2).and.thor.le.
     + xs_edge(3).and.thor.le.xs_edge(4).and.w.le.0.01.and.w
     + .ge.0.0) then
	      xsc=xs+xsc
              avasy=xs*asy+avasy
	      cacc=cacc+1
C	      print*,'xs= ',xsc,' asym= ',avasy
	   endif
	     
!	print*, 'XS=',xs,' Ei=',Ei,' Theta=',thetat
!	print*, 'asy=',asy
**********quasi-elastic*****************************	
	else 
	   if (choice.eq.2) then ! compute quasi-elastic
	      if((Ei/(1.+2.*Ei/2.808*(sin(ang/2.))**2.)-eps).le.0.0054)
     +              then 
				! below quasi-elastic threshold
		                ! no need to compute it		 
		 xs=0.0
		 asy=0.
		 !print*,'trop petit'
	      else ! compute cross section and asymetry

              ! initialyze values for interpolation.
		 rewind(24)	! come back at the beginning of the file
		 read(24,*)xom,sigma,qsq,atp,atlp,atpp,
     +           atpn,atlpp,atlpn,al,app
		 xomprev=0.	!we don't want it to be xom as initial value !
		 if (aspin.eq.0.) then ! spin orientation
		    asyprev=al	! longitudinal asymetries
		    xsprev=sigma
		 else
		    asyprev=app !transverse
		    xsprev=sigma
		 endif
	      ! end of initialization

		 rewind(24)	! come back at the beginning of the file
		 do i=1,500
		    read(24,*,end=997)xom,sigma,qsq,atp,atlp,atpp,
     +              atpn,atlpp,atlpn,al,app
		    xom=xom/1000. ! now xom in GeV (it was in MeV in the file)
	                    ! xom ordered in file by increasing values
		    if(xom-(Ei-Eps).ge.0)then ! take closest sup energy lost value 
				! and go out of the loop 
		 !print*, 'going out'
		       if (aspin.eq.0.) then ! spin orientation
			  asy=al ! longitudinal asymetries
			  xs=sigma
		       else
			  if(aspin.eq.90.) then
			     asy=app !transverse
			     xs=sigma
			  else 
			     print*,'no good angle for q-elastic (must be 0 or 90)'
			     goto 996
			  endif 
		       endif
		       goto 995
		    endif
		    xomprev=xom
		    if (aspin.eq.0.) then ! spin orientation
		       asyprev=al ! longitudinal asymetries
		       xsprev=sigma
		    else
		       asyprev=app !transverse
		       xsprev=sigma
		    endif
 997		    continue
		 enddo

 995		 continue	  
	    !print*, 'qtties:',asy,xs,xom
		 xs=(xs-xsprev)/(xom-xomprev)*(Ei-Eps-xomprev)+xsprev !interpolate
		 asy=(asy-asyprev)/(xom-xomprev)*(Ei-Eps-xomprev)+asyprev
	    !print*,'prev qtties:',asyprev,xsprev,xomprev
	      endif
 	    !print*, Ei-Eps,xom,asy,xs  
	   else   
*********phase space*********************************************
	      xs=0.
	      asy=0.	   
	   endif !if elastic
	endif	!if quasi-elastic

C Monte-Carlo is finished, copy result to Raw Wise Ntuple
	ztup(1)=zreact
	ztup(2)=xs
	ztup(3)=asy
	ztup(4)=w
	ztup(5)=mott
	ztup(6)=rec_y
	ztup(7)=rec_ph
	ztup(8)=rec_th
	ztup(9)=deltat
	ztup(10)=xfoc
	ztup(11)=yfoc
	ztup(12)=phfoc
	ztup(13)=thfoc
	ztup(14)=Eps
	ztup(15)=Epor
	ztup(16)=Etemp
	ztup(17)=ang
	ztup(18)=xbeam
	ztup(19)=ybeam
	ztup(20)=angstrag
	ztup(21)=z
	ztup(22)=wor
	ztup(23)=dpor
	ztup(24)=thor
	ztup(25)=phor
	ztup(26)=yt
	call hfn(1,ztup)
 500	continue

	enddo  !end of trials loop
	call endpaw

	close(20)
	close(21)
	close(22)
	close(23)
	close(24)
	close(25)
	close(30)
	close(31)
	close(32)
	close(33)
C output physics choice.
	if (choice.eq.0) then
	   print*, 'phase space'
	   print*,'acceptance in the chosen cut: ', cacc/tacc
	   print*,'total events in the full SA: ', tacc
	   print*,'Events that are within our SA: ', cacc
	   print*,'Available SA: ', cacc/tacc*full_sa,' Gev st.rad'
	endif
	if (choice.eq.1) then
	   print*, 'elastic'
	   print*,'Xsection= ',xsc/tacc,' nbarn'
           print*,'<asy>= ',avasy/xsc*100.,' %'
	   print*,'acceptance in the chosen cut: ', cacc/tacc
	   print*,'total events in the full SA: ', tacc
	   print*,'Events that are within our SA: ', cacc
	   print*,'Available SA: ', cacc/tacc*full_sa,' st.rad'
	   print*,'number of counts in acceptance ',xsc
	endif
	if (rc.eq.1) then
	      print*, 'radiative  correction on'
	      else
	      print*, 'radiative  correction off'	 
	endif
	if (choice.eq.2) then
	   print*, 'quasi-elastic'
	endif
	if (sieve.eq.1) print*, 'sieve slit in'
 996	continue
      	end





*******************************************************************************

        subroutine musc(rad_len,dth,dph,Ep,th_spec,z)
C+_____________________________________________________________________
!
! MUSC - Simulate multiple scattering of an electron in the HRS
!   spectrometer. This subroutine is a repaired and simplified version
!   of the old SLAC subroutine (used in MONTP) by the same name.
!   This subroutine will only calculate multiple scattering for electrons
!   in the HRS spectrometer.
!
! ASSUMPTIONS: DTH and DPH given in milli-radians, RAD_LEN in radiation
!   lengths. The formula used is due to Rossi and Greisen (See the book
!   by Segre, NUCLEI AND PARTICLES, 1982, p. 48.) The formula assumes a
!   gaussian distribution for the scattering angle. It is further assumed
!   that the angles DTH and DPH are the delta-theta and delta-phi angles
!   of an electron in the HRS, located at an angle with respect to
!   the beam direction that is large enough so that DTH-DPH space can
!   be approximated as a cartesian space.
!   Again, not the same subroutine as in uniform since we are in the HRS right
!   handed frame, not in the uniform left handed frame!
! from D. Potterveld 
C-_____________________________________________________________________

        implicit none
	real es,rndm,t1,pi,phi_sigma,phi_scat,theta_scat,u1,u2,rand
        real*4 rad_len,dth,dph,z
        real*4 Ep,th_spec
        pi = 3.141592654
        es = 13.6E-03 ! so energy has to be GeV

C Compute scattering angles, PHI_SCAT from a gaussian distribution,
C THETA_SCAT from uniform distribution.
        u1=rand()!rndm(1)
	u2=rand()!rndm(2)
	t1=sqrt(-2*log(u1))*cos(2*pi*u2)
C       theta_sigma = es*sqrt(rad_len)/Ep
C    CHANGED TO HAVE THE CORRECTION TERM ( FROM THE BOOKLET) 2/13/89
        phi_sigma = es*sqrt(rad_len)/Ep*z*(1.+0.038*LOG(RAD_LEN))
        phi_scat = t1*phi_sigma/2.  
        theta_scat = 2*pi*rand()!rndm(3)

C Compute new trajectory angles (units are mr)
        dth = dth + phi_scat*cos(theta_scat)*1000.
        dph = dph + phi_scat*sin(theta_scat)*1000.
        return
        end
*****************************************************************************

       subroutine project(xtr,ytr,z_drift,dth,dph)

C+___________________________________________________________________________
!
!   Calculate new HRS coordinates after drifting in a field
!   free region for a distance z_drift. It is assumed that all X, and Y
!   variables have the same (linear) dimensions, and that the angles DTH and
!   DPH are in units of milliradians and are close enough to zero that a
!   paraxial approximation may be made.
!
!   not the same subroutine as in uniform since we are in the HRS right
!   handed frame, not in the uniform left handed frame
C-___________________________________________________________________________

        implicit  none

        real*4 z_drift,xtr,ytr,dth,dph

        xtr = xtr + (dth/1000.) * z_drift
        ytr = ytr + (dph/1000.) * z_drift

        return
        end
*******************************************************************************
c
	subroutine sieve_slit(x,y,flag)
c
c added 050305,  cuts on 49 sieve holes (thin sieve)
c v.sulkosky 050405
c cut rays at the sieve slit to simulate sieve slit data
c
c x, xsieve
c y, ysieve
	real*4 x
	real*4 y
	real*4 flag
	real*4 ans
	integer i,j

c     
c     r-arm sieve slit:
c Old survey A870
c	real ysvpos(7)/0.0151401,0.0103649,0.0055897,0.0008145,
c     >  -0.0053069,-0.0114283,-0.0175497/ ! in m
c New survey A870r
	real xsvpos(7)/0.0380788,0.0247692,0.0114596,-0.00185,
     >  -0.0151596,-0.0284692,-0.0417788/ ! in m
	real ysvpos(7)/0.0155925,0.0108173,0.0060421,0.0012669,
     >  -0.0048545,-0.0109759,-0.0170973/ ! in m
c without offset
c	real xsvpos(7)/0.0399288,0.0266192,0.0133096,0.0,
c     >  -0.0133096,-0.0266192,-0.0399288/ ! in m
c	real ysvpos(7)/0.0143256,0.0095504,0.0047752,0.0,
c     >  -0.0061219,-0.0122438,-0.0183657/ ! in m

	real holesize(2)/0.000698,0.001346/ ! radius of hole in m
      
	i=1
	j=1
	flag=0.0
	do while ((i.lt.8) .and. (flag.eq.0.0)) ! sieve row
	   do while ((j.lt.8) .and. (flag.eq.0.0)) ! sieve column
				! check for large holes
	      if((i.eq.4 .and. j.eq.4) .or. (i.eq.6 .and. j.eq.5)) then
		 ans=sqrt((x-xsvpos(i))**2+(y-ysvpos(j))**2)
		 if(ans.gt.holesize(2)) then
		    flag = 0.0
		 else 
		    flag = 1.0
		 endif
	      else
		 ans=sqrt((x-xsvpos(i))**2+(y-ysvpos(j))**2)
		 if(ans.gt.holesize(1)) then
		    flag = 0.0
		 else 
		    flag = 1.0
c       write(*,*)' ans = ',ans, ' flag = ',flag
		 endif
	      endif		! sieve row loop
	      j=j+1
	   enddo		! sieve column loop
	   i=i+1
	   j=1
	enddo
	
	return
	end
c
c       end of sieve slit


*******************************************************************************
!!!!!!!!function to compute He 3 elastic cross section and assymetries
	real function cross_section(e,theta)

	implicit none

	real alpha,e,ef,fc,fm,Ge,Gm,mt,mott,nbarn,nu
	real q,qfm,recoil,s2,tau,t2,theta,w1,w2
	mt=2.809 !3He mass = 2.809 GeV/c**2
	alpha=7.3e-3 !	alpha = 1/137
	nbarn=0.389e6 ! barn: (1 GeV)**-2 = 0.389e-3 barn

	s2=sin(abs(theta)/2.)**2
	ef = e/(1+2.*e*(sin(abs(theta)/2.))**2/mt)
	q = (e-ef)*2.*mt
c	Calculate energy transfer
	nu = e - ef
	tau = nu/(2.*mt)

c	Calculate Mott cross section in barn/(GeV-sr)
	mott = ((alpha * cos(theta/2.) / (2. * e*s2) )**2)*nbarn

c	Calculate recoil factor
	recoil = ef/e

c	Estimate elastic form factors for 3He.
	qfm = q / 0.197328**2. ! h-bar c = 0.197328 GeV-fm
	Call Form2(qfm,fc,fm)
 	Ge = 2. * fc
	Gm = -6.3663 * fm
	t2 = tan(abs(theta)/2.)**2

c	Calculate cross sections at elastic peak (unpolarized)
	w1 = Gm**2 * tau
	w2 = (Ge**2 + tau * Gm**2)/(1 + tau)
	cross_section = recoil * mott * (w2 + 2 * w1 * t2)
	end
*******************************************************************************
      Real function asym_calc3(e,a,aspin)
	implicit none
	real a,ae,am,anum,anume,anumm,adif,aspin,astar
	real denom,e,ef,fc,fm,Ge,Gm,mt,nu,phistar,q,q3
	real q3sq,qfm,t,t2,tau,thetaq,vl,vt,vtlp,vtp

	mt=2.80793 !3He mass = 2.80793 GeV/c**2

	ef = e/(1+2.*e*(sin(abs(a)/2.))**2/mt)
	q = (e-ef)*2.*mt
	nu = e - ef
	tau = nu/(2.*mt)
	qfm = q / 0.197328**2. ! h-bar c = 0.197328 GeV-fm
	Call Form2(qfm,fc,fm)
 	Ge = 2. * fc
	Gm = -6.3663 * fm
c	Calculate the angle between the spin direction and
c	the momentum transfer.
	q3sq = q + nu**2
	q3 = sqrt(q3sq)
	thetaq = -asin(ef*sin(a)/q3)

	astar = abs( aspin*3.14159/180. - thetaq )
	if (astar.ge.0.0.and.astar.le.3.14159) then
	   phistar = 0.0
	else
	   phistar = 3.14159
	endif
c	astar = abs(astar)

C$$$	astar = abs(aspin*3.14159/180. - thetaq)
C$$$	adif = thetaq - aspin*3.14159/180.
C$$$	if ( astar .gt. 3.14159 ) astar = 6.28318 - astar
C$$$	if ( (adif .le. 0.) .and.
C$$$	1	(adif .gt. -3.14159) ) then
C$$$		phistar = 3.14159
C$$$	elseif (adif .gt. 3.14159) then
C$$$		phistar = 3.14159
C$$$	else
C$$$		phistar = 0
C$$$	endif
c
c to test if phistar has a wrong sign:
c this results an asym of +7.4% instead of -8.4%
c
c	phistar=3.14159-phistar
c
c
c now fix phistar and astar to E99117 elastic conditions
c
	phistar = 3.14159
	if (astar.ge.0.0.and.astar.le.3.14159/2) then
	   astar=3.14159-astar
	endif
	phistar = 0
	if (astar.ge.3.14159/2.and.astar.le.3.14159) then
	   astar=3.14159-astar
	endif

c to study the sys err due to the uncertainty in target spin direction: \pm 1 deg
c	astar=astar-3.1415926/360.

	t = tan(abs(a)/2.)
	t2 = t**2

	vl = (q/q3sq)**2
	vt = q/(2*q3sq) + t2
	vtp = t*sqrt(t2 + q/q3sq)
	vtlp = q*t/(1.4142*q3sq)
c
c this is to test if sign for astar is wrong, (like should be pi-aster
c  or sth. ) this results in an asym of
c  -7.4% instead of -8.4%.
c
c	anum = - 2.*tau*GM**2*cos(astar)*vtp + 2*sqrt(2*tau*
c
	anum = 2.*tau*GM**2*cos(astar)*vtp + 2*sqrt(2*tau*
     1         (1.+tau))*GM*GE*sin(astar)*cos(phistar)*vtlp
	denom = vl*(1.+tau)*GE**2 + vt*2*tau*GM**2


	anumm = 2.*tau*GM**2*cos(astar)*vtp 
	anume=2*sqrt(2*tau*(1.+tau))*GM*GE*sin(astar)*cos(phistar)*vtlp
	am = -anumm/denom
	ae = -anume/denom

c	Calculate Asymmetry
	asym_calc3 = -(anum/denom)
	end
*********************************************************************
	subroutine form2(q2,fc,fm)
c
c	Calculates the normalized charge and magnetic form factors
c	for 3He using the parametrization of Amround et al (1994)
c
	real q2,fc,fm

	fc = formc_1994(q2)
	fm = formm_1994(q2)
	return
	end
******for 1994 data****************
      FUNCTION FORMC_1994(Q2)
c     Function to compute a SOG form factor
c     Q2 : momentum transfer squared (fm-2)
C     NR : number of Gaussians
C     GA : Gaussians rms (usually 0.8)
C     RR : Gaussians central positions
C     QQ : Gaussians amplitudes
      real RR(12),QQ(12)
      real  Q2,GA,A,B
      DATA NR/12/,GA/0.65/
      DATA RR/0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.0,4.6,5.2/
      DATA QQ/.027614,.170847,.219805,.170486,.134453,.100953,
     +        .074310,.053970,.023689,.017502,.002034,.004338/
      Q=SQRT(Q2)
      G2=GA*GA
      A=EXP(-Q2*G2/4.)
      S=0.
      DO I=1,NR
          B=2.*RR(I)**2/G2
          QR=Q*RR(I)
          IF (QR.EQ.0.) THEN
              SS=1.+B
          ELSE
              SS=COS(QR)+B*SIN(QR)/QR
          END IF
          SS=QQ(I)/(1.+B)*SS
          S=S+SS
      END DO
      FORMC_1994=A*S
      RETURN
      END
******************************
      FUNCTION FORMM_1994(Q2)
c     Function to compute a SOG form factor
c     Q2 : momentum transfer squared (fm-2)
c     NR : number of Gaussians
c     GA : Gaussians rms (usually 0.8)
C     RR : Gaussians central positions
C     QQ : Gaussians amplitudes
      real RR(12),QQ(12)
      real  Q2,GA,A,B
      DATA NR/12/,GA/0.65/
      DATA RR/0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.0,4.6,5.2/
      DATA QQ/.059785,.138368,.281326,.000037,.289808,.019056,
     +        .114825,.042296,.028345,.018312,.007843,.000000/
      Q=SQRT(Q2)
      G2=GA*GA
      A=EXP(-Q2*G2/4.)
      S=0.
      DO I=1,NR
          B=2.*RR(I)**2/G2
          QR=Q*RR(I)
          IF (QR.EQ.0.) THEN
              SS=1.+B
          ELSE
              SS=COS(QR)+B*SIN(QR)/QR
          END IF
          SS=QQ(I)/(1.+B)*SS
          S=S+SS
      END DO
      FORMM_1994=A*S
      RETURN
      END

****************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to get the energy loss by external with the right proba due 
!Ionisation and bremstralung. Formula are from Xiaodong Jiang Thesis
! and Claude Marchand Thesis
* CEBAF, A. D 10/25/1999
!Note that the function does not depend much on a and z. Their purpose is 
! just to get a more accurate b.
! input and output Units are GeV and radiation length
      real function incl(E0,z,t)
 
      implicit none
 
      real b,E0,rndm,t,z,rand

      b=4./3.*(1.+(z+1.)/9./(log(1194.*z**(-2./3.))+z*log(
     +184.15*z**(-1./3.))))
      incl=(E0)*(rndm(1)*0.999)**(1./(b*t)) ! Bremstralung

 	if (.not.((incl.ge.-100.).and.(incl.le.10**20.))) then
 	   print*,'Arg, Something wrong with external RC for this event'
           print*,'external:',E0*1000.,z,t,incl
	   incl=10**20.
 	endif		
      end
!routine to compute the probability to lose energy by internal
! Formula are from C Marchand thesis
* CEBAF, A. D 10/25/1999
* Input Units are GeV  and radian
      real function intrad(E0,thp)

      implicit none

      real alpha,E,E0,me,Mhe3,nu,pi,rndm,qsq,thp,rand

      alpha=1/137.
      pi=3.141592654
      Mhe3=2.808 ! masse nucleus in GeV   
      me=0.000511 ! mass e- in GeV   
 
      E=E0/(1.+2.*E0/Mhe3*(sin(thp/2.))**2.)
      qsq=2.*E0*E*(1.-cos(thp))
      nu=2.*alpha/pi*(log(qsq/me**2.)-1.) 
      nu=nu/2.
 	intrad=E*(rndm(1)*0.999)**(1./(nu)) 
	if (.not.((intrad.ge.-100000.).and.(intrad.le.10000000.))) then
	   print*,'Arg, Something wrong with internal RC for this event'
      print*,'internal:',' thp:',thp,' Q^2:',qsq,' rcw:',intrad
	endif	
      end

!function to get energy loss by ionization
!Formula are from the Leo
* CEBAF, A. D 11/27/1999
! input and output Units are GeV and radiation length
      real function ioni(E0,a,z,xd)
 
      implicit none
 
      real a,betas,betasmo,d0,de,E0,ksi
      real lneps,pi,phil,rndm,u,xd,z,rand
      integer j
      parameter (pi=3.1415926)
      common/landau/ld(2,68)
      real ld !landau integral ld(1,i)=lambda, ld(2,i)=int(phi(lambda))
      E0=E0*1000.!convert E0 in MeV
      betas=1-(0.511/E0)**2.
      betasmo=(0.511/E0)**2. ! beta^2 of the electron
      ksi=0.154/betas*z/a*xd ! warning, ksi in MeV
      lneps=log((betasmo)*(0.0135*z)**2./2./0.511/betas)+betas
      d0=ksi*(log(ksi)-lneps+0.37)  ! ionisation energy loss (MeV) most probable      
      u=rndm(1)
      do j=1,68 ! look at the inverse of the integrated landau tail
         if(u.le.ld(2,j)) then ! interpolate
            phil=(ld(1,j)-ld(1,j-1))/(ld(2,j)-ld(2,j-1))
            phil=phil*(u-ld(2,j-1))+ld(1,j-1)
         goto 2
         endif
      enddo
 2    continue
      de=d0+0.06+phil*ksi!express the Landau shape in term of energy loss (MeV)
      E0=E0/1000.
      ioni=de/1000.
     
      end 
!function to get most probable energy loss by ionization
* CEBAF, A. D 11/27/1999
! input and output Units are GeV and radiation length
!      real function deltazero(E0,a,z,xd)

!      implicit none
 
!      real a,betas,betasmo,d0,E0,ksi
!      real lneps,xd,z
 
!      E0=E0*1000.!convert E0 in MeV
!      betas=1-(0.511/E0)**2.
!      betasmo=(0.511/E0)**2. ! beta^2 of the electron
!      ksi=0.154/betas*z/a*xd	! warning, ksi in MeV
!      lneps=log((betasmo)*(0.0135*z)**2./2./0.511/betas)+betas
!      d0=ksi*(log(ksi)-lneps+0.37)  ! ionisation energy loss (MeV) most probable 
!      E0=E0/1000. 	
!      deltazero=d0/1000. !(GeV)
    
!      end 

!! subroutine to fill the common block with landau dist
      subroutine stuffit
	implicit none
	integer i
      common/landau/ld(2,68)
      real lambda(68),intphi(68),ld
      data lambda/-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-.5,.0 ,.5 ,1.0,1.5 ,2. 
     +,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5 ,9.0 ,9.5,
     +10.0,10.5,11.0 ,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,
     +16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5,21.,21.5,22.,22.5,23.,
     +23.5,24.0,24.5,25.0,25.5,26.,26.5,27.,27.5,28.,28.5,29.,29.5,30. /
      data intphi/.0,5.00000E-04,6.50000E-03,2.90000E-02,8.15000E-02,
     +.1515,.239,.329,.414,.4875,.55,.6025,.6475,.6875001,.7195001,
     +.7455001,.7670001,.7850001,.8005001,.8145001,.8270001,.8385001,
     +.8482501,.8570001,.8650001,.8722501,.8787501,.8847501,.8902501,
     +.8952501,.8998752,.9041252,.9080002,.9115002,.9147502,.9177502,
     +.9205002,.9230002,.9253752,.9276252,.9297502,.9317502,.9336252,
     +.9353752,.9370502,.9386502,.9402002,.9417002,.9430752,.9443253,
     +.9454502,.9464502,.9474002,.9482752,.9491102,.9498602,.9505252,
     +.9511502,.9517627,.9523627,.9529527,.9535328,.9541028,.9546627,
     +.9552028,.9557328,.9562528,.9567627 /
	do i=1,68
	   ld(1,i)=lambda(i)
	   ld(2,i)=intphi(i)
	enddo
      end
!!subroutines for row wise ntuples
      subroutine initpaw
      common/pawc/blanc(150000)
      character*8 tags(26)
*      data tags/'zlab','deltat','thetat','phit','yt','xs','asy'
*     +,'wmm','mott'/
      data tags/'zreact','xs','asy','wmm','mott','rec_y','rec_ph'
     +,'rec_th','deltat','xfoc','yfoc','phfoc','thfoc'
     +,'Eps','Epor','Etemp','ang','xbeam','ybeam','angstrag'
     +,'z','wor','dpor','thor','phor','yt'/
*	data tags/'zlab','xs','asy','wmm','mott','yor','phor'
*     +,'thor','dpp','xfoc','yfoc','phfoc','thfoc'/
      call hlimit (150000)
      call hlimit (150000)
      call hropen (55,'bidon','mchc12ntp.hbook','n',8190,irc)
      call hbookn(1,'truc',26,'bidon',1000,tags)
      end
c
      subroutine endpaw
      call hrout(0,icycle,' ')
      call hrend('bidon')
      end






