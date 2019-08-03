cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     File: simplejet.f -- simplified jet model works best with synchrotron only (set comsw = 1.d0 
c                          for compton) there is no accretion disk or BB component
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine xrbjet(ear,ne,param,ifl,photar,photeng)
      include 'param_dbl.inc'
      integer ifl,ne,nebin
      double precision  ear(0:ne),param(15),photar(ne),photeng(ne)
      parameter(nebin=10000)
      double precision fsc,asum,esum
      double precision tbb2,bbf2,bbf1,reff,reff2
c      double precision mdot
      double precision toten,bbdil,tst,uphdil2
      double precision bbirad,bbnrm,bbtrm

      integer nz,nelec,nsyn,ncom,nzdum,ncarr,nsarr,njet
      double precision ebin(nebin),dele(nebin),ntot0,ntot,mbh,inclin
      double precision jetrat,mxsw,plotsw,plfrac,disksw
      double precision zmin,zmax,zinc,zcut,rnumin,rnumax,xemin
      DOUBLE PRECISION xemax
csm   next are parameters that used to be changeable that are now fixed
      double precision outfac,zfrac,dopsw,comsw,bbsw
      parameter(zinc=.1D0,nzdum=200)
      parameter(xemin=12.D0,xemax=23.D0,njet=2)
      parameter(nsyn=100,ncom=60)
      parameter(ncarr=nzdum*ncom*njet,nsarr=nzdum*njet*nsyn)
      parameter(nelec=200,zfrac=200.d0)
      parameter(dopsw=1.d0,bbsw=1.d0)
      integer i,j,k,nw,m,check
      DOUBLE PRECISION alpha,zsh,r0,eltemp,eddlum,rin,r_g,betas0
      double precision gam0,rvel0,vel0,b_en,b0,h0,gamax0
      DOUBLE PRECISION ajet,sum1,sum2,tm1,tm2
      DOUBLE PRECISION endens,betat,einc,game,rvel,gamax,r,gamv2,gamv
      double precision beta,gshift,renorm,totlos,tmlos
      DOUBLE PRECISION oldnum,gshock,ub,ucom
      DOUBLE PRECISION uphdil,accon,syncon,comcon,escom,qutrmb,qutrmc
      double precision emax,fplot(nebin),complot(nebin),presyn(nebin)
      double precision postsyn(nebin),bbplot(nebin),visco
      DOUBLE PRECISION emin,cnorm,rnuinc,absd,phoden
      double precision energ(nebin)
      DOUBLE PRECISION bbcon,area,vol,xeinc,com
      DOUBLE PRECISION sflx,cflx,frq,bbflx,hbb,rout,gemax,bemax
      DOUBLE PRECISION z,delz,rdlgen(nelec),rden(nelec),rdedn(nelec)
      DOUBLE PRECISION rdend(nelec),lelec,ledens,eled(nelec)
      DOUBLE PRECISION etemp(nelec),dtemp(nelec)
      double precision ytb(nelec),ytc(nelec),ytd(nelec)
      double precision yelb,yelc,yeld
      DOUBLE PRECISION tstrm(nelec),nurad(nsyn)
      double precision elen,drtrm,yderb,yderc,yderd
      DOUBLE PRECISION ephxr(ncom),comspc(njet,nzdum,ncom),nutot(nebin)
      DOUBLE PRECISION snu(nsyn),sdump(nsyn),cnu(ncom),cdump(ncom)
      DOUBLE PRECISION ephot,synabs(njet,nzdum,nsyn)
      double precision nphot(nsyn),phodis,yphot
      DOUBLE PRECISION nusyn(njet,nzdum,nsyn),nucom(njet,nzdum,ncom)
      double precision ysyb(nsyn),ysyc(nsyn),ysyd(nsyn)
      double precision ycob(nsyn),ycoc(nsyn),ycod(nsyn)
      DOUBLE PRECISION bfield,arg,val,elenmn,elenmx,ephpass,ephmax
      double precision ephmin
      DOUBLE PRECISION hratio,eddrat,tin,thrlum
      DOUBLE PRECISION phofrq(nsyn),phoint(nsyn),ebot,edtrm
      DOUBLE PRECISION pltrm,k2,enorm,bete,reled,elmin,elmax
      DOUBLE PRECISION thmfrac,dopfac(njet,nzdum)
      double precision dkpc,dist
      DOUBLE PRECISION shelen(nelec),shden(nelec),zpass,shocknd
      DOUBLE PRECISION shockr,shockgb,ulim
      DOUBLE PRECISION blim,bbnumax
      double precision nubb(nsyn)
      double precision ynub(nsyn),ynuc(nsyn),ynud(nsyn)
      DOUBLE PRECISION nubot,nutop,bbint
      DOUBLE PRECISION theff,tbbeff,equip,photmx,peaksw
      double precision jetu(nsyn),disku(nsyn),nprot
      double precision avphot,avelec,typmax
      double precision gamfac,gamin
      double precision cb,cc,cd
      double precision synb,sync,synd
      logical initialized
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cDM Playing with parallelization
c      INTEGER NTHREADS, OMP_GET_NUM_THREADS, CHUNK
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Playing with irradiation by the jet.
c      double precision irradsw
c      integer nrad,nradbig,nbig,nface
c      parameter(nradbig=100,nface=1,nbig=nebin*nradbig)
c      double precision rmin,rmax,rinc,direct,scatt(nsyn)
c      double precision radius(nradbig),rbinu(nradbig),rbind(nradbig)
c      double precision doppler,thetd(nradbig),diska(nradbig),annuarea
c      double precision nusynd(nzdum,nradbig,nsyn)
c      double precision nucomd(nzdum,nradbig,ncom)
c      double precision syndsk(nzdum,nradbig,nsyn)
c      double precision comdsk(nzdum,nradbig,ncom)
c      double precision nuplt(nebin),refemis(nradbig,nebin)
c      double precision total(nebin),T_irrad(nradbig)
c      double precision disk_albedo,tout
c      data refemis/nbig*0.d0/,total/nebin*0.d0/,T_irrad/nradbig*0.d0/
c      data bbirad/0.d0/, totlos/0.d0/, gshock/0.d0/, cnorm/0.d0/
c      data bbcon/0.d0/, shocknd/0.d0/, shockr/0.d0/, shockgb/0.d0/
c      data reff/0.d0/,reff2/0.d0/,nrad/0/, bbf1/0.0d0/
c      data disk_albedo/0.0/, tout/0.0/
c      external mcdobs1,mcdjet1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      external comint
      common /syn/arg(47),val(47),synb(47),sync(47),synd(47)
      common /eldis1/lelec(nelec),ledens(nelec),yelb(nelec),
     *yelc(nelec),yeld(nelec)
      common /eldis2/elenmn,elenmx,bfield,dist
      common /derivs/elen(nelec),drtrm(nelec),yderb(nelec),yderc(nelec),
     &yderd(nelec)

      common /com/ephot(nsyn),phodis(nsyn),yphot(nsyn)
      common /com2/ephpass,ephmax,ephmin
      common /jpars/visco,mbh,zmin
      common /bbpars/rin,rout,hbb,tin,zpass
      common /erdpar/frq,inclin
      common /jetg/gamv
      common /comsplinpar/ cb(nsyn),cc(nsyn),cd(nsyn)

      data   initialized/.FALSE./
      data   etemp,dtemp/nelec*0.D0,nelec*0.D0/      
      data   synabs,comspc/nsarr*0.D0,ncarr*0.D0/
      data arg/0.0001d0,0.0002d0,0.0005d0, 0.001d0, 0.002d0, 0.005d0,  
     *         0.01d0,  0.03d0,  0.05d0,   0.07d0,  0.1d0,   0.2d0,   
     *         0.3d0,   0.4d0,   0.5d0,    0.6d0,   0.7d0,   0.8d0,   
     *         0.9d0,   1.d0,    1.5d0,    2.d0,    2.5d0,   3.d0,   
     *         3.5d0,   4.d0,    4.5d0,    5.d0,    5.5d0,   6.d0,   
     *         6.5d0,   7.d0,    7.5d0,    8.d0,    8.5d0,   9.d0,   
     *         9.5d0,   10.d0,   12.d0,    14.d0,   16.d0,   18.d0,
     *         20.d0,   25.d0,   30.d0,    40.d0,   50.d0 /
      data val/ -1.002d+00, -9.031d-01, -7.696d-01, -6.716d-01,
     *          -5.686d-01, -4.461d-01, -3.516d-01, -2.125d-01,
     *          -1.537d-01, -1.192d-01, -8.725d-02, -4.383d-02,
     *          -3.716d-02, -4.528d-02, -5.948d-02, -7.988d-02,
     *          -1.035d-01, -1.296d-01, -1.586d-01, -1.838d-01,
     *          -3.507d-01, -5.214d-01, -6.990d-01, -8.861d-01,
     *          -1.073d+00, -1.267d+00, -1.470d+00, -1.670d+00,
     *          -1.870d+00, -2.073d+00, -2.279d+00, -2.483d+00,
     *          -2.686d+00, -2.893d+00, -3.097d+00, -3.303d+00,
     *          -3.510d+00, -3.717d+00, -4.550d+00, -5.388d+00,
     *          -6.230d+00, -7.075d+00, -7.921d+00, -1.005d+01,
     *          -1.218d+01, -1.646d+01, -2.076d+01/


       mbh=param(1)
       r_g=gconst*mbh*msun/clite**2
       eddlum=1.25D38*mbh
       jetrat=param(2)
       jetrat=eddlum*jetrat
       alpha=param(3)
       zsh=param(4)*r_g
       r0=param(5)*r_g
       hratio=param(6)
       h0=hratio*r0
       inclin=pi*param(7)/180.0d0
       eltemp=param(8)
       plfrac=param(9)
       dkpc=param(10)
       dist=dkpc*1.D3*pc
       comsw=0.d0
       plotsw=param(12)
       fsc=param(13)
       zmax=param(14)
       equip=param(15)

       print*,'Running jet model ...'
       check=0

       disksw=1.d0
       if (plotsw.eq.1.d0) then
          call getlun(ilun12)
          open(ilun12,file='total.dat' )
          call getlun(ilun13)
          open(ilun13,file='com.dat')
          call getlun(ilun14)
          open(ilun14,file='presyn.dat')
          call getlun(ilun16)
          open(ilun16,file='postsyn.dat')
       endif

       do i=1,nebin
          fplot(i) = 0.0d0
          complot(i)=0.d0
          presyn(i)=0.d0
          postsyn(i)=0.d0
       enddo

c ** define energy array from ear array, and then convert into frequencies
c ** to be compatible with old routine.  ear, ebin, dele in keV

       do i=1,ne
          ebin(i)=ear(i-1)+ (ear(i)-ear(i-1))/2.d0
          dele(i)=ear(i)-ear(i-1)
          nutot(i)=log10(ebin(i)/hkev)
       end do    

       rnumin=log10(0.5d0*ear(0)/hkev)
       rnumax=log10(ear(ne)/hkev)
       rnuinc=(rnumax-rnumin)/nsyn
      do i=1,nsyn
        nubb(i)=rnumin+(i-.5D0)*rnuinc
        nurad(i)=10.D0**nubb(i)
        energ(i)=herg*nurad(i)
        ephot(i)=log10(energ(i))
      end do 

      xeinc=(xemax-xemin)/ncom
      do i=1,ncom
         ephxr(i)=10.D0**(xemin+(i-.5D0)*xeinc)
      end do

      zmin=int(log10(r0*0.3d0)/zinc)*zinc      
      nz=(zmax-zmin)/zinc +1
      if (nz.gt.nzdum) then
         print*,nz,'nz greater than array size!'
         stop
      end if      

c     initialization of the synchrotron interpolation
      call initspline(47,arg,val,synb,sync,synd)

      thmfrac=1.D0-plfrac


      betas0=sqrt((gad4_3-1.D0)/(gad4_3+1.D0))
      gam0=1.D0/sqrt(1.D0-betas0**2)
      rvel0=gam0*betas0
      vel0=clite*rvel0
      zcut=zfrac*r0

c    see pg 291, iv (assump 2, 4 or 5) 


      emin=2.23d0*kboltz*eltemp
      gamin=emin/emerg
      emax=gamfac*gamin*emerg

      nprot=(jetrat/4.D0)/(rvel0*clite*pmgm*clite**2*pi*r0**2)

      if (equip.eq.1.d0) then
         ntot0=nprot
         b_en=ntot0*3.28173d0*kboltz*eltemp
      else if (equip.eq.0.d0) then
         ntot0=nprot
         b_en=ntot0*(pmgm*clite**2-3.28173d0*kboltz*eltemp)
      else
         ntot0=nprot
         b_en=equip*ntot0*3.28173d0*kboltz*eltemp
      endif      

      b0=sqrt(8.D0*pi*b_en)
      

c  Freqs in Hz, fluxes in mJy
c  set gamax0 to unity so gamax/gamax0 gives correct mach number factor
c  for energy shift

      gamax0=1.D0

      endens=3.28173D0*kboltz*eltemp*ntot0
      betat=emerg/(kboltz*eltemp)
      enorm=ntot0*betat/k2(betat)
      elmin=emerg
      elmax=50.D0*kboltz*eltemp
      einc=log10(elmax/elmin)/nelec
      do i=1,nelec
         rdlgen(i)=log10(elmin)+(i-.5D0)*einc
         rden(i)=10.D0**rdlgen(i)
         game=rden(i)/emerg
         if (game.le.1.D3)then
            bete=sqrt(game**2-1.D0)/game
         else
            bete=1.D0-1.D0/(2.D0*game**2)-1.D0/(8.D0*game**4)
         endif
         reled=enorm*(bete/emerg)**3*rden(i)**2*
     *        exp(-betat*sqrt(1.D0+(bete*rden(i)/emerg)**2))
         rdend(i)=reled*rden(i)
         rdedn(i)=log10(reled)
      end do

c      print*, 'EQUIP in Noz', b_en/endens
      if (.not. initialized) call initjet

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  START THE BIG LOOP OVER Z HERE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cDM For reffrac_base maybe do k=1,5 or 10???
      do k=1,nz
         z=zmin+zinc*(k-1)
         delz=10.D0**z-10**(z-zinc)
         z=10.D0**z
         zpass=z
         call jetpars(z,b0,ntot0,gamax0,r0,h0,rvel,bfield,ntot,gamax,r)

         area=pi*r**2
         vol=delz*area

         gamv2=1.D0+rvel**2
         gamv=sqrt(gamv2)

         if (gamv2.gt.1.D5)then
            beta=1.D0-1.D0/(2.D0*gamv2)-1.D0/(8.D0*gamv2**2)
         else
            beta=sqrt(gamv2-1.D0)/gamv
         endif


         if (dopsw.eq.1)then
            do l=1,njet
               dopfac(l,k)=1.D0/(gamv*(1.D0-beta*cos(inclin)*(-1.d0)**
     *                   (l-1)))      
            end do
         else
            do l=1,njet
               dopfac(l,k)=1.D0
            end do
         endif
         
         gshift=gamax/gamax0
         nw=0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if past the shock, can read stored distribution directly
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (z .gt. zsh .and. check.eq.2) goto 21
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  take read in distrib, Shift the energies down by mach**-1/3, 
c  keep track of lost ones by renormalizing first and then losing
c  the first bins...and renorm'ing

         do i=1,nelec
            lelec(i)=rdlgen(i)
            elen(i)=10.D0**lelec(i)
            ledens(i)=rdedn(i)
            ledens(i)=log10(ntot/ntot0*10.D0**ledens(i))
         end do

         do i=1,nelec
            etemp(i)=log10(gshift*elen(i))
            if (gshift*elen(i)/emerg.le.1.D0)then
               nw=nw+1
            endif
            dtemp(i)=ledens(i)
         end do

         sum1=0.D0
         do i=1,nelec-1
            tm1=(10.D0**etemp(i+1)-10.D0**etemp(i))*.5D0*((10.D0**dtemp
     *          (i))+(10.D0**dtemp(i+1)))
            sum1=sum1+tm1
         end do

cjw         renorm=ntot/sum1
         renorm=log10(ntot/sum1)
         do i=1,nelec
cjw            dtemp(i)=log10(renorm*10.D0**dtemp(i))
            dtemp(i)=renorm+dtemp(i)
         end do

c   spline_dbl and resize the array into the old lelec,ledens if you
c   lose particles at low energy
         if (nw.gt.0) then
            if (nw.eq.1) then
               totlos=(10.D0**etemp(2)-10.D0**etemp(1))*10.D0**dtemp(1)
            else if (nw.gt.1) then
               totlos=0.D0
               do i=1,nw-1
                  tmlos=(10.D0**etemp(i+1)-10.D0**etemp(i))*.5D0*(10.D0
     *                  **dtemp(i)+10.D0**dtemp(i+1))
                  totlos=totlos+tmlos
               end do
            endif
            oldnum=ntot
            ntot=ntot-totlos
         
cjw            renorm=ntot/oldnum
            renorm=log10(ntot/oldnum)
            do i=1,nelec
cjw               dtemp(i)=log10(renorm*10.D0**dtemp(i))
               dtemp(i)=renorm+dtemp(i)
            end do

            call initspline(nelec,etemp,dtemp,ytb,ytc,ytd)

            einc=(etemp(nelec)-etemp(nw+1))/nelec
            do i=1,nelec
               lelec(i)=etemp(nw+1)+(i-.5D0)*einc

               call interspline(nelec,etemp,dtemp,ytb,ytc,ytd,lelec(i),
     *               ledens(i))

               elen(i)=10.D0**lelec(i)
               eled(i)=elen(i)*10.D0**dtemp(i)
            end do
         else
            do i=1,nelec
               lelec(i)=etemp(i)
               elen(i)=10.D0**lelec(i)
               eled(i)=elen(i)*10.D0**dtemp(i)
               ledens(i)=dtemp(i)
            end do
         endif
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if you're only in nozzle...go directly to S & C, otherwise be accel'd
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (z.lt.zsh) goto 23
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  start commands for onetime thing
c  skip over 21 commands to energy shift and synchro calc
c         print*, 'THIS SHOULD ONLY APPEAR ONCE FOR Z,seg=', z,k

         gshock=gshift

         ajet=pi*r**2

c to get idea of SSC rate...use flux from previous component at distance
c delz...will be overestimate

c         print*,'ACCELPARS',r,z,dopfac(1,k),dopfac(2,k),bfield,ntot
         ub=bfield**2/(8.D0*pi)

         sum1=0
         do m=1,njet
            do i=1,nsyn
               phofrq(i)=10.D0**nusyn(m,k-1,i)/dopfac(m,k-1)
c 3 powers of dopfac...one was inside synint for angle aberration
               phoint(i)=mjy*10.D0**synabs(m,k-1,i)/dopfac(m,k-1)**3

c     handle differently:  not flux but assuming that the synchrotron
c                          was produced in this segment...overestimate
c                          but hopefully not too much because r,delz
c                          also increased
c old:              phoint(i)=phoint(i)*4.*dist**2/(clite*delz**2)
c new:  without cancelled units:  flux * D^2/ area of slab * time to cross
c        phoint(i)*4*pi*D^2 * delz/c /(pi r^2 delz)
               phoint(i)=phoint(i)*4.D0*dist**2/(r**2*clite)
            end do
            
c rough integration to get approx energy density 
            do i=1,nsyn-1
               tm1=(phofrq(i+1)-phofrq(i))*.5D0*(phoint(i)+phoint(i+1))
               sum1=sum1+tm1
            end do
         end do            

         ucom=sum1


c  see pg 274 in ntbk III to explain rates 
c see pg 32 ntbk VI for the syncon=escom version


c         accon=3.D0/4.D0*((vjet/clite)**2)*clite*charg*bfield/fsc
         accon=3.d0/4.d0*fsc*clite*charg*bfield
         syncon=4.D0/3.D0*sigtom*ub/(emgm**2*clite**3)
         comcon=syncon*ucom/ub
         escom=beta*clite/z
         qutrmB=escom/(syncon+comcon)
         qutrmC=accon/(syncon+comcon)
         emax=(-qutrmB+sqrt(qutrmB**2+4.D0*qutrmC))/2.D0
         gemax=emax/emerg
         if (gemax.ge.100.D0) then
            bemax=1.D0-1.D0/(2.D0*gemax**2)-1.D0/(8.D0*gemax**4)
         else
            bemax=sqrt(gemax**2-1.D0)/gemax
         endif

c  from this emax can generate powerlaw between min and emax
c  ...use emin normalized to scaled thermal peak at shock, and after shock
c  then scale the two down together (see pg 14, ntbk 4)

         emin=2.23D0*kboltz*eltemp*gshift 

c   create powerlaw of normalization plfrac*ntot

         cnorm=plfrac*ntot*(1.D0-alpha)/(emax**(1.D0-alpha)-emin**
     *         (1.D0-alpha))

c   make energy array between lelec(1) and log10(emax), add power
c   law and particle distribution together in appropriate ratios,
c   then write into array  and after shock, read it in
c   and shift appropriately according to shock

c   in theory this should be correctly normalized, but it's not.  renorm to 
c   ntot at end to make exacter

         call initspline(nelec,lelec,ledens,ytb,ytc,ytd)
         
         ebot=lelec(1)
         einc=(log10(emax)-ebot)/nelec
         do i=1,nelec
            etemp(i)=ebot+(i-.5D0)*einc
            if (etemp(i).le.lelec(nelec)) then

               call interspline(nelec,lelec,ledens,ytb,ytc,ytd,etemp(i),
     *                          edtrm)

               edtrm=thmfrac*10.D0**edtrm
            else
               edtrm=0.D0
            endif
            if (etemp(i).ge.log10(emin)) then
               pltrm=cnorm*(10.D0**etemp(i))**(-alpha)
            else
               pltrm=0.D0
            endif
            dtemp(i)=log10(edtrm+pltrm)
            lelec(i)=etemp(i)
            elen(i)=10.D0**lelec(i)
         enddo

         sum1=0.D0
         do i=1,nelec-1
            tm1=(elen(i+1)-elen(i))*.5D0*(10.D0**dtemp(i)+10.D0**dtemp
     *          (i+1))
            sum1=sum1+tm1
         end do
         renorm=ntot/sum1


         shocknd=ntot
         shockr=r
         shockgb=rvel
         
         do i=1,nelec
            ledens(i)=log10(renorm*10.D0**dtemp(i))
            eled(i)=elen(i)*(10.D0**ledens(i))
            shelen(i)=etemp(i)
            shden(i)=dtemp(i)
         end do

         sum1=0.D0
         sum2=0.D0
         do i=1,nelec-1
            tm1=(elen(i+1)-elen(i))*.5D0*(eled(i)+eled(i+1))
            tm2=(elen(i+1)-elen(i))*.5D0*(10.D0**ledens(i)+10.D0**
     *          ledens(i+1))
            sum1=sum1+tm1
            sum2=sum2+tm2
         end do
c         print*,'EQUIP at Shock',bfield**2/(8.D0*pi*sum1)
c         print*,'NTOT AFTER SHOCK',sum2
         check=2
        

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c    skip reading this file over for shock r,z
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         goto 23
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


 21      do i=1,nelec
            etemp(i)=shelen(i)
            dtemp(i)=shden(i)
         end do

c    shift number density down from shock, not from ntot0 anymore

         ntot=shocknd*(shockr/r)**2*(shockgb/rvel)
c         print*,'compare',shocknd,ntot

c  Now shift down in energies and lose the lost guys
c  NOTE: the shift is from the shock=gshift/gshock

         nw=0
         do i=1,nelec
            etemp(i)=log10(gshift/gshock*10.D0**etemp(i))
            if (10.D0**etemp(i)/emerg.le.1.D0) nw=nw+1
         end do
         
         sum1=0.D0
         do i=1,nelec-1
            tm1=(10.D0**etemp(i+1)-10.D0**etemp(i))*.5D0*((10.D0**
     *          dtemp(i))+(10.D0**dtemp(i+1)))
            sum1=sum1+tm1
         end do

cjw         renorm=ntot/sum1
         renorm=log10(ntot/sum1)
         do i=1,nelec
cjw            dtemp(i)=log10(renorm*10.D0**dtemp(i))
            dtemp(i)=renorm+dtemp(i)
         end do


c   spline_dbl and resize the array into the old lelec,ledens if you
c   lose particles at low energy
         if (nw.gt.0) then
            if (nw.eq.1) then
               totlos=(10.D0**etemp(2)-10.D0**etemp(1))*10.D0**dtemp(1)
            else if (nw.gt.1) then
               totlos=0.D0
               do i=1,nw-1
                  tmlos=(10.D0**etemp(i+1)-10.D0**etemp(i))*.5D0*(10.D0
     *                  **dtemp(i)+10.D0**dtemp(i+1))
                  totlos=totlos+tmlos
               end do
            endif
            oldnum=ntot
            ntot=ntot-totlos
         
cjw            renorm=ntot/oldnum
            renorm=log10(ntot/oldnum)
            do i=1,nelec
cjw               dtemp(i)=log10(renorm*10.D0**dtemp(i))
               dtemp(i)=renorm+dtemp(i)
            end do
            
            call initspline(nelec,etemp,dtemp,ytb,ytc,ytd)

            einc=(etemp(nelec)-etemp(nw+1))/nelec
            do i=1,nelec
               lelec(i)=etemp(nw+1)+(i-.5D0)*einc
               call interspline(nelec,etemp,dtemp,ytb,ytc,ytd,lelec(i),
     *               ledens(i))

               elen(i)=10.D0**lelec(i)
               eled(i)=elen(i)*(10.D0**ledens(i))
            end do
         else
            do i=1,nelec
               lelec(i)=etemp(i)
               elen(i)=10.D0**lelec(i)
               ledens(i)=dtemp(i)
               eled(i)=elen(i)*(10.D0**ledens(i))
            end do
         endif


c  Double check equip
         sum1=0.D0
         do i=1,nelec-1
            tm1=(elen(i+1)-elen(i))*.5D0*(eled(i)+eled(i+1))
            sum1=sum1+tm1
         end do
         ub=bfield**2/(8.D0*pi)

 23      continue

c get average electron energy
         sum1=0.d0
         sum2=0.d0
         do i=1,nelec-1
            tm1=(elen(i+1)-elen(i))*.5D0*(eled(i)+eled(i+1))
            tm2=(elen(i+1)-elen(i))*.5D0*(10.d0**ledens(i+1)+10.d0**
     *             ledens(i))
            sum1=sum1+tm1
            sum2=sum2+tm2
         end do
         avelec=sum1/sum2

         if (ntot.le.0.D0) then
            do i=1,nsyn      
               do m=1,njet
                  nusyn(m,k,i)=log10(nurad(i)*dopfac(m,k))
                  synabs(m,k,i)=-100.D0
               end do
            end do
            print*,'ran out of particles'
            stop
         endif
         call initspline(nelec,lelec,ledens,yelb,yelc,yeld)
         
         elenmn=elen(1)
         elenmx=elen(nelec)
         do i=1,nelec
            tstrm(i)=(10.D0**ledens(i))/elen(i)**2
         end do
         do i=1,nelec-1
            drtrm(i)=(tstrm(i+1)-tstrm(i))/(elen(i+1)-elen(i))
         end do
         drtrm(nelec)=drtrm(nelec-1)
         call initspline(nelec,elen,drtrm,yderb,yderc,yderd)



c  4.*pi in synabs is for assumed isotropic source
c  synint returns units of ergs/cm^2/s/Hz for absorbed 
c  convert to Jansky=1.e-23 erg/cm^2/s/Hz


         do i=1,nsyn
            do m=1,njet
               nusyn(m,k,i)=log10(nurad(i)*dopfac(m,k))
            enddo
         enddo

         do i=1,nsyn
cjw         ...first perform integrations for sychrotron emissivity and
cjw         ...absorption at current frequency
            call synintegrals(nurad(i),esum,asum)
            do m=1,njet
               call synint(nurad(i),r,inclin,dopfac(m,k),delz,
     &                esum,asum,absd,phoden)
               nphot(i)=phoden
               if (absd.lt.1.D-100) then
                  do j=i,nsyn
                     do l=1,njet
                        synabs(l,k,j)=-100.D0
                     end do
                     nphot(j)=0.D0
                  end do
                  goto 8
               endif
               synabs(m,k,i)=log10(absd*dopfac(m,k)**2/mjy)
            end do
         end do
         

 8       continue

         if (z.gt.zcut) then
            goto 999
         endif


c  spline_dbl photon density for passing to compton routine, divide by herg
c  *clite*energy(erg)*pi*r**2 to get #/cm^3/erg (see pg 12, ntbk 4)
c  NEW: add in the diluted blackbody from the disk, see pg 266/293 ntbk III
c  Newer:  make the max energy based on tin...
         

         do i=1,nsyn
            if (nphot(i).eq.0)then 
               nphot(i)=nphot(i-1)-4.D0
            else
               nphot(i)=log10(nphot(i))
            endif
            jetu(i)=10.D0**nphot(i)/(clite*herg*energ(i)*pi*r**2)
            phodis(i)=log10(jetu(i))
         end do
         
c CHECK SSC ENERGY DENSITY TO COMPARE TO UCOM1 AT SHOCK

         sum1=0.d0
         do i=1,nsyn-1
            tm1=(10.d0**nubb(i+1)-10.d0**nubb(i))*.5d0*
     *           (10.d0**nphot(i+1)+10.d0**nphot(i))
            sum1=sum1+tm1
         end do
         toten=sum1/(clite*pi*r**2)
c         print*,'ENDEN FOR SSC',toten

c CHECK ENERGY DENSITIES

         sum1=0.d0
         sum2=0.d0
         do i=1,nsyn-1
            tm1=(10.d0**(phodis(i+1)+ephot(i+1))+10.d0**(phodis(i)+
     *          ephot(i)))*0.5d0*(10.d0**ephot(i+1)-10.d0**ephot(i))
            tm2=(10.d0**phodis(i+1)+10.d0**phodis(i))*0.5d0*(
     *         10.d0**ephot(i+1)-10.d0**ephot(i))  
            sum1=sum1+tm1
            sum2=sum2+tm2
         end do

         avphot=sum1/sum2

         ephmax=10.D0**ephot(nsyn)
         ephmin=10.D0**ephot(1)
         call initspline(nsyn,ephot,phodis,cb,cc,cd)


c  call integrate compton function over electrons, returns #/cm^3/s/erg

         do i=1,ncom
c           (ephxr is in Hz)         
            ephpass=ephxr(i)*herg
            do m=1,njet
               nucom(m,k,i)=log10(ephpass*dopfac(m,k)/herg)
            end do
cjw         avoid out of bounds error
            if (i.ge.3.and. comspc(1,k,i-1).le.-20.D0.and.
     *              comspc(1,k,i-1)-comspc(1,k,i-2).le.0.D0) then
               do m=1,njet
                  comspc(m,k,i)=comspc(m,k,i-1)-4.D0
               end do
               goto 988
            endif

c   see (224,iv) for why this line is important:  limits integration to right
c   range
            typmax=avphot*(avelec/emerg)**2
            if (ephpass.ge.typmax) then
               blim=max(log(elenmn/emerg),log(sqrt(ephpass/
     *              avphot)/30.d0),log(ephpass/emerg))+1.d-15
            else
               blim=max(log(elenmn/emerg),log(ephpass/emerg))+1.d-15
            endif

            if (blim.ge.log(elenmx/emerg)) then
               com=0.d0
            else
               call qrombe_dbl(comint,5.D-4,blim,log(elenmx/emerg),com)
            endif 
c            print*,'fullcom success',i,k

c           convert to ergs/cm^2/s/Hz and then mjy
c           write in Hz and mJy
            do m=1,njet
               comspc(m,k,i)=com*ephxr(i)*herg**2*vol*dopfac(m,k)**2
     *              /(mjy*4.D0*pi*dist**2)
               if (comspc(m,k,i).eq.0.d0) then
                  comspc(m,k,i)=comspc(m,k,i-1)-2.D0
               else
                  comspc(m,k,i)=log10(comspc(m,k,i))
               endif
            end do
 988     end do  

 999     continue

      end do

c  spline_dbl arrays for adding to big total array for z < zmax, and add
c  BB flux, nutot in Hz

      do m=1,njet
          do k=1,nz
cDM        do k=1,5 or 10 as per reffrac_base.f ????
            z=zmin+zinc*(k-1)
            do i=1,nsyn
               snu(i)=nusyn(m,k,i)
               sdump(i)=synabs(m,k,i)
            end do
            call initspline(nsyn,snu,sdump,ysyb,ysyc,ysyd)

            if (z.le.log10(zcut)) then
               do i=1,ncom
                  cnu(i)=nucom(m,k,i)
                  cdump(i)=comspc(m,k,i)
c                  write(ilun21,*) cnu(i),cdump(i)
               end do
               call initspline(ncom,cnu,cdump,ycob,ycoc,ycod)
            endif

            do i=1,ne
               if (nutot(i).ge.snu(1).and.nutot(i).le.
     *              snu(nsyn)) then
                  call interspline(nsyn,snu,sdump,ysyb,ysyc,ysyd,
     *                             nutot(i),sflx)
               else
                  sflx=-200.D0
               endif
               if (z.gt.log10(zcut)) then
                  cflx=-200.D0
               else 
                  if (nutot(i).ge.cnu(1).and.nutot(i).le.cnu(ncom))then
                     call interspline(ncom,cnu,cdump,ycob,ycoc,ycod,
     *                                nutot(i),cflx)
                  else
                     cflx=-200.d0
                  endif
               end if
               if (plotsw.eq.1.d0) then
                  complot(i)=complot(i)+10.d0**cflx
                  if (z.lt.log10(zsh)) then
                     presyn(i)=presyn(i)+10.d0**sflx
                  else
                     postsyn(i)=postsyn(i)+10.d0**sflx
                  endif
               endif
               fplot(i)=fplot(i)+10.D0**sflx+10.D0**cflx
            end do
         end do
      end do



      do i=1,ne
         frq=10.D0**nutot(i)
         if (plotsw.eq.1.d0) then
            if (fplot(i).eq.0) then
               write(ilun12,*) nutot(i), -200.D0
            else
               write(ilun12,*) nutot(i),log10(fplot(i))
            endif
            if (complot(i).eq.0) then
               write(ilun13,*) nutot(i), -200.D0
            else
               write(ilun13,*) nutot(i),log10(complot(i))
            endif
            if (presyn(i).eq.0) then
               write(ilun14,*) nutot(i), -200.D0
            else
               write(ilun14,*) nutot(i),log10(presyn(i))
            endif
            if (postsyn(i).eq.0) then
               write(ilun16,*) nutot(i), -200.D0
            else
               write(ilun16,*) nutot(i),log10(postsyn(i))
            endif
         endif
         if (fplot(i).eq.0) then
            photeng(i) = nutot(i)
            photar(i) = -200.D0
         else
            photeng(i) =  nutot(i)
            photar(i) = log10(fplot(i))
         endif
      end do


      if (plotsw.eq.1.d0) then
         close(ilun12)
         close(ilun13)
         close(ilun14)
         close(ilun16)
c         close(ilun21)
         call frelun(-1) 
      endif


 88   format (5(3x,G14.7))
 89   format(a,5(3x,g14.7))
 100  format(a,f3.1,a)
 101  format(a,f4.1,a)
 103  format(4(2x,g14.7))


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine synintegrals:  does frequency and electron spectrum 
c     dependent calc, single component returns kernals necessary for 
c     emissivity and self-absorption for use in synint, which calculates
c     the actual spectrum
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine synintegrals(freq,esum,asum)
      double precision esum,asum,freq
      double precision blim,ulim,elenmn,elenmx,bfield,dist
      double precision synemis,absfnc,nupass

      common /eldis2/elenmn,elenmx,bfield,dist
      common /spass/nupass
      external synemis,absfnc

      nupass=freq
      blim=log(elenmn)
      ulim=log(elenmx)
      call qrombe_dbl(synemis,1.D-5,blim,ulim,esum)

c   see pg 172-179 in ntbk III.  using absorbed spectrum for scattering
c   updated pg 12, iv

      call qrombe_dbl(absfnc,1.D-5,blim,ulim,asum)

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine synint...for selfabsorption, single component
c      returns specific intensity in erg/s/cm^2/st/Hz
c     integrates distributions over electron spectrum
c     sin(pitch)=2/3 equiv to angle averaging over isotropic
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine synint(freq,r,angle,dopfac,h,esum,asum,fluxa,fluxa2)
      include 'param_dbl.inc'
      DOUBLE PRECISION freq,angle,fluxa,dopfac,fluxa2
      double precision pitch
      parameter(pitch=0.73D0)
      DOUBLE PRECISION h,area,esum,acons,asum,tsyn
      DOUBLE PRECISION absfac,tsyn2,absfac2
      DOUBLE PRECISION r,asyn
      double precision elcons,elenmn,elenmx,bfield,dist
      double precision epsasyn
      common /eldis2/elenmn,elenmx,bfield,dist

      elcons=sqrt(3.D0)*(charg**3)*bfield*sin(pitch)/emerg
      acons=-clite**2/(8.D0*pi*freq**2)
cjw      eps=elcons*esum
      asyn=acons*elcons*asum
      epsasyn=esum/(acons*asum)
      area=4.D0*pi*dist**2

c  this first term is for what is seen externally...includes skin depth
c  and angle effects...see pg 12 ntbk 4 for consideration..

      tsyn=pi/2.D0*asyn*r/(dopfac*sin(angle))

      if (tsyn.ge.1.D0) then
         absfac=(1.D0-exp(-tsyn))
      else
cjw         absfac=tsyn-tsyn**2/2.D0+tsyn**3/6.D0-tsyn**4/24.D0+tsyn**5/
cjw     *          120.D0
         absfac=tsyn*(1.D0+tsyn*(-0.5D0+tsyn*(1.D0/6.D0+
     &        tsyn*(-1.D0/24.D0+tsyn/120.D0))))

      endif
c  to test what it looks like with minimal self absorption:
c      tsyn=.1
c      absfac=tsyn
c      print*,tsyn,absfac

cjw      fluxa=(2.D0*r*h*sin(angle)*dopfac*absfac*eps/asyn)/area
      fluxa=(2.D0*r*h*sin(angle)*dopfac*absfac*epsasyn)/area

c  this second term is the same as above, but assuming what is "seen" locally
c  by the particles for compton scattering.  pass 'flux', erg/s/Hz to 
c  main code
      
c  erg/s/cm^2/hz:

      tsyn2=pi/2.D0*asyn*r
      if (tsyn2.ge.1.D0) then
         absfac2=(1.D0-exp(-tsyn2))
      else
cjw         absfac2=tsyn2-tsyn2**2/2.D0+tsyn2**3/6.D0-tsyn2**4/24.D0+tsyn2
cjw     *           **5/120.D0

         absfac2=tsyn2*(1.D0+tsyn2*(-0.5D0+tsyn2*(1.D0/6.D0+
     &        tsyn2*(-1.D0/24.D0+tsyn2/120.D0))))

      endif

cjw      fluxa2=2.D0*r*h*absfac2*eps/asyn
      fluxa2=2.D0*r*h*absfac2*epsasyn

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     synemis--function for synchrotron emissivity to be integrated over 
c              particle distribution in synint
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function synemis(e1e)
      include 'param_dbl.inc'
      integer nelec
      parameter(nelec=200)
      DOUBLE PRECISION synemis,e1e,ele,game,x,func,funcl,ellog,eden
      DOUBLE PRECISION arg,val
      DOUBLE PRECISION elenmn,elenmx,bfield,nupass,dist
      double precision lelec,ledens,yelb,yelc,yeld
      double precision synb,sync,synd
      common /syn/arg(47),val(47),synb(47),sync(47),synd(47)
      common /eldis1/lelec(nelec),ledens(nelec),yelb(nelec),
     *yelc(nelec),yeld(nelec)
      common /eldis2/elenmn,elenmx,bfield,dist
      common /spass/nupass

      ele=exp(e1e)
      game=ele/emerg
      x=nupass*4.D0*pi*emgm*clite/(3.D0*charg*bfield*game**2)

      if (x.le.1.D-4) then
         func=4.D0*pi*(x/2.D0)**(1.D0/3.D0)/(sqrt(3.D0)*2.68D0)
      else if (x.gt.50.D0)then
         func=sqrt(pi*x/2.D0)*exp(-x)
      else
         call interspline(47,arg,val,synb,sync,synd,x,funcl)
         func=10.D0**funcl
      endif     

cjw   ellog=log10(ele)
cjw   0.434... is log10(exp(1d0))  
      ellog=e1e*0.434294481903252D0
      call interspline(nelec,lelec,ledens,yelb,yelc,yeld,ellog,eden)
      eden=10.D0**eden
      synemis=ele*eden*func

      return
      end
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     absfnc--function for self absorption, to be integrated over the derived
c             function of the particle distribution in synint
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function absfnc(e2)
      include 'param_dbl.inc'
      integer nelec
      parameter(nelec=200)
      DOUBLE PRECISION absfnc,e2,elec,game,x,func,funcl,deden
      DOUBLE PRECISION arg,val
      DOUBLE PRECISION elenmn,elenmx,bfield,nupass,dist
      double precision elen,drtrm,yderb,yderc,yderd
      double precision synb,sync,synd
      common /syn/arg(47),val(47),synb(47),sync(47),synd(47)
      common /eldis2/elenmn,elenmx,bfield,dist
      common /derivs/elen(nelec),drtrm(nelec),yderb(nelec),yderc(nelec),
     &yderd(nelec)
      common /spass/nupass

      elec=exp(e2)
      game=elec/emerg
      x=nupass*4.D0*pi*emgm*clite/(3.D0*charg*bfield*game**2)
      if (x.le.1.D-4) then
         func=4.D0*pi*(x/2.D0)**(1.D0/3.D0)/(sqrt(3.D0)*2.68D0)
      else if (x.gt.50.D0)then
         func=sqrt(pi*x/2.D0)*exp(-x)
      else
         call interspline(47,arg,val,synb,sync,synd,x,funcl)
         func=10.D0**funcl
      endif     

c   see functional form, pg 164 ntbk 3

      call interspline(nelec,elen,drtrm,yderb,yderc,yderd,elec,deden)
cjw      absfnc=elec*func*(elec**2)*deden
      absfnc=func*(elec**3)*deden
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function comint for SSC of synchro flux--numerically integrates 
c     over photon distribution and cross section
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function comint(gam)
      include 'param_dbl.inc'
      integer nelec
      parameter(nelec=200)
      DOUBLE PRECISION comint,gam,game,e1,gnorm,valu,ellog,den
      DOUBLE PRECISION ephpass,gpass,ephmax,ephmin,econst
      double precision lelec,ledens,yelb,yelc,yeld
      double precision blim,ulim
      common /com2/ephpass,ephmax,ephmin
      common /eldis1/lelec(nelec),ledens(nelec),yelb(nelec),
     *yelc(nelec),yeld(nelec)
      common /cpass/gpass,e1
      external comfnc
         
      game=exp(gam)
      gpass=game
      e1=ephpass/(game*emerg)
c      print*,e1
c  convert N(E)dE to N(gam)dgam....see pg 170 ntbk III
      
      econst=2.d0*pi*(re0**2)*clite
      gnorm=emerg

c DON'T CHANGE!  OTHERWISE GOES OUT OF ARRAY BOUNDS
c added incremental to avoid hitting blim/utst boundary
      blim=max(log(ephpass/(4.d0*game*(game-ephpass/emerg))),
     *         log(ephmin))+1.d-14
      ulim=log(min(ephpass,ephmax))
c      eg4b=4.d0*exp(blim)*game
c      eg4u=4.d0*exp(ulim)*game
c      btst=eg4b/(emerg+eg4b)
c      utst=eg4u/(emerg+eg4u)


      if (ulim.le.blim) then 
         comint=0.d0
         return
      endif

c     DON'T CHANGE THIS EPS!!!

c      print*,'blim,ulim',blim,ulim
c      print*,'before 1'
      call qrombe_dbl2(comfnc,5.D-6,blim,ulim,valu)
c      print*,'after 1'

      ellog=log10(game*emerg)
      call interspline(nelec,lelec,ledens,yelb,yelc,yeld,ellog,den)

      den=10.D0**den
c fullexpress:   comint=game*gnorm*econst*den*valu/game**2
      comint=gnorm*econst*den*valu/game

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function comfnc--kernel of eq 2.48 in Blumenthal & Gould 1970
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function comfnc(ein)
      include 'param_dbl.inc'
      DOUBLE PRECISION ein,einit,utst,biggam,q,elg,phonum
      double precision tm1,tm2,tm3
cjw   note: nbin has to equal parameter nsyn in the main program
cjw   (should we put this common block in an include file?)
      integer nbin
      parameter(nbin=100)
      DOUBLE PRECISION ephot,phodis,yphot,gpass,e1,btst
      common /com/ephot(nbin),phodis(nbin),yphot(nbin)
      common /cpass/gpass,e1

      double precision eg4
      double precision cb,cc,cd
      common /comsplinpar/ cb(nbin),cc(nbin),cd(nbin)

      comfnc=0.0d0
      einit=exp(ein)
      btst=einit/(gpass*emerg)
cjw      utst=4.D0*einit*gpass/(emerg*(1.D0+4.D0*einit*gpass/emerg))
      eg4=4.D0*einit*gpass
      utst=eg4/(emerg+eg4)

cjw
cjw the following costs 3 seconds in a typical run
cjw are these sanity checks really needed?
cjw

      if ( (e1.lt.btst.and.(btst-e1)/btst.ge.3.d-8).or.(e1.gt.utst.and.
     * (e1-utst)/utst.ge.3.d-8)) then
         print*,'gotfailures',btst,e1,utst
c         print*,btst-e1,(btst-e1)/btst,(e1-utst)/utst
          comfnc=0.D0
          return
      endif

      biggam=eg4/emerg
      q=e1/(biggam*(1.D0-e1))
cjw   elg=log10(einit)
cjw   0.434... is log10(exp(1d0))  
      elg=ein*0.434294481903252D0

      call interspline(nbin,ephot,phodis,cb,cc,cd,elg,phonum)

      tm1=2.D0*q*log(q)
      tm2=(1.D0+2.D0*q)*(1.D0-q)
      tm3=0.5D0*((biggam*q)**2)*(1.D0-q)/(1.D0+biggam*q)
c      print*,tm1,tm2,tm3,phonum
      comfnc=(tm1+tm2+tm3)*10.d0**phonum
      
c     took out einit: here's full expression
c      comfnc=einit*(2.D0*q*log(q)+(1.D0+2.D0*q)*(1.D0-q)+.5D0*
c     *        ((biggam*q)**2)*
c     *       (1.D0-q)/(1.D0+biggam*q))*(10.D0**phonum)/einit


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jetpars(z,b0,n0,g0,r0,h0,gb,b,n,g,r)
      integer maxnv
      DOUBLE PRECISION clite
      parameter (maxnv=120,clite=3.D10)
      integer nv
      DOUBLE PRECISION x,y,gbx(maxnv),gby(maxnv),gbs0
      double precision gbyb(maxnv),gbyc(maxnv),gbyd(maxnv)
      DOUBLE PRECISION visco,mbh,bff,gbin,zmin,zff
      common /jpars/visco,mbh,zmin
      common /cbgb/gbx,gby,gbyb,gbyc,gbyd,gbs0,nv

      DOUBLE PRECISION z,b0,n0,g0,r0,h0,gb,b,n,g,r,mj

c     Input:
c     z  = distance from origin(!) along jet axis (same units as r0,h0)
c     b0 = magnetic field in nozzle
c     n0 = density in nozzle
c     g0 = electron Lorentz factor in nozzle
c     r0 = width of nozzle
c     h0 = height of nozzle

c     Output:
c     gb= bulk velocity of jet (gamma*beta)
c     b = magnetic field at z
c     n = density at z
c     g = electron Lorentz factor at z
c     r = width at z


c     mj = Mach number, but relative to soundspeed in nozzle!
c            local Mach number will be slightly higher due to cooling

      DOUBLE PRECISION betas0,Gammas
      Gammas=4.d0/3.D0
      betas0=sqrt((Gammas-1.D0)/(Gammas+1.D0))
cjw      gbs0=betas0/sqrt(1.D0-betas0**2)
      gbs0=betas0/sqrt((1.D0-betas0)*(1.d0+betas0))

      x=log10((max(z-h0,0.0D0)+r0)/r0)



c      Write (*,*) 'call splint_dbl with x=',x
      
      zff=10.d0**zmin
      if (z.lt.h0) then
         if (visco.le.0.D0)then
            y=gbs0
         else
            bff=sqrt(2.D0*6.67D-8*2.D33*mbh/r0)/clite*visco
            gbin=bff/sqrt(1.D0-bff**2)
            y=gbin+(gbs0-gbin)/(h0-zff)*(z-zff)
         endif
      else
        call interspline(nv,gbx,gby,gbyb,gbyc,gbyd,x,y)
      endif


      gb=y
      mj=gb/gbs0
      r=r0+max(z-h0,0.0D0)/mj
      n=n0*(r0/r)**2/mj
      if (z.lt.h0) then
c      if (z.lt.0) then 
         b=b0*(r0/r)/mj**(0.5D0)
         g=g0
      else
         b=b0*(r0/r)/mj**(0.5D0+1.D0/6.D0)
         g=g0/mj**(1.D0/3.D0)
      endif

c      write (*,10) z,r,mj,gb,b,n,g
 10   format (7(3x,G12.5))
      return
      end

      subroutine initjet
      integer maxnv
      parameter (maxnv=120)
      integer nv
      DOUBLE PRECISION gbx(maxnv),gby(maxnv),gbs0
      double precision gbyb(maxnv),gbyc(maxnv),gbyd(maxnv)
      common /cbgb/gbx,gby,gbyb,gbyc,gbyd,gbs0,nv

c     Looks a bit ugly, but better explicitly define the arrays
c     gbx(nv) and gby(nv) in the code so that the file v.dat is
c     no longer needed. Note v.dat had only 101 rows with data!
      data gby/0.436646d0, 0.862060d0, 1.059490d0, 1.216690d0,
     *         1.352370d0, 1.473980d0, 1.585380d0, 1.688930d0,
     *         1.786170d0, 1.878200d0, 1.965820d0, 2.049640d0,
     *         2.130140d0, 2.207710d0, 2.282660d0, 2.355250d0,
     *         2.425710d0, 2.494210d0, 2.560930d0, 2.626000d0,
     *         2.689540d0, 2.751660d0, 2.812450d0, 2.872000d0,
     *         2.930390d0, 2.987690d0, 3.043940d0, 3.099230d0,
     *         3.153580d0, 3.207060d0, 3.259700d0, 3.311540d0,
     *         3.362620d0, 3.412980d0, 3.462640d0, 3.511630d0,
     *         3.559990d0, 3.607730d0, 3.654880d0, 3.701460d0,
     *         3.747490d0, 3.793000d0, 3.838000d0, 3.882500d0,
     *         3.926530d0, 3.970100d0, 4.013220d0, 4.055910d0,
     *         4.098180d0, 4.140040d0, 4.181520d0, 4.222600d0,
     *         4.263320d0, 4.303670d0, 4.343670d0, 4.383330d0,
     *         4.422650d0, 4.461650d0, 4.500330d0, 4.538700d0,
     *         4.576770d0, 4.614540d0, 4.652030d0, 4.689230d0,
     *         4.726160d0, 4.762830d0, 4.799230d0, 4.835370d0,
     *         4.871260d0, 4.906910d0, 4.942310d0, 4.977480d0,
     *         5.012420d0, 5.047130d0, 5.081620d0, 5.115890d0,
     *         5.149950d0, 5.183800d0, 5.217440d0, 5.250880d0,
     *         5.284120d0, 5.317170d0, 5.350030d0, 5.382700d0,
     *         5.415180d0, 5.447490d0, 5.479610d0, 5.511570d0,
     *         5.543340d0, 5.574950d0, 5.606400d0, 5.637680d0,
     *         5.668790d0, 5.699750d0, 5.730550d0, 5.761200d0,
     *         5.791700d0, 5.822050d0, 5.852250d0, 5.882300d0,
     *         5.912220d0,
     *         0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     *         0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     *         0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /
      data gbx/120*0.d0/
      do nv=1,101
       gbx(nv)=(nv-1.d0)*0.2d0
      enddo

c      Write (*,*) nv,' Points read.'
c      call spline_dbl(gbx,gby,nv,1.0D0,0.0D0,gby2)
c      call splint_dbl(gbx,gby,gby2,nv,0.d0,gbs0)

      call initspline(nv,gbx,gby,gbyb,gbyc,gbyd)
      call interspline(nv,gbx,gby,gbyb,gbyc,gbyd,0.d0,gbs0)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION function k0(x)
      implicit none
      DOUBLE PRECISION x,y,t,i0,t2,y2,y1
c
c  k0 from A&S good to 2.d-7
c
      y = x/2.d0
      y2=y*y
      y1=1.d0/y

      if (y.lt.1.D0) then
        t = x/3.75d0
        t2= t*t
c        i0 = (1.D0+3.5156229d0*t**2+3.0899424d0*t**4+
c     +       1.2067492d0*t**6+0.2659732d0*t**8+0.0360768d0*t**10+
c     +       0.0045813d0*t**12)

        i0 = 1.D0+t2*(3.5156229d0+t2*(3.0899424d0+
     +       t2*(1.2067492d0*t**6+t2*(0.2659732d0*t**8+t2* (0.0360768d0+
     +       t2*0.0045813d0)))))

c        k0 = log(2.D0/x)*i0-0.5772156649d0+
c     +       0.42278420d0*y**2+0.23069756d0*y**4+
c     +       0.03488590d0*y**6+0.00262698d0*y**8+
c     +       0.00010750d0*y**10+0.00000740d0*y**12


        k0 = log(2.D0/x)*i0-0.5772156649d0+
     +       y2*(0.42278420d0+y2*(0.23069756d0+
     +       y2*(0.03488590d0+y2*(0.00262698d0+
     +       y2*(0.00010750d0+y2*0.00000740d0)))))
      else

c        k0 = exp(-x)/sqrt(x)*(1.25331414d0-0.07832358d0/y+
c     +          0.02189568d0/y**2-0.01062446d0/y**3+
c     +          0.00587872d0/y**4-0.00251540d0/y**5+
c     +          0.00053208d0/y**6)

         k0 = exp(-x)/sqrt(x)*(1.25331414d0+y1*(-0.07832358d0+
     +        y1*(0.02189568d0+y1*(-0.01062446d0+
     +        y1*(0.00587872d0+y1*(-0.00251540d0+
     +        y1*0.00053208d0))))))
      end if
      return
      end
c----------------------------------------------------------------
      DOUBLE PRECISION function k1(x)
      implicit none
      DOUBLE PRECISION x,y,t,i1,t2,y2,Y1
c
c  k1 from A&S good to 2.d-7
c
      y = x/2.d0
      y2=y*y
      y1=1.d0/y

      if (y.lt.1.D0) then
        t = x/3.75d0
        t2=t*t
c        i1 = x*(0.5d0+0.87890594d0*t**2+
c     +          0.51498869d0*t**4+
c     +          0.15084934d0*t**6+0.02658733d0*t**8+
c     +          0.00301532d0*t**10+
c     +          0.00032411d0*t**12)

        i1=x*(0.5d0+t2*(0.87890594d0 + t2*(0.51498869d0 +
     +     t2*(0.15084934d0 + t2*(0.02658733d0+
     +     t2*(0.00301532d0 + t2*0.00032411d0))))))

c        k1 = 1.D0/x*(x*log(y)*i1+1.d0+
c     +          0.15443144d0*y**2-0.67278579d0*y**4-
c     +          0.18156897d0*y**6-0.01919402d0*y**8-
c     +          0.00110404d0*y**10-0.00004686d0*y**12)


        k1 = 1.D0/x*(x*log(y)*i1+1.d0+
     +          y2*(0.15443144d0-y2*(0.67278579d0-
     +          y2*(0.18156897d0-y2*(0.01919402d0-
     +          y2*(0.00110404d0-y2*0.00004686d0))))))

      else
c        k1 = exp(-x)/sqrt(x)*(1.25331414d0+0.23498618d0/y-
c     +          0.03655620d0/y**2+0.01504268d0/y**3-
c     +          0.00780353d0/y**4+0.00325614d0/y**5-
c     +          0.00068245d0/y**6)

        k1 = exp(-x)/sqrt(x)*(1.25331414d0+y1*(0.23498618d0-
     +          y1*(0.03655620d0+y1*(0.01504268d0-
     +          y1*(0.00780353d0+y1*(0.00325614d0-
     +          y1*0.00068245d0))))))
      end if
      return
      end
c----------------------------------------------------------------
      DOUBLE PRECISION function k2(x)
      implicit none
      DOUBLE PRECISION k0,k1,x
c
c  k2 from A&S is just:
c
      k2 = 2.d0*k1(x)/x+k0(x)
      return
      end
c----------------------------------------------------------------
      DOUBLE PRECISION function k3(x)
      implicit none
      DOUBLE PRECISION k1,k2,x
c
c  k3 from A&S is just:
c
      k3 = 4.d0*k2(x)/x+k1(x)
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE trapzd_dbl(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x

      if (n.eq.1) then
        s=0.5D0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5D0*del
        sum=0.D0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5D0*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software jA##!+)2+.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE trapzd_dbl2(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x

      if (n.eq.1) then
        s=0.5D0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5D0*del
        sum=0.D0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5D0*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software jA##!+)2+.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE polint_dbl(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.D0) then
             print*,'crashed in polint'
             return
          endif
c             pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software jA##!+)2+.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE polint_dbl2(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.D0) then
             print*,'crashed in polint2'
             return
          endif
c          pause 'failure in polint2'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software jA##!+)2+.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     qrombe_dbl--qromb modified so that eps is chosen within main to fit needs
c                 and double precision for codes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE qrombe_dbl(func,EPS,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (JMAX=35, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
      
      h(1)=1.D0
      do 11 j=1,JMAX
        call trapzd_dbl(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint_dbl(h(j-KM),s(j-KM),K,0.D0,ss,dss)
c          print*,j,abs(dss),eps*abs(ss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)

        h(j+1)=0.25D0*h(j)
11    continue
c      pause 'too many steps in qromb'
      print*,'too many steps in qromb',abs(dss),eps*abs(ss)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software jA##!+)2+.


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     qrombe_dbl2--qrombe_dbl with new name for nested loops
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE qrombe_dbl2(func,EPS,a,b,ss)
      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (JMAX=35, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint2,trapzd2
      INTEGER j
      DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
      
      h(1)=1.D0
      do 11 j=1,JMAX
        call trapzd_dbl2(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint_dbl2(h(j-KM),s(j-KM),K,0.D0,ss,dss)
c          print*,j,abs(dss),eps*abs(ss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)

        h(j+1)=0.25D0*h(j)
11    continue
c      pause 'too many steps in qromb'
      print*,'too many steps in qromb2',abs(dss),eps*abs(ss)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software jA##!+)2+.




C------------------------------------------------------------------------------
        subroutine getlun(iounit)

C       get an unallocated logical unit number
C
C        O  (i) iounit - An unopened logical unit number
C       
C       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996
C       Copied with minor changes from the FITSIO routine ftgiou.

        integer iounit

        iounit=0
        call lunlst(iounit)
        end
C------------------------------------------------------------------------------
        subroutine frelun(iounit)

C       free specified logical unit number; if iounit=-1, then free all units
C       
C        I  (i) iounit - The logical unit number to be freed
C       
C       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996
C       Copied with minor changes from the FITSIO routine ftfiou.

        integer iounit

        call lunlst(iounit)
        end
C------------------------------------------------------------------------------
        subroutine lunlst(iounit)

C       generic routine to manage logical unit numbers in the range 10-49
C       
C       I/O  (i) iounit - The logical unit number to be allocated/freed
C       
C       James Peachey, HEASARC/GSFC/NASA  Hughes STX, November, 1996
C       Copied with minor changes from the FITSIO routine ftxiou.

        integer iounit,i
        integer array(40)
        save array
        data array/40*0/

        if (iounit .eq. 0)then
C           get an unused logical unit number
            do 10 i=40,1,-1

C        The following would be a more robust way of testing for
C        an available unit number, however, this cannot work
C        when building XANLIB using the IRAF/SPP version, because
C        IRAF does not use Fortran I/O.
C
C                inquire(unit=iounit, exist=exists, opened=open)
C                if(exists .and. .not. open)then
C                    array(iounit-9)=1
C                    return
C                end if

                 if (array(i) .eq. 0)then
                     array(i)=1
                     iounit=i+9
                     return
                 end if
10          continue
C           error: all units are allocated
            iounit=-1
c            call xaerror(
c     &           'GETLUN has no more available unit numbers.', 1)
c            print*,'something wrong with units'

        else if (iounit .eq. -1)then
C           deallocate all the unit numbers
            do 20 i=1,40
                 array(i)=0
20          continue

        else
C            deallocat a specific unit number
             if (iounit .ge. 10 .and. iounit .le. 49)then
                 array(iounit-9)=0
             end if
        endif
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initspline (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)

cjw
cjw   source: http://www.netlib.org/fmm/spline.f, and modified
cjw
cjw   note: boundary conditions here are not fully what we might
cjw         need, but I think this can easily be patched in the
cjw         production version

c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return

      if ( n .eq. 2 ) then
         b(1) = (y(2)-y(1))/(x(2)-x(1))
         c(1) = 0.
         d(1) = 0.
         b(2) = b(1)
         c(2) = 0.
         d(2) = 0.
         return
      endif
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
      enddo
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if (n.gt.3) then 
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
         c(1) = c(1)*d(1)**2/(x(4)-x(1))
         c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      endif 
c
c  forward elimination
c
      do i=2,n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
      enddo
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      enddo
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
      enddo
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interspline(n, xa, y, cb,  cc,  cd, xp, yres)
      integer n
      double precision xa(n),y(n),cb(n),cc(n),cd(n),yres,xp,dx
      INTEGER k,khi,klo

      integer enn
      logical dosearch

      save klo,enn
      data enn/-1/

c     
c     caching of last search interval
c
      dosearch=.true.
      if ( n.eq.enn.and. klo.lt.n ) then
         khi=klo+1
c        ... is it still valid?
         if (xa(klo).le.xp .and. xa(khi).gt.xp ) then
            dosearch=.false.
         else 
c           ... well, since we're usually searching with the search
c           ... parameter in ascending order, perhaps the neighboring
c           ... interval is the one we want?
            if (khi.lt.n) then
               klo=klo+1
               khi=khi+1
               if (xa(klo).le.xp .and. xa(khi).gt.xp ) then
                  dosearch=.false.
               endif
            endif
         endif
      endif
c
c     ... if we weren't successful above, do a binary search
c
      if (dosearch) then 
         klo=1
         khi=n
 1       if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if(xa(k).gt.xp)then
               khi=k
            else
               klo=k
            endif
            goto 1
         endif
c        ... remember the length of the array (good indicator whether
c        ... we're still using the same xa array)
         enn=n

      endif 

c     ... now use the parameters of klo for the evaluation
c     ... using a Horner scheme
      dx=xp-xa(klo)
      yres=(( cd(klo)*dx+cc(klo) )*dx+cb(klo) )*dx + y(klo)

      return
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xrbinterp(ear,energ,phot,photar,ne,newne)
      implicit none
c     Interpolate from grid sized with ne bins to one of newne bins
c     Adapted from S-Lang version, for greater speed (Mike Noble)
      integer ne, newne, i, iplus1, j, jplus1
      double precision ear(ne), energ(ne), phot(ne)
      real*4 photar(newne)
      double precision emid, phflux

      j = 1
      do i=1,newne

c        middle of bin

         iplus1 = i + 1
         emid = (ear(i) + ear(iplus1))/2.0;

c        linear search if we don't bracket yet
         if (j .eq. -1) j = 1

         do while (j .le. ne .and. energ(j) .lt. emid)
            j = j + 1
         enddo
 

         jplus1 = j
         j = j - 1

         if (j .lt. 1 .or . j .gt. ne) then
            photar(i) = 0.0
         else
c               ... ph/cm^2/s/keV
            phflux = phot(j) + (phot(jplus1) - phot(j)) *
     +                   (emid - energ(j)) / (energ(jplus1) - energ(j))

            photar(i) = phflux * (ear(iplus1) - ear(i))
         endif

      enddo

      end

