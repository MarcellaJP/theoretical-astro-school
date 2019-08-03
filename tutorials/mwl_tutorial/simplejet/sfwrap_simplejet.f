c     This file was generated by SLIRP (version pre2.0.0-34)
c     (c) 2005-2008 Massachusetts Institute of Technology
c     It contains subroutine wrappers for the Fortran functions
c     (and possibly COMMON blocks) wrapped by the simplejet module,
c     to maximize portability.
c     mnoble@space.mit.edu

      subroutine synemissfwrap (OUTPUT,e1e)
      external synemis
      double precision e1e
      double precision OUTPUT,synemis
      OUTPUT = synemis(e1e)
      end

      subroutine absfncsfwrap (OUTPUT,e2)
      external absfnc
      double precision e2
      double precision OUTPUT,absfnc
      OUTPUT = absfnc(e2)
      end

      subroutine comintsfwrap (OUTPUT,gam)
      external comint
      double precision gam
      double precision OUTPUT,comint
      OUTPUT = comint(gam)
      end

      subroutine comfncsfwrap (OUTPUT,ein)
      external comfnc
      double precision ein
      double precision OUTPUT,comfnc
      OUTPUT = comfnc(ein)
      end

      subroutine k0sfwrap (OUTPUT,x)
      external k0
      double precision x
      double precision OUTPUT,k0
      OUTPUT = k0(x)
      end

      subroutine k1sfwrap (OUTPUT,x)
      external k1
      double precision x
      double precision OUTPUT,k1
      OUTPUT = k1(x)
      end

      subroutine k2sfwrap (OUTPUT,x)
      external k2
      double precision x
      double precision OUTPUT,k2
      OUTPUT = k2(x)
      end

      subroutine k3sfwrap (OUTPUT,x)
      external k3
      double precision x
      double precision OUTPUT,k3
      OUTPUT = k3(x)
      end

      subroutine sfwrapderivs (initfunc)
      external initfunc
      double precision elen(200)
      double precision drtrm(200)
      double precision yderb(200)
      double precision yderc(200)
      double precision yderd(200)
      common        /derivs/ elen,drtrm,yderb,yderc,yderd
      call initfunc(0,5,elen,drtrm,yderb,yderc,yderd)
      end

      subroutine sfwrapeldis2 (initfunc)
      external initfunc
      double precision elenmn
      double precision elenmx
      double precision bfield
      double precision dist
      common        /eldis2/ elenmn,elenmx,bfield,dist
      call initfunc(1,4,elenmn,elenmx,bfield,dist)
      end

      subroutine sfwrapcom (initfunc)
      external initfunc
      double precision ephot(100)
      double precision phodis(100)
      double precision yphot(100)
      common        /com/ ephot,phodis,yphot
      call initfunc(2,3,ephot,phodis,yphot)
      end

      subroutine sfwrapspass (initfunc)
      external initfunc
      double precision nupass
      common        /spass/ nupass
      call initfunc(3,1,nupass)
      end

      subroutine sfwrapcpass (initfunc)
      external initfunc
      double precision gpass
      double precision e1
      common        /cpass/ gpass,e1
      call initfunc(4,2,gpass,e1)
      end

      subroutine sfwrapbbpars (initfunc)
      external initfunc
      double precision rin
      double precision rout
      double precision hbb
      double precision tin
      double precision zpass
      common        /bbpars/ rin,rout,hbb,tin,zpass
      call initfunc(5,5,rin,rout,hbb,tin,zpass)
      end

      subroutine sfwraperdpar (initfunc)
      external initfunc
      double precision frq
      double precision inclin
      common        /erdpar/ frq,inclin
      call initfunc(6,2,frq,inclin)
      end

      subroutine sfwrapjpars (initfunc)
      external initfunc
      double precision visco
      double precision mbh
      double precision zmin
      common        /jpars/ visco,mbh,zmin
      call initfunc(7,3,visco,mbh,zmin)
      end

      subroutine sfwrapcom2 (initfunc)
      external initfunc
      double precision ephpass
      double precision ephmax
      double precision ephmin
      common        /com2/ ephpass,ephmax,ephmin
      call initfunc(8,3,ephpass,ephmax,ephmin)
      end

      subroutine sfwrapcomsplinpar (initfunc)
      external initfunc
      double precision cb(100)
      double precision cc(100)
      double precision cd(100)
      common        /comsplinpar/ cb,cc,cd
      call initfunc(9,3,cb,cc,cd)
      end

      subroutine sfwrapcbgb (initfunc)
      external initfunc
      double precision gbx(120)
      double precision gby(120)
      double precision gbyb(120)
      double precision gbyc(120)
      double precision gbyd(120)
      double precision gbs0
      integer nv
      common        /cbgb/ gbx,gby,gbyb,gbyc,gbyd,gbs0,nv
      call initfunc(10,7,gbx,gby,gbyb,gbyc,gbyd,gbs0,nv)
      end

      subroutine sfwrapeldis1 (initfunc)
      external initfunc
      double precision lelec(200)
      double precision ledens(200)
      double precision yelb(200)
      double precision yelc(200)
      double precision yeld(200)
      common        /eldis1/ lelec,ledens,yelb,yelc,yeld
      call initfunc(11,5,lelec,ledens,yelb,yelc,yeld)
      end

      subroutine sfwrapjetg (initfunc)
      external initfunc
      double precision gamv
      common        /jetg/ gamv
      call initfunc(12,1,gamv)
      end

      subroutine sfwrapsyn (initfunc)
      external initfunc
      double precision arg(47)
      double precision val(47)
      double precision synb(47)
      double precision sync(47)
      double precision synd(47)
      common        /syn/ arg,val,synb,sync,synd
      call initfunc(13,5,arg,val,synb,sync,synd)
      end
