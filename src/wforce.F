!-----------------------------------------------------------------------------
!
!   /////////////////////         BEGIN         \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\     MODULE FORCE      ////////////////////
!     
!-----------------------------------------------------------------------------

  MODULE FORCE_MODULE
  
!-----------------------------------------------------------------------------

  implicit none

  integer, parameter :: nfmax = 10 ! Max number of forcing areas

  integer ::  iwshap(nfmax) = 1 ! for W forcing
  integer ::  ivshap(nfmax) = 1
  integer ::  iushap(nfmax) = 1
  integer ::  iqshap(nfmax) = 1
  integer ::  itshap(nfmax) = 1
  integer ::  isshap(nfmax) = 1 ! for QV/SS forcing
  integer ::  isslcl(nfmax) = 1 ! for QV/SS forcing (use LCL as lower boundary of forcing)
  
  integer ::  iforce = 0
  
  integer ::  iforcetyp = 1
  integer ::  nwfor     = 0  ! number of W forcing regions (max of nfmax)
  integer ::  nssfor    = 0  ! number of QV/SS forcing regions (max of nfmax)
  integer ::  iwforce   = 0  ! forcing flag for W
  integer ::  issforce  = 0  ! forcing for QV (supersaturation)
  integer ::  iqrforce  = 0  ! forcing for QR
  integer ::  itopforce = 0
  integer ::  itopforcetype = 1 ! 1 = rain, 2 = graupel, etc., 3 = graupel as 100% rain
  integer ::  itopforcerainopt = 1 ! rain options (binforce) 1 = default, 2 = MP using topforcen0 and topforcerainrate
  real    ::  topforceq = 0.0
  real    ::  topforceN = 0.0
  real    ::  topforceX = 0.0 ! xlocation
  real    ::  topforceY = 0.0 ! ylocation
  real    ::  topforceR = 4000.0 ! radius of forcing region (meters)
  real    ::  topforcedia = 1.5e-3 ! mean particle diameter (meters)
  real    ::  topforcen0  = 8.0e6 ! Initial intercept parameter of distribution (-1 to use topforcedia instead)
  real    ::  topforcerainrate  = 54.0 ! rainrate to use with MP rain distribution
  real    ::  topforcealpha = 0    ! Initial shape parameter

  real    ::  qrfmeso(nfmax) = 0.0 ! value of QR for forcing
  real    ::  qhfmeso(nfmax) = 0.0 ! value of QH for forcing
  real    ::  ssfmeso(nfmax) = 0.0 ! value of QV/SS for forcing
  real    ::  xssfcen(nfmax) = 0.0 ! x-center of forcing
  real    ::  yssfcen(nfmax) = 0.0 ! y-center of forcing

  real    ::  rssfrad(nfmax) = 0.0 ! radius of cylindrical forcing region
  
  real    ::  wfmeso(nfmax) = 0.0 ! value of W for forcing
  real    ::  ufmeso(nfmax) = 0.1 ! value of U for forcing 
  real    ::  umaxmeso(nfmax) = 60.0 ! max value of U when forcing is turned off
  real    ::  xwfcen(nfmax) = 0.0 ! x-center of forcing
  real    ::  ywfcen(nfmax) = 0.0 ! y-center of forcing
  real    ::  zwfcen(nfmax) = 0.0 ! z-center of forcing

  real    ::  xwfcennew(nfmax) = 0.0 ! x-center of forcing
  real    ::  ywfcennew(nfmax) = 0.0 ! y-center of forcing
  real    ::  zwfcennew(nfmax) = 0.0 ! z-center of forcing

  real    ::  xwfrad(nfmax) = 0.0 ! radii of ellipsoidal forcing region
  real    ::  ywfrad(nfmax) = 0.0 
  real    ::  zwfrad(nfmax) = 0.0 

  real    ::  xwfmov(nfmax) = 0.0 ! movement rate of forcing region
  real    ::  ywfmov(nfmax) = 0.0 
  real    ::  zwfmov(nfmax) = 0.0 

  integer ::  twstrt(nfmax) = 0   ! starting time of W forcing
  integer ::  twstop(nfmax) = 0   ! ending time of W forcing

  integer ::  tsstrt(nfmax) = 0   ! starting time of QV/SS forcing
  integer ::  tsstop(nfmax) = 0   ! ending time of QV/SS forcing

  real ::  tsslow(nfmax) = 253.15   ! low (i.e., top) temperature of QV/SS forcing
  real ::  tsshigh(nfmax) = 285.15   ! high (i.e., bottom) temperature of QV/SS forcing
  
  integer :: nphys = 1   ! number of substeps for microphysics
  integer :: igamrain = 1  ! 1 = gamma of volume, 2 = gamma of diameter
  integer :: igamsnow = 1  ! 1 = gamma of volume, 2 = gamma of diameter
  
  real :: chargeperparticle = 100.e-15
  


! Parker (2008) low level cooling
  integer              :: low_level_cooling_flag = 0 ! 0 = off; 1 = on
  real                 :: low_level_cooling_depth = 1000. ! meters
  real                 :: low_level_cooling_rate  = 3.0 ! degrees per hour
  integer              :: low_level_cooling_start = 10800
  integer              :: low_level_cooling_end   = 21600

 
  real, allocatable :: tt13sl(:,:) ! U flux
  real, allocatable :: tt23sl(:,:) ! V flux
  real, allocatable :: td13sl(:,:) ! deformation term
  real, allocatable :: td23sl(:,:) ! deformation term
  real, allocatable :: twt3sl(:,:) ! TKE buoyancy term
  real, allocatable :: thf3sl(:,:) ! TH flux
  real, allocatable :: tmf3sl(:,:) ! QV flux
  real, allocatable :: tdh3sl(:,:) ! not used


  CONTAINS
  
!-----------------------------------------------------------------------------
   SUBROUTINE FORCE_INIT
   RETURN
   END SUBROUTINE FORCE_INIT
   
!-----------------------------------------------------------------------------
  
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
!   /////////////////////          BEGIN           \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\    SUBROUTINE WFORCE      ////////////////////
!     general forcing functions
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
   SUBROUTINE WFORCE(xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,  &
                     nx,ny,nz,dt,w,wt,fw,u,ut,fu,t0,             &
                     gxt,gyt,gzt,loop,time,time_real,    &
                     uinit,vinit,ugrid,vgrid)

!   USE PARAM_MODULE, only: RKSCHEME,ng,pii,luno,bcy
#ifdef MPI
!   USE COMMASMPI_MODULE
#endif
      use input ! ib,ie, etc.
      use constants

!-----------------------------------------------------------------------------
   implicit none

   real, dimension(ib:ie)   :: xh,rxh,uh,ruh
   real, dimension(ib:ie+1) :: xf,rxf,uf,ruf
   real, dimension(jb:je)   :: yh,vh,rvh
   real, dimension(jb:je+1) :: yf,vf,rvf

   integer :: nx,ny,nz
   real    :: dt
   real    :: w (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! for Naylor forcing
   real    :: wt(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! for Naylor forcing
   real    :: fw(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! w tendency (for updraft acceleration term)
   real    :: u (-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! for u forcing
   real    :: ut(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng)
   real    :: fu(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! U tendency
   real    :: t0(-ng+1:nx+ng,-ng+1:ny+ng,-ng+1:nz+ng) ! work array
   real    :: gxt(-ng+1:nx+ng,4), gyt(-ng+1:ny+ng,4), gzt(-ng+1:nz+ng,4) ! axis grid lengths
   real    :: uinit(-ng+1:nz+ng), vinit(-ng+1:nz+ng) ! base state u and v 
   real    :: ugrid, vgrid ! grid motion vector
   integer :: loop
   integer :: time
   real    :: time_real
   
   integer :: n
   integer :: i,j,k
   real    :: x,y,z,radius
   real    :: wtot, u1, v1
   logical :: work_to_do
   integer :: idofor(50)
   real    :: gamma

!-----------------------------------------------------------------------------
! MPI LOCAL VARIABLES
#ifdef MPI
   integer :: ixb, jyb, kzb
   integer :: ixe, jye, kze
#endif

!-----------------------------------------------------------------------------

   work_to_do = .false.
   idofor(:)  = 0
   
   
   DO n = 1,nwfor


    IF ( time_real .ge. twstrt(n) .and. time_real .le. twstop(n) ) THEN
      work_to_do = .true.
      idofor(n) = 1
    ENDIF

! set movement to mean layer wind
      IF ( xwfmov(n) .le. -999. ) THEN
       wtot = 0.0
       u1 = 0.0
       v1 = 0.0

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg
    
       do k = kzb,kze
#else
       DO k = 1,nz-1
#endif
        z = gzt(k,1)

         radius = sqrt(((z-zwfcen(n))/zwfrad(n))**2)
      
       IF ( radius .lt. 1.0 ) THEN
       
        IF ( xwfmov(n) .eq. -999. ) THEN       ! weighted mean wind in forcing depth
         u1 = u1 + (uinit(k) - ugrid)*(cos(pii*radius)+1.)/2./gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)*(cos(pii*radius)+1.)/2./gzt(k,3)
         wtot = wtot + (cos(pii*radius)+1.)/2./gzt(k,3)
        ELSEIF ( xwfmov(n) .eq. -1000. ) THEN  ! mean wind in forcing depth
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
        ENDIF
       ELSEIF ( xwfmov(n) .eq. -1001. ) THEN ! mean wind through depth of domain
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1002. .and. z .lt. 10.e3 ) THEN ! mean wind through 10 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1002. .and. z .lt. 6.e3 ) THEN ! mean wind through 6 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ELSEIF ( xwfmov(n) .eq. -1003. .and. z .gt. 5.e3 .and. z .lt. 10.e3 ) THEN ! mean wind 5-10 km
         u1 = u1 + (uinit(k) - ugrid)/gzt(k,3)
         v1 = v1 + (vinit(k) - vgrid)/gzt(k,3)
         wtot = wtot + 1./gzt(k,3)
       ENDIF
       
      ENDDO

       xwfmov(n) = u1/wtot
       ywfmov(n) = v1/wtot
       
       write(luno,*) 'wforce: n, xwfmov, ywfmov = ',n,xwfmov(n),ywfmov(n)
      
      ENDIF

       IF ( Abs(xwfmov(n)) .lt. 900. ) THEN
        xwfcennew(n) = xwfcen(n) + Max(0.0, time_real - twstrt(n))*xwfmov(n)
!        print*, 'Update xwfcen(',n,') = ',xwfcen(n),dt,xwfmov(n)
       ELSE
         xwfcennew(n) = xwfcen(n) 
       ENDIF
       IF ( Abs(ywfmov(n)) .lt. 900. ) THEN
        ywfcennew(n) = ywfcen(n) + Max(0.0, time_real - twstrt(n))*ywfmov(n)
       ELSE
        ywfcennew(n) = ywfcen(n) 
       ENDIF
       IF ( Abs(zwfmov(n)) .lt. 900. ) THEN
        zwfcennew(n) = zwfcen(n) ! + Max(0.0, time_real - twstrt(n))*zwfmov(n)
       ELSE
        zwfcennew(n) = zwfcen(n) 
       ENDIF
      
   ENDDO
   

!    IF ( .not. work_to_do )  RETURN
    IF ( work_to_do ) THEN ! RETURN
    
    t0(:,:,:) = 0.0
!
!  compute forcing based on shape I (cosine squared function)
!
   DO n = 1,nwfor
     IF ( abs(iwshap(n)) .eq. 1 ) THEN
      
!     IF ( time .ge. twstrt(n) .and. time .le. twstop(n) ) THEN
      IF ( idofor(n) .eq. 1 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        z = gzt(k,1)
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)
#else
       DO k = 1,nz-1
        z = gzt(k,1)
        DO j = 1,ny-1
         y = gyt(j,1)
         DO i = 1,nx-1
          x = gxt(i,1)
#endif
         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((y-ywfcennew(n))/ywfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ELSE
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ENDIF

       IF ( radius .lt. 1.0 ) THEN
!       IF ( iwshap(n) .lt. 0 .and. dist .le. 1.0 ) THEN
        IF ( iwshap(n) .eq. 1 ) THEN
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
         IF ( wfmeso(n) .gt. 0.0 ) THEN
           t0(i,j,k) = Max(t0(i,j,k), wfmeso(n)*(cos(pii*radius)+1.)/2.)
         ELSE
           t0(i,j,k) = wfmeso(n)*(cos(pii*radius)+1.)/2.
         ENDIF
        ELSEIF ( iwshap(n) .eq. -1 .and. w(i,j,k) .ge. -0.1 ) THEN
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
          t0(i,j,k) = Max(t0(i,j,k), wfmeso(n)*(cos(pii*radius)+1.)/2.)
        ENDIF
!       ENDIF
       ENDIF

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

     ELSEIF ( abs(iwshap(n)) .eq. 2 ) THEN ! nudge toward 10 m/s (Naylor et al. 2012, MWR)
      
!     IF ( time .ge. twstrt(n) .and. time .le. twstop(n) ) THEN
      IF ( idofor(n) .eq. 1 ) THEN

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        z = gzt(k,1)
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)
#else
       DO k = 1,nz-1
        z = gzt(k,1)
        DO j = 1,ny-1
         y = gyt(j,1)
         DO i = 1,nx-1
          x = gxt(i,1)
#endif
         IF ( ny .gt. 2 ) THEN
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((y-ywfcennew(n))/ywfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ELSE
         radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
         ENDIF

       IF ( radius .lt. 1.0 ) THEN
!       IF ( iwshap(n) .lt. 0 .and. dist .le. 1.0 ) THEN
        IF ( iwshap(n) .eq. 2 ) THEN
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
           t0(i,j,k) = wfmeso(n)*Max(0.0, 10.0*(cos(pii*radius)+1.)/2. - wt(i,j,k) )
        ELSEIF ( iwshap(n) .eq. -2  ) THEN
           t0(i,j,k) =(wfmeso(n))*Max(0.0, 10.0*(cos(pii*radius)+1.)/2. - w(i,j,k) )
!         fw(i,j,k) = fw(i,j,k) + wfmeso(n)*(cos(pii*radius)+1.)/2.
!          t0(i,j,k) = Max(t0(i,j,k), wfmeso(n)*(cos(pii*radius)+1.)/2.)
        ENDIF
!       ENDIF
       ENDIF

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF
 
      ELSEIF ( abs(iwshap(n)) .eq. 3 ) THEN ! nudge U (for workshop squall line case)
      
!     IF ( time .ge. twstrt(n) .and. time .le. twstop(n) ) THEN
      IF ( idofor(n) .eq. 1 ) THEN

      gamma = 1.0
      IF (time_real >= twstop(n)-300.) THEN
        gamma = 1.0+(0.0-1.0)*( time_real - (twstop(n)-300.) )/(300.)
      ENDIF

#ifdef MPI
       kzb = -ng+1
       kze = ktile+ng
       if (kzbeg .eq. nzbeg) kzb = 1
       if (kzend .eq. nzend) kze = kzend-kzbeg

       jyb = -ng+1
       jye = jtile+ng
       if (jybeg .eq. nybeg) jyb = 1
       if (jyend .eq. nyend) jye = jyend-jybeg

       ixb = -ng+1
       ixe = itile+ng
       if (ixbeg .eq. nxbeg) ixb = 1
       if (ixend .eq. nxend) ixe = ixend-ixbeg

       do k = kzb,kze
        z = gzt(k,1)
        
        do j = jyb,jye
        y = gyt(j,1)
         
         do i = ixb,ixe
         x = gxt(i,1)
#else
       DO k = 1,nz-1
        z = gzt(k,1)
        DO j = 1,ny-1
         y = gyt(j,1)
         DO i = 1,nx-1
          x = gxt(i,1)
#endif
!         IF ( ny .gt. 2 ) THEN
!         radius = sqrt(((x-xwfcen(n))/xwfrad(n))**2  &
!                      +((y-ywfcen(n))/ywfrad(n))**2  &
!                      +((z-zwfcen(n))/zwfrad(n))**2)
!         ELSE
! 2D forcing -- assumes periodic domain or channel simulation
          radius = sqrt(((x-xwfcennew(n))/xwfrad(n))**2  &
                      +((z-zwfcen(n))/zwfrad(n))**2)
!         ENDIF

!        IF ( Abs(x-xwfcennew(n)) < xwfrad(n) .and. Abs(z - zwfcen(n)) < zwfrad(n) ) THEN
        IF ( ( bcy == 2 .and. Abs(x-xwfcennew(n)) < xwfrad(n) .and.  z < zwfrad(n) ) .or. &
     &       ( Abs(y-ywfcennew(n))/ywfrad(n) <= 1. .and. Abs(x-xwfcennew(n)) < xwfrad(n) .and.  z < zwfrad(n) )   ) THEN
!        IF ( radius < 1. ) THEN
!          uten1(i,j,k)=uten1(i,j,k)+0.10*gamma           &
!                       *cos(0.5*pi*(xh(i)-0.0)/10000.0)  &
!                       *((cosh(2.5*(zh(i,j,k)-0.0)/10000.0))**(-2))

         IF ( abs( u(i,j,k) ) <= umaxmeso(n) ) THEN
         fu(i,j,k) = fu(i,j,k) + ufmeso(n)*gamma           &
                       *cos(0.5*pii*(x-xwfcennew(n))/xwfrad(n))  &
                       *((cosh(2.5*z/zwfrad(n)))**(-2))
!                       *((cosh(2.5*(z-zwfcen(n))/zwfrad(n)))**(-2))
         ENDIF
        ENDIF

         ENDDO
        ENDDO
       ENDDO
      
      ENDIF

     ENDIF ! iwshap
   
   ENDDO

#ifdef MPI
    kzb = -ng+1
    kze = ktile+ng
    if (kzbeg .eq. nzbeg) kzb = 1
    if (kzend .eq. nzend) kze = kzend-kzbeg

    jyb = -ng+1
    jye = jtile+ng
    if (jybeg .eq. nybeg) jyb = 1
    if (jyend .eq. nyend) jye = jyend-jybeg

    ixb = -ng+1
    ixe = itile+ng
    if (ixbeg .eq. nxbeg) ixb = 1
    if (ixend .eq. nxend) ixe = ixend-ixbeg

    do k = kzb,kze ; do j = jyb,jye ; do i = ixb,ixe
#else
    DO k = 1,nz-1
     DO j = 1,ny-1
      DO i = 1,nx-1
#endif
        fw(i,j,k) = fw(i,j,k) + t0(i,j,k)
        
      ENDDO
     ENDDO
    ENDDO
    
    ENDIF ! work_to_do

!   DO n = 1,nwfor
!    IF ( loop .eq. RKSCHEME ) THEN
!       IF ( Abs(xwfmov(n)) .lt. 900. ) THEN
!        xwfcen(n) = xwfcen(n) + dt*xwfmov(n)
!!        print*, 'Update xwfcen(',n,') = ',xwfcen(n),dt,xwfmov(n)
!       ENDIF
!       IF ( Abs(ywfmov(n)) .lt. 900. ) THEN
!        ywfcen(n) = ywfcen(n) + dt*ywfmov(n)
!       ENDIF
!       IF ( Abs(zwfmov(n)) .lt. 900. ) THEN
!        zwfcen(n) = zwfcen(n) + dt*zwfmov(n)
!       ENDIF
!     ENDIF
!   ENDDO
   
   RETURN
   
   END SUBROUTINE WFORCE

   
  END MODULE FORCE_MODULE
