module stellar_evolution
  use parameters, only: dp
  use constants , only:  Myr2sec, sigma_boltz, L_sun, M_sun, yr2sec, pi, factG_in_cgs
  use feedbackMOD, only: mdot_const,wind_terminal_const, wind_model, path_to_stellar_evolutionary_tracks
  !! now we first load all the data from file with a subroutine
  implicit none
  public :: interp
  interface interp
    module procedure load_stellar_data, evolve_star, SN_time!, choose_model!, load_stellar_data!choose_model, star_status, load_stellar_data, 
  end interface
  interface winds
    module procedure C_TEMP,BISTABILITY_MDOT 
  end interface
  interface utils
    module procedure random_stduniform,random_uniform
  end interface
  integer, parameter, private:: nmodels = 24, nlines=400 
  real(dp),dimension(nmodels,nlines) , private, save:: star_age,star_mass,star_logL,star_logTe,star_logMdot, star_gedd
  real(dp),dimension(nmodels,nlines) , private, save:: star_X,star_Y,star_C12,star_N14,star_O16
  real(dp),dimension(nmodels,nlines) , private, save:: star_Al26
  real(dp),dimension(nmodels)        , private, save:: star_mass_init,star_final_age
  real(dp),dimension(nlines)         , private, save:: normage
  character(13),private:: type_dummy
  logical, private,save:: first=.true., rotation=.true.
  integer, dimension(nlines), private, save:: ok_lines 
  integer, dimension(nmodels), private, save:: models_indx
  real(dp), parameter, public:: stellar_scale_t = Myr2sec
  real(dp), parameter:: kms = 1.0d5




contains 
!!!! LEGACY SUBROUTINES FOR STELLAR EVOLUTION INTERPOLATION
!  subroutine interp_evo(m_ini, tt, temp3, lum3, mdot3, gamma_ed3,m3)
!   use parameters, only: dp

!   implicit none
!   real(dp):: m_ini
  
!   integer:: closer1, closer2, dummy_int
!   real(dp):: tt!mass=0., teff=0., lum=0.,mdot=0., gamma_ed=0.,
!   real(dp):: m1=0., m2=0., temp1=0,temp2=0, lum1=0,lum2=0, mdot1=0,mdot2=0, gamma_ed1=0. ,gamma_ed2=0.
!   real(dp):: m3, temp3, lum3, mdot3, gamma_ed3


!   if (first) then 
!        !load data from file 
       
!        call load_stellar_data
!   endif

  
! endsubroutine interp_evo

!> Return the supernova time for a star of initial mass m_ini by
!> interpolating the pre-tabulated final ages.
subroutine SN_time(m_ini, t_sn)
  use parameters, only: dp
  use interpolation

  implicit none
  real(dp), intent(in) ::  m_ini
  real(dp), intent(out)::  t_sn
  integer:: i
  if (first) then 
      ! Load stellar evolution tracks on the first invocation.
      call load_stellar_data
  endif

  ! Mass-age interpolation delivers the supernova time in Myr.
  t_sn      = lin_interp_bound(star_mass_init,star_final_age, m_ini )
  endsubroutine

!> Interpolate stellar properties (mass, Teff, luminosity, etc.)
!> at a given fractional age tt/t_sn for the requested initial mass.
subroutine evolve_star(m_ini, tt, mass, teff, lum,mdot, gamma_ed, v_term)
  use parameters, only: dp
  use interpolation

  implicit none
  real(dp), intent(in):: m_ini, tt
  real(dp), intent(out):: mass, teff, lum, mdot, gamma_ed, v_term
  integer:: idx1, idx2
  real(dp):: t_sn, m1=0.0, m2=0.0, teff1=0.0,teff2=0.0, lum1=0,lum2=0, mdot1=0.0,mdot2=0.0, gamma_ed1=0.0&
            & ,gamma_ed2=0.0,stellar_radius=0.0,t_interp=0.0, m_interp=0.0

  !mass to stellar evo units 
  
  ! cap age betweem 1.0 and 0  
  !t_interp=tt
  if (first) then 
      ! Lazy-load tracks so the interpolation tables are available.
      call load_stellar_data
  endif

  ! Work with normalized age and clamp outputs before interpolation.
  call SN_time(m_ini, t_sn)
 ! print*, tt, t_sn
  t_interp=tt/t_sn
  mass=0.
  teff=0.
  mdot=0.
  gamma_ed=0.
  lum=0.
if (trim(wind_model) == "constant") then 
    mdot=mdot_const ! !m_sun/yr
    v_term = wind_terminal_const * kms ! cgs
    !some padding values for relevant stuff
    teff=20000.
    lum=1d36
    
else if (trim(wind_model)=="stellar_evo" .or.trim(wind_model)=="bistability_jump" ) then 
  idx1=locate(star_mass_init,m_ini)

  
  if ((idx1==1 .and. m_ini<=star_mass_init(idx1)) .or. (idx1 == nmodels .and. m_ini>=star_mass_init(idx1))) then 
    idx2=idx1
    
    mass      = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_mass   (idx1 , 1:ok_lines(idx1)), t_interp )
    lum       = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_logL   (idx1 , 1:ok_lines(idx1)), t_interp )
    teff      = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_logTe  (idx1 , 1:ok_lines(idx1)), t_interp )
    mdot      = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_logMdot(idx1 , 1:ok_lines(idx1)), t_interp )
    gamma_ed  = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_Gedd   (idx1 , 1:ok_lines(idx1)), t_interp )

  else
    idx2=idx1+1
    
    m1         = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_mass   (idx1 , 1:ok_lines(idx1)), t_interp )
    lum1       = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_logL   (idx1 , 1:ok_lines(idx1)), t_interp )
    teff1      = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_logTe  (idx1 , 1:ok_lines(idx1)), t_interp )
    mdot1      = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_logMdot(idx1 , 1:ok_lines(idx1)), t_interp )
    gamma_ed1  = lin_interp_bound( star_age(idx1 , 1:ok_lines(idx1)), star_Gedd   (idx1 , 1:ok_lines(idx1)), t_interp )
    
    m2         = lin_interp_bound( star_age(idx2 , 1:ok_lines(idx2)), star_mass   (idx2 , 1:ok_lines(idx2)), t_interp )
    lum2       = lin_interp_bound( star_age(idx2 , 1:ok_lines(idx2)), star_logL   (idx2 , 1:ok_lines(idx2)), t_interp )
    teff2      = lin_interp_bound( star_age(idx2 , 1:ok_lines(idx2)), star_logTe  (idx2 , 1:ok_lines(idx2)), t_interp )
    mdot2      = lin_interp_bound( star_age(idx2 , 1:ok_lines(idx2)), star_logMdot(idx2 , 1:ok_lines(idx2)), t_interp )
    gamma_ed2  = lin_interp_bound( star_age(idx2 , 1:ok_lines(idx2)), star_Gedd   (idx2 , 1:ok_lines(idx2)), t_interp )

    mass     = m1        + (m_ini-star_mass_init(idx1))* (m2        - m1       ) / (star_mass_init(idx2) - star_mass_init(idx1))
    lum      = lum1      + (mass-m1)                   * (lum2      - lum1     ) / (m2-m1)
    teff     = teff1     + (mass-m1)                   * (teff2     - teff1    ) / (m2-m1)
    mdot     = mdot1     + (mass-m1)                   * (mdot2     - mdot1    ) / (m2-m1)
    gamma_ed = gamma_ed1 + (mass-m1)                   * (gamma_ed2 - gamma_ed1) / (m2-m1)

  endif

  ! Convert back from log-space and derive terminal velocity from the
  ! stellar radius and escape speed scaling.
  lum = 10**lum
  teff= 10**teff
  mdot= 10**mdot
  stellar_radius=(1./teff**2.)*sqrt(lum*L_sun/(4.*pi*sigma_boltz))! cm 
  v_term =sqrt(2.*factG_in_cgs*(mass*M_sun)/stellar_radius*(1.0-gamma_ed))*C_TEMP(teff) !!! cm/s

  if (trim(wind_model)=="bistability_jump") then             
    mdot = BISTABILITY_MDOT(mass, lum, teff, gamma_ed, stellar_radius, v_term)
  endif
  ! to user_units
  !v_term=v_term/kms
endif

end subroutine evolve_star


!> Populate the stellar evolution tables in memory by parsing the
!> Geneva tracks for all available models. Subsequent calls no-op.
subroutine load_stellar_data
  use parameters, only: dp 
  use particleMOD
  use feedbackMOD
   implicit none
   character(400)::filename
   character(1):: type_checker ="r",Rot_dummy
   integer :: i, j,  ilun, imodel, nlines_counter
   integer ::  Line_dummy
   real(dp):: m_ini_dummy,Zini_dummy,Time_dummy,Mass_dummy,logL_dummy,logTe_dummy,final_age_dummy
   real(dp):: X_dummy,Y_dummy,C12_dummy,C13_dummy,N14_dummy,O16_dummy,O17_dummy,O18_dummy,Ne20_dumm,Ne22_dummy,Al26_dummy
   real(dp):: QCC_dummy,logTe_u_dummy,logdM_dt_dummy,log_rhoc_dummy,logTc_dummmy
   real(dp):: Xc_dummy ,Yc_dummy ,C12c_dummy,C13c_dumm,N14c_dummy,O16c_dummy,O17c_dummy,O18c_dummy,Ne20c_dummy, Ne22c_dummy,Al26c_dummy
   real(dp):: Omegas_dummy,Omegac_dummy,oblat_dummy,dM_dtR_dummy,vcrit1_dummy,vcrit2_dummy,veq_dummy,OOc_dummy,Gedd_dummy,dM_dtm_dummy,Ltot_dummy
   real(dp):: ttest(500), mtest(100)
   real(dp):: dummy 
   integer:: dummy_int
   ilun=283485
   if (.not.first ) return
   do dummy_int = 1 ,nmodels
    models_indx(dummy_int)=dummy_int
   enddo
   first=.False.
   filename=TRIM(path_to_stellar_evolutionary_tracks)//'/J_A+A_537_A146/tables.dat' 
   !if (myid==1) write(*,*) "Loading stellar data from "//filename
   open(ilun,file=filename, action='read')

   imodel=0
   i=1
   ! Read the ASCII table sequentially, storing each evolutionary line.
   do while (i<= nlines*nmodels*2) 
     read(ilun,*, end=9999) m_ini_dummy,Zini_dummy,Rot_dummy,Line_dummy,Time_dummy,Mass_dummy,logL_dummy,logTe_dummy,&
     & X_dummy,Y_dummy,C12_dummy,C13_dummy,N14_dummy,O16_dummy,O17_dummy,O18_dummy,Ne20_dumm,Ne22_dummy,Al26_dummy,&
     & QCC_dummy,logTe_u_dummy,logdM_dt_dummy,log_rhoc_dummy,logTc_dummmy,&
     & Xc_dummy ,Yc_dummy ,C12c_dummy,C13c_dumm,N14c_dummy,O16c_dummy,O17c_dummy,O18c_dummy,Ne20c_dummy, Ne22c_dummy,Al26c_dummy,&
     & Omegas_dummy,Omegac_dummy,oblat_dummy,dM_dtR_dummy,vcrit1_dummy,vcrit2_dummy,veq_dummy,OOc_dummy,Gedd_dummy,dM_dtm_dummy,Ltot_dummy
     if (type_checker == trim(Rot_dummy)) then 
       j= Line_dummy
       if (Line_dummy==1) imodel=imodel+1
       star_mass_init(imodel) = m_ini_dummy
       star_age    (imodel,j) = Time_dummy/1.0d6
       star_mass   (imodel,j) = Mass_dummy
       star_logL   (imodel,j) = logL_dummy
       star_logTe  (imodel,j) = logTe_dummy
       !! set the zeros in the data appropriately for log conversion, otherwise they are converted to ones. 
       if (logdM_dt_dummy>-1.0) then 
         star_logMdot(imodel,j) = -13.522878745 ! this is just the solar wind. 
       else 
        star_logMdot(imodel,j) = logdM_dt_dummy
       endif 
       star_Gedd   (imodel,j) = Gedd_dummy
       star_X      (imodel,j) = X_dummy
       star_Y      (imodel,j) = Y_dummy
       star_C12    (imodel,j) = C12_dummy
       star_N14    (imodel,j) = N14_dummy
       star_O16    (imodel,j) = O16_dummy
       star_Al26   (imodel,j) = Al26_dummy
       !if (Line_dummy==400) print*, imodel, i,nlines * nmodels, Line_dummy, m_ini_dummy, Rot_dummy,Gedd_dummy, myid
       if (Line_dummy==400) star_final_age(imodel)=Time_dummy/1.0d6
     endif
     i=i+1
   enddo
 9999 continue!if (myid==1)  write(*,*) "Stellar evolution data loaded"
 !777  format(F6.2,13X,E22.15,1X,F11.6,F10.6,1X,F9.6,3(1X,E14.7),16X,E14.7,1X,E14.7,61X,E10.3,19X,F8.3)
 !888  format(F6.2,F5.3,A1,I3,E22.15,F11.6,F10.6,F9.6,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,&
 !& E14.7,E14.7,E10.3,F7.4,F9.6,F8.3,F9.6,F9.6,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,E14.7,&
 !& E10.3,E10.3,E10.3,E10.3,E10.3,E9.2,E9.2,E9.2,F9.6,F9.6,E14.7,E17.10)
  
 ! normalize ages of the models to the final age

  ! Normalize the age axis so all tracks span 0-1 across their lifetimes.
  do imodel=1,nmodels
    star_age(imodel,1:nlines) = star_age(imodel,1:nlines) / star_age(imodel,nlines)
  enddo

 ! keep only the part of the array where the star keeps evolving

  do imodel=1,nmodels
    ok_lines(imodel) = 1 
    do j=1,nlines

        if(star_age(imodel,j) == star_age(imodel,nlines)   ) then !j+1) ) then
           
          !print*, star_age(imodel,j), star_age(imodel,nlines), imodel,ok_lines(imodel)
          exit
        endif 
        ok_lines(imodel) = ok_lines(imodel) + 1
    end do
    !print*, "last and current t and mdot", star_mass_init(imodel),star_mass(imodel,ok_lines(imodel)),&
    !& star_age(imodel,ok_lines(imodel)), star_logMdot(imodel,ok_lines(imodel)), ok_lines(imodel)
  end do
  !do i = 1, nmodels
   ! write(*,*) "jjh"
    !print*, ok_lines(i)
    !do j=1,ok_lines(i) !(star_mass(i,j), j=1,nlines)
    !print*, star_logMdot(i,ok_lines(i)),star_final_age, star_mass_init(i)
    !enddo
 ! enddo 
  ! plot the fitted tracks
  !do i= 1,size(mtest) 
  !  mtest(i)=i
    !call random_uniform(0.9_dp,119.0_dp, mtest(i))
  !enddo
  !ttest(1)=1.0/size(ttest)
  !do i=2,size(ttest)
  !  ttest(i)=ttest(i-1)+1.0/size(ttest)
  !enddo
  !do i = 1,size(mtest) 
  !  do j= 1,size(ttest)
  !    write(*,*),i
      !    call evolve_star(mtest(i),ttest(j), Mass_dummy, logTe_dummy, logL_dummy,logdM_dt_dummy, Gedd_dummy,logTe_u_dummy)
      
  !    write(i+1000,*)ttest(j),mtest(i),Mass_dummy, logTe_dummy, logL_dummy,logdM_dt_dummy, Gedd_dummy, logTe_u_dummy/kms
  !  enddo
    
  !enddo 
    
endsubroutine
!> Produce a reproducible double-precision uniform deviate in (0,1].
subroutine random_stduniform(u)
  implicit none
  real(dp),intent(out) :: u
  real :: r
  call random_number(r)
  u = 1 - r
end subroutine random_stduniform
!> Sample a uniform deviate x in [a,b] using random_stduniform.
subroutine random_uniform(a,b,x)
  implicit none
  real(dp),intent(in) :: a,b
  real(dp),intent(out) :: x
  real(dp) :: u
  call random_stduniform(u)
  x = (b-a)*u + a
end subroutine random_uniform
!> Piecewise mapping between temperature and v_term / v_esc scaling.
function C_TEMP(temp) result(v_term_to_v_esc)
  use parameters, only: dp
  implicit none
  real(dp) , intent(in)  :: temp
  real(dp)               :: v_term_to_v_esc
  !print*, temp
 if (temp .le. 10000.) then 

     v_term_to_v_esc=1.

 else if ((temp .gt. 10000.) .and. (temp .le. 21000.) ) then 

     v_term_to_v_esc=1.4

 else if (temp .gt. 21000) then 
     v_term_to_v_esc=2.65

 endif
end function C_TEMP

!> Empirical mass-loss recipe across the bistability jump; switches
!> between the cool and hot wind branches based on effective temperature.
function BISTABILITY_MDOT(mass, lum, temp, gamma_ed, stellar_radius,v_term) result(mdot)
  !!!!!ATTENTION, THIS METOD DO NOT GIVE RESULTS FOR T<12kK
  use constants, only: pi, M_sun, R_sun,yr2sec
  use parameters, only: dp
  implicit none

  real(dp) , intent(in)  :: mass, lum, temp, gamma_ed, stellar_radius, v_term
  real(dp)             :: mdot
  real(dp)             :: rho_jump, temp_jump, jump1, jump2
  jump1=22500.
  jump2=27500.
  !print*, temp
  if (temp .gt. 22500. .and. temp .le. 27500 ) then
    rho_jump=10.**(-14.94+3.2*gamma_ed)
    temp_jump=1000*(49.1+1.64*log10(rho_jump))
    if (temp .le. temp_jump) then 
        jump1=temp_jump
    else if (temp .gt. temp_jump) then 
        jump2=temp_jump
    endif
    endif
  if (temp .le. jump1) then !!actually woorks only until 12000 K 
    mdot= 10.**(                           & 
          -6.688                           &
          +2.210*LOG10(lum/1.0d5)          &
          -1.339*log10(mass/30.)           &
          -1.601*log10(C_TEMP(temp)/2.)    &
          +1.070*log10(temp/20000.))

  else if (temp .gt. jump2) then !.and. (temp .le. 50000.) ) then !!!lets allow it... 
    mdot= 10.**(                           &
          -6.697                           &
          +2.194*LOG10(lum/1.0d5)          &
          -1.313*log10(mass/30.)           &
          -1.226*log10(C_TEMP(temp)/2.)    &
          +0.933*log10(temp/40000.)        &
          -10.92*(log10(temp/40000.))**2)

  endif 
  !print*,mass,  mdot

end function BISTABILITY_MDOT
endmodule