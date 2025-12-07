
module particleMOD
    use generalMOD,only: dp

    use constants
    implicit None


    ! particle scales
    real(dp), parameter:: scale_mp = M_sun      ! Msun
    real(dp), parameter:: scale_tp = Myr2sec    ! Myr
    real(dp), parameter:: scale_vp = kms      ! kms
    real(dp)           :: pscale_m = M_sun      ! Msun
    real(dp)           :: pscale_t = Myr2sec    ! Myr
    real(dp)           :: pscale_v = kms        ! kms
    real(dp)           :: pscale_l = pc2cm

endmodule


module feedbackMOD
    use generalMOD,only: dp
    use constants, only: Myr2sec, yr2sec, kms
    implicit None
    
    ! feedback channels
    character(1000),save:: path_to_stellar_evolutionary_tracks = "/home/mattia/codes/Stevia/data/"

    character(20)     :: wind_model      =  "stellar_evo"      !!!               chose modelling of wind
    real(dp),save     :: mdot_const          =  1.0d-7     !!! [M_sun/yr]    constant mass loss rate M_sun/yr
    real(dp),save     :: wind_terminal_const = 2.0d3  !!! [km/s]        constant terminal velocity of the stellar wind km/s
    
endmodule


module parameters
    use feedbackMOD
    use particleMOD
    use generalMOD
    implicit None

endmodule parameters

