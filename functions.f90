module functions
    Implicit none
contains

function Pk(kk)    ! Power spectrum in terms of wavenumber kk, Amp controls the amplitude 
    use fiducial 
    Implicit none
    Real*8 :: kk, Pk

    Pk = Amp*kk**n*dexp(-kk**2*lambda**2/2.d0)

End function Pk

function P1G(A_k)    ! One-point distribution of a Gaussian random field. 
    use fiducial
    Implicit none
    Real*8 :: P1G,A_k

    P1G = 2.d0*A_k*exp(-A_k**2)

end function P1G

function P4(k_index)    ! Normalised polyspectra for N=4. Equation (45) in published version of Matsubara's paper for k,k,-k,-k
    use arrays
    use fiducial
    Implicit none
    Double Complex :: P4,term1,term2,term3,term4
    Integer*4 :: k_index,i

    term1 = dcmplx(0.d0,0.d0) 

    term2 = dcmplx(0.d0,0.d0)

    term3 = dcmplx(0.d0,0.d0)

    term4 = dcmplx(0.d0,0.d0)

    If ((k_index .eq. 1) .or. (k_index .eq. Nyquist)) then

        Do i=1,number_of_realizations

            term1 = All_alphak(i,k_index)**4 + term1

            term2 = All_alphak(i,k_index) + term2 

            term4 = All_alphak(i,k_index)**2 + term4

        End Do

        term1 = term1/dble(number_of_realizations)

        term2 = term2/dble(number_of_realizations)

        term4 = term4/dble(number_of_realizations)

        P4 = Volume*(term1 + 3.d0*term2**4 - 4.d0*term2**2*term4)

    Else

        Do i=1,number_of_realizations

            term1 = All_alphak(i,k_index)**2*All_alphak(i,number_of_sites - k_index + 2)**2 + term1

            term2 = All_alphak(i,k_index) + term2 

            term3 = All_alphak(i,number_of_sites - k_index + 2) + term3

            term4 = All_alphak(i,k_index)*All_alphak(i,number_of_sites - k_index + 2) + term4

        End Do

        term1 = term1/dble(number_of_realizations)

        term2 = term2/dble(number_of_realizations)

        term3 = term3/dble(number_of_realizations)

        term4 = term4/dble(number_of_realizations)

        P4 = Volume*(term1 + 3.d0*term2**2*term3**2 - 4.d0*term2*term3*term4)
    
    End If

end function P4

function P1NG(A_k,k_index)    ! One-point distribution of a non-Gaussian random field. 
    use fiducial
    use arrays
    Implicit none
    Real*8 :: P1NG,A_k
    Integer*4 :: k_index!,realization

!    P1NG = (1.d0 + (All_Ak(realization,k_index)**4/4.d0 - All_Ak(realization,k_index)**2 &
!    + 1.d0/2.d0)*real(P4(k_index))/Volume)*P1G(All_Ak(realization,k_index))

    P1NG = (1.d0 + (A_k**4/4.d0 - A_k**2 + 1.d0/2.d0)*real(P4(k_index))/Volume)*P1G(A_k)

end function P1NG


!function B(k1, k2, Amp)
    ! Returns the bispectrum for local quadratic nongaussianity
    ! Inputs: ki = wavenumbers, Amp=power spectrum amplitude.
!    Use Stuff
!    Implicit none
!    Real*8 :: k1, k2, P1, P2, P3, Amp, B

!    P1 = Pk(k1,Amp)
!    P2 = Pk(k2,Amp)
!    P3 = Pk(-k1-k2,Amp)
!    B = 2d0*(P1*P2+P1*P3+P2*P3)/sqrt(P1*P2*P3)

!end function B

!function Pk4(k1,k2,k3,k4,Amp)
    ! This is an auxiliary function for the local quadratic nongaussian trispectrum.
    ! See Smith & Kamionkowski (2012) arxiv:1203.6654 equation 9.
    ! Inputs: ki = wavenumbers, Amp=power spectrum amplitude.
!    Use Stuff
!    Implicit none
!    Real*8 :: k1,k2,k3,k4,Amp,P1,P2,P3,P4,Pk4

!    P1 = Pk(k1,Amp)
!    P2 = Pk(k2,Amp)
!    P3 = Pk(k3,Amp)
!    P4 = Pk(k4,Amp)
!    Pk4 = 4d0*Pk(k1+k2,Amp)*(P1*P3+P1*P4+P2*P3+P2*P4)

!end function Pk4

!function Trisp(k1,k2,k3,k4,Amp)
    ! This is the local quadratic nongaussian trispectrum.
    ! See Smith & Kamionkowski (2012) arxiv:1203.6654 equation 8.
    ! Inputs: ki = wavenumbers, Amp=power spectrum amplitude.
!    Use Stuff
!    Implicit none
!    Real*8 :: k1,k2,k3,k4,Amp,Trisp

!    Trisp = (Pk4(k1,k2,k3,k4,Amp) + Pk4(k1,k3,k2,k4,Amp) + Pk4(k1,k4,k2,k3,Amp))/sqrt(Pk(k1,Amp)*Pk(k2,Amp)*Pk(k3,Amp)*Pk(k4,Amp))

!end function Trisp

function rndgauss() ! Generates a random gaussian variate, zero mean, unit variance
    Implicit none
!    Complex ( kind = 8 ) c8_normal_01
    Real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
!    Real ( kind = 8 ) r8_uniform_01
!    Integer ( kind = 4 ) seed
    Real ( kind = 8 ) v1
    Real ( kind = 8 ) v2
    Real ( kind = 8 ) x_r
    Real*4 x
    Real*8 :: rndgauss
    call random_number(x)
    v1 = dble(x)
    call random_number(x)
    v2 = dble(x)
    x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
    rndgauss = x_r
    return
end function rndgauss

end module functions 
