Module fiducial

Implicit none

save 

Real*8,parameter :: Amp = 1.d0 ! Amplitude of the power spectrum 
Real*8,parameter :: n = 0.d0 ! Spectral index in the power spectrum 
Real*8,parameter :: lambda_factor = 3.d-2 ! constant in smoothing scale in the power spectrum  
Real*8,parameter :: Pi = 3.141592653589793d0
Real*8,parameter :: fundamental_wavenumber = 1.d0 !2.d0*Pi/L ! lower virial mass fiducial value in Hill and Spergel paper 1312.4525
Real*8,parameter :: L = 2.d0*Pi/fundamental_wavenumber    !1.0d0 ! fundamental length in real space
Real*8,parameter :: lambda = lambda_factor*L ! smoothing scale in the power spectrum 
Real*8,parameter :: fNL = 2.d-1 ! 5.d-1 ! 1.d0 ! Non-Gaussianity parameter for the quadratic mapping. Matsubara uses three cited values
!Real*8,parameter :: H_0 = 69.7d3 ! Hubble parameter at present time in m s-1 Mpc-1
!Real*8,parameter :: h = 69.7d-2 ! Dimensionless 
!Real*8,parameter :: Omega_m0 = 0.282d0 ! Matter parameter density at present time
!Real*8,parameter :: Omega_L0 = 0.7181d0 ! 1 - Omega_m0 : Dark energy parameter density (LCDM model) at present time in a flat universe
!Real*8,parameter :: Omega_b_h2 = 0.02240d0 ! Baryon density today 
!Real*8,parameter :: ns = 0.9646d0 ! scalar spectrum power-law index 
!Real*8,parameter :: z_dec = 1100d0 ! decoupling red-shift
!Real*8,parameter :: A_s = 2.43d-9 ! dimensionless curvature power spectrum at k_0 = 0.05 Mpc-1
!Real*8,parameter :: k_0 = 0.002 ! Mpc-1
!Real*8,parameter :: kmin = 1.d-4 !H_0*1.e-2 # minimum wavenumber to compute angular power spectrum of lensing potential
!Real*8,parameter :: kmax = 1.d5 ! #*H_0*1.e-2 # maximum wavenumber to compute angular power spectrum of lensing potential 
!Real*8,parameter :: zmax = 1.d1 ! upper red-shift fiducial value in Hill and Spergel paper 1312.4525
!Real*8,parameter :: zmin = 5.d-3 ! lower red-shift fiducial value in Hill and Spergel paper 1312.4525 
!Real*8,parameter :: sigma8 = 0.817d0  !#0.8285 # RMS matter fluctuations today in linear theory 
!Real*8,parameter :: Mmax = 5.d15 ! upper virial mass fiducial value in Hill and Spergel paper 1312.4525
Integer*4,parameter :: number_of_sites = 256 ! Number of sites to generate Gaussian field (of the form 2**N)
Integer*4,parameter :: Nyquist = number_of_sites/2 + 1 ! Nyquist index 
Integer*4,parameter :: number_of_realizations = 1d4 ! number of realizations of field
Integer*4,parameter :: number_of_bins = 1d3 ! number of bins for A_k distribution
!Integer*4,parameter :: lmax = 1d4 ! highest multipole
!Integer*4,parameter :: lmin = 1d1 ! lowest multipole 
Real*8,parameter :: Volume = dble(number_of_sites)*L ! Volume in real space (1 dimension) 
!Real*8 :: Normalization ! Normalization constant for matter power spectrum (Equation (A7) in Eisenstein and Hu)


End Module fiducial
