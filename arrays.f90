Module arrays
!    Real*8                              :: AveragePower,keps, sigmasq
!    Real*8, dimension(:,:), allocatable :: Bisp
    Double Complex, dimension(:), allocatable :: fkG,GaussianField, NGField, fkNG!,alphak!, calphak
    Real*8,dimension(:),allocatable :: k,PowerSpectrum,Ak,Ak_edges!, Trisp0, c1, c2, Power, PowerRaw, maxL, error, fNLest, fNLsqTest
    Real*8, dimension(:,:), allocatable :: All_Ak
    Double Complex, dimension(:,:), allocatable :: All_alphak 
    Integer :: status1,status2,status3
End Module arrays
