Program ngpdf
!########
! Modules 
!########

    Use omp_lib
    Use MKL_DFTI
    Use fiducial 
    Use functions
    Use arrays 
    Implicit none

!#############
! Declarations 
!#############

Integer*4 :: i,j,index,counter
type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
Double Complex,dimension(number_of_sites) :: GaussianField2,MeanGaussianField,VarianceGaussianField,NGField2,CovMatrixAlpha
Double Complex,dimension(number_of_sites) :: MeanAlphak
Integer :: isign,Status
Real*8 :: mean_Ak

!############
! Assignments
!############


!#############################
! Allocating memory for arrays 
!#############################

allocate(k(number_of_sites),PowerSpectrum(number_of_sites),fkG(number_of_sites),GaussianField(number_of_realizations,number_of_sites),&
NGField(number_of_realizations,number_of_sites),fkNG(number_of_sites),Ak(number_of_realizations),&
All_alphak(number_of_realizations,number_of_sites),All_Ak(number_of_realizations,number_of_sites),&
Ak_edges(1:number_of_bins),stat=status1) 

!#########################################
! Setting seed for random number generator
!#########################################

call random_seed

!############################
! Filling wavevector array in
!############################

If (number_of_sites .ge. 2) then 

    Do i=1, Nyquist-1

        k(i)        =  dble(i-1)*fundamental_wavenumber

        k(number_of_sites-i+1) = -dble(i)*fundamental_wavenumber

    End do

Else 

    k(number_of_sites) = 0.d0

End If

k(Nyquist) = abs(k(Nyquist))

!###########################################
! Filling (Gaussian) Power Spectrum array in
!###########################################

Do i=1,number_of_sites

    PowerSpectrum(i) = Pk(k(i))

End do

!###############################################################################################
! Main loop starts here; it is splitted in the available processors at the current node (16 or 12)
! which means that every processor computes a realization.
! Main calculations are term(n,m); n labels the trial, and m the term in Matsubara eq 58, in the
! order that they appear in his paper.Note that the p^(3) and p^(4) terms have Dirac delta 
! functions in them, and this collapses many of the sums; this has been done mathematically; the
! code reflects the smaller dimension of the summations.
!###############################################################################################

open(15,file='./output/one-point-distribution-functions.dat') ! Open file to test Gaussian Field in real space 

open(16,file='./output/Ak-distribution-simulations.dat')    ! File for Ak distribution from simulations

!$omp Parallel 
!$omp Do Private(fKG,GaussianField2,My_Desc1_Handle,My_Desc2_Handle,isign)
Do j=1,number_of_realizations

    call random_seed ! Seed for random number generator. It should be called independently by each processor

    ! Set FT of gaussian field. Treat zero wavenumber and Nyquist mode separately, and ensure the field 
    ! is real by assigning complex conjugates to k < 0 modes:
 
    fkG(1) = dcmplx(0.d0,0.d0)*sqrt(Volume)    ! k = 0 component  

    Do i=2,Nyquist-1

        fkG(i)        = sqrt(PowerSpectrum(i)/2.d0)*dcmplx(rndgauss(),rndgauss())*sqrt(Volume)

        fkG(number_of_sites-i+2) = conjg(fkG(i))

    End Do

    fkG(Nyquist) = sqrt(PowerSpectrum(Nyquist))*dcmplx(rndgauss(),0.d0)*sqrt(Volume)     ! Nyquist component is real; treat separately

    GaussianField2  = dcmplx(0.d0,0.d0)    ! Set auxiliary Gaussian field array to zero 

    ! Perform IFFT using intel MKL library:

    isign = +1
    Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, number_of_sites )
    Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiCommitDescriptor( My_Desc1_Handle )
    Status = DftiComputeBackward( My_Desc1_Handle, fkG, GaussianField2)
    Status = DftiFreeDescriptor(My_Desc1_Handle)

    Do i=1,number_of_sites

        GaussianField(j,i) = GaussianField2(i)

    End Do

End Do

!$omp end do
!$omp end parallel

Do i=1,number_of_sites    ! It computes ensemble average of initial Gaussian fields

    MeanGaussianField(i) = sum(GaussianField(:,i))/dble(number_of_realizations) 

End Do

Do j=1,number_of_realizations    ! Redefine Gaussian field to be zero mean field (ensemble average)

    Do i=1,number_of_sites

        GaussianField(j,i) = GaussianField(j,i) - MeanGaussianField(i)

    End Do

End Do    ! By now the ensemble average of the Gaussian field is zero

Do i=1,number_of_sites    ! It computes ensemble variance of (ensemble) zero mean Gaussian fields

    VarianceGaussianField(i) = sum( ( GaussianField(:,i) - MeanGaussianField(i) )**2 )/dble(number_of_realizations) 

End Do

Do j=1,number_of_realizations    ! Redefine Gaussian field to be unit variance (ensemble)

    Do i=1,number_of_sites

        GaussianField(j,i) = GaussianField(j,i)/sqrt(VarianceGaussianField(i))

    End Do

End Do    ! By now the ensemble variance of the Gaussian field should approximately one

Do j=1,number_of_realizations    ! Define non-Gaussian field 

    Do i=1,number_of_sites

        NGField(j,i) = GaussianField(j,i) + fNL*(GaussianField(j,i)**2 - VarianceGaussianField(i))

    End Do

End Do    

!$omp Parallel 
!$omp Do Private(fkNG,NGField2,My_Desc1_Handle,My_Desc2_Handle,isign)
Do j=1,number_of_realizations

    Do i=1,number_of_sites

        NGField2(i) = NGField(j,i)

    End Do

! FFT using intel MKL libraries:

    isign = -1
    Status = DftiCreateDescriptor( My_Desc2_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 1, number_of_sites )
    Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiCommitDescriptor( My_Desc2_Handle )
    Status = DftiComputeForward( My_Desc2_Handle, NGField2, fkNG )
    Status = DftiFreeDescriptor(My_Desc2_Handle)
    

    !fkNG = fkNG/dble(number_of_sites)    ! This is needed so the coefficients agree with inputs when fNL=0.

    fkNG = fKNG/sqrt(Volume)    ! Volume-normalized Fourier coefficients (Equation (12) in Matsubara's paper)

    Do i=1,number_of_sites

        All_alphak(j,i) = fkNG(i)/sqrt(PowerSpectrum(i))    ! Define normalized variables (Equation (40) in Matsubara's paper)

        All_Ak(j,i) = abs(All_alphak(j,i))    ! Define Fourier amplitude (Equation (56) in Matsubara's paper)

        !thetha(i) =  ! Define Fourier phase

    End Do

!####################
! Main loop ends here
!#################### 

End do

!$omp end do
!$omp end parallel

Do i=1,number_of_sites    ! It computes ensemble average of normalized Fourier coefficients of non-Gaussian field

    MeanAlphak(i) = sum(All_alphak(:,i))/dble(number_of_realizations) 

End Do

print *,All_alphak(10,2)
print *,All_alphak(10,number_of_sites)
print *,conjg(All_alphak(10,2))
print *,abs(All_alphak(10,2))**2
print *,All_alphak(10,2)*All_alphak(10,number_of_sites)
print *,MeanAlphak(2),MeanAlphak(number_of_sites)
stop


Do j=1,number_of_realizations    ! Redefine Gaussian field to be zero mean field (ensemble average)

    Do i=1,number_of_sites

        All_alphak(j,i) = All_alphak(j,i) - MeanAlphak(i)

    End Do

End Do    ! By now the ensemble average of the Gaussian field is zero

!CovMatrixAlpha(50) = dcmplx(0.d0,0.d0)

!Do j=1,number_of_realizations    ! It computes ensemble average of initial Gaussian fields

!    CovMatrixAlpha(50) = All_alphak(j,50)*All_alphak(j,number_of_sites) +  CovMatrixAlpha(50)

!End Do

!print *, CovMatrixAlpha(50)/dble(number_of_realizations)
!stop
!Do i=2,Nyquist-1

!    print *, All_alphak(1,i),All_alphak(1,number_of_sites - i + 2)
!    print *, All_alphak(1,i)*All_alphak(1,number_of_sites - i )!+ 2)
!    CovMatrixAlpha(50) = All_alphak(1,i)*All_alphak(1,number_of_sites - i ) + CovMatrixAlpha(50)
!    stop
!End Do

!print *, CovMatrixAlpha(50)
stop
Do i=1,number_of_sites

    If (((maxval(All_Ak(:,i)) .lt. 1.d1) .and. (maxval(All_Ak(:,i)) .gt. 2.d0) ).and. (minval(All_Ak(:,i)) .gt. 1.d-3)) then

        index = i

        print *,'Found!'

        exit

    End If

End Do

!print *,All_Ak(1,index),All_Ak(1,number_of_sites - index + 2)

!stop

Do i=1,number_of_realizations

    Ak(i) = 10**(log10(minval(All_Ak(:,index)) ) + real(i-1)*(log10(maxval(All_Ak(:,index))) -&
    log10(minval(All_Ak(:,index))))/real(number_of_realizations-1))   

End Do

Ak_edges(1) = minval(All_Ak(:,index))

Do i=2,number_of_bins 

    Ak_edges(i) = minval(All_Ak(:,index)) + dble(i)*(maxval(All_Ak(:,index)) - minval(All_Ak(:,index)))/dble(number_of_bins)

End Do

Do i=1,number_of_bins - 1

    mean_Ak = 0.d0

    counter = 0

    Do j=1,number_of_realizations

        If (i .eq. number_of_bins-1) then

            If ( (All_Ak(j,index) .ge. Ak_edges(i) ) .and. (All_Ak(j,index) .le. Ak_edges(i+1) ) ) then

                mean_Ak = mean_Ak + All_Ak(j,index)

                counter = counter + 1


!                If ( (All_Ak(j,number_of_sites - index + 2) .ge. Ak_edges(i) ) .and. (All_Ak(j,number_of_sites -&
!                index + 2) .le. Ak_edges(i+1) ) ) then

!                    mean_Ak = mean_Ak + All_Ak(j,number_of_sites - index + 2)

!                    counter = counter + 1

!                End If

            End If

        Else

            If ( (All_Ak(j,index) .ge. Ak_edges(i) ) .and. (All_Ak(j,index) .lt. Ak_edges(i+1) ) ) then

                mean_Ak = mean_Ak + All_Ak(j,index)

                counter = counter + 1


!                If ( (All_Ak(j,number_of_sites - index + 2) .ge. Ak_edges(i) ) .and. (All_Ak(j,number_of_sites - &
!                index + 2) .lt. Ak_edges(i+1) ) ) then

!                    mean_Ak = mean_Ak + All_Ak(j,number_of_sites - index +2)

!                    counter = counter + 1

!                End If

            End If

        End If
    
    End Do

    write(16,'(2es18.10)') mean_Ak/dble(counter), dble(counter)/dble(number_of_realizations)/(Ak_edges(i+1) - Ak_edges(i))/Pi/2.d0

End Do

Do i=1,number_of_realizations 
   
    write(15,'(3es18.10)') Ak(i),P1G(Ak(i)),P1NG(Ak(i),index)

!    write(16,'(es18.10)') All_Ak(i,index)

!    write(16,'(es18.10)') All_Ak(i,number_of_sites - index + 2)

End Do

close(15)

close(16)

!print *, 'Main loop has ended '
!stop

!If (npts .ge. 2) then 

!    do n=1,min(10,ntrials)
!        write(6,600)n,maxL(n),error(n),term(n,7)/V,term(n,8)/V
!        600 format(i4,es17.10,es17.10,es17.10,es17.10)
!    end do

!End If

!coeff1 = sum(c1(1:ntrials))                           ! What do these terms mean ? 
!coeff2 = sum(c2(1:ntrials)) -sum(c1(1:ntrials)**2)/2d0 ! 

!write(6,6900)coeff1,coeff2
!6900 format('c1,c2 :  ',es17.10,es17.10)

!If (npts .ge. 2) then 

!    fNLhat      = -coeff1/(2d0*coeff2)
!    errorhat    =  1d0/sqrt(2d0*abs(coeff2))

!    write(6,6000)fNLhat,errorhat
!    6000 format('Posterior: fNL     ',es17.10,' +/-',es17.10)

!    write(6,61)sum(maxL)/dble(ntrials),sqrt(1d0/sum(1d0/error**2))

!    61 format('Mean of estimators:',es17.10,' +/-',es17.10)
!Else

!    fNLhat = 0.d0
!    write(6,6023)fNLhat
!    6023 format('Posterior: fNL     ',es17.10)

!End If

!write(6,63)fNL,ntrials,npts
!63 format('True value : ',es17.10,/,' Trials: ',i5,' Number of points :',i6)

!write(6,29)sum(fNLest)/dble(ntrials)
!write(6,30)sum(fNLsqTest)/dble(ntrials)

!If (npts .ge. 2) then 
!    write(2,6000)fNLhat,errorhat
!    write(2,61)sum(maxL)/dble(ntrials),sqrt(1d0/sum(1d0/error**2))
!End If

!write(2,28)fNL
!28 format('True fNL:',es17.10)
!write(2,29)sum(fNLest)/dble(ntrials)
!29 format('Crude estimate from B:',es17.10)

!write(2,30)sum(fNLsqTest)/dble(ntrials)
!30 format('Crude estimate of fNL^2 from T:',es17.10)

!If (npts .ge. 2) then 
!    do n=1,ntrials
!        write(2,20)maxL(n),error(n),fNLest(n)
!    20  format(es17.10,es17.10,es17.10)
!    enddo
!End If

!close(2)

!####################
! Deallocating memory 
!####################

deallocate(k,PowerSpectrum,fkG,GaussianField,NGField,fkNG,Ak)

end Program ngpdf
