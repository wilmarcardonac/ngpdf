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
Double Complex :: VarianceGaussianField,MeanGaussianField
Integer :: isign,Status
Real*8 :: mean_Ak

!############
! Assignments
!############


!#############################
! Allocating memory for arrays 
!#############################

allocate(k(number_of_sites),PowerSpectrum(number_of_sites),fkG(number_of_sites),GaussianField(number_of_sites),&
NGField(number_of_sites),fkNG(number_of_sites),Ak(number_of_realizations),All_alphak(number_of_realizations,number_of_sites),&
All_Ak(number_of_realizations,number_of_sites),Ak_edges(1:number_of_bins),stat=status1) !, , , calphak(npts), Aksq(npts)
!allocate(Trisp0(npts), term(ntrials,8), c1(ntrials), c2(ntrials), Power(npts), PowerRaw(npts), Bisp(iNyq,iNyq))
!allocate(maxL(ntrials), error(ntrials), fNLest(ntrials), fNLsqTest(ntrials))

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

!###################################################################
! Precompute the bispectrum, for speed (unnormalised at this stage):
!###################################################################

!If (npts .ge. 2) then  
!    Do i1=2,iNyq
!        Do i2=2,iNyq
!            Bisp(i1,i2) = B(k(i1),k(i2),Amp)
!            !if(npts<17) print *,i1,i2,Bisp(i1,i2)
!        End do
!    End do
!End If
!print *, 'Bispectrum ', Bisp
!stop

!##########################
! Open file to write output
!##########################

!open(unit=2,file=fileout)
!write(2,21)NN,npts,fNL,ntrials
!21 format(i3,i10,f8.5,i6)

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
!!$omp do private( calphak,Power, Amp, Trisp0, Ntriangles, NT, temp, ak, k1, rescale_T, AveragePower)

!$omp Parallel 
!$omp Do Private(fKG,GaussianField,MeanGaussianField,VarianceGaussianField,NGField,My_Desc1_Handle,My_Desc2_Handle,isign,fkNG)
Do j=1,number_of_realizations

!    rescale_T = dble(npts)
    
    call random_seed ! Seed for random number generator. It should be called independently by each processor

    ! Set FT of gaussian field. Treat zero wavenumber and Nyquist mode separately, and ensure the field 
    ! is real by assigning complex conjugates to k < 0 modes:
 
    fkG(1) = dcmplx(0.d0,0.d0)    ! k = 0 component  

    Do i=2,Nyquist-1
        fkG(i)        = sqrt(PowerSpectrum(i)/2.d0)*dcmplx(rndgauss(),rndgauss())
        fkG(number_of_sites-i+2) = conjg(fkG(i))
    End Do

    fkG(Nyquist) = sqrt(PowerSpectrum(Nyquist))*dcmplx(rndgauss(),0.d0)     ! Nyquist component is real; treat separately

!    Do i=1,number_of_sites
!        write(15,'(2es18.10)') real(fKG(i)),imag(fKG(i))
!    End Do

    GaussianField   = 0.d0    ! Set Gaussian field array to zero 

    ! Perform IFFT using intel MKL library:

    isign = +1
    Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, number_of_sites )
    Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiCommitDescriptor( My_Desc1_Handle )
    Status = DftiComputeBackward( My_Desc1_Handle, fkG, GaussianField )
    Status = DftiFreeDescriptor(My_Desc1_Handle)

!    Do i=1,number_of_sites
!        write(15,'(2es18.10)') real(GaussianField(i)),imag(GaussianField(i))
!    End Do

    MeanGaussianField = sum(GaussianField)/dble(number_of_sites)    ! Compute mean of Gaussian Field

    GaussianField = GaussianField - MeanGaussianField    ! Make the GaussianField zero mean

!    Do i=1,number_of_sites
!        write(15,'(2es18.10)') real(GaussianField(i)),imag(GaussianField(i))
!    End Do

    MeanGaussianField = sum(GaussianField)/dble(number_of_sites)    ! Compute mean of transformed Gaussian field

    VarianceGaussianField = sum((GaussianField - MeanGaussianField)**2)/dble(number_of_sites)    ! Compute variance of Gaussian field 

!        Amp             = AmpRaw/Variance   ! Rescale amplitude of the power spectrum with variance of Gaussian field  
!        Power           = PowerRaw/Variance                  ! Rescale power spectrum with variance of Gaussian field 
    GaussianField = GaussianField/sqrt(VarianceGaussianField)       ! Rescale Gaussian field to have unit variance 
!        newVariance = real(sum(GaussianField**2)/dble(npts)) ! Variance of rescaled Gaussian field (should be unit)

    MeanGaussianField = sum(GaussianField)/dble(number_of_sites)    ! Mean of Gaussian Field (zero by now)

    VarianceGaussianField = sum((GaussianField - MeanGaussianField)**2)/dble(number_of_sites)    ! Compute variance of Gaussian field (one by now) 

    NGField = GaussianField + fNL*(GaussianField**2 - VarianceGaussianField)    ! Apply local NG

!    Do i=1,number_of_sites
!        write(15,'(2es18.10)') real(NGField(i)),imag(NGField(i))
!    End Do

! FFT using intel MKL libraries:

    isign = -1
    Status = DftiCreateDescriptor( My_Desc2_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 1, number_of_sites )
    Status = DftiSetValue( My_Desc2_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiCommitDescriptor( My_Desc2_Handle )
    Status = DftiComputeForward( My_Desc2_Handle, NGField, fkNG )
    Status = DftiFreeDescriptor(My_Desc2_Handle)

    fkNG = fkNG/dble(number_of_sites)    ! This is needed so the coefficients agree with inputs when fNL=0.

    fkNG = fKNG/sqrt(Volume)    ! Volume-normalized Fourier coefficients (Equation (12) in Matsubara's paper)

!    Do i=1,number_of_sites
!        write(15,'(2es18.10)') real(fkNG(i)),imag(fkNG(i))
!    End Do

    Do i=1,number_of_sites

        All_alphak(j,i) = fkNG(i)/sqrt(PowerSpectrum(i))    ! Define normalized variables (Equation (40) in Matsubara's paper)

        All_Ak(j,i) = abs(sqrt( real(All_alphak(j,i))**2 + imag(All_alphak(j,i))**2 ))    ! Define Fourier amplitude (Equation (56) in Matsubara's paper)

        !thetha(i) =  ! Define Fourier phase

    End Do

! Start with some crude estimates (equal weight) for fNL:

!    fNLest(n)   = 0d0
!    Ntriangles  = 0
!    Var         = 0d0

!    If (npts .ge. 4) then
!        do i1 = 2, iNyq-1
!            do i2 = 2, iNyq-i1+1

!                P1 = Pk(k(i1),Amp)
!                P2 = Pk(k(i2),Amp)
!                P3 = Pk(k(i1+i2-1),Amp)

!                xx = dble(alphak(i1)*alphak(i2)*calphak(i1 + i2 - 1))*sqrt(P1*P2*P3)/(2d0*(P1*P2+P1*P3+P2*P3))

!                fNLest(n)   = fNLest(n) + xx
!                Var         = Var + xx**2
!                Ntriangles  = Ntriangles + 1
!            end do
!        end do
    
!        fNLest(n) = fNLest(n)/dble(Ntriangles)

!        write(6,50)fNLest(n),sqrt((Var/dble(Ntriangles)-fNLest(n)**2)/dble(Ntriangles))
!        50  format('Simple fNL   estimate from B:',es17.10,' +/- ',es17.10)
    
!    Else 

!        write(6,53)fNLest(n)
!        53  format('Bispectrum is zero for 1-,2-point function. Simple fNL estimate from B:',es17.10)

!    End If 
!    stop
! Trispectrum test:

! Simple one first T(k,k,-k,-k):

!    fNLsqTest(n)    = 0d0
!    NT              = 0
!    Var             = 0d0

!    If (npts .ge. 16) then    

!        Do i1 = 2, iNyq/2-1
!            ak = alphak(i1)
!            k1 = k(i1)
!        print *,i1,ak,dble(abs(ak)**4),Trisp(k1,k1,-k1,-k1,Amp),dble(abs(ak)**4)/Trisp(k1,k1,-k1,-k1,Amp)
!            xx              = (dble(abs(ak)**4)-2d0*AveragePower**2)/(Trisp(k1,k1,-k1,-k1,Amp)*dble(npts))
!            fNLsqTest(n)    = fNLsqTest(n) + xx
!            Var             = Var + xx**2
!            NT              = NT  + 1
!        End do

!       Var             = Var/dble(NT)-(fNLsqTest(n)/dble(NT))**2

!        write(6,51)fNLsqTest(n)/dble(NT),sqrt(Var/dble(NT))
!        51  format('Simple fNL^2 estimate from T:',es17.10,' +/- ',es17.10)
!        write(6,52)sqrt(fNLsqTest(n)/dble(NT))
!        52  format('Simple fNL   estimate from T:',es17.10)

!    Else 

!        write(6,54) sqrt(fNLsqTest(n))
!        54  format('Trispectrum T(k,k,-k,-k) is zero for 1-,..,15-point function. Simple fNL estimate from T:',es17.10)
    
!    End If
!stop
!    fNLsqTest(n)    = 0d0
!    NT              = 0

! These modes have no disconnected part [ -> -> -> <- ]:

!    If (npts .ge. 8) then 

!        do i1 = 2, iNyq-1
!            do i2 = 2, iNyq-i1
!                do i3 = 2, iNyq - i1 - i2 + 2

!                    temp = Trisp(k(i1),k(i2),k(i3),-k(i1+i2+i3-2),Amp)
!                    fNLsqTest(n) = (dble(alphak(i1)*alphak(i2)*alphak(i3)*calphak(i1 + i2 + i3 - 2)))/temp
!                    NT = NT + 1

!                    if(abs(k(i1)+k(i2)+k(i3)-k(i1+i2+i3-2))>keps) print *,'Problem'
!                end do
!            end do
!        end do

!        fNLsqTest(n) = fNLsqTest(n)*dble(npts**2)/dble(NT)

!        print *,'Crude fNL^2 estimate from T:',fNLsqTest(n)

!    Else 

!        print *,'Trispectrum -> -> -> <- is zero for # points < 8. Crude fNL^2 estimate from T:',fNLsqTest(n)

!    End If
!stop 
! This is really where the Matsubara calculation starts.
!=======================================================

!    If (npts .ge. 4) then 

!        do i1 = 2, iNyq-1
!            do i2 = 2, iNyq-i1+1
!                term(n, 1) = term(n, 1) + dble(alphak(i1)*alphak(i2)*calphak(i1 + i2 - 1))*Bisp(i1,i2)
!            end do
!        end do
    
        ! Correct for renormalisation of Amp:

!        term(n, 1) = term(n, 1)/sqrt(Variance)

!    End If

! Terms 2,6 (slowest, with 3,7):

!    If ( npts .ge. 8) then 

!        do i1 = 2, iNyq-1
!            do i2 = 2, iNyq-i1
!                do i3=2, iNyq - i1 - i2 + 2
!                    temp = dble(alphak(i1)*alphak(i2)*alphak(i3)*calphak(i1 + i2 + i3 - 2))
!                    term(n, 2) = term(n, 2) + temp * Trisp(k(i1), k(i2), k(i3), -(k(i1) + k(i2) + k(i3)), Amp)
!                    term(n, 6) = term(n, 6) + temp * Bisp(i1,i2) * Bisp(i3,i1+i2+i3-2)
!                end do
!            end do
!        end do

        ! Correct Bisp^2 terms for renormalisation of Amp

!        term(n, 2) =   term(n, 2) / 3d0
!        term(n, 6) = - term(n, 6) / Variance

!    End If

! Terms 3,7 (slowest, with 2,6):

!    If (npts .ge. 4) then 

!        do i1 = 2, iNyq-1
!            do i2 = 2, iNyq-i1+1
!                do i3 = 2, i1 + i2 - 2
!                    temp = dble(alphak(i1)*alphak(i2)*calphak(i3)*calphak(i1 + i2 - i3))
!                    term(n, 3) = term(n, 3) + temp * Trisp(k(i1), k(i2),-k(i3), -(k(i1) + k(i2) - k(i3)),Amp)
!                    term(n, 7) = term(n, 7) + temp * (Bisp(i1,i2)*Bisp(i3,i1+i2-i3)+4d0*Bisp(i1,i3)*Bisp(i2,i1+i2-i3))
!                end do
!            end do
!        end do

!        term(n, 3) =   term(n, 3) / 4d0
!        term(n, 7) = - term(n, 7) / (4d0 * Variance)

!    End If

! Terms 4,5,8:

!    If ( npts .ge. 2) then 

!        do i1 = 2, iNyq
!            do i2 = 2, iNyq
!                term(n, 4) = term(n, 4) + (Aksq(i1) + Aksq(i2) - 1d0)* &
!                Trisp(k(i1), k(i2), -k(i1), -k(i2), Amp)
!            end do
!        end do

!        term(n, 4) = -term(n, 4) / 2d0
!        term(n, 5) =  (term(n, 1)**2) / 2d0

!    End If

!    If (npts .ge. 4) then 

!        do i1 = 2, iNyq-1
!            do i2 = 2, iNyq-i1+1
!                term(n, 8) = term(n, 8) + (Aksq(i1) + Aksq(i2) + Aksq(i1 + i2 - 1) - 1d0) * Bisp(i1,i2)**2
!            end do
!        end do

!        term(n, 8) = term(n, 8) / (2d0 * Variance)

!    End If

! Rescale trispectrum terms:

!    term(n,2) = term(n,2) * rescale_T
!    term(n,3) = term(n,3) * rescale_T
!    term(n,4) = term(n,4) * rescale_T

!    do i=1,8
!        if(i==1) then
!            print *,i,term(n,1)/sqrt(V)
!        else
!            print *,i,term(n,i)/V
!        endif
!    end do

!    coeff1 = term(n, 1)/sqrt(V)
!    coeff2 = sum(term(n, 2:8))/V

!    c1(n) = coeff1
!    c2(n) = coeff2
 
!    If (npts .ge. 2) then 
!        maxL(n)     = -coeff1/(2d0*coeff2)
!        error(n)    = 1d0/sqrt(2d0*abs(coeff2))
!    End If 

!    print *,'Amp:',Amp
!    print *,'Coefficients:', coeff1, coeff2
!    print *,'Perturbations at true fNL:',coeff1*fNL,coeff2*fNL**2

!    tid = omp_get_thread_num()
    
!    If (npts .ge. 2) then
!        write(6,60)n,maxL(n),error(n),tid
!        60  format('Trial ',i5,' MaxL fNL:  ',es17.10,' +/-',es17.10,' Thread ',i4)
!        print *,' '
!    End If

!####################
! Main loop ends here
!#################### 

End do

!$omp end do
!$omp end parallel

Do i=2,Nyquist-1

    print *, All_alphak(1,i),All_alphak(1,number_of_sites - i + 2)
    stop
End Do
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
