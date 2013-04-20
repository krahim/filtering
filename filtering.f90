
subroutine  filterTimeSeriesDirect(timeSeries, &
     theFilter, nTimeSeries, nFilter, nOutput, res) 
  implicit none
   
  integer :: i, j, nTimeSeries, nFilter, nOutput, nFilterM1
  double precision :: timeSeries(nTimeSeries), theFilter(nFilter), &
       res(nOutput)
  
  nFilterM1 = nFilter -1
  call revDouble(theFilter, nFilter)

  do i = 1, nOutput
     res(i) =  timeSeries(i) * theFilter(1)
     do j = 1, nFilterM1
        res(i) = res(i) + timeSeries(i+j) * theFilter(j +1)
     enddo
  enddo

end subroutine filterTimeSeriesDirect


subroutine filterWfft(timeSeries, padFilter, ndata, nfilter, nfft, nRes, &
     blockLen, res, fftFilter, res1, tsfft) 
!                K-1
!          y_t = Sum g_k x_{t+K-1-k},  t = 0, ..., N-K+1
!                k=0

!      http://en.wikipedia.org/wiki/Cyclic_convolution
    
  implicit none

  integer :: ndata, nfilter, nfft, nRes, blockLen, &
       i, lenLeft
  real(kind(0.0d0)) :: timeSeries(ndata), res(nRes), padFilter(nfft), &
       res1(nfft)
  complex(kind(0.0d0)) :: fftFilter(nfft), tsfft(nfft)
  integer*8 :: p1, p2

  call fft_r2c(p1, nfft, padFilter, fftFilter)
  
  i = 0

  while1: do while ((i*blockLen + nfft) < ndata)
     call fft_r2c(p1, nfft, &
          timeSeries((i*blockLen +1) : (i*blockLen + nfft)), &
          tsfft)

     tsfft = (tsfft * fftFilter)/dble(nfft)
       
     call fft_c2r(p2, nfft, tsfft, res1)
     
     res((i*blockLen +1) : (i*blockLen + blockLen)) = &
          res1(nfilter:nfft)
     
     i = i + 1
     
  end do while1
  
  lenLeft = ndata - i*blockLen  
  
  if (lenLeft /= 0) then 
     res1 = 0.0d0
     res1(1:lenLeft) = timeSeries((i*blockLen +1) : ndata)
    
     call fft_r2c(p1, nfft, &
          res1, &
          tsfft)
     tsfft = (tsfft * fftFilter)/dble(nfft)
     
     call fft_c2r(p2, nfft, tsfft, res1)

     res((i*blockLen +1) : nRes) = &
          res1(nfilter:(nfilter + lenLeft - nfilter))
  end if

  call dfftw_destroy_plan(p1)
  call dfftw_destroy_plan(p2)
  
end subroutine filterWfft
