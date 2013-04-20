program tester1
  real(kind(0.0d0)) :: tmp(10)
  integer :: i
  complex(kind(0.0d0)) :: res(10)
  integer*8 p

  i=10
  tmp = (/ ( dble(i), i=1,10 ) /)

  call fft_r2c(p, i, tmp, res)
  
  write(*,*) 'tester1', res, size(tmp) , shape(tmp)

  res = (0.0d0, 1.0d0)

  write(*,*) res

  call dfftw_destroy_plan(p)

end program tester1
