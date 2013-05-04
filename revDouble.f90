subroutine revDouble(dat, n)
  integer i, n, ndiv2, nMiP1
  double precision :: dat(n), tmp

  ndiv2 = n/2

  do i = 1, ndiv2
     nMiP1 = n - i + 1
     tmp = dat(i)
     dat(i) = dat(nMiP1)
     dat(nMiP1) = tmp
  end do

end subroutine revDouble
