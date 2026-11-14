subroutine fcfe(rho, sig, cs, cfe) bind(C)
	use CFElocalmod
	implicit none
	real (kind=8), intent(in):: rho, sig, cs ! in cgs
	real (kind=8), intent(out):: cfe(5)
	integer i
	cfe = f_cfelocal(rho*1.e3, sig/100.0, cs/100.0)
!	do i=1,5
!		write(*,*) cfe(i)
!	end do
	return
end subroutine fcfe
