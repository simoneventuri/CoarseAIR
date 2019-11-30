!-----------------------------------------------------------------------------!
!> This module provides an interface for C functions
Module C_interface

  implicit none

  interface

    subroutine fit(n, x, y, a, b, c, err) bind(C,name="fit")
      use iso_C_binding
      integer(C_int), intent(in) :: n
      real(C_double), dimension(n), intent(in) :: x
      real(C_double), dimension(n), intent(in) :: y
      real(C_double), intent(out) :: a
      real(C_double), intent(out) :: b
      real(C_double), intent(out) :: c
      real(C_double), intent(out) :: err
    end subroutine fit

  end interface

end module C_interface
!---------------------------------------------------------------------------!
