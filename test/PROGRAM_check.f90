program check
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use test_cgrid, only : collect_cgrid
  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    new_testsuite("suite1", collect_cgrid) &
    ]
  WRITE(error_unit, *) ""
  WRITE(error_unit, *) ""
  WRITE(error_unit, *) "-----------------------------------------------------------------"
  WRITE(error_unit, *) "<<<<<<<<<<<<<<<<<<<<<<<<< STARTING TEST >>>>>>>>>>>>>>>>>>>>>>>>>"
  WRITE(error_unit, *) "-----------------------------------------------------------------"
  WRITE(error_unit, *) ""
  WRITE(error_unit, *) ""


  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program check