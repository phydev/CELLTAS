program main

  use run_cells_m
  use global_m

  implicit none

  ! arg(1): sim_id
  ! arg(2): iseed
  ! arg(3): porosity
  ! arg(4,5): cell position y,z
  do iarg = 1, 6
    call get_command_argument(iarg,length=arglen)
    allocate(character(arglen) :: args(iarg)%arg)
    call get_command_argument(iarg,args(iarg)%arg)
  end do

  read (args(1)%arg,*) sim_id
  read (args(2)%arg,*) iseed
  read (args(3)%arg,*) porosity
  read (args(4)%arg,*) rinit(1)
  read (args(5)%arg,*) rinit(2)
  read (args(6)%arg,*) rinit(3)

  call run_cells(sim_id, iseed,  porosity, rinit)



end program main
