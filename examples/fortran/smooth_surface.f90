! File-based native example; forwards to the unchanged smoother entry points.
! smooth_surface scaffold.gidmsh cad.txt output_root [order] [refinements] [two_stage]
program smooth_surface
  implicit none
  character(4096) :: mesh, cad, root, argument
  integer*8 :: order, refinements, two_stage, filetype, ier
  if (command_argument_count()<3) then
    print *, 'Usage: smooth_surface scaffold cad output_root [order] [refinements] [two_stage]'
    stop 1
  endif
  call get_command_argument(1,mesh)
  call get_command_argument(2,cad)
  call get_command_argument(3,root)
  order=8
  refinements=0
  two_stage=1
  call get_command_argument(4,argument)
  if (len_trim(argument)>0) read(argument,*) order
  call get_command_argument(5,argument)
  if (len_trim(argument)>0) read(argument,*) refinements
  call get_command_argument(6,argument)
  if (len_trim(argument)>0) read(argument,*) two_stage
  call get_filetype(trim(mesh),filetype,ier)
  if (ier/=0) stop 2
  if (two_stage==0) then
    call multiscale_mesher_unif_refine(trim(mesh),filetype,1_8,trim(cad),16_8, &
        order,refinements,1_8,5d0,trim(root),ier)
  else
    call multiscale_mesher_unif_refine_two_stage(trim(mesh),filetype,1_8,trim(cad),16_8, &
        order,refinements,1_8,5d0,trim(root),ier)
  endif
  if (ier/=0) then
    print *, 'Smoother error: ',ier
    stop 3
  endif
end program smooth_surface
