program bie_read_go3
  ! Link ONLY to the existing BIE library, never to SurfSmooth3D.
  implicit none
  character(len=2000) :: filename, output
  integer*8 :: npatches, npoints
  integer*8, allocatable :: orders(:), offsets(:), types(:)
  real*8, allocatable :: values(:,:), coefs(:,:), weights(:)
  call get_command_argument(1,filename)
  call get_command_argument(2,output)
  call open_gov3_geometry_mem(trim(filename),npatches,npoints)
  allocate(orders(npatches),offsets(npatches+1),types(npatches))
  allocate(values(12,npoints),coefs(9,npoints),weights(npoints))
  call open_gov3_geometry(trim(filename),npatches,orders,offsets,types, &
    npoints,values,coefs,weights)
  if (any(weights.le.0)) stop 1
  open(unit=33,file=trim(output),access='stream',form='unformatted',status='replace')
  write(33) npatches,npoints,values,coefs,weights
  close(33)
  print *, 'BIE_GO3_READ_PASSED',npatches,npoints,sum(weights)
end program
