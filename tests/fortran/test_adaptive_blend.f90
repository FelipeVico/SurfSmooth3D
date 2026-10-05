program test_adaptive_blend
  use Mod_Adaptive_Blend_Surface
  use iso_c_binding
  implicit none
  character(1024) :: scaffold,cad,root
  integer(c_int64_t) :: id,ier,n,np,nt,j
  real(c_double) :: sphere(4),targets(3,2),before(4,2),after(4,2)
  real(c_double), allocatable :: original(:,:),values(:,:),meta(:,:),uv(:,:),projected(:,:)
  integer(c_int64_t), allocatable :: ids(:)
  call get_command_argument(1,scaffold)
  call get_command_argument(2,cad)
  call get_command_argument(3,root)
  if (len_trim(root)==0) stop 1
  call blend_open(trim(scaffold)//c_null_char,trim(cad)//c_null_char,trim(root)//c_null_char, &
      4_c_int64_t,4_c_int64_t,1_c_int64_t,2d0,500000_c_int64_t,id,ier)
  if (ier/=0) stop 2
  call blend_size(id,n,np,sphere,ier)
  if (ier/=0 .or. n<=0 .or. np/=15) stop 3
  allocate(original(12,n*np),meta(9,n))
  call blend_get(id,n,np,original,meta)
  targets=original(1:3,1:2)
  call blend_sigma(id,2_c_int64_t,targets,before,ier)
  if (ier/=0) stop 4
  allocate(ids(2))
  ids=[1_c_int64_t,n]
  call blend_refine(id,2_c_int64_t,ids,n*np,ier)
  if (ier/=8) stop 5
  call blend_refine(id,2_c_int64_t,ids,500000_c_int64_t,ier)
  if (ier/=0) stop 6
  j=n
  call blend_size(id,n,np,sphere,ier)
  if (n/=j+6) stop 7
  deallocate(meta)
  allocate(values(12,n*np),meta(9,n))
  call blend_get(id,n,np,values,meta)
  if (maxval(abs(values(:,4*np+1:5*np)-original(:,np+1:2*np)))/=0) stop 8
  if (any(abs(meta(3:8,2)-[.5d0,.5d0,-.5d0,0d0,0d0,-.5d0])>1d-14)) stop 9
  call blend_sigma(id,2_c_int64_t,targets,after,ier)
  if (ier/=0 .or. maxval(abs(before-after))/=0) stop 10
  nt=2
  allocate(uv(2,nt),projected(12,nt))
  ids=[1_c_int64_t,2_c_int64_t]
  uv(:,1)=[.2d0,.3d0]
  uv(:,2)=[.1d0,.1d0]
  call blend_project(id,nt,ids,uv,projected,ier)
  if (ier/=0 .or. any(abs(projected)>1d100)) stop 11
  call blend_close(id,ier)
  if (ier/=0) stop 12
  call blend_size(id,n,np,sphere,ier)
  if (ier/=10) stop 13
  call blend_close_all()
  print *, 'PASS: adaptive blend session, selective children, fixed sigma, projection and cleanup.'
end program
