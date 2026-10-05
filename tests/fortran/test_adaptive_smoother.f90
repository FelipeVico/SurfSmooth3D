program test_adaptive_smoother
  use Mod_Adaptive_Smooth_Surface
  implicit none
  type(adaptive_basis) :: b
  type(adaptive_leaf), allocatable :: leaves(:)
  type(Geometry) :: g
  type(adaptive_leaf) :: checkleaf
  real*8, allocatable :: rv(:,:),u(:,:),v(:,:),w(:)
  real*8 :: a(2,3), original(12,15), error, jac
  integer*8 :: p,n,i,j,ids(1),ier
  character(1024) :: scaffold,cad,root
  do p=1,6
    call setup_adaptive_basis(p,b)
    n=b%np
    allocate(rv(2,n),u(n,n),v(n,n),w(n),leaves(2))
    call vioreanu_simplex_quad(p,n,rv,u,v,w)
    do i=1,2
      leaves(i)%parent=i
      leaves(i)%map=0
      leaves(i)%map(1,2)=1
      leaves(i)%map(2,3)=1
      leaves(i)%points(1,:)=[0d0,1d0,0d0,.5d0,.5d0,0d0]+10*i
      leaves(i)%points(2,:)=[0d0,0d0,1d0,0d0,.5d0,.5d0]
      leaves(i)%points(3,:)=0
      leaves(i)%normals=0
      leaves(i)%normals(3,:)=1
      leaves(i)%diameter=sqrt(2d0)
      allocate(leaves(i)%values(12,n),leaves(i)%height(n))
      leaves(i)%height=0
      leaves(i)%values=0
      leaves(i)%values(1,:)=10*i+rv(1,:)
      leaves(i)%values(2,:)=rv(2,:)
      leaves(i)%values(4,:)=1
      leaves(i)%values(8,:)=1
      leaves(i)%values(12,:)=1
      call adaptive_tail(leaves(i),b)
      if (maxval(leaves(i)%indicator(1:2))>1d-12) stop 1
    enddo
    if (p>=2) then
      ! Inject one normalized highest-degree position mode with known RMS.
      leaves(1)%values(3,:)=.01d0*v(:,n)
      call adaptive_tail(leaves(1),b)
      if (abs(leaves(1)%indicator(1)-.01d0)>1d-12) stop 2
      leaves(1)%values(3,:)=0
      ! Translation and isotropic scaling leave the dimensionless tails fixed.
      leaves(1)%values(1:3,:)=7*leaves(1)%values(1:3,:)+123
      leaves(1)%diameter=7*leaves(1)%diameter
      leaves(1)%values(3,:)=leaves(1)%values(3,:)+.07d0*v(:,n)
      call adaptive_tail(leaves(1),b)
      if (abs(leaves(1)%indicator(1)-.01d0)>1d-11) stop 8
      leaves(1)%values(1:3,:)=leaves(2)%values(1:3,:)
      leaves(1)%values(1,:)=leaves(1)%values(1,:)-10
      leaves(1)%diameter=sqrt(2d0)
      leaves(1)%values(10,:)=.1d0*v(:,n)
      call adaptive_tail(leaves(1),b)
      if (leaves(1)%indicator(2)<.01d0) stop 9
      leaves(1)%values(10,:)=0
    endif
    if (p==4) then
      original=leaves(2)%values
      g%norder_smooth=p
      g%n_order_sf=n
      allocate(g%Points(3,1),g%Tri(6,1),g%Normal_Vert(3,1))
      ids=[1_8]
      call pack_adaptive_parents(g,leaves,ids)
      call refine_geometry_smart(g)
      call funcion_Base_Points(g)
      allocate(g%S_smooth(3,4*n),g%N_smooth(3,4*n),g%du_smooth(3,4*n),g%dv_smooth(3,4*n))
      g%S_smooth=g%Base_Points
      g%N_smooth=g%Base_Points_N
      g%du_smooth=g%Base_Points_U
      g%dv_smooth=g%Base_Points_V
      call merge_adaptive_children(leaves,g,ids)
      if (size(leaves)/=5) stop 3
      if (any(leaves(5)%values/=original)) stop 4
      do j=1,4
        call adaptive_child_map(j,a)
        if (any(leaves(j)%map/=a)) stop 5
        error=maxval(abs(leaves(j)%values(1:2,:)- &
            spread([10d0,0d0]+a(:,1),2,n)-matmul(a(:,2:3),rv)))
        if (error>1d-12) stop 6
        jac=norm2(cross_product_adaptive(leaves(j)%values(4:6,1),leaves(j)%values(7:9,1)))
        if (abs(jac-.25d0)>1d-12) stop 7
      enddo
    endif
    deallocate(rv,u,v,w,leaves)
  enddo
  ! An aliased field can have zero tails and a nonzero off-node residual.
  checkleaf%diameter=2
  checkleaf%indicator=0
  call update_adaptive_checks(checkleaf,1d0,.02d0,[0d0,0d0,-2d0],[0d0,0d0,1d0],[0d0,0d0,1d0])
  if (abs(checkleaf%indicator(3)-.005d0)>1d-15 .or. checkleaf%indicator(4)/=0) stop 10
  checkleaf%indicator=0
  call update_adaptive_checks(checkleaf,1d0,0d0,[0d0,0d0,-2d0],[0d0,0d0,1d0],[1d0,0d0,0d0])
  if (checkleaf%indicator(4)<1d90) stop 11
  checkleaf%indicator=0
  call update_adaptive_checks(checkleaf,1d0,0d0,[0d0,0d0,0d0],[0d0,0d0,1d0],[0d0,0d0,1d0])
  if (minval(checkleaf%indicator(3:4))<1d90) stop 12
  print *, 'PASS: planar and modal tails, selective packing, central child, unchanged leaves, Jacobian scaling.'
  print *, 'PASS: scaled/translated tails, normal tails, independent residual, orientation and gradient checks.'
  if (command_argument_count()>=2) then
    call get_command_argument(1,scaffold)
    call get_command_argument(2,cad)
    call get_command_argument(3,root)
    if (len_trim(root)==0) root = 'adaptive_bounds_result'
    call run_adaptive_smoother(trim(scaffold),trim(cad),trim(root), &
        3_8,8_8,4_8,1_8,5d0,1d-3,2_8,200000_8,ier)
    if (ier/=0) stop 13
    print *, 'PASS: full adaptive solve with bounds-checked adaptive module.'
  endif
end program
