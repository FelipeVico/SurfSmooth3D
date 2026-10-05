program test_newton_radius_guard
  use Mod_Smooth_Surface
  use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf
  implicit none
  type(Geometry) :: geom
  type(Feval_stuff) :: feval
  type(projection_radius_guard) :: guard, disabled
  double precision :: base(3,6), normals(3,6), h(6), F(6), dF(6), grad(3,6), targets(3,6)
  double precision :: center(3), scale
  integer *8 :: flags(6), ier, i
  integer :: variant
  logical :: exists
  logical :: failure_mask(6)
  character(len=1000) :: root

  call get_command_argument(1,root)
  if (len_trim(root).eq.0) root = 'test_projection_guard'
  allocate(geom%Points(3,6),geom%skeleton_Points(3,1))
  allocate(geom%Base_Points(3,6),geom%Base_Points_N(3,6))
  geom%npoints = 6
  geom%n_Sk_points = 1
  geom%n_Sf_points = 6
  flags = 0
  normals = 0.0d0
  normals(3,:) = 2.0d0

  do variant=1,2
    center = 0.0d0
    scale = 1.0d0
    if (variant.eq.2) then
      center = [30.0d0,-20.0d0,10.0d0]
      scale = 4.0d0
    endif
    geom%Points = spread(center,2,6)
    do i=1,3
      geom%Points(i,2*i-1) = center(i)-scale
      geom%Points(i,2*i) = center(i)+scale
    enddo
    geom%skeleton_Points(:,1) = center
    call initialize_projection_guard(guard,geom,trim(root)//'_boundary',ier)
    call require(ier.eq.0,'valid enclosing sphere')
    call require(norm2(guard%center-center).lt.1.0d-14,'translation invariant center')
    call require(abs(guard%radius-scale).lt.1.0d-14,'scale invariant radius')
    call begin_projection_guard(guard,1_8,6_8)
    base = spread(center,2,6)
    h = 0.0d0
    h(1) = 5.0d0*scale
    call check_projection_targets(guard,6_8,base,normals,h,flags,1_8,ier)
    call require(ier.eq.0,'point exactly on 10R is accepted')
    h(1) = h(1)*(1.0d0+1.0d-10)
    h(2) = h(1)
    call check_projection_targets(guard,6_8,base,normals,h,flags,1_8,ier)
    call require(ier.eq.7.and.guard%report_written,'all escaped points produce a report')
    call require(.not.any(guard%step_available),'failure before first step has unavailable metrics')
    call initialize_projection_guard(guard,geom,trim(root)//'_boundary',ier)
    inquire(file=trim(root)//'_boundary_newton_failure.txt',exist=exists)
    call require(.not.exists,'initialization removes stale failure report')
  enddo

  call initialize_projection_guard(guard,geom,trim(root)//'_snapshot',ier)
  guard%refinement = 2
  call begin_projection_guard(guard,2_8,6_8)
  base(1,2) = center(1)+scale/2.0d0
  base(2,3) = center(2)+scale/2.0d0
  base(3,4) = center(3)-scale/2.0d0
  base(1,5) = center(1)-scale/2.0d0
  base(2,6) = center(2)-scale/2.0d0
  h = [6.0d0,1.0d0,0.001d0,0.0d0,0.0d0,0.0d0]*scale
  do i=1,3
    call record_projection_step(guard,i,h(i),normals(:,i),3_8)
  enddo
  call record_projection_step(guard,4_8,0.0001d0,normals(:,4),2_8)
  call record_projection_step(guard,6_8,ieee_value(0.0d0,ieee_positive_inf),normals(:,6),3_8)
  flags(4) = 1
  call check_projection_targets(guard,6_8,base,normals,h,flags,3_8,ier)
  call require(ier.eq.7.and.guard%report_written,'full snapshot with one escape')
  call require(abs(guard%step_length(2)-2.0d0*scale).lt.1.0d-14,'physical step includes normal norm')
  call require(guard%step_iteration(4).eq.2,'converged point retains its own latest iteration')
  call require(.not.guard%step_available(5),'never-updated point stays unavailable')
  call require(.not.ieee_is_finite(guard%step_length(6)),'nonfinite step remains explicit')
  guard%filename = trim(root)//'_newton_newton_failure.txt'
  h(1) = scale
  failure_mask = .false.
  failure_mask(1:2) = .true.
  call report_newton_failure(guard,6_8,base,normals,h,flags,3_8, &
      'Newton correction exceeded the existing 1.1 stopping threshold',failure_mask)
  call require(guard%failed.and.guard%report_written,'ordinary Newton failure also writes a full snapshot')
  call require(all(abs(h)*2.0d0.lt.10.0d0*scale),'Newton stop does not require a radius escape')
  call report_newton_failure(disabled,6_8,base,normals,h,flags,3_8,'Disabled',failure_mask)
  call require(.not.disabled%failed,'one-stage Newton reporting remains disabled')
  h(1) = 6.0d0*scale
  guard%filename = trim(root)//'_pre_fmm_newton_failure.txt'

  geom%Base_Points = base
  geom%Base_Points_N = normals
  call begin_projection_guard(guard,1_8,6_8)
  call project_points_eval_roots(6_8,h,base,normals,geom,F,dF,feval,1_8,flags,grad,targets,guard,1_8,ier)
  call require(ier.eq.7,'vertex evaluator returns before touching uninitialized FMM setup')
  call begin_projection_guard(guard,2_8,6_8)
  call fun_roots_derivative(h,geom,F,dF,feval,1_8,flags,grad,targets,guard,1_8,ier)
  call require(ier.eq.7,'RV evaluator returns before touching uninitialized FMM setup')

  call initialize_projection_guard(guard,geom,trim(root)//'_nonfinite',ier)
  call begin_projection_guard(guard,1_8,6_8)
  h(2) = ieee_value(0.0d0,ieee_quiet_nan)
  h(3) = ieee_value(0.0d0,ieee_positive_inf)
  call check_projection_targets(guard,6_8,base,normals,h,flags,1_8,ier)
  call require(ier.eq.7.and.guard%report_written,'NaN and infinity reject without FMM')
  guard%filename = trim(root)//'/missing_directory/newton_failure.txt'
  call check_projection_targets(guard,6_8,base,normals,h,flags,1_8,ier)
  call require(ier.eq.7.and..not.guard%report_written,'diagnostic I/O failure still stops the solve')
  call check_projection_targets(disabled,6_8,base,normals,h,flags,1_8,ier)
  call require(ier.eq.0,'disabled one-stage guard is a no-op')
  call begin_projection_guard(guard,2_8,3_8)
  call require(size(guard%step_length).eq.3.and..not.any(guard%step_available),'new stage resets snapshot')
  print *, 'PASS: radius/Newton failures, physical steps, snapshots, resets and pre-FMM rejection'

contains
  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ',message
      stop 1
    endif
  end subroutine require
end program test_newton_radius_guard
