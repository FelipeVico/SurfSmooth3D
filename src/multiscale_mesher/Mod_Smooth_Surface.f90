Module Mod_Smooth_Surface
  
  use Mod_Feval
  use ModType_Smooth_Surface
  use Mod_Newton_Recovery
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  
  implicit none

  ! Per-solve diagnostic state; the original one-stage path leaves it disabled.
  type projection_radius_guard
    logical :: enabled = .false., failed = .false., report_written = .false.
    logical :: recovery_enabled = .true., recovery_report_started = .false.
    integer*8 :: recovery_deferred = 0, recovery_recovered = 0, recovery_unresolved = 0
    double precision :: center(3) = 0.0d0, radius = 0.0d0
    integer *8 :: stage = 0, refinement = 0, iteration = 0
    character(:), allocatable :: filename
    character(:), allocatable :: recovery_filename
    logical, allocatable :: recovered_mask(:)
    double precision, allocatable :: step_length(:)
    integer *8, allocatable :: step_iteration(:)
    logical, allocatable :: step_available(:)
  end type projection_radius_guard
  
contains

  logical function projection_recovery_enabled(guard) result(enabled)
    type(projection_radius_guard), optional, intent(in) :: guard
    enabled = .false.
    if (.not.present(guard)) return
    enabled = guard%enabled .and. guard%recovery_enabled
  end function

  subroutine defer_projection_target(guard,index,iteration,base,normal,initial,last,candidate,f,grad,event,deferred)
    type(projection_radius_guard), intent(in) :: guard
    integer*8, intent(in) :: index,iteration
    real*8, intent(in) :: base(3),normal(3),initial,last,candidate,f,grad(3)
    type(recovery_event), intent(inout) :: event
    logical, intent(out) :: deferred
    real*8 :: ratio,point(3)
    integer :: reason
    point = base + normal*candidate
    call projection_target_status(guard,point,ratio,reason)
    deferred = reason /= 0
    if (.not.deferred) return
    event=recovery_event()
    event%last_residual=ieee_value(0d0,ieee_quiet_nan)
    event%last_gradient=ieee_value(0d0,ieee_quiet_nan)
    event%base_point=base
    event%normal=normal
    event%initial_point=base+normal*initial
    event%last_point=base+normal*last
    event%index = index
    event%reason = reason
    event%iteration = iteration
    event%initial_height = initial
    event%last_height = last
    event%rejected_height = candidate
    event%rejected_point = point
    if (ieee_is_finite(f)) event%last_residual = f
    if (all(ieee_is_finite(grad))) event%last_gradient = grad
  end subroutine

  subroutine append_recovery_event(events,event,used)
    type(recovery_event), allocatable, intent(inout) :: events(:)
    type(recovery_event), intent(in) :: event
    integer*8, intent(inout) :: used
    type(recovery_event), allocatable :: expanded(:)
    if (.not.allocated(events)) allocate(events(16))
    if (used==size(events,kind=8)) then
      allocate(expanded(2*size(events,kind=8)))
      expanded(1:used)=events(1:used)
      call move_alloc(expanded,events)
    endif
    used=used+1
    events(used)=event
  end subroutine

  subroutine evaluate_projection_subset(n,heights,base,normals,active,g,feval,mode,f,df,grad,r,guard,iteration,ier)
    integer*8, intent(in) :: n,mode,iteration
    real*8, intent(in) :: heights(n),base(3,n),normals(3,n)
    logical, intent(in) :: active(n)
    type(Geometry), intent(in) :: g
    type(Feval_stuff), intent(inout) :: feval
    type(projection_radius_guard), intent(inout) :: guard
    real*8, intent(out) :: f(n),df(n)
    real*8, intent(inout) :: grad(3,n),r(3,n)
    integer*8, intent(out) :: ier
    real*8, allocatable :: targets(:,:),directions(:,:),values(:),derivatives(:,:)
    integer*8 :: i,j,nt,flags(n)
    flags=0
    ! Safe retained positions are checked independently of the active mask.
    call check_projection_targets(guard,n,base,normals,heights,flags,iteration,ier)
    if (ier/=0) return
    nt=count(active,kind=8)
    if (nt==0) return
    allocate(targets(3,nt),directions(3,nt),values(nt),derivatives(3,nt))
    j=0
    do i=1,n
      if (.not.active(i)) cycle
      r(:,i)=base(:,i)+normals(:,i)*heights(i)
      j=j+1
      targets(:,j)=r(:,i)
      directions(:,j)=normals(:,i)
    enddo
    call eval_density_grad_FMM(g,targets,directions,nt,values,derivatives,feval,mode)
    j=0
    do i=1,n
      if (.not.active(i)) cycle
      j=j+1
      f(i)=values(j)
      ! Match My_Newton's existing arithmetic for the unaffected active group.
      df(i)=derivatives(1,j)*normals(1,i)+derivatives(2,j)*normals(2,i)+derivatives(3,j)*normals(3,i)
      grad(:,i)=derivatives(:,j)
      r(:,i)=targets(:,j)
    enddo
  end subroutine

  subroutine finish_projection_recovery(n,base,normals,initial,heights,deferred,events,grad,r,g,feval,mode, &
      flags,guard,ier)
    integer*8, intent(in) :: n,mode,flags(n)
    real*8, intent(in) :: base(3,n),normals(3,n),initial(n)
    real*8, intent(inout) :: heights(n),grad(3,n),r(3,n)
    logical, intent(in) :: deferred(n)
    type(recovery_event), intent(inout) :: events(:)
    type(Geometry), intent(in) :: g
    type(Feval_stuff), intent(inout) :: feval
    type(projection_radius_guard), intent(inout) :: guard
    integer*8, intent(out) :: ier
    logical :: failures(n)
    real*8 :: snapshot(n)
    integer*8 :: i,j
    ier=0
    if (.not.any(deferred)) return
    call recover_levelset_targets(n,base,normals,initial,heights,deferred,grad,r,g,feval,mode, &
        guard%center,guard%radius,events,ier)
    do j=1,size(events,kind=8)
      i=events(j)%index
      guard%recovered_mask(i)=events(j)%recovered
    enddo
    call write_recovery_report(guard,n,deferred,events)
    if (ier==0) return
    failures=.false.
    snapshot=heights
    do j=1,size(events,kind=8)
      i=events(j)%index
      failures(i)=.not.events(j)%recovered
      if (failures(i)) snapshot(i)=events(j)%rejected_height
    enddo
    call report_newton_failure(guard,n,base,normals,snapshot,flags,guard%iteration, &
        'Local downward-crossing recovery failed; see Newton recovery diagnostic',failures)
    ier=7
  end subroutine

  subroutine write_recovery_report(guard,n,deferred,events)
    type(projection_radius_guard), intent(inout) :: guard
    integer*8, intent(in) :: n
    logical, intent(in) :: deferred(n)
    type(recovery_event), intent(in) :: events(:)
    integer :: unit,ios,close_ios
    integer*8 :: i,j,nd,nr
    nd=count(deferred,kind=8)
    nr=0
    do i=1,size(events,kind=8)
      if (events(i)%recovered) nr=nr+1
    enddo
    guard%recovery_deferred=guard%recovery_deferred+nd
    guard%recovery_recovered=guard%recovery_recovered+nr
    guard%recovery_unresolved=guard%recovery_unresolved+nd-nr
    print *, 'Newton recovery stage/refinement:',guard%stage,guard%refinement
    print *, '  deferred/recovered/unresolved:',nd,nr,nd-nr
    if (.not.allocated(guard%recovery_filename)) then
      print *, 'WARNING: no Newton recovery diagnostic filename'
      return
    endif
    if (guard%recovery_report_started) then
      open(newunit=unit,file=guard%recovery_filename,status='old',position='append',iostat=ios)
    else
      open(newunit=unit,file=guard%recovery_filename,status='replace',iostat=ios)
    endif
    if (ios/=0) then
      print *, 'WARNING: unable to write Newton recovery diagnostic: ',guard%recovery_filename
      return
    endif
    if (.not.guard%recovery_report_started) write(unit,'(a)',iostat=ios) 'NEWTON_RECOVERY_V1'
    guard%recovery_report_started=.true.
    if (ios==0) write(unit,*,iostat=ios) 'BEGIN ',guard%stage,guard%refinement,nd
    do i=1,size(events,kind=8)
      if (ios/=0) cycle
      write(unit,*,iostat=ios) 'TARGET ',events(i)%index,events(i)%reason,events(i)%iteration
      if (ios==0) write(unit,'(a)',iostat=ios) 'OUTCOME '//trim(events(i)%outcome)
      if (ios==0) write(unit,*,iostat=ios) 'BASE_POINT ',events(i)%base_point
      if (ios==0) write(unit,*,iostat=ios) 'NORMAL ',events(i)%normal
      if (ios==0) write(unit,*,iostat=ios) 'INITIAL_POINT ',events(i)%initial_point
      if (ios==0) write(unit,*,iostat=ios) 'LAST_POINT ',events(i)%last_point
      if (ios==0) write(unit,*,iostat=ios) 'INITIAL_HEIGHT ',events(i)%initial_height
      if (ios==0) write(unit,*,iostat=ios) 'LAST_HEIGHT ',events(i)%last_height
      if (ios==0) write(unit,*,iostat=ios) 'REJECTED_HEIGHT ',events(i)%rejected_height
      if (ios==0) write(unit,*,iostat=ios) 'REJECTED_POINT ',events(i)%rejected_point
      if (ios==0) write(unit,*,iostat=ios) 'LAST_RESIDUAL ',events(i)%last_residual
      if (ios==0) write(unit,*,iostat=ios) 'LAST_GRADIENT ',events(i)%last_gradient
      if (ios==0) write(unit,*,iostat=ios) 'SIGMA ',events(i)%sigma
      if (ios==0) write(unit,*,iostat=ios) 'INTERVAL ',events(i)%lo,events(i)%hi
      if (ios==0) write(unit,*,iostat=ios) 'ROOT ',events(i)%root
      if (ios==0) write(unit,*,iostat=ios) 'RESIDUAL ',events(i)%residual
      if (ios==0) write(unit,*,iostat=ios) 'SLOPE ',events(i)%slope
      if (ios==0) write(unit,*,iostat=ios) 'POSITION_TOLERANCE ',events(i)%position_tolerance
      if (ios==0) write(unit,*,iostat=ios) 'COUNTS ',events(i)%evaluations,events(i)%brackets,events(i)%iterations
      if (allocated(events(i)%bracket_lower)) then
        do j=1,size(events(i)%bracket_lower,kind=8)
          if (ios==0) write(unit,*,iostat=ios) 'BRACKET ',events(i)%bracket_lower(j),events(i)%bracket_upper(j), &
              events(i)%bracket_lower_residual(j),events(i)%bracket_upper_residual(j)
          if (ios==0) write(unit,*,iostat=ios) 'CANDIDATE ',events(i)%candidate_root(j), &
              events(i)%candidate_residual(j),events(i)%candidate_slope(j),events(i)%candidate_iterations(j)
          if (ios==0) write(unit,'(a)',iostat=ios) 'CANDIDATE_OUTCOME '//trim(events(i)%candidate_outcome(j))
        enddo
      endif
      if (ios==0) write(unit,'(a)',iostat=ios) 'END_TARGET'
    enddo
    if (ios==0) write(unit,*,iostat=ios) 'END ',nr,nd-nr
    close(unit,iostat=close_ios)
    if (ios/=0.or.close_ios/=0) print *, 'WARNING: incomplete Newton recovery diagnostic: ',guard%recovery_filename
  end subroutine

  subroutine validate_recovered_patches(g,feval,mode,guard,ier)
    type(Geometry), intent(in) :: g
    type(Feval_stuff), intent(inout) :: feval
    type(projection_radius_guard), intent(inout) :: guard
    integer*8, intent(in) :: mode
    integer*8, intent(out) :: ier
    integer*8 :: p,np,m,nc,k,i,j,l,first,last,nt,nb,nbatch,owner,patches(g%ntri)
    integer*8, allocatable :: owners(:)
    real*8, allocatable :: rv(:,:),umat(:,:),vmat(:,:),w(:),uv(:,:),pol(:,:),du(:,:),dv(:,:)
    real*8, allocatable :: coefficients(:,:),targets(:,:),normals(:,:),grad(:,:),f(:),crosses(:,:),jac0(:)
    real*8 :: d(2,g%n_order_sf),anchor(3),r_u(3),r_v(3),b_u(3),b_v(3),cr(3),cr0(3),jac,ng
    real*8 :: u,v,t,lu(6),lv(6),ratio,dotn,radial
    integer :: reason
    ier=0
    if (.not.allocated(guard%recovered_mask)) return
    if (.not.any(guard%recovered_mask)) return
    np=g%n_order_sf
    if (size(guard%recovered_mask,kind=8)/=g%ntri*np) then
      call reject_recovered_geometry(guard,'Invalid recovered patch layout',ier)
      return
    endif
    nb=0
    do k=1,g%ntri
      first=(k-1)*np+1
      last=k*np
      if (.not.any(guard%recovered_mask(first:last))) cycle
      nb=nb+1
      patches(nb)=k
      do i=first,last
        if (.not.all(ieee_is_finite(g%du_smooth(:,i))).or. &
            .not.all(ieee_is_finite(g%dv_smooth(:,i))).or..not.ieee_is_finite(g%ds_smooth(i))) then
          call reject_recovered_geometry(guard,'Nonfinite recovered patch derivatives',ier)
          return
        endif
        call crossproduct(g%du_smooth(:,i),g%dv_smooth(:,i),cr)
        call crossproduct(g%Base_Points_U(:,i),g%Base_Points_V(:,i),cr0)
        jac=norm2(cr)
        ng=norm2(g%N_smooth(:,i))
        if (.not.all(ieee_is_finite(cr)).or..not.all(ieee_is_finite(cr0)).or. &
            .not.all(ieee_is_finite(g%N_smooth(:,i))).or..not.ieee_is_finite(jac).or. &
            .not.ieee_is_finite(norm2(cr0)).or..not.ieee_is_finite(ng).or. &
            jac<=0d0.or.norm2(cr0)<=0d0.or.ng<=0d0) then
          call reject_recovered_geometry(guard,'Degenerate recovered patch node',ier)
          return
        endif
        ratio=jac/norm2(cr0)
        dotn=dot_product(cr/jac,g%N_smooth(:,i)/ng)
        if (.not.ieee_is_finite(ratio).or..not.ieee_is_finite(dotn).or.ratio<1d-3.or.dotn<.99d0) then
          call reject_recovered_geometry(guard,'Recovered patch node orientation or Jacobian failed',ier)
          return
        endif
      enddo
    enddo
    p=g%norder_smooth
    m=2*p+4
    nc=(m+1)*(m+2)/2
    allocate(rv(2,np),umat(np,np),vmat(np,np),w(np),uv(2,nc),pol(np,nc),du(np,nc),dv(np,nc))
    allocate(coefficients(3,np))
    call vioreanu_simplex_quad(p,np,rv,umat,vmat,w)
    l=0
    do i=0,m
      do j=0,m-i
        l=l+1
        uv(:,l)=[real(i,8)/m,real(j,8)/m]
        call koorn_ders(uv(:,l),p,np,pol(:,l),d)
        du(:,l)=d(1,:)
        dv(:,l)=d(2,:)
      enddo
    enddo
    ! Oversampled patch checks are grouped to avoid one source-tree build per patch.
    nbatch=max(1_8,32768_8/nc)
    do owner=1,nb,nbatch
      nt=min(nbatch,nb-owner+1)*nc
      allocate(targets(3,nt),normals(3,nt),crosses(3,nt),jac0(nt),f(nt),grad(3,nt),owners(nt))
      l=0
      do k=owner,min(nb,owner+nbatch-1)
        first=(patches(k)-1)*np+1
        last=patches(k)*np
        anchor=g%S_smooth(:,first)
        coefficients=matmul(g%S_smooth(:,first:last)-spread(anchor,2,np),transpose(umat))
        do j=1,nc
          l=l+1
          owners(l)=patches(k)
          targets(:,l)=anchor+matmul(coefficients,pol(:,j))
          r_u=matmul(coefficients,du(:,j))
          r_v=matmul(coefficients,dv(:,j))
          call crossproduct(r_u,r_v,cr)
          crosses(:,l)=cr
          u=uv(1,j)
          v=uv(2,j)
          t=1-u-v
          lu=[1-4*t,4*u-1,0d0,4*(t-u),4*v,-4*v]
          lv=[1-4*t,0d0,4*v-1,-4*u,4*u,4*(t-v)]
          ! Subtract an anchor for translation-stable launch derivatives.
          b_u=matmul(g%Points(:,g%Tri(:,patches(k)))- &
              spread(g%Points(:,g%Tri(1,patches(k))),2,6),lu)
          b_v=matmul(g%Points(:,g%Tri(:,patches(k)))- &
              spread(g%Points(:,g%Tri(1,patches(k))),2,6),lv)
          call crossproduct(b_u,b_v,cr0)
          jac0(l)=norm2(cr0)
          normals(:,l)=t*g%Normal_Vert(:,g%Tri(1,patches(k)))+ &
              u*g%Normal_Vert(:,g%Tri(2,patches(k)))+v*g%Normal_Vert(:,g%Tri(3,patches(k)))
          call projection_target_status(guard,targets(:,l),radial,reason)
          if (reason/=0.or..not.all(ieee_is_finite(cr)).or..not.all(ieee_is_finite(cr0)).or. &
              .not.ieee_is_finite(jac0(l)).or..not.ieee_is_finite(norm2(cr)).or. &
              jac0(l)<=0d0.or.norm2(cr)<=0d0) then
            call reject_recovered_geometry(guard,'Invalid oversampled recovered patch',ier)
            return
          endif
          if (norm2(cr)/jac0(l)<1d-3) then
            call reject_recovered_geometry(guard,'Degenerate oversampled recovered patch Jacobian',ier)
            return
          endif
        enddo
      enddo
      call eval_density_grad_FMM(g,targets,normals,nt,f,grad,feval,mode)
      do i=1,nt
        jac=norm2(crosses(:,i))
        ng=norm2(grad(:,i))
        if (.not.all(ieee_is_finite(grad(:,i))).or..not.ieee_is_finite(ng).or.ng*guard%radius<=1d-14) then
          call reject_recovered_geometry(guard,'Invalid oversampled field gradient',ier)
          return
        endif
        ratio=jac/jac0(i)
        dotn=-dot_product(crosses(:,i)/jac,grad(:,i)/ng)
        if (.not.ieee_is_finite(ratio).or..not.ieee_is_finite(dotn).or.ratio<1d-3.or.dotn<.99d0) then
          call reject_recovered_geometry(guard,'Oversampled recovered patch orientation or Jacobian failed',ier)
          return
        endif
      enddo
      deallocate(targets,normals,crosses,jac0,f,grad,owners)
    enddo
  end subroutine

  subroutine reject_recovered_geometry(guard,message,ier)
    type(projection_radius_guard), intent(inout) :: guard
    character(*), intent(in) :: message
    integer*8, intent(out) :: ier
    integer :: unit,ios
    ier=6
    guard%failed=.true.
    print *, 'Recovered surface validation failed: ',message
    if (.not.allocated(guard%recovery_filename)) return
    open(newunit=unit,file=guard%recovery_filename,status='old',position='append',iostat=ios)
    if (ios/=0) then
      print *, 'WARNING: cannot append recovery validation diagnostic'
      return
    endif
    write(unit,'(a)',iostat=ios) 'VALIDATION_FAILED '//trim(message)
    close(unit)
  end subroutine


  double precision function newton_correction_limit(guard) result(limit)
    implicit none
    type(projection_radius_guard), optional, intent(in) :: guard

    ! Preserve the original one-stage Newton stopping rule. Opt-in two-stage
    ! and adaptive solves instead bound physical targets with their guard.
    limit = 1.1d0
    if (present(guard)) then
      if (guard%enabled) limit = huge(1.0d0)
    endif
  end function newton_correction_limit

  subroutine initialize_projection_guard(guard, Geometry1, output_root, ier, newton_recovery)
    implicit none
    type(projection_radius_guard), intent(out) :: guard
    type(Geometry), intent(in) :: Geometry1
    character(len=*), intent(in) :: output_root
    integer *8, intent(out) :: ier
    logical, optional, intent(in) :: newton_recovery
    double precision :: lower(3), upper(3)
    integer :: unit, ios
    integer *8 :: i

    ier = 0
    guard%enabled = .true.
    if (present(newton_recovery)) guard%recovery_enabled = newton_recovery
    guard%filename = trim(output_root)//'_newton_failure.txt'
    guard%recovery_filename = trim(output_root)//'_newton_recovery.txt'
    open(newunit=unit, file=guard%recovery_filename, status='old', iostat=ios)
    if (ios.eq.0) close(unit, status='delete', iostat=ios)
    open(newunit=unit, file=guard%filename, status='old', iostat=ios)
    if (ios.eq.0) close(unit, status='delete', iostat=ios)

    if (.not.all(ieee_is_finite(Geometry1%Points)) .or. &
        .not.all(ieee_is_finite(Geometry1%skeleton_Points))) then
      ier = 7
      print *, 'Newton radius guard: nonfinite original geometry'
      return
    endif
    lower = min(minval(Geometry1%Points,dim=2), minval(Geometry1%skeleton_Points,dim=2))
    upper = max(maxval(Geometry1%Points,dim=2), maxval(Geometry1%skeleton_Points,dim=2))
    guard%center = lower + (upper-lower)/2.0d0
    do i=1,Geometry1%npoints
      guard%radius = max(guard%radius, norm2(Geometry1%Points(:,i)-guard%center))
    enddo
    do i=1,Geometry1%n_Sk_points
      guard%radius = max(guard%radius, norm2(Geometry1%skeleton_Points(:,i)-guard%center))
    enddo
    if (.not.ieee_is_finite(guard%radius) .or. guard%radius.le.0.0d0) then
      ier = 7
      print *, 'Newton radius guard: invalid original enclosing sphere'
    endif
  end subroutine initialize_projection_guard

  subroutine begin_projection_guard(guard, stage, npoints)
    implicit none
    type(projection_radius_guard), optional, intent(inout) :: guard
    integer *8, intent(in) :: stage, npoints
    if (.not.present(guard)) return
    if (.not.guard%enabled) return
    if (allocated(guard%step_length)) deallocate(guard%step_length,guard%step_iteration,guard%step_available)
    allocate(guard%step_length(npoints),guard%step_iteration(npoints),guard%step_available(npoints))
    if (allocated(guard%recovered_mask)) deallocate(guard%recovered_mask)
    allocate(guard%recovered_mask(npoints))
    guard%recovered_mask = .false.
    guard%stage = stage
    guard%iteration = 0
    guard%failed = .false.
    guard%report_written = .false.
    guard%step_length = ieee_value(0.0d0,ieee_quiet_nan)
    guard%step_iteration = 0
    guard%step_available = .false.
  end subroutine begin_projection_guard

  subroutine record_projection_step(guard, index, step, normal, iteration)
    implicit none
    type(projection_radius_guard), optional, intent(inout) :: guard
    integer *8, intent(in) :: index, iteration
    double precision, intent(in) :: step, normal(3)
    if (.not.present(guard)) return
    if (.not.guard%enabled) return
    guard%step_length(index) = abs(step)*norm2(normal)
    guard%step_iteration(index) = iteration
    guard%step_available(index) = .true.
  end subroutine record_projection_step

  subroutine check_projection_targets(guard, npoints, base, normals, heights, flag_con, iteration, ier)
    implicit none
    type(projection_radius_guard), optional, intent(inout) :: guard
    integer *8, intent(in) :: npoints, flag_con(npoints), iteration
    double precision, intent(in) :: base(3,npoints), normals(3,npoints), heights(npoints)
    integer *8, intent(out) :: ier
    integer *8 :: i, nbad
    double precision :: point(3), radial_ratio
    integer :: reason

    ier = 0
    if (.not.present(guard)) return
    if (.not.guard%enabled) return
    guard%iteration = iteration
    nbad = 0
    do i=1,npoints
      point = base(:,i) + normals(:,i)*heights(i)
      call projection_target_status(guard,point,radial_ratio,reason)
      if (reason.ne.0) nbad = nbad + 1
    enddo
    if (nbad.eq.0) return
    guard%failed = .true.
    ier = 7
    print *, 'Newton radius guard stopped stage/refinement/iteration:', guard%stage,guard%refinement,iteration
    print *, '  rejected points:', nbad, '  enclosing radius:', guard%radius, '  limit: 10*R'
    call write_projection_failure(guard,npoints,base,normals,heights,flag_con,nbad)
  end subroutine check_projection_targets

  subroutine projection_target_status(guard, point, radial_ratio, reason)
    implicit none
    type(projection_radius_guard), intent(in) :: guard
    double precision, intent(in) :: point(3)
    double precision, intent(out) :: radial_ratio
    integer, intent(out) :: reason
    reason = 0
    radial_ratio = ieee_value(0.0d0,ieee_quiet_nan)
    if (.not.all(ieee_is_finite(point))) then
      reason = 2
      return
    endif
    radial_ratio = norm2((point-guard%center)/guard%radius)
    if (.not.ieee_is_finite(radial_ratio) .or. radial_ratio.gt.10.0d0) reason = 1
  end subroutine projection_target_status

  subroutine report_newton_failure(guard,npoints,base,normals,heights,flag_con,iteration,message,failure_mask)
    implicit none
    type(projection_radius_guard), optional, intent(inout) :: guard
    integer *8, intent(in) :: npoints, flag_con(npoints), iteration
    double precision, intent(in) :: base(3,npoints), normals(3,npoints), heights(npoints)
    character(len=*), intent(in) :: message
    logical, intent(in) :: failure_mask(npoints)
    integer *8 :: i, nbad
    integer :: reason
    double precision :: point(3), radial_ratio

    if (.not.present(guard)) return
    if (.not.guard%enabled) return
    guard%iteration = iteration
    guard%failed = .true.
    nbad = 0
    do i=1,npoints
      point = base(:,i) + normals(:,i)*heights(i)
      call projection_target_status(guard,point,radial_ratio,reason)
      if (failure_mask(i) .or. reason.ne.0) nbad = nbad + 1
    enddo
    print *, 'Two-stage Newton stopped: ',message
    print *, '  stage/refinement/iteration:',guard%stage,guard%refinement,iteration,'  failed points:',nbad
    call write_projection_failure(guard,npoints,base,normals,heights,flag_con,nbad,failure_mask,message)
  end subroutine report_newton_failure

  subroutine write_projection_failure(guard,npoints,base,normals,heights,flag_con,nbad,failure_mask,message)
    implicit none
    intrinsic :: merge
    type(projection_radius_guard), intent(inout) :: guard
    integer *8, intent(in) :: npoints, flag_con(npoints), nbad
    double precision, intent(in) :: base(3,npoints), normals(3,npoints), heights(npoints)
    logical, optional, intent(in) :: failure_mask(npoints)
    character(len=*), optional, intent(in) :: message
    integer :: unit, ios, close_ios, reason
    integer *8 :: i
    double precision :: point(3), radial_ratio
    character(len=8) :: stage_name
    character(len=*), parameter :: columns = &
        'columns point_index base_x base_y base_z target_x target_y target_z '// &
        'rejected reason converged step_available step_length step_iteration distance_over_radius'

    guard%report_written = .false.
    stage_name = 'rv_nodes'
    if (guard%stage.eq.1) stage_name = 'vertices'
    open(newunit=unit,file=guard%filename,status='replace',action='write',iostat=ios)
    if (ios.ne.0) then
      print *, 'WARNING: unable to write Newton failure diagnostic: ',guard%filename
      return
    endif
    if (present(message)) then
      write(unit,'(a,/,a,1x,a)',iostat=ios) 'NEWTON_PROJECTION_FAILURE_V2','failure_reason',message
    else
      write(unit,'(a)',iostat=ios) 'NEWTON_RADIUS_FAILURE_V1'
    endif
    if (ios.eq.0) write(unit,'(a,1x,a,/,a,1x,i0,/,a,1x,i0,/,a,3(1x,es26.17e3),/,a,1x,es26.17e3,'// &
        '/,a,1x,i0,/,a,1x,i0,/,a,1x,i0,/,a)',iostat=ios) &
        'stage',stage_name,'refinement',guard%refinement, &
        'iteration',guard%iteration,'center',guard%center,'radius',guard%radius, &
        'limit_multiplier',10,'num_points',npoints,'num_rejected',nbad,columns
    if (ios.eq.0) then
      do i=1,npoints
        point = base(:,i) + normals(:,i)*heights(i)
        call projection_target_status(guard,point,radial_ratio,reason)
        if (present(failure_mask)) then
          if (failure_mask(i) .and. reason.eq.0) reason = 3
        endif
        write(unit,'(i0,6(1x,es26.17e3),4(1x,i0),1x,es26.17e3,1x,i0,1x,es26.17e3)',iostat=ios) &
            i,base(:,i),point,merge(1,0,reason.ne.0),reason,flag_con(i), &
            merge(1,0,guard%step_available(i)),guard%step_length(i),guard%step_iteration(i),radial_ratio
        if (ios.ne.0) exit
      enddo
    endif
    if (ios.eq.0) then
      if (present(message)) then
        write(unit,'(a)',iostat=ios) 'END_NEWTON_PROJECTION_FAILURE'
      else
        write(unit,'(a)',iostat=ios) 'END_NEWTON_RADIUS_FAILURE'
      endif
    endif
    close(unit,iostat=close_ios)
    guard%report_written = ios.eq.0 .and. close_ios.eq.0
    if (.not.guard%report_written) print *, 'WARNING: incomplete Newton failure diagnostic: ',guard%filename
  end subroutine write_projection_failure
  
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_cos(v,All_Normals,coef,n_Normals)
    implicit none

    !List of calling arguments
    integer *8, intent(in) :: n_Normals
    double precision, intent(out) :: coef
    double precision, intent(in) :: All_Normals(3,n_Normals), v(3)

    !List of local variables
    integer *8 count
    double precision n_times, tol, my_dot, norm1, norm2
    n_times=0
    tol=1.0d-10
    norm1=sqrt(v(1)**2+v(2)**2+v(3)**2)
    do count=1,n_Normals
      norm2=sqrt(All_Normals(1,count)**2+All_Normals(2,count)**2+All_Normals(3,count)**2)
      my_dot=dot_product(v,All_Normals(:,count))/(norm1*norm2)
      if (my_dot<0) then
        my_dot=0
      endif
      n_times=n_times+my_dot**2
    enddo
    coef=1/n_times
    return
  end subroutine find_cos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine refine_geometry(Geometry1,n_order_sf)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1
    integer *8 , intent(in) :: n_order_sf

    !List of local variables
    integer *8 count,contador_indices
    double precision, allocatable :: Points(:,:), Normal_Vert(:,:)
    integer *8, allocatable :: Tri(:,:)
    double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
    double precision Pa(3),Pb(3),Pc(3),Pd(3),Pe(3),Pf(3),Pg(3),Ph(3),Pi(3)
    double precision Nor1(3),Nor2(3),Nor3(3),Nor4(3),Nor5(3),Nor6(3)
    double precision Nor_a(3),Nor_b(3),Nor_c(3),Nor_d(3),Nor_e(3),Nor_f(3),Nor_g(3),Nor_h(3),Nor_i(3)
    double precision U(9),V(9)
    double precision F_x(9),F_y(9),F_z(9),dS(9)
    double precision nP_x(9),nP_y(9),nP_z(9)
    double precision U_x(9),U_y(9),U_z(9),V_x(9),V_y(9),V_z(9)
    integer *8 m,N,n_order_aux

    allocate(Points(3,Geometry1%ntri*15))
    allocate(Normal_Vert(3,Geometry1%ntri*15))
    allocate(Tri(6,Geometry1%ntri*4))
    contador_indices=1
    n_order_aux=9
    U=[0.2500d0,0.7500d0,0d0,0.2500d0,0.5000d0,0.7500d0,0.2500d0,0d0,0.2500d0]
    V=[0d0,0d0,0.2500d0,0.2500d0,0.2500d0,0.2500d0,0.5000d0,0.7500d0,0.7500d0]
    do count=1,Geometry1%ntri
      P1=Geometry1%Points(:,Geometry1%Tri(1,count))
      P2=Geometry1%Points(:,Geometry1%Tri(2,count))
      P3=Geometry1%Points(:,Geometry1%Tri(3,count))
      P4=Geometry1%Points(:,Geometry1%Tri(4,count))
      P5=Geometry1%Points(:,Geometry1%Tri(5,count))
      P6=Geometry1%Points(:,Geometry1%Tri(6,count))
      call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_aux)
      Pa=[F_x(1),F_y(1),F_z(1)]
      Pb=[F_x(2),F_y(2),F_z(2)]
      Pc=[F_x(3),F_y(3),F_z(3)]
      Pd=[F_x(4),F_y(4),F_z(4)]
      Pe=[F_x(5),F_y(5),F_z(5)]
      Pf=[F_x(6),F_y(6),F_z(6)]
      Pg=[F_x(7),F_y(7),F_z(7)]
      Ph=[F_x(8),F_y(8),F_z(8)]
      Pi=[F_x(9),F_y(9),F_z(9)]
      Nor1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
      Nor2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
      Nor3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
      Nor4=(Nor1+Nor2)/2.0d0
      Nor5=(Nor2+Nor3)/2.0d0
      Nor6=(Nor3+Nor1)/2.0d0
      Nor_a=(Nor1+Nor4)/2.0d0
      Nor_b=(Nor4+Nor2)/2.0d0
      Nor_c=(Nor1+Nor6)/2.0d0
      Nor_d=(Nor4+Nor6)/2.0d0
      Nor_e=(Nor4+Nor5)/2.0d0
      Nor_f=(Nor5+Nor2)/2.0d0
      Nor_g=(Nor5+Nor6)/2.0d0
      Nor_h=(Nor3+Nor6)/2.0d0
      Nor_i=(Nor3+Nor5)/2.0d0
      Points(:,contador_indices)=P1
      Points(:,contador_indices+1)=Pa
      Points(:,contador_indices+2)=P4
      Points(:,contador_indices+3)=Pb
      Points(:,contador_indices+4)=P2
      Points(:,contador_indices+5)=Pc
      Points(:,contador_indices+6)=Pd
      Points(:,contador_indices+7)=Pe
      Points(:,contador_indices+8)=Pf
      Points(:,contador_indices+9)=P6
      Points(:,contador_indices+10)=Pg
      Points(:,contador_indices+11)=P5
      Points(:,contador_indices+12)=Ph
      Points(:,contador_indices+13)=Pi
      Points(:,contador_indices+14)=P3

      Normal_Vert(:,contador_indices)=Nor1
      Normal_Vert(:,contador_indices+1)=Nor_a
      Normal_Vert(:,contador_indices+2)=Nor4
      Normal_Vert(:,contador_indices+3)=Nor_b
      Normal_Vert(:,contador_indices+4)=Nor2
      Normal_Vert(:,contador_indices+5)=Nor_c
      Normal_Vert(:,contador_indices+6)=Nor_d
      Normal_Vert(:,contador_indices+7)=Nor_e
      Normal_Vert(:,contador_indices+8)=Nor_f
      Normal_Vert(:,contador_indices+9)=Nor6
      Normal_Vert(:,contador_indices+10)=Nor_g
      Normal_Vert(:,contador_indices+11)=Nor5
      Normal_Vert(:,contador_indices+12)=Nor_h
      Normal_Vert(:,contador_indices+13)=Nor_i
      Normal_Vert(:,contador_indices+14)=Nor3
      Tri(:,(count-1)*4+1)=contador_indices*[1, 1, 1, 1, 1, 1]+[0, 2, 9, 1, 6, 5]
      Tri(:,(count-1)*4+2)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 11, 9, 7, 10, 6]
      Tri(:,(count-1)*4+3)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 4, 11, 3, 8, 7]
      Tri(:,(count-1)*4+4)=contador_indices*[1, 1, 1, 1, 1, 1]+[9, 11, 14, 10, 13, 12]
      contador_indices=contador_indices+15
    enddo
    Geometry1%npoints=Geometry1%ntri*15
    Geometry1%ntri=Geometry1%ntri*4
    m=Geometry1%npoints
    N=Geometry1%ntri
    Geometry1%n_Sf_points=N*n_order_sf

    if (allocated(Geometry1%Points)) then
      deallocate(Geometry1%Points)
    endif
    allocate(Geometry1%Points(3,m))
    if (allocated(Geometry1%Tri)) then
      deallocate(Geometry1%Tri)
    endif
    allocate(Geometry1%Tri(6,N))
    if (allocated(Geometry1%Normal_Vert)) then
      deallocate(Geometry1%Normal_Vert)
    endif
    allocate(Geometry1%Normal_Vert(3,m))

    Geometry1%Points=Points
    Geometry1%Tri=Tri
    Geometry1%Normal_Vert=Normal_Vert

    deallocate(Points)
    deallocate(Normal_Vert)
    deallocate(Tri)

    return
  end subroutine refine_geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine funcion_Base_Points(Geometry1)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1

    !List of local variables
    double precision U(Geometry1%n_order_sf),V(Geometry1%n_order_sf),w(Geometry1%n_order_sf)
    double precision UV(2,Geometry1%n_order_sf)
    double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),N1(3),N2(3),N3(3)
    double precision F_x(Geometry1%n_order_sf),F_y(Geometry1%n_order_sf)
    double precision F_z(Geometry1%n_order_sf),dS(Geometry1%n_order_sf)
    double precision nP_x(Geometry1%n_order_sf),nP_y(Geometry1%n_order_sf),nP_z(Geometry1%n_order_sf)
    double precision U_x(Geometry1%n_order_sf),U_y(Geometry1%n_order_sf),U_z(Geometry1%n_order_sf)
    double precision V_x(Geometry1%n_order_sf),V_y(Geometry1%n_order_sf),V_z(Geometry1%n_order_sf)

    integer *8 :: count,n_order_sf, itype, npols, norder_smooth,npol_smooth
    double precision, allocatable :: umatr(:,:),vmatr(:,:)


    
    norder_smooth = Geometry1%norder_smooth
    n_order_sf=Geometry1%n_order_sf

    npol_smooth = (norder_smooth+1)*(norder_smooth+2)/2
    allocate(umatr(npol_smooth,npol_smooth))
    allocate(vmatr(npol_smooth,npol_smooth))
    

    !
    ! get the smooth nodes and quadrature weights
    !
     npols = npol_smooth
     call vioreanu_simplex_quad(norder_smooth,npols,UV,umatr,vmatr,w)
     U = UV(1,:)
     V = UV(2,:)
    
    


    
    if (allocated(Geometry1%Base_Points)) then
      deallocate(Geometry1%Base_Points)
    endif
    if (allocated(Geometry1%Base_Points_N)) then
      deallocate(Geometry1%Base_Points_N)
    endif
    if (allocated(Geometry1%Base_Points_U)) then
      deallocate(Geometry1%Base_Points_U)
    endif
    if (allocated(Geometry1%Base_Points_V)) then
      deallocate(Geometry1%Base_Points_V)
    endif
    if (allocated(Geometry1%w_smooth)) then
      deallocate(Geometry1%w_smooth)
    endif
    allocate(Geometry1%Base_Points(3,Geometry1%n_Sf_points))
    allocate(Geometry1%Base_Points_N(3,Geometry1%n_Sf_points))
    allocate(Geometry1%Base_Points_U(3,Geometry1%n_Sf_points))
    allocate(Geometry1%Base_Points_V(3,Geometry1%n_Sf_points))
    allocate(Geometry1%w_smooth(Geometry1%n_Sf_points))
    do count=1,Geometry1%ntri
      P1=Geometry1%Points(:,Geometry1%Tri(1,count))
      P2=Geometry1%Points(:,Geometry1%Tri(2,count))
      P3=Geometry1%Points(:,Geometry1%Tri(3,count))
      P4=Geometry1%Points(:,Geometry1%Tri(4,count))
      P5=Geometry1%Points(:,Geometry1%Tri(5,count))
      P6=Geometry1%Points(:,Geometry1%Tri(6,count))
      N1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
      N2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
      N3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
      call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_sf)
      nP_x=N1(1)*(1.0d0-U-V)+N2(1)*U+N3(1)*V
      nP_y=N1(2)*(1.0d0-U-V)+N2(2)*U+N3(2)*V
      nP_z=N1(3)*(1.0d0-U-V)+N2(3)*U+N3(3)*V
      Geometry1%Base_Points(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_x
      Geometry1%Base_Points(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_y
      Geometry1%Base_Points(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_z
      Geometry1%Base_Points_N(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_x
      Geometry1%Base_Points_N(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_y
      Geometry1%Base_Points_N(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_z
      Geometry1%Base_Points_U(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_x
      Geometry1%Base_Points_U(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_y
      Geometry1%Base_Points_U(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_z

      Geometry1%Base_Points_V(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_x
      Geometry1%Base_Points_V(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_y
      Geometry1%Base_Points_V(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_z

      Geometry1%w_smooth((count-1)*n_order_sf+1:(count)*n_order_sf)=w
    enddo
    return
  end subroutine funcion_Base_Points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  
subroutine find_smooth_surface(Geometry1, Feval_stuff_1, adapt_flag, ier, guard)
  implicit none

  !
  ! this is the main driver routine that finds the smooth surface via
  ! convolution
  !
  double precision :: my_cross(3)

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1
  type ( Feval_stuff ), intent(inout) :: Feval_stuff_1
  integer *8, intent(in) :: adapt_flag

  !List of local variables
  integer *8 :: flag,maxiter,ipoint,count,n_order_sf,itri
  double precision :: h(Geometry1%n_Sf_points),tol
  !double precision :: Norm_N
  double precision :: r_t(3,Geometry1%n_Sf_points)
  double precision :: F(Geometry1%n_Sf_points)
  double precision :: grad_F(3,Geometry1%n_Sf_points)
  double precision :: dVdu(3),dVdv(3),dhdu,dhdv,h_var
  integer *8 ier
  type(projection_radius_guard), optional, intent(inout) :: guard

  ier = 0

  !
  ! check the memory allocation
  !
  
  if (allocated(Geometry1%du_smooth)) then
    deallocate(Geometry1%du_smooth)
  endif
  allocate(Geometry1%du_smooth(3,Geometry1%n_Sf_points))
        
  if (allocated(Geometry1%dv_smooth)) then
    deallocate(Geometry1%dv_smooth)
  endif
  allocate(Geometry1%dv_smooth(3,Geometry1%n_Sf_points))
        
  if (allocated(Geometry1%ds_smooth)) then
    deallocate(Geometry1%ds_smooth)
  endif
  allocate(Geometry1%ds_smooth(Geometry1%n_Sf_points))

  if (allocated(Geometry1%S_smooth)) then
    deallocate(Geometry1%S_smooth)
  endif
  allocate(Geometry1%S_smooth(3,Geometry1%n_Sf_points))

  if (allocated(Geometry1%N_smooth)) then
    deallocate(Geometry1%N_smooth)
  endif
  allocate(Geometry1%N_smooth(3,Geometry1%n_Sf_points))

  if (allocated(Geometry1%ru_smooth)) then
    deallocate(Geometry1%ru_smooth)
  endif
  allocate(Geometry1%ru_smooth(3,Geometry1%n_Sf_points))

  if (allocated(Geometry1%rv_smooth)) then
    deallocate(Geometry1%rv_smooth)
  endif
  allocate(Geometry1%rv_smooth(3,Geometry1%n_Sf_points))



   if (.not.allocated(Geometry1%height)) then
     allocate (Geometry1%height(Geometry1%n_Sf_points))

     do count=1,Geometry1%n_Sf_points
       Geometry1%height(count)=0.0d0
       h(count)=0.0d0
     enddo
   else

     do count=1,Geometry1%n_Sf_points
       h(count)=Geometry1%height(count)
     enddo
   endif

  !
  ! set some parameters for the newton routine for finding the surface
  !
  tol = 1.0d-9
  maxiter = 14
  flag = 0


  call My_Newton(h, tol, maxiter, Geometry1, flag, Feval_stuff_1, &
      adapt_flag, grad_F, r_t, guard)

  if (flag.eq.7) then
    ier = 7
    return
  endif
  
  if (flag .eq. 1) then
    ier = 3
    return
  end if

  if (flag.eq. 3) then
     ier = 4
     return
  endif

  !
  ! copy over the points on the surface
  !
  do count = 1,Geometry1%n_Sf_points
    Geometry1%S_smooth(:,count) = r_t(:,count)
  enddo

  !
  ! compute various partial derivatives
  !
  n_order_sf = Geometry1%n_order_sf
  
  do count=1,Geometry1%n_Sf_points

    Geometry1%N_smooth(:,count) = -1.0d0*grad_F(:,count) &
        /(sqrt(grad_F(1,count)**2 + grad_F(2,count)**2 &
        + grad_F(3,count)**2))

    itri = (count-1)/n_order_sf+1
    h_var=h(count)

    dVdu = Geometry1%Normal_Vert(:,Geometry1%Tri(2,itri)) &
        - Geometry1%Normal_Vert(:,Geometry1%Tri(1,itri))

    dVdv = Geometry1%Normal_Vert(:,Geometry1%Tri(3,itri)) &
        - Geometry1%Normal_Vert(:,Geometry1%Tri(1,itri))

    dhdu = -dot_product(Geometry1%N_smooth(:,count), &
        dVdu*h_var+Geometry1%Base_Points_U(:,count)) &
        /dot_product(Geometry1%N_smooth(:,count) &
        ,Geometry1%Base_Points_N(:,count))

    dhdv = -dot_product(Geometry1%N_smooth(:,count), &
        dVdv*h_var+Geometry1%Base_Points_V(:,count)) &
        /dot_product(Geometry1%N_smooth(:,count), &
        Geometry1%Base_Points_N(:,count))

    Geometry1%ru_smooth(:,count) = dVdu*h_var &
        + Geometry1%Base_Points_N(:,count)*dhdu &
        + Geometry1%Base_Points_U(:,count)

    Geometry1%rv_smooth(:,count) = dVdv*h_var &
        + Geometry1%Base_Points_N(:,count)*dhdv &
        + Geometry1%Base_Points_V(:,count)

    Geometry1%du_smooth(:,count)=Geometry1%ru_smooth(:,count)
    Geometry1%dv_smooth(:,count)=Geometry1%rv_smooth(:,count)
    call crossproduct(Geometry1%ru_smooth(:,count)&
      &,Geometry1%rv_smooth(:,count),my_cross)
    Geometry1%ds_smooth(count)=norm2(my_cross)

    call crossproduct(Geometry1%ru_smooth(:,count), &
        Geometry1%rv_smooth(:,count), my_cross)

    Geometry1%w_smooth(count) = &
        Geometry1%w_smooth(count)*norm2(my_cross)

    Geometry1%ru_smooth(:,count) = &
        Geometry1%ru_smooth(:,count)/norm2(Geometry1%ru_smooth(:,count))

    call crossproduct(Geometry1%N_smooth(:,count), &
        Geometry1%ru_smooth(:,count), my_cross)

    Geometry1%rv_smooth(:,count) = my_cross

  enddo

  if (present(guard)) then
    if (allocated(guard%recovered_mask)) then
      if (any(guard%recovered_mask)) then
        call validate_recovered_patches(Geometry1,Feval_stuff_1,adapt_flag,guard,ier)
        if (ier/=0) return
      endif
    endif
  endif

  !
  ! copy over the h height function
  !
  do count=1,Geometry1%n_Sf_points
    Geometry1%height(count)=h(count)
  enddo

  return
end subroutine find_smooth_surface






subroutine My_Newton(x,tol,maxiter,Geometry1,flag, &
    Feval_stuff_1 ,adapt_flag, grad_F, r_t, guard)
  implicit none
  intrinsic :: merge

  !List of calling arguments
  type (Geometry), intent(in) :: Geometry1
  integer *8, intent(in) :: maxiter
  double precision, intent(inout) :: x(Geometry1%n_Sf_points)
  double precision, intent(in) :: tol
  integer *8, intent(out) :: flag
  type ( Feval_stuff ), intent(inout) :: Feval_stuff_1      !! data type that
  integer *8, intent(in) :: adapt_flag
  double precision, intent(inout) :: r_t(3,Geometry1%n_Sf_points)
  double precision, intent(inout) :: grad_F(3,Geometry1%n_Sf_points)
  type(projection_radius_guard), optional, intent(inout) :: guard


  !
  ! this routine runs a newton iteration to find the smooth surface
  !
  ! Input:
  !   x - initial guess at the h function
  !   tol - pointwise convergence criterion for Newton
  !   maxiter - max iteration count
  !   Geometry1 - all the geo info
  !   Feval_stuff_1 -
  !   adapt_flag - determines what kind of smoothing kernel to use
  !
  ! Output:
  !   x - the h function from the surface
  !   flag - whether or not Newton succeeded
  !   grad_F - gradient of SOMETHING
  !   r_t - points on the smooth surface
  !   
  
  !List of local variables
  integer *8 count,count2
  double precision F(Geometry1%n_Sf_points)
  double precision :: dF(Geometry1%n_Sf_points), err(Geometry1%n_Sf_points)
  double precision :: correction_limit
  integer *8  flag_con(Geometry1%n_Sf_points)
  integer *8 guard_error
  logical :: recovery, deferred(Geometry1%n_Sf_points)
  real*8 :: candidate
  real*8, allocatable :: initial(:),last_finite_f(:),last_finite_grad(:,:)
  type(recovery_event) :: candidate_event
  integer*8 :: events_used
  type(recovery_event), allocatable :: events(:)

  correction_limit = newton_correction_limit(guard)
  call begin_projection_guard(guard,2_8,Geometry1%n_Sf_points)

  recovery = projection_recovery_enabled(guard)
  deferred = .false.
  events_used=0
  if (recovery) then
    initial=x
    allocate(last_finite_f(Geometry1%n_Sf_points),last_finite_grad(3,Geometry1%n_Sf_points))
    last_finite_f=ieee_value(0d0,ieee_quiet_nan)
    last_finite_grad=ieee_value(0d0,ieee_quiet_nan)
    do count2=1,Geometry1%n_Sf_points
      if (.not.all(ieee_is_finite(Geometry1%Base_Points_N(:,count2))) .or. &
          .not.ieee_is_finite(norm2(Geometry1%Base_Points_N(:,count2))).or. &
          norm2(Geometry1%Base_Points_N(:,count2))<=0d0) then
        flag_con=0
        call report_newton_failure(guard,Geometry1%n_Sf_points,Geometry1%Base_Points, &
            Geometry1%Base_Points_N,x,flag_con,0_8,'Zero or nonfinite pseudonormal', &
            [(count==count2,count=1,Geometry1%n_Sf_points)])
        flag=3
        return
      endif
    enddo
  endif

  ! initialize the errors and convergence flags
  do count2 = 1,Geometry1%n_Sf_points
    err(count2) = tol+1.0d0
    flag_con(count2) = 0
  enddo

  count=1
  flag=0

 4999 format(i10,i9,5x,e10.2)

  ! print out convergence information
  write (*,*) 'iteration,  #targets,   err'  
  write (*,4999) count, &
        Geometry1%n_Sf_points-sum(flag_con)-sum(merge(1_8,0_8,deferred)), 1.0d0

  ! the x function has been initialized to 0
  do while ( (maxval(err)>tol) .and. (count<maxiter) )

    ! call a step of Newton

    !print *, "Entering fun_roots_derivative"

    if (any(deferred)) then
      call evaluate_projection_subset(Geometry1%n_Sf_points,x,Geometry1%Base_Points,Geometry1%Base_Points_N, &
          flag_con==0.and..not.deferred,Geometry1,Feval_stuff_1,adapt_flag,F,dF,grad_F,r_t,guard,count,guard_error)
    else
      call fun_roots_derivative(x, Geometry1, F, dF, Feval_stuff_1, &
          adapt_flag, flag_con, grad_F, r_t, guard, count, guard_error)
    endif
    if (guard_error.ne.0) then
      flag = guard_error
      return
    endif

    count = count+1
    do count2 = 1,Geometry1%n_Sf_points

      if (flag_con(count2) .eq. 0 .and. .not.deferred(count2)) then

        if (recovery) then
          if (ieee_is_finite(F(count2))) last_finite_f(count2)=F(count2)
          if (all(ieee_is_finite(grad_F(:,count2)))) last_finite_grad(:,count2)=grad_F(:,count2)
        endif
        err(count2) = F(count2)/dF(count2)
        call record_projection_step(guard,count2,err(count2),Geometry1%Base_Points_N(:,count2),count-1)
        if (recovery) then
          candidate=x(count2)-err(count2)
          call defer_projection_target(guard,count2,count-1,Geometry1%Base_Points(:,count2), &
              Geometry1%Base_Points_N(:,count2),initial(count2),x(count2),candidate,last_finite_f(count2), &
              last_finite_grad(:,count2),candidate_event,deferred(count2))
          if (deferred(count2)) then
            call append_recovery_event(events,candidate_event,events_used)
            err(count2)=0d0
            cycle
          endif
        endif
        x(count2) = x(count2) - err(count2)

        err(count2) = abs(err(count2))
        
        if ( abs(x(count2)) .lt. 1.0d0) then
          err(count2) = abs(err(count2))
        else
          err(count2) = abs(err(count2)/x(count2))
        end if
        
        
        if (abs(err(count2))<tol) then
          flag_con(count2)=1
        endif
      endif
    enddo

    call check_projection_targets(guard,Geometry1%n_Sf_points,Geometry1%Base_Points, &
        Geometry1%Base_Points_N,x,flag_con,count-1,guard_error)
    if (guard_error.ne.0) then
      flag = guard_error
      return
    endif
    write (*,4999) count, &
        Geometry1%n_Sf_points-sum(flag_con)-sum(merge(1_8,0_8,deferred)), maxval(err)
    if (maxval(err) > correction_limit) then
       call report_newton_failure(guard,Geometry1%n_Sf_points,Geometry1%Base_Points, &
           Geometry1%Base_Points_N,x,flag_con,count-1, &
           'Newton correction exceeded the solver limit',err.gt.correction_limit)
       flag = 3
       if (any(deferred)) call write_recovery_report(guard,Geometry1%n_Sf_points,deferred,events(1:events_used))
       return
    endif
        
  end do

!  stop

  
  if ( maxval(err)>tol ) then
    call report_newton_failure(guard,Geometry1%n_Sf_points,Geometry1%Base_Points, &
        Geometry1%Base_Points_N,x,flag_con,count-1,'Newton iteration limit reached',flag_con.eq.0)
    flag=1
    if (any(deferred)) call write_recovery_report(guard,Geometry1%n_Sf_points,deferred,events(1:events_used))
    return
  endif

  if (any(deferred)) then
    call finish_projection_recovery(Geometry1%n_Sf_points,Geometry1%Base_Points,Geometry1%Base_Points_N, &
        initial,x,deferred,events(1:events_used),grad_F,r_t,Geometry1,Feval_stuff_1,adapt_flag,flag_con,guard,flag)
  endif

  return
end subroutine My_Newton






  
subroutine project_points_eval_roots(nproj, h, base_points, base_normals, &
    Geometry1, F, dF, Feval_stuff_1, adapt_flag, flag_con, grad_F, r_t, guard, iteration, guard_error)
  implicit none

  integer *8, intent(in) :: nproj
  double precision, intent(in) :: h(nproj)
  double precision, intent(in) :: base_points(3,nproj)
  double precision, intent(in) :: base_normals(3,nproj)
  type (Geometry), intent(in) :: Geometry1
  double precision, intent(out) :: F(nproj), dF(nproj)
  type (Feval_stuff), intent(inout) :: Feval_stuff_1
  integer *8, intent(in) :: adapt_flag
  integer *8, intent(in) :: flag_con(nproj)
  double precision, intent(inout) :: grad_F(3,nproj), r_t(3,nproj)
  type(projection_radius_guard), optional, intent(inout) :: guard
  integer *8, optional, intent(in) :: iteration
  integer *8, optional, intent(out) :: guard_error

  integer *8 :: count, ipointer, ntarg
  integer *8 :: current_iteration, local_error
  double precision, allocatable :: r_t2(:,:), v_norm(:,:)
  double precision, allocatable :: grad_F2(:,:), F2(:)

  current_iteration = 0
  if (present(iteration)) current_iteration = iteration
  call check_projection_targets(guard,nproj,base_points,base_normals,h,flag_con,current_iteration,local_error)
  if (present(guard_error)) guard_error = local_error
  if (local_error.ne.0) return

  do count=1,nproj
    r_t(:,count) = base_points(:,count) + base_normals(:,count)*h(count)
  enddo

  ntarg = nproj - sum(flag_con)
  if (ntarg .eq. 0) return

  allocate(r_t2(3,ntarg))
  allocate(v_norm(3,ntarg))
  allocate(grad_F2(3,ntarg))
  allocate(F2(ntarg))

  ipointer = 1
  do count=1,nproj
    if (flag_con(count) .eq. 0) then
      r_t2(:,ipointer) = r_t(:,count)
      v_norm(:,ipointer) = base_normals(:,count)
      ipointer = ipointer + 1
    endif
  enddo

  call eval_density_grad_FMM(Geometry1, r_t2, v_norm, ntarg, F2, &
      grad_F2, Feval_stuff_1, adapt_flag)

  ipointer = 1
  do count=1,nproj
    if (flag_con(count) .eq. 0) then
      dF(count) = dot_product(grad_F2(:,ipointer), base_normals(:,count))
      F(count) = F2(ipointer)
      r_t(:,count) = r_t2(:,ipointer)
      grad_F(:,count) = grad_F2(:,ipointer)
      ipointer = ipointer + 1
    endif
  enddo

  deallocate(r_t2)
  deallocate(v_norm)
  deallocate(grad_F2)
  deallocate(F2)

  return
end subroutine project_points_eval_roots




subroutine project_points_to_levelset(Geometry1, Feval_stuff_1, adapt_flag, &
    nproj, base_points, base_normals, projected_points, heights, grad_F, ier, guard)
  implicit none
  intrinsic :: merge

  type (Geometry), intent(in) :: Geometry1
  type (Feval_stuff), intent(inout) :: Feval_stuff_1
  integer *8, intent(in) :: adapt_flag, nproj
  double precision, intent(in) :: base_points(3,nproj), base_normals(3,nproj)
  double precision, intent(out) :: projected_points(3,nproj), heights(nproj)
  double precision, intent(out) :: grad_F(3,nproj)
  integer *8, intent(out) :: ier
  type(projection_radius_guard), optional, intent(inout) :: guard

  integer *8 :: count, iter, maxiter
  integer *8 :: flag_con(nproj)
  double precision :: F(nproj), dF(nproj), err(nproj)
  double precision :: r_t(3,nproj), maxerr, normn, step, tol, correction_limit
  logical :: failure_mask(nproj), recovery, deferred(nproj)
  real*8 :: candidate
  real*8, allocatable :: initial(:),last_finite_f(:),last_finite_grad(:,:)
  type(recovery_event) :: candidate_event
  integer*8 :: events_used
  type(recovery_event), allocatable :: events(:)

  ier = 0
  correction_limit = newton_correction_limit(guard)
  tol = 1.0d-9
  maxiter = 14
  heights = 0.0d0
  projected_points = base_points
  grad_F = 0.0d0
  flag_con = 0
  call begin_projection_guard(guard,1_8,nproj)
  recovery=projection_recovery_enabled(guard)
  deferred=.false.
  events_used=0
  if (recovery) then
    initial=heights
    allocate(last_finite_f(nproj),last_finite_grad(3,nproj))
    last_finite_f=ieee_value(0d0,ieee_quiet_nan)
    last_finite_grad=ieee_value(0d0,ieee_quiet_nan)
  endif

  do count=1,nproj
    normn = sqrt(sum(base_normals(:,count)**2))
    if (recovery) normn = norm2(base_normals(:,count))
    if (.not.(normn.gt.0.0d0).or..not.ieee_is_finite(normn).or. &
        .not.all(ieee_is_finite(base_normals(:,count)))) then
      failure_mask = .false.
      failure_mask(count) = .true.
      call report_newton_failure(guard,nproj,base_points,base_normals,heights,flag_con,0_8, &
          'Zero or nonfinite vertex pseudonormal',failure_mask)
      ier = 3
      return
    endif
    flag_con(count) = 0
    err(count) = tol + 1.0d0
  enddo

  iter = 1
  write (*,*) 'two-stage vertex Newton iteration,  #targets,   err'
  write (*,'(i10,i9,5x,e10.2)') iter, nproj-sum(flag_con)-sum(merge(1_8,0_8,deferred)), 1.0d0

  do while ((maxval(err).gt.tol) .and. (iter.lt.maxiter))
    if (any(deferred)) then
      call evaluate_projection_subset(nproj,heights,base_points,base_normals,flag_con==0.and..not.deferred, &
          Geometry1,Feval_stuff_1,adapt_flag,F,dF,grad_F,r_t,guard,iter,ier)
    else
      call project_points_eval_roots(nproj, heights, base_points, &
          base_normals, Geometry1, F, dF, Feval_stuff_1, adapt_flag, &
          flag_con, grad_F, r_t, guard, iter, ier)
    endif
    if (ier.ne.0) return

    iter = iter + 1
    maxerr = 0.0d0
    do count=1,nproj
      if (flag_con(count) .eq. 0 .and. .not.deferred(count)) then
        if (recovery) then
          if (ieee_is_finite(F(count))) last_finite_f(count)=F(count)
          if (all(ieee_is_finite(grad_F(:,count)))) last_finite_grad(:,count)=grad_F(:,count)
        endif
        if (abs(dF(count)) .le. 1.0d-300) then
          failure_mask = .false.
          failure_mask(count) = .true.
          call report_newton_failure(guard,nproj,base_points,base_normals,heights,flag_con,iter-1, &
              'Directional level-set derivative is too small for Newton',failure_mask)
          ier = 4
          if (any(deferred)) call write_recovery_report(guard,nproj,deferred,events(1:events_used))
          return
        endif
        step = F(count)/dF(count)
        call record_projection_step(guard,count,step,base_normals(:,count),iter-1)
        if (recovery) then
          candidate=heights(count)-step
          call defer_projection_target(guard,count,iter-1,base_points(:,count),base_normals(:,count), &
              initial(count),heights(count),candidate,last_finite_f(count),last_finite_grad(:,count), &
              candidate_event,deferred(count))
          if (deferred(count)) then
            call append_recovery_event(events,candidate_event,events_used)
            err(count)=0d0
            cycle
          endif
        endif
        heights(count) = heights(count) - step
        err(count) = abs(step)
        if (abs(heights(count)) .ge. 1.0d0) then
          err(count) = abs(step/heights(count))
        endif
        if (err(count) .lt. tol) flag_con(count) = 1
        maxerr = max(maxerr, err(count))
      endif
    enddo

    call check_projection_targets(guard,nproj,base_points,base_normals,heights,flag_con,iter-1,ier)
    if (ier.ne.0) return
    write (*,'(i10,i9,5x,e10.2)') iter, nproj-sum(flag_con)-sum(merge(1_8,0_8,deferred)), &
        maxval(err)
    if (maxerr .gt. correction_limit) then
      call report_newton_failure(guard,nproj,base_points,base_normals,heights,flag_con,iter-1, &
          'Newton correction exceeded the solver limit',err.gt.correction_limit)
      ier = 4
      if (any(deferred)) call write_recovery_report(guard,nproj,deferred,events(1:events_used))
      return
    endif
  enddo

  if (maxval(err) .gt. tol) then
    call report_newton_failure(guard,nproj,base_points,base_normals,heights,flag_con,iter-1, &
        'Newton iteration limit reached',flag_con.eq.0)
    ier = 3
    if (any(deferred)) call write_recovery_report(guard,nproj,deferred,events(1:events_used))
    return
  endif

  if (any(deferred)) then
    call finish_projection_recovery(nproj,base_points,base_normals,initial,heights,deferred,events(1:events_used), &
        grad_F,r_t,Geometry1,Feval_stuff_1,adapt_flag,flag_con,guard,ier)
    if (ier/=0) return
  endif

  do count=1,nproj
    projected_points(:,count) = base_points(:,count) + &
        base_normals(:,count)*heights(count)
  enddo

  return
end subroutine project_points_to_levelset




double precision function triangle_area_from_points(P1, P2, P3)
  implicit none

  double precision, intent(in) :: P1(3), P2(3), P3(3)
  double precision :: u(3), v(3), cr(3)

  u = P2 - P1
  v = P3 - P1
  cr(1) = u(2)*v(3) - u(3)*v(2)
  cr(2) = u(3)*v(1) - u(1)*v(3)
  cr(3) = u(1)*v(2) - u(2)*v(1)
  triangle_area_from_points = 0.5d0*sqrt(sum(cr**2))

  return
end function triangle_area_from_points




subroutine clear_smooth_surface_work_arrays(Geometry1)
  implicit none

  type (Geometry), intent(inout) :: Geometry1

  if (allocated(Geometry1%height)) deallocate(Geometry1%height)
  if (allocated(Geometry1%Base_Points)) deallocate(Geometry1%Base_Points)
  if (allocated(Geometry1%Base_Points_N)) deallocate(Geometry1%Base_Points_N)
  if (allocated(Geometry1%Base_Points_U)) deallocate(Geometry1%Base_Points_U)
  if (allocated(Geometry1%Base_Points_V)) deallocate(Geometry1%Base_Points_V)
  if (allocated(Geometry1%Dummy_targ)) deallocate(Geometry1%Dummy_targ)
  if (allocated(Geometry1%w_smooth)) deallocate(Geometry1%w_smooth)
  if (allocated(Geometry1%du_smooth)) deallocate(Geometry1%du_smooth)
  if (allocated(Geometry1%dv_smooth)) deallocate(Geometry1%dv_smooth)
  if (allocated(Geometry1%ds_smooth)) deallocate(Geometry1%ds_smooth)
  if (allocated(Geometry1%S_smooth)) deallocate(Geometry1%S_smooth)
  if (allocated(Geometry1%N_smooth)) deallocate(Geometry1%N_smooth)
  if (allocated(Geometry1%ru_smooth)) deallocate(Geometry1%ru_smooth)
  if (allocated(Geometry1%rv_smooth)) deallocate(Geometry1%rv_smooth)

  return
end subroutine clear_smooth_surface_work_arrays




subroutine project_scaffold_vertices_to_levelset(Geometry1, Feval_stuff_1, &
    adapt_flag, ier, guard)
  implicit none

  type (Geometry), intent(inout) :: Geometry1
  type (Feval_stuff), intent(inout) :: Feval_stuff_1
  integer *8, intent(in) :: adapt_flag
  integer *8, intent(out) :: ier
  type(projection_radius_guard), optional, intent(inout) :: guard

  logical, allocatable :: is_corner(:),recovered_vertices(:)
  real*8, allocatable :: candidate_points(:,:)
  real*8 :: old_cross(3),new_cross(3)
  integer *8, allocatable :: vertex_ids(:)
  double precision, allocatable :: base_points(:,:), base_normals(:,:)
  double precision, allocatable :: projected_points(:,:), heights(:)
  double precision, allocatable :: grad_F(:,:), area0(:)
  integer *8 :: i, j, itri, nproj, id
  double precision :: disp, mean_disp, max_disp
  double precision :: area1, ratio, min_ratio
  double precision :: p1(3), p2(3), p3(3)

  ier = 0
  allocate(is_corner(Geometry1%npoints))
  is_corner = .false.

  do itri=1,Geometry1%ntri
    do j=1,3
      id = Geometry1%Tri(j,itri)
      if (id.ge.1 .and. id.le.Geometry1%npoints) is_corner(id) = .true.
    enddo
  enddo

  nproj = count(is_corner)
  if (nproj .le. 0) then
    ier = 6
    return
  endif

  allocate(vertex_ids(nproj))
  allocate(base_points(3,nproj))
  allocate(base_normals(3,nproj))
  allocate(projected_points(3,nproj))
  allocate(heights(nproj))
  allocate(grad_F(3,nproj))
  allocate(area0(Geometry1%ntri))

  j = 0
  do i=1,Geometry1%npoints
    if (is_corner(i)) then
      j = j + 1
      vertex_ids(j) = i
      base_points(:,j) = Geometry1%Points(:,i)
      base_normals(:,j) = Geometry1%Normal_Vert(:,i)
    endif
  enddo

  do itri=1,Geometry1%ntri
    p1 = Geometry1%Points(:,Geometry1%Tri(1,itri))
    p2 = Geometry1%Points(:,Geometry1%Tri(2,itri))
    p3 = Geometry1%Points(:,Geometry1%Tri(3,itri))
    area0(itri) = triangle_area_from_points(p1, p2, p3)
    if (area0(itri) .le. 1.0d-300) then
      ier = 6
      return
    endif
  enddo

  call project_points_to_levelset(Geometry1, Feval_stuff_1, adapt_flag, &
      nproj, base_points, base_normals, projected_points, heights, &
      grad_F, ier, guard)
  if (ier .ne. 0) then
    print *, "two-stage scaffold vertex projection failed, ier=", ier
    return
  endif

  if (present(guard)) then
    if (allocated(guard%recovered_mask)) then
      if (any(guard%recovered_mask)) then
        candidate_points=Geometry1%Points
        allocate(recovered_vertices(Geometry1%npoints))
        recovered_vertices=.false.
        do i=1,nproj
          candidate_points(:,vertex_ids(i))=projected_points(:,i)
          recovered_vertices(vertex_ids(i))=guard%recovered_mask(i)
        enddo
        do itri=1,Geometry1%ntri
          p1=candidate_points(:,Geometry1%Tri(1,itri))
          p2=candidate_points(:,Geometry1%Tri(2,itri))
          p3=candidate_points(:,Geometry1%Tri(3,itri))
          ratio=triangle_area_from_points(p1,p2,p3)/area0(itri)
          if (.not.all(ieee_is_finite([p1,p2,p3])).or..not.ieee_is_finite(ratio).or.ratio<1d-10) then
            call reject_recovered_geometry(guard,'Recovered scaffold area or coordinates failed',ier)
            return
          endif
          if (.not.any(recovered_vertices(Geometry1%Tri(1:3,itri)))) cycle
          call crossproduct(p2-p1,p3-p1,new_cross)
          call crossproduct(Geometry1%Points(:,Geometry1%Tri(2,itri))- &
              Geometry1%Points(:,Geometry1%Tri(1,itri)),Geometry1%Points(:,Geometry1%Tri(3,itri))- &
              Geometry1%Points(:,Geometry1%Tri(1,itri)),old_cross)
          if (dot_product(old_cross/norm2(old_cross),new_cross/norm2(new_cross))<=0d0) then
            call reject_recovered_geometry(guard,'Recovered scaffold triangle reversed orientation',ier)
            return
          endif
        enddo
      endif
    endif
  endif

  mean_disp = 0.0d0
  max_disp = 0.0d0
  do i=1,nproj
    disp = sqrt(sum((projected_points(:,i)-base_points(:,i))**2))
    mean_disp = mean_disp + disp
    max_disp = max(max_disp, disp)
    Geometry1%Points(:,vertex_ids(i)) = projected_points(:,i)
  enddo
  mean_disp = mean_disp/dble(nproj)

  do itri=1,Geometry1%ntri
    p1 = Geometry1%Points(:,Geometry1%Tri(1,itri))
    p2 = Geometry1%Points(:,Geometry1%Tri(2,itri))
    p3 = Geometry1%Points(:,Geometry1%Tri(3,itri))
    Geometry1%Points(:,Geometry1%Tri(4,itri)) = (p1 + p2)/2.0d0
    Geometry1%Points(:,Geometry1%Tri(5,itri)) = (p2 + p3)/2.0d0
    Geometry1%Points(:,Geometry1%Tri(6,itri)) = (p3 + p1)/2.0d0
  enddo

  min_ratio = huge(1.0d0)
  do itri=1,Geometry1%ntri
    p1 = Geometry1%Points(:,Geometry1%Tri(1,itri))
    p2 = Geometry1%Points(:,Geometry1%Tri(2,itri))
    p3 = Geometry1%Points(:,Geometry1%Tri(3,itri))
    area1 = triangle_area_from_points(p1, p2, p3)
    ratio = area1/area0(itri)
    min_ratio = min(min_ratio, ratio)
  enddo

  print *, "two-stage scaffold projection summary"
  print *, "  projected vertices:", nproj
  print *, "  mean/max vertex displacement:", mean_disp, max_disp
  print *, "  min projected/original triangle area ratio:", min_ratio

  if (min_ratio .lt. 1.0d-10) then
    print *, "two-stage projected scaffold is degenerate"
    ier = 6
    return
  elseif (min_ratio .lt. 1.0d-6) then
    print *, "WARNING: two-stage projected scaffold has tiny triangle area ratio"
  endif

  call clear_smooth_surface_work_arrays(Geometry1)

  return
end subroutine project_scaffold_vertices_to_levelset




  subroutine check_Gauss(Geometry1,x0,y0,z0,err_rel)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1
    double precision, intent(in) :: x0,y0,z0
    double precision, intent(out) :: err_rel

    !List of local variables
    integer *8 umio,count1,count2,flag,n_order_sf
    double precision  F,Ex,Ey,Ez,R,x,y,z,pi,w,nx,ny,nz

    pi=3.141592653589793238462643383d0
    F=0.0d0

    do count1=1,Geometry1%n_Sf_points
      x=Geometry1%S_smooth(1,count1)
      y=Geometry1%S_smooth(2,count1)
      z=Geometry1%S_smooth(3,count1)
      w=Geometry1%w_smooth(count1)
      nx=Geometry1%N_smooth(1,count1)
      ny=Geometry1%N_smooth(2,count1)
      nz=Geometry1%N_smooth(3,count1)

      R=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      Ex=(x-x0)/(4*pi*R**3)
      Ey=(y-y0)/(4*pi*R**3)
      Ez=(z-z0)/(4*pi*R**3)
      F=F+(Ex*nx+Ey*ny+Ez*nz)*w
    enddo
    err_rel=abs(F-1)
    
    return
  end subroutine check_Gauss





  
  subroutine fun_roots_derivative(h, Geometry1, F, dF, &
      Feval_stuff_1, &
      adapt_flag, flag_con, grad_F, r_t, guard, iteration, guard_error)
    implicit none

    !List of calling arguments
    type (Geometry), intent(in) :: Geometry1
    double precision, intent(in) :: h(Geometry1%n_Sf_points)
    double precision, intent(out) ::  F(Geometry1%n_Sf_points), &
        dF(Geometry1%n_Sf_points)
    type ( Feval_stuff ), intent(inout) :: Feval_stuff_1
    
    integer *8, intent(in) :: adapt_flag
    integer *8, intent(in) :: flag_con(Geometry1%n_Sf_points)
    double precision, intent(inout) :: r_t(3,Geometry1%n_Sf_points), &
        grad_F(3,Geometry1%n_Sf_points)
    type(projection_radius_guard), optional, intent(inout) :: guard
    integer *8, optional, intent(in) :: iteration
    integer *8, optional, intent(out) :: guard_error



    !List of local variables
    !double precision Norm_N
    integer *8 :: count, ipointer, ntarg
    integer *8 :: current_iteration, local_error
    double precision, allocatable :: r_t2(:,:),v_norm(:,:),grad_F2(:,:),F2(:)

    current_iteration = 0
    if (present(iteration)) current_iteration = iteration
    call check_projection_targets(guard,Geometry1%n_Sf_points,Geometry1%Base_Points, &
        Geometry1%Base_Points_N,h,flag_con,current_iteration,local_error)
    if (present(guard_error)) guard_error = local_error
    if (local_error.ne.0) return

    allocate(r_t2(3,Geometry1%n_Sf_points-sum(flag_con)))
    allocate(v_norm(3,Geometry1%n_Sf_points-sum(flag_con)))
    allocate(grad_F2(3,Geometry1%n_Sf_points-sum(flag_con)))
    allocate(F2(Geometry1%n_Sf_points-sum(flag_con)))


    ipointer=1
    do count=1,Geometry1%n_Sf_points

      ! compute the current estimate for the point on the surface
      r_t(:,count) = Geometry1%Base_Points(:,count) &
          + Geometry1%Base_Points_N(:,count)*h(count)

      ! if not converged, add it to the list of targets
      if (flag_con(count) .eq. 0) then
        r_t2(:,ipointer) = r_t(:,count)
        v_norm(:,ipointer) = Geometry1%Base_Points_N(:,count)
        ipointer = ipointer+1
      endif
    enddo

    !
    ! get values and gradients
    !
    ntarg = Geometry1%n_Sf_points-sum(flag_con)

    call eval_density_grad_FMM(Geometry1, r_t2, v_norm, &
        ntarg, F2, grad_F2, Feval_stuff_1, adapt_flag)

    !
    ! update the points
    !
    ipointer=1
    do count=1,Geometry1%n_Sf_points

      if ( flag_con(count) .eq. 0) then

        dF(count)=(grad_F2(1,ipointer)*Geometry1%Base_Points_N(1,count)+grad_F2(2,ipointer)*&
            &Geometry1%Base_Points_N(2,count)+grad_F2(3,ipointer)*Geometry1%Base_Points_N(3,count))
        F(count)=F2(ipointer)
        r_t(:,count) = r_t2(:,ipointer)
        grad_F(:,count) = grad_F2(:,ipointer)
        ipointer=ipointer+1
      endif

    enddo

    deallocate(r_t2)
    deallocate(v_norm)
    deallocate(grad_F2)
    deallocate(F2)

    return
  end subroutine fun_roots_derivative





  
  subroutine refine_geometry_smart(Geometry1)
    implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1

    !List of local variables
    character ( len=100 ) plot_name
    integer *8 count,contador_indices
    double precision, allocatable :: Points(:,:), Normal_Vert(:,:)
    integer *8, allocatable :: Tri(:,:)
    double precision, allocatable :: h_new(:),h_tri(:),h_1(:),h_2(:),h_3(:)
    double precision, allocatable :: h_4(:)

    double precision P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
    double precision Pa(3),Pb(3),Pc(3),Pd(3),Pe(3),Pf(3),Pg(3),Ph(3),Pi(3)
    double precision Nor1(3),Nor2(3),Nor3(3),Nor4(3),Nor5(3),Nor6(3)
    double precision Nor_a(3),Nor_b(3),Nor_c(3),Nor_d(3),Nor_e(3),Nor_f(3),Nor_g(3),Nor_h(3),Nor_i(3)
    double precision U(9),V(9)
    double precision F_x(9),F_y(9),F_z(9),dS(9)
    double precision nP_x(9),nP_y(9),nP_z(9)
    double precision U_x(9),U_y(9),U_z(9),V_x(9),V_y(9),V_z(9)
    double precision, allocatable :: ximat(:,:)
    integer *8 m,N,n_order_aux,n_order_sf,istart
    integer *8 nover,norder, itmp

    n_order_sf=Geometry1%n_order_sf

    allocate(h_tri(n_order_sf))
    allocate(h_1(n_order_sf))
    allocate(h_2(n_order_sf))
    allocate(h_3(n_order_sf))
    allocate(h_4(n_order_sf))

    nover = n_order_sf*4
    allocate(ximat(nover,n_order_sf))

    norder = Geometry1%norder_smooth
    call get_refine_interp_mat(norder,n_order_sf,nover,ximat)

    allocate(Points(3,Geometry1%ntri*15))
    allocate(Normal_Vert(3,Geometry1%ntri*15))
    allocate(Tri(6,Geometry1%ntri*4))
    allocate(h_new(Geometry1%n_Sf_points*4))

    contador_indices=1
    n_order_aux=9
    U=[0.2500d0,0.7500d0,0d0,0.2500d0,0.5000d0,0.7500d0,0.2500d0,0d0,0.2500d0]
    V=[0d0,0d0,0.2500d0,0.2500d0,0.2500d0,0.2500d0,0.5000d0,0.7500d0,0.7500d0]

    do count=1,Geometry1%ntri

      P1=Geometry1%Points(:,Geometry1%Tri(1,count))
      P2=Geometry1%Points(:,Geometry1%Tri(2,count))
      P3=Geometry1%Points(:,Geometry1%Tri(3,count))
      P4=Geometry1%Points(:,Geometry1%Tri(4,count))
      P5=Geometry1%Points(:,Geometry1%Tri(5,count))
      P6=Geometry1%Points(:,Geometry1%Tri(6,count))

      call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_aux)
      Pa=[F_x(1),F_y(1),F_z(1)]
      Pb=[F_x(2),F_y(2),F_z(2)]
      Pc=[F_x(3),F_y(3),F_z(3)]
      Pd=[F_x(4),F_y(4),F_z(4)]
      Pe=[F_x(5),F_y(5),F_z(5)]
      Pf=[F_x(6),F_y(6),F_z(6)]
      Pg=[F_x(7),F_y(7),F_z(7)]
      Ph=[F_x(8),F_y(8),F_z(8)]
      Pi=[F_x(9),F_y(9),F_z(9)]
      Nor1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
      Nor2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
      Nor3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
      Nor4=(Nor1+Nor2)/2.0d0
      Nor5=(Nor2+Nor3)/2.0d0
      Nor6=(Nor3+Nor1)/2.0d0
      Nor_a=(Nor1+Nor4)/2.0d0
      Nor_b=(Nor4+Nor2)/2.0d0
      Nor_c=(Nor1+Nor6)/2.0d0
      Nor_d=(Nor4+Nor6)/2.0d0
      Nor_e=(Nor4+Nor5)/2.0d0
      Nor_f=(Nor5+Nor2)/2.0d0
      Nor_g=(Nor5+Nor6)/2.0d0
      Nor_h=(Nor3+Nor6)/2.0d0
      Nor_i=(Nor3+Nor5)/2.0d0
      Points(:,contador_indices)=P1
      Points(:,contador_indices+1)=Pa
      Points(:,contador_indices+2)=P4
      Points(:,contador_indices+3)=Pb
      Points(:,contador_indices+4)=P2
      Points(:,contador_indices+5)=Pc
      Points(:,contador_indices+6)=Pd
      Points(:,contador_indices+7)=Pe
      Points(:,contador_indices+8)=Pf
      Points(:,contador_indices+9)=P6
      Points(:,contador_indices+10)=Pg
      Points(:,contador_indices+11)=P5
      Points(:,contador_indices+12)=Ph
      Points(:,contador_indices+13)=Pi
      Points(:,contador_indices+14)=P3

      Normal_Vert(:,contador_indices)=Nor1
      Normal_Vert(:,contador_indices+1)=Nor_a
      Normal_Vert(:,contador_indices+2)=Nor4
      Normal_Vert(:,contador_indices+3)=Nor_b
      Normal_Vert(:,contador_indices+4)=Nor2
      Normal_Vert(:,contador_indices+5)=Nor_c
      Normal_Vert(:,contador_indices+6)=Nor_d
      Normal_Vert(:,contador_indices+7)=Nor_e
      Normal_Vert(:,contador_indices+8)=Nor_f
      Normal_Vert(:,contador_indices+9)=Nor6
      Normal_Vert(:,contador_indices+10)=Nor_g
      Normal_Vert(:,contador_indices+11)=Nor5
      Normal_Vert(:,contador_indices+12)=Nor_h
      Normal_Vert(:,contador_indices+13)=Nor_i
      Normal_Vert(:,contador_indices+14)=Nor3
      Tri(:,(count-1)*4+1)=contador_indices*[1, 1, 1, 1, 1, 1]+[0, 2, 9, 1, 6, 5]
      Tri(:,(count-1)*4+2)=contador_indices*[1, 1, 1, 1, 1, 1]+[11, 9, 2, 10, 6, 7]![2, 11, 9, 7, 10, 6]
      Tri(:,(count-1)*4+3)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 4, 11, 3, 8, 7]
      Tri(:,(count-1)*4+4)=contador_indices*[1, 1, 1, 1, 1, 1]+[9, 11, 14, 10, 13, 12]

      h_tri=Geometry1%height((count-1)*n_order_sf+1:(count)*n_order_sf)

      istart = (count-1)*nover+1
      call dmatvec(nover,n_order_sf,ximat,h_tri,h_new(istart))


      contador_indices=contador_indices+15

    enddo

    Geometry1%npoints=Geometry1%ntri*15
    Geometry1%ntri=Geometry1%ntri*4
    m=Geometry1%npoints
    N=Geometry1%ntri
    Geometry1%n_Sf_points=N*n_order_sf

    if (allocated(Geometry1%Points)) then
      deallocate(Geometry1%Points)
    endif
    allocate(Geometry1%Points(3,m))

    if (allocated(Geometry1%Tri)) then
      deallocate(Geometry1%Tri)
    endif
    allocate(Geometry1%Tri(6,N))

    if (allocated(Geometry1%Normal_Vert)) then
      deallocate(Geometry1%Normal_Vert)
    endif
    allocate(Geometry1%Normal_Vert(3,m))

    Geometry1%Points = Points
    Geometry1%Tri = Tri
    Geometry1%Normal_Vert = Normal_Vert

    deallocate(Points)
    deallocate(Normal_Vert)
    deallocate(Tri)


    deallocate(Geometry1%height)

    allocate(Geometry1%height(Geometry1%n_Sf_points))
    Geometry1%height=h_new

    return
  end subroutine refine_geometry_smart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module Mod_Smooth_Surface
