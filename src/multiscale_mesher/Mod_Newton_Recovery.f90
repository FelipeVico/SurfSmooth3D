! Local, derivative-free recovery of guarded Newton targets.  The CAD sources
! and sigma definition are borrowed read-only from the ordinary field evaluator.
module Mod_Newton_Recovery
  use, intrinsic :: ieee_arithmetic
  use ModType_Smooth_Surface, only: Geometry
  use Mod_Feval, only: Feval_stuff, eval_density_grad_FMM
  use Mod_Fast_Sigma, only: function_eval_sigma
  implicit none
  private

  integer, parameter, public :: recovery_scan_half = 128
  integer, parameter, public :: recovery_max_iterations = 80
  integer, parameter, public :: recovery_batch_limit = 32768
  double precision, parameter, public :: recovery_sigma_extent = 8.0d0
  double precision, parameter, public :: recovery_residual_tolerance = 1.0d-8
  double precision, parameter, public :: recovery_relative_position_tolerance = 1.0d-9
  integer, parameter :: scan_count = 2*recovery_scan_half+1

  type, public :: recovery_event
    integer*8 :: index=0, reason=0, iteration=0, evaluations=0, brackets=0, iterations=0
    double precision :: initial_height=0, last_height=0, rejected_height=0
    double precision :: base_point(3)=0, normal(3)=0, initial_point(3)=0, last_point(3)=0
    double precision :: rejected_point(3)=0, last_residual=0, last_gradient(3)=0
    double precision :: sigma=0, lo=0, hi=0, root=0, residual=0, slope=0, position_tolerance=0
    logical :: recovered=.false.
    character(96) :: outcome='not_attempted'
    ! All detected downward brackets, including candidates rejected in validation.
    double precision, allocatable :: bracket_lower(:), bracket_upper(:)
    double precision, allocatable :: bracket_lower_residual(:), bracket_upper_residual(:)
    double precision, allocatable :: candidate_root(:), candidate_residual(:), candidate_slope(:)
    integer*8, allocatable :: candidate_iterations(:)
    character(96), allocatable :: candidate_outcome(:)
  end type

  ! Reverse-communication Brent-Dekker state.  status: 0 running, 1 success,
  ! 2 invalid input/evaluation, 3 iteration limit.  No evaluator is embedded.
  type, public :: brent_state
    double precision :: a=0, b=0, c=0, fa=0, fb=0, fc=0, d=0, e=0
    double precision :: position_tolerance=0, residual_tolerance=0, trial=0
    double precision :: lower=0, upper=0, lower_value=0, upper_value=0
    integer :: iterations=0, max_iterations=recovery_max_iterations, status=2
    logical :: waiting=.false., exact_mode=.false.
  end type

  public :: brent_initialize, brent_request, brent_submit
  public :: find_downward_brackets, recovery_search_interval, recovery_root_is_valid
  public :: select_recovery_root, recover_levelset_targets

contains

  pure subroutine brent_initialize(state,lo,hi,flo,fhi,position_tol,residual_tol,max_iterations)
    type(brent_state), intent(out) :: state
    double precision, intent(in) :: lo,hi,flo,fhi,position_tol,residual_tol
    integer, optional, intent(in) :: max_iterations
    state=brent_state()
    if (present(max_iterations)) state%max_iterations=max_iterations
    if (.not.all(ieee_is_finite([lo,hi,flo,fhi,position_tol,residual_tol]))) return
    if (lo.ge.hi.or.position_tol.le.0.or.residual_tol.lt.0.or.state%max_iterations.lt.1) return
    if (.not.((flo.gt.0.and.fhi.lt.0).or.(flo.lt.0.and.fhi.gt.0))) return
    state%a=lo; state%b=hi; state%c=lo
    state%fa=flo; state%fb=fhi; state%fc=flo
    state%d=hi-lo; state%e=state%d
    state%lower=lo; state%upper=hi
    state%lower_value=flo; state%upper_value=fhi
    state%position_tolerance=position_tol
    state%residual_tolerance=residual_tol
    state%status=0
  end subroutine

  pure subroutine brent_request(state,trial,requested)
    type(brent_state), intent(inout) :: state
    double precision, intent(out) :: trial
    logical, intent(out) :: requested
    double precision :: xm,tol,p,q,r,s,oldb,oldfb,min1,min2
    requested=.false.; trial=state%b
    if (state%status.ne.0) return
    if (state%waiting) then
      trial=state%trial; requested=.true.; return
    endif
    if (state%fb.eq.0.0d0) state%exact_mode=.true.
    if (state%exact_mode) then
      ! Even an exactly zero rounded function value must establish physical
      ! position certainty.  Probe on both sides until its sign bracket is
      ! narrow enough; a zero plateau cannot provide that certificate.
      if (state%upper-state%lower.le.state%position_tolerance) then
        state%c=state%upper
        if (state%b-state%lower.gt.state%upper-state%b) state%c=state%lower
        state%status=1; trial=state%b; return
      endif
      if (state%iterations.ge.state%max_iterations) then
        state%status=3; return
      endif
      if (state%b-state%lower.ge.state%upper-state%b) then
        trial=state%lower+0.5d0*(state%b-state%lower)
      else
        trial=state%b+0.5d0*(state%upper-state%b)
      endif
      if (trial.eq.state%b.or.trial.le.state%lower.or.trial.ge.state%upper) then
        state%status=2; return
      endif
      state%trial=trial; state%waiting=.true.; requested=.true.; return
    endif
    if ((state%fb.gt.0.and.state%fc.gt.0).or.(state%fb.lt.0.and.state%fc.lt.0)) then
      state%c=state%a; state%fc=state%fa
      state%d=state%b-state%a; state%e=state%d
    endif
    if (abs(state%fc).lt.abs(state%fb)) then
      oldb=state%b; oldfb=state%fb
      state%a=oldb; state%fa=oldfb
      state%b=state%c; state%fb=state%fc
      state%c=oldb; state%fc=oldfb
    endif
    xm=0.5d0*(state%c-state%b)
    tol=0.5d0*state%position_tolerance
    if (abs(xm).le.tol.and.abs(state%fb).le.state%residual_tolerance) then
      state%status=1; trial=state%b; return
    endif
    if (state%iterations.ge.state%max_iterations) then
      state%status=3; return
    endif
    if (abs(xm).le.tol) then
      ! The function has not met its residual tolerance despite a tiny bracket.
      ! Repeated rounded evaluations cannot establish a trustworthy root.
      state%status=2; return
    endif
    ! Nearly multiple roots can stall the usual interpolation criterion.
    ! Alternate interpolation with bisection to bound bracket contraction.
    ! Forty bisections fit within the fixed 80-evaluation recovery budget.
    if (mod(state%iterations,2).eq.1) then
      state%d=xm; state%e=xm
    elseif (abs(state%e).ge.tol.and.abs(state%fa).gt.abs(state%fb)) then
      s=state%fb/state%fa
      if (state%a.eq.state%c) then
        p=2.0d0*xm*s; q=1.0d0-s
      else
        q=state%fa/state%fc; r=state%fb/state%fc
        p=s*(2.0d0*xm*q*(q-r)-(state%b-state%a)*(r-1.0d0))
        q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
      endif
      if (p.gt.0) q=-q
      p=abs(p)
      min1=3.0d0*xm*q-abs(tol*q); min2=abs(state%e*q)
      if (ieee_is_finite(p).and.ieee_is_finite(q).and.q.ne.0.and.2.0d0*p.lt.min(min1,min2)) then
        state%e=state%d; state%d=p/q
      else
        state%d=xm; state%e=xm
      endif
    else
      state%d=xm; state%e=xm
    endif
    state%a=state%b; state%fa=state%fb
    if (abs(state%d).gt.tol) then
      trial=state%b+state%d
    else
      trial=state%b+sign(tol,xm)
    endif
    ! Keep every requested evaluation strictly inside the current bracket.
    if (.not.ieee_is_finite(trial).or.trial.le.min(state%b,state%c).or.trial.ge.max(state%b,state%c)) &
      trial=state%b+xm
    if (trial.eq.state%b.or.trial.eq.state%c) then
      state%status=2; return
    endif
    state%trial=trial; state%waiting=.true.; requested=.true.
  end subroutine

  pure subroutine brent_submit(state,value)
    type(brent_state), intent(inout) :: state
    double precision, intent(in) :: value
    if (state%status.ne.0.or..not.state%waiting) return
    state%waiting=.false.; state%iterations=state%iterations+1
    if (.not.ieee_is_finite(value)) then
      state%status=2; return
    endif
    if (state%exact_mode) then
      if (value.eq.0.0d0) then
        state%status=2; return
      endif
      if (state%trial.lt.state%b) then
        if ((value.gt.0).neqv.(state%lower_value.gt.0)) then
          state%status=2; return
        endif
        state%lower=state%trial; state%lower_value=value
      else
        if ((value.gt.0).neqv.(state%upper_value.gt.0)) then
          state%status=2; return
        endif
        state%upper=state%trial; state%upper_value=value
      endif
      return
    endif
    if (value.ne.0.0d0) then
      if ((value.gt.0).eqv.(state%lower_value.gt.0)) then
        state%lower=state%trial; state%lower_value=value
      else
        state%upper=state%trial; state%upper_value=value
      endif
    endif
    state%b=state%trial; state%fb=value
  end subroutine

  ! Ignore runs of near-zero values, requiring reliable opposite signs on
  ! their two sides.  This produces one bracket for a near-zero cluster.
  pure subroutine find_downward_brackets(x,values,zero_tol,left,right,nb)
    double precision, intent(in) :: x(:),values(:),zero_tol
    integer, intent(out) :: left(:),right(:),nb
    integer :: i,positive
    nb=0; positive=0; left=0; right=0
    if (size(x).ne.size(values)) return
    do i=1,size(x)
      if (.not.ieee_is_finite(values(i)).or..not.ieee_is_finite(x(i))) then
        positive=0; cycle
      endif
      if (values(i).gt.zero_tol) then
        positive=i
      elseif (values(i).lt.-zero_tol) then
        if (positive.gt.0) then
          if (x(i).gt.x(positive).and.nb.lt.min(size(left),size(right))) then
            nb=nb+1; left(nb)=positive; right(nb)=i
          endif
        endif
        positive=0
      endif
    enddo
  end subroutine

  pure subroutine recovery_search_interval(base,normal,h_initial,anchor_sigma,center,radius, &
      unit,s_anchor,lo,hi,position_tol,valid)
    double precision, intent(in) :: base(3),normal(3),h_initial,anchor_sigma,center(3),radius
    double precision, intent(out) :: unit(3),s_anchor,lo,hi,position_tol
    logical, intent(out) :: valid
    double precision :: normal_length,anchor(3),q(3),along,perp(3),remaining,width,scale,tolerance
    valid=.false.; unit=0; s_anchor=0; lo=0; hi=0; position_tol=0
    if (.not.all(ieee_is_finite(base)).or..not.all(ieee_is_finite(normal))) return
    if (.not.all(ieee_is_finite([h_initial,anchor_sigma,radius])).or..not.all(ieee_is_finite(center))) return
    if (radius.le.0.or.anchor_sigma.le.0) return
    normal_length=norm2(normal)
    if (.not.ieee_is_finite(normal_length).or.normal_length.le.0) return
    unit=normal/normal_length; s_anchor=h_initial*normal_length
    anchor=base+s_anchor*unit
    if (.not.all(ieee_is_finite(anchor))) return
    if (norm2((anchor-center)/radius).gt.10.0d0) return
    ! Work in radius-scaled coordinates to avoid overflowing the sphere's
    ! quadratic discriminant on otherwise representable inputs.
    q=(base-center)/radius; along=dot_product(q,unit); perp=q-along*unit
    remaining=100.0d0-dot_product(perp,perp)
    if (.not.ieee_is_finite(remaining).or.remaining.lt.0) return
    width=sqrt(remaining)
    lo=max(s_anchor-recovery_sigma_extent*anchor_sigma,radius*(-along-width))
    hi=min(s_anchor+recovery_sigma_extent*anchor_sigma,radius*(-along+width))
    if (.not.all(ieee_is_finite([lo,hi,s_anchor])).or.lo.ge.hi) return
    ! Analytic sphere endpoints can round a few ulps outside 10R.  Move the
    ! clipped endpoint inward by only the floating-point accuracy floor.
    scale=max(radius,maxval(abs(base)),maxval(abs(center)),abs(s_anchor),anchor_sigma)
    tolerance=64.0d0*epsilon(1.0d0)*scale
    if (norm2((base+lo*unit-center)/radius).gt.10.0d0) lo=min(s_anchor,lo+tolerance)
    if (norm2((base+hi*unit-center)/radius).gt.10.0d0) hi=max(s_anchor,hi-tolerance)
    position_tol=max(recovery_relative_position_tolerance*anchor_sigma,tolerance)
    if (lo.ge.hi.or.position_tol.le.0) return
    valid=.true.
  end subroutine

  pure logical function recovery_root_is_valid(point,residual,gradient,unit,s,lo,hi,center,radius, &
      lower_residual,upper_residual)
    double precision, intent(in) :: point(3),residual,gradient(3),unit(3),s,lo,hi,center(3),radius
    double precision, intent(in) :: lower_residual,upper_residual
    double precision :: gradient_length
    recovery_root_is_valid=.false.
    if (.not.all(ieee_is_finite(point)).or..not.all(ieee_is_finite(gradient))) return
    if (.not.all(ieee_is_finite(unit)).or..not.all(ieee_is_finite(center))) return
    if (.not.all(ieee_is_finite([residual,s,lo,hi,radius])).or.radius.le.0) return
    if (s.lt.lo.or.s.gt.hi.or.norm2((point-center)/radius).gt.10.0d0) return
    if (abs(residual).gt.recovery_residual_tolerance) return
    gradient_length=norm2(gradient)
    if (.not.ieee_is_finite(gradient_length).or.gradient_length.le.0) return
    if (gradient_length.le.1.0d-14/radius) return
    if (dot_product(gradient,unit).ge.-1.0d-12*gradient_length) return
    if (.not.(lower_residual.gt.0.and.upper_residual.lt.0)) return
    recovery_root_is_valid=.true.
  end function

  pure subroutine select_recovery_root(roots,valid,s_anchor,position_tol,selected,ambiguous)
    double precision, intent(in) :: roots(:),s_anchor,position_tol
    logical, intent(in) :: valid(:)
    integer, intent(out) :: selected
    logical, intent(out) :: ambiguous
    double precision :: distance,best
    integer :: i
    selected=0; ambiguous=.false.; best=huge(1.0d0)
    do i=1,min(size(roots),size(valid))
      if (.not.valid(i)) cycle
      if (.not.ieee_is_finite(roots(i))) cycle
      distance=abs(roots(i)-s_anchor)
      if (distance.lt.best-position_tol) then
        selected=i; best=distance; ambiguous=.false.
      elseif (abs(distance-best).le.position_tol.and.selected.ne.0) then
        if (abs(roots(i)-roots(selected)).gt.position_tol) ambiguous=.true.
      endif
    enddo
    if (ambiguous) selected=0
  end subroutine

  subroutine recover_levelset_targets(n,base,normals,h_initial,heights,deferred,grad,r_t, &
      geometry1,feval,mode,center,radius,events,ier)
    integer*8, intent(in) :: n,mode
    double precision, intent(in) :: base(3,n),normals(3,n),h_initial(n),center(3),radius
    logical, intent(in) :: deferred(n)
    double precision, intent(inout) :: heights(n),grad(3,n),r_t(3,n)
    type(Geometry), intent(in) :: geometry1
    type(Feval_stuff), intent(inout) :: feval
    type(recovery_event), intent(inout) :: events(:)
    integer*8, intent(out) :: ier
    integer*8, allocatable :: indices(:)
    logical, allocatable :: seen(:)
    integer :: group_start,group_end,group_size,j,ndeferred
    ier=0; ndeferred=count(deferred)
    if (size(events).ne.ndeferred) then
      ier=7
      events%outcome='invalid_event_count'
      return
    endif
    if (ndeferred.eq.0) return
    ! Events are compact: one record per deferred target in deterministic
    ! deferral order, rather than a large record for every ordinary target.
    allocate(indices(ndeferred),seen(n)); seen=.false.
    do j=1,ndeferred
      indices(j)=events(j)%index
      if (indices(j).lt.1.or.indices(j).gt.n) then
        events(j)%outcome='invalid_event_target_index'; ier=7; cycle
      endif
      if (.not.deferred(indices(j)).or.seen(indices(j))) then
        events(j)%outcome='invalid_or_duplicate_deferred_target'; ier=7; cycle
      endif
      seen(indices(j))=.true.
    enddo
    if (ier.ne.0) return
    deallocate(seen)
    ! Each group's entire scan and all its Brent trials fit the batch ceiling.
    group_size=recovery_batch_limit/scan_count
    do group_start=1,ndeferred,group_size
      group_end=min(ndeferred,group_start+group_size-1)
      call recover_group(indices(group_start:group_end),n,base,normals,h_initial,heights,grad,r_t, &
        geometry1,feval,mode,center,radius,events,group_start-1)
    enddo
    do j=1,ndeferred
      if (.not.events(j)%recovered) ier=7
    enddo
  end subroutine

  subroutine recover_group(indices,n,base,normals,h_initial,heights,grad,r_t,geometry1,feval,mode, &
      center,radius,events,event_offset)
    integer*8, intent(in) :: indices(:),n,mode
    double precision, intent(in) :: base(3,n),normals(3,n),h_initial(n),center(3),radius
    double precision, intent(inout) :: heights(n),grad(3,n),r_t(3,n)
    type(Geometry), intent(in) :: geometry1
    type(Feval_stuff), intent(inout) :: feval
    type(recovery_event), intent(inout) :: events(:)
    integer, intent(in) :: event_offset
    integer :: ng,j,k,nt,nb,q,selected,total_brackets,first,last,requested_count
    integer*8 :: idx,eidx
    double precision, allocatable :: points(:,:),directions(:,:),values(:),field_grad(:,:),sigma(:)
    double precision, allocatable :: unit(:,:),anchor_s(:),xs(:,:),scan_values(:,:)
    double precision, allocatable :: root_grad(:,:),root_values(:)
    logical, allocatable :: searchable(:),root_valid(:)
    integer, allocatable :: scan_owner(:),scan_slot(:),owner(:),local_slot(:),pending(:),begin_bracket(:),end_bracket(:)
    type(brent_state), allocatable :: states(:)
    integer :: left(scan_count),right(scan_count)
    logical :: valid,requested,ambiguous
    double precision :: lo,hi,tol,trial,gradient_length
    ng=size(indices)
    allocate(points(3,recovery_batch_limit),directions(3,recovery_batch_limit))
    allocate(values(recovery_batch_limit),field_grad(3,recovery_batch_limit))
    allocate(sigma(ng),unit(3,ng),anchor_s(ng),xs(scan_count,ng),scan_values(scan_count,ng))
    allocate(searchable(ng),scan_owner(recovery_batch_limit),scan_slot(recovery_batch_limit))
    allocate(begin_bracket(ng),end_bracket(ng))
    searchable=.false.; scan_values=ieee_value(0.0d0,ieee_quiet_nan)
    begin_bracket=0; end_bracket=-1; nt=0

    ! Reject invalid anchors before even the sigma evaluator sees them.
    do j=1,ng
      idx=indices(j); eidx=event_offset+j; events(eidx)%outcome='invalid_anchor'; events(eidx)%recovered=.false.
      if (.not.all(ieee_is_finite(base(:,idx))).or..not.all(ieee_is_finite(normals(:,idx)))) cycle
      if (.not.ieee_is_finite(h_initial(idx))) cycle
      gradient_length=norm2(normals(:,idx))
      if (.not.ieee_is_finite(gradient_length).or.gradient_length.le.0) cycle
      points(:,nt+1)=base(:,idx)+h_initial(idx)*normals(:,idx)
      if (.not.all(ieee_is_finite(points(:,nt+1)))) cycle
      if (.not.ieee_is_finite(radius).or.radius.le.0) cycle
      if (.not.all(ieee_is_finite(center))) cycle
      if (norm2((points(:,nt+1)-center)/radius).gt.10.0d0) cycle
      nt=nt+1; scan_owner(nt)=j; searchable(j)=.true.
    enddo
    if (nt.eq.0) return
    call function_eval_sigma(feval%FSS_1,points(:,1:nt),int(nt,8),values(1:nt), &
      field_grad(1,1:nt),field_grad(2,1:nt),field_grad(3,1:nt),mode)
    do k=1,nt
      j=scan_owner(k); sigma(j)=values(k)
    enddo
    nt=0
    do j=1,ng
      if (.not.searchable(j)) cycle
      idx=indices(j); eidx=event_offset+j
      events(eidx)%sigma=sigma(j); events(eidx)%outcome='invalid_sigma_or_interval'
      call recovery_search_interval(base(:,idx),normals(:,idx),h_initial(idx),sigma(j),center,radius, &
        unit(:,j),anchor_s(j),lo,hi,tol,valid)
      searchable(j)=valid
      if (.not.valid) cycle
      events(eidx)%lo=lo; events(eidx)%hi=hi; events(eidx)%position_tolerance=tol
      events(eidx)%outcome='no_downward_bracket'
      do k=1,scan_count
        if (k.le.recovery_scan_half+1) then
          xs(k,j)=lo+(anchor_s(j)-lo)*dble(k-1)/recovery_scan_half
        else
          xs(k,j)=anchor_s(j)+(hi-anchor_s(j))*dble(k-recovery_scan_half-1)/recovery_scan_half
        endif
        ! Preserve the anchor and endpoints exactly, avoiding interpolation ulps.
        if (k.eq.1) xs(k,j)=lo
        if (k.eq.recovery_scan_half+1) xs(k,j)=anchor_s(j)
        if (k.eq.scan_count) xs(k,j)=hi
        points(:,nt+1)=base(:,idx)+xs(k,j)*unit(:,j)
        if (.not.all(ieee_is_finite(points(:,nt+1)))) cycle
        if (norm2((points(:,nt+1)-center)/radius).gt.10.0d0) cycle
        nt=nt+1; directions(:,nt)=unit(:,j); scan_owner(nt)=j; scan_slot(nt)=k
      enddo
    enddo
    if (nt.eq.0) return
    call eval_density_grad_FMM(geometry1,points(:,1:nt),directions(:,1:nt),int(nt,8), &
      values(1:nt),field_grad(:,1:nt),feval,mode)
    do k=1,nt
      j=scan_owner(k); scan_values(scan_slot(k),j)=values(k)
      idx=indices(j); eidx=event_offset+j; events(eidx)%evaluations=events(eidx)%evaluations+1
    enddo

    total_brackets=0
    do j=1,ng
      if (.not.searchable(j)) cycle
      idx=indices(j); eidx=event_offset+j
      if (.not.all(ieee_is_finite(scan_values(:,j)))) then
        events(eidx)%outcome='nonfinite_scan'; searchable(j)=.false.; cycle
      endif
      call find_downward_brackets(xs(:,j),scan_values(:,j),recovery_residual_tolerance,left,right,nb)
      events(eidx)%brackets=nb
      call allocate_candidate_report(events(eidx),nb)
      if (nb.eq.0) cycle
      begin_bracket(j)=total_brackets+1; end_bracket(j)=total_brackets+nb
      total_brackets=total_brackets+nb
      do k=1,nb
        events(eidx)%bracket_lower(k)=xs(left(k),j); events(eidx)%bracket_upper(k)=xs(right(k),j)
        events(eidx)%bracket_lower_residual(k)=scan_values(left(k),j)
        events(eidx)%bracket_upper_residual(k)=scan_values(right(k),j)
      enddo
      events(eidx)%outcome='no_valid_downward_root'
    enddo
    if (total_brackets.eq.0) return
    allocate(states(total_brackets),owner(total_brackets),local_slot(total_brackets),pending(total_brackets))
    allocate(root_grad(3,total_brackets),root_values(total_brackets),root_valid(total_brackets))
    root_valid=.false.; root_grad=0; root_values=huge(1.0d0)
    do j=1,ng
      idx=indices(j); eidx=event_offset+j
      do q=begin_bracket(j),end_bracket(j)
        k=q-begin_bracket(j)+1; owner(q)=j; local_slot(q)=k
        call brent_initialize(states(q),events(eidx)%bracket_lower(k),events(eidx)%bracket_upper(k), &
          events(eidx)%bracket_lower_residual(k),events(eidx)%bracket_upper_residual(k), &
          events(eidx)%position_tolerance,recovery_residual_tolerance)
      enddo
    enddo
    do
      requested_count=0
      do q=1,total_brackets
        call brent_request(states(q),trial,requested)
        if (.not.requested) cycle
        j=owner(q); idx=indices(j); eidx=event_offset+j
        points(:,requested_count+1)=base(:,idx)+trial*unit(:,j)
        if (.not.all(ieee_is_finite(points(:,requested_count+1)))) then
          states(q)%status=2; cycle
        endif
        if (norm2((points(:,requested_count+1)-center)/radius).gt.10.0d0) then
          states(q)%status=2; cycle
        endif
        requested_count=requested_count+1; pending(requested_count)=q
        directions(:,requested_count)=unit(:,j)
      enddo
      if (requested_count.eq.0) exit
      call eval_density_grad_FMM(geometry1,points(:,1:requested_count),directions(:,1:requested_count), &
        int(requested_count,8),values(1:requested_count),field_grad(:,1:requested_count),feval,mode)
      do k=1,requested_count
        q=pending(k); idx=indices(owner(q)); eidx=event_offset+owner(q)
        call brent_submit(states(q),values(k))
        events(eidx)%evaluations=events(eidx)%evaluations+1
      enddo
    enddo

    ! Freshly evaluate every converged candidate, not stale interpolation data.
    nt=0
    do q=1,total_brackets
      j=owner(q); idx=indices(j); eidx=event_offset+j; k=local_slot(q)
      events(eidx)%iterations=events(eidx)%iterations+states(q)%iterations
      events(eidx)%candidate_iterations(k)=states(q)%iterations
      events(eidx)%candidate_root(k)=states(q)%b
      events(eidx)%candidate_residual(k)=states(q)%fb
      events(eidx)%candidate_outcome(k)='brent_invalid_or_unresolved'
      if (states(q)%status.eq.3) events(eidx)%candidate_outcome(k)='brent_iteration_limit'
      if (states(q)%status.ne.1) cycle
      ! Evaluate the precise original-parameterization coordinates returned
      ! to callers, rather than the algebraically equivalent unit-line form.
      trial=states(q)%b/norm2(normals(:,idx))
      if (.not.ieee_is_finite(trial)) then
        events(eidx)%candidate_outcome(k)='nonfinite_recovered_height'; cycle
      endif
      points(:,nt+1)=base(:,idx)+normals(:,idx)*trial
      if (.not.all(ieee_is_finite(points(:,nt+1)))) then
        events(eidx)%candidate_outcome(k)='nonfinite_recovered_point'; cycle
      endif
      if (norm2((points(:,nt+1)-center)/radius).gt.10.0d0) then
        events(eidx)%candidate_outcome(k)='recovered_point_outside_guard'; cycle
      endif
      nt=nt+1; pending(nt)=q
      directions(:,nt)=unit(:,j)
    enddo
    if (nt.gt.0) then
      call eval_density_grad_FMM(geometry1,points(:,1:nt),directions(:,1:nt),int(nt,8), &
        values(1:nt),field_grad(:,1:nt),feval,mode)
      do k=1,nt
        q=pending(k); j=owner(q); idx=indices(j); eidx=event_offset+j; nb=local_slot(q)
        root_values(q)=values(k); root_grad(:,q)=field_grad(:,k)
        events(eidx)%evaluations=events(eidx)%evaluations+1
        events(eidx)%candidate_residual(nb)=values(k)
        gradient_length=norm2(field_grad(:,k))
        if (ieee_is_finite(gradient_length).and.gradient_length.gt.0) &
          events(eidx)%candidate_slope(nb)=dot_product(field_grad(:,k),unit(:,j))/gradient_length
        root_valid(q)=recovery_root_is_valid(points(:,k),values(k),field_grad(:,k),unit(:,j),states(q)%b, &
          events(eidx)%lo,events(eidx)%hi,center,radius,states(q)%lower_value,states(q)%upper_value)
        events(eidx)%candidate_outcome(nb)='root_validation_failed'
        trial=states(q)%b/norm2(normals(:,idx))
        if (.not.ieee_is_finite(trial)) then
          root_valid(q)=.false.; events(eidx)%candidate_outcome(nb)='nonfinite_recovered_height'
        endif
        if (root_valid(q)) events(eidx)%candidate_outcome(nb)='valid_downward_root'
      enddo
    endif
    do j=1,ng
      first=begin_bracket(j); last=end_bracket(j)
      if (first.le.0.or.last.lt.first) cycle
      idx=indices(j); eidx=event_offset+j
      call select_recovery_root(events(eidx)%candidate_root,root_valid(first:last),anchor_s(j), &
        events(eidx)%position_tolerance,selected,ambiguous)
      if (ambiguous) then
        events(eidx)%outcome='ambiguous_equidistant_roots'; cycle
      endif
      if (selected.eq.0) cycle
      q=first+selected-1
      trial=states(q)%b/norm2(normals(:,idx))
      if (.not.ieee_is_finite(trial)) then
        events(eidx)%outcome='nonfinite_recovered_height'; cycle
      endif
      events(eidx)%root=states(q)%b; events(eidx)%residual=root_values(q)
      events(eidx)%slope=events(eidx)%candidate_slope(selected)
      heights(idx)=trial
      grad(:,idx)=root_grad(:,q); r_t(:,idx)=base(:,idx)+normals(:,idx)*heights(idx)
      events(eidx)%recovered=.true.; events(eidx)%outcome='recovered'
    enddo
  end subroutine

  subroutine allocate_candidate_report(event,nb)
    type(recovery_event), intent(inout) :: event
    integer, intent(in) :: nb
    if (allocated(event%bracket_lower)) then
      deallocate(event%bracket_lower,event%bracket_upper,event%bracket_lower_residual,event%bracket_upper_residual)
      deallocate(event%candidate_root,event%candidate_residual,event%candidate_slope)
      deallocate(event%candidate_iterations,event%candidate_outcome)
    endif
    allocate(event%bracket_lower(nb),event%bracket_upper(nb))
    allocate(event%bracket_lower_residual(nb),event%bracket_upper_residual(nb))
    allocate(event%candidate_root(nb),event%candidate_residual(nb),event%candidate_slope(nb))
    allocate(event%candidate_iterations(nb),event%candidate_outcome(nb))
    event%candidate_root=0; event%candidate_residual=0; event%candidate_slope=0
    event%candidate_iterations=0; event%candidate_outcome='not_attempted'
  end subroutine

end module Mod_Newton_Recovery
