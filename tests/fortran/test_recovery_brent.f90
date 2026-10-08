program test_recovery_brent
  use Mod_Newton_Recovery
  use ModType_Smooth_Surface, only: Geometry
  use Mod_Feval, only: Feval_stuff
  use, intrinsic :: ieee_arithmetic
  implicit none
  type(brent_state) :: state
  type(Geometry) :: geometry1
  type(Feval_stuff) :: feval
  type(recovery_event) :: events(2)
  double precision :: x(257),y(257),trial,unit(3),anchor,lo,hi,tol
  double precision :: base(3),normal(3),center(3),roots(4),nan
  double precision :: bases(3,4),normals(3,4),initial(4),heights(4),grad(3,4),targets(3,4)
  integer*8 :: ier
  integer :: left(257),right(257),nb,i,selected,calls
  logical :: requested,valid,ambiguous,eligible(4)
  logical :: deferred(4)

  ! Reverse communication must stay bracketed and find roots for several
  ! conditioning patterns without reading a derivative.
  do i=1,8
    call solve_case(i)
  enddo
  call brent_initialize(state,0.0d0,1.0d0,1.0d0,2.0d0,1.0d-10,1.0d-8)
  call require(state%status.eq.2,'a sign-changing bracket is mandatory')
  call brent_initialize(state,0.0d0,1.0d0,1.0d0,-1.0d0,0.0d0,1.0d-8)
  call require(state%status.eq.2,'zero position tolerance is invalid')
  call brent_initialize(state,0.0d0,1.0d0,1.0d0,-2.0d0,1.0d-12,1.0d-8,1)
  call brent_request(state,trial,requested)
  call require(requested,'one iteration requests a value')
  call brent_submit(state,0.25d0)
  call brent_request(state,trial,requested)
  call require(state%status.eq.3.and..not.requested,'iteration exhaustion is explicit')
  nan=ieee_value(0.0d0,ieee_quiet_nan)
  call brent_initialize(state,0.0d0,1.0d0,1.0d0,-1.0d0,1.0d-10,1.0d-8)
  call brent_request(state,trial,requested)
  call brent_submit(state,nan)
  call require(state%status.eq.2,'nonfinite evaluations fail closed')
  call brent_initialize(state,-1.0d0,1.0d0,1.0d0,-1.0d0,1.0d-10,1.0d-8)
  do
    call brent_request(state,trial,requested)
    if (.not.requested) exit
    if (abs(trial).le.0.1d0) then
      call brent_submit(state,0.0d0)
    else
      call brent_submit(state,-trial)
    endif
  enddo
  call require(state%status.eq.2,'a rounded zero plateau cannot establish position certainty')

  do i=1,257
    x(i)=-4.0d0+dble(i-1)/32.0d0
    y(i)=-(x(i)+2.0d0)*x(i)*(x(i)-2.0d0)
  enddo
  call find_downward_brackets(x,y,1.0d-8,left,right,nb)
  call require(nb.eq.2,'two downward crossings detected; upward crossing excluded')
  call require(x(left(1)).lt.-2.and.x(right(1)).gt.-2,'first bracket crosses left root')
  call require(x(left(2)).lt.2.and.x(right(2)).gt.2,'second bracket crosses right root')
  y=-x*x
  call find_downward_brackets(x,y,1.0d-8,left,right,nb)
  call require(nb.eq.0,'tangent contact is not a downward bracket')
  y=x
  call find_downward_brackets(x,y,1.0d-8,left,right,nb)
  call require(nb.eq.0,'upward crossing rejected')
  y=1.0d0
  y(100:110)=0.5d-8
  y(111:)=-1.0d0
  call find_downward_brackets(x,y,1.0d-8,left,right,nb)
  call require(nb.eq.1.and.left(1).eq.99.and.right(1).eq.111,'near-zero run forms one bracket')
  y(105)=nan
  call find_downward_brackets(x,y,1.0d-8,left,right,nb)
  call require(nb.eq.0,'nonfinite gap cannot form a trustworthy bracket')
  y=1.0d0
  call find_downward_brackets(x,y,1.0d-8,left,right,nb)
  call require(nb.eq.0,'no crossing remains unresolved')

  roots=[-2.0d0,2.0d0,0.75d0,0.25d0]
  eligible=[.true.,.true.,.true.,.false.]
  call select_recovery_root(roots,eligible,0.0d0,1.0d-9,selected,ambiguous)
  call require(selected.eq.3.and..not.ambiguous,'nearest eligible root selected')
  eligible=[.true.,.true.,.false.,.false.]
  call select_recovery_root(roots,eligible,0.0d0,1.0d-9,selected,ambiguous)
  call require(selected.eq.0.and.ambiguous,'equidistant distinct roots fail as ambiguous')
  roots=[1.0d0,1.0d0+0.25d-9,3.0d0,4.0d0]
  call select_recovery_root(roots,eligible,0.0d0,1.0d-9,selected,ambiguous)
  call require(selected.eq.1.and..not.ambiguous,'duplicate root coordinates are not ambiguous')
  eligible=.false.
  call select_recovery_root(roots,eligible,0.0d0,1.0d-9,selected,ambiguous)
  call require(selected.eq.0.and..not.ambiguous,'no valid candidate fails closed')

  base=[0.0d0,0.0d0,0.0d0]; normal=[2.0d0,0.0d0,0.0d0]; center=0
  call recovery_search_interval(base,normal,0.25d0,0.1d0,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(valid.and.abs(anchor-0.5d0).lt.1.0d-15,'physical anchor uses normal length')
  call require(abs(lo+0.3d0).lt.1.0d-15.and.abs(hi-1.3d0).lt.1.0d-15,'local eight-sigma interval')
  call require(abs(tol-1.0d-10).lt.1.0d-20,'local position tolerance')
  call require(norm2(unit-[1.0d0,0.0d0,0.0d0]).lt.1.0d-15,'unit direction')
  normal=[1.0d0,0.0d0,0.0d0]
  call recovery_search_interval(base,normal,0.0d0,100.0d0,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(valid.and.abs(lo+10).lt.1.0d-12.and.abs(hi-10).lt.1.0d-12,'window clips sphere')
  base=[0.0d0,6.0d0,0.0d0]
  call recovery_search_interval(base,normal,0.0d0,100.0d0,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(valid.and.abs(lo+8).lt.1.0d-12.and.abs(hi-8).lt.1.0d-12,'off-center sphere intersection')
  base=0
  call recovery_search_interval(base,normal,10.0d0,0.1d0,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(valid.and.anchor.eq.hi,'anchor exactly at 10R leaves one-sided interval')
  do i=1,257
    if (i.le.129) then
      x(i)=lo+(anchor-lo)*dble(i-1)/128.0d0
    else
      x(i)=anchor
    endif
    y(i)=9.6d0-x(i)
  enddo
  call find_downward_brackets(x,y,1.0d-8,left,right,nb)
  call require(nb.eq.1,'repeated clipped endpoints do not duplicate a bracket')
  base=[0.0d0,11.0d0,0.0d0]
  call recovery_search_interval(base,normal,0.0d0,1.0d0,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(.not.valid,'outside anchor rejected')
  base=0; normal=0
  call recovery_search_interval(base,normal,0.0d0,1.0d0,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(.not.valid,'zero pseudonormal rejected')
  normal=[1.0d0,0.0d0,0.0d0]
  call recovery_search_interval(base,normal,0.0d0,0.0d0,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(.not.valid,'nonpositive sigma rejected')
  call recovery_search_interval(base,normal,0.0d0,nan,center,1.0d0,unit,anchor,lo,hi,tol,valid)
  call require(.not.valid,'nonfinite sigma rejected')
  ! Translation and scaling preserve the physical line and search radius.
  base=[30.0d0,-20.0d0,10.0d0]; center=base; normal=[3.0d0,0.0d0,0.0d0]
  call recovery_search_interval(base,normal,2.0d0/3.0d0,0.4d0,center,4.0d0, &
    unit,anchor,lo,hi,tol,valid)
  call require(valid.and.abs(anchor-2).lt.1.0d-14,'translated physical anchor')
  call require(abs(lo+1.2d0).lt.1.0d-14.and.abs(hi-5.2d0).lt.1.0d-14,'scaled search interval')
  base=0; center=0; normal=[1.0d-150,0.0d0,0.0d0]
  call recovery_search_interval(base,normal,1.0d150,0.05d0,center,1.0d0, &
    unit,anchor,lo,hi,tol,valid)
  call require(valid.and.abs(anchor-1).lt.1.0d-14,'tiny nonunit normal retains physical anchor')
  call require(abs(lo-0.6d0).lt.1.0d-14.and.abs(hi-1.4d0).lt.1.0d-14,'tiny normal search unchanged')
  normal=[1.0d150,0.0d0,0.0d0]
  call recovery_search_interval(base,normal,1.0d-150,0.05d0,center,1.0d0, &
    unit,anchor,lo,hi,tol,valid)
  call require(valid.and.abs(anchor-1).lt.1.0d-14,'huge nonunit normal retains physical anchor')
  base=[1.0d12,-2.0d12,3.0d12]; center=base; normal=[1.0d0,0.0d0,0.0d0]
  call recovery_search_interval(base,normal,0.5d0,0.1d0,center,1.0d0, &
    unit,anchor,lo,hi,tol,valid)
  call require(valid.and.abs(anchor-0.5d0).lt.1.0d-14,'large translation retains physical anchor')
  call require(tol.ge.64.0d0*epsilon(1.0d0)*3.0d12,'large translation uses representability floor')

  center=0; base=0; unit=[1.0d0,0.0d0,0.0d0]
  valid=recovery_root_is_valid(base,1.0d-10,[-1.0d0,0.0d0,0.0d0],unit,0.0d0,-1.0d0,1.0d0, &
    center,1.0d0,1.0d0,-1.0d0)
  call require(valid,'valid downward transverse root')
  valid=recovery_root_is_valid(base,0.0d0,[1.0d0,0.0d0,0.0d0],unit,0.0d0,-1.0d0,1.0d0, &
    center,1.0d0,1.0d0,-1.0d0)
  call require(.not.valid,'upward root rejected even with zero residual')
  valid=recovery_root_is_valid(base,0.0d0,[0.0d0,1.0d0,0.0d0],unit,0.0d0,-1.0d0,1.0d0, &
    center,1.0d0,1.0d0,-1.0d0)
  call require(.not.valid,'tangent root rejected')
  valid=recovery_root_is_valid(base,0.0d0,[-1.0d-15,0.0d0,0.0d0],unit,0.0d0,-1.0d0,1.0d0, &
    center,1.0d0,1.0d0,-1.0d0)
  call require(.not.valid,'near-zero gradient rejected')
  valid=recovery_root_is_valid(base,2.0d-8,[-1.0d0,0.0d0,0.0d0],unit,0.0d0,-1.0d0,1.0d0, &
    center,1.0d0,1.0d0,-1.0d0)
  call require(.not.valid,'fresh residual required')
  valid=recovery_root_is_valid(base,0.0d0,[-1.0d0,0.0d0,0.0d0],unit,0.0d0,-1.0d0,1.0d0, &
    center,1.0d0,-1.0d0,1.0d0)
  call require(.not.valid,'original bracket must establish downward crossing')

  ! Compact event records map back to original target IDs.  These deliberately
  ! invalid anchors must fail before the uninitialized FMM/sigma state is used.
  deferred=[.false.,.true.,.false.,.true.]
  bases=nan; normals=1.0d0; initial=0; heights=7.0d0; grad=8.0d0; targets=9.0d0
  events(1)%index=4; events(2)%index=2
  call recover_levelset_targets(4_8,bases,normals,initial,heights,deferred,grad,targets, &
    geometry1,feval,1_8,center,1.0d0,events,ier)
  call require(ier.eq.7.and.all(events%outcome.eq.'invalid_anchor'),'compact original indices fail before FMM')
  call require(all(heights.eq.7).and.all(grad.eq.8).and.all(targets.eq.9),'failed targets preserve retained data')
  events(1)%index=2; events(2)%index=2
  call recover_levelset_targets(4_8,bases,normals,initial,heights,deferred,grad,targets, &
    geometry1,feval,1_8,center,1.0d0,events,ier)
  call require(ier.eq.7.and.events(2)%outcome.eq.'invalid_or_duplicate_deferred_target','duplicate event rejected')
  events(1)%index=5; events(2)%index=2
  call recover_levelset_targets(4_8,bases,normals,initial,heights,deferred,grad,targets, &
    geometry1,feval,1_8,center,1.0d0,events,ier)
  call require(ier.eq.7.and.events(1)%outcome.eq.'invalid_event_target_index','out-of-range event rejected')
  call recover_levelset_targets(4_8,bases,normals,initial,heights,deferred,grad,targets, &
    geometry1,feval,1_8,center,1.0d0,events(1:1),ier)
  call require(ier.eq.7.and.events(1)%outcome.eq.'invalid_event_count','event count must match deferred targets')
  deferred=.false.
  call recover_levelset_targets(4_8,bases,normals,initial,heights,deferred,grad,targets, &
    geometry1,feval,1_8,center,1.0d0,events(1:0),ier)
  call require(ier.eq.0,'empty recovery returns without FMM')
  print *, 'PASS: derivative-free Brent, bracketing, sphere clipping, selection, and root validation'

contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ',message
      stop 1
    endif
  end subroutine

  double precision function value(case_number,s)
    integer, intent(in) :: case_number
    double precision, intent(in) :: s
    select case(case_number)
    case(1)
      value=0.25d0-s
    case(2)
      value=0.3d0-exp(s)
    case(3)
      value=cos(s)-s
    case(4)
      value=0.3d0-s*s*s
    case(5)
      value=-(s-0.17325d0)**3
    case(6)
      value=1.0d0/(1.0d0+exp(20.0d0*s))-0.314159d0
    case(7)
      value=1.0d-10*(0.125d0-s)
    case(8)
      value=(0.43125d0-s)*(1.0d0+s*s)
    end select
  end function

  subroutine solve_case(case_number)
    integer, intent(in) :: case_number
    double precision :: a,b,expected
    a=-2.0d0; b=2.0d0
    select case(case_number)
    case(1)
      expected=0.25d0
    case(2)
      expected=log(0.3d0)
    case(3)
      expected=0.7390851332151607d0
    case(4)
      expected=0.3d0**(1.0d0/3.0d0)
    case(5)
      expected=0.17325d0
    case(6)
      expected=log(1.0d0/0.314159d0-1.0d0)/20.0d0
    case(7)
      expected=0.125d0
    case(8)
      expected=0.43125d0
    end select
    call brent_initialize(state,a,b,value(case_number,a),value(case_number,b),1.0d-10,1.0d-8)
    calls=0
    do
      call brent_request(state,trial,requested)
      if (.not.requested) exit
      call require(trial.gt.a.and.trial.lt.b,'every trial lies inside input bracket')
      call brent_submit(state,value(case_number,trial))
      calls=calls+1
      call require(calls.le.80,'Brent respects evaluation limit')
    enddo
    if (state%status.ne.1.or.abs(state%b-expected).gt.2.0d-9) then
      print *, 'case/state/root/expected: ',case_number,state%status,state%b,expected
    endif
    call require(state%status.eq.1,'Brent converged')
    call require(abs(state%c-state%b).le.state%position_tolerance,'success requires narrow position bracket')
    call require(abs(state%b-expected).le.2.0d-9,'root has requested physical accuracy')
    call require(abs(value(case_number,state%b)).le.1.0d-8,'root meets residual tolerance')
  end subroutine
end program
