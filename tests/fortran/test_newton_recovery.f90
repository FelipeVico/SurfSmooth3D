program test_newton_recovery
  use Mod_Smooth_Surface
  implicit none
  type(Geometry) :: g
  type(Feval_stuff) :: fev
  type(projection_radius_guard) :: guard
  character(2048) :: mesh,cad,root
  integer*8 :: ier,calls0,deferred0,recovered0
  double precision :: bases(3,2),normals(3,2),projected(3,2),h(2),grad(3,2),f(2)
  double precision :: x(1),r(3,1),df(3,1),initial_grad(3),reference,sigma(1),ds(3,1)
  double precision :: original_sources(3,1),unsafe_proposal(3)
  logical :: exists

  call get_command_argument(1,mesh)
  call get_command_argument(2,cad)
  call get_command_argument(3,root)
  call require(len_trim(mesh)>0.and.len_trim(cad)>0.and.len_trim(root)>0,'three paths required')
  call readgeometry(g,trim(mesh),3_8,8_8,4_8,ier)
  call require(ier==0,'read cylinder scaffold')
  call load_cad_skeleton(g,trim(cad))
  call funcion_normal_vert(g)
  call initialize_projection_guard(guard,g,trim(root),ier)
  call require(ier==0,'initialize fixed radius')
  call start_Feval_tree(fev,g,5d0,1_8)
  fev%audit_targets=.true.
  fev%audit_center=guard%center
  fev%audit_radius=guard%radius
  original_sources=g%skeleton_Points(:,1:1)
  call check_candidate_guard()

  ! The first point is just inside the smoothed barrel. Its pseudonormal is
  ! almost circumferential, with a tiny outward gradient component. Newton's
  ! first proposal leaves 10R, although a local outward crossing exists.
  bases=0d0
  bases(1,:)=[9.9d0,9.5d0]
  normals=0d0
  normals(1,:)=1d0
  call eval_density_grad_FMM(g,bases,normals,2_8,f,grad,fev,1_8)
  initial_grad=grad(:,1)
  call tangent_direction(initial_grad,normals(:,1))
  call require(f(1)>0d0,'recovery anchor is inside the fixed level set')
  unsafe_proposal=bases(:,1)-f(1)/dot_product(initial_grad,normals(:,1))*normals(:,1)
  call require(norm2((unsafe_proposal-guard%center)/guard%radius)>10d0,'first proposal exceeds 10R')
  call function_eval_sigma(fev%FSS_1,bases(:,1:1),1_8,sigma,ds(1,:),ds(2,:),ds(3,:),1_8)
  reference=independent_root(bases(:,1),normals(:,1),8d0*sigma(1))

  deferred0=guard%recovery_deferred
  recovered0=guard%recovery_recovered
  call project_points_to_levelset(g,fev,1_8,2_8,bases,normals,projected,h,grad,ier,guard)
  call require(ier==0,'mixed ordinary/deferred vertex targets succeed')
  call require(guard%recovery_deferred==deferred0+1,'exactly one target deferred')
  call require(guard%recovery_recovered==recovered0+1,'exactly one target recovered')
  call require(guard%recovered_mask(1).and..not.guard%recovered_mask(2),'stable target scattering')
  call require(abs(h(1)-reference)<1d-6,'recovery agrees with independent direct-kernel bisection')
  call require(h(1)>0d0.and.h(1)<8d0*sigma(1),'closest outward crossing stays local')
  call require(.not.fev%unsafe_target_seen,'unsafe Newton proposal never reached FMM')
  call eval_density_grad_FMM(g,projected,normals,2_8,f,grad,fev,1_8)
  call require(maxval(abs(f))<=1d-8,'all accepted roots have fresh small residuals')
  call require(dot_product(grad(:,1),normals(:,1))<0d0,'accepted root decreases outward')
  inquire(file=trim(root)//'_newton_recovery.txt',exist=exists)
  call require(exists,'recovery report persists')
  print *, 'PASS: mixed targets recover after ordinary Newton; all FMM targets stay safe.'
  call check_nonleading_target_indices()

  ! Exercise an all-deferred solve, which must never issue an empty FMM call.
  deferred0=guard%recovery_deferred
  call project_points_to_levelset(g,fev,1_8,1_8,bases(:,1:1),normals(:,1:1), &
      r,x,df,ier,guard)
  call require(ier==0.and.guard%recovery_deferred==deferred0+1,'all-deferred vertex solve succeeds')
  call require(abs(x(1)-reference)<1d-6,'all-deferred solve selects the same root')

  ! A deferred point must not start recovery when another active point hits
  ! the unchanged Newton iteration limit.
  if(allocated(g%Base_Points)) deallocate(g%Base_Points)
  if(allocated(g%Base_Points_N)) deallocate(g%Base_Points_N)
  allocate(g%Base_Points(3,2),g%Base_Points_N(3,2))
  g%n_Sf_points=2
  g%Base_Points=bases;g%Base_Points_N=normals
  h=0d0;projected=bases;grad=0d0
  recovered0=guard%recovery_recovered
  calls0=fev%evaluation_calls
  call My_Newton(h,1d-9,2_8,g,ier,fev,1_8,grad,projected,guard)
  call require(ier==1,'remaining active targets retain the existing iteration limit')
  call require(guard%recovery_recovered==recovered0,'active-group failure prevents deferred recovery')
  call require(h(1)==0d0,'unsafe proposal does not overwrite the deferred safe height')
  call require(fev%evaluation_calls==calls0+1,'active-group failure issues no recovery field evaluations')

  ! My_Newton also serves arbitrary independent/adaptive-blend targets. It
  ! must not assume that n_Sf_points is a multiple of the patch node count.
  if(allocated(g%Base_Points)) deallocate(g%Base_Points)
  if(allocated(g%Base_Points_N)) deallocate(g%Base_Points_N)
  allocate(g%Base_Points(3,1),g%Base_Points_N(3,1))
  g%n_Sf_points=1
  g%Base_Points=bases(:,1:1)
  g%Base_Points_N=normals(:,1:1)
  x=0d0;r=g%Base_Points;df=0d0
  deferred0=guard%recovery_deferred
  call My_Newton(x,1d-9,14_8,g,ier,fev,1_8,df,r,guard)
  call require(ier==0.and.guard%recovery_deferred==deferred0+1,'arbitrary-target Newton recovers')
  call require(abs(x(1)-reference)<1d-6,'arbitrary-target root agrees with reference')
  call require(maxval(abs(r(:,1)-g%Base_Points(:,1)-x(1)*g%Base_Points_N(:,1)))<1d-13, &
      'recovered coordinate uses original base and direction')
  call require(.not.fev%unsafe_target_seen,'all recovery evaluations stay inside 10R')
  call check_inherited_height()
  print *, 'PASS: all-deferred and arbitrary-target solves use the original projection line.'

  ! Disabled recovery retains the original guarded abort, before a bad target
  ! can reach the next FMM evaluation.
  guard%recovery_enabled=.false.
  deferred0=guard%recovery_deferred
  x=0d0;r=g%Base_Points;df=0d0
  call My_Newton(x,1d-9,14_8,g,ier,fev,1_8,df,r,guard)
  call require(ier==7,'disabled recovery retains radius failure')
  call require(guard%recovery_deferred==deferred0,'disabled recovery does not queue targets')
  call require(.not.fev%unsafe_target_seen,'disabled solve also blocks unsafe evaluations')
  guard%recovery_enabled=.true.

  ! Existing iteration-limit failures do not trigger recovery.
  x=0d0;r=g%Base_Points;df=0d0
  deferred0=guard%recovery_deferred
  calls0=fev%evaluation_calls
  call My_Newton(x,1d-9,1_8,g,ier,fev,1_8,df,r,guard)
  call require(ier==1,'ordinary iteration limit remains unchanged')
  call require(guard%recovery_deferred==deferred0,'iteration limit does not defer')
  call require(fev%evaluation_calls==calls0,'zero available iterations issue no FMM call')

  ! A circumferential line entirely outside the smoothed barrel has no local
  ! eligible crossing. Fail closed without publishing a bogus root.
  g%Base_Points(1,1)=9.99d0
  call eval_density_grad_FMM(g,g%Base_Points,g%Base_Points_N,1_8,f(1:1),df,fev,1_8)
  call tangent_direction(df(:,1),g%Base_Points_N(:,1))
  x=0d0;r=g%Base_Points
  recovered0=guard%recovery_recovered
  call My_Newton(x,1d-9,14_8,g,ier,fev,1_8,df,r,guard)
  call require(ier/=0,'missing local crossing fails closed')
  call require(guard%recovery_recovered==recovered0,'failed recovery publishes no successful target')
  call require(.not.fev%unsafe_target_seen,'failed recovery also evaluates only safe points')
  call require(all(g%skeleton_Points(:,1:1)==original_sources),'CAD sources unchanged by recovery')
  call check_invalid_patches()
  call check_transactional_scaffold()
  call destroy_Feval_tree(fev)
  print *, 'PASS: disabled recovery, ordinary iteration limits, and absent crossings fail safely.'
contains
  subroutine check_inherited_height()
    double precision :: anchor(3),expected(3),recorded_height,recorded_point(3)
    integer :: unit,ios
    character(1024) :: line
    anchor=g%Base_Points(:,1)
    expected=anchor+reference*g%Base_Points_N(:,1)
    g%Base_Points(:,1)=anchor-.1d0*g%Base_Points_N(:,1)
    x=.1d0;r(:,1)=anchor;df=0d0
    call My_Newton(x,1d-9,14_8,g,ier,fev,1_8,df,r,guard)
    call require(ier==0,'inherited nonzero height recovers')
    call require(abs(x(1)-(.1d0+reference))<1d-6,'inherited height retains the original parameterization')
    call require(maxval(abs(r(:,1)-expected))<1d-6,'inherited height recovers the same physical root')
    recorded_height=-1d0;recorded_point=huge(1d0)
    open(newunit=unit,file=trim(root)//'_newton_recovery.txt',status='old',iostat=ios)
    call require(ios==0,'open inherited recovery diagnostic')
    do
      read(unit,'(a)',iostat=ios) line
      if(ios/=0) exit
      line=adjustl(line)
      if(index(line,'INITIAL_HEIGHT ')==1) read(line(15:),*) recorded_height
      if(index(line,'INITIAL_POINT ')==1) read(line(14:),*) recorded_point
    enddo
    close(unit)
    call require(recorded_height==.1d0,'diagnostic preserves inherited initial height')
    call require(maxval(abs(recorded_point-anchor))<1d-13,'recovery window is centered on inherited initial point')
    g%Base_Points(:,1)=anchor
    print *, 'PASS: inherited nonzero height keeps its initial point and recovers the same physical root.'
  end subroutine

  subroutine check_transactional_scaffold()
    type(Geometry) :: scaffold
    type(projection_radius_guard) :: scaffold_guard
    double precision :: original(3,6)
    integer*8 :: status
    ! The almost tangential corner recovers across the opposite edge while
    ! its two radial neighbors converge normally. Each point has a valid
    ! root, but publishing the corners would reverse the triangle.
    scaffold=g
    deallocate(scaffold%Points,scaffold%Tri,scaffold%Normal_Vert)
    allocate(scaffold%Points(3,6),scaffold%Tri(6,1),scaffold%Normal_Vert(3,6))
    scaffold%ntri=1;scaffold%npoints=6
    scaffold%Tri(:,1)=[1_8,2_8,3_8,4_8,5_8,6_8]
    scaffold%Points(:,1)=bases(:,1)
    scaffold%Points(:,2)=bases(:,1)+[0d0,-.01d0,0d0]
    scaffold%Points(:,3)=bases(:,1)+[0d0,0d0,-.01d0]
    scaffold%Points(:,4)=(scaffold%Points(:,1)+scaffold%Points(:,2))/2d0
    scaffold%Points(:,5)=(scaffold%Points(:,2)+scaffold%Points(:,3))/2d0
    scaffold%Points(:,6)=(scaffold%Points(:,3)+scaffold%Points(:,1))/2d0
    scaffold%Normal_Vert=0d0;scaffold%Normal_Vert(1,:)=1d0
    scaffold%Normal_Vert(:,1)=normals(:,1)
    original=scaffold%Points
    call initialize_projection_guard(scaffold_guard,scaffold,trim(root)//'_transactional',status)
    call require(status==0,'initialize transactional scaffold guard')
    call project_scaffold_vertices_to_levelset(scaffold,fev,1_8,status,scaffold_guard)
    call require(status==6,'recovered scaffold orientation reversal is rejected')
    call require(scaffold_guard%recovery_recovered==1,'transactional test exercised actual corner recovery')
    call require(all(scaffold%Points==original),'rejected scaffold leaves original corners and midpoints unchanged')
    call require(.not.fev%unsafe_target_seen,'transactional validation also protects all FMM targets')
    print *, 'PASS: recovered scaffold orientation failure is transactional.'
  end subroutine

  subroutine check_nonleading_target_indices()
    double precision :: base(3,3),normal(3,3),points(3,3),height(3),gradient(3,3)
    integer*8 :: status,before_deferred,before_recovered
    base(:,1)=bases(:,2)
    base(:,2)=bases(:,1)
    base(:,3)=bases(:,1)
    normal(:,1)=normals(:,2)
    normal(:,2)=normals(:,1)
    normal(:,3)=2d0*normals(:,1)
    before_deferred=guard%recovery_deferred
    before_recovered=guard%recovery_recovered
    call project_points_to_levelset(g,fev,1_8,3_8,base,normal,points,height,gradient,status,guard)
    call require(status==0,'multiple nonleading deferred targets recover together')
    call require(guard%recovery_deferred==before_deferred+2.and. &
        guard%recovery_recovered==before_recovered+2,'multiple recovery counters are accurate')
    call require(all(guard%recovered_mask.eqv.[.false.,.true.,.true.]),'compact recovery scatters original indices')
    call require(abs(height(2)-reference)<1d-6.and.abs(2d0*height(3)-reference)<1d-6, &
        'physical displacement is independent of pseudonormal length')
    call require(maxval(abs(points(:,2)-points(:,3)))<1d-10,'scaled pseudonormals recover identical positions')
    print *, 'PASS: multiple nonleading target indices scatter correctly with unnormalized pseudonormals.'
  end subroutine

  subroutine check_candidate_guard()
    type(recovery_event) :: event
    type(recovery_event),allocatable :: events(:)
    double precision :: base(3),normal(3),gradient(3),candidate
    logical :: deferred
    integer*8 :: before,used,target
    base=guard%center;normal=[1d0,0d0,0d0];gradient=[-1d0,0d0,0d0]
    before=fev%evaluation_calls
    candidate=10d0*guard%radius
    call defer_projection_target(guard,1_8,2_8,base,normal,0d0,0.25d0,candidate,1d0,gradient,event,deferred)
    call require(.not.deferred,'exactly 10R candidate remains accepted')
    candidate=10d0*guard%radius*(1d0+4d0*epsilon(1d0))
    call defer_projection_target(guard,1_8,2_8,base,normal,0d0,0.25d0,candidate,1d0,gradient,event,deferred)
    call require(deferred.and.event%reason==1,'beyond-10R candidate is deferred')
    call require(event%initial_height==0d0.and.event%last_height==0.25d0,'initial and last safe heights retained')
    candidate=ieee_value(0d0,ieee_quiet_nan)
    call defer_projection_target(guard,1_8,3_8,base,normal,0d0,0.25d0,candidate,1d0,gradient,event,deferred)
    call require(deferred.and.event%reason==2,'nonfinite candidate is deferred')
    call require(event%last_residual==1d0.and.all(event%last_gradient==gradient),'last finite field data retained')
    call require(fev%evaluation_calls==before,'candidate rejection performs no FMM evaluation')
    used=0
    allocate(event%bracket_lower(1))
    do target=35,1,-1
      event%index=target
      event%bracket_lower=real(target,8)
      call append_recovery_event(events,event,used)
    enddo
    call require(used==35.and.size(events)>=35,'compact event storage grows beyond its initial capacity')
    event%bracket_lower=0d0
    do target=1,used
      call require(events(target)%index==36-target,'compact events retain original target identities')
      call require(events(target)%bracket_lower(1)==real(36-target,8),'event growth preserves owned bracket data')
    enddo
    print *, 'PASS: candidate checks preserve safe state and enforce exact 10R/nonfinite triggers.'
  end subroutine

  subroutine check_invalid_patches()
    type(Geometry) :: patch
    type(projection_radius_guard) :: patch_guard
    integer*8,parameter :: np=15
    integer*8 :: status,node,before
    double precision :: rv(2,np),umat(np,np),vmat(np,np),weights(np)
    call vioreanu_simplex_quad(4_8,np,rv,umat,vmat,weights)
    patch%ntri=1;patch%n_order_sf=np;patch%n_Sf_points=np;patch%norder_smooth=4
    patch%npoints=6
    allocate(patch%Tri(6,1),patch%Points(3,6),patch%Normal_Vert(3,6))
    allocate(patch%S_smooth(3,np),patch%N_smooth(3,np),patch%du_smooth(3,np), &
        patch%dv_smooth(3,np),patch%ds_smooth(np),patch%Base_Points_U(3,np),patch%Base_Points_V(3,np))
    patch%Tri(:,1)=[1_8,2_8,3_8,4_8,5_8,6_8]
    patch%Points(1,:)=[0d0,1d0,0d0,.5d0,.5d0,0d0]
    patch%Points(2,:)=[0d0,0d0,1d0,0d0,.5d0,.5d0]
    patch%Points(3,:)=0d0
    patch%Normal_Vert=0d0;patch%Normal_Vert(3,:)=1d0
    patch%Base_Points_U=0d0;patch%Base_Points_U(1,:)=1d0
    patch%Base_Points_V=0d0;patch%Base_Points_V(2,:)=1d0
    patch%N_smooth=0d0;patch%N_smooth(3,:)=1d0
    patch%du_smooth=patch%Base_Points_U
    patch%dv_smooth=-patch%Base_Points_V
    patch%ds_smooth=1d0
    patch%S_smooth=0d0
    patch%S_smooth(1:2,:)=rv
    patch_guard%enabled=.true.;patch_guard%radius=20d0
    patch_guard%filename=trim(root)//'_patch_failure.txt'
    call begin_projection_guard(patch_guard,2_8,np)
    patch_guard%recovered_mask=.true.
    before=fev%evaluation_calls
    call validate_recovered_patches(patch,fev,1_8,patch_guard,status)
    call require(status==6,'flipped recovered patch is rejected')
    call require(fev%evaluation_calls==before,'nodal orientation failure avoids field evaluation')

    ! r=(u,u*v,0) has positive Jacobians at the interior RV nodes but a
    ! collapsed edge u=0. The oversampling lattice must inspect that edge.
    do node=1,np
      patch%S_smooth(:,node)=[rv(1,node),rv(1,node)*rv(2,node),0d0]
      patch%du_smooth(:,node)=[1d0,rv(2,node),0d0]
      patch%dv_smooth(:,node)=[0d0,rv(1,node),0d0]
      patch%ds_smooth(node)=rv(1,node)
    enddo
    call require(minval(rv(1,:))>1d-3,'synthetic collapsed edge passes nodal Jacobian threshold')
    call begin_projection_guard(patch_guard,2_8,np)
    patch_guard%recovered_mask=.true.
    call validate_recovered_patches(patch,fev,1_8,patch_guard,status)
    call require(status==6,'collapsed boundary between RV nodes is rejected')
    call require(fev%evaluation_calls==before,'invalid oversampled geometry avoids field evaluation')

    ! Finite derivative components can still overflow their cross product;
    ! NaN/Inf comparisons must not let that nodal geometry reach the field.
    patch%S_smooth=0d0;patch%S_smooth(1:2,:)=rv
    patch%du_smooth=0d0;patch%du_smooth(1,:)=1d200
    patch%dv_smooth=0d0;patch%dv_smooth(2,:)=1d200
    patch%ds_smooth=1d0
    call begin_projection_guard(patch_guard,2_8,np)
    patch_guard%recovered_mask=.true.
    call validate_recovered_patches(patch,fev,1_8,patch_guard,status)
    call require(status==6,'finite derivative cross-product overflow is rejected')
    call require(fev%evaluation_calls==before,'overflowing nodal geometry avoids field evaluation')

    patch%du_smooth(1,1)=ieee_value(0d0,ieee_quiet_nan)
    call begin_projection_guard(patch_guard,2_8,np)
    patch_guard%recovered_mask=.true.
    call validate_recovered_patches(patch,fev,1_8,patch_guard,status)
    call require(status==6,'nonfinite implicit recovered derivatives are rejected')
    print *, 'PASS: flipped nodes, collapsed oversampled edges, and nonfinite derivatives are rejected.'
  end subroutine

  subroutine tangent_direction(gradient,direction)
    double precision,intent(in) :: gradient(3)
    double precision,intent(out) :: direction(3)
    double precision :: tangent(3)
    tangent=[-gradient(2),gradient(1),0d0]
    tangent=tangent/norm2(tangent)
    direction=(tangent-1d-6*gradient/norm2(gradient))/sqrt(1d0+1d-12)
  end subroutine

  function independent_root(base,direction,upper) result(value)
    double precision,intent(in) :: base(3),direction(3),upper
    double precision :: value,left,right,middle,fl,fr,fm
    integer :: iteration
    left=0d0;right=upper
    fl=direct_residual(base+left*direction)
    fr=direct_residual(base+right*direction)
    call require(fl>0d0.and.fr<0d0,'independent downward bracket exists')
    do iteration=1,40
      middle=(left+right)/2d0
      fm=direct_residual(base+middle*direction)
      if(fm>0d0) then
        left=middle
      else
        right=middle
      endif
    enddo
    value=(left+right)/2d0
  end function

  function direct_residual(point) result(value)
    ! An independent O(ns) evaluation avoids using the Brent implementation
    ! or FMM target-tree behavior to validate the accepted root.
    double precision,intent(in) :: point(3)
    double precision :: value,local_sigma(1),local_ds(3,1),target(3,1)
    double precision :: displacement(3),distance,kernel,d1,d2,pi
    integer*8 :: source
    target(:,1)=point
    call function_eval_sigma(fev%FSS_1,target,1_8,local_sigma,local_ds(1,:), &
        local_ds(2,:),local_ds(3,:),1_8)
    value=0d0;pi=4d0*atan(1d0)
    do source=1,g%n_Sk_points
      displacement=point-g%skeleton_Points(:,source)
      distance=norm2(displacement)
      call surf_smooth_ker(distance,local_sigma(1),kernel,d1,d2)
      value=value-dot_product(displacement,g%skeleton_N(:,source))*g%skeleton_w(source)*kernel/(4d0*pi)
    enddo
    value=value-0.5d0
  end function

  subroutine require(condition,label)
    logical,intent(in) :: condition
    character(*),intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      stop 1
    endif
  end subroutine
end program
