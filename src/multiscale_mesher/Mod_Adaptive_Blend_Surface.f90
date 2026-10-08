module Mod_Adaptive_Blend_Surface
  use Mod_Adaptive_Smooth_Surface
  implicit none
  private
  public :: blend_open, blend_open_recovery, blend_close, blend_close_all, blend_size, &
      blend_get, blend_refine, blend_project, blend_sigma

  type blend_context
    type(Geometry) :: g
    type(Feval_stuff) :: feval
    type(projection_radius_guard) :: guard
    type(adaptive_leaf), allocatable :: leaves(:)
    type(adaptive_basis) :: basis
    integer*8 :: id=0, mode=1
    logical :: tree_ready=.false.
  end type
  type context_slot
    type(blend_context), pointer :: ctx=>null()
  end type
  ! Only ownership lives in this registry; all numerical state belongs to a run.
  type(context_slot), save :: slots(8)
  integer*8, save :: next_id=1
contains
  function lookup(id) result(c)
    integer*8, intent(in) :: id
    type(blend_context), pointer :: c
    integer :: j
    nullify(c)
    do j=1,size(slots)
      if (.not.associated(slots(j)%ctx)) cycle
      if (slots(j)%ctx%id==id) then
        c=>slots(j)%ctx
        return
      endif
    enddo
  end function

  subroutine blend_close(id,ier) bind(C,name='adaptive_blend_close')
    integer(c_int64_t), intent(in) :: id
    integer(c_int64_t), intent(out) :: ier
    integer :: j
    ier=10
    do j=1,size(slots)
      if (.not.associated(slots(j)%ctx)) cycle
      if (slots(j)%ctx%id/=id) cycle
      if (slots(j)%ctx%tree_ready) call destroy_Feval_tree(slots(j)%ctx%feval)
      deallocate(slots(j)%ctx)
      nullify(slots(j)%ctx)
      ier=0
      return
    enddo
  end subroutine

  subroutine blend_close_all() bind(C,name='adaptive_blend_close_all')
    integer*8 :: id,ier
    integer :: j
    do j=1,size(slots)
      if (.not.associated(slots(j)%ctx)) cycle
      id=slots(j)%ctx%id
      call blend_close(id,ier)
    enddo
  end subroutine

  subroutine blend_open(fname,cad,root,nquad,p,mode,rlam,budget,id,ier) bind(C,name='adaptive_blend_open')
    character(c_char), intent(in) :: fname(*),cad(*),root(*)
    integer(c_int64_t), intent(in) :: nquad,p,mode,budget
    real(c_double), intent(in) :: rlam
    integer(c_int64_t), intent(out) :: id,ier
    call blend_open_recovery(fname,cad,root,nquad,p,mode,rlam,budget,1_c_int64_t,id,ier)
  end subroutine

  subroutine blend_open_recovery(fname,cad,root,nquad,p,mode,rlam,budget,newton_recovery,id,ier) &
      bind(C,name='adaptive_blend_open_recovery')
    character(c_char), intent(in) :: fname(*),cad(*),root(*)
    integer(c_int64_t), intent(in) :: nquad,p,mode,budget,newton_recovery
    real(c_double), intent(in) :: rlam
    integer(c_int64_t), intent(out) :: id,ier
    type(blend_context), pointer :: c
    integer*8 :: j,close_error
    id=0
    ier=11
    nullify(c)
    do j=1,size(slots)
      if (associated(slots(j)%ctx)) cycle
      allocate(slots(j)%ctx)
      c=>slots(j)%ctx
      exit
    enddo
    if (.not.associated(c)) return
    id=next_id
    next_id=next_id+1
    c%id=id
    c%mode=mode
    call readgeometry(c%g,c_string(fname),3_8,nquad,p,ier)
    if (ier/=0) goto 900
    if (c%g%ntri>budget/((p+1)*(p+2)/2)) then
      ier=8
      goto 900
    endif
    call load_cad_skeleton(c%g,c_string(cad))
    call funcion_normal_vert(c%g)
    call initialize_projection_guard(c%guard,c%g,c_string(root),ier,newton_recovery/=0)
    if (ier/=0) goto 900
    call start_Feval_tree(c%feval,c%g,rlam,mode)
    c%tree_ready=.true.
    call project_scaffold_vertices_to_levelset(c%g,c%feval,mode,ier,c%guard)
    if (ier/=0) goto 900
    call funcion_normal_vert(c%g)
    call funcion_Base_Points(c%g)
    call find_smooth_surface(c%g,c%feval,mode,ier,c%guard)
    if (ier/=0) goto 900
    allocate(c%leaves(c%g%ntri))
    do j=1,c%g%ntri
      call gather_adaptive_leaf(c%g,j,c%leaves(j))
      c%leaves(j)%parent=j
      c%leaves(j)%map=0
      c%leaves(j)%map(1,2)=1
      c%leaves(j)%map(2,3)=1
    enddo
    call setup_adaptive_basis(p,c%basis)
    return
900 call blend_close(id,close_error)
    id=0
  end subroutine

  subroutine blend_size(id,n,np,sphere,ier) bind(C,name='adaptive_blend_size')
    integer(c_int64_t), intent(in) :: id
    integer(c_int64_t), intent(out) :: n,np,ier
    real(c_double), intent(out) :: sphere(4)
    type(blend_context), pointer :: c
    c=>lookup(id)
    ier=10
    n=0
    np=0
    if (.not.associated(c)) return
    n=size(c%leaves,kind=8)
    np=c%basis%np
    sphere=[c%guard%center,c%guard%radius]
    ier=0
  end subroutine

  subroutine blend_get(id,n,np,values,meta) bind(C,name='adaptive_blend_get')
    integer(c_int64_t), intent(in) :: id,n,np
    real(c_double), intent(out) :: values(12,np,n),meta(9,n)
    type(blend_context), pointer :: c
    integer*8 :: j
    c=>lookup(id)
    do j=1,n
      values(:,:,j)=c%leaves(j)%values
      meta(1:2,j)=[real(c%leaves(j)%parent,8),real(c%leaves(j)%depth,8)]
      meta(3:8,j)=reshape(c%leaves(j)%map,[6])
      meta(9,j)=c%leaves(j)%diameter
    enddo
  end subroutine

  subroutine blend_refine(id,n,ids,budget,ier) bind(C,name='adaptive_blend_refine')
    integer(c_int64_t), intent(in) :: id,n,ids(n),budget
    integer(c_int64_t), intent(out) :: ier
    type(blend_context), pointer :: c
    integer*8 :: j
    c=>lookup(id)
    ier=10
    if (.not.associated(c)) return
    ier=12
    if (n<1 .or. any(ids<1) .or. any(ids>size(c%leaves))) return
    do j=2,n
      if (ids(j)<=ids(j-1)) return
    enddo
    ier=8
    if (n>(budget/c%basis%np-size(c%leaves,kind=8))/3) return
    call pack_adaptive_parents(c%g,c%leaves,ids)
    call refine_geometry_smart(c%g)
    call funcion_Base_Points(c%g)
    c%guard%refinement=maxval(c%leaves(ids)%depth)+1
    call find_smooth_surface(c%g,c%feval,c%mode,ier,c%guard)
    if (ier/=0) return
    call merge_adaptive_children(c%leaves,c%g,ids)
  end subroutine

  subroutine launch_at(leaf,uv,b,bu,bv,n,nu,nv)
    type(adaptive_leaf), intent(in) :: leaf
    real*8, intent(in) :: uv(2)
    real*8, intent(out) :: b(3),bu(3),bv(3),n(3),nu(3),nv(3)
    real*8 :: u,v,w,l(6),lu(6),lv(6)
    u=uv(1)
    v=uv(2)
    w=1-u-v
    l=[w*(2*w-1),u*(2*u-1),v*(2*v-1),4*u*w,4*u*v,4*v*w]
    lu=[1-4*w,4*u-1,0d0,4*(w-u),4*v,-4*v]
    lv=[1-4*w,0d0,4*v-1,-4*u,4*u,4*(w-v)]
    b=matmul(leaf%points,l)
    bu=matmul(leaf%points,lu)
    bv=matmul(leaf%points,lv)
    n=w*leaf%normals(:,1)+u*leaf%normals(:,2)+v*leaf%normals(:,3)
    nu=leaf%normals(:,2)-leaf%normals(:,1)
    nv=leaf%normals(:,3)-leaf%normals(:,1)
  end subroutine

  subroutine blend_project(id,nt,ids,uv,values,ier) bind(C,name='adaptive_blend_project')
    integer(c_int64_t), intent(in) :: id,nt,ids(nt)
    real(c_double), intent(in) :: uv(2,nt)
    real(c_double), intent(out) :: values(12,nt)
    integer(c_int64_t), intent(out) :: ier
    type(blend_context), pointer :: c
    real*8, allocatable :: saved_b(:,:),saved_n(:,:),bu(:,:),bv(:,:),nu(:,:),nv(:,:)
    real*8, allocatable :: h(:),grad(:,:),r(:,:),pol(:),coeff(:),f(:)
    real*8 :: denominator,hu,hv,cross(3),jacobian
    logical :: recovered
    integer*8 :: j,old_nt,flag
    c=>lookup(id)
    ier=10
    if (.not.associated(c)) return
    ier=12
    if (nt<1 .or. any(ids<1) .or. any(ids>size(c%leaves))) return
    if (any(.not.ieee_is_finite(uv)) .or. any(uv< -1d-12)) return
    if (any(sum(uv,dim=1)>1+1d-12)) return
    values=ieee_value(0d0,ieee_quiet_nan)
    allocate(bu(3,nt),bv(3,nt),nu(3,nt),nv(3,nt),h(nt),grad(3,nt),r(3,nt),f(nt))
    allocate(pol(c%basis%np),coeff(c%basis%np))
    ! My_Newton already accepts arbitrary targets through these launch arrays.
    old_nt=c%g%n_Sf_points
    call move_alloc(c%g%Base_Points,saved_b)
    call move_alloc(c%g%Base_Points_N,saved_n)
    allocate(c%g%Base_Points(3,nt),c%g%Base_Points_N(3,nt))
    c%g%n_Sf_points=nt
    do j=1,nt
      call launch_at(c%leaves(ids(j)),uv(:,j),c%g%Base_Points(:,j),bu(:,j),bv(:,j), &
          c%g%Base_Points_N(:,j),nu(:,j),nv(:,j))
      call koorn_pols(uv(:,j),c%basis%order,c%basis%np,pol)
      coeff=matmul(c%basis%transform,c%leaves(ids(j))%height)
      h(j)=dot_product(coeff,pol)
    enddo
    c%guard%refinement=maxval(c%leaves(ids)%depth)
    call My_Newton(h,1d-9,14_8,c%g,flag,c%feval,c%mode,grad,r,c%guard)
    ier=0
    if (flag/=0) then
      ier=4
      if (flag==1) ier=3
      if (flag==7) ier=7
      goto 900
    endif
    do j=1,nt
      r(:,j)=c%g%Base_Points(:,j)+h(j)*c%g%Base_Points_N(:,j)
    enddo
    call eval_density_grad_FMM(c%g,r,c%g%Base_Points_N,nt,f,grad,c%feval,c%mode)
    do j=1,nt
      recovered=.false.
      if (allocated(c%guard%recovered_mask)) then
        if (size(c%guard%recovered_mask)==nt) recovered=c%guard%recovered_mask(j)
      endif
      denominator=dot_product(grad(:,j),c%g%Base_Points_N(:,j))
      if (.not.ieee_is_finite(denominator)) then
        if (recovered) then
          call reject_recovered_geometry(c%guard,'nonfinite directional derivative at recovered blend target',ier)
          goto 900
        endif
        cycle
      endif
      if (abs(denominator)<=1d-12*norm2(grad(:,j))*norm2(c%g%Base_Points_N(:,j))) then
        if (recovered) then
          call reject_recovered_geometry(c%guard,'singular directional derivative at recovered blend target',ier)
          goto 900
        endif
        cycle
      endif
      if (norm2(grad(:,j))*c%guard%radius<=1d-14) then
        if (recovered) then
          call reject_recovered_geometry(c%guard,'vanishing field gradient at recovered blend target',ier)
          goto 900
        endif
        cycle
      endif
      if (recovered) then
        if (denominator>=-1d-12*norm2(grad(:,j))*norm2(c%g%Base_Points_N(:,j))) then
          call reject_recovered_geometry(c%guard,'outward crossing lost at recovered blend target',ier)
          goto 900
        endif
      endif
      hu=-dot_product(grad(:,j),bu(:,j)+h(j)*nu(:,j))/denominator
      hv=-dot_product(grad(:,j),bv(:,j)+h(j)*nv(:,j))/denominator
      values(1:3,j)=r(:,j)
      values(4:6,j)=bu(:,j)+h(j)*nu(:,j)+hu*c%g%Base_Points_N(:,j)
      values(7:9,j)=bv(:,j)+h(j)*nv(:,j)+hv*c%g%Base_Points_N(:,j)
      values(10:12,j)=-grad(:,j)/norm2(grad(:,j))
      if (recovered) then
        cross(1)=values(5,j)*values(9,j)-values(6,j)*values(8,j)
        cross(2)=values(6,j)*values(7,j)-values(4,j)*values(9,j)
        cross(3)=values(4,j)*values(8,j)-values(5,j)*values(7,j)
        jacobian=norm2(cross)
        if (any(.not.ieee_is_finite(values(:,j))) .or. .not.ieee_is_finite(jacobian)) then
          call reject_recovered_geometry(c%guard,'nonfinite implicit derivatives at recovered blend target',ier)
          goto 900
        endif
        if (jacobian<=0d0 .or. dot_product(cross,values(10:12,j))<=0d0) then
          call reject_recovered_geometry(c%guard,'invalid implicit Jacobian at recovered blend target',ier)
          goto 900
        endif
      endif
    enddo
900 deallocate(c%g%Base_Points,c%g%Base_Points_N)
    call move_alloc(saved_b,c%g%Base_Points)
    call move_alloc(saved_n,c%g%Base_Points_N)
    c%g%n_Sf_points=old_nt
  end subroutine

  subroutine blend_sigma(id,nt,r,values,ier) bind(C,name='adaptive_blend_sigma')
    integer(c_int64_t), intent(in) :: id,nt
    real(c_double), intent(in) :: r(3,nt)
    real(c_double), intent(out) :: values(4,nt)
    integer(c_int64_t), intent(out) :: ier
    type(blend_context), pointer :: c
    c=>lookup(id)
    ier=10
    if (.not.associated(c)) return
    ier=12
    if (any(.not.ieee_is_finite(r))) return
    call function_eval_sigma(c%feval%FSS_1,r,nt,values(1,:),values(2,:),values(3,:),values(4,:),c%mode)
    ier=0
  end subroutine
end module Mod_Adaptive_Blend_Surface
