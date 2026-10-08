module Mod_Adaptive_Smooth_Surface
  use Mod_Smooth_Surface
  use iso_c_binding, only: c_char, c_null_char, c_int64_t, c_double, c_funptr, c_f_procpointer
  implicit none

  ! A leaf owns its launch chart, not a new level-set source.
  type adaptive_leaf
    integer*8 :: parent = 0, depth = 0
    real*8 :: map(2,3), points(3,6), normals(3,6), diameter
    real*8, allocatable :: values(:,:), height(:)
    real*8 :: indicator(4) = -1d0
    logical :: checked = .false.
  end type

  type adaptive_basis
    integer*8 :: order, np, nc
    real*8, allocatable :: uv(:,:), transform(:,:), pol(:,:), du(:,:), dv(:,:)
  end type

  abstract interface
    subroutine adaptive_progress(pass,npatches,unresolved,depth,indicators) bind(C)
      import c_int64_t,c_double
      integer(c_int64_t), intent(in) :: pass,npatches,unresolved,depth
      real(c_double), intent(in) :: indicators(4)
    end subroutine
  end interface

contains

  subroutine adaptive_child_map(child, a)
    integer*8, intent(in) :: child
    real*8, intent(out) :: a(2,3)
    a = 0
    a(1,2) = .5d0
    a(2,3) = .5d0
    select case(child)
    case(2)
      a(:,1) = .5d0
      a(1,2) = -.5d0
      a(2,3) = -.5d0
    case(3)
      a(1,1) = .5d0
    case(4)
      a(2,1) = .5d0
    end select
  end subroutine

  subroutine setup_adaptive_basis(p, b)
    integer*8, intent(in) :: p
    type(adaptive_basis), intent(out) :: b
    real*8, allocatable :: rv(:,:), v(:,:), w(:), x(:), gu(:,:), gv(:,:), gw(:), d(:,:)
    real*8 :: a(2,3), t
    integer*8 :: n, ng, j, k, offset
    n = (p+1)*(p+2)/2
    ng = p+2
    b%order = p
    b%np = n
    b%nc = 4*n+3*ng+3
    allocate(rv(2,n),v(n,n),w(n),b%transform(n,n))
    call vioreanu_simplex_quad(p,n,rv,b%transform,v,w)
    allocate(b%uv(2,b%nc),b%pol(n,b%nc),b%du(n,b%nc),b%dv(n,b%nc),d(2,n))
    do k=1,4
      call adaptive_child_map(k,a)
      do j=1,n
        b%uv(:,(k-1)*n+j) = a(:,1)+matmul(a(:,2:3),rv(:,j))
      enddo
    enddo
    allocate(x(ng),gu(ng,ng),gv(ng,ng),gw(ng))
    call legeexps(1_8,ng,x,gu,gv,gw)
    offset = 4*n
    do j=1,ng
      t = (1+x(j))/2
      b%uv(:,offset+j) = [t,0d0]
      b%uv(:,offset+ng+j) = [1-t,t]
      b%uv(:,offset+2*ng+j) = [0d0,1-t]
    enddo
    b%uv(:,b%nc-2) = [0d0,0d0]
    b%uv(:,b%nc-1) = [1d0,0d0]
    b%uv(:,b%nc) = [0d0,1d0]
    do j=1,b%nc
      call koorn_ders(b%uv(:,j),p,n,b%pol(:,j),d)
      b%du(:,j) = d(1,:)
      b%dv(:,j) = d(2,:)
    enddo
  end subroutine

  subroutine gather_adaptive_leaf(g, k, leaf)
    type(Geometry), intent(in) :: g
    integer*8, intent(in) :: k
    type(adaptive_leaf), intent(inout) :: leaf
    integer*8 :: i, first, last
    leaf%points = g%Points(:,g%Tri(:,k))
    leaf%normals = g%Normal_Vert(:,g%Tri(:,k))
    leaf%diameter = max(norm2(leaf%points(:,2)-leaf%points(:,1)), &
        norm2(leaf%points(:,3)-leaf%points(:,1)),norm2(leaf%points(:,3)-leaf%points(:,2)))
    first = (k-1)*g%n_order_sf+1
    last = k*g%n_order_sf
    if (allocated(leaf%values)) deallocate(leaf%values,leaf%height)
    allocate(leaf%values(12,g%n_order_sf),leaf%height(g%n_order_sf))
    leaf%values(1:3,:) = g%S_smooth(:,first:last)
    leaf%values(4:6,:) = g%du_smooth(:,first:last)
    leaf%values(7:9,:) = g%dv_smooth(:,first:last)
    leaf%values(10:12,:) = g%N_smooth(:,first:last)
    leaf%height = g%height(first:last)
    leaf%checked = .false.
    leaf%indicator = -1
  end subroutine

  subroutine pack_adaptive_parents(g, leaves, ids)
    type(Geometry), intent(inout) :: g
    type(adaptive_leaf), intent(in) :: leaves(:)
    integer*8, intent(in) :: ids(:)
    integer*8 :: j, k, n, np
    n = size(ids,kind=8)
    np = g%n_order_sf
    call clear_smooth_surface_work_arrays(g)
    deallocate(g%Points,g%Tri,g%Normal_Vert)
    g%ntri = n
    g%npoints = 6*n
    g%n_Sf_points = n*np
    allocate(g%Points(3,6*n),g%Normal_Vert(3,6*n),g%Tri(6,n),g%height(n*np))
    do j=1,n
      k = ids(j)
      g%Points(:,6*j-5:6*j) = leaves(k)%points
      g%Normal_Vert(:,6*j-5:6*j) = leaves(k)%normals
      g%Tri(:,j) = [(6*j-6+k,k=1,6)]
      g%height((j-1)*np+1:j*np) = leaves(ids(j))%height
    enddo
  end subroutine

  subroutine adaptive_coefficients(leaf,b,c,nc)
    type(adaptive_leaf), intent(in) :: leaf
    type(adaptive_basis), intent(in) :: b
    real*8, intent(out) :: c(3,b%np), nc(3,b%np)
    real*8 :: r(3,b%np), normals(3,b%np), center(3), length
    integer*8 :: j
    center = leaf%values(1:3,1)
    do j=1,b%np
      r(:,j) = leaf%values(1:3,j)-center
      length = norm2(leaf%values(10:12,j))
      normals(:,j) = leaf%values(10:12,j)/max(length,tiny(1d0))
    enddo
    c = matmul(r,transpose(b%transform))
    nc = matmul(normals,transpose(b%transform))
  end subroutine

  subroutine adaptive_tail(leaf,b)
    type(adaptive_leaf), intent(inout) :: leaf
    type(adaptive_basis), intent(in) :: b
    real*8 :: c(3,b%np), nc(3,b%np)
    integer*8 :: degree, first
    if (.not.all(ieee_is_finite(leaf%values)) .or. leaf%diameter<=0) then
      leaf%indicator = 1d100
      leaf%checked = .true.
      return
    endif
    call adaptive_coefficients(leaf,b,c,nc)
    degree = max(2_8,b%order-1)
    first = degree*(degree+1)/2+1
    leaf%indicator(1) = sqrt(2*sum(c(:,first:b%np)**2))/leaf%diameter
    degree = max(1_8,b%order-1)
    first = degree*(degree+1)/2+1
    leaf%indicator(2) = sqrt(2*sum(nc(:,first:b%np)**2))
  end subroutine

  subroutine assess_adaptive_leaves(leaves,b,g,feval,guard,mode,tol)
    type(adaptive_leaf), intent(inout) :: leaves(:)
    type(adaptive_basis), intent(in) :: b
    type(Geometry), intent(inout) :: g
    type(Feval_stuff), intent(inout) :: feval
    type(projection_radius_guard), intent(in) :: guard
    integer*8, intent(in) :: mode
    real*8, intent(in) :: tol
    integer*8, allocatable :: owner(:), queue(:)
    real*8, allocatable :: targets(:,:), normals(:,:), pnormals(:,:), f(:), grad(:,:)
    real*8 :: c(3,b%np), nc(3,b%np), r(3,b%nc), du(3,b%nc), dv(3,b%nc), nn(3,b%nc)
    real*8 :: crossn(3), lengthn, rad, jac
    integer*8 :: i,j,k,nq,start,last,nt,capacity
    integer :: reason
    allocate(queue(size(leaves)))
    nq = 0
    do i=1,size(leaves,kind=8)
      if (leaves(i)%checked) cycle
      call adaptive_tail(leaves(i),b)
      leaves(i)%checked = .true.
      if (maxval(leaves(i)%indicator(1:2))>tol) cycle
      nq = nq+1
      queue(nq) = i
      leaves(i)%indicator(3:4) = 0
    enddo
    ! Bound independent-check workspace; each FMM call still batches many patches.
    if (nq==0) return
    capacity = min(nq,max(1_8,100000_8/b%nc))
    allocate(targets(3,capacity*b%nc),normals(3,capacity*b%nc), &
        pnormals(3,capacity*b%nc),f(capacity*b%nc),grad(3,capacity*b%nc),owner(capacity*b%nc))
    do start=1,nq,capacity
      last = min(nq,start+capacity-1)
      nt = 0
      do k=start,last
        i = queue(k)
        call adaptive_coefficients(leaves(i),b,c,nc)
        r = matmul(c,b%pol)
        du = matmul(c,b%du)
        dv = matmul(c,b%dv)
        nn = matmul(nc,b%pol)
        do j=1,b%nc
          r(:,j) = r(:,j)+leaves(i)%values(1:3,1)
          call projection_target_status(guard,r(:,j),rad,reason)
          crossn = cross_product_adaptive(du(:,j)/guard%radius,dv(:,j)/guard%radius)
          jac = norm2(crossn)
          lengthn = norm2(nn(:,j))
          if (reason/=0 .or. .not.ieee_is_finite(jac) .or. &
              .not.ieee_is_finite(lengthn) .or. jac<=1d-28 .or. lengthn<=1d-14) then
            leaves(i)%indicator(3:4) = 1d100
            cycle
          endif
          nt = nt+1
          targets(:,nt) = r(:,j)
          normals(:,nt) = nn(:,j)/lengthn
          pnormals(:,nt) = crossn/jac
          owner(nt) = i
        enddo
      enddo
      if (nt==0) cycle
      call eval_density_grad_FMM(g,targets(:,1:nt),normals(:,1:nt),nt,f(1:nt),grad(:,1:nt),feval,mode)
      do j=1,nt
        i = owner(j)
        call update_adaptive_checks(leaves(i),guard%radius,f(j),grad(:,j),normals(:,j),pnormals(:,j))
      enddo
    enddo
  end subroutine

  subroutine update_adaptive_checks(leaf,radius,residual,gradient,normal,polynomial_normal)
    type(adaptive_leaf), intent(inout) :: leaf
    real*8, intent(in) :: radius,residual,gradient(3),normal(3),polynomial_normal(3)
    real*8 :: length,gn(3),errorn
    length=norm2(gradient)
    if (.not.ieee_is_finite(residual) .or. .not.ieee_is_finite(length) .or. length*radius<=1d-14) then
      leaf%indicator(3:4)=1d100
      return
    endif
    gn=-gradient/length
    errorn=max(norm2(gn-normal),norm2(gn-polynomial_normal))
    if (dot_product(gn,polynomial_normal)<=0 .or. dot_product(gn,normal)<=0) errorn=1d100
    leaf%indicator(3)=max(leaf%indicator(3),abs(residual)/length/leaf%diameter)
    leaf%indicator(4)=max(leaf%indicator(4),errorn)
  end subroutine

  function cross_product_adaptive(a,b) result(c)
    real*8, intent(in) :: a(3),b(3)
    real*8 :: c(3)
    c = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  subroutine merge_adaptive_children(leaves,g,ids)
    type(adaptive_leaf), allocatable, intent(inout) :: leaves(:)
    type(Geometry), intent(in) :: g
    integer*8, intent(in) :: ids(:)
    type(adaptive_leaf), allocatable :: next(:)
    integer*8 :: i,j,k,q,m,n
    real*8 :: a(2,3)
    n = size(leaves,kind=8)
    allocate(next(n+3*size(ids)))
    m=1
    q=0
    do i=1,n
      if (m<=size(ids)) then
        if (ids(m)==i) then
          do j=1,4
            q=q+1
            k=4*(m-1)+j
            call gather_adaptive_leaf(g,k,next(q))
            next(q)%parent=leaves(i)%parent
            next(q)%depth=leaves(i)%depth+1
            call adaptive_child_map(j,a)
            next(q)%map(:,1)=leaves(i)%map(:,1)+matmul(leaves(i)%map(:,2:3),a(:,1))
            next(q)%map(:,2:3)=matmul(leaves(i)%map(:,2:3),a(:,2:3))
          enddo
          m=m+1
          cycle
        endif
      endif
      q=q+1
      next(q)=leaves(i)
    enddo
    call move_alloc(next,leaves)
  end subroutine

  subroutine run_adaptive_smoother(fname,cad,outroot,filetype,nquad,p,mode,rlam,tol,maxdepth,maxpoints,ier,progress)
    character(*), intent(in) :: fname,cad,outroot
    integer*8, intent(in) :: filetype,nquad,p,mode,maxdepth,maxpoints
    real*8, intent(in) :: rlam,tol
    integer*8, intent(out) :: ier
    procedure(adaptive_progress), optional :: progress
    call run_adaptive_smoother_recovery(fname,cad,outroot,filetype,nquad,p,mode,rlam,tol, &
        maxdepth,maxpoints,ier,progress,.true.)
  end subroutine

  subroutine run_adaptive_smoother_recovery(fname,cad,outroot,filetype,nquad,p,mode,rlam,tol, &
      maxdepth,maxpoints,ier,progress,newton_recovery)
    character(*), intent(in) :: fname,cad,outroot
    integer*8, intent(in) :: filetype,nquad,p,mode,maxdepth,maxpoints
    real*8, intent(in) :: rlam,tol
    integer*8, intent(out) :: ier
    procedure(adaptive_progress), optional :: progress
    logical, intent(in) :: newton_recovery
    type(Geometry) :: g
    type(Feval_stuff) :: feval
    type(projection_radius_guard) :: guard
    type(adaptive_leaf), allocatable :: leaves(:)
    type(adaptive_basis) :: basis
    integer*8, allocatable :: ids(:)
    real*8 :: history(11,maxdepth+1), maxima(4)
    integer*8 :: np,pass,i,n,marked,unresolved,deepest,stopcode
    ier=0
    call readgeometry(g,trim(fname),filetype,nquad,p,ier)
    if (ier/=0) return
    np=(p+1)*(p+2)/2
    if (g%ntri>maxpoints/np) then
      ier=8
      return
    endif
    call load_cad_skeleton(g,cad)
    call funcion_normal_vert(g)
    call initialize_projection_guard(guard,g,outroot,ier,newton_recovery)
    if (ier/=0) return
    call start_Feval_tree(feval,g,rlam,mode)
    call project_scaffold_vertices_to_levelset(g,feval,mode,ier,guard)
    if (ier/=0) goto 900
    call funcion_normal_vert(g)
    call funcion_Base_Points(g)
    call find_smooth_surface(g,feval,mode,ier,guard)
    if (ier/=0) goto 900
    allocate(leaves(g%ntri))
    do i=1,g%ntri
      call gather_adaptive_leaf(g,i,leaves(i))
      leaves(i)%parent=i
      leaves(i)%map=0
      leaves(i)%map(1,2)=1
      leaves(i)%map(2,3)=1
    enddo
    call setup_adaptive_basis(p,basis)
    history=0
    do pass=0,maxdepth
      call assess_adaptive_leaves(leaves,basis,g,feval,guard,mode,tol)
      n=size(leaves,kind=8)
      allocate(ids(n))
      marked=0
      unresolved=0
      deepest=0
      maxima=-1
      do i=1,n
        maxima=max(maxima,leaves(i)%indicator)
        deepest=max(deepest,leaves(i)%depth)
        if (maxval(leaves(i)%indicator)<=tol) cycle
        unresolved=unresolved+1
        if (leaves(i)%depth>=maxdepth) cycle
        marked=marked+1
        ids(marked)=i
      enddo
      history(:,pass+1)=[real(pass,8),real(n,8),real(n*np,8),real(marked,8), &
          real(unresolved,8),real(deepest,8),maxima,real(count([(leaves(i)%indicator(3)>=0,i=1,n)]),8)]
      print *, 'Adaptive pass, patches, unresolved, depth:',pass,n,unresolved,deepest
      print *, '  maximum indicators:',maxima
      if (present(progress)) call progress(pass,n,unresolved,deepest,maxima)
      stopcode=0
      if (unresolved==0) exit
      stopcode=1
      if (marked==0) exit
      stopcode=2
      if (marked>(maxpoints/np-n)/3) exit
      call pack_adaptive_parents(g,leaves,ids(:marked))
      call refine_geometry_smart(g)
      call funcion_Base_Points(g)
      guard%refinement=pass+1
      call find_smooth_surface(g,feval,mode,ier,guard)
      if (ier/=0) goto 900
      call merge_adaptive_children(leaves,g,ids(:marked))
      deallocate(ids)
    enddo
    call write_adaptive_result(outroot,p,leaves,history(:,:pass+1),stopcode,guard,ier)
900 continue
    call destroy_Feval_tree(feval)
  end subroutine

  subroutine write_adaptive_result(root,p,leaves,history,stopcode,guard,ier)
    character(*), intent(in) :: root
    integer*8, intent(in) :: p,stopcode
    type(adaptive_leaf), intent(in) :: leaves(:)
    real*8, intent(in) :: history(:,:)
    type(projection_radius_guard), intent(in) :: guard
    integer*8, intent(out) :: ier
    integer :: unit,ios
    integer*8 :: i,j,k,n
    ier=9
    n=size(leaves,kind=8)
    open(newunit=unit,file=root//'.go3',status='replace',iostat=ios)
    if (ios/=0) return
    write(unit,*,iostat=ios) p
    if (ios==0) write(unit,*,iostat=ios) n
    ! GO3 stores one component at a time, across all patches.
    do j=1,12
      do i=1,n
        if (ios==0) write(unit,'(es26.17e3)',iostat=ios) leaves(i)%values(j,:)
      enddo
    enddo
    close(unit)
    if (ios/=0) return
    open(newunit=unit,file=root//'_leaves.txt',status='replace',iostat=ios)
    if (ios/=0) return
    do i=1,n
      if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) real(leaves(i)%parent,8), &
          real(leaves(i)%depth,8),leaves(i)%map,leaves(i)%indicator,leaves(i)%diameter
    enddo
    close(unit)
    if (ios/=0) return
    open(newunit=unit,file=root//'_history.txt',status='replace',iostat=ios)
    if (ios/=0) return
    do i=1,size(history,2)
      if (ios==0) write(unit,'(*(es26.17e3,1x))',iostat=ios) history(:,i)
    enddo
    close(unit)
    if (ios/=0) return
    open(newunit=unit,file=root//'_status.txt',status='replace',iostat=ios)
    if (ios/=0) return
    write(unit,*,iostat=ios) stopcode,guard%center,guard%radius
    close(unit)
    if (ios==0) ier=0
  end subroutine

  subroutine adaptive_c(fname,cad,outroot,filetype,nquad,p,mode,rlam,tol,maxdepth,maxpoints,ier,progress_function) &
      bind(C,name='multiscale_mesher_adaptive_c')
    character(c_char), intent(in) :: fname(*),cad(*),outroot(*)
    integer(c_int64_t), intent(in) :: filetype,nquad,p,mode,maxdepth,maxpoints
    real(c_double), intent(in) :: rlam,tol
    integer(c_int64_t), intent(out) :: ier
    type(c_funptr), value :: progress_function
    procedure(adaptive_progress), pointer :: progress
    call c_f_procpointer(progress_function,progress)
    call run_adaptive_smoother(c_string(fname),c_string(cad),c_string(outroot), &
        filetype,nquad,p,mode,rlam,tol,maxdepth,maxpoints,ier,progress)
  end subroutine

  subroutine adaptive_recovery_c(fname,cad,outroot,filetype,nquad,p,mode,rlam,tol,maxdepth,maxpoints, &
      newton_recovery,ier,progress_function) &
      bind(C,name='multiscale_mesher_adaptive_recovery_c')
    character(c_char), intent(in) :: fname(*),cad(*),outroot(*)
    integer(c_int64_t), intent(in) :: filetype,nquad,p,mode,maxdepth,maxpoints,newton_recovery
    real(c_double), intent(in) :: rlam,tol
    integer(c_int64_t), intent(out) :: ier
    type(c_funptr), value :: progress_function
    procedure(adaptive_progress), pointer :: progress
    call c_f_procpointer(progress_function,progress)
    call run_adaptive_smoother_recovery(c_string(fname),c_string(cad),c_string(outroot), &
        filetype,nquad,p,mode,rlam,tol,maxdepth,maxpoints,ier,progress,newton_recovery/=0)
  end subroutine

  function c_string(chars) result(str)
    character(c_char), intent(in) :: chars(*)
    character(:), allocatable :: str
    integer :: n,j
    n=0
    do while(chars(n+1)/=c_null_char)
      n=n+1
    enddo
    allocate(character(n)::str)
    do j=1,n
      str(j:j)=chars(j)
    enddo
  end function
end module Mod_Adaptive_Smooth_Surface
