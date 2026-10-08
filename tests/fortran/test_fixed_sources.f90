program test_fixed_sources
  use Mod_Adaptive_Smooth_Surface
  implicit none
  integer*8, parameter :: nt=32
  type(Geometry) :: g
  type(Feval_stuff) :: fev
  type(projection_radius_guard) :: guard
  type(adaptive_leaf), allocatable :: leaves(:)
  integer*8 :: ier, i, j, nsrc, nseed, nleaf, mode, ids(2)
  character(len=1024) :: meshfile, cadfile, root, argument
  double precision :: rlam, alpha0, mean0
  double precision :: targets(3,nt), normals(3,nt), f0(nt), f1(nt)
  double precision :: grad0(3,nt), grad1(3,nt)
  double precision :: s0(nt), s1(nt), ds0(3,nt), ds1(3,nt)
  double precision, allocatable :: cadp(:,:), cadn(:,:), cadw(:), launch0(:,:)
  double precision, allocatable :: seedp(:,:), seeds(:), unchanged(:,:)
  double precision, allocatable :: first(:), different(:), repeated(:)

  ! Usage: test_fixed_sources mesh.gidmsh cad.txt output_root [mode] [rlam]
  ! Link with the native library and matching surface objects/LP64 BLAS, or
  ! with the MATLAB archive and its matching ILP64 BLAS. Do not mix ABIs.
  call get_command_argument(1,meshfile)
  call get_command_argument(2,cadfile)
  call get_command_argument(3,root)
  call require(len_trim(meshfile)>0.and.len_trim(cadfile)>0.and.len_trim(root)>0,'three paths required')
  mode=1
  rlam=5.0d0
  call get_command_argument(4,argument)
  if (len_trim(argument)>0) read(argument,*) mode
  call get_command_argument(5,argument)
  if (len_trim(argument)>0) read(argument,*) rlam
  call readgeometry(g,trim(meshfile),3_8,16_8,4_8,ier)
  call require(ier==0,'read scaffold')
  call load_cad_skeleton(g,trim(cadfile))
  call funcion_normal_vert(g)
  call initialize_projection_guard(guard,g,trim(root)//'_invariance',ier)
  call require(ier==0,'initialize radius guard')
  call start_Feval_tree(fev,g,rlam,mode)
  nsrc=g%n_Sk_points
  nseed=g%ntri_sk
  cadp=g%skeleton_Points
  cadn=g%skeleton_N
  cadw=g%skeleton_w
  launch0=g%Points
  seedp=fev%FSS_1%TreeLRD_1%W_Pts_mem
  seeds=fev%FSS_1%TreeLRD_1%W_sgmas_mem
  alpha0=fev%FSS_1%alpha
  mean0=fev%FSS_1%sgma_mean
  do i=1,nt
    j=1+(i-1)*(nsrc-1)/(nt-1)
    targets(:,i)=cadp(:,j)+0.03d0*cadn(:,j)
  enddo
  normals=0.0d0
  normals(3,:)=1.0d0
  call eval_density_grad_FMM(g,targets,normals,nt,f0,grad0,fev,mode)
  call function_eval_sigma(fev%FSS_1,targets,nt,s0,ds0(1,:),ds0(2,:),ds0(3,:),mode)
  call require(all(ieee_is_finite(f0)).and.all(ieee_is_finite(grad0)),'finite reference field')
  call require(all(ieee_is_finite(s0)).and.all(ieee_is_finite(ds0)),'finite reference sigma')

  ! The included cylinder admits an almost tangential projection line with
  ! an unsafe Newton proposal and a valid nearby downward crossing. Verify
  ! the fixed-source invariants across that real recovery event as well.
  if (mode==1.and.abs(rlam-5d0)<1d-14.and.index(trim(meshfile),'cylinder')>0) then
    call verify_recovery_event()
    call verify('after deferred Newton recovery')
  endif

  call project_scaffold_vertices_to_levelset(g,fev,mode,ier,guard)
  call require(ier==0,'stage-one projection')
  print *, 'Max scaffold coordinate change:',maxval(abs(g%Points-launch0))
  call require(maxval(abs(g%Points-launch0))>1.0d-10,'nontrivial scaffold movement')
  call verify('after scaffold projection')
  call funcion_normal_vert(g)
  call funcion_Base_Points(g)
  call verify('after rebuilding normals and base points')
  call find_smooth_surface(g,fev,mode,ier,guard)
  call require(ier==0,'stage-two projection')
  call verify('after full smooth-surface projection')

  call refine_geometry_smart(g)
  call funcion_Base_Points(g)
  guard%refinement=1
  call find_smooth_surface(g,fev,mode,ier,guard)
  call require(ier==0,'uniform refinement projection')
  call verify('after uniform refinement and projection')

  nleaf=g%ntri
  call require(nleaf>2,'enough leaves for a selective refinement')
  allocate(leaves(nleaf))
  do i=1,nleaf
    call gather_adaptive_leaf(g,i,leaves(i))
    leaves(i)%parent=i
    leaves(i)%map=0.0d0
    leaves(i)%map(1,2)=1.0d0
    leaves(i)%map(2,3)=1.0d0
  enddo
  unchanged=leaves(2)%values
  ids=[1_8,nleaf]
  call pack_adaptive_parents(g,leaves,ids)
  call verify('after packing selected parents')
  call refine_geometry_smart(g)
  call funcion_Base_Points(g)
  guard%refinement=2
  call find_smooth_surface(g,fev,mode,ier,guard)
  call require(ier==0,'selective refinement projection')
  call merge_adaptive_children(leaves,g,ids)
  call require(size(leaves,kind=8)==nleaf+6,'only selected leaves split')
  call require(all(leaves(5)%values==unchanged),'unselected leaf is unchanged')
  call verify('after selective adaptive refinement and projection')
  call destroy_Feval_tree(fev)
  call require(.not.associated(fev%FSS_1),'tree owner cleared')
  call destroy_Feval_tree(fev)
  print *, 'FIXED_PHI_ADAPTIVE_INVARIANCE_PASS mode=',mode

  ! Exercise the real legacy public entry point A, B, A in one process.
  ! Changing rlam must change B but must not contaminate the second A.
  call legacy_run('a',5.0d0,first)
  call legacy_run('b',4.0d0,different)
  call legacy_run('a_repeat',5.0d0,repeated)
  call require(size(first)==size(repeated).and.size(first)==size(different),'legacy output dimensions')
  call require(maxval(abs(first-repeated))<1.0d-12,'repeated legacy result is independent of previous call')
  call require(maxval(abs(first-different))>1.0d-10,'legacy changed parameter takes effect')
  print *, 'LEGACY_REPEATED_CALL_PASS max delta A/A=',maxval(abs(first-repeated))

contains
  subroutine verify_recovery_event()
    double precision :: base(3,1),normal(3,1),projected(3,1),height(1),gradient(3,1),residual(1)
    double precision :: tangent(3),unit_gradient(3)
    integer*8 :: before,local_status
    base(:,1)=[9.9d0,0d0,0d0]
    normal(:,1)=[1d0,0d0,0d0]
    call eval_density_grad_FMM(g,base,normal,1_8,residual,gradient,fev,mode)
    unit_gradient=gradient(:,1)/norm2(gradient(:,1))
    tangent=[-gradient(2,1),gradient(1,1),0d0]
    tangent=tangent/norm2(tangent)
    normal(:,1)=(tangent-1d-6*unit_gradient)/sqrt(1d0+1d-12)
    before=guard%recovery_recovered
    call project_points_to_levelset(g,fev,mode,1_8,base,normal,projected,height,gradient,local_status,guard)
    call require(local_status==0,'fixed-source recovery succeeds')
    call require(guard%recovery_recovered==before+1,'fixed-source test exercised actual recovery')
  end subroutine verify_recovery_event

  subroutine verify(label)
    character(len=*), intent(in) :: label
    call require(g%n_Sk_points==nsrc.and.g%ntri_sk==nseed,'source counts unchanged')
    call require(all(g%skeleton_Points==cadp),'CAD positions unchanged')
    call require(all(g%skeleton_N==cadn),'CAD normals unchanged')
    call require(all(g%skeleton_w==cadw),'CAD weights unchanged')
    call require(all(fev%FSS_1%TreeLRD_1%W_Pts_mem==seedp),'sigma centers unchanged')
    call require(all(fev%FSS_1%TreeLRD_1%W_sgmas_mem==seeds),'sigma seeds unchanged')
    call require(fev%FSS_1%alpha==alpha0.and.fev%FSS_1%sgma_mean==mean0,'sigma parameters unchanged')
    call eval_density_grad_FMM(g,targets,normals,nt,f1,grad1,fev,mode)
    call function_eval_sigma(fev%FSS_1,targets,nt,s1,ds1(1,:),ds1(2,:),ds1(3,:),mode)
    call require(all(ieee_is_finite(f1)).and.all(ieee_is_finite(grad1)),'finite current field')
    call require(all(ieee_is_finite(s1)).and.all(ieee_is_finite(ds1)),'finite current sigma')
    print *, trim(label)
    print *, '  max delta phi, grad_phi:',maxval(abs(f1-f0)),maxval(abs(grad1-grad0))
    print *, '  max delta sigma, grad_sigma:',maxval(abs(s1-s0)),maxval(abs(ds1-ds0))
    call require(maxval(abs(f1-f0))<1.0d-12,'fixed phi invariant')
    call require(maxval(abs(grad1-grad0))<1.0d-12,'fixed phi gradient invariant')
    call require(maxval(abs(s1-s0))<1.0d-14,'fixed sigma invariant')
    call require(maxval(abs(ds1-ds0))<1.0d-14,'fixed sigma gradient invariant')
  end subroutine verify

  subroutine legacy_run(suffix,lambda,values)
    character(len=*), intent(in) :: suffix
    double precision, intent(in) :: lambda
    double precision, allocatable, intent(out) :: values(:)
    character(len=:), allocatable :: output
    integer*8 :: status,p,npatch,n
    integer :: unit,ios
    output=trim(root)//'_legacy_'//suffix
    call multiscale_mesher_unif_refine(trim(meshfile),3_8,1_8,trim(cadfile), &
        16_8,4_8,0_8,1_8,lambda,output,status)
    call require(status==0,'legacy public driver succeeds')
    open(newunit=unit,file=output//'_o04_r00.go3',status='old',iostat=ios)
    call require(ios==0,'legacy output exists')
    read(unit,*,iostat=ios) p
    call require(ios==0,'legacy order header')
    read(unit,*,iostat=ios) npatch
    call require(ios==0,'legacy patch header')
    n=12*npatch*(p+1)*(p+2)/2
    allocate(values(n))
    read(unit,*,iostat=ios) values
    close(unit)
    call require(ios==0.and.all(ieee_is_finite(values)),'legacy finite output')
  end subroutine legacy_run

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      stop 1
    endif
  end subroutine require
end program test_fixed_sources
