
program smoother

  implicit none

  integer *8 :: N, count,nrefine, ifplot
  integer *8 :: adapt_flag, ifflatten
  integer *8 :: interp_flag
  integer *8 :: norder_skel, norder_smooth
  integer *8 :: ifcad

  character (len=4096) :: fnamein, fnameout_root,name_aux, fcad
  character (len=21) :: plot_name
  character (len=8) :: istr1,istr2
  character (len=2) :: arg_comm

  real *8 rlam
  integer *8 i, ier, ifiletype
  


  call prini(6,13)

  call test_surf_smooth_ker_small_u()


  ! order with which to discretize the skeleton patches (pick
  ! something high-order)
  norder_skel = 16

  ! order with which to discretize the smooth patches, choose
  ! something reasonable: 4, 6, 8, 10, etc.
  norder_smooth = 8

  ! Define number of refinements of smooth surface to be output in
  ! the go3 format
  ! 
  !
  nrefine = 1
  ! nrefine=1

  ! this is to enable adaptativity (otherwise sigma is constant)
  ! adapt_flag = 0  ->  no adaptivity, mean triangle size
  ! adapt_flag = 1  ->  some adaptivity, alpha form
  ! adapt_flag = 2  ->  full recursive definition, slightly slower
  adapt_flag = 1


  !
  !
  ! rlam flag decides the proportionality value for \sigma_{j}
  ! in relation to triangle diameter
  ! \sigma_{j} = D_{j}/rlam
  !
!  rlam = 10.0d0 !(usual value)
  rlam = 5.0d0

  !rlam = .5d0
  !rlam = 1
  !rlam = 2.5d0


!
! specify the msh file to read in
!

!    fnamein='../../geometries/meshes/cuboid_a1_b2_c1p3.tri'
!    fnameout_root='../../geometries/cuboid_a1_b2_c1p3'
!    ifiletype = 2    
    fcad = 'tmp'

!    fnamein='../../geometries/meshes/prism_50.gidmsh'
!    fnameout_root='../../geometries/prism_50'
!    ifiletype = 3

!    fnamein='../../geometries/meshes/sphere.msh'
!    fnameout_root='../../geometries/sphere'
 
!    fnamein = '../../geometries/meshes/cow_new.msh'
!    fnameout_root = '../../geometries/cow_new'

!    fnamein = '../../geometries/meshes/cow_new_gmshv4.msh'
!    fnameout_root = '../../geometries/cow_new'

!    fnamein = '../../geometries/meshes/lens_r00.msh'
!    fnameout_root = '../../geometries/lens_r00'

!    fnamein = '../../geometries/meshes/lens_r00_gmshv4.msh'
!    fnameout_root = '../../geometries/lens_r00'
!
     fnamein = '../../geometries/meshes/cylinder.gidmsh'
     fnameout_root = '../../geometries/cylinder'
     ifcad = 1
     fcad = '../../geometries/meshes/cylinder_skeleton.txt'

! The extracted test receives fixture/output paths from its build driver.
     if (command_argument_count()/=3) stop 'Expected mesh, CAD and output paths'
     call get_command_argument(1,fnamein)
     call get_command_argument(2,fcad)
     call get_command_argument(3,fnameout_root)

!    Determine file type
!        * ifiletype = 1, for .msh
!        * ifiletype = 2, for .tri
!        * ifiletype = 3, for .gidmsh
!        * ifiletype = 4, for .msh gmsh v2
!        * ifiletype = 5, for .msh gmsh v4
    ier = 0
    call get_filetype(fnamein, ifiletype, ier)

    if(ier.ne.0) then
      print *, "File type not recognized"
      stop
    endif



    ier = 0
    call multiscale_mesher_unif_refine(fnamein, ifiletype, ifcad, fcad, &
       norder_skel, norder_smooth, nrefine, adapt_flag, rlam, &
       fnameout_root, ier)
    if (ier/=0) stop 1

contains

subroutine test_surf_smooth_ker_small_u()
  implicit none

  integer :: i, j
  double precision :: sigma, r, h, hh, hhh, tol
  double precision, parameter :: sq2 = &
      1.414213562373095048801688724209698079d0
  double precision, parameter :: sigma_vals(5) = &
      (/1.0d-6, 1.0d-4, 1.0d0, 1.0d4, 1.0d6/)
  double precision, parameter :: u_vals(12) = &
      (/0.0d0, 1.0d-6, 0.0099d0, 0.01d0, 0.0101d0, 0.02d0, &
        0.05d0, 0.099999d0, 0.1d0, 0.100001d0, 0.2d0, 1.0d0/)

  ! Independent references from 100-digit quadrature on 0 <= t <= 1:
  ! h*sigma**3 = sqrt(2/pi) * integral t**2*exp(-u**2*t**2) dt
  ! hh*sigma**5 = -sqrt(2/pi) * integral t**4*exp(-u**2*t**2) dt
  ! Include zero, both sides of the old and new cutoffs, and larger u.
  double precision, parameter :: h_ref(12) = (/ &
      0.26596152026762178529d0, 0.26596152026746220838d0, &
      0.26594588068190753508d0, 0.26594556314630849608d0, &
      0.26594524241985508199d0, 0.26589769862049268736d0, &
      0.26556293395493362979d0, 0.26437146726412553417d0, &
      0.26437143557598217326d0, 0.26437140388752645338d0, &
      0.25966869263217634284d0, 0.15117705942927216154d0/)
  double precision, parameter :: hh_ref(12) = (/ &
      -0.15957691216057307118d0, -0.15957691216045908767d0, &
      -0.15956574106267798804d0, -0.15956551425296157036d0, &
      -0.15956528516410687459d0, -0.15953132584863218809d0, &
      -0.15929223024331476783d0, -0.15844152032209126768d0, &
      -0.15844149770208524444d0, -0.15844147508185653729d0, &
      -0.15508772768057870723d0, -0.080002925970168342412d0/)

  do i=1,size(sigma_vals)
    sigma = sigma_vals(i)
    do j=1,size(u_vals)
      r = sq2*sigma*u_vals(j)
      call surf_smooth_ker(r, sigma, h, hh, hhh)
      call check_kernel_value('h', sigma, u_vals(j), &
          h*sigma**3, h_ref(j), 1.0d-12)
      tol = 1.0d-12
      ! The direct formula above the cutoff loses some precision.
      if (u_vals(j).ge.0.1d0) tol = 5.0d-11
      call check_kernel_value('hh', sigma, u_vals(j), &
          hh*sigma**5, hh_ref(j), tol)
    enddo
  enddo

  print *, 'surf_smooth_ker small-u regression: PASS'

end subroutine test_surf_smooth_ker_small_u


subroutine check_kernel_value(label, sigma, u, value, reference, tolerance)
  implicit none

  character (len=*), intent(in) :: label
  double precision, intent(in) :: sigma, u, value, reference, tolerance
  double precision :: relerr

  relerr = abs(value-reference)/abs(reference)
  ! This form also rejects NaN results.
  if (.not.(relerr.le.tolerance)) then
    print *, 'surf_smooth_ker regression failure: ', trim(label)
    print *, 'sigma=', sigma, ' u=', u
    print *, 'computed=', value, ' reference=', reference
    print *, 'relative error=', relerr, ' tolerance=', tolerance
    stop 1
  endif

end subroutine check_kernel_value

end program smoother
