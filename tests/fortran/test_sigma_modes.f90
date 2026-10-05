program test_sigma_modes
  use Mod_Feval
  implicit none
  type(Geometry) :: geom
  type(Fast_Sigma_stuff), pointer :: tree => null(), original => null()
  type(Feval_stuff) :: feval
  integer *8, parameter :: ntri=8, ntargets=4
  integer *8 :: i, j, mode
  double precision, parameter :: rlam=3.0d0
  double precision :: centers(3,ntri), longest(ntri), shortest(ntri), lengths(3)
  double precision :: targets(3,ntargets), sigma(ntargets), gradient(3,ntargets)
  double precision :: reference(ntargets), reference_gradient(3,ntargets)
  double precision :: origin(3), value, mean_seed
  double precision :: original_sigma(ntargets), original_gradient(3,ntargets)
  double precision :: saved_points(3,3*ntri)

  geom%ntri = ntri
  geom%ntri_sk = ntri
  geom%npoints = 3*ntri
  allocate(geom%Points(3,3*ntri),geom%Tri(3,ntri))
  do i=1,ntri
    origin = [2.0d0*mod(i-1,2_8),2.0d0*mod((i-1)/2,2_8),2.0d0*((i-1)/4)]
    geom%Tri(:,i) = [3*i-2,3*i-1,3*i]
    geom%Points(:,3*i-2) = origin
    geom%Points(:,3*i-1) = origin + [0.5d0+0.08d0*i,0.0d0,0.0d0]
    geom%Points(:,3*i) = origin + [0.0d0,1.2d0+0.13d0*i,0.0d0]
    centers(:,i) = sum(geom%Points(:,3*i-2:3*i),dim=2)/3.0d0
    lengths = [norm2(geom%Points(:,3*i-1)-geom%Points(:,3*i-2)), &
        norm2(geom%Points(:,3*i)-geom%Points(:,3*i-1)), &
        norm2(geom%Points(:,3*i)-geom%Points(:,3*i-2))]
    longest(i) = maxval(lengths)/rlam
    shortest(i) = minval(lengths)/rlam
  enddo
  saved_points = geom%Points
  targets(:,1) = [0.1d0,0.2d0,0.3d0]
  targets(:,2) = [2.6d0,0.8d0,0.4d0]
  targets(:,3) = [1.1d0,2.4d0,1.6d0]
  targets(:,4) = [2.2d0,2.7d0,2.3d0]

  call setup_tree_sigma_geometry(original,geom,rlam)
  call function_eval_sigma(original,targets,ntargets,original_sigma, &
      original_gradient(1,:),original_gradient(2,:),original_gradient(3,:),1_8)
  do mode=0,3
    call setup_tree_sigma_geometry(tree,geom,rlam,mode)
    do i=1,ntri
      value = huge(1.0d0)
      do j=1,ntri
        if (norm2(tree%TreeLRD_1%W_Pts_mem(:,i)-centers(:,j)).lt.1.0d-13) then
          value = longest(j)
          if (mode.eq.3) value = shortest(j)
          exit
        endif
      enddo
      call require(abs(tree%TreeLRD_1%W_sgmas_mem(i)-value).lt.1.0d-14,'triangle seeds')
    enddo
    call function_eval_sigma(tree,targets,ntargets,sigma, &
        gradient(1,:),gradient(2,:),gradient(3,:),mode)
    call require(all(sigma.gt.0.0d0),'positive sigma')
    if (mode.eq.0) then
      mean_seed = sum(longest)/ntri
      call require(maxval(abs(sigma-mean_seed)).lt.1.0d-14,'constant mean longest-side sigma')
      call require(all(gradient.eq.0.0d0),'zero constant gradient')
    else if (mode.eq.1.or.mode.eq.3) then
      if (mode.eq.1) then
        call direct_gaussian(longest,reference,reference_gradient)
        call require(all(sigma.eq.original_sigma),'omitted mode preserves longest-side sigma')
        call require(all(gradient.eq.original_gradient),'omitted mode preserves gradient')
      else
        call direct_gaussian(shortest,reference,reference_gradient)
      endif
      call require(maxval(abs(sigma-reference)).lt.1.0d-12,'direct Gaussian sigma')
      call require(maxval(abs(gradient-reference_gradient)).lt.1.0d-12,'direct Gaussian gradient')
    else
      call function_eval_sigma(original,targets,ntargets,reference, &
          reference_gradient(1,:),reference_gradient(2,:),reference_gradient(3,:),2_8)
      call require(all(sigma.eq.reference),'recursive sigma unchanged')
      call require(all(gradient.eq.reference_gradient),'recursive gradient unchanged')
    endif
    call require(all(geom%Points.eq.saved_points),'sigma initialization does not change scaffold')
    call destroy_fast_sigma(tree)
  enddo

  call start_Feval_tree(feval,geom,rlam)
  call function_eval_sigma(feval%FSS_1,targets,ntargets,sigma, &
      gradient(1,:),gradient(2,:),gradient(3,:),1_8)
  call require(all(sigma.eq.original_sigma),'legacy Feval initializer')
  call destroy_Feval_tree(feval)
  call start_Feval_tree(feval,geom,rlam,3_8)
  call function_eval_sigma(feval%FSS_1,targets,ntargets,sigma, &
      gradient(1,:),gradient(2,:),gradient(3,:),3_8)
  call direct_gaussian(shortest,reference,reference_gradient)
  call require(maxval(abs(sigma-reference)).lt.1.0d-12,'shortest-side Feval initializer')
  call destroy_Feval_tree(feval)
  call destroy_fast_sigma(original)
  print *, 'PASS: constant, longest, recursive and shortest sigma; direct gradients; optional initializers'

contains
  subroutine direct_gaussian(seeds,values,derivatives)
    double precision, intent(in) :: seeds(ntri)
    double precision, intent(out) :: values(ntargets), derivatives(3,ntargets)
    double precision :: weights(ntri), dweights(3,ntri), delta(3), width, denominator
    integer :: k, t
    width = sqrt(5.0d0/2.0d0)*maxval(seeds)
    do t=1,ntargets
      do k=1,ntri
        delta = targets(:,t)-centers(:,k)
        weights(k) = exp(-sum(delta**2)/(2.0d0*width**2))
        dweights(:,k) = -delta/width**2*weights(k)
      enddo
      denominator = sum(weights)
      values(t) = sum(weights*seeds)/denominator
      do k=1,3
        derivatives(k,t) = sum(dweights(k,:)*(seeds-values(t)))/denominator
      enddo
    enddo
  end subroutine direct_gaussian

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ',message
      stop 1
    endif
  end subroutine require
end program test_sigma_modes
