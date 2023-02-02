module global

	use iso_fortran_env, only : int64
	implicit none
	integer, parameter :: dp = SELECTED_REAL_KIND(15,99)

	! --- global variables
	character(len=2) :: typ_1, typ_2
	character(len=2), allocatable, dimension (:) :: t_atm
	integer :: max_frame, n_frame, dt, n_map,rdf_switch
	integer :: n_atm,n_bin,n_ow,n_i,n_k,i_trj,max_clust,i_clust,n_k_av,n_k_tot
	real(kind=dp), allocatable, dimension (:) :: lx,ly,lz
	real(kind=dp), allocatable, dimension (:,:) :: x,y,z,x_clu,y_clu,z_clu
	real(kind=dp), allocatable, dimension (:) :: x_ow,y_ow,z_ow
	real(kind=dp), allocatable, dimension (:) :: x_new,y_new,z_new,x_old,y_old,z_old
	real(kind=dp) :: vol_av, init_cut
	real(kind=dp), parameter :: pi = 3.14159265d0
	integer, allocatable, dimension (:) :: clust_dist, n_bnd
	integer, allocatable, dimension (:,:) :: clust_size
	integer :: cluster_switch, n_c

	namelist /input_deck/ max_frame, dt, typ_1, typ_2, n_bin, n_map, init_cut, cluster_switch

end module global

