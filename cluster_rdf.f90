program cluster_rdf
	use global
	implicit none

	! --- read input file
	open(unit = 1, file ='input.dat', status = 'old')
		read(unit = 1, nml = input_deck)
	close(unit = 1) 

	! --- initialise random number generator
	call rand_init

	! --- calculate number of frames to be used
	n_frame = max_frame / dt

	! --- read in trajectory file (.gro format)
	call read_trj

        ! --- assign clusters to each frame
        call cluster

	! --- calculate 1:1 atomistic RDF
	rdf_switch = 1
	call calc_rdf

	! --- calculate 4:1 atomistic RDF
	rdf_switch = 2
	call calc_rdf

end program

subroutine read_trj
        use global
	implicit none
        ! --- variables
        character(len=13) :: A, nam
        integer :: i, i_atm, atm_i, i_frame
	real(kind=dp), allocatable, dimension (:) :: vol
	logical :: first_pass 

        ! --- read in initial coordinates for 3D periodic membrane model
        100 format (a13, a2, i5, f8.3, f8.3, f8.3)
        101 format (f10.5, f10.5, f10.5)

        ! --- allocate frame arrays
        allocate(vol(max_frame),lx(max_frame),ly(max_frame),lz(max_frame))

        open(unit=10,file='output.gro')
	first_pass = .true.
        do i_frame = 1, max_frame
                read(10,*) A
                read(10,*) n_atm
                ! --- allocate atom arrays in first frame
                if (first_pass.eqv..true.) then
			allocate(n_bnd(n_atm),x(max_frame,n_atm),y(max_frame,n_atm),z(max_frame,n_atm),t_atm(n_atm))
        		if ((typ_1.eq.'CA').or.(typ_1.eq.'CB')) then
                		open(unit=50,file='bonding.dat')
                		read(50,*) n_c
                		n_bnd(:) = 0
                		do i = 1, n_c
                       			read(50,*) n_bnd(i)
                		end do
                		close(unit=50)
        		end if
		end if

                do i_atm = 1, n_atm
                        read(10,100) nam, t_atm(i_atm), atm_i, x(i_frame, i_atm), y(i_frame, i_atm), z(i_frame, i_atm)
			if (t_atm(i_atm).eq.'CA') then
				if (n_bnd(i_atm).eq.2) then
					t_atm(i_atm) = 'CB'
				end if
			end if
                        if (first_pass.eqv..true.) then
                                if (t_atm(i_atm).eq.typ_1) then
                                      	n_i = n_i + 1
				else if (t_atm(i_atm).eq.typ_2) then
					n_ow = n_ow + 1
				end if
			end if
                end do
                read(10,*) lx(i_frame), ly(i_frame), lz(i_frame)
		vol(i_frame) = lx(i_frame) * ly(i_frame) * lz(i_frame)
		first_pass = .false.
        end do
	vol_av = sum(vol) / real(max_frame)
        close(unit=10)
end subroutine

! needs to made more general so same subroutine can be used for clusters
subroutine calc_rdf
	use global
	implicit none
	! --- variables
	real(kind=dp) :: bin_width, rx, ry, rz, rxyz, bin_vol, max_bin, rcut
	integer :: i_k, i_frame, i_atm, j_atm, i_bin
	real(kind=dp), allocatable, dimension (:) :: g

	200 format (f8.3, f8.3)	

	allocate(g(n_bin))

	rcut = 2.0d0
	bin_width = rcut / real(n_bin)

	! --- write histograms for Dirac delta function
	do i_frame = 1, n_frame
		i_trj = i_frame * dt
		! --- atomistic RDF
		if (rdf_switch.eq.1) then
			do i_atm = 1, n_atm - 1
				if(t_atm(i_atm).eq.typ_1) then
					do j_atm = i_atm + 1, n_atm
						if(t_atm(j_atm).eq.typ_2) then
							rx = x(i_trj,i_atm) - x(i_trj,j_atm)
							ry = y(i_trj,i_atm) - y(i_trj,j_atm)
							rz = z(i_trj,i_atm) - z(i_trj,j_atm)
							call pbc(rx,ry,rz,lx(i_trj),ly(i_trj),lz(i_trj))
							rxyz = sqrt(rx*rx+ry*ry+rz*rz)
							if (rxyz.lt.rcut) then
								i_bin = ceiling(rxyz/bin_width)
								g(i_bin) = g(i_bin) + 1
							end if
						end if
					end do
				end if
			end do
		! --- CG water RDF - only those with 4 water molecules
		else if (rdf_switch.eq.2) then
                        do i_atm = 1, n_atm
                                if(t_atm(i_atm).eq.typ_1) then
                                        do i_k = 1, n_k
						! --- always include 4 water clusters
						!if (clust_size(i_trj,i_k).eq.4) then
                                        		rx = x(i_trj,i_atm) - x_clu(i_trj,i_k)
                                                	ry = y(i_trj,i_atm) - y_clu(i_trj,i_k)
                                                	rz = z(i_trj,i_atm) - z_clu(i_trj,i_k)
                                                	call pbc(rx,ry,rz,lx(i_trj),ly(i_trj),lz(i_trj))
                                                	rxyz = sqrt(rx*rx+ry*ry+rz*rz)
                                                	if (rxyz.lt.rcut) then
                                                		i_bin = ceiling(rxyz/bin_width)
                                                		g(i_bin) = g(i_bin) + 1
                                                	end if	
						!end if
						! --- include 3 water clusters
						!if ((cluster_switch.eq.34).or.(cluster_switch.eq.345)) then
                                                !	if (clust_size(i_trj,i_k).eq.3) then
                                            	 !        	rx = x(i_trj,i_atm) - x_clu(i_trj,i_k)
                                                !	        ry = y(i_trj,i_atm) - y_clu(i_trj,i_k)
                                                 !      		rz = z(i_trj,i_atm) - z_clu(i_trj,i_k)
                                                  !      	call pbc(rx,ry,rz,lx(i_trj),ly(i_trj),lz(i_trj))
                                                   !    		rxyz = sqrt(rx*rx+ry*ry+rz*rz)
                                                    !    	if (rxyz.lt.rcut) then
                                                     !           	i_bin = ceiling(rxyz/bin_width)
                                                      !          	g(i_bin) = g(i_bin) + 1
                                                       ! 	end if
						!	end if
						!end if
						! --- include 5 water clusters
                                                !if ((cluster_switch.eq.35).or.(cluster_switch.eq.345)) then
                                                !        if (clust_size(i_trj,i_k).eq.5) then
                                                !    		rx = x(i_trj,i_atm) - x_clu(i_trj,i_k)
                                                !       	 	ry = y(i_trj,i_atm) - y_clu(i_trj,i_k)
                                                !       	 	rz = z(i_trj,i_atm) - z_clu(i_trj,i_k)
                                                !        	call pbc(rx,ry,rz,lx(i_trj),ly(i_trj),lz(i_trj))
                                                !        	rxyz = sqrt(rx*rx+ry*ry+rz*rz)
                                                !        	if (rxyz.lt.rcut) then
                                                !                	i_bin = ceiling(rxyz/bin_width)
                                                !                	g(i_bin) = g(i_bin) + 1
                                                !        	end if
						!	end if
                                                !end if
                                        end do
                                end if
                        end do
		end if
	end do

	! --- write radial distribution function
	if (rdf_switch.eq.1) then
		open(unit=20,file='rdf_atom.xvg')
	else if (rdf_switch.eq.2) then
		open(unit=20,file='rdf_CG.xvg')	
	end if
	do i_bin = 1, n_bin
		! --- calculate bin volume
		bin_vol = (4.0/3.0) * pi * (i_bin**3 - (i_bin-1)**3) * bin_width**3
	        open(unit=10,file='clust_size.dat')	! --- normalise and plot radial distribution function
		if (rdf_switch.eq.1) then
			g(i_bin) = vol_av * g(i_bin) / bin_vol / real(n_i) / real(n_ow) / real(n_frame)
		else if (rdf_switch.eq.2) then	
			g(i_bin) = vol_av * g(i_bin) / bin_vol / real(n_i) / real(n_k_av) / real(n_frame)
		end if
		write(20,200) (i_bin - 0.5d0) * bin_width, g(i_bin)
	end do
	close(unit=20)
end subroutine

subroutine cluster
	use global
	implicit none
	! --- variables
	integer :: i_k, i_frame,i_ow, count, count_k, i_atm,j_ow
	real(kind=dp) :: ran, rx, ry, rz, rxyz
	logical :: first_pass
	logical :: assigned(n_ow)
	integer :: max_cycle, cycles

	first_pass = .true.	
	max_cycle = 100000

	! --- initialise k_means algorithm
	n_k = n_ow / n_map
	max_clust = n_map + 10

	allocate(clust_size(n_frame,n_k))
	allocate(clust_dist(0:max_clust))
	clust_size(:,:) = 0
	allocate(x_old(n_k),y_old(n_k),z_old(n_k),x_new(n_k),y_new(n_k),z_new(n_k))
	allocate(x_clu(n_frame,n_k),y_clu(n_frame,n_k),z_clu(n_frame,n_k))
	allocate(x_ow(n_atm),y_ow(n_atm),z_ow(n_atm))

	! --- loop over all frames
	do i_frame = 1, n_frame
                i_trj = i_frame * dt
		count = 0
		! --- write new atom position array for this frame
		do i_atm = 1, n_atm
			if (t_atm(i_atm).eq.typ_2) then
                                count = count + 1
                                x_ow(count) = x(i_trj,i_atm)
                                y_ow(count) = y(i_trj,i_atm)
                                z_ow(count) = z(i_trj,i_atm)
				assigned(count) = .false.
			end if
                end do
		! --- on first frame use random OW coordinates to initialise clusters
		if (first_pass.eqv..true.) then
			count_k = 0
			cycles = 0
			do
				cycles = cycles + 1
                                if (cycles.ge.max_cycle) then
                                        print*, "Number of cycles exceeds limit. Reduce init_cut."
                                        stop
                                end if
                		call random_number(ran)
                		i_ow = ceiling(ran * n_ow)
				! --- if water molecule has already been selecte find a new one
                		if (assigned(i_ow)) cycle
				count_k = count_k + 1
                                assigned(i_ow) = .true.
                                x_old(count_k) = x_ow(i_ow)
                                y_old(count_k) = y_ow(i_ow)
                                z_old(count_k) = z_ow(i_ow)
				! --- prevent any neighbouring waters (within cutoff) from future selection
				do j_ow = 1, n_ow
					if (i_ow.ne.j_ow) then
						rx = x_ow(i_ow) - x_ow(j_ow)
                                                ry = y_ow(i_ow) - y_ow(j_ow)
                                                rz = z_ow(i_ow) - z_ow(j_ow)
                                		call pbc(rx,ry,rz,lx(i_trj),ly(i_trj),lz(i_trj))
                                		rxyz = sqrt(rx*rx+ry*ry+rz*rz)
						if (rxyz.lt.init_cut) then
							assigned(j_ow) = .true.
						end if
					end if					
				end do
                		if (count_k.eq.n_k) exit
			end do
		end if
		first_pass = .false.

		call k_means

		! --- use centres of mass of clusters from previous frame to initialise
		do i_k = 1, n_k
			x_clu(i_trj,i_k) = x_new(i_k)
               		y_clu(i_trj,i_k) = y_new(i_k)
                	z_clu(i_trj,i_k) = z_new(i_k)
			x_old(i_k) = x_new(i_k)
			y_old(i_k) = y_new(i_k)
			z_old(i_k) = z_new(i_k)
		end do

	end do

	! --- calculate cluster size distribution
	n_k_tot = 0	
	clust_dist(:) = 0
	open(unit=10,file='clust_size.dat')
        do i_clust = 0, max_clust
		do i_frame = 1, n_frame
                	i_trj = i_frame * dt
                        do i_k = 1, n_k
				if (i_clust.eq.clust_size(i_trj,i_k)) then
	                        	clust_dist(i_clust) = clust_dist(i_clust) + 1
			!		if (i_clust.eq.4) then
						n_k_tot = n_k_tot + 1
			!		end if
			!		if ((cluster_switch.eq.3).or.(cluster_switch.eq.345)) then
			!			if (i_clust.eq.3) then
                         !                       	n_k_tot = n_k_tot + 1
                          !              	end if
			!		end if
                         !               if ((cluster_switch.eq.5).or.(cluster_switch.eq.345)) then
                          !                      if (i_clust.eq.5) then
                           !                             n_k_tot = n_k_tot + 1
                            !                    end if
                             !           end if
				end if
                        end do
		end do
                write(10,*) i_clust, real(clust_dist(i_clust)) / real(n_k) / real(n_frame)
	end do
	n_k_av = real(n_k_tot) / real(n_frame)
	close(unit=10)

end subroutine
		
subroutine k_means
	use global
	implicit none
	! --- implicit none
	real(kind=dp) :: k_tol, converge, rxyz2_tot, twopi
	real(kind=dp), allocatable, dimension (:) :: rxyz2
	real(kind=dp) :: x_xi, x_zeta, y_xi, y_zeta, z_xi, z_zeta
	integer, allocatable, dimension (:) :: min, clust, k_size
	integer :: max_cycle, cycles, i_ow, i_k, count, tot_count	
	real(kind=dp) :: rx, ry, rz

	max_cycle = 100
	k_tol = 1.0d-3
	twopi = 2.0d0*pi

        allocate (rxyz2(n_k),min(n_k),clust(n_ow),k_size(n_k))

	open(unit=10,file='kmeans.dat')
	write(10,*) "*******"
        write(10,*) "Frame:", i_trj
	cycles = 0
	! --- repeat until convergence criteria met
	do
		cycles = cycles + 1
		k_size(:) = 0
		! --- assign each atom to a cluster
		do i_ow = 1, n_ow 
			do i_k = 1, n_k
                                rx = x_ow(i_ow) - x_old(i_k)
                                ry = y_ow(i_ow) - y_old(i_k)
                                rz = z_ow(i_ow) - z_old(i_k)
				call pbc(rx,ry,rz,lx(i_trj),ly(i_trj),lz(i_trj))
				rxyz2(i_k) = rx*rx+ry*ry+rz*rz
			end do
			min = minloc(rxyz2)
			clust(i_ow) = min(1)
		end do

		! --- count the number of atoms assigned to each cluster and determine new cluster CoM
		do i_k = 1, n_k
			x_xi = 0.0d0
			x_zeta = 0.0d0
			y_xi = 0.0d0
			y_zeta = 0.0d0
			z_xi = 0.0d0
			z_zeta = 0.0d0
			count = 0
			do i_ow = 1, n_ow
				if (clust(i_ow).eq.i_k) then
					count = count + 1
					! --- map particle coordinates from a line to a circle
					x_xi = x_xi + cos(x_ow(i_ow)*twopi/lx(i_trj))
					x_zeta = x_zeta + sin(x_ow(i_ow)*twopi/lx(i_trj))
                                        y_xi = y_xi + cos(y_ow(i_ow)*twopi/ly(i_trj))
                                        y_zeta = y_zeta + sin(y_ow(i_ow)*twopi/ly(i_trj))
                                        z_xi = z_xi + cos(z_ow(i_ow)*twopi/lz(i_trj))
                                        z_zeta = z_zeta + sin(z_ow(i_ow)*twopi/lz(i_trj))
				end if
			end do	
			if (count.ne.0) then
				! --- map coordinates back from circle to a line
				x_xi = x_xi / count
				x_zeta = x_zeta / count
				x_new(i_k) = (atan2(-1.0d0*x_zeta,-1.0d0*x_xi)+pi)*lx(i_trj)/twopi
                        	y_xi = y_xi / count
                        	y_zeta = y_zeta / count
                        	y_new(i_k) = (atan2(-1.0d0*y_zeta,-1.0d0*y_xi)+pi)*ly(i_trj)/twopi
                        	z_xi = z_xi / count
                        	z_zeta = z_zeta / count
                        	z_new(i_k) = (atan2(-1.0d0*z_zeta,-1.0d0*z_xi)+pi)*lz(i_trj)/twopi
			end if
			! --- determine centre of mass change since last iteration
			rx = x_new(i_k) - x_old(i_k)
			ry = y_new(i_k) - y_old(i_k)
			rz = z_new(i_k) - z_old(i_k)
			call pbc(rx,ry,rz,lx(i_trj),ly(i_trj),lz(i_trj))
			rxyz2(i_k) = rx*rx+ry*ry+rz*rz
                        x_old(i_k) = x_new(i_k)
                        y_old(i_k) = y_new(i_k)
                        z_old(i_k) = z_new(i_k)
			k_size(i_k) = count
		end do
		rxyz2_tot = sum(rxyz2)
		write(10,*) cycles, rxyz2_tot
		! --- convergence test
		if (rxyz2_tot.le.k_tol) then
			! --- save all clusters with 4 water molecules
			do i_k = 1, n_k
				clust_size(i_trj,i_k) = k_size(i_k)
			end do
			exit
		end if
		if (cycles.gt.max_cycle) then
			write(10,*) "Could not converge K-means algorithm"
			write(10,*) "Increase number of cycles"
			stop
		end if
	end do

	! some sort of test to check that there are n_map waters in each cluster
	! WHY should we end up with n_map waters in each cluster - easy to get stuck in local minima if initial waters not chosen sensibly

	write(10,*) "K-means algorithm converged in", cycles, "cycles"
	write(10,*)
	
end subroutine

subroutine pbc(rx,ry,rz,lx,ly,lz)
	use global, only: dp
	implicit none
	real(kind=dp) :: rx,ry,rz,lx,ly,lz
        rx = rx - lx * nint(rx/lx)
        ry = ry - ly * nint(ry/ly)
        rz = rz - lz * nint(rz/lz)
end subroutine

subroutine rand_init
        implicit none
        integer :: n_seed, i_seed
        integer, allocatable, dimension (:) :: seed
        integer :: count
        call system_clock(count)
        call random_seed(size=n_seed)
        allocate (seed(n_seed))
        do i_seed = 1, n_seed
                seed(i_seed) = lcg(count)
        end do
        call random_seed(put=seed)
        Contains
                function lcg(s)
                use iso_fortran_env, only: int64
                integer :: lcg
                integer :: s
                if (s.eq.0) then
                        s = 104729
                else
                        s = mod(s,4294967296_int64)
                end if
                s = mod(s * 279470273_int64, 4294967291_int64)
                lcg = int(mod(s,int(huge(0),int64)), kind(0))
                end function lcg
end subroutine rand_init  

