module sor

    !This module uses iterative methods to solve the Poisson equation
    ! Successive Over Relaxation

    implicit none

    contains

        subroutine sor_cgrid(f, e1u, e2u, e1v, e2v, e1t, e2t, phi_zero, bd_iminvals, &
            bd_imaxvals, bd_jminvals, bd_jmaxvals, bd_iminnm, bd_imaxnm, bd_jminnm, &
            bd_jmaxnm, bd_ipd, bd_jpd, tmask, omega, thresh, niterations_max, phi)
            !
            !sor_cgrid
            !
            !Subroutine that solves the poisson equation on a C grid. This function has been made to solve
            !for the compressible component of the barotropic flow where
            !
            ! U = k x grad(psi) + grad(phi)
            !
            ! div(U) = lap(phi)
            !
            !The subroutine uses Successive Over-Relaxation to solve the poisson equation. This is quite old fashioned
            !and not the fastesr but it makes up for that in simplicity.
            !
            !Variable declaration below
            !
            real, intent(in) :: f(:,:)    !Forcing array of poisson equation (y,x)
            real, intent(in) :: e1u(size(f,1), size(f,2))    !ArrayS of u cell widths      (y,x)
            real, intent(in) :: e2u(size(f,1), size(f,2))
            real, intent(in) :: e1v(size(f,1), size(f,2))    !Arrays of v cell widths      (y,x)
            real, intent(in) :: e2v(size(f,1), size(f,2))
            real, intent(in) :: e1t(size(f,1), size(f,2))    !Arrays of t cell widths      (y,x)
            real, intent(in) :: e2t(size(f,1), size(f,2))
            real, intent(in) :: phi_zero(size(f,1), size(f,2))   !Initial solution for phi (y,x)
            real, intent(in) :: bd_jminvals(size(f,1))       !Boundary values: - If using a dirichlet boundary condition
            real, intent(in) :: bd_jmaxvals(size(f,1))       !                   the vector describes the constant values 
                                                             !                   at the boundary.
            real, intent(in) :: bd_iminvals(size(f,2))       !                 - If using a Neumann boundary condition the 
            real, intent(in) :: bd_imaxvals(size(f,2))       !                   vector describes the gradient at the boundary.
            !
            !                                                                   !Boundary Classes >>>>
            logical, intent(in) :: bd_iminnm, bd_imaxnm,  bd_jminnm, bd_jmaxnm  !If True, applies a Neumann 
                                                                                !boundary condition.

            logical, intent(in) :: bd_ipd, bd_jpd                               !If True, applies a periodic boundary condition in the 
                                                                                !i direction (bd_ipd==True) and/or j direction (bd_ipd==True)
                                                                                !Periodic boundary conditions will override Neumann boundary conditions
            !
                                                                                !Any boundary that is not Neumann or periodic will have
                                                                                !a Dirichlet boundary condition

            logical, intent(in) :: tmask(size(f,1), size(f,2))   !The mask for t poins on the c grid
            real, intent(in) :: omega                            !Relaxation parameter for Successive Over-Relaxation
            real, intent(in) :: thresh                           !Threshold for satisfactory solution. Iterations teminate 
                                                                 !when sum(res)/sum(f) < thresh 
            integer, intent(in) :: niterations_max
            !
            integer :: i,j,imax,jmax, ip, jp
            integer :: niterations
            logical :: correction_log
            logical :: land_mask_mu(size(f,1), size(f,2))        
            logical :: land_mask_nu(size(f,1), size(f,2))
            real :: mu(size(f,1), size(f,2))
            real :: nu(size(f,1), size(f,2))
            real :: delta(size(f,1), size(f,2))
            real :: f_mod(size(f,1), size(f,2))
            real :: e1v_col(size(f,2))
            real :: e2u_col(size(f,1))
            double precision :: sum_f
            real :: res, res_norm
            logical :: kbd_iminnm, kbd_imaxnm,  kbd_jminnm, kbd_jmaxnm
            !
            real, intent(out) :: phi(size(f,1), size(f,2))
            !
            imax = size(f,2)
            jmax = size(f,1)
            !
            !Force field will need to be modified if Neumann BC's are used or there is a tmask
            f_mod = f
            PRINT *, "jminval(69): ", bd_jminvals(69)
            PRINT *, "SUM_F = ", sum(abs(f_mod))
            !
            !Construct mu and nu.
            mu = e2u/e1u
            nu = e1v/e2v
            !
            !Periodic boundaries override Neumann boundary conditions
            kbd_iminnm = (bd_iminnm.AND.(.NOT.(bd_ipd))) !--> Set to false if periodic
            kbd_imaxnm = (bd_imaxnm.AND.(.NOT.(bd_ipd)))
            kbd_jminnm = (bd_imaxnm.AND.(.NOT.(bd_jpd)))
            kbd_jmaxnm = (bd_imaxnm.AND.(.NOT.(bd_jpd)))

            !If using Neuann boundary conditions: modify mu, nu and f accordingly
            if (kbd_iminnm) then
                !
                do j=1,jmax
                    mu(j,1) = 0.0          ! VV Non-zero Neumann BC's are currently broken VV
                    f_mod(j,2) = f_mod(j,2)! + e2u(j,1)*bd_iminvals(j)
                end do
                !
            end if
            !
            PRINT *, "SUM_F = ", sum(abs(f_mod))
            if (kbd_imaxnm) then
                !
                do j=1,jmax
                    mu(j,imax-1) = 0.0               ! VV Non-zero Neumann BC's are currently broken VV
                    f_mod(j,imax-1) = f_mod(j,imax-1)! + e2u(j,imax-1)*bd_imaxvals(j)
                end do
                !
            end if
            !
            PRINT *, "SUM_F = ", sum(abs(f_mod))
            if (kbd_jminnm) then

                do i=1,imax
                    nu(1,i) = 0.0          ! VV Non-zero Neumann BC's are currently broken VV
                    f_mod(2,i) = f_mod(2,i)! + e1v(1,i)*bd_jminvals(i)
                end do
                !
            end if
            !
            PRINT *, "SUM_F = ", sum(abs(f_mod))
            if (kbd_jmaxnm) then
                !
                do i=1,imax
                    nu(jmax-1,i) = 0.0
                    f_mod(jmax-1,i) = f_mod(jmax-1,i)! + e1v(jmax-1,i)*bd_jminvals(i)
                end do
                !
            end if

            !Treat any tmask points like land and update mu and nu for no normal flow
            PRINT *, "SUM_F = ", sum(abs(f_mod))
            where (tmask)
                !
                f_mod = 0
                !
            end where
            PRINT *, "SUM_F = ", sum(abs(f_mod))
            land_mask_mu = tmask
            land_mask_nu = tmask

            do  j = 1,jmax-1
                do i = 1,imax-1
                    !
                    land_mask_mu(j,i) = (tmask(j,i) .or. tmask(j,i+1))
                    land_mask_nu(j,i) = (tmask(j,i) .or. tmask(j+1,i))
                    !
                end do
            end do

            where (land_mask_mu)
                !
                mu = 0.0
                !
            end where 

            where (land_mask_nu)
                !
                nu = 0.0
                !
            end where

            !Scale f according to the t cell area
            PRINT *, "SUM_F = ", sum(abs(f_mod))
            f_mod = f_mod*e1t*e2t

            !If all boundaries have Neumann boundary conditions or periodic, set one point on the boundary to
            ! Dirichlet. This is required for stability in the solver as the grid needs a connection to
            ! the real world
            correction_log =  ( ((kbd_iminnm.AND.kbd_imaxnm).OR.(bd_ipd))    &
                           .AND.((kbd_jminnm.AND.kbd_jmaxnm).OR.(bd_jpd)) ) 

            if ( correction_log ) then
                !Find an unmasked point and save it's index.
                !This is needed to keep the potential bounded to small values
                ip=1
                jp=1
                do i = 1,imax
                    do j = 1,jmax
                        if (.NOT.tmask(j,i)) then !If the point is unmasked
                            jp = j !(jp, ip) is the location of an unmasked point
                            ip = i
                            exit
                        end if
                    end do
                end do
                print *, "jp = ", jp, " ip = ", ip
            end if


            !Now that mu and nu have been properly constructed, construct delta
            delta = -(mu + nu)

            do i = 2,imax
                do j = 2,jmax

                    delta(j,i) = delta(j,i) - mu(j,i-1) - nu(j-1,i)

                end do
            end do

            if (bd_ipd) then
                do j = 2,jmax
                    delta(j,1) = delta(j,1) - mu(j,imax) - nu(j-1,1)
                end do
            end if

            if (bd_jpd) then
                do i = 2,imax
                    delta(1,i) = delta(1,i) - mu(1,i-1) - nu(jmax,i)
                end do
            end if

            if (bd_ipd.AND.bd_jpd) then
                delta(1,1) = delta(1,1) - mu(1,imax) - nu(jmax,1)
            end if

            where (abs(delta) < 1e-3)
                delta = 1
            end where 

            phi = phi_zero


            if (.NOT.(bd_jpd)) then
                if(.NOT.(kbd_jminnm))  phi(1,:) = bd_jminvals
                if(.NOT.(kbd_jmaxnm)) phi(jmax,:) = bd_jmaxvals
            end if

            if (.NOT.(bd_ipd)) then
                if(.NOT.(kbd_iminnm)) phi(:,1) = bd_iminvals
                if(.NOT.(kbd_imaxnm)) phi(:,imax) = bd_imaxvals
            end if

            sum_f = sum(abs(f_mod))
            print *, "SUM_F: ", sum_f

            res_norm = thresh + 100

            print *, "Res_norm (before SOR) : ", res_norm
            niterations = 0

            do while ((res_norm > thresh).AND.(niterations<niterations_max))
                !
                res_norm = 0.0
                !
                ! Calculate residual for interior points
                do j = 2,jmax-1
                    do i = 2,imax-1

                        res = mu(j,i)*phi(j,i+1) + mu(j,i-1)*phi(j,i-1) &
                        + nu(j,i)*phi(j+1,i) + nu(j-1,i)*phi(j-1,i) &
                        + delta(j,i)*phi(j,i) - f_mod(j,i)

                        phi(j,i) = phi(j,i) - omega*res/delta(j,i)

                        res_norm = res_norm + abs(res)

                    end do
                end do

                !If periodic calculate residual for border values
                if(bd_jpd) then
                    do i = 2,imax-1
                        ! j=1 border
                        res = mu(1,i)*phi(1,i+1) + mu(1,i-1)*phi(1,i-1) &
                        + nu(1,i)*phi(2,i) + nu(jmax,i)*phi(jmax,i) &
                        + delta(1,i)*phi(1,i) - f_mod(1,i)

                        phi(1,i) = phi(1,i) - omega*res/delta(1,i)

                        res_norm = res_norm + abs(res)

                        ! j=jmax border
                        res = mu(jmax,i)*phi(jmax,i+1) + mu(jmax,i-1)*phi(jmax,i-1) &
                        + nu(jmax,i)*phi(1,i) + nu(jmax-1,i)*phi(jmax-1,i) &
                        + delta(jmax,i)*phi(jmax,i) - f_mod(jmax,i)

                        phi(jmax,i) = phi(jmax,i) - omega*res/delta(jmax,i)

                        res_norm = res_norm + abs(res)
                    end do
                end if

                if(bd_ipd) then
                    do j = 2,jmax-1
    
                        ! i=1 border
                        res = mu(j,1)*phi(j,2) + mu(j,imax)*phi(j,imax) &
                        + nu(j,1)*phi(j+1,1) + nu(j-1,1)*phi(j-1,1) &
                        + delta(j,1)*phi(j,1) - f_mod(j,1)

                        phi(j,1) = phi(j,1) - omega*res/delta(j,1)

                        res_norm = res_norm + abs(res)

                        ! i=imax border
                        res = mu(j,imax)*phi(j,1) + mu(j,imax-1)*phi(j,imax-1) &
                        + nu(j,imax)*phi(j+1,imax) + nu(j-1,imax)*phi(j-1,imax) &
                        + delta(j,imax)*phi(j,imax) - f_mod(j,imax)

                        phi(j,imax) = phi(j,imax) - omega*res/delta(j,imax)

                        res_norm = res_norm + abs(res)

    
                    end do

                end if

                if(bd_ipd.AND.bd_jpd) then
                    !Account for corners if both boundaries are periodic
                    !(j=1,i=1) corner
                    res = mu(1,1)*phi(1,2) + mu(1,imax)*phi(1,imax) &
                        + nu(1,1)*phi(2,1) + nu(jmax,1)*phi(jmax,1) &
                        + delta(1,1)*phi(1,1) - f_mod(1,1)

                        phi(1,1) = phi(1,1) - omega*res/delta(1,1)

                        res_norm = res_norm + abs(res)

                    !(j=1,i=imax) corner
                    res = mu(1,imax)*phi(1,1) + mu(1,imax-1)*phi(1,imax-1) &
                        + nu(1,imax)*phi(2,imax) + nu(jmax,imax)*phi(jmax,imax) &
                        + delta(1,imax)*phi(1,imax) - f_mod(1,imax)

                        phi(1,imax) = phi(1,imax) - omega*res/delta(1,imax)

                        res_norm = res_norm + abs(res)

                    !(j=jmax, i=1)
                    res = mu(jmax,1)*phi(jmax,2) + mu(jmax,imax)*phi(jmax,imax) &
                        + nu(jmax,1)*phi(1,1) + nu(jmax-1,1)*phi(jmax-1,1) &
                        + delta(jmax,1)*phi(jmax,1) - f_mod(jmax,1)

                        phi(jmax,1) = phi(jmax,1) - omega*res/delta(jmax,1)

                        res_norm = res_norm + abs(res)
                    
                    !(j=jmax, i=imax)
                    res = mu(jmax,imax)*phi(jmax,1) + mu(jmax,imax-1)*phi(jmax,imax-1) &
                        + nu(jmax,imax)*phi(1,imax) + nu(jmax-1,imax)*phi(jmax-1,imax) &
                        + delta(jmax,imax)*phi(jmax,imax) - f_mod(jmax,imax)

                        phi(jmax,imax) = phi(jmax,imax) - omega*res/delta(jmax,imax)

                        res_norm = res_norm + abs(res)
                end if
                !
                !If the absolute values of the potential are getting too large for single precision
                !subtract a constant amount from them to reset
                if ( correction_log) then
                    if(phi(jp,ip).GT.1e10) then
                        !Subtract the absolute value of an unmasked point
                        phi = phi - phi(jp,ip)
                    end if
                end if
                
                res_norm = res_norm/sum_f
                niterations = niterations + 1
                if ((MOD(niterations, 200000).EQ.0).OR.(niterations.EQ.1)) then
                    print *, "Iteration: ", niterations
                    print *, "Res_norm : ", res_norm
                end if

                !
            end do

            if (niterations >= niterations_max) then
                print *, "Number of maximum iterations exceeded... oh dear"
            else
                print *, "Completed in ", niterations, " iterations"
            end if 


            !Now the interior points of phi have been found, we set the boundary values

            if (kbd_iminnm) then
                phi(:,1) = phi(:,2) - e1u(:,1)*bd_iminvals
            else if (.NOT.bd_ipd) then
                phi(:,1) = bd_iminvals
            end if

            if (kbd_imaxnm) then
                phi(:,imax) = phi(:,imax-1) + e1u(:,imax-1)*bd_imaxvals
            else if (.NOT.bd_ipd) then
                phi(:,imax) = bd_imaxvals
            end if

            if (kbd_jminnm) then
                phi(1,:) = phi(2,:) - e2v(1,:)*bd_jminvals
            else if (.NOT.bd_jpd) then
                phi(1,:) = bd_jminvals
            end if

            if (kbd_jmaxnm) then
                phi(jmax,:) = phi(jmax-1,:) + e2v(jmax-1,:)*bd_jmaxvals
            else if (.NOT.bd_jpd) then
                phi(jmax,:) = bd_jmaxvals
            end if

        end subroutine sor_cgrid



        subroutine sor_while_fast(a, b, x_0, omega, thresh, sy, x )
            !
            real, intent(in) :: a(:,:)
            real, intent(in) :: b(SIZE(a,2))
            real, intent(in) :: x_0(size(a,2))
            real, intent(in) :: omega
            real, intent(in) :: thresh
            integer, intent(in) :: sy
            real :: res(size(a,2))
            real :: x_i(size(a,2))
            integer :: count, n, i
            real, intent(out) :: x(SIZE(a,2))
            real :: res_norm
            !
            n = size(x_0)

            DO i = 1,n

                res(i) = a(i,i)*x_0(i) - b(i)

                if (i .ne. 1) then
                    res(i) =  res(i) + a(i,i-1)*x_0(i-1)
                end if

                if (i .ne. n) then
                    res(i) = res(i) + a(i,i+1)*x_0(i+1)
                end if

                if ( i-sy >= 1 ) then
                    res(i) = res(i) + a(i,i-sy)*x_0(i-sy)
                end if

                if (i+sy <= n ) then
                    res(i) = res(i) + a(i,i+sy)*x_0(i+sy)
                end if

            END DO

            res_norm = sqrt(dot_product(res,res))/sqrt(dot_product(b,b))

            !print *, "INITIAL RES_NORM: ", res_norm

            x_i = x_0

            count = 0

            DO WHILE (res_norm > thresh)
                !
                !Calulate new iteration of x

                DO i = 1,n
                    x(i) = x_i(i) - (omega/a(i,i))*res(i)
                END DO

                !Calculate new residual from new x

                DO i = 1,n

                    res(i) = a(i,i)*x(i) - b(i)
    
                    if (i .ne. 1) then
                        res(i) = res(i) + a(i,i-1)*x(i-1)
                    end if
    
                    if (i .ne. n) then
                        res(i) = res(i) + a(i,i+1)*x(i+1)
                    end if
    
                    if ( i-sy >= 1 ) then
                        res(i) = res(i) + a(i,i-sy)*x(i-sy)
                    end if
    
                    if (i+sy <= n ) then
                        res(i) = res(i) + a(i,i+sy)*x(i+sy)
                    end if
                    !
                END DO
                !
                !
                res_norm = sqrt(dot_product(res,res))/sqrt(dot_product(b,b))
                x_i = x
                !
                count = count + 1
                !
                print *, "res_norm = ", res_norm
                !
            END DO
            !
            print *, "Complete after ", count, " iterations"
            
        end subroutine sor_while_fast



        subroutine sor_while(a, b, ld_inv, x_0, omega, thresh, x )
            !
            real, intent(in) :: a(:,:)
            real, intent(in) :: b(SIZE(a,2))
            real, intent(in) :: ld_inv(SIZE(a,1),SIZE(a,2))
            real, intent(in) :: x_0(size(a,2))
            real, intent(in) :: omega
            real, intent(in) :: thresh
            real :: res(size(a,2))
            real :: x_i(size(a,2))
            integer :: count
            real, intent(out) :: x(SIZE(a,2))
            real :: res_norm
            !
            res = matmul(a,x_0) - b
            res_norm = sqrt(dot_product(res,res))/sqrt(dot_product(b,b))

            x_i = x_0

            count = 0

            DO WHILE (res_norm > thresh)
                !
                x = x_i - omega*matmul(ld_inv,res)

                res = matmul(a,x) - b

                res_norm = sqrt(dot_product(res,res))/sqrt(dot_product(b,b))

                x_i = x
                !
                count = count + 1
                !
                print *, res_norm

            END DO
            !
            !print *, "Complete after ", count, " iterations"
        end subroutine sor_while

        subroutine neumann_mask_proc(a, mask, output)
            !
            ! Applies no normal gradient conditions to points in an array A
            ! near masked values.
            !
            !PROCEDURE
            !
            ! 1) Run through the elements of A
            ! 
            ! For each element...
            ! 2) identify masked neighbouring points (if any)
            ! 3) If masked neighbouring points are found, copy neighbouring unmasked values (if any) to set zero gradient
            !    The direction of copying is always from the lower value of index to the highest. For example
            !
            !    o - Unmasked value, x - Masked value
            !
            !                            o
            !                            |
            !                            v
            !                            o
            !    Case 1:       o --> o   x   o --> o  
            !                            o
            !                            |
            !                            v

            real, intent(in) :: a(:,:) !Array to be treated
            logical, intent(in) :: mask(size(a,1), size(a,2)) !Mask for array a (1 = masked, 0 = unmasked)
            real, intent(out)   :: output(size(a,1), size(a,2)) !Treated array (mask is unchanged)

            integer :: i, j !For loop indices
            integer :: imax, jmax
            logical :: mask_im1, mask_ip1, mask_jm1, mask_jp1, mask_ij

            imax = size(a,1)
            jmax = size(a,2)

            output = a

            do i = 1,imax
                do j = 1,jmax

                    mask_ij = mask(i,j)

                    if (mask_ij .eqv. .true.) THEN
                        ! If masked, set to zero and move to next element
                        output(i,j) = 0
                        CYCLE
                    end if

                    if (i .eq. 1) THEN 
                        mask_im1 = .true.
                    else 
                        mask_im1 = mask(i-1,j)
                    end if

                    if (i .eq. imax) THEN 
                        mask_ip1 = .true.
                    else 
                        mask_ip1 = mask(i+1,j)
                    end if

                    if (j .eq. 1) THEN 
                        mask_jm1 = .true.
                    else 
                        mask_jm1 = mask(i,j-1)
                    end if

                    if (j .eq. jmax) THEN 
                        mask_jp1 = .true.
                    else 
                        mask_jp1 = mask(i,j+1)
                    end if

                    !First resolve problem in x direction

                    if (( mask_im1 .eqv. .false. ).and.(mask_ip1 .eqv. .true.)) THEN
                        !Masked neighbour on RHS
                        output(i,j) = output(i-1,j)
                        !
                    else if (( mask_im1 .eqv. .true. ).and.(mask_ip1 .eqv. .false.)) THEN
                        !Masked neighbour on LHS
                        output(i+1,j) = output(i,j) 

                    end if

                    !Then resolve problem in y direction

                    if (( mask_jm1 .eqv. .false. ).and.(mask_jp1 .eqv. .true.)) THEN
                        !Masked neighbour above (larger j)
                        output(i,j) = output(i,j-1)
                        !
                    else if (( mask_im1 .eqv. .true. ).and.(mask_ip1 .eqv. .false.)) THEN
                        !Masked neighbour below (smaller j)
                        output(i,j+1) = output(i,j)
                        
                    end if


                end do
            end do

        end subroutine neumann_mask_proc

        subroutine neumann_mask_proc2(a, mask, output)
            !
            real, intent(in) :: a(:,:)
            logical, intent(in) :: mask(size(a,1), size(a,2))

            integer :: mark_array(size(a,1), size(a,2))
            integer :: mark_count

            integer :: i, j, imax, jmax, nmasked_marked, nmasked, im
            

            logical :: mask_im1, mask_ip1, mask_jm1, mask_jp1, mask_list(4)
            integer :: mark_im1, mark_ip1, mark_jm1, mark_jp1, mark_max, mark_list(4), mark_tmp


            real,intent(out) :: output(size(a,1), size(a,2))

            mark_array = 0

            mark_count = 1

            output = a

            imax = size(a,1)
            jmax = size(a,2)

            do i = 1, imax
                do j = 1, jmax

                    print * , "OUTPUT>>>>>>>>>>"
                    print *, output
                    print * , ">>>>>>>>>>>>>>>>>"

                    print * , "MARKS>>>>>>>>>>"
                    print *, mark_array
                    print * , ">>>>>>>>>>>>>>>>>"

                    !If element is masked, move on to next element
                    if (mask(i,j) .eqv. .true.) THEN
                        cycle
                    end if

                    if (i .eq. 1) then
                        mask_im1 = .true.
                        mark_im1 = 0
                    else
                        mask_im1 = mask(i-1,j)
                        mark_im1 = mark_array(i-1,j)
                    end if

                    if (i .eq. imax) then
                        mask_ip1 = .true.
                        mark_ip1 = 0
                    else
                        mask_ip1 = mask(i+1,j)
                        mark_ip1 = mark_array(i+1,j)
                    end if

                    if (j .eq. 1) then
                        mask_jm1 = .true.
                        mark_jm1 = 0
                    else
                        mask_jm1 = mask(i,j-1)
                        mark_jm1 = mark_array(i,j-1)
                    end if

                    if (j .eq. jmax) then
                        mask_jp1 = .true.
                        mark_jp1 = 0
                    else
                        mask_jp1 = mask(i,j+1)
                        mark_jp1 = mark_array(i,j+1)
                    end if
                    !
                    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    !
                    ! Count the number of masked and marked points
                    !
                    mask_list(1) = mask_im1
                    mask_list(2) = mask_ip1
                    mask_list(3) = mask_jm1
                    mask_list(4) = mask_jp1
                    !
                    mark_list(1) = mark_im1
                    mark_list(2) = mark_ip1
                    mark_list(3) = mark_jm1
                    mark_list(4) = mark_jp1
                    !
                    nmasked_marked = count((mask_list .and. (mark_list .ne. 0)))
                    nmasked = count(mask_list)

                    print *, "i,j", i,j
                    print *, "nmasked_marked: ", nmasked_marked

                    if (nmasked .eq. 0) then
                        cycle
                    
                    else if (nmasked_marked .eq. 0) then

                        if (mask_im1.and.(i .ne. 1)) then
                            output(i-1,j) = output(i,j)
                            mark_array(i-1,j) = mark_count
                        end if

                        if (mask_ip1.and.(i .ne. imax)) then
                            output(i+1,j) = output(i,j)
                            mark_array(i+1,j) = mark_count
                        end if

                        if (mask_jm1.and.(j .ne. 1)) then
                            output(i,j-1) = output(i,j)
                            mark_array(i,j-1) = mark_count
                        end if

                        if (mask_jp1.and.(j .ne. jmax)) then
                            output(i,j+1) = output(i,j)
                            mark_array(i,j+1) = mark_count
                        end if

                        mark_array(i,j) = mark_count

                        mark_count = mark_count + 1

                    else
                        !
                        mark_max = maxval(mark_list,mask = mask_list)
                        print *, "mark_max: ", mark_max
                        print *, "mark_list: ", mark_list
                        print *, "mask_list: ", mask_list
                        !
                        if (( mark_max .eq. mark_im1 ) .and. (mask_im1)) then
                            output(i,j) = output(i-1,j)
                        else if (( mark_max .eq. mark_ip1 ) .and. (mask_ip1)) then
                            output(i,j) = output(i+1,j)
                        else if (( mark_max .eq. mark_jm1 ) .and. (mask_jm1)) then
                            output(i,j) = output(i,j-1)
                        else if (( mark_max .eq. mark_jp1 ) .and. (mask_jp1)) then
                            output(i,j) = output(i,j+1)
                        end if

                        do im = 1,4

                            if ((mark_list(im).ne.mark_max).and.(mark_list(im).ne.0) .and. (mask_list(im))) then
                            !
                                where (mark_array .eq. mark_list(im) )
                                !
                                    output = output(i,j)
                                    mark_array = mark_max
                                !
                                end where
                            !
                            end if

                        end do


                    end if




                end do
            end do



            !
        end subroutine neumann_mask_proc2

end module SOR

