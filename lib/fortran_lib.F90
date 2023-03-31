module fortran_lib

    !Ths module contains various grid point calculations that are significnatly faster if coded in Fortran

    ! Remember to compile this code when any changes are made using the command
    ! f2py -c -m fortran_lib fortran_lib.F90 ; mv fortran_lib.*.so fortran_lib.so

    implicit none

    contains

    subroutine residual_overturning(v_array, rhop_array, rhop_coord, e3, e1, vmask, res_ov, rhop_depth, output_mask )
        !
        real, intent(in) :: v_array(:,:,:,:) ! Array of meridional velocities (t,z,y,x)
        !
        real, intent(in) :: rhop_array(SIZE(v_array,1),SIZE(v_array,2),SIZE(v_array,3),SIZE(v_array,4))
        !
        real, intent(in) :: rhop_coord(:)    !Density coordinates to interpolate to (rho)
        !
        real, intent(in) :: e3(SIZE(v_array,2),SIZE(v_array,3),SIZE(v_array,4)) ! Cell thicknesses for meridional velocities (z,y,x)
        !
        real, intent(in) :: e1(SIZE(v_array,3),SIZE(v_array,4)) ! Zonal cell widths for meridional velocities (y,x)
        !
        logical, intent(in) :: vmask(SIZE(v_array,2),SIZE(v_array,3),SIZE(v_array,4)) ! Mask for v points (True = Masked)
        !
        real, intent(out) :: res_ov(SIZE(v_array,1),SIZE(rhop_coord,1),SIZE(v_array,3),SIZE(v_array,4)) !Residual overturning stream function (t,rho,y,x)
        !
        real, intent(out) :: rhop_depth(SIZE(v_array,1),SIZE(rhop_coord,1),SIZE(v_array,3),SIZE(v_array,4)) !Depth of isopycnals
        !
        logical, intent(out) :: output_mask(SIZE(v_array,1),SIZE(rhop_coord,1),SIZE(v_array,3),SIZE(v_array,4)) !Mask for both output arrays (True = Masked)
        !
        integer :: ii, ij, ik, ip, it, ni, nj, nk, nt, np, nk_avail, ip_min, kind, ik_min, iip
        real :: v_single, e3_single, rhop_single, irhop, drhop, rhop_close, res_ov_ipm1, rhop_depth_ipm1
        real :: v_column(SIZE(v_array,2))
        logical :: vmask_column(SIZE(v_array,2))
        real :: rhop_column(SIZE(v_array,2))
        real :: e3_column(SIZE(v_array,2))
        real :: rhop_bounds(SIZE(v_array,2)+1)
        real :: v_zint(SIZE(v_array,1),SIZE(v_array,3),SIZE(v_array,4))
        logical :: first_bath_rho_log,ext_rho_log

        ni = SIZE(v_array,4)
        nj = SIZE(v_array,3)
        nk = SIZE(v_array,2)
        nt = SIZE(v_array,1)
        np = SIZE(rhop_coord,1)

        output_mask (:,:,:,:) = .false.
        rhop_depth(:,:,:,:) = 0.
        res_ov(:,:,:,:) = 0.
        v_zint(:,:,:) = 0.

        do ii = 1,ni
            do ij = 1,nj
                !
                !Look at the geometry of the column and determine if there are more than one data points
                vmask_column = vmask(:,ij,ii)

                !Count the number of unmasked points
                nk_avail = COUNT(.NOT.vmask_column)

                select case(nk_avail)
                    !
                    case(0) !Land column - No unmasked points
                        output_mask(:,:,ij,ii) = .true.
                        res_ov(:,:,ij,ii) = 0.
                        rhop_depth(:,:,ij,ii) = 0.
                    !
                    case(1) !Single point column - One unmasked point
                        output_mask(:,:,ij,ii) = .true.
                        res_ov(:,:,ij,ii) = 0.
                        rhop_depth(:,:,ij,ii) = 0.
                        !
                        !Find the location of the single unmasked point
                        do ik = 1,nk
                            if (.NOT.vmask_column(ik)) then
                                kind = ik
                                exit
                            end if
                        end do
                        !
                        e3_single = e3(kind, ij, ii)
                        !
                        do it = 1,nt
                            v_single = v_array(it,kind,ij,ii)
                            rhop_single = rhop_array(it,kind,ij,ii)
                            !
                            ip_min = MINLOC( ABS(rhop_coord - rhop_single), 1 ) !Find the closest isopycnal coordinate
                            rhop_close = rhop_coord(ip_min)
                            !
                            output_mask(it,ip_min,ij,ii) = .false.
                            res_ov(it,ip_min,ij,ii) = v_single
                            rhop_depth(it,ip_min,ij,ii) = e3_single / 2

                        end do
                    
                    case default !More than one unmasked point in the column
                        !
                        e3_column = e3(:,ij,ii)
                        !
                        do it = 1,nt
                            v_column = v_array(it,:,ij,ii)
                            rhop_column = rhop_array(it,:,ij,ii)
                            !
                            !Estimate the values of density between the v points
                            rhop_bounds(:) = 0.
                            !
                            !Extrapolate at the top to find the surface density
                            drhop = (rhop_column(2)-rhop_column(1))*e3_column(1)/(e3_column(2)+e3_column(1))
                            rhop_bounds(1) = rhop_column(1) - drhop
                            !
                            !Extrapolate at the bottom to find the bottom density
                            drhop = (rhop_column(nk_avail)-rhop_column(nk_avail-1))                                   &
                                    *e3_column(nk_avail)/(e3_column(nk_avail)+e3_column(nk_avail-1))
                            rhop_bounds(nk_avail+1) = rhop_column(nk_avail) + drhop
                            !
                            do ik = 2,nk_avail
                                !
                                drhop = (rhop_column(ik)-rhop_column(ik-1))                                           &
                                    * e3_column(ik-1)/(e3_column(ik-1)+e3_column(ik))
                                rhop_bounds(ik) = rhop_column(ik-1) + drhop
                                !
                            end do
                            !
                            !
                            ! Set temporary counting variables to zero before searching for density surfaces
                            ik_min = 1
                            res_ov_ipm1 = 0.
                            rhop_depth_ipm1 = 0.
                            first_bath_rho_log = .false.
                            !
                            do ip = 1,np
                                !
                                irhop = rhop_coord(ip)
                                !
                                !Determine if isopycnal lies outside the discrete domain below the previous isopycnal
                                ext_rho_log = ( (irhop < rhop_bounds(ik_min)))!.OR.(irhop > rhop_bounds(nk_avail+1)))
                                !

                                
                                ! Determine if isopycnal is the first one that lies in the bathymetry
                                ! Extrapolation needed in this case
                                ! if (ip /= 1) then
                                !    first_bath_rho_log = ((rhop_coord(ip-1) >= rhop_bounds(1))   &
                                !                         .AND.(rhop_coord(ip-1) <= rhop_bounds(nk_avail+1)))
                                !    first_bath_rho_log = (first_bath_rho_log.AND.(irhop > rhop_bounds(nk_avail+1)))
                                ! else
                                !    first_bath_rho_log = .false.
                                ! end if
                               
                                ! If this isopycnal is the first to hit the sea floor, set to full depth integral
                                ! and mask all values for higher density isopycnals
                                ! if (ext_rho_log.AND.first_bath_rho_log) then
                                !     res_ov(it,ip,ij,ii) = SUM(v_column*e3_column, MASK=.NOT.vmask_column)
                                !     rhop_depth(it,ip,ij,ii) = SUM(e3_column, MASK=.NOT.vmask_column)
                                !     if (ip /= np) then
                                !         do iip = ip+1, np
                                !             res_ov(it,iip, ij, ii) = 0.
                                !             rhop_depth(it,iip, ij, ii) = 0.
                                !             output_mask(it,iip, ij, ii) = .true.
                                !         end do
                                !     end if

                                !     exit
                                
                                ! end if

                                ! else

                                ! if (ext_rho_log.AND.(.NOT.first_bath_rho_log) ) then
                                if (ext_rho_log.OR.first_bath_rho_log) then
                                    output_mask(it,ip,ij,ii) = .true.
                                    res_ov(it,ip,ij,ii) = 0.
                                    rhop_depth(it,ip,ij,ii) = 0.

                                !if ((irhop < rhop_bounds(1)).OR.(irhop > rhop_bounds(nk_avail+1))) then
                                !    output_mask(it,ip,ij,ii) = .true.
                                !    res_ov(it,ip,ij,ii) = 0.
                                !    rhop_depth(it,ip,ij,ii) = 0.
                                
                                else
                                    !
                                    do ik = ik_min,nk_avail
                                        
                                        if ( ( irhop >= rhop_bounds(ik) ).AND.(irhop < rhop_bounds(ik+1)) ) then
                                            rhop_depth(it,ip,ij,ii) = rhop_depth_ipm1                   &
                                            + e3_column(ik) * (irhop - rhop_bounds(ik))                        &
                                            /(rhop_bounds(ik+1)-rhop_bounds(ik))

                                            !If density surface is not deeper than previous surface, mask and ignore
                                            if ( (ip > 1).AND.(rhop_depth(it,ip,ij,ii) <  rhop_depth(it,ip-1,ij,ii)) ) then
                                                output_mask(it,ip,ij,ii) = .true.
                                                res_ov(it,ip,ij,ii) = 0.

                                            !Otherwise calculate final part of residual overturning
                                            else
                                                res_ov(it,ip,ij,ii) = res_ov_ipm1                          &
                                                + v_column(ik)*e3_column(ik)*(irhop - rhop_bounds(ik))             &
                                                /(rhop_bounds(ik+1)-rhop_bounds(ik))

                                            end if

                                            exit
                                        else

                                            res_ov_ipm1 = res_ov_ipm1 + v_column(ik)*e3_column(ik)
                                            res_ov(it,ip,ij,ii) = res_ov_ipm1

                                            rhop_depth_ipm1 = rhop_depth_ipm1 + e3_column(ik)
                                            rhop_depth(it,ip,ij,ii) = rhop_depth_ipm1
                                            ik_min = ik_min + 1

                                            if (ik == nk_avail) then
                                                first_bath_rho_log = .true.
                                            end if


                                        end if


                                        
                                    end do
                                    !
                                end if


                                !
                            end do
                            !
                        end do

                end select

            end do
        end do

        !Convert to an integral from the bottom by subtracting the bottom value of res_ov
        ! do it = 1,nt
        !     do ik = 1,nk
        !         where(.NOT.vmask(ik,:,:))
        !             v_zint(it,:,:) = v_zint(it,:,:) + v_array(it,ik,:,:) * e3(ik,:,:)
        !         end where
        !     end do
        ! end do
        
        ! do it = 1,nt
        !     res_ov(it,np,:,:) = SUM(v_array(it,:,:,:)*e3(:,:,:), DIM=1, MASK=.NOT.vmask)
        !     rhop_depth(it,np,:,:) = SUM(e3, DIM=1, MASK=.NOT.vmask)
        !     output_mask(it,np,:,:) = ALL(vmask, DIM=1)
        ! end do
        

        do it = 1,nt
            do ip = 1,np
                res_ov(it,ip,:,:) = res_ov(it,ip,:,:)*e1(:,:)  !Multiply by cell width to calculate volume flux
            end do
        end do

    end subroutine residual_overturning


    ! VVVVV z2rhop subroutine is currently not used VVVVV !

    subroutine z2rhop( array, rhop_array, mask, rhop_coord, boundary_mode, output_array, output_mask )
        !Transforma z-coordinate array into density coordinates using linear interpolation
        !Note that this requires the density to NOT increase with depth at any point

        !Array containing the data we want to interpolate (t,z,y,x)
        real, intent(in) :: array(:,:,:,:)

        !Array containing the density information (t,z,y,x)
        real, intent(in) :: rhop_array(SIZE(array,1),SIZE(array,2),SIZE(array,3),SIZE(array,4))


        !mask for array and rhop_array (1=unmasked, 0=masked) !! NEMO convention !!
        logical, intent(in) :: mask(SIZE(array,2),SIZE(array,3),SIZE(array,4)) 


        !Density levels to interpolate onto (rho)
        real, intent(in) :: rhop_coord(:)

        !Treatment of values beyond the density range in the column
        character(len=20), intent(in) :: boundary_mode
        ! boundary_mode = 'mask' --> Mask values that lie beyond the range of densities in the column
        !               = 'extend' --> Values that lie beyond the range of densities in the column are equal to the first/last value in the range
        !               = 'extrapolate --> Values that lie beyond the range of densities in the column are equal to the extrapolated value.
        !

        !Array interpolated onto density coordinates (t,rho,y,x)
        real,intent(out) :: output_array(SIZE(array,1),SIZE(rhop_coord,1),SIZE(array,3),SIZE(array,4))

        !Mask for output_array (1=masked, 0=unmasked) !! Python convention !! (t,rho,y,x)
        logical,intent(out) :: output_mask(SIZE(array,1),SIZE(rhop_coord,1),SIZE(array,3),SIZE(array,4))

        integer :: ii, ij, it, ip, ik, ni, nj, nt, nk, nk_avail, np, kmin, kind, pmin, k1, k2, ind
        real :: column(SIZE(array,2))
        logical :: column_mask(SIZE(array,2))
        real :: rhop_column(SIZE(array,2))
        real :: single_rhop, single_point, rhop1, rhop2, point1, point2, irhop, rhop_column_max, rhop_column_min

        ni = SIZE(array,4)
        nj = SIZE(array,3)
        nk = SIZE(array,2)
        nt = SIZE(array,1)
        np = SIZE(rhop_coord)

        output_array(:,:,:,:) = 0.
        output_mask(:,:,:,:) = .FALSE.

        !Operate on each column
        do ii = 1,ni
            do ij = 1,nj
                do it = 1,nt

                    column_mask = mask(:,ij,ii)
                    column = array(it,:,ij,ii)
                    rhop_column = rhop_array(it,:,ij,ii)
                    !
                    !Filter the columns for any decreases in density with depth. >>>>
                    ind = 1
                    do ik = 2,nk
                        if ((rhop_column(ik) > rhop_column(ind)).AND.(column_mask(ik))) then
                            rhop_column(ind + 1) = rhop_column(ik)   
                            column(ind + 1) = column(ik)
                            column_mask(ind + 1) = column_mask(ik)

                            ind = ind + 1

                        end if
                    end do
                    !
                    do ik = ind + 1, nk
                        rhop_column(ik) = 0.
                        column(ik) = 0.
                        column_mask(ik) = .false.
                    end do

                    nk_avail = COUNT(column_mask)

                                            
                    if (nk_avail == 0) THEN
                        !
                        !If all points are masked, set all output values to zero and mask
                        output_array(:,:,ij,ii) = 0.
                        output_mask(:,:,ij,ii) = .TRUE.

                    else if (nk_avail == 1) THEN
                        !
                        kind = 0
                        do ik = 1,nk
                            if (column_mask(ik)) THEN 
                                kind = ik
                                exit
                            end if
                        end do

                        single_rhop = rhop_array(it,kind,ij,ii) 
                        single_point = array(it, kind, ij, ii)

                        !Find the nearest density level
                        pmin = MINLOC(ABS(rhop_coord-single_rhop), 1)

                        !Single data value is assigned to the nearest density level
                        output_array(it, pmin, ij, ii) = single_point

                        !All other density levels are masked 
                        output_mask(it,:,ij, ii) = .true.
                        output_mask(it, pmin, ij, ii) = .false.
                    
                    else

                        do ip = 1,np
                            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            irhop = rhop_coord(ip)
                            !
                            !
                            rhop_column_max = maxval(rhop_column, mask=column_mask)
                            rhop_column_min = minval(rhop_column, mask=column_mask)

                            if ( (irhop < rhop_column_max).AND.(irhop > rhop_column_min) ) THEN

                                kmin = MINLOC(ABS(irhop-rhop_column), 1, MASK=column_mask)
                                !
                                if (irhop - rhop_column(kmin) > 0) THEN
                                    k1 = kmin
                                    k2 = kmin + 1
                                else
                                    k1 = kmin - 1
                                    k2 = kmin
                                end if
                                !
                                !If k1 and k2 lie within the model domain we do not need to worry about boundaries
                                ! if ((k1 >= 1).and.(k2<=nk_avail)) THENls 
                                rhop1 = rhop_column(k1)
                                rhop2 = rhop_column(k2)

                                point1 = column(k1)
                                point2 = column(k2)
                                
                                output_array(it,ip,ij,ii) = point1 + (point2-point1)*(irhop - rhop1)/(rhop2-rhop1)

                            else 

                                !Values near the top and bottom boundaries depend on the boundary_mode
                                select case (boundary_mode)
                                    !
                                    case("mask")
                                        !
                                        output_array(it,ip,ij,ii) = 0.
                                        output_mask(it,ip,ij,ii) = .TRUE.
                                        !
                                    case("extend")
                                        !
                                        !Extend 
                                        if ( irhop <= rhop_column_min  ) THEN
                                            output_array(it,ip,ij,ii) = array(it,1,ij,ii)
                                        else
                                            output_array(it,ip,ij,ii) = array(it,nk_avail,ij,ii)
                                        end if
                                        !
                                    case("extrapolate")
                                        !
                                        if ( irhop <= rhop_column_min  ) THEN
                                            point1 = column(1)
                                            point2 = column(2)
                                            rhop1 = rhop_column(1)
                                            rhop2 = rhop_column(2)
                                        else
                                            point1 = column(nk_avail - 1)
                                            point2 = column(nk_avail)
                                            rhop1 = rhop_column(nk_avail-1)
                                            rhop2 = rhop_column(nk_avail)
                                        end if
                                        !
                                        output_array(it,ip,ij,ii) = point1 + (point2 - point1)*(irhop-rhop1)/(rhop2-rhop1)
                                    !
                                end select

                            end if
                        end do
                    end if
                end do




            end do
        end do
    end subroutine z2rhop

end module fortran_lib

