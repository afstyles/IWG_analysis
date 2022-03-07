module fortran_lib

    !Ths module contains various grid point calculations that are significnatly faster if coded in Fortran

    implicit none

    contains

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

        integer :: ii, ij, it, ip, ik, ni, nj, nt, nk, nk_avail, np, kmin, kind, pmin, k1, k2
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
                column_mask = mask(:,ij,ii)
                nk_avail = COUNT(column_mask)
                                         
                if (nk_avail == 0) THEN
                    !
                    !If all points are masked, set all output values to zero and mask
                    output_array(:,:,ij,ii) = 0.
                    output_mask(:,:,ij,ii) = .TRUE.
                    continue

                else if (nk_avail == 1) THEN
                    !
                    do it = 1,nt

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
                    end do
                    !
                    continue
                end if


                do it = 1,nt
                    do ip = 1,np
                        !
                        column = array(it,:,ij,ii)
                        rhop_column = rhop_array(it,:,ij,ii)
                        irhop = rhop_coord(ip)
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
                end do




            end do
        end do
    end subroutine z2rhop

end module fortran_lib