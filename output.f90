module output_netcdf
USE netcdf
use shallow_water
integer :: ncid
integer :: x_id,xu_id,y_id,yv_id,t_id
integer :: eta_id,u_id,v_id
Integer :: ierr
integer :: nb

logical :: cg_netcdf = .false.

contains

subroutine change_output(option)

logical :: option
cg_netcdf = option

end subroutine change_output


subroutine CreatenetCDF


! creating a file named « output.nc »
! ncid is the identifier of the NetCDF file (to be used subsequently in all NetCDF calls)

if (cg_netcdf) then
  ierr=nf90_create("cg_output.nc",NF90_CLOBBER,ncid)
else
  ierr=nf90_create("output.nc",NF90_CLOBBER,ncid)
end if

! To create NetCDF dimensions, provide ncid, dimension name, the corresponding size, and an identifying integer

ierr = nf90_def_dim(ncid,'x',curgrid%Nx+2,x_id)
ierr = nf90_def_dim(ncid,'xu',curgrid%Nx+1,xu_id)
ierr = nf90_def_dim(ncid,'y',curgrid%Ny+2,y_id)
ierr = nf90_def_dim(ncid,'yv',curgrid%Ny+1,yv_id)

! For time, specify that the dimension is not predefined -> NF90_UNLIMITED

ierr = nf90_def_dim(ncid,'t',NF90_UNLIMITED,t_id)

! Variable definitions: provide ncid, variable name, variable type, variable dimensions (+ time dimension)

ierr=nf90_def_var(ncid,'eta',NF90_DOUBLE,(/x_id,y_id,t_id/),eta_id)
ierr=nf90_def_var(ncid,'u',NF90_DOUBLE,(/xu_id,y_id,t_id/),u_id)
ierr=nf90_def_var(ncid,'v',NF90_DOUBLE,(/x_id,yv_id,t_id/),v_id)

! close the <<header>> of the NetCDF file
ierr=nf90_enddef(ncid)

nb = 0
end subroutine

subroutine OutPutNetCDF


! writing data values

nb = nb + 1
ierr=nf90_put_var(ncid,eta_id,curgrid%eta(0:curgrid%Nx+1,0:curgrid%Ny+1),start=(/1,1,nb/))
ierr=nf90_put_var(ncid,u_id,curgrid%u(1:curgrid%Nx+1,0:curgrid%Ny+1),start=(/1,1,nb/))
ierr=nf90_put_var(ncid,v_id,curgrid%v(1:curgrid%Nx+1,1:curgrid%Ny+1),start=(/1,1,nb/))

ierr=nf90_sync(ncid)


end subroutine

subroutine CloseNetCDF

ierr=nf90_close(ncid)

end subroutine

end module output_netcdf
