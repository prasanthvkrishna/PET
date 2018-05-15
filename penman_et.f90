program Penman_et 
 ! Written By Prasanth Valayamkunnath, Virginia Tech
 ! This fortran program estimates Penman Potential Evapotranspiration from WRF-ARW output (NetCDF Format).
 ! Load nessasary Libraries NetCDF and OpenMP (shared memory parallel processing)
  use netcdf
  use omp_lib

 
  implicit none
 !-----------Variable Names ---------------------!
 ! nlats   --> dimention along latitude
 ! nlong   --> dimension along longitude
 ! ndims   --> number of dimension
 ! sigma   --> Stefan-Boltzmann constant=5.67051D-8 W/m^2K^4    
 ! swin    --> shortwave in (Wm-2)
 ! lwin    --> longwave in (Wm-2)
 ! albedo  --> albedo from MODIS monthly values
 ! emiss   --> surface emissivity, epsilon
 ! tskin   --> surface temperature (K)
 ! delta   --> slope of saturation vapor pressure curve (FAO 56)
 ! qv_2m   --> humidity ratio at 2m (kg/kg)
 ! v_2m    --> wind speed at 2m (m/s)
 ! E_A     --> penamn original Rome wind function (mm)
 ! ea      --> sataurated vapor pressure 
 ! es      --> actual vapor pressure (kPa)
 ! gam     --> gamma- psychrometric constant (kPa/degreeC)
 ! cp      --> specific heat capaicity of moist air (0.001013 MJ/kg/degreeC)
 ! press   --> surface pressure (kPa)
 ! lambda  --> latent heat of vaporization-temperature depended formula (MJ/kg)
 ! pet     --> penman pet (mm)

  interface
    subroutine list_files(path,wrfFileNames,NFiles)
      character(len=100), dimension(:),allocatable, intent(inout) :: wrfFileNames
      character(*), intent(in) :: path
      integer,intent(out) :: NFiles
    end subroutine list_files
  end interface  

 !Define variables and Arrays: All variables are 3D (time,lat,lon)
  Integer, parameter 		 :: nlats  =324 , nlons =485, ndims =1 ,nmonths=12
  integer, parameter             :: styear = 2001, enyear = 2015, stmn=6, enmn=8
  integer, parameter,dimension(3):: date   = (/30,31,31/) 
  real,parameter		 :: sigma = 5.67051D-8 !W/m^2K^4 
  real,parameter		 :: cp = 0.001013 
  integer 			 :: lon_dimid, lat_dimid, rec_dimid,i,NF
  real    			 :: lats(NLATS), lons(NLONS)
  real    			 :: delta, E_A, es, ea, gam,alb,Rnet
  integer 			 :: lon_varid, lat_varid
  real, dimension(1:nlons,1:nlats)  :: swin,lwin,emiss,tskin,&
                                    qv_2m,v_2m, press, gh,t_2m ,pet
  real, dimension(1:nlons,1:nlats,1:nmonths)  ::albedo !from 12 month data
  character(len=100), dimension(:), allocatable :: wrfFileNames
  character(*), parameter :: filepath = "./"
  integer :: ncid, dimid, varid, is_err,ii,jj,num_threads
  


  ! Get all wrfoutput files in the present working directory
  call list_files(filepath,wrfFileNames,NF)

  ! Read MODIS Monthly albedo
  ! Open NetCDF file
  print*,"Program is reading MODIS ALBEDO"
  is_err = nf90_open("geo_em.d01.nc", NF90_NOWRITE, ncid)
  if (is_err /= 0) then
	write(*,'("Penman PET:  Problem opening geo_em.d01.nc file: ''", A, "''")')
	stop
  endif
  is_err = nf90_inq_varid(ncid,  "ALBEDO12M",  varid)
  is_err = nf90_get_var(ncid, varid, values=albedo, start=(/1,1,1,1/), count=(/nlons,nlats,nmonths,1/)) ! dimensions extracted from netcdf file
  if (is_err /= 0) then
        print*, 'ncid = ', ncid, "Error Reading Variable albedo in input NetCDF file."
  else
        print*,"Albedo loaded"
  endif
  !print*,albedo(106,106,1)
	

  do i = 1,NF ! loop through files

	      !Open NetCDF file
	      print*,"Program is reading ",wrfFileNames(i)
	      is_err = nf90_open(wrfFileNames(i), NF90_NOWRITE, ncid)
	      !print*,is_err," ",ncid


	      if (is_err /= 0) then
	         write(*,'("Penman PET:  Problem opening NetCDF file: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      ! Read all variables from NetCDF files

	      !--------- swin----------------!
	      call rd_2d_var("SW", ncid, swin, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  SW not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      !--------- lwin----------------!
	      call rd_2d_var("LW", ncid, lwin, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  LW not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      !--------- emiss----------------!
	      call rd_2d_var("EMI", ncid, emiss, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  EMI not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      !--------- tskin----------------!
	      call rd_2d_var("TSK", ncid, tskin, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  TSK not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      !--------- read qv_2m----------------!
	      call rd_2d_var("Q2", ncid, qv_2m, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  Q2 not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif
              
	      !--------- read v_2m----------------!
	      call rd_2d_var("WS", ncid, v_2m, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  WS not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      !--------- read press----------------!
	      call rd_2d_var("PSFC", ncid, press, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  PSFC not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      !--------- read gh----------------!
	      call rd_2d_var("GH", ncid, gh, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  GH not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif

	      !--------- read t_2m----------------!
	      call rd_2d_var("T2", ncid, t_2m, 1, nlats, 1, nlons, is_err)
	      
	      if (is_err /= 0) then
	         write(*,'("Penman PET:  T2 not found: ''", A, "''")') trim(wrfFileNames(i))
	         stop
	      endif
              !print*,gh(100,100),tskin(100,100)

	      
	      !~~~~~~~~ Loops i,j~~~~~~! Estimates PET
		  num_threads=1
		  !$ call omp_set_num_threads(num_threads)
		  !$omp parallel do
	      do ii=1,nlons
		 do jj=1,nlats
			! Suuumer average albedo
			alb  = albedo(ii,jj,6)+albedo(ii,jj,7)+albedo(ii,jj,8) ! percentage should be multiply with 0.01
			! Net Radiation similat to Noah LSM physics (1-alb)swin+lwin-emiss*sigma*tskin**4
		        Rnet = ((1-alb*0.01)*swin(ii,jj))+lwin(ii,jj)-(emiss(ii,jj)*sigma*(tskin(ii,jj)**4))
			delta= 4098.0*( 0.6108*exp( 17.27*( t_2m(ii,jj)-273.15 ) /( t_2m(ii,jj)-35.85 )))/(t_2m(ii,jj)-35.85)**2
			gam  = (0.665D-3)*press(ii,jj)/1000
			es   = 0.6108*exp( 17.27*( t_2m(ii,jj)-273.15 ) /( t_2m(ii,jj)-35.85 )) !kPa
			ea   = qv_2m(ii,jj)*press(ii,jj)*0.001/(0.622*qv_2m(ii,jj))/100
                        !print*,ea
  			ea   = max (ea, 0.001) !kPa
			E_A  = 0.26*(1+0.54*v_2m(ii,jj)*(es-ea))/24 !mm h-1

                        pet(ii,jj)=(delta*((Rnet-gh(ii,jj))*3600D-6/(2.501-(0.002361*(t_2m(ii,jj)-273.15)))))+(gam*E_A)/(delta+gam)
                        !print*,pet(ii,jj)
                        
		 end do
              end do
			  !$omp end parallel do
	      !~~~~~~~~~~~~~~~Add PET to output netcdf and write a netcdf file~~~~~~~~~!
	     
              call wrt_pet_var(pet,wrfFileNames(i),nlons,nlats)
              
  end do

end program Penman_et



!----------------------------------------------------------------------------------
  subroutine rd_2d_var(vname, ncid, d_array, st_lat, en_lat, st_lon, en_lon, iserr)
    use netcdf    
    character(len=*), intent(in) :: vname
    integer, intent(in)  :: ncid
    integer, intent(in)  :: st_lat, en_lat, st_lon, en_lon
    real, dimension(st_lon:en_lon,st_lat:en_lat), intent(out) :: d_array
    integer              :: ierr, varid
    integer, intent(out) :: iserr

    ierr = nf90_inq_varid(ncid,  vname,  varid)
    !print*,varid
    if (ierr /= 0) then
          print*, 'ncid = ', ncid, "Error Reading Variable "//trim(vname)//" in input NetCDF file."
    else
          iserr = ierr
    endif

    ierr = nf90_get_var(ncid, varid, values=d_array, start=(/st_lon,st_lat,1/), count=(/en_lon-st_lon+1,en_lat-st_lat+1,1/))
    if (ierr /= 0) then
          print*, 'ncid = ', ncid, "Error Reading Variable "//trim(vname)//" in input NetCDF file."
    else
          iserr = ierr
          return
    endif

  end subroutine rd_2d_var
!----------------------------------------------------------------------------------
subroutine list_files(path,wrfFileNames,NFiles)
   character(len=100), dimension(:),allocatable, intent(inout) :: wrfFileNames
   character(*), intent(in) :: path
   integer :: cnt,iserr,cntwrf
   integer,intent(out) :: NFiles
   real    :: irer
   ! get the files
   call system("ls "//trim(path)//"OBJ1_* >filelist.txt")
   open(41,FILE='filelist.txt',action="read")
   !counts of file
   cnt = 0
   do
     read(41,FMT='(a)',iostat=iserr) irer
     if (iserr/=0) EXIT
     cnt = cnt+1
   end do
   NFiles = cnt
   print*, "Number of NetCDF files: " , NFiles
   allocate(wrfFileNames(NFiles))
   rewind(41)
   do cntwrf = 1,NFiles
    read(41,'(a)') wrfFileNames(cntwrf)
   !print*,wrfFileNames(cntwrf)
   end do
   close(41)
   return
end subroutine list_files

!---------------------------------------------------------------------------------------

 subroutine wrt_pet_var(pet_out,outfname,nlons,nlats)
  use netcdf 
  character (len = *), intent(in) ::outfname   
  character (len = *), parameter :: FILE_NAME = "PET_"
  integer, intent(in):: nlons,nlats
  integer, parameter :: nrecs=1,ndims=3

  integer    :: lats(nlats), lons(nlons)
  integer :: lon_varid  , lat_varid

  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: REC_NAME = "time"
  integer                        :: lon_dimid, lat_dimid, rec_dimid
  integer, dimension(3)          :: start, counts
  character (len = *), parameter :: PET_NAME="Penman_PET"
  integer :: pet_varid
  integer :: dimids(2)!dimids(3)

  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PET_UNITS = "mmh-1"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  real,dimension(1:nlons, 1:nlats), intent(in):: pet_out
  integer, parameter :: st_lat = 1, st_lon = 1
  integer :: lat, lon, rec, i

  do lat = 1, nlats
     lats(lat) = lat
  end do
  do lon = 1, nlons
     lons(lon) = lon
  end do
  !print*,FILE_NAME//outfname(16:)

  call iserror( nf90_create(FILE_NAME//outfname(16:), nf90_clobber, ncid) )

  call iserror( nf90_def_dim(ncid, LAT_NAME, nlats, lat_dimid) )
  call iserror( nf90_def_dim(ncid, LON_NAME, nlons, lon_dimid) )
  !call iserror( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

  call iserror( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
  call iserror( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )

  call iserror( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
  call iserror( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )

  dimids = (/ lon_dimid, lat_dimid /)

  call iserror( nf90_def_var(ncid, PET_NAME, NF90_REAL, dimids, pet_varid) )

  call iserror( nf90_put_att(ncid, pet_varid, UNITS, PET_UNITS) )

  call iserror( nf90_enddef(ncid) )

  call iserror( nf90_put_var(ncid, lat_varid, lats) )
  call iserror( nf90_put_var(ncid, lon_varid, lons) )
  counts = (/ nlons, nlats, 1 /)
  start = (/ 1, 1, 1 /)

  !do rec = 1, nrecs
     !start(3) = rec
     call iserror( nf90_put_var(ncid, pet_varid, pet_out, start = (/1,1/)) )
  !end do


  contains
    subroutine iserror(status)
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
    end subroutine iserror  
end subroutine wrt_pet_var














 
