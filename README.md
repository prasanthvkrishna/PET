# PET
Estimate Potential EvapoTranspiration from WRF-Noah Output

This Fortran 90 script estimates Penman (1948) potential evapotranspitation from wrfoutput.
Write the PET into a netcdf file
~Prasant valayamkunnath, Virginia Tech (pvk03@vt.edu)

To compile:
	Libraries, NETCDF OpenMP
        
	gfortran -o Penman -fopenmp penman_et.f90 -I$NETCDF/include -L$NETCDF/lib -lnetcdf -lnetcdff
To run,
	Extract data.tar.gz to current working directory
	./Penman

Questions?: pvk03@vt.edu
