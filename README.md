# Eco3M-CarbOx
Biogeochemical Model

Eco3M is freely available under CeCILL license agreement (a French equivalent to the L-GPL license; 
http://cecill.info/licences/Licence_CeCILL_V1.1-US.html; last access: 10 February 2020). 
The Eco3M-CarbOx model is written in Fortran-90/95 and the plotting code is written in Matlab®. 
It is available on a Git Hub repository at https://github.com/klajaunie/Eco3M-CarbOx (last access: 10 February 2020).

To run Eco3M-CarbOx:
>> make !two executable will be created : eco3M_ini.exe and eco3M.exe
- the file config.ini allow to define:
    the time, time step, and save time of simulation
    variables 
    biogeochemical process 
- Results files are stocked in "SORTIES" directory
- Boundary conditions and forcings data are stocked in "DATA" directory
- All subroutines of biogeochemical processes are stocked in "F_PROCESS" directory 
