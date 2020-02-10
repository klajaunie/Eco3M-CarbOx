# Eco3M-CarbOx
Biogeochemical Model

Eco3M is freely available under CeCILL license agreement (a French equivalent to the L-GPL license; 
http://cecill.info/licences/Licence_CeCILL_V1.1-US.html; last access: 10 February 2020). 
The Eco3M-CarbOx model is written in Fortran-90/95 and the plotting code is written in matlab. 
It is available on a Git Hub repository at https://github.com/klajaunie/Eco3M-CarbOx (last access: 10 February 2020).

To run Eco3M-CarbOx:
>> make !two executable will be created : eco3M_ini.exe and eco3M.exe
- the files config.ini allow to definie:
    the time, time step, and save time of simulation
    variables 
    biogeochemical process 
- Results files are stocked on "SORTIES" repertory
- Boundary conditions and forcings data are stocked on "DATA" repertory
- All subroutine of biogeochemical processes are stocked on "F_PROCESS" repertory 
