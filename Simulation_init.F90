!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the wind tunnel with a step problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_tAmbient', sim_tAmbient)

  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('eos_singleSpeciesA', eos_singleSpeciesA)
  call RuntimeParameters_get('eos_singleSpeciesZ', eos_singleSpeciesZ)

  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallp', sim_smallP)

  call RuntimeParameters_get('sim_nozHeight', sim_nozHeight)
  call RuntimeParameters_get('sim_nozExit', sim_nozExit)
  call RuntimeParameters_get('sim_nozIn', sim_nozIn)
  call RuntimeParameters_get('sim_pBacking', sim_pBacking)
  call RuntimeParameters_get('sim_chamHeight', sim_chamHeight)
  call RuntimeParameters_get('sim_chamWidth', sim_chamWidth)

  call RuntimeParameters_get('xmin', xmin)
  call RuntimeParameters_get('xmax', xmax)
  call RuntimeParameters_get('ymin', ymin)
  call RuntimeParameters_get('ymax', ymax)


  call Logfile_stamp("initializing for windtunnel + step", 'run_init')
  write (*,*) "flash:  initializing for wind tunnel + step"

end subroutine Simulation_init
