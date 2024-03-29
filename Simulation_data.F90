!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the wind tunnel problem with a step
!!
!! ARGUMENTS
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
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_pAmbient, sim_TAmbient, sim_rhoAmbient
  real, save :: sim_gamma, sim_smallP, sim_smallX, eos_singleSpeciesA, eos_singleSpeciesZ
  real, save :: sim_nozHeight, sim_nozExit, sim_nozIn, sim_chamHeight, sim_chamWidth
  real, save :: sim_pBacking
  real, save :: xmin, xmax, ymin, ymax

integer, save :: sim_meshMe
end module Simulation_data


