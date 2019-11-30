add_subdirectory ( InitialStates )
add_subdirectory ( MolecularState )
add_subdirectory ( ODE )
add_subdirectory ( Species )
add_subdirectory ( Trajectories )

add_sources (
  Collision_Class.F90
  ImpactParameter_Class.F90
  Integrator_Class.F90
  Matrix_Inversion_Module.F90
  Transformation_Class.F90
)
