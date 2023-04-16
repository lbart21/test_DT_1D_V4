"""
Testing file for V4 of 1D Design Tool

Contents:
1. Regular Sod's Shock Tube With SinglePhaseMultiSpeciesNonReactivePipe Components
2. Sod's Shock Tube With SinglePhaseMultiSpeciesMultiInletNonReactivePipe Components
    With 1 Inlet On Each Inlet Face
3. Sod's Shock Tube With SinglePhaseMultiSpeciesMultiInletNonReactivePipe Components
    With 3 & 8 Inlets On The Inlet Faces
4. 
"""
import numpy as np

from copy import deepcopy

from gdtk.gas import GasModel, GasState, ThermochemicalReactor

from Algorithms.DT_1D_V4.fluid_models.flow_state import FlowState

from Algorithms.DT_1D_V4.extras.join_blocks import JointBlock

from Algorithms.DT_1D_V4.models.inlet_block_to_mimic_multi_inlet_flow.\
        inlet_block_to_mimic_multi_inlet_flow \
            import InletBlockToMimicMultiInletFlow

from Algorithms.DT_1D_V4.models.single_phase_multi_species_finite_rate_reactive_pipe.\
        single_phase_multi_species_finite_rate_reactive_pipe \
            import SinglePhaseMultiSpeciesFiniteRateReactivePipe

from Algorithms.DT_1D_V4.models.single_phase_multi_species_multi_inlet_nonreactive_pipe.\
        single_phase_multi_species_multi_inlet_nonreactive_pipe \
            import SinglePhaseMultiSpeciesMultiInletNonReactivePipe

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_conical_nozzle.\
        single_phase_multi_species_nonreactive_conical_nozzle \
            import SinglePhaseMultiSpeciesNonReactiveConicalNozzle

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_quadratic_nozzle.\
        single_phase_multi_species_nonreactive_quadratic_nozzle \
            import SinglePhaseMultiSpeciesNonReactiveQuadraticNozzle

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_elliptic_nozzle.\
        single_phase_multi_species_nonreactive_elliptic_nozzle \
            import SinglePhaseMultiSpeciesNonReactiveEllipticNozzle

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_mixing_pipe.\
        single_phase_multi_species_nonreactive_mixing_pipe \
            import SinglePhaseMultiSpeciesNonReactiveMixingPipe

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_mixing_elliptic_nozzle.\
        single_phase_multi_species_nonreactive_mixing_elliptic_nozzle \
            import SinglePhaseMultiSpeciesNonReactiveMixingEllipticNozzle

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_mixing_quadratic_nozzle.\
        single_phase_multi_species_nonreactive_mixing_quadratic_nozzle \
            import SinglePhaseMultiSpeciesNonReactiveMixingQuadraticNozzle

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_pipe_with_heat_flux.\
        single_phase_multi_species_nonreactive_pipe_with_heat_flux \
            import SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_pipe_with_heat_generation.\
        single_phase_multi_species_nonreactive_pipe_with_heat_generation \
            import SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen

from Algorithms.DT_1D_V4.models.single_phase_multi_species_nonreactive_pipe.\
        single_phase_multi_species_nonreactive_pipe \
            import SinglePhaseMultiSpeciesNonReactivePipe

from Algorithms.DT_1D_V4.models.\
        single_phase_multi_species_reactive_pipe_with_uniform_combustion_fraction.\
        single_phase_multi_species_reactive_pipe_with_uniform_combustion_fraction \
            import SinglePhaseMultiSpeciesReactivePipeWithUniformCombustionFraction

from Algorithms.DT_1D_V4.models.single_phase_multi_species_reactive_pipe_with_fast_chemistry.\
        single_phase_multi_species_reactive_pipe_with_fast_chemistry \
            import SinglePhaseMultiSpeciesReactivePipeWithUniformFastChemistry

from Algorithms.DT_1D_V4.models.single_phase_multi_species_reactive_pipe_with_partial_equilibrium_chemistry.\
        single_phase_multi_species_reactive_pipe_with_partial_equilibrium_chemistry \
            import SinglePhaseMultiSpeciesReactivePipeWithPartialEquilibriumChemistry

from Algorithms.DT_1D_V4.models.single_phase_multi_species_reactive_pipe_with_fast_chemistry_for_goal_massf.\
        single_phase_multi_species_reactive_pipe_with_fast_chemistry_for_goal_massf \
            import SinglePhaseMultiSpeciesReactivePipeWithFastChemistryForGoalMassf

from Algorithms.DT_1D_V4.models.single_phase_multi_species_reactive_pipe_with_tanh_fast_chemistry_for_goal_massf_profile.\
        single_phase_multi_species_reactive_pipe_with_tanh_fast_chemsitry_for_goal_massf_profile \
            import SinglePhaseMultiSpeciesReactivePipeWithTanhFastChemistryForGoalMassfProfile

from Algorithms.StoichiometricRatioSolver.determine_stoichiometric_coefficients \
        import SystemOfReactions
from Algorithms.StoichiometricRatioSolver.molecule import Molecule
from Algorithms.StoichiometricRatioSolver.reaction import Reaction

from Algorithms.DT_1D_V4.solver.solver import Solver

########################################## COMMON PARAMETERS #######################################
############################################## cfl_flag ############################################
#   [False, fixed_dt]
#   [True, "no_ramping", fixed_cfl]
#   [True, "cfl_ramp_from_tCurrent", [[cfl_1, t_1], [cfl_2, t_2], ..., [cfl_n, t_n]]]
#   [True, "cfl_ramp_from_step", [[cfl_1, step_1], [cfl_2, step_2], ..., [cfl_n, step_n]]]
#   [True, "dt_ramp_then_cfl_ramp_from_tCurrent", [dt_init, scale], [[cfl_1, t_1], [cfl_2, t_2], ... , [cfl_n, t_n]]]
#   [True, "dt_ramp_then_cfl_ramp_from_step", [dt_init, scale], [[cfl_1, step_1], [cfl_2, step_2], ... , [cfl_n, step_n]]]

########################################## boundary conditions #####################################
# [interface_id, ["WallNoSlip_BC", n_layers]]
# [interface_id, ["WallWithSlip_BC", n_layers]]
# [interface_id, ["SupersonicInFlow_BC", n_layers, flow_state]]
# [interface_id, ["FromStagnationInFlow_BC", n_layers, stag_fluid_state]]
# [interface_id, ["FromStagnationWithMassFlowRateInFlow_BC", n_layers, [stag_fluid_state, mass_flux]]]
# [interface_id, ["SimpleOutFlow_BC", n_layers]]
# [interface_id, ["SimpleExtrapolateOutFlow_BC", n_layers]]
# [interface_id, ["FixedPOutFlow_BC", n_layers, fluid_state]]
# [interface_id, ["FixedPTOutFlow_BC", n_layers, fluid_state]]

############################################# flux_scheme ##########################################
# "AUSMPlusUPOriginal"
# "AUSMPlusUPPaper"

############################################ recon_scheme ##########################################
# ["eilmer_L1R1", [1, 1]]
# ["eilmer_L2R2", [2, 2]]
# ["eilmer_L3R3", [3, 3]]
# ["eilmer_L2R1", [2, 1]]
# ["eilmer_L1R2", [1, 2]]
# ["MUSCL", [1, 1]] # [1, 1] is wrong, don't use yet
# ["THINC", [1, 1]] # [1, 1] is wrong, don't use yet
# ["MUSCLTHINCBVD", [1, 1]] # [1, 1] is wrong, don't use yet

############################################# limiter ##############################################
# "eilmer_vanAlbada_L3R3"
# "eilmer_vanAlbada_L2R2"
# "nonuniform_vanAlbada"
# "nonuniform_vanLeer"
# "nonuniform_minmod"
# "nonuniform_superbee"
# "nonuniform_MC"
# "uniform_Koren"
# "uniform_minmod"
# "uniform_MC"
# "uniform_Osher"
# "uniform_ospre"
# "uniform_superbee"
# "uniform_Sweby"
# "uniform_UMIST"
# "uniform_vanAlbada1"
# "uniform_vanAlbada2"
# "uniform_vanLeer"
# "uniform_Gregori"

###################### Quantities available for spatial distribution extraction ####################
# p, T, vel_x, rho, Ma, p_t, T_t, massf, molef, geometry data (automatically done)

####################### Quantities available for transient cell data extraction ####################
# p, T, vel_x, rho, Ma, p_t, T_t, massf, molef

#################### Quantites available for transient interface data extraction ###################
# p, mass_flux, momentum_flux, enthalpy_flux





####################################################################################################
######### Regular Sod's Shock Tube With SinglePhaseMultiSpeciesNonReactivePipe Components ##########
"""
# Left component is SinglePhaseMultiSpeciesNonReactivePipe
# Right component is SinglePhaseMultiSpeciesNonReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS, geometry = [D, L], \
            init_flow_state = init_fs_lft, recon_props = ["p", "T", "vel_x"], \
            recon_scheme = ["eilmer_L1R1", [1, 1]], limiter = "eilmer_vanAlbada_L2R2", update_from = "pT", \
            comp_label = "P1", flux_scheme = "AUSMPlusUPOriginal", \
            reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS, geometry = [D, L], \
            init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x"], \
            recon_scheme = ["eilmer_L1R1", [1, 1]], limiter = "eilmer_vanAlbada_L2R2", update_from = "pT", \
            comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
            reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 1, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
######## Sod's Shock Tube With SinglePhaseMultiSpeciesMultiInletNonReactivePipe Components #########
################################# With 1 Inlet On Each Inlet Face ##################################
"""
# Left component is SinglePhaseMultiSpeciesMultiInletNonReactivePipe
# Right component is reversed SinglePhaseMultiSpeciesMultiInletNonReactivePipe
# Both components have 1 inlet on the specified inlet boundary
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

A_inlet = 0.25 * np.pi * D ** 2.0

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesMultiInletNonReactivePipe(n_cells = NCELLS, n_inlets = 1, \
            bulk_geometry = [D, L], inlet_areas = [A_inlet], init_flow_state = init_fs_lft, \
            comp_label = "P1", limiter = "eilmer_vanAlbada", recon_scheme = ["copy", [1, 1]], \
            recon_props = ["p", "T", "vel_x"], update_from = "pT", \
            flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_mirrored_flow = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesMultiInletNonReactivePipe(n_cells = NCELLS, n_inlets = 1, \
            bulk_geometry = [D, L], inlet_areas = [A_inlet], init_flow_state = init_fs_rght, \
            comp_label = "P2", limiter = "eilmer_vanAlbada", recon_scheme = ["copy", [1, 1]], \
            recon_props = ["p", "T", "vel_x"], update_from = "pT", \
            flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_mirrored_flow = True)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 2, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
######## Sod's Shock Tube With SinglePhaseMultiSpeciesMultiInletNonReactivePipe Components #########
############################## With 3 & 8 Inlets On The Inlet Faces ################################
"""
# Left component is SinglePhaseMultiSpeciesMultiInletNonReactivePipe
# Right component is reversed SinglePhaseMultiSpeciesMultiInletNonReactivePipe
# Left component has 3 identical inlets
# Right component has 8 identical inlets
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

N_INLETS_LFT = 3
N_INLETS_RGHT = 8

A_inlet = 0.25 * np.pi * D ** 2.0

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesMultiInletNonReactivePipe(n_cells = NCELLS, n_inlets = N_INLETS_LFT, \
            bulk_geometry = [D, L], \
            inlet_areas = [A_inlet / N_INLETS_LFT] * N_INLETS_LFT, init_flow_state = init_fs_lft, \
            comp_label = "P1", limiter = "eilmer_vanAlbada", recon_scheme = ["copy", [1, 1]], \
            recon_props = ["p", "T", "vel_x"], update_from = "pT", \
            flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_mirrored_flow = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])
P1.add_boundary_conditions(BC = [1, ["SimpleOutFlow_BC", 2]])
P1.add_boundary_conditions(BC = [2, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesMultiInletNonReactivePipe(n_cells = NCELLS, n_inlets = N_INLETS_RGHT, \
            bulk_geometry = [D, L], inlet_areas = [A_inlet / N_INLETS_RGHT] * N_INLETS_RGHT, \
            init_flow_state = init_fs_rght, \
            comp_label = "P2", limiter = "eilmer_vanAlbada", recon_scheme = ["copy", [1, 1]], \
            recon_props = ["p", "T", "vel_x"], update_from = "pT", \
            flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_mirrored_flow = True)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])
P2.add_boundary_conditions(BC = [NCELLS + 1, ["SimpleOutFlow_BC", 2]])
P2.add_boundary_conditions(BC = [NCELLS + 2, ["SimpleOutFlow_BC", 2]])
P2.add_boundary_conditions(BC = [NCELLS + 3, ["SimpleOutFlow_BC", 2]])
P2.add_boundary_conditions(BC = [NCELLS + 4, ["SimpleOutFlow_BC", 2]])
P2.add_boundary_conditions(BC = [NCELLS + 5, ["SimpleOutFlow_BC", 2]])
P2.add_boundary_conditions(BC = [NCELLS + 6, ["SimpleOutFlow_BC", 2]])
P2.add_boundary_conditions(BC = [NCELLS + 7, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS + 2, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 3, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""     
####################################################################################################
######## Sod's Shock Tube With SinglePhaseMultiSpeciesNonReactiveConicalNozzle Components ##########
##################################### With Constant Diameter #######################################
"""
# Left component is SinglePhaseMultiSpeciesNonReactiveConicalNozzle
# Right component is SinglePhaseMultiSpeciesNonReactiveConicalNozzle
# Both components have uniform diameter
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesNonReactiveConicalNozzle(n_cells = NCELLS, geometry = [D, D, L], \
                init_flow_state = init_fs_lft, recon_props = ["p", "T", "vel_x"], \
                recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
                comp_label = "P1", flux_scheme = "AUSMPlusUPOriginal", \
                reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactiveConicalNozzle(n_cells = NCELLS, geometry = [D, D, L], \
                init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x"], \
                recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
                comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
                reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 4, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
####### Sod's Shock Tube With SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux Components ########
##################################### With Heat Flux = 0.0 #########################################
"""
# Left component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux
# Right component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux
# Both components have heat flux = 0.0
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux(n_cells = NCELLS, geometry = [D, L], \
            init_flow_state = init_fs_lft, recon_props = ["p", "T", "vel_x"], \
            recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
            comp_label = "P1", flux_scheme = "AUSMPlusUPOriginal", heat_flux = 0.0, \
            reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux(n_cells = NCELLS, geometry = [D, L], \
            init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x"], \
            recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
            comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", heat_flux = 0.0, \
            reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 5, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
####### Sod's Shock Tube With SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen Components #########
##################################### With Heat Gen = 0.0 ##########################################
"""
# Left component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen
# Right component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen
# Both components have heat gen = 0.0
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen(n_cells = NCELLS, geometry = [D, L], \
                init_flow_state = init_fs_lft, recon_props = ["p", "T", "vel_x"], \
                recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
                comp_label = "P1", flux_scheme = "AUSMPlusUPOriginal", heat_gen = 0.0, \
                reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen(n_cells = NCELLS, geometry = [D, L], \
                init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x"], \
                recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
                comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", heat_gen = 0.0, \
                reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 6, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
######### Sod's Shock Tube With SinglePhaseMultiSpeciesNonReactiveMixingPipe Components ############
"""
# Left component is SinglePhaseMultiSpeciesNonReactiveMixingPipe
# Right component is SinglePhaseMultiSpeciesNonReactiveMixingPipe
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
print(init_gm_lft.species_names)
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesNonReactiveMixingPipe(n_cells = NCELLS, geometry = [D, L], \
            init_flow_state = init_fs_lft, recon_props = ["p", "T", "vel_x", "massf"], \
            recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
            comp_label = "P1", flux_scheme = "AUSMPlusUPOriginal", \
            reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactiveMixingPipe(n_cells = NCELLS, geometry = [D, L], \
            init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x", "massf"], \
            recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", update_from = "pT", \
            comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
            reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 7, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
#### Sod's Shock Tube With SinglePhaseMultiSpeciesReactivePipeWithCombustionFraction Components ####
##################################### With Combustion Fraction = 0.0 ###############################
"""
# Left component is SinglePhaseMultiSpeciesReactivePipeWithCombustionFraction
# Right component is SinglePhaseMultiSpeciesReactivePipeWithCombustionFraction
# In both components, combustion fraction = 0.0
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 1MPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

air = Molecule("air")
reaction = Reaction(reactants = [air], products = [air])
s = SystemOfReactions()
s.add_reaction(reaction = reaction)
s.solve_for_coefficients()

NCELLS = 100
L = 0.5
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
print(init_gm_lft.species_names)
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e6
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesReactivePipeWithUniformCombustionFraction(n_cells = NCELLS, \
                geometry = [D, L], init_flow_state = init_fs_lft, comp_label = "P1", \
                combustion_fraction = [[0.0, "air", reaction]], limiter = "eilmer_vanAlbada", \
                recon_scheme = ["copy", [1, 1]], recon_props = ["p", "T", "vel_x", "massf"], \
                update_from = "pT", flux_scheme = "AUSMPlusUPOriginal", \
                reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesReactivePipeWithUniformCombustionFraction(n_cells = NCELLS, \
                geometry = [D, L], init_flow_state = init_fs_rght, comp_label = "P2", \
                combustion_fraction = [[0.0, "air", reaction]], limiter = "eilmer_vanAlbada", \
                recon_scheme = ["copy", [1, 1]], recon_props = ["p", "T", "vel_x", "massf"], \
                update_from = "pT", flux_scheme = "AUSMPlusUPOriginal", \
                reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 8, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
####### Stationary Flow With SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen Components ##########
#################################### With Non-Zero Heat Gen ########################################
"""
# Left component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen
# Right component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen
# Both components have heat gen = 1e6 W/M^3
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
print(init_gm_lft.species_names)
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e5
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen(n_cells = NCELLS, geometry = [D, L], \
                        init_flow_state = init_fs_lft, recon_props = ["p", "T", "vel_x", "massf"], \
                        recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                        update_from = "pT", comp_label = "P1", flux_scheme = "AUSMPlusUPOriginal", \
                        heat_gen = 1e6, reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatGen(n_cells = NCELLS, geometry = [D, L], \
                        init_flow_state = init_fs_rght, recon_scheme = ["copy", [1, 1]], \
                        recon_props = ["p", "T", "vel_x", "massf"], limiter = "eilmer_vanAlbada", \
                        update_from = "pT", comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
                        heat_gen = 1e6, reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 9, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
####### Stationary Flow With SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux Components #########
################################### With Non-Zero Heat Flux ########################################
"""
# Left component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux
# Right component is SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux
# Heat flux is calculated to have the equivalent energy source term as the heat gen case
# Total length = 1.0m
# Diameter = 0.45m
# Both boundary conditions are simple outflow BCs
# Fluid: ideal air
# Left state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 100
L = 0.5
D = 0.45

heat_flux = 0.25 * 1e6 * D

init_gm_lft = GasModel('ideal-air-gas-model.lua')
print(init_gm_lft.species_names)
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e5
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

P1 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux(n_cells = NCELLS, geometry = [D, L], \
                        init_flow_state = init_fs_lft, recon_props = ["p", "T", "vel_x", "massf"], \
                        recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                        update_from = "pT", comp_label = "P1", flux_scheme = "AUSMPlusUPOriginal", \
                        heat_flux = heat_flux, reverse_direction_for_ghost_cells = False)

P1.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipeWithHeatFlux(n_cells = NCELLS, geometry = [D, L], \
                        init_flow_state = init_fs_rght, recon_scheme = ["copy", [1, 1]], \
                        recon_props = ["p", "T", "vel_x", "massf"], limiter = "eilmer_vanAlbada", \
                        update_from = "pT", comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
                        heat_flux = heat_flux, reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = NCELLS, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 10, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
################ Supersonic Inlet Test Case Using InletBlockToMimicMultiInletFlow ##################
###################### and SinglePhaseMultiSpeciesNonReactivePipe Components #######################
################################### Not Run To Steady State ########################################
"""
# Left component is InletBlockToMimicMultiInletFlow -> Single Cell
# Right component is SinglePhaseMultiSpeciesNonReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Left boundary is a wall BC
# Right boundary is a simple outflow BC
# Fluid: ideal air
# Left state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Inflow state: Ma = 2, T = 300K, p = 500kPa
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 200
L = 1.0
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e5
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

inlet_gm = GasModel('ideal-air-gas-model.lua')
inlet_gs = GasState(inlet_gm)
inlet_gs.p = 5e5
inlet_gs.T = 300.0
inlet_gs.update_thermo_from_pT()
Ma_inlet = 2.0
inlet_fs = FlowState(state = inlet_gs, vel_x = Ma_inlet * inlet_gs.a)

P1 = InletBlockToMimicMultiInletFlow(inlet_flow_state = inlet_fs, init_flow_state = init_fs_lft, \
                                    geometry = [D, L / NCELLS], \
                                    inlet_area = 0.25 * np.pi * D ** 2.0, comp_label = "P1", \
                                    recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                                    recon_props = ["p", "T", "vel_x"], update_from = "pT", \
                                    flux_scheme = "AUSMPlusUPOriginal")

P1.add_boundary_conditions(BC = [0, ["WallNoSlip_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS - 1, \
                        geometry = [D, L * (1.0 - 1.0 / NCELLS)], \
                        init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x"], \
                        recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                        update_from = "pT", comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
                        reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS - 1, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = 1, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 11, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
######## Supersonic Inlet Test Case Using SinglePhaseMultiSpeciesNonReactivePipe Component #########
################################### Not Run To Steady State ########################################
"""
# Component is SinglePhaseMultiSpeciesNonReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Left boundary is a supersonic inlet BC
# Right boundary is a simple outflow BC
# Fluid: ideal air
# Initial state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Inflow state: Ma = 2, T = 300K, p = 500kPa
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 200
L = 1.0
D = 0.45

init_gm = GasModel('ideal-air-gas-model.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 300.0
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

inlet_gm = GasModel('ideal-air-gas-model.lua')
inlet_gs = GasState(inlet_gm)
inlet_gs.p = 5e5
inlet_gs.T = 300.0
inlet_gs.update_thermo_from_pT()
Ma_inlet = 2.0
inlet_fs = FlowState(state = inlet_gs, vel_x = Ma_inlet * inlet_gs.a)

Pipe = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS, geometry = [D, L], \
                                init_flow_state = init_fs, recon_props = ["p", "T", "vel_x"], \
                                recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                                update_from = "pT", comp_label = "Pipe", \
                                flux_scheme = "AUSMPlusUPOriginal", \
                                reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SupersonicInFlow_BC", 2, inlet_fs]])
Pipe.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-4, data_save_dt = 1e-4, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 12, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
################ Supersonic Inlet Test Case Using InletBlockToMimicMultiInletFlow ##################
###################### and SinglePhaseMultiSpeciesNonReactivePipe Components #######################
#################################### Run To Steady State ###########################################
"""
# Left component is InletBlockToMimicMultiInletFlow -> Single Cell
# Right component is SinglePhaseMultiSpeciesNonReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Left boundary is a wall BC
# Right boundary is a simple outflow BC
# Fluid: ideal air
# Left state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Inflow state: Ma = 2, T = 300K, p = 500kPa, Area = 0.25 * pi * D ^ 2
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 200
L = 1.0
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e5
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

inlet_gm = GasModel('ideal-air-gas-model.lua')
inlet_gs = GasState(inlet_gm)
inlet_gs.p = 5e5
inlet_gs.T = 300.0
inlet_gs.update_thermo_from_pT()
Ma_inlet = 2.0
inlet_fs = FlowState(state = inlet_gs, vel_x = Ma_inlet * inlet_gs.a)

P1 = InletBlockToMimicMultiInletFlow(inlet_flow_state = inlet_fs, init_flow_state = init_fs_lft, \
                                    geometry = [D, L / NCELLS], \
                                    inlet_area = 0.25 * np.pi * D ** 2.0, comp_label = "P1", \
                                    recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                                    recon_props = ["p", "T", "vel_x"], update_from = "pT", \
                                    flux_scheme = "AUSMPlusUPOriginal")

P1.add_boundary_conditions(BC = [0, ["WallNoSlip_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS - 1, \
                        geometry = [D, L * (1.0 - 1.0 / NCELLS)], \
                        init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x"], \
                        recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                        update_from = "pT", comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
                        reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS - 1, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = 1, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-3, data_save_dt = 1e-3, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 13, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
######## Supersonic Inlet Test Case Using SinglePhaseMultiSpeciesNonReactivePipe Component #########
##################################### Run To Steady State ##########################################
"""
# Component is SinglePhaseMultiSpeciesNonReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Left boundary is a supersonic inlet BC
# Right boundary is a simple outflow BC
# Fluid: ideal air
# Initial state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Inflow state: Ma = 2, T = 300K, p = 500kPa
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 200
L = 1.0
D = 0.45

init_gm = GasModel('ideal-air-gas-model.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 300.0
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

inlet_gm = GasModel('ideal-air-gas-model.lua')
inlet_gs = GasState(inlet_gm)
inlet_gs.p = 5e5
inlet_gs.T = 300.0
inlet_gs.update_thermo_from_pT()
Ma_inlet = 2.0
inlet_fs = FlowState(state = inlet_gs, vel_x = Ma_inlet * inlet_gs.a)

Pipe = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS, geometry = [D, L], \
                                init_flow_state = init_fs, recon_props = ["p", "T", "vel_x"], \
                                recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                                update_from = "pT", comp_label = "Pipe", \
                                flux_scheme = "AUSMPlusUPOriginal", \
                                reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SupersonicInFlow_BC", 2, inlet_fs]])
Pipe.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-3, data_save_dt = 1e-3, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 14, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
############# InletBlockToMimicMultiInletFlow Test Case With Non-Zero Wall Boundary ################
##################### SinglePhaseMultiSpeciesNonReactivePipe To Complete Pipe ######################
################################### Not Run To Steady State ########################################
"""
# Left component is InletBlockToMimicMultiInletFlow -> Single Cell
# Right component is SinglePhaseMultiSpeciesNonReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Left boundary is a wall BC
# Right boundary is a simple outflow BC
# Fluid: ideal air
# Left state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Right state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Inflow state: Ma = 2, T = 300K, p = 500kPa, Area = 3 / 4 * 0.25 * pi * D ^ 2
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 200
L = 1.0
D = 0.45

init_gm_lft = GasModel('ideal-air-gas-model.lua')
init_gs_lft = GasState(init_gm_lft)
init_gs_lft.p = 1e5
init_gs_lft.T = 300.0
init_gs_lft.update_thermo_from_pT()
init_fs_lft = FlowState(state = init_gs_lft, vel_x = 0.0)

init_gm_rght = GasModel('ideal-air-gas-model.lua')
init_gs_rght = GasState(init_gm_rght)
init_gs_rght.p = 1e5
init_gs_rght.T = 300.0
init_gs_rght.update_thermo_from_pT()
init_fs_rght = FlowState(state = init_gs_rght, vel_x = 0.0)

inlet_gm = GasModel('ideal-air-gas-model.lua')
inlet_gs = GasState(inlet_gm)
inlet_gs.p = 5e5
inlet_gs.T = 300.0
inlet_gs.update_thermo_from_pT()
Ma_inlet = 2.0
inlet_area_fraction = 3/4
inlet_fs = FlowState(state = inlet_gs, vel_x = Ma_inlet * inlet_gs.a)

P1 = InletBlockToMimicMultiInletFlow(inlet_flow_state = inlet_fs, init_flow_state = init_fs_lft, \
                                    geometry = [D, L / NCELLS], \
                                    inlet_area = inlet_area_fraction * 0.25 * np.pi * D ** 2.0, \
                                    comp_label = "P1", \
                                    recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                                    recon_props = ["p", "T", "vel_x"], update_from = "pT", \
                                    flux_scheme = "AUSMPlusUPOriginal")

P1.add_boundary_conditions(BC = [0, ["WallNoSlip_BC", 2]])

P2 = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS - 1, \
                        geometry = [D, L * (1.0 - 1.0 / NCELLS)], \
                        init_flow_state = init_fs_rght, recon_props = ["p", "T", "vel_x"], \
                        recon_scheme = ["copy", [1, 1]], limiter = "eilmer_vanAlbada", \
                        update_from = "pT", comp_label = "P2", flux_scheme = "AUSMPlusUPOriginal", \
                        reverse_direction_for_ghost_cells = False)

P2.add_boundary_conditions(BC = [NCELLS - 1, ["SimpleOutFlow_BC", 2]])

Pipe = JointBlock(mesh_object_1 = P1, mesh_object_2 = P2, \
            block_1_interface_id_being_replaced = 1, block_2_interface_id_being_replaced = 0, \
            new_interface = deepcopy(P2.interface_array[0]), adding_ghost_cells_bool = False, \
            new_component_label = "Pipe")

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-3, data_save_dt = 1e-3, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 15, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
#### SinglePhaseMultiSpeciesMultiInletNonReactivePipe Test Case With Wall and SuperSonic Inlets ####
################################### Not Run To Steady State ########################################
"""
# Component is SinglePhaseMultiSpeciesMultiInletNonReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Left boundary:
#   - Supersonic inlet on 3/4 of total area
#   - Wall BC on 1/4 of total area
# Right boundary is a simple outflow BC
# Fluid: ideal air
# Initial state: p = 100kPa, T = 300.0K, vel_x = 0.0ms^-1
# Inflow state: Ma = 2, T = 300K, p = 500kPa
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada

NCELLS = 200
L = 1.0
D = 0.45
total_area = 0.25 * np.pi * D ** 2.0
inlet_area_fraction = 3/4

init_gm = GasModel('ideal-air-gas-model.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 300.0
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

inlet_gm = GasModel('ideal-air-gas-model.lua')
inlet_gs = GasState(inlet_gm)
inlet_gs.p = 5e5
inlet_gs.T = 300.0
inlet_gs.update_thermo_from_pT()
Ma_inlet = 2.0
inlet_fs = FlowState(state = inlet_gs, vel_x = Ma_inlet * inlet_gs.a)

Pipe = SinglePhaseMultiSpeciesMultiInletNonReactivePipe(n_cells = NCELLS, n_inlets = 2, \
                    bulk_geometry = [D, L], \
                    inlet_areas = [inlet_area_fraction * total_area, (1.0 - inlet_area_fraction) * total_area], \
                    init_flow_state = init_fs, comp_label = "Pipe", limiter = "eilmer_vanAlbada", \
                    recon_scheme = ["copy", [1, 1]], recon_props = ["p", "T", "vel_x"], \
                    update_from = "pT", flux_scheme = "AUSMPlusUPOriginal", \
                    reverse_direction_for_mirrored_flow = False)

Pipe.add_boundary_conditions(BC = [0, ["SupersonicInFlow_BC", 2, inlet_fs]])
Pipe.add_boundary_conditions(BC = [1, ["WallNoSlip_BC", 2]])

Pipe.add_boundary_conditions(BC = [NCELLS + 1, ["SimpleOutFlow_BC", 2]])

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-3, data_save_dt = 1e-3, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 16, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
################ Multiple Reaction Test Case To Mimic Multiple Reversible Reactions ################
########## Using SinglePhaseMultiSpeciesReactivePipeWithUniformCombustionFraction Component ########
"""
# Component is SinglePhaseMultiSpeciesReactivePipeWithUniformCombustionFraction
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 2 species nitrogen
# Initial state: p = 100kPa, T = 4000K, vel_x = 0.0ms^-1
# Reconstruction: Copy
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada
# Trying to mimic 2 reversible reactions, so 4 reactions are defined.
#   - Rates for each reaction must be tuned to approach finite rate model.

N2_r1 = Molecule(chemical_form = "N2")
N_r1 = Molecule(chemical_form = "N")

r1 = Reaction(reactants = [N2_r1, N2_r1], products = [N_r1, N_r1, N2_r1])

N2_r2 = Molecule(chemical_form = "N2")
N_r2 = Molecule(chemical_form = "N")

r2 = Reaction(reactants = [N_r2, N_r2, N2_r2], products = [N2_r2, N2_r2])

N2_r3 = Molecule(chemical_form = "N2")
N_r3 = Molecule(chemical_form = "N")

r3 = Reaction(reactants = [N2_r3, N_r3], products = [N_r3, N_r3, N_r3])

N2_r4 = Molecule(chemical_form = "N2")
N_r4 = Molecule(chemical_form = "N")

r4 = Reaction(reactants = [N_r4, N_r4, N_r4], products = [N2_r4, N_r4])

s = SystemOfReactions()
s.add_reaction(reaction = r1)
s.add_reaction(reaction = r2)
s.add_reaction(reaction = r3)
s.add_reaction(reaction = r4)

s.solve_for_coefficients()

L = 1.0
D = 0.45

init_gm = GasModel('nitrogen-2sp.lua')
print(init_gm.species_names)
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 4000.0
init_gs.massf = [0.8, 0.2]
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)
reaction_speed_up_coefficient = 3.0
x_N2 = 0.05000029 * reaction_speed_up_coefficient
x_N = 0.231459 * reaction_speed_up_coefficient

Pipe = SinglePhaseMultiSpeciesReactivePipeWithUniformCombustionFraction(n_cells = 1, \
        geometry = [D, L], init_flow_state = init_fs, comp_label = "Pipe", \
        combustion_fraction = [[x_N2, "N2", r1], [x_N2, "N2", r2], [x_N2, "N2", r3], [x_N, "N", r4]], \
        limiter = "eilmer_vanAlbada_L2R2", recon_scheme = ["eilmer_L3R3", [3, 3]], \
        recon_props = ["p", "T", "vel_x", "massf"], update_from = "pT", \
        flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 1]])
Pipe.add_boundary_conditions(BC = [1, ["SimpleOutFlow_BC", 1]])
Pipe.cell_idx_to_track = [0]

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-3, data_save_dt = 1e-3, \
            transient_cell_flow_property_variables_to_write = ["p", "T", "rho", "u", "massf"], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 17, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "u", "massf"], \
            rapid_data_save_steps = 1)
"""
####################################################################################################
#### Nitrogen Reactor Test Case Using SinglePhaseMultiSpeciesFiniteRateReactivePipe Component ######
"""
# Component is SinglePhaseMultiSpeciesFiniteRateReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 2 species nitrogen
# Initial state: p = 100kPa, T = 4000K, vel_x = 0.0ms^-1
# Reconstruction: Starting at L3R3, testing that the hierarchy system works
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2

D = 0.45
L = 1.0

init_gm = GasModel('nitrogen-2sp.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 4000.0
init_gs.massf = [0.8, 0.2]
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

reactor = ThermochemicalReactor(gmodel = init_gm, filename1 = "chem.lua")

Pipe = SinglePhaseMultiSpeciesFiniteRateReactivePipe(n_cells = 1, geometry = [D, L], \
            limiter = "eilmer_vanAlbada_L2R2", comp_label = "Pipe", init_flow_state = init_fs, \
            recon_scheme = ["eilmer_L3R3", [3, 3]], recon_props = ["p", "T", "vel_x", "massf"], \
            update_from = "pT", flux_scheme = "AUSMPlusUPOriginal", \
            shared_reactor_object = reactor, reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 1]])
Pipe.add_boundary_conditions(BC = [1, ["SimpleOutFlow_BC", 1]])

Pipe.cell_idx_to_track = [0]

S = Solver(mesh_object = Pipe, \
        cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
        t_final = 5e-3, data_save_dt = 1e-3, \
        transient_cell_flow_property_variables_to_write = ["p", "T", "rho", "u", "massf"], \
        transient_interface_flow_property_variables_to_write = [], \
        sim_number = 18, \
        spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "u", "massf"], \
        rapid_data_save_steps = 1)
"""
####################################################################################################
### Nitrogen Reactor Test Case Using SinglePhaseMultiSpeciesReactivePipeWithUniformFastChemistry ###
######################################## Component #################################################
"""
# Component is SinglePhaseMultiSpeciesReactivePipeWithUniformFastChemistry
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 2 species nitrogen
# Initial state: p = 100kPa, T = 4000K, vel_x = 0.0ms^-1
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 1e-3

D = 0.45
L = 1.0

NCELLS = 1

init_gm = GasModel('nitrogen-2sp.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 4000.0
init_gs.massf = [0.8, 0.2]
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

equil_gm = GasModel('cea-n2-gas-model.lua')
equil_gs = GasState(equil_gm)
fast_chemistry_fraction = 0.3

Pipe = SinglePhaseMultiSpeciesReactivePipeWithUniformFastChemistry(n_cells = NCELLS, \
            comp_label = "Pipe", geometry = [D, L], init_flow_state = init_fs, \
            fast_chemistry_fraction = fast_chemistry_fraction, fast_chemistry_gs = equil_gs, \
            limiter = "eilmer_vanAlbada_L2R2", recon_scheme = ["eilmer_L1R1", [1, 1]], \
            recon_props = ["p", "T", "vel_x", "massf"], update_from = "pT", \
            flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 1]])
Pipe.add_boundary_conditions(BC = [1, ["SimpleOutFlow_BC", 1]])

Pipe.cell_idx_to_track = [0]

S = Solver(mesh_object = Pipe, \
        cfl_flag = [False, 1e-6], \
        t_final = 1e-3, data_save_dt = 2e-4, \
        transient_cell_flow_property_variables_to_write = ["p", "T", "rho", "u", "massf"], \
        transient_interface_flow_property_variables_to_write = [], \
        sim_number = 19, \
        spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "u", "massf"], \
        rapid_data_save_steps = 1)
"""
####################################################################################################
# Nitrogen Reactor Test Case Using 
# SinglePhaseMultiSpeciesReactivePipeWithPartialEquilibriumChemistry ###
######################################## Component #################################################
"""
# Component is SinglePhaseMultiSpeciesReactivePipeWithPartialEquilibriumChemistry
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 2 species nitrogen
# Initial state: p = 100kPa, T = 4000K, vel_x = 0.0ms^-1
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 1e-3

D = 0.45
L = 1.0

NCELLS = 1

init_gm = GasModel('nitrogen-2sp.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 4000.0
init_gs.massf = [0.8, 0.2]
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

equil_gm = GasModel('cea-n2-gas-model.lua')
equil_gs = GasState(equil_gm)
fast_chemistry_fraction = 0.3
partial_equilibrium_fraction = 0.5

Pipe = SinglePhaseMultiSpeciesReactivePipeWithPartialEquilibriumChemistry(n_cells = NCELLS, \
            comp_label = "Pipe", geometry = [D, L], init_flow_state = init_fs, \
            partial_equilibrium_fraction= partial_equilibrium_fraction, \
            fast_chemistry_fraction = fast_chemistry_fraction, fast_chemistry_gs = equil_gs, \
            limiter = "eilmer_vanAlbada_L2R2", recon_scheme = ["eilmer_L1R1", [1, 1]], \
            recon_props = ["p", "T", "vel_x", "massf"], update_from = "pT", \
            flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 1]])
Pipe.add_boundary_conditions(BC = [1, ["SimpleOutFlow_BC", 1]])

Pipe.cell_idx_to_track = [0]

S = Solver(mesh_object = Pipe, \
        cfl_flag = [False, 1e-6], \
        t_final = 1e-3, data_save_dt = 2e-4, \
        transient_cell_flow_property_variables_to_write = ["p", "T", "rho", "u", "massf"], \
        transient_interface_flow_property_variables_to_write = [], \
        sim_number = 20, \
        spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "u", "massf"], \
        rapid_data_save_steps = 1)
"""
####################################################################################################
# H2 + O2 Reactor Test Case Using 
# SinglePhaseMultiSpeciesFiniteRateReactivePipe ###
######################################## Component #################################################
"""
# Component is SinglePhaseMultiSpeciesFiniteRateReactivePipe
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 6 species O, O2, H, H2, OH, H2O
# Initial state: p = 100kPa, T = 4000K, vel_x = 0.0ms^-1
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 1e-3

D = 0.45
L = 1.0

NCELLS = 1

init_gm = GasModel('H2_O2_6sp.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 4000.0
print(init_gm.species_names)
#init_gs.massf = [0.8, 0.2]
#init_gs.update_thermo_from_pT()
#init_fs = FlowState(state = init_gs, vel_x = 0.0)
"""
####################################################################################################
# N2 + N Reactor Test Case Using 
# SinglePhaseMultiSpeciesReactivePipeWithFastChemistryForGoalMassf ###
######################################## Component #################################################
"""
# Component is SinglePhaseMultiSpeciesReactivePipeWithFastChemistryForGoalMassf
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 6 species O, O2, H, H2, OH, H2O
# Initial state: p = 100kPa, T = 4000K, vel_x = 0.0ms^-1
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 1e-3

D = 0.45
L = 1.0

NCELLS = 1

init_gm = GasModel('nitrogen-2sp.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 4000.0
init_gs.massf = [0.8, 0.2]
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

N = Molecule(chemical_form = "N")
N2 = Molecule(chemical_form = "N2")

r = Reaction(reactants = [N], products = [N2])
s = SystemOfReactions()
s.add_reaction(reaction = r)
s.solve_for_coefficients()
print("reactants_stoichiometric_coefficients", r.reactants_stoichiometric_coefficients)

combustion_parameters = ["N2", 0.86, 0.03, r]

Pipe = SinglePhaseMultiSpeciesReactivePipeWithFastChemistryForGoalMassf(n_cells = NCELLS, \
            geometry = [D, L], init_flow_state = init_fs, comp_label = "Pipe", \
            combustion_parameters = combustion_parameters, limiter = "eilmer_vanAlbada_L2R2", \
            recon_scheme = ["eilmer_L1R1", [1, 1]], recon_props = ["p", "T", "vel_x", "massf"], \
            update_from = "pT", flux_scheme = "AUSMPlusUPOriginal", \
            reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SimpleOutFlow_BC", 1]])
Pipe.add_boundary_conditions(BC = [1, ["SimpleOutFlow_BC", 1]])

Pipe.cell_idx_to_track = [0]

S = Solver(mesh_object = Pipe, \
        cfl_flag = [False, 1e-6], \
        t_final = 1e-3, data_save_dt = 2e-4, \
        transient_cell_flow_property_variables_to_write = ["p", "T", "rho", "u", "massf"], \
        transient_interface_flow_property_variables_to_write = [], \
        sim_number = 22, \
        spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "u", "massf"], \
        rapid_data_save_steps = 1)
"""
####################################################################################################
### Pushing species out of pipe with SuperSonic Inlet Using  ###
#################### SinglePhaseMultiSpeciesNonReactiveMixingPipe Component ########################
"""
# Component is SinglePhaseMultiSpeciesNonReactiveMixingPipe
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 2 species N2, Ar
# Initial state: p = 100kPa, T = 300K, vel_x = 0.0ms^-1, All Ar
# Inlet: p = 200kPa, T = 300K, Ma = 2, All N2
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 1e-3

D = 0.45
L = 1.0

NCELLS = 100

init_gm = GasModel('N2-Ar-thermally-perfect-gas-model.lua')
init_gs = GasState(init_gm)
init_gs.p = 1e5
init_gs.T = 300.0
init_gs.massf = [0.0, 1.0]
init_gs.update_thermo_from_pT()
init_fs = FlowState(state = init_gs, vel_x = 0.0)

inlet_Ma = 2.0
inlet_gm = GasModel('N2-Ar-thermally-perfect-gas-model.lua')
inlet_gs = GasState(inlet_gm)
inlet_gs.p = 1e5
inlet_gs.T = 300.0
inlet_gs.massf = [1.0, 0.0]
inlet_gs.update_thermo_from_pT()
inlet_vel_x = inlet_Ma * inlet_gs.a
inlet_fs = FlowState(state = inlet_gs, vel_x = inlet_vel_x)

Pipe = SinglePhaseMultiSpeciesNonReactiveMixingPipe(n_cells = NCELLS, geometry = [D, L], \
                init_flow_state = init_fs, recon_props = ["p", "T", "vel_x", "massf"], \
                recon_scheme = ["eilmer_L1R1", [1, 1]], limiter = "eilmer_vanAlbada_L2R2", \
                update_from = "pT", comp_label = "Pipe", flux_scheme = "AUSMPlusUPOriginal", \
                reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["SupersonicInFlow_BC", 2, inlet_fs]])
Pipe.add_boundary_conditions(BC = [NCELLS, ["SimpleOutFlow_BC", 2]])

S = Solver(mesh_object = Pipe, \
        cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
        t_final = 1e-3, data_save_dt = 2e-4, \
        transient_cell_flow_property_variables_to_write = ["p", "T", "rho", "u", "massf"], \
        transient_interface_flow_property_variables_to_write = [], \
        sim_number = 23, \
        spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "u", "massf"], \
        rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
### Testing CEA model
"""
gm = GasModel("cea-H2-O2-H2O-gas-model.lua")
gs = GasState(gm)
gs.p = 2.1312e7
gs.T = 859.26

gs.update_thermo_from_pT()
print(gs.ceaSavedData)
"""
####################################################################################################
### SSME model Configuration 1 using goal mass fraction checmistry approach
"""
# Component is SinglePhaseMultiSpeciesNonReactiveMixingPipe
# Total length = 1.0m
# Diameter = 0.45m
# Both boundaries are simple outflow BCs
# Fluid: 2 species N2, Ar
# Initial state: p = 100kPa, T = 300K, vel_x = 0.0ms^-1, All Ar
# Inlet: p = 200kPa, T = 300K, Ma = 2, All N2
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 1e-3

D1 = 0.45
L1 = 0.25

NCELLS1 = 100

comb_chamber_init_gm = GasModel("H2-O2-H2O-thermally-perfect-gas-model.lua")
print(comb_chamber_init_gm.species_names)
comb_chamber_init_gs = GasState(comb_chamber_init_gm)
#comb_chamber_init_massf = (np.array([0.0567333, 0.87496, 0.068436])/sum([0.87496, 0.0567333, 0.068436])).tolist()
#comb_chamber_init_gs.massf = comb_chamber_init_massf
#print(comb_chamber_init_massf)
#print(sum(comb_chamber_init_massf))
comb_chamber_init_molef = [2/3, 1/3, 0.0]
comb_chamber_init_gs.molef = comb_chamber_init_molef
comb_chamber_init_gs.p = 1.9795e7
#comb_chamber_init_gs.T = 3588.706
comb_chamber_init_gs.T = 300.0
comb_chamber_init_gs.update_thermo_from_pT()
comb_chamber_init_fs = FlowState(state = comb_chamber_init_gs, vel_x = 0.0)

H2 = Molecule(chemical_form = "H2")
O2 = Molecule(chemical_form = "O2")
H2O = Molecule(chemical_form = "H2O")

r = Reaction(reactants = [H2, O2], products = [H2O])

s = SystemOfReactions()
s.add_reaction(reaction = r)
s.solve_for_coefficients()


comb_chamber_inlet_gm = GasModel("H2-O2-H2O-thermally-perfect-gas-model.lua")
comb_chamber_inlet_gs = GasState(comb_chamber_inlet_gm)
inlet_BC_molef = comb_chamber_init_molef
comb_chamber_inlet_gs.molef = inlet_BC_molef
comb_chamber_inlet_gs.T = 300.0
comb_chamber_inlet_gs.p = 2.1339e7
comb_chamber_inlet_gs.update_thermo_from_pT()
print(comb_chamber_inlet_gs.massf)


Pipe = SinglePhaseMultiSpeciesReactivePipeWithTanhFastChemistryForGoalMassfProfile(n_cells = NCELLS1, \
            comp_label = "Pipe", geometry = [D1, L1], init_flow_state = comb_chamber_init_fs, \
            goal_massf_profile = ["H2O", comb_chamber_inlet_gs.massf[2], 1.0, 0.01, r, 6.0, L1/2.0], \
            limiter = "eilmer_vanAlbada_L2R2", recon_scheme = ["eilmer_L1R1", [1, 1]], \
            recon_props = ["p", "T", "vel_x", "massf"], update_from = "pT", \
            flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

#Pipe.add_boundary_conditions(BC = [0, ["FromStagnationWithMassFlowRateInFlow_BC", 2, [comb_chamber_inlet_gs, 476.090552]]])
Pipe.add_boundary_conditions(BC = [0, ["WallNoSlip_BC", 2]])
Pipe.add_boundary_conditions(BC = [NCELLS1, ["WallNoSlip_BC", 2]])

S = Solver(mesh_object = Pipe, \
        cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
        t_final = 8e-6, data_save_dt = 2e-6, \
        transient_cell_flow_property_variables_to_write = ["p", "T", "rho", "u", "massf"], \
        transient_interface_flow_property_variables_to_write = [], \
        sim_number = 25, \
        spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "u", "massf"], \
        rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
##### Nonreactive Pipe and Nozzle Testing
"""
# Pipe component is SinglePhaseMultiSpeciesNonReactivePipe
# Converging nozzle component is SinglePhaseMultiSpeciesNonReactiveQuadraticNozzle
# Diverging nozzle components is SinglePhaseMultiSpeciesNonReactiveEllipticNozzle
# Pipe length = 0.25m
# Pipe Diameter = 0.45m
# Converging nozzle length = 0.12m
# Throat Diameter = 0.276m
# Diverging nozzle length = 1.0m
# Exit Diameter = 0.618m
# Inlet BC = stagnation with mass flow rate
# Outlet BC = OutFlow
# Fluid: 3 species H2, O2, H2O
# Initial state: p = 100kPa, T = 300K, vel_x = 0.0ms^-1
# Inlet: p = 200kPa, T = 300K, mdot = 400kgs^-1
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 1e-1

NCELLS1 = 50
D_1 = 0.45
L_1 = 0.25

NCELLS2 = 20
D_t = 0.276
L_2 = 0.12

NCELLS3 = 60
D_e = 0.618
L_3 = 0.43


gm_init = GasModel('ideal-air-gas-model.lua')
gs_init = GasState(gm_init)
gs_init.p = 100e3
gs_init.T = 300.0
gs_init.update_thermo_from_pT()
fs_init = FlowState(state = gs_init)

gm_inlet = GasModel('ideal-air-gas-model.lua')
gs_inlet = GasState(gm_inlet)
gs_inlet.p = 500e3
gs_inlet.T = 300.0
gs_inlet.update_thermo_from_pT()

Pipe = SinglePhaseMultiSpeciesNonReactivePipe(n_cells = NCELLS1, geometry = [D_1, L_1], \
                init_flow_state = fs_init, recon_props = ["p", "T", "vel_x"], \
                recon_scheme = ["eilmer_L1R1", [1, 1]], limiter = "eilmer_vanAlbada_L2R2", \
                update_from = "pT", comp_label = "Pipe", flux_scheme = "AUSMPlusUPOriginal", \
                reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["FromStagnationWithMassFlowRateInFlow_BC", 2, [gs_inlet, 400.0]]])

Conv_nozzle = SinglePhaseMultiSpeciesNonReactiveQuadraticNozzle(n_cells = NCELLS2, \
                geometry = [D_1, D_t, L_2], init_flow_state = fs_init, \
                recon_props = ["p", "T", "vel_x"], recon_scheme = ["eilmer_L1R1", [1, 1]], \
                limiter = "eilmer_vanAlbada_L2R2", update_from = "pT", comp_label = "Conv_nozzle", \
                flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

Div_nozzle = SinglePhaseMultiSpeciesNonReactiveEllipticNozzle(n_cells = NCELLS3, \
                geometry = [D_t, D_e, L_3], init_flow_state = fs_init, \
                recon_props = ["p", "T", "vel_x"], recon_scheme = ["eilmer_L1R1", [1, 1]], \
                limiter = "eilmer_vanAlbada_L2R2", update_from = "pT", comp_label = "Div_nozzle", \
                flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

Div_nozzle.add_boundary_conditions(BC = [NCELLS3, ["SimpleOutFlow_BC", 2]])

Pipe_and_conv_nozzle = JointBlock(mesh_object_1 = Pipe, mesh_object_2 = Conv_nozzle, \
                                    block_1_interface_id_being_replaced = NCELLS1, \
                                    block_2_interface_id_being_replaced = 0, \
                                    new_interface = deepcopy(Conv_nozzle.interface_array[0]), \
                                    adding_ghost_cells_bool = False, new_component_label = "Pipe")

Full_component = JointBlock(mesh_object_1 = Pipe_and_conv_nozzle, mesh_object_2 = Div_nozzle, \
                            block_1_interface_id_being_replaced = NCELLS1 + NCELLS2, \
                            block_2_interface_id_being_replaced = 0, \
                            new_interface = deepcopy(Div_nozzle.interface_array[0]), \
                            adding_ghost_cells_bool = False, new_component_label = "Rocket")

S = Solver(mesh_object = Full_component, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-5, data_save_dt = 2e-2, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 1, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "Ma"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
##### Nonreactive Pipe and Nozzle Testing With 5 species air
"""
# Pipe component is SinglePhaseMultiSpeciesNonReactiveMixingPipe
# Converging nozzle component is SinglePhaseMultiSpeciesNonReactiveMixingQuadraticNozzle
# Diverging nozzle components is SinglePhaseMultiSpeciesNonReactiveMixingEllipticNozzle
# Pipe length = 0.25m
# Pipe Diameter = 0.45m
# Converging nozzle length = 0.12m
# Throat Diameter = 0.276m
# Diverging nozzle length = 1.0m
# Exit Diameter = 0.618m
# Inlet BC = stagnation with mass flow rate
# Outlet BC = OutFlow
# Fluid: 3 species H2, O2, H2O
# Initial state: p = 100kPa, T = 600K, vel_x = 0.0ms^-1, [H2, O2, H2O] = [0.036, 0.006, 0.958] massf
# Inlet: p = 200kPa, T = 600K, mdot = , [H2, O2, H2O] = [0.036, 0.006, 0.958] massf
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 5e-2

NCELLS1 = 50
D_1 = 0.45
L_1 = 0.25

NCELLS2 = 20
D_t = 0.276
L_2 = 0.12

NCELLS3 = 60
D_e = 0.618
L_3 = 0.43

gm_init = GasModel('five-species-air-N2-O2-Ar-CO2-Ne.lua')
gs_init = GasState(gm_init)
gs_init.p = 100e3
gs_init.T = 300.0
massf_list = [0.75511, 0.2314, 0.0129, 0.00063, 0.000013]
massf_list = (np.array(massf_list) / sum(massf_list)).tolist()
gs_init.massf = massf_list
gs_init.update_thermo_from_pT()
fs_init = FlowState(state = gs_init)

gm_inlet = GasModel('five-species-air-N2-O2-Ar-CO2-Ne.lua')
gs_inlet = GasState(gm_inlet)
gs_inlet.p = 500e3
gs_inlet.T = 300.0
gs_inlet.massf = massf_list
gs_inlet.update_thermo_from_pT()

Pipe = SinglePhaseMultiSpeciesNonReactiveMixingPipe(n_cells = NCELLS1, geometry = [D_1, L_1], \
                    init_flow_state = fs_init, recon_props = ["p", "T", "vel_x", "massf"], \
                    recon_scheme = ["eilmer_L1R1", [1, 1]], limiter = "eilmer_vanAlbada_L2R2", \
                    update_from = "pT", comp_label = "Pipe", flux_scheme = "AUSMPlusUPOriginal", \
                    reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["FromStagnationWithMassFlowRateInFlow_BC", 2, [gs_inlet, 400.0]]])

Conv_nozzle = SinglePhaseMultiSpeciesNonReactiveMixingQuadraticNozzle(n_cells = NCELLS2, \
                    geometry = [D_1, D_t, L_2], init_flow_state = fs_init, \
                    recon_props = ["p", "T", "vel_x", "massf"], \
                    recon_scheme = ["eilmer_L1R1", [1, 1]], limiter = "eilmer_vanAlbada_L2R2", \
                    update_from = "pT", comp_label = "Conv_nozzle", \
                    flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

Div_nozzle = SinglePhaseMultiSpeciesNonReactiveMixingEllipticNozzle(n_cells = NCELLS3, \
                    geometry = [D_t, D_e, L_3], init_flow_state = fs_init, \
                    recon_props = ["p", "T", "vel_x", "massf"], \
                    recon_scheme = ["eilmer_L1R1", [1, 1]], limiter = "eilmer_vanAlbada_L2R2", \
                    update_from = "pT", comp_label = "Div_nozzle", \
                    flux_scheme = "AUSMPlusUPOriginal", reverse_direction_for_ghost_cells = False)

Div_nozzle.add_boundary_conditions(BC = [NCELLS3, ["SimpleOutFlow_BC", 2]])

Pipe_and_conv_nozzle = JointBlock(mesh_object_1 = Pipe, mesh_object_2 = Conv_nozzle, \
                                    block_1_interface_id_being_replaced = NCELLS1, \
                                    block_2_interface_id_being_replaced = 0, \
                                    new_interface = deepcopy(Conv_nozzle.interface_array[0]), \
                                    adding_ghost_cells_bool = False, new_component_label = "Pipe")

Full_component = JointBlock(mesh_object_1 = Pipe_and_conv_nozzle, mesh_object_2 = Div_nozzle, \
                            block_1_interface_id_being_replaced = NCELLS1 + NCELLS2, \
                            block_2_interface_id_being_replaced = 0, \
                            new_interface = deepcopy(Div_nozzle.interface_array[0]), \
                            adding_ghost_cells_bool = False, new_component_label = "Rocket")

S = Solver(mesh_object = Full_component, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 5e-2, data_save_dt = 1e-2, \
            transient_cell_flow_property_variables_to_write = [], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 5, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "Ma", "massf"], \
            rapid_data_save_steps = 1_000_000)
"""
####################################################################################################
##### FV reactor with 19 reactions, 11sp
"""
# Pipe component is SinglePhaseMultiSpeciesFiniteRateReactivePipe
# Pipe length = 0.1m
# Pipe Diameter = 0.1m
# Inlet BC = Wall 
# Outlet BC = Wall
# Fluid: 11sp, H2, O2, O, H, OH, H2O, HO2, H2O2, He, Ar, N2
# Initial state: p = 0.3atm, T = 880K, 99% N2, 0.5% O2, 0.5% H2
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 0.035

D = 0.1
L = 0.1
NCELLS = 1

gm_init = GasModel('H2-O2-11sp-thermally-perfect-gas-model.lua')
gs_init = GasState(gm_init)
gs_init.p = 30397.5
gs_init.T = 880.0
# species names:    ['H',   'H2',   'O',    'O2',   'H2O2', 'H2O',  'HO2',  'OH',   'Ar',   'He',   'N2']
molef =             [0.0,   0.005,  0.0,    0.005,  0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.99]
gs_init.molef = molef

gs_init.update_thermo_from_pT()
fs_init = FlowState(state = gs_init, vel_x = 0.0)

reactor = ThermochemicalReactor(gmodel = gm_init, \
            filename1 = "hydrogen-oxygen-19r-chain-dissociation-recombination-and-formation-chemistry.lua")

Pipe = SinglePhaseMultiSpeciesFiniteRateReactivePipe(n_cells = NCELLS, geometry = [D, L], \
            limiter = "eilmer_van_albada_L2R2", comp_label = "Pipe", init_flow_state = fs_init, \
            recon_scheme = ["eilmer_L1R1", [1, 1]], recon_props = ["p", "T", "vel_x", "massf"], \
            update_from = "pT", flux_scheme = "AUSMPlusUPOriginal", \
            shared_reactor_object = reactor, \
            reverse_direction_for_ghost_cells = False)

Pipe.add_boundary_conditions(BC = [0, ["WallWithSlip_BC", 1]])
Pipe.add_boundary_conditions(BC = [1, ["WallWithSlip_BC", 1]])

Pipe.cell_idx_to_track = [0]

S = Solver(mesh_object = Pipe, \
            cfl_flag = [True, "cfl_ramp_from_tCurrent", [[0.01, 1e-5], [0.05, 1e-4]]], \
            t_final = 0.1, data_save_dt = 0.1, \
            transient_cell_flow_property_variables_to_write = ["massf", "molef", "conc", "p", "T", "rho", "Ma", "vel_x"], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 8, \
            spatial_cell_flow_property_variables_to_write = ["p", "T", "rho", "vel_x", "Ma", "massf"], \
            rapid_data_save_steps = 50)
"""
####################################################################################################
##### FV reactor with 18 reactions, 9sp

# Pipe component is SinglePhaseMultiSpeciesFiniteRateReactivePipe
# Pipe length = 0.1m
# Pipe Diameter = 0.1m
# Inlet BC = Wall 
# Outlet BC = Wall
# Fluid: 11sp, H2, O2, O, H, OH, H2O, HO2, H2O2, He, Ar, N2
# Initial state: p = 0.3atm, T = 880K, 99% N2, 0.5% O2, 0.5% H2
# Reconstruction: eilmer_L1R1
# Flux Scheme: AUSMPlusUPOriginal
# Limiter: eilmer_van_albada_L2R2
# t_final = 0.1

D = 0.1
L = 0.1
NCELLS = 1

gm_init = GasModel('h2-o2-n2-9sp.lua')
gs_init = GasState(gm_init)
#print(gm_init.species_names)
# species names:    ['O',   'O2',   'N2',   'H',    'H2',   'H2O',  'HO2',  'OH',   'H2O2']
molef =             [0.0,   0.005,  0.99,   0.0,    0.005,  0.0,    0.0,    0.0,    0.0]
gs_init.molef = molef
gs_init.p = 100e3
gs_init.T = 2000.0

gs_init.update_thermo_from_pT()
fs_init = FlowState(state = gs_init, vel_x = 0.0)

gm_reactor = GasModel('h2-o2-n2-9sp.lua')

reactor = ThermochemicalReactor(gmodel = gm_reactor, \
            filename1 = "h2-o2-n2-9sp-18r.lua")

Pipe = SinglePhaseMultiSpeciesFiniteRateReactivePipe(n_cells = NCELLS, geometry = [D, L], \
            limiter = "eilmer_van_albada_L2R2", comp_label = "Pipe", init_flow_state = fs_init, \
            recon_scheme = ["eilmer_L1R1", [1, 1]], recon_props = ["p", "T", "vel_x", "massf"], \
            update_from = "pT", flux_scheme = "AUSMPlusUPOriginal", \
            shared_reactor_object = reactor, reverse_direction_for_ghost_cells = False)
print("Starting species mass cq:", Pipe.cell_array[0].cqs["spcs_mass"]) 
print("Starting mass fractions: ", Pipe.cell_array[0].flow_state.fluid_state.massf)       

Pipe.add_boundary_conditions(BC = [0, ["WallWithSlip_BC", 1]])
Pipe.add_boundary_conditions(BC = [1, ["WallWithSlip_BC", 1]])

Pipe.cell_idx_to_track = [0]

S = Solver(mesh_object = Pipe, \
            cfl_flag = [False, 1e-7], \
            t_final = 1e-7, data_save_dt = 0.1, \
            transient_cell_flow_property_variables_to_write = ["massf", "molef", "conc", "T", "p", "rho", "u", "Ma", "vel_x"], \
            transient_interface_flow_property_variables_to_write = [], \
            sim_number = 1, \
            spatial_cell_flow_property_variables_to_write = ["T", "p", "rho", "u", "vel_x", "Ma", "massf"], \
            rapid_data_save_steps = 1)