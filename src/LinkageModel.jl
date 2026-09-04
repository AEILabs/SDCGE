# Rewritten by John-Mary Matovu
#             Research Fellow
#             Research for Transformation and Development
#             Kampala, Uganda
#             www.rtdug.org

module LinkageModel

using JuMP
using LinearAlgebra
using DelimitedFiles
using XLSX
using Complementarity
using PATHSolver
using DataFrames
# PATH licence: use PATH_LICENSE_STRING from the environment if set, otherwise
# fall back to the public courtesy licence (valid until 31 Dec 2035).
PATHSolver.c_api_License_SetString(get(ENV, "PATH_LICENSE_STRING",
    "1259252040&Courtesy&&&USR&GEN2035&5_1_2026&1000&PATH&GEN&31_12_2035&0_0_0&6000&0_0"))

include("Types.jl")
include("SAM.jl")
include("Calibration.jl")
include("ParameterTables.jl")
include("Functions.jl")
include("Initialization.jl")
include("Variables.jl")
include("Production.jl")
include("Income.jl")
include("Demand.jl")
include("Trade.jl")
include("Equilibrium.jl")
include("Closure.jl")
include("Factors.jl")
include("Other.jl")
# Dynamics.jl and Accounting.jl removed: their equations referenced legacy
# variables (RGDPMP, LambdaL/K, AIDADS, EVFA, ...) that were dropped during
# the static-model cleanup. Recursive dynamics is implemented in
# RecursiveDynamic.jl / PolicyScenarios.jl instead.
include("Results.jl")
# Plotting.jl is the only file that needs Plots.jl (and, transitively, GR and ~200 further
# packages).  Callers that render their own charts — a web backend, a headless batch run —
# can set LCGE_NO_PLOTS=1 before loading this file to skip it, in which case plot_results,
# plot_dynamic_results, plot_all_scenarios and export_results_and_plots! are simply not
# defined (the names below are still exported, which is legal and harmless).
const LCGE_WITH_PLOTS = !(lowercase(get(ENV, "LCGE_NO_PLOTS", "")) in ("1", "true", "yes"))
if LCGE_WITH_PLOTS
    include("Plotting.jl")
end
include("Diagnostics.jl")
include("ModelBuilder.jl")
include("RecursiveDynamic.jl")
include("PolicyScenarios.jl")

export LinkageData, init_data, default_sets!, read_sets_csv!, setup_sam_accounts!, build_default_large_sam!,
       read_sam_csv!, read_sam_excel!, validate_sam!, balance_sam_ras!, sam_balance_table, sam_balance_summary, assert_balanced_sam!, export_sam_balance_report!, export_balanced_sam!, calibrate_from_sam!,
       precompute_parameters, sanitize_parameters!, parameter_diagnostics, parameters, prepare_data!, add_variables!, initialize_from_sam!, enforce_nlp_safe_bounds_and_starts!, initialization_diagnostics, check_initialization!,
       results_dataframe, export_results!, solve_and_save!, plot_results, export_results_and_plots!,
       ces_demand, ces_price, ces_alpha, cet_supply, cet_price, cet_gamma_from_primal, cet_primal_from_gamma,
       armington_market_shares, armington_beta, armington_alpha_update, armington_long_run_elasticity, build_linkage_model!,
       build_model, model, solve_model!, solve_linkage!, run_linkage!, run_model!, print_model_diagnostics, print_solver_status,
       diagnose_model, print_equation_diagnostics, diagnostic_variables, diagnostic_constraints,
       diagnostic_equation_matches, diagnostic_variables_without_equations,
       diagnostic_variables_with_multiple_equations, diagnostic_redundant_equations,
       run_recursive_dynamic!, update_period_data!, plot_dynamic_results,
       _collect_period_values, _period_summary,
       Scenario, write_policy_template, read_policy_scenarios, run_policy_experiments!,
       run_scenario!, update_period_data_scenario!, plot_all_scenarios

end
