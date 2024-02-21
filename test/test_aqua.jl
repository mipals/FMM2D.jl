using Test
using Aqua
using FMM2D

Aqua.test_undefined_exports(FMM2D)
Aqua.test_project_extras(FMM2D)
Aqua.test_unbound_args(FMM2D)
Aqua.test_ambiguities(FMM2D)
Aqua.test_deps_compat(FMM2D)
Aqua.test_stale_deps(FMM2D)
Aqua.test_piracies(FMM2D)
Aqua.test_persistent_tasks(FMM2D)
