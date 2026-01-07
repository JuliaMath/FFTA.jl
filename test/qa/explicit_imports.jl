using FFTA, Test
import ExplicitImports

@testset "ExplicitImports" begin
    # No implicit imports in FFTA (ie. no `using MyPkg`)
    @test ExplicitImports.check_no_implicit_imports(FFTA) === nothing

    # No non-owning imports in FFTA (ie. no `using LinearAlgebra: map`)
    @test ExplicitImports.check_all_explicit_imports_via_owners(FFTA) === nothing

    # No non-public imports in FFTA (ie. no `using MyPkg: _non_public_internal_func`)
    @test ExplicitImports.check_all_explicit_imports_are_public(FFTA) === nothing

    # No stale imports in FFTA (ie. no `using MyPkg: func` where `func` is not used in FFTA)
    @test ExplicitImports.check_no_stale_explicit_imports(FFTA) === nothing

    # No non-owning accesses in FFTA (ie. no `... LinearAlgebra.map(...)`)
    @test ExplicitImports.check_all_qualified_accesses_via_owners(FFTA) === nothing

    # No non-public accesses in FFTA (ie. no `... MyPkg._non_public_internal_func(...)`)
    @test ExplicitImports.check_all_qualified_accesses_are_public(FFTA; ignore = (:require_one_based_indexing, :Fix1)) === nothing

    # No self-qualified accesses in FFTA (ie. no `... FFTA.func(...)`)
    @test ExplicitImports.check_no_self_qualified_accesses(FFTA) === nothing
end
