#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#
# Build a single docs/src page for fast debugging, instead of the whole
# manual. Runs the page's @example blocks and renders just that page.
#
# Usage:
#   julia docs/make_debug.jl DataProcessing.md
#
# Output goes to docs/build_debug/ (separate from the real docs/build/,
# so this doesn't clobber a full build you might already have).
#
# Note: pages are built in isolation, so any @ref links to other manual
# pages won't resolve — that's expected and downgraded to a warning here.
# To achieve that isolation, the target page is copied into a scratch
# source directory containing nothing else: Documenter renders every .md
# file it finds under `source` regardless of what's listed in `pages` (that
# keyword only controls the navigation menu), so pointing it at the full
# docs/src would silently rebuild the whole site every time.

ENV["GKSwstype"] = "100"

using Pkg
Pkg.activate(@__DIR__)   # docs/Project.toml: Documenter, Plots, NMRflux via [sources]
Pkg.instantiate()

using NMRflux
using Documenter

page = length(ARGS) >= 1 ? ARGS[1] : "DataProcessing.md"

# The scratch source dir must live inside the git repo (not the system temp
# dir): Documenter computes "Edit on GitHub" links relative to the repo
# root, and throws MissingRemoteError if the page it's rendering lives
# outside it.
mktempdir(@__DIR__) do tmpsrc
    cp(joinpath(@__DIR__, "src", page), joinpath(tmpsrc, page))

    makedocs(
        sitename  = "NMRflux.jl (debug)",
        modules   = [NMRflux],
        clean     = true,
        source    = tmpsrc,
        build     = joinpath(@__DIR__, "build_debug"),
        format    = Documenter.HTML(prettyurls = false, collapselevel = 2, assets = String[]),
        pages     = [page],
        checkdocs = :none,
        warnonly  = true,
    )
end

println("\nBuilt to docs/build_debug/$(splitext(page)[1]).html")
