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

ENV["GKSwstype"] = "100"

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NMRflux
using Documenter

page = length(ARGS) >= 1 ? ARGS[1] : "DataProcessing.md"

makedocs(
    sitename  = "NMRflux.jl (debug)",
    modules   = [NMRflux],
    clean     = true,
    build     = joinpath(@__DIR__, "build_debug"),
    format    = Documenter.HTML(prettyurls = false, collapselevel = 2, assets = String[]),
    pages     = [page],
    checkdocs = :none,
    warnonly  = true,
)

println("\nBuilt to docs/build_debug/$(splitext(page)[1]).html")
