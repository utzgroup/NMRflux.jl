#!/usr/bin/env julia
# -*- coding: utf-8 -*-

###########################################################################
# File:        make.jl
# Project:     NMRflux.jl
#
# Description:
#   Build the NMRflux.jl documentation using Documenter.jl.
#   Configures the documentation pages and HTML output and runs the
#   documentation build.
###########################################################################

ENV["GKSwstype"] = "100"   # force GR (Plots.jl backend) into headless mode, so
                            # rendering plots for the docs doesn't pop up a GUI window

using Pkg
Pkg.activate(@__DIR__)   # activate docs/Project.toml, which has Documenter, Plots (for
                         # the @example blocks that plot), and NMRflux itself via
                         # [sources]. Lets `julia docs/make.jl` work from any cwd.
Pkg.instantiate()

using NMRflux
using Documenter
using Documenter.Remotes

# -------------------------------------------------------------------------
# Configuration: define project metadata and documentation settings
# -------------------------------------------------------------------------

makedocs(
    sitename = "NMRflux.jl",
    authors  = "Manaz Kaleel & Marcel Utz",
    clean    = true,

    repo = Remotes.GitHub("utzgroup", "NMRflux.jl"),

    modules = [NMRflux],

    format = Documenter.HTML(
        prettyurls   = false,   # set to true later for GitHub Pages
        collapselevel = 2,
        #assets       = String[],        # Documenter default but uses only half a page
        assets = ["custom.css"],
    ),

    pages = [
        "Home" => "index.md",

        "Quick Start" => [
            "Getting Started" => "QuickStart.md",
        ],

        "Manual" => [
            "User Manual" => "Manual.md",
        ],

        "Advanced topics" => [
            "Data Loading" => "DataLoading.md",
            "SpectData" => "SpectData.md",
            "Data Processing" => "DataProcessing.md",
            "Spin Dynamics and FID generation" => "SpinDynamics.md",
            # "RINSE" => "RINSE.md",
            "Sinusoidal Decomposition (Craft)" => "Craft.md",
        ],

        "Development" => [
            "Roadmap to 1.0" => "NMRflux-1.0.md",
        ],

        "Reference" => [
            "API" => "Reference.md",
        ],

        "Machine learning demonstrator" => [
            "RINSE" => "RINSE.md",
        ],
        "Release Notes" => "ReleaseNotes.md",
    ],

    # For now, don't fail if some docstrings are not yet included in the manual
    checkdocs = :none,
)

# ------------------------------------------------------------------------------
# Deploy docs needed when pushing to GitHub Pages
# Uncomment when GitHub CI is set up:
# deploydocs(
#     repo   = "https://github.com/utzgroup/NMRflux.jl.git",
#     target = "build",
# )
# ------------------------------------------------------------------------------
