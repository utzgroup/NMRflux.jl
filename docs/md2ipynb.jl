#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#
# Convert a Documenter.jl markdown page into a Jupyter notebook.
#
# Prose becomes markdown cells. Code fences containing actual Julia code
# (```@example ...```, ```julia```, ```@julia```) become Julia code cells,
# in document order, so they execute top-to-bottom in a single kernel just
# like a Documenter `@example` session does. ```math``` fences are folded
# into the surrounding markdown as `$$...$$` so they render as LaTeX in
# Jupyter. Everything else (```@docs```, ```@meta```, ```toml```, ...) is
# left untouched as a literal fenced block inside its markdown cell, since
# those are Documenter build directives, not runnable Julia.
#
# Usage:
#   julia docs/md2ipynb.jl docs/src/DataProcessing.md
#   julia docs/md2ipynb.jl docs/src/DataProcessing.md /tmp/out.ipynb
#
# Limitations (this is a simple, line-based converter, not a full CommonMark
# parser): fences must start at column 1 with exactly three backticks and no
# indentation, and fenced blocks may not contain nested triple-backtick
# fences. `![](file.svg)` references from `@example` blocks are left as-is;
# they will only resolve once the notebook is actually run and regenerates
# those files (e.g. via `savefig(...)`) next to itself.

using JSON

const FENCE = r"^```(.*)$"

function fence_kind(info::AbstractString)
    tok = isempty(strip(info)) ? "" : first(split(strip(info)))
    if tok == "@example" || tok == "julia" || tok == "@julia"
        return :code
    elseif tok == "math"
        return :math
    else
        return :verbatim
    end
end

"Split text into lines, each keeping its trailing newline (nbformat convention), except the last."
function nb_source_lines(text::AbstractString)
    isempty(text) && return String[]
    lines = split(text, '\n'; keepempty=true)
    # split on '\n' introduces one trailing empty string if text ends with '\n'; drop it
    if lines[end] == ""
        lines = lines[1:end-1]
    end
    return [i < length(lines) ? lines[i] * "\n" : lines[i] for i in eachindex(lines)]
end

markdown_cell(text::AbstractString, id::AbstractString) = Dict(
    "cell_type" => "markdown",
    "id"        => id,
    "metadata"  => Dict(),
    "source"    => nb_source_lines(rstrip(text)),
)

code_cell(text::AbstractString, id::AbstractString) = Dict(
    "cell_type"         => "code",
    "id"                => id,
    "metadata"          => Dict(),
    "execution_count"   => nothing,
    "outputs"           => [],
    "source"            => nb_source_lines(rstrip(text)),
)

function parse_md_to_cells(text::AbstractString)
    cells = Any[]
    md_buf = IOBuffer()
    counter = Ref(0)
    next_id() = (counter[] += 1; string(counter[]; base=16, pad=8))

    function flush_markdown!()
        s = String(take!(md_buf))
        if !isempty(strip(s))
            push!(cells, markdown_cell(s, next_id()))
        end
        md_buf = IOBuffer()
    end

    lines = split(text, '\n')
    i = 1
    n = length(lines)
    while i <= n
        line = lines[i]
        m = match(FENCE, line)
        if m === nothing
            println(md_buf, line)
            i += 1
            continue
        end

        info = m.captures[1]
        kind = fence_kind(info)

        # collect body lines up to the closing fence (a bare ``` on its own line)
        body = String[]
        j = i + 1
        while j <= n && strip(lines[j]) != "```"
            push!(body, lines[j])
            j += 1
        end
        body_text = join(body, "\n")

        if kind == :code
            flush_markdown!()
            push!(cells, code_cell(body_text, next_id()))
        elseif kind == :math
            println(md_buf, "\$\$")
            println(md_buf, body_text)
            println(md_buf, "\$\$")
        else
            # verbatim: keep the whole fence (open/body/close) as literal markdown text
            println(md_buf, line)
            for bl in body
                println(md_buf, bl)
            end
            println(md_buf, "```")
        end

        i = j + 1   # skip past the closing fence line
    end

    flush_markdown!()
    return cells
end

function notebook(cells)
    return Dict(
        "cells"    => cells,
        "metadata" => Dict(
            "kernelspec" => Dict(
                "display_name" => "Julia $(VERSION)",
                "language"     => "julia",
                "name"         => "julia-$(VERSION.major).$(VERSION.minor)",
            ),
            "language_info" => Dict(
                "name"    => "julia",
                "version" => string(VERSION),
            ),
        ),
        "nbformat"        => 4,
        "nbformat_minor"  => 5,
    )
end

function md2ipynb(inpath::AbstractString, outpath::AbstractString)
    text = read(inpath, String)
    cells = parse_md_to_cells(text)
    nb = notebook(cells)
    open(outpath, "w") do io
        JSON.print(io, nb, 1)
    end
    return outpath
end

function main()
    if isempty(ARGS)
        println("Usage: julia docs/md2ipynb.jl <input.md> [output.ipynb]")
        exit(1)
    end
    inpath = ARGS[1]
    outpath = length(ARGS) >= 2 ? ARGS[2] : replace(inpath, r"\.md$" => ".ipynb")
    md2ipynb(inpath, outpath)
    println("Wrote $outpath")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
