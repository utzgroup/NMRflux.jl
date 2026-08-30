#!/usr/bin/env julia
# -*- coding: utf-8 -*-
#
# Convert a Jupyter notebook back into a Documenter.jl markdown page — the
# inverse of md2ipynb.jl.
#
# Markdown cells are emitted as prose. `$$...$$` blocks that stand alone on
# their own lines are turned back into ```math``` fences. Code cells are
# wrapped in ```@example <name>``` fences, all sharing one session name (so
# Documenter runs them in one continuous session, matching how a notebook's
# cells share one kernel). Cell `outputs` are ignored: Documenter regenerates
# those itself when it builds the docs, so any outputs recorded in the
# notebook would be redundant at best, stale at worst.
#
# Usage:
#   julia docs/ipynb2md.jl docs/src/DataProcessing.ipynb
#   julia docs/ipynb2md.jl notebook.ipynb docs/src/Out.md
#
# Limitations (simple, line-based converter): raw cells are passed through
# as-is; `$$...$$` blocks are only recognised when each `$$` sits alone on
# its own line (the form md2ipynb.jl produces) — inline math and other
# LaTeX delimiters are left untouched.

using JSON

"Cell `source` is either a JSON array of lines or a single string; normalise to one string."
cell_text(c) = (s = c["source"]; s isa AbstractVector ? join(s) : s)

function example_name(outpath::AbstractString)
    base = first(splitext(basename(outpath)))
    slug = replace(base, r"[^A-Za-z0-9]" => "")
    isempty(slug) && (slug = "Doc")
    occursin(r"^[0-9]", slug) && (slug = "Doc" * slug)
    return slug * "Eg"
end

"Turn `\$\$...\$\$` blocks that stand alone on their own lines back into ```math``` fences."
function restore_math_fences(text::AbstractString)
    lines = split(text, '\n')
    out = IOBuffer()
    i, n = 1, length(lines)
    while i <= n
        if strip(lines[i]) == "\$\$"
            j = i + 1
            body = String[]
            while j <= n && strip(lines[j]) != "\$\$"
                push!(body, lines[j])
                j += 1
            end
            if j <= n   # found a closing $$
                println(out, "```math")
                for bl in body
                    println(out, bl)
                end
                println(out, "```")
                i = j + 1
                continue
            end
        end
        println(out, lines[i])
        i += 1
    end
    return String(take!(out))
end

function cell_markdown(c::AbstractDict, egname::AbstractString)
    text = rstrip(cell_text(c))
    isempty(text) && return ""
    ct = c["cell_type"]
    if ct == "markdown" || ct == "raw"
        return rstrip(restore_math_fences(text))
    elseif ct == "code"
        return "```@example $egname\n$text\n```"
    else
        error("unknown cell_type: $ct")
    end
end

function notebook_to_md(nb::AbstractDict, egname::AbstractString)
    pieces = [cell_markdown(c, egname) for c in nb["cells"]]
    filter!(!isempty, pieces)
    return join(pieces, "\n\n") * "\n"
end

function ipynb2md(inpath::AbstractString, outpath::AbstractString)
    nb = JSON.parsefile(inpath)
    md = notebook_to_md(nb, example_name(outpath))
    write(outpath, md)
    return outpath
end

function main()
    if isempty(ARGS)
        println("Usage: julia docs/ipynb2md.jl <input.ipynb> [output.md]")
        exit(1)
    end
    inpath = ARGS[1]
    outpath = length(ARGS) >= 2 ? ARGS[2] : replace(inpath, r"\.ipynb$" => ".md")
    ipynb2md(inpath, outpath)
    println("Wrote $outpath")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
