# this script plot pairplots of posteriors predicted by Neural likelihood estimator 
#
# Kai Li
# 24 May 2025

using CSV, Tables
using DataFrames

using CairoMakie
using PairPlots

galaxy = CSV.read("samples_30May2025_galaxy_final_time_v2.csv", DataFrame; delim = ',', header = false)

rename!(galaxy, ["nStar", "Pps", "Psp", "gamma", "pa","angle_prolif"]) 

fig = Figure()

fig = pairplot(
    PairPlots.Series( galaxy[10000:100:end,:], color=Makie.RGBA(0.12, 0.47, 0.71, 0.5)) => (
        PairPlots.HexBin(colormap=Makie.cgrad([:transparent, Makie.RGB(0.12, 0.47, 0.71)]),),
        PairPlots.Scatter(filtersigma=2),
        PairPlots.Contour(color=Makie.RGB(0.12, 0.47, 0.71)),
        PairPlots.MarginDensity(
            color=Makie.RGB(0.12, 0.47, 0.71),
            linewidth=1.5f0
        ),
        PairPlots.MarginQuantileText()
    ),
    labels = Dict(
    # basic string
    :nStar => L"n^*",
    # Makie rich text
    :Pps => L"P_{ps}",
    # LaTeX String
    :Psp => L"P_{sp}",
    :gamma => L"\gamma",
    :pa => L"p_a",
    :angle_prolif => L"\theta",
    ),
    axis=(;
        nStar = (; lims = (;low=0, high=1)),
        Pps = (; lims = (;low=0, high=1)),
        Psp = (; lims = (;low=0, high=1)),
        gamma = (; lims = (;low=0, high=1)),
        pa = (; lims = (;low=0, high=1)),
        angle_prolif = (; lims = (;low=0, high=pi/20))
    ),
)

# save("pairplot_galaxy_30May2025.pdf",fig)


fig = pairplot(galaxy[10000:100:end,:],
    labels = Dict(
    # basic string
    :nStar => L"n^*",
    # Makie rich text
    :Pps => L"P_{ps}",
    # LaTeX String
    :Psp => L"P_{sp}",
    :gamma => L"\gamma",
    :pa => L"p_a",
    :angle_prolif => L"\theta",
    ),
    axis=(;
        nStar = (; lims = (;low=0, high=1)),
        Pps = (; lims = (;low=0, high=1)),
        Psp = (; lims = (;low=0, high=1)),
        gamma = (; lims = (;low=0, high=1)),
        pa = (; lims = (;low=0, high=1)),
        angle_prolif = (; lims = (;low=0, high=pi/20))
    ),
)
# save("pairplot_galaxy_30May2025.pdf",fig)