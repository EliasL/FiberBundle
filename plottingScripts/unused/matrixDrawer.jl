using Luxor
using JLD2
using MathTeXEngine
using LaTeXStrings
using Plots
using Random
using DataStructures
include("../support/dataManager.jl")

function drawmatrix(L=3, pixel_size = 1)
    tiles = Tiler(L*pixel_size*4, L*pixel_size*4, L, L, margin=0)
    for (pos, n) in tiles
        c_pos = pos .+ pixel_size/2
        box(pos, tiles.tilewidth*0.6, tiles.tileheight*0.6, :clip)
        background(0.8,0.8,0.8)
        clipreset()
        Luxor.text("$n", c_pos, halign=:right, valign=:center)
    end
end

function plotMatrix()
    L=3
    ps = 100 #Pixel size
    ly, lx = (2,1)
        
    pixel_L = L*ps
    title_space = ceil(Int, pixel_L/2)
    subtitle_space = ceil(Int, pixel_L/5)
    spacing_x = ceil(Int, pixel_L/20)
    spacing_y = spacing_x + subtitle_space
    image_size_x = pixel_L*lx+(lx-1)*spacing_x
    image_size_y = pixel_L*ly+(ly)*spacing_y + title_space
    font_size = 4*pixel_L/32
    title_font_size = font_size * 4/3

    Drawing(image_size_x, image_size_y, "plots/Visualizations/LocalLoadSharing.png")
    background(1,1,1) #White background
    fontface("Computer Modern")
    fontsize(title_font_size)
    Luxor.text(latexstring("Matrix"), image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
    fontsize(font_size)
    for j=1:ly, i=1:lx
        origin((pixel_L+spacing_x)*(i-1), (pixel_L+spacing_y)*(j-1) + title_space)
        drawmatrix(L, pixel_L/20)
        Luxor.text("test$j$i", pixel_L/2, pixel_L+subtitle_space, halign=:center, valign=:bottom)
    end

    finish()
end


plotMatrix()