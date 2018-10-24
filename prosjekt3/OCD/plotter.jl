import PyPlot
const plt = PyPlot

function plottify(poss, dims = 3, labels = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptun", "Pluto"])
    labels = split(labels, r"'|[|]|,| ")
    filter!(e->e≠"",labels)
    filter!(e->e≠"[",labels)
    filter!(e->e≠"]",labels)
    for i = 1:Int64(length(labels))
        if dims == 3
            plt.plot3D(poss[i, 1, :], poss[i, 2, :], poss[i, 3, :], label = labels[i])
        elseif dims == 2
            plt.plot(poss[i, 1, :], poss[i, 2, :], label = labels[i])
        end
    end
    if dims == 3
        plt.zlabel("pos z [AU]")
    end
    plt.xlabel("pos x [AU]")
    plt.ylabel("pos y [AU]")
    plt.axis("equal")
    plt.legend()
    plt.show()
end