module Utils

using GLMakie

export set_project_theme!

function set_project_theme!()
    set_theme!(theme_black())
    update_theme!(
        fontsize = 20,
        Axis = (
            xgridvisible = false,
            ygridvisible = false,
        )
    )
end

end
