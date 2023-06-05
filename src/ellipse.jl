function conic_to_axes(coeffs, normalise = true)
    # https://www.wikiwand.com/en/Ellipse#/General_ellipse
    # Extract  conic parameters
    if normalise && coeffs[1] != 0
        coeffs = coeffs/coeffs[1]
    end

    A, B, C, D, E, F = coeffs

    # Usefult intermediates
    disc = B * B - 4 * A * C

    # Center
    x0 = (2 * C * D - B * E) / disc
    y0 = (2 * A * E - B * D) / disc

    # Axis lengths
    term1 = 2 * (A * E * E + C * D * D - B * D * E + disc * F)
    a = - sqrt(term1 * (A + C + sqrt((A - C) ^ 2 + B ^ 2))) / disc
    b = - sqrt(term1 * (A + C - sqrt((A - C) ^ 2 + B ^ 2))) / disc

    if A < C
        if B == 0
            phi_b = 0
        else
            phi_b = atan(B/(A-C)) / 2
        end
    else
        if B == 0
            phi_b = π/2
        else
            phi_b = atan(B/(A-C)) / 2 - π/2
        end
    end


    return x0, y0, a, b, phi_b

end

"""
    fit_ellipse(positions::Vector{Point2f})

Document this function
"""
function fit_ellipse(positions)
    zzz = point_conic_eq.(to_value(positions))
    b= -ones(length(to_value(positions)))
    A = hcat(zzz...)'
    coeffs =  A  \ b
    return vcat(coeffs,[1])
end  # function fit_ellipse

function point_conic_eq(p::Point2f)
    x,y  = p.data
    return [x^2, x*y, y^2, x, y]
end


function getellipsepoints(cx, cy, rx, ry, θ)
	t = range(0, 2*pi, length=100)
	ellipse_x_r = @. rx * cos(t)
	ellipse_y_r = @. ry * sin(t)
	R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
	r_ellipse = [ellipse_x_r ellipse_y_r] * R
	x = @. cx + r_ellipse[:,1]
	y = @. cy + r_ellipse[:,2]
	(x,y)
end

