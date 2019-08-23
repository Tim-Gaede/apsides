#───────────────────────────────────────────────────────────────────────────────
# Timothy Gaede 2019-08-23
#
# NOTE: GM (μ for short) of a celestial body is typically known with more
# relative ("percent-wise") precision than M.  The GM of a planetary
# body is inferred more directly than M.
#
#
# Returns periapsis and apoapsis altitudes (i.e. distances from the surface).
function apsides(pos, vel, GM, R)

    peri, apo = apsides(pos, vel, GM)

    peri - R,    apo - R
end


# Returns periapsis and apoapsis distances from center of the celestial body.
function apsides(pos, vel, GM)

    function mag(a)
        ∑aᵢ² = 0.0
        for aᵢ in a;    ∑aᵢ² += aᵢ^2; end
        √∑aᵢ²
    end

    # WARNING: length(vec) would return the number of items in the array.
    # A "3D vector" in physics terminology can be represented as a 1D vector
    # with three items in computer terminology.
    r = mag(pos)
    v = mag(vel)

    h = r*v # angular momentum per unit mass at periapsis and apoapsis
    ϵ = (v^2)/2 - GM/r # orbital energy per unit mass

    # Substituting v = h/r into ϵ = (v^2)/2 - GM/r yields:
    # ϵ = (h^2)/2r^2 - GM/r  
    # r ≠ 0.  Rearrange into a quadratic equation to solve for r:
    # ϵr^2 + GMr - (h^2)/2 = 0   
    a = ϵ
    b = GM
    c = -(h^2)/2


    ((-b + √(b^2 - 4a*c)) / 2a),    ((-b - √(b^2 - 4a*c)) / 2a)
end
#───────────────────────────────────────────────────────────────────────────────



using Formatting

#═══════════════════════════════════════════════════════════════════════════════
function main()
    println("-"^50, "\n"^2)
    G = 6.673 * 10.0^-11
    M = 5.972 * 10.0^24 # Integer rollover if using 10^24
    k = 1000
    R = 6378k

    pos = [R+200k, 0, 0]
    vel = [0, 7790, 0]

    peri, apo = apsides(pos, vel, G*M)

    peri_km_text = format(peri / k, precision = 2) * "km"
    apo_km_text  = format(apo / k, precision = 2) * "km"
    println("Distances from the center of the celestial body:")
    println("periapsis: ", peri_km_text,  "     apoapsis: ", apo_km_text)


    perialt, apoalt = apsides(pos, vel, G*M, R)
    perialt_km_text = format(perialt / k, precision = 2) * "km"
    apoalt_km_text  = format(apoalt / k, precision = 2) * "km"
    println("\n"^2)
    println("altitudes (distance above the surface):")
    println("periapsis: ", perialt_km_text,  "     apoapsis: ", apoalt_km_text)

end
#═══════════════════════════════════════════════════════════════════════════════
main()
