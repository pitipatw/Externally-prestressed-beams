function get_fps(fₚₑ::Number, fc::Number , ρ::Number , fₚᵤ::Number) ::Number
    temp1 = fₚₑ + 70 + fc/(100*ρ)
    temp2 = fₚₑ + 420
    temp3 = 1300
    fₚₛ = minimum([temp1,temp2,temp3])
    if fₚₑ>=0.5*fₚᵤ && fₚₛ==temp1
        println("Invalid Assumption for Eq.1")
        return 0
    end
return fₚₛ
end

