module EdesVersions

    module Maas15
        include("Maas15.jl")
    end 
    EDES(::Type{Val{:Maas15}}) = Maas15.EDES_Maas15

    # module Macsbio
    # include("Macsbio.jl")
    # end
    # EDES(::Type{Val{:Macsbio}})= Macsbio.EDES_Macsbio

    function EDES(::Type{Val{version}}) where {version}
        println("No version $version known, obtaining Maas15 settings")
        EDES(:Maas15)
    end

    EDES(version::Symbol) = EDES(Val{version})

    export EDES
end

