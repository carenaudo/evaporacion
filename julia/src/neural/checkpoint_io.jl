"""
    Evaporacion.Neural.CheckpointIO

Serialization and weight conversion utilities between PyTorch NPZ dumps and Lux parameter structures.
"""
module CheckpointIO

using NPZ
using Lux
using ComponentArrays

export load_option3b_weights_from_npz, save_lux_weights_npz

"""
    load_option3b_weights_from_npz(npz_path::String) -> ComponentArray

Loads weights exported from PyTorch option3b checkpoint and maps them into Lux named tuples.
"""
function load_option3b_weights_from_npz(npz_path::String)
    dict = npzread(npz_path)

    # PyTorch linear weights are (out_features, in_features). Lux Dense weights are (out_features, in_features).
    ps = (
        encoder = (
            layer_1 = (weight = dict["encoder.net.0.weight"], bias = dict["encoder.net.0.bias"]),
            layer_2 = (weight = dict["encoder.net.2.weight"], bias = dict["encoder.net.2.bias"]),
            layer_3 = (weight = dict["encoder.net.4.weight"], bias = dict["encoder.net.4.bias"])
        ),
        rhs = (
            layer_1 = (weight = dict["rhs.net.0.weight"], bias = dict["rhs.net.0.bias"]),
            layer_2 = (weight = dict["rhs.net.2.weight"], bias = dict["rhs.net.2.bias"]),
            layer_3 = (weight = dict["rhs.net.4.weight"], bias = dict["rhs.net.4.bias"]),
            layer_4 = (weight = dict["rhs.net.6.weight"], bias = dict["rhs.net.6.bias"])
        )
    )

    return ps
end

"""
    save_lux_weights_npz(ps, npz_path::String)

Saves Lux parameter hierarchy into an NPZ archive.
"""
function save_lux_weights_npz(ps, npz_path::String)
    dict = Dict{String, Any}()
    # Flatten parameter names
    function _flatten(prefix, obj)
        for (k, v) in pairs(obj)
            key = isempty(prefix) ? string(k) : "$(prefix).$(k)"
            if v isa NamedTuple
                _flatten(key, v)
            elseif v isa AbstractArray
                dict[key] = Array(v)
            end
        end
    end
    _flatten("", ps)
    npzwrite(npz_path, dict)
    return nothing
end

end # module CheckpointIO
