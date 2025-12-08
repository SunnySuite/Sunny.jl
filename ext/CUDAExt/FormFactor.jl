struct FormFactorDevice
    j0 :: Sunny.ExpandedBesselIntegral
    j2 :: Sunny.ExpandedBesselIntegral
    g  :: Float64
end

FormFactorDevice(host::Sunny.FormFactor) = FormFactorDevice(host.j0, host.j2, host.g)

function Adapt.adapt_structure(to, data::FormFactorDevice)
    j0 = Adapt.adapt_structure(to, data.j0)
    j2 = Adapt.adapt_structure(to, data.j2)
    g = Adapt.adapt_structure(to, data.g)
    FormFactorDevice(j0, j2, g)
end
