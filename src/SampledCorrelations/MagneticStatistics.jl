mutable struct StaticMagneticStatistics
    n      :: Int64
    μ      :: Vec3
    Caccum :: Mat3
end

StaticMagneticStatistics() = StaticMagneticStatistics(0, zero(Vec3), zero(Mat3))

function add_magnetic_sample!(sms::StaticMagneticStatistics, sys)
    m = magnetic_moment_per_site(sys) # Could use q=0 information, but this guarantees correct observables and application of g-factor
    sms.n += 1
    δ₁ = m - sms.μ
    sms.μ += δ₁/sms.n
    δ₂ = m - sms.μ
    sms.Caccum += δ₁ * δ₂'

    return nothing
end

magnetic_moment_per_site(sms::StaticMagneticStatistics) = sms.μ

function magnetic_susceptibility_per_site(sms::StaticMagneticStatistics, kT)
    (; Caccum, n) = sms 
    return (1/kT) * (Caccum + Caccum') / 2*(n-1)
end