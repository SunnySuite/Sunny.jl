struct BroadeningDevice{F} <: Sunny.AbstractBroadening
    kernel :: F   # Function mapping x = (ω - ϵ) to an intensity scaling factor
    fwhm :: Float64
end

BroadeningDevice(host::Sunny.Broadening) = BroadeningDevice(host.kernel, host.fwhm)

function Adapt.adapt_structure(to, data::BroadeningDevice)
    kernel = Adapt.adapt_structure(to, data.kernel)
    fwhm = Adapt.adapt_structure(to, data.fwhm)
    BroadeningDevice(kernel, fwhm)
end

function (b::BroadeningDevice)(ϵ, ω)
    b.kernel(ω - ϵ)
end

function _broaden(data, bands_data, disp, energies, kernel)
    iq = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if iq > size(disp, 2)
        return
    end
    iω = threadIdx().y + (blockIdx().y - 1) * blockDim().y
    if iω > length(energies)
        return
    end
    ω = energies[iω]
    total = 0.
    for (ib, b) in enumerate(view(disp, :, iq))
        #norm(bands.data[ib, iq]) < cutoff && continue
        total += kernel(b, ω) * bands_data[ib, iq]
    end
    data[iω, iq] = total 
    return
end

function broaden!(data::CuArray{Ret}, bands::BandIntensitiesDevice{Ret}; energies, kernel) where Ret
    energies_d = CuArray(collect(Float64, energies))
    #issorted(energies) || error("energies must be sorted")

    nω = length(energies)
    nq = size(bands.qpts.qs,1)
    (nω, nq...) == size(data) || error("Argument data must have size ($nω×$(sizestr(bands.qpts)))")

    #asdf = norm.(vec(bands.data))
    #cutoff = 1e-12 * Statistics.quantile(asdf, 0.95)

    kernel_d = BroadeningDevice(kernel)
    gpu_kernel = CUDA.@cuda launch=false _broaden(data, bands.data, bands.disp, energies_d, kernel_d)
    config = launch_configuration(gpu_kernel.fun)
    optimal_threads_1d = config.threads
    
    threads_x = floor(Int, sqrt(optimal_threads_1d))
    threads_x = Base.min(nq, threads_x)

    threads_y = optimal_threads_1d ÷ threads_x
    threads_y = Base.min(nω, threads_y)

    threads = (threads_x, threads_y) # e.g., (16, 32) or similar, max product is 1024

    blocks_x = cld(nq, threads_x)
    blocks_y = cld(nω, threads_y)
    blocks = (blocks_x, blocks_y)
    gpu_kernel(data, bands.data, bands.disp, energies_d, kernel_d; threads=threads, blocks=blocks)

    return data
end

function Sunny.broaden(bands::BandIntensitiesDevice; energies, kernel)
    data = CUDA.zeros(eltype(bands.data), length(energies), size(bands.qpts.qs)...)
    broaden!(data, bands; energies, kernel)
    return IntensitiesDevice(bands.crystal, bands.qpts, collect(Float64, energies), data)
end
