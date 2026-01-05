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
    iω = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    if iω > length(energies)
        return
    end
    iq = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    if iq > size(disp, 2)
        return
    end

    bands_buf = CuDynamicSharedArray(Float64, (size(bands_data, 1), blockDim().y))
    bands_bufq = view(bands_buf, :, threadIdx().y)
    k = threadIdx().x
    while k <= size(bands_data, 1)
        bands_bufq[k] = bands_data[k, iq]
        k += blockDim().x
    end

    disp_buf = CuDynamicSharedArray(Float64, (size(disp, 1), blockDim().y), size(bands_data, 1) * blockDim().y* sizeof(Float64))
    disp_bufq = view(disp_buf, :, threadIdx().y)
    k = threadIdx().x
    while k <= size(disp, 1)
        disp_bufq[k] = disp[k, iq]
        k += blockDim().x
    end
    CUDA.sync_threads()

    ω = energies[iω]
    total = 0.
    for ib in eachindex(disp_bufq)
        total += kernel(disp_bufq[ib], ω) * bands_bufq[ib]
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

    function get_shmem(threads; rows=size(bands.data,1))
        if length(threads) == 2
            return 2 * threads[2] * rows * sizeof(Float64)
        else
            return 2 * threads * rows * sizeof(Float64)
        end
    end

    config = launch_configuration(gpu_kernel.fun, shmem=threads->get_shmem(threads))
    optimal_threads_1d = config.threads

    threads_x = Base.min(nω, optimal_threads_1d)
    threads_y = Base.min(nq, optimal_threads_1d ÷ threads_x)
    threads = (threads_x, threads_y) # e.g., (16, 32) or similar, max product is 1024

    blocks_x = cld(nω, threads_x)
    blocks_y = cld(nq, threads_y)
    blocks = (blocks_x, blocks_y)

    gpu_kernel(data, bands.data, bands.disp, energies_d, kernel_d; threads=threads, blocks=blocks, shmem=get_shmem(threads))
    return data
end

function Sunny.broaden(bands::BandIntensitiesDevice; energies, kernel)
    data = CUDA.zeros(eltype(bands.data), length(energies), size(bands.qpts.qs)...)
    broaden!(data, bands; energies, kernel)
    return IntensitiesDevice(bands.crystal, bands.qpts, collect(Float64, energies), data)
end
