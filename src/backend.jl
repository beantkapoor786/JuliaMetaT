using Preferences

const _BACKEND       = Ref{Symbol}(:unset)
const _CUDA_LOADED   = Ref(false)
const _AMDGPU_LOADED = Ref(false)
const _METAL_LOADED  = Ref(false)

function __init__()
    try; @eval Main using CUDA;   _CUDA_LOADED[]   = true; catch; end
    try; @eval Main using AMDGPU; _AMDGPU_LOADED[] = true; catch; end
    try; @eval Main using Metal;  _METAL_LOADED[]  = true; catch; end

    _BACKEND[] = _detect_active_backend()
end

"""
Read JACC's configured backend. Priority order:
1. `JACC_BACKEND` environment variable (allows SLURM scripts to set the backend
   without modifying LocalPreferences.toml in the project directory)
2. JACC's LocalPreferences.toml entry (`default_backend`)
3. Fall back to `threads`

Returns one of `:cuda`, `:amdgpu`, `:metal`, `:oneapi`, `:threads`.
"""
function _detect_active_backend()
    env_val = get(ENV, "JACC_BACKEND", "")
    if !isempty(env_val)
        return Symbol(lowercase(strip(env_val)))
    end
    s = try
        load_preference(Base.UUID("0979c8fe-16a4-4796-9b82-89a9f10403ea"),
                        "default_backend", "threads")
    catch
        "threads"
    end
    return Symbol(s)
end

"""
    set_backend(; force=:auto) -> Symbol

Reports the currently active JACC backend. JACC's backend is configured
at the project level via Preferences.jl in `LocalPreferences.toml` and
is fixed for the duration of the Julia session.

To change the backend, run (in a fresh Julia session):

    julia> using JACC
    julia> JACC.set_backend("threads")   # or "cuda", "amdgpu", "metal"

then restart Julia.

The `force` keyword is accepted for API stability but ignored — the
active backend cannot be changed at runtime.
"""
function set_backend(; force::Symbol = :auto)
    chosen = _BACKEND[]

    env_src = !isempty(get(ENV, "JACC_BACKEND", "")) ? " (from JACC_BACKEND env)" : ""
    msg = chosen === :cuda    ? "Backend: NVIDIA CUDA$env_src" :
          chosen === :amdgpu  ? "Backend: AMD ROCm$env_src"    :
          chosen === :metal   ? "Backend: Apple Metal$env_src" :
          chosen === :threads ? "Backend: CPU Threads$env_src" :
                                "Backend: $chosen$env_src"

    if force !== :auto && force !== chosen
        @warn "Requested backend `$force` differs from active backend `$chosen`. " *
              "JACC backend cannot be changed at runtime — edit LocalPreferences.toml " *
              "and restart Julia."
    end

    println(msg)
    return chosen
end

current_backend() = _BACKEND[]
