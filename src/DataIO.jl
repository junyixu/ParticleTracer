"""
The DataIO module handles input/output operations for particle data.
Main functionalities include:
- Batch saving of large-scale particle data
- Creating data indices
- Reading particle data
- Managing data storage structure
"""
module DataIO

using HDF5
using Dates
using Distributed
using ..PtcStruct
using ..UserInputs
using ..UserInputs: TotalSteps, SavePerNSteps
using TimeZones  # 需要添加这个包

export save_particle_data, read_particle_data, create_index_file, merge_process_files, SaveConfig

"""
    SaveConfigType

Configuration type for data saving operations.

Fields:
- batch_size::Int: Number of particles in each batch
- output_dir::String: Output directory path
- timestamp::String: Timestamp of data saving operation
"""
const SaveConfigType = NamedTuple{(:batch_size, :output_dir, :timestamp),
                            Tuple{Int, String, String}}

"""
    SaveConfig(; batch_size=1000)

Constructor for SaveConfig instance.

Parameters:
- batch_size: Number of particles per batch, defaults to 1000

Returns:
- SaveConfig instance
"""
function SaveConfig(; batch_size=1000)
    output_dir = if isdefined(UserInputs, :output_dir) 
        UserInputs.output_dir 
    else 
        "./DataAnalysis"
    end
    mkpath(output_dir)
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    
    # Create a SaveConfigType configuration object
    # batch_size: Number of particles per batch
    # output_dir: Output directory path 
    # timestamp: Timestamp of data saving operation
    return (batch_size=batch_size, output_dir=output_dir, timestamp=timestamp)
end

"""
    save_particle_data(ptc_data, n, x0, p0, B0, config)

Save particle data using process-wise batch strategy. Each process writes to its own
temporary file, which will be merged later.

Parameters:
- ptc_data: Particle trajectory data
- n: Global particle number
- x0: Initial position
- p0: Initial momentum
- B0: Initial magnetic field
- config: Save configuration
"""
function save_particle_data(ptc_data::ParticleData, n::Int, x0, p0, B0, config::SaveConfigType)
    batch_number = ceil(Int, n/config.batch_size)
    # Each process writes to its own temporary file
    filename = joinpath(config.output_dir, 
                       "particles_$(config.timestamp)_batch$(batch_number)_proc$(myid()).h5")
    
    h5open(filename, "cw") do file
        particle_id = mod1(n, config.batch_size)
        save_single_particle(file, particle_id, ptc_data, n, x0, p0, B0)
        update_batch_metadata!(file, batch_number, config)
    end
end

"""
    save_single_particle(file, particle_id, ptc_data, global_n, x0, p0, B0)

Save single particle data to HDF5 file.

Parameters:
- file: Open HDF5 file
- particle_id: Particle ID within current batch
- ptc_data: Particle trajectory data
- global_n: Global particle number
- x0, p0, B0: Initial conditions
"""
function save_single_particle(file, particle_id, ptc_data, global_n, x0, p0, B0)
    g_ptc = create_group(file, "particle_$particle_id")
    
    # Save trajectory data
    g_traj = create_group(g_ptc, "trajectory")
    g_traj["position"] = ptc_data.X
    g_traj["momentum"] = ptc_data.P
    g_traj["magnetic_field"] = ptc_data.B
    
    # Save metadata
    g_meta = create_group(g_ptc, "metadata")
    attrs(g_meta)["total_steps"] = TotalSteps
    attrs(g_meta)["save_per_n_steps"] = SavePerNSteps
    attrs(g_meta)["global_particle_number"] = global_n
    attrs(g_meta)["creation_date"] = string(now(localzone()))
    attrs(g_meta)["dt"] = UserInputs.Δt
    
    # Save initial conditions
    g_meta["initial_position"] = x0
    g_meta["initial_momentum"] = p0
    g_meta["initial_B_field"] = B0
end

"""
    update_batch_metadata!(file, batch_number, config)

Update or create global metadata for batch file.

Parameters:
- file: Open HDF5 file
- batch_number: Current batch number
- config: Save configuration
"""
function update_batch_metadata!(file, batch_number, config::SaveConfigType)
    if !haskey(file, "global_metadata")
        g_global = create_group(file, "global_metadata")
        attrs(g_global)["batch_number"] = batch_number
        attrs(g_global)["batch_size"] = config.batch_size
        attrs(g_global)["total_particles"] = UserInputs.N
        attrs(g_global)["creation_date"] = string(Dates.now())
    end
end

"""
    create_index_file(config)

Create index file recording information for all batches.

Parameters:
- config: Save configuration

Index file contains:
- Total number of particles
- Batch size
- Total number of batches
- Creation time
- Simulation timestamp
"""
function create_index_file(config::SaveConfigType)
    filename = joinpath(config.output_dir, "index_$(config.timestamp).h5")
    h5open(filename, "w") do file
        attrs(file)["total_particles"] = UserInputs.N
        attrs(file)["batch_size"] = config.batch_size
        attrs(file)["total_batches"] = ceil(Int, UserInputs.N/config.batch_size)
        attrs(file)["creation_date"] = string(Dates.now())
        attrs(file)["simulation_timestamp"] = config.timestamp
        attrs(file)["dt"] = UserInputs.Δt
    end
end

"""
    read_particle_data(n, timestamp)

Read data for specified particle.

Parameters:
- n: Global particle number
- timestamp: Data saving timestamp

Returns:
Dictionary containing position, momentum, magnetic field, and metadata
"""
function read_particle_data(n::Int, timestamp::String)
    config = SaveConfig()
    batch_number = ceil(Int, n/config.batch_size)
    particle_id = mod1(n, config.batch_size)
    
    filename = joinpath(config.output_dir, 
                       "particles_$(timestamp)_batch$(batch_number).h5")
    
    h5open(filename, "r") do file
        g_ptc = file["particle_$particle_id"]
        return Dict(
            "position" => read(g_ptc["trajectory/position"]),
            "momentum" => read(g_ptc["trajectory/momentum"]),
            "magnetic_field" => read(g_ptc["trajectory/magnetic_field"]),
            "metadata" => read_metadata(g_ptc["metadata"])
        )
    end
end

"""
    read_metadata(g_meta)

Read all metadata attributes from HDF5 group.

Parameters:
- g_meta: HDF5 metadata group

Returns:
Dictionary containing all metadata
"""
function read_metadata(g_meta)
    Dict(String(name) => read(attr) for (name, attr) in attrs(g_meta))
end

"""
    merge_process_files(config)

Merge all process temporary data files into final batch files.
Call this function after all computations are complete.

Parameters:
- config: Save configuration
"""
function merge_process_files(config::SaveConfigType)
    total_batches = ceil(Int, UserInputs.N/config.batch_size)
    
    for batch_number in 1:total_batches
        # Create final batch file
        final_filename = joinpath(config.output_dir, 
                                "particles_$(config.timestamp)_batch$(batch_number).h5")
        
        h5open(final_filename, "w") do final_file
            # Merge all process files
            for proc_id in workers()
                proc_filename = joinpath(config.output_dir, 
                                      "particles_$(config.timestamp)_batch$(batch_number)_proc$(proc_id).h5")
                
                if isfile(proc_filename)
                    h5open(proc_filename, "r") do proc_file
                        # 只复制 particle_i 的数据
                        filter(particle_id -> startswith(particle_id, "particle_"), keys(proc_file)) |>
                            particle_ids -> foreach(particle_id -> copy_object(proc_file[particle_id], final_file, particle_id), particle_ids)
                    end
                    # Delete temporary file
                    rm(proc_filename)
                end
            end
            
            # Update final file metadata
            update_batch_metadata!(final_file, batch_number, config)
        end
    end
end

end # module 
