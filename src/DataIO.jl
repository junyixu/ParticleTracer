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
export save_particle_data, read_particle_data, create_index_file, merge_process_files, SaveConfig, SaveConfigType

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

保存粒子数据，使用进程级批处理策略。每个进程写入自己的临时文件，之后合并。

Parameters:
- ptc_data: 粒子轨迹数据
- n: 全局粒子编号
- x0: 初始位置
- p0: 初始动量
- B0: 初始磁场
- config: 保存配置
"""
function save_particle_data(ptc_data::AbstractParticleData, n::Int, x0, p0, B0, config::SaveConfigType)
    batch_number = ceil(Int, n/config.batch_size)
    # Each process writes to its own temporary file
    filename = joinpath(config.output_dir, 
                       "particles_$(config.timestamp)_batch$(batch_number)_proc$(myid()).h5")
    
    # Open file in create/write mode:
    # 'c' - create file if it doesn't exist
    # 'w' - open for writing
    h5open(filename, "cw") do file
        # Using mod1 instead of mod because:
        # 1. mod1 returns values in range 1 to m (while mod returns 0 to m-1)
        # 2. We need particle IDs to start from 1 for HDF5 group names
        # 3. This maintains consistency with 1-based indexing convention
        # eg:
        # julia> mod1(5,5)
        # 5
        # julia> mod(5,5)
        # 0
        particle_id = mod1(n, config.batch_size)
        save_single_particle(file, particle_id, ptc_data, n, x0, p0, B0)
        update_batch_metadata!(file, batch_number, config)
    end
end

"""
    save_single_particle(file, particle_id, ptc_data, global_n, x0, p0, B0)

保存单个粒子数据到 HDF5 文件。
"""
function save_single_particle(file, particle_id, ptc_data::MagneticParticleData, global_n, x0, p0, B0)
    g_ptc = create_group(file, "particle_$particle_id")
    
    # 保存轨迹数据
    g_traj = create_group(g_ptc, "trajectory")
    g_traj["position"] = ptc_data.X
    g_traj["momentum"] = ptc_data.P
    g_traj["magnetic_field"] = ptc_data.B
    
    # 保存元数据
    save_metadata(g_ptc, global_n, x0, p0, B0)
end

function save_single_particle(file, particle_id, ptc_data::EMParticleData, global_n, x0, p0, B0)
    g_ptc = create_group(file, "particle_$particle_id")
    
    # 保存轨迹数据
    g_traj = create_group(g_ptc, "trajectory")
    g_traj["position"] = ptc_data.X
    g_traj["momentum"] = ptc_data.P
    g_traj["magnetic_field"] = ptc_data.B
    g_traj["electric_field"] = ptc_data.E
    
    # 保存元数据
    save_metadata(g_ptc, global_n, x0, p0, B0)
end

"""
    save_metadata(g_ptc, global_n, x0, p0, B0)

保存粒子的元数据。
"""
function save_metadata(g_ptc, global_n, x0, p0, B0)
    g_meta = create_group(g_ptc, "metadata")
    attrs(g_meta)["total_steps"] = TotalSteps
    attrs(g_meta)["save_per_n_steps"] = SavePerNSteps
    attrs(g_meta)["global_particle_number"] = global_n
    attrs(g_meta)["creation_date"] = string(now(localzone()))
    attrs(g_meta)["Δt"] = UserInputs.Δt
    
    # 保存初始条件
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
        attrs(file)["Δt"] = UserInputs.Δt
    end
end

"""
    read_particle_data(n, timestamp)

读取指定粒子的数据。

返回：
包含位置、动量、磁场（和电场，如果存在）以及元数据的字典
"""
function read_particle_data(n::Int, timestamp::String)
    config = SaveConfig()
    batch_number = ceil(Int, n/config.batch_size)
    particle_id = mod1(n, config.batch_size)
    
    filename = joinpath(config.output_dir, 
                       "particles_$(timestamp)_batch$(batch_number).h5")
    
    h5open(filename, "r") do file
        g_ptc = file["particle_$particle_id"]
        data = Dict(
            "position" => read(g_ptc["trajectory/position"]),
            "momentum" => read(g_ptc["trajectory/momentum"]),
            "magnetic_field" => read(g_ptc["trajectory/magnetic_field"]),
            "metadata" => read_metadata(g_ptc["metadata"])
        )
        
        # 如果存在电场数据，则添加到返回结果中
        if haskey(g_ptc["trajectory"], "electric_field")
            data["electric_field"] = read(g_ptc["trajectory/electric_field"])
        end
        
        return data
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
