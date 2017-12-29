#=============================#
include("atoms.jl");
include("cell.jl");


function mergefiles(infile::AbstractArray{String}, outfile::String)
  out=open(outfile,"w")
  # Loop over the files
  for i=1:size(infile)[1]
    # Reading file
    open(infile[i],"r") do f
      for line in eachline(f)
        # Writting line into merger
        write(out,line)
      end
    end
  end
  close(out)
end

include("PDB.jl");
