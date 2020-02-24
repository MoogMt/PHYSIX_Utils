# Loading file
include("contactmatrix.jl")


Volumes=[8.6,8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500]
Cut_Off=[1.75]

folder_base=string("/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/")

for cut_off in Cut_Off
    for V in Volumes
        for T in Temperatures
            folder_local=string(folder_base,V,"/",T,"K/Data/")
            if isfile(string(folder_local,"dimer-",cut_off,".dat")) && isfile(string(folder_local,"molecules_allinfo_cutoff-",cut_off,".dat"))
                #----------------
                file_molecules=open(string(folder_local,"molecules_allinfo_cutoff-",cut_off,".dat"))
                lines_molecules=readlines(file_molecules)
                close(file_molecules)
                #-------------------
                file_dimers=open(string(folder_local,"dimer-",cut_off,".dat"))
                lines_dimers=readlines(file_dimers)
                close(file_dimers)
                #--------------------
                file_dimers_alone=open(string(folder_local,"dimer_alone-",cut_off,".dat"),"w")
                file_dimers_chain=open(string(folder_local,"dimer_chain-",cut_off,".dat"),"w")
                alone=0
                step_init=1
                for dimer=1:size(lines_dimers)[1]
                    print("V=",V," T=",T,"K Progress: ", dimer/size(lines_dimers)[1]*100,"% \n")
                    line=split(lines_dimers[dimer])
                    step=parse(Int64,line[1])
                    C1=parse(Int64,line[2])
                    C2=parse(Int64,line[3])
                    for molecule=step_init:size(lines_molecules)[1]
                        line2=split(lines_molecules[molecule])
                        if step != parse(Int64, line2[1])
                            continue
                        end
                        # Check only with carbon because too hard otherwise
                        if C1 == parse(Int64, line2[3+1]) && C2 == parse(Int64, line2[3+2] )
                            if parse(Int64,line2[size(line2)[1]]) == 6
                                for k=1:size(line)[1]
                                    write(file_dimers_alone,string(line[k]," "))
                                end
                                write(file_dimers_alone,string("\n"))
                                alone +=1
                            else
                                for k=1:size(line)[1]
                                    write(file_dimers_chain,string(line[k]," "))
                                end
                                write(file_dimers_chain,string("\n"))
                            end
                            # No point looking backward
                            step_init=molecule+1
                        end
                    end
                end
                file_prop=open(string(folder_local,"dimer_in-out_prop.dat"),"w")
                write(file_prop,string(alone/size(lines_dimers)[1]," ",1-alone/size(lines_dimers)[1],"\n"))
                close(file_prop)
                close(file_dimers_alone)
                close(file_dimers_chain)
            end
        end
    end
end
