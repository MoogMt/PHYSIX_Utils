GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2
max_coord=5

restart=false

for V in Volumes
    for T in Temperatures
        for cut_off in Cut_Off
            #for cut_off in Cut_Off

            folder_in=string(folder_base,V,"/",T,"K/")
            file=string(folder_in,"TRAJEC_wrapped.xyz")

            if ! isfile(file)
                continue
            end

            folder_out=string(folder_in,"Data/")

            traj = filexyz.readFastFile(file)
            cell=cell_mod.Cell_param(V,V,V)

            nb_atoms=size(traj[1].names)[1]
            nb_steps=size(traj)[1]

            restart=true

            coord_matrix=ones(nb_steps,nbC,8)*(-1)

            for step_sim=1:nb_steps
                print("Progress: ",step_sim/nb_steps*100,"%\n")
                bond_matrix=zeros(nb_atoms,nb_atoms)
                for atom1=1:nb_atoms
                    for atom2=atom1+1:nb_atoms
                        if cell_mod.distance(traj[step_sim],cell,atom1,atom2) < cut_off
                            bond_matrix[atom1,atom2]=1
                            bond_matrix[atom2,atom1]=1
                        end
                    end
                end

                for carbon=1:nbC
                    count_coord=1
                    for carbon2=1:nbC
                        if carbon == carbon2
                            continue
                        end
                        if bond_matrix[carbon,carbon2] > 0
                            coord_matrix[step_sim,carbon,count_coord]=sum(bond_matrix[carbon2,:])
                            count_coord += 1
                        end
                    end
                    count_coord=1
                    for oxygen=1:nbO
                        if bond_matrix[carbon,nbC+oxygen] > 0
                            if count_coord > 4
                                continue
                            end
                            coord_matrix[step_sim,carbon,4+count_coord]=sum(bond_matrix[nbC+oxygen,:])
                            count_coord += 1
                        end
                    end
                    # sort
                    for i=1:4
                        for j=i+1:4
                            if coord_matrix[step_sim,carbon,i] < coord_matrix[step_sim,carbon,j]
                                stock=coord_matrix[step_sim,carbon,i]
                                coord_matrix[step_sim,carbon,i]=coord_matrix[step_sim,carbon,j]
                                coord_matrix[step_sim,carbon,j]=stock
                            end
                        end
                    end
                    for i=5
                        for j=i+1:8
                            if coord_matrix[step_sim,carbon,i] < coord_matrix[step_sim,carbon,j]
                                stock=coord_matrix[step_sim,carbon,i]
                                coord_matrix[step_sim,carbon,i]=coord_matrix[step_sim,carbon,j]
                                coord_matrix[step_sim,carbon,j]=stock
                            end
                        end
                    end
                end
            end

            guess_cases=false
            case_matrix=zeros(Int,nb_steps,nbC)

            cases=zeros(39,8)
            cases[1,:]=[-1,-1,-1,-1,1,1,-1,-1]
            cases[2,:]=[-1,-1,-1,-1,2,1,-1,-1]
            cases[3,:]=[-1,-1,-1,-1,2,2,-1,-1]
            cases[4,:]=[-1,-1,-1,-1,2,2,2,-1]
            cases[5,:]=[-1,-1,-1,-1,2,2,1,-1]
            cases[6,:]=[-1,-1,-1,-1,2,1,1,-1]
            cases[7,:]=[-1,-1,-1,-1,1,1,1,-1]
            cases[8,:]=[-1,-1,-1,-1,2,2,2,2]
            cases[9,:]=[-1,-1,-1,-1,2,2,2,1]
            cases[10,:]=[-1,-1,-1,-1,2,2,1,1]
            cases[11,:]=[-1,-1,-1,-1,2,1,1,1]
            cases[12,:]=[-1,-1,-1,-1,1,1,1,1]
            cases[13,:]=[2,-1,-1,-1,2,-1,-1,-1]
            cases[14,:]=[2,-1,-1,-1,1,-1,-1,-1]
            cases[15,:]=[3,-1,-1,-1,2,-1,-1,-1]
            cases[16,:]=[3,-1,-1,-1,1,-1,-1,-1]
            cases[17,:]=[4,-1,-1,-1,2,-1,-1,-1]
            cases[18,:]=[4,-1,-1,-1,1,-1,-1,-1]
            cases[19,:]=[2,-1,-1,-1,2,2,-1,-1]
            cases[20,:]=[2,-1,-1,-1,2,1,-1,-1]
            cases[21,:]=[2,-1,-1,-1,1,1,-1,-1]
            cases[22,:]=[3,-1,-1,-1,2,2,-1,-1]
            cases[23,:]=[3,-1,-1,-1,2,1,-1,-1]
            cases[24,:]=[3,-1,-1,-1,1,1,-1,-1]
            cases[25,:]=[4,-1,-1,-1,2,2,-1,-1]
            cases[26,:]=[4,-1,-1,-1,2,1,-1,-1]
            cases[27,:]=[4,-1,-1,-1,1,1,-1,-1]
            cases[28,:]=[2,-1,-1,-1,2,2,2,-1]
            cases[29,:]=[2,-1,-1,-1,2,2,1,-1]
            cases[30,:]=[2,-1,-1,-1,2,1,1,-1]
            cases[31,:]=[2,-1,-1,-1,1,1,1,-1]
            cases[32,:]=[3,-1,-1,-1,2,2,2,-1]
            cases[33,:]=[3,-1,-1,-1,2,2,1,-1]
            cases[34,:]=[3,-1,-1,-1,2,1,1,-1]
            cases[35,:]=[3,-1,-1,-1,1,1,1,-1]
            cases[36,:]=[4,-1,-1,-1,2,2,2,-1]
            cases[37,:]=[4,-1,-1,-1,2,2,1,-1]
            cases[38,:]=[4,-1,-1,-1,2,1,1,-1]
            cases[39,:]=[4,-1,-1,-1,1,1,1,-1]

            count_cases=zeros(39)

            for step_sim=1:nb_steps
                print("Computing cases - Progress: ",step_sim/nb_steps*100,"%\n")
                count_cases
                for carbon=1:nbC
                    index=1
                    i=1
                    d=0
                    for k=1:size(cases)[2]
                        d+= (coord_matrix[step_sim,carbon,k] - cases[i,k])*(coord_matrix[step_sim,carbon,k] - cases[i,k])
                    end
                    d_min=d
                    for i=2:size(cases)[1]
                        d=0
                        for k=1:size(cases)[2]
                            d+= (coord_matrix[step_sim,carbon,k] - cases[i,k])*(coord_matrix[step_sim,carbon,k] - cases[i,k])
                        end
                        if d < d_min
                            index=i
                            d_min=d
                        end
                    end
                    count_cases[index] += 1
                    case_matrix[step_sim,carbon]=index
                end
            end

            percent=count_cases/sum(count_cases)*100

            new_size=0
            index_keep=[]
            cases_keep=zeros(0,8)
            for i=1:size(count_cases)[1]
                if percent[i] > 0.5
                    new_size += 1
                    push!(index_keep, i)
                    cases_keep=[ cases_keep ; transpose(cases[i,:]) ]
                end
            end

            count_cases_keep=zeros(size(cases_keep)[1])
            case_matrix_keep=zeros(Int,nb_steps,nbC)
            for step_sim=1:nb_steps
                print("Computing cases - 2 - Progress: ",step_sim/nb_steps*100,"%\n")
                count_cases_keep
                for carbon=1:nbC
                    index=1
                    i=1
                    d=0
                    for k=1:size(cases_keep)[2]
                        d+= (coord_matrix[step_sim,carbon,k] - cases_keep[i,k])*(coord_matrix[step_sim,carbon,k] - cases_keep[i,k])
                    end
                    d_min=d
                    for i=2:size(cases_keep)[1]
                        d=0
                        for k=1:size(cases_keep)[2]
                            d+= (coord_matrix[step_sim,carbon,k] - cases_keep[i,k])*(coord_matrix[step_sim,carbon,k] - cases_keep[i,k])
                        end
                        if d < d_min
                            index=i
                            d_min=d
                        end
                    end
                    count_cases_keep[index] += 1
                    case_matrix_keep[step_sim,carbon]=index
                end
            end

            file_cases_keep=open(string(folder_out,"cases-",cut_off,".dat"),"w")
            if size(count_cases_keep)[1] == 1
                for i=1:8
                    write(file_cases_keep,string(cases_keep[i]," "))
                end
                write(file_cases_keep,string(count_cases_keep,"\n"))
            else
                for i=size(count_cases_keep)[1]
                    for j=1:8
                        write(file_cases_keep,string(cases_keep[i,j]," "))
                    end
                    write(file_cases_keep,string(count_cases_keep,"\n"))
                end
            end
            close(file_cases_keep)

            file_cases=open(string(folder_out,"cases_occurences-",cut_off,".dat"),"w")
            for step_sim=1:nb_steps
                for carbon=1:nbC
                    write(file_cases,string(case_matrix_keep[step_sim,carbon]," "))
                end
                write(file_cases,string("\n"))
            end
            close(file_cases)

            lag_min=1
            d_lag=10
            lag_max=5001
            case_transition=zeros(Float64,size(cases_keep)[1],size(cases_keep)[1],Int(trunc((lag_max-lag_min)/d_lag)))
            count_lag=1
            for lag=lag_min:d_lag:lag_max-1
                print("Lag Compute - Progress: ",lag/lag_max*100,"%\n")
                for carbon=1:nbC
                    for step_sim=lag+1:nb_steps
                        case_transition[ case_matrix_keep[step_sim-lag,carbon], case_matrix_keep[step_sim,carbon], count_lag ] += 1
                    end
                end
                count_lag += 1
            end

            case_proba=case_transition[:,:,:]
            for lag=1:size(case_proba)[3]
                for i=1:size(cases_keep)[1]
                    case_proba[:,i,lag] /= sum( case_proba[:,i,lag] )
                end
            end
            file_proba=open(string(folder_out,"proba_transition-",cut_off,".dat"),"w")
            for lag=1:size(case_proba)[3]
                if size(count_cases_keep)[1] == 1
                    write(file_proba,string(case_proba[1,1,lag]," "))
                    continue
                end
                for i=1:size(case_proba)[1]
                    for j=1:size(case_proba)[2]
                        write(file_proba,string(case_proba[i,j,lag]," "))
                    end
                end
                write(file_proba,string("\n"))
            end
            close(file_proba)


            file_lag=open(string(folder_out,"markovian_evolution_test-",cut_off,".dat"),"w")
            file_out_2=open(string(folder_out,"markovian_evolution_test2-",cut_off,".dat"),"w")
            for lag=1:count_lag-1
                d_test=0
                d_test2=0
                count2=0
                for i=1:size(cases_keep)[1]
                    for j=1:size(cases_keep)[1]
                        # Specific i-j transition
                        p_truth=case_proba[i,j,lag]
                        p_test=0
                        p_test2=0
                        for k=1:size(cases_keep)[1]
                            p_test += case_proba[i,k,Int(trunc(lag/2))+1]*case_proba[k,j,Int(trunc(lag/2))+1]
                        end
                        d_test += abs(p_truth-p_test)
                        d_test2 += d_test*d_test
                        count2 += 1
                        write(file_lag,string(p_truth," ",p_test," "))
                    end
                end
                write(file_lag,string("\n"))
                write(file_out_2,string(lag*d_lag*0.005," ",d_test/count2," ",d_test2/count2-(d_test/count2)^2,"\n"))
            end
            close(file_lag)
            close(file_out_2)

            file_lag=open(string(folder_out,"markovian_evolution_test-3_",cut_off,".dat"),"w")
            for lag=2:count_lag-1
                case_proba_test=case_proba[:,:,1]^lag
                write(file_lag,string(lag*d_lag*0.005," "))
                for i=1:size(cases_keep)[1]
                    for j=1:size(cases_keep)[1]
                        # Specific i-j transition
                        write(file_lag,string(case_proba[i,j,lag]," ",case_proba_test[i,j]," "))
                    end
                end
                write(file_lag,string("\n"))
            end

            #
        end
    end
end
