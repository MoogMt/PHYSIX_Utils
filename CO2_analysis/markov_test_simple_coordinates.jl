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

V=8.82
T=3000
cut_off=1.75

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")

traj = filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

nb_atoms=size(traj[1].names)[1]
nb_steps=size(traj)[1]

restart=true

coord_matrix=zeros(nb_steps,nbC,2)

file_out=open(string(folder_out,"coordinance-C-",cut_off,".dat"),"w")
for step_sim=1:nb_steps

    print("Progress: ",step_sim/nb_steps*100,"%\n")
    write(file_out,string(step_sim," "))
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
        coord_matrix[step_sim,carbon,1]=sum(bond_matrix[carbon,1:nbC])
        coord_matrix[step_sim,carbon,2]=sum(bond_matrix[carbon,nbC+1:nbC+nbO])
        write(file_out,string(sum(bond_matrix[carbon,1:nbC])," ",sum(bond_matrix[carbon,nbC+1:nbC+nbO])," ") )
    end
    write(file_out,string("\n"))
end
close(file_out)

guess_cases=false
case_matrix=zeros(Int,nb_steps,nbC)

global cases=ones(0,1)

global count_cases=[1]

if guess_cases
    for step_sim=1:nb_steps
        print("Computing cases - Progress: ",step_sim/nb_steps*100,"%\n")
        global cases
        global count_cases
        for carbon=1:nbC
            global cases
            global count_cases
            if size(cases)[1] == 0
                global cases=ones(1,2)
                cases[1,1]= coord_matrix[step_sim,carbon,1]
                cases[1,2]= coord_matrix[step_sim,carbon,2]
            else
                global new=true
                if size(cases)[2] == 1
                    if cases[1] == coord_matrix[step_sim,carbon,1] && cases[2] == coord_matrix[step_sim,carbon,2]
                        global count[i] += 1
                        global new=false
                    end
                else
                    for i=1:size(cases)[1]
                        if cases[i,1] == coord_matrix[step_sim,carbon,1] && cases[i,2] == coord_matrix[step_sim,carbon,2]
                            global count_cases[i] += 1
                            global new=false
                        end
                    end
                end
                if new
                    global cases=[cases; transpose(coord_matrix[step_sim,carbon,:]) ]
                    push!(count_cases,1)
                end
            end
        end
    end
    file_out=open(string(folder_out,"cases_simple-",cut_off,".dat"),"w")
    for step_sim=1:nb_steps
        print("Case affectation - Progress:",step_sim/nb_steps*100,"%\n")
        for carbon=1:nbC
            for i=1:size(cases)[1]
                if cases[i,1] == coord_matrix[step_sim,carbon,1] && cases[i,2] == coord_matrix[step_sim,carbon,2]
                    case_matrix[step_sim,carbon]=i
                end
            end
            write(file_out,string(case_matrix[step_sim,carbon]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
else
    cases_CO2=[0,2]
    cases_CO3=[0,3]
    cases_CO4=[0,4]
    cases_C1O2=[1,2]
    cases_C1O3=[1,3]
    cases=zeros(5,2)
    cases[1,:]=cases_CO2
    cases[2,:]=cases_CO3
    cases[3,:]=cases_CO4
    cases[4,:]=cases_C1O2
    cases[5,:]=cases_C1O3
    count_cases=zeros(5)

    file_out=open(string(folder_out,"cases_simple-",cut_off,"-imposed.dat"),"w")
    for step_sim=1:nb_steps
        print("Computing cases - Progress: ",step_sim/nb_steps*100,"%\n")
        global count_cases
        for carbon=1:nbC
            global index=1
            i=1
            global d=0
            for k=1:2
                global d+= (coord_matrix[step_sim,carbon,k] - cases[i,k])*(coord_matrix[step_sim,carbon,k] - cases[i,k])
            end
            global d_min=d
            for i=2:size(cases)[1]
                global d=0
                for k=1:2
                    global d+= (coord_matrix[step_sim,carbon,k] - cases[i,k])*(coord_matrix[step_sim,carbon,k] - cases[i,k])
                end
                if d < d_min
                    global index=i
                    global d_min=d
                end
            end
            count_cases[index] += 1
            case_matrix[step_sim,carbon]=index
        end
    end
    close(file_out)
end


lag_min=1
d_lag=10
lag_max=5001
case_transition=zeros(Float64,size(cases)[1],size(cases)[1],Int(trunc((lag_max-lag_min)/d_lag)))
global count_lag=1
for lag=lag_min:d_lag:lag_max-1
    print("Lag Compute - Progress: ",lag/lag_max*100,"%\n")
    for carbon=1:nbC
        for step_sim=lag+1:nb_steps
            case_transition[ case_matrix[step_sim-lag,carbon], case_matrix[step_sim,carbon], count_lag ] += 1
        end
    end
    global count_lag += 1
end

global max_case=5
if ! guess_cases
    global max_cases=size(count_cases)[1]
end

case_proba=case_transition[1:max_case,1:max_case,:]
for lag=1:size(case_proba)[3]
    for i=1:max_case
        case_proba[:,i,lag] /= sum( case_proba[:,i,lag] )
    end
end

file_lag=open(string(folder_out,"markovian_evolution_test.dat"),"w")
for lag=1:count_lag-1
    write(file_lag,string(lag*d_lag*0.005," "))
    for i=1:max_case
        for j=1:max_case
            # Specific i-j transition
            global p_truth=case_proba[i,j,lag]
            global p_test=0
            global p_test2=0
            for k=1:max_case
                global p_test += case_proba[i,k,Int(trunc(lag/2))+1]*case_proba[k,j,Int(trunc(lag/2))+1]
            end
            write(file_lag,string(p_truth," ",p_test," "))
        end
    end
    write(file_lag,string("\n"))
end
close(file_lag)

file_lag=open(string(folder_out,"markovian_evolution_test-2.dat"),"w")
for lag=2:count_lag-1
    case_proba_test=case_proba[:,:,1]^lag
    write(file_lag,string(lag*d_lag*0.005," "))
    for i=1:max_case
        for j=1:max_case
            # Specific i-j transition
            write(file_lag,string(case_proba[i,j,lag]," ",case_proba_test[i,j]," "))
        end
    end
    write(file_lag,string("\n"))
end
