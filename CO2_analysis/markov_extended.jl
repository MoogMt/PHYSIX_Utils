GPfolder=string("/home/moogmt/PHYSIX_Utils/GPlib/Julia/")

include(string(GPfolder,"contactmatrix.jl"))
include(string(GPfolder,"geom.jl"))

function buildCoordinationMatrix( traj::Vector{T1}, cell::T2, cut_off_bond::T3 ) where { T1 <: atom_mod.AtomList, T2 <: cell_mod.Cell_param , T3 <: Real }
    nb_atoms=size(traj[1].names)[1]
    nb_steps=size(traj)[1]
    n_dim=9
    coord_matrix=ones(nb_steps,nbC,n_dim)*(-1)
    for step_sim=1:nb_steps

        print("Building Coordination Signal - Progress: ",step_sim/nb_steps*100,"%\n")

        # Bond Matrix
        bond_matrix=zeros(nb_atoms,nb_atoms)
        for carbon=1:32
            for atom=1:96
                if carbon == atom
                    continue
                end
                if cell_mod.distance( traj[step_sim], cell, carbon, atom) < cut_off_bond
                    bond_matrix[carbon,atom]=1
                    bond_matrix[atom,carbon]=1
                end
            end
        end
        for oxygen1=1:64-1
            for oxygen2=oxygen1+1:64
                if cell_mod.distance( traj[step_sim], cell, nbC+oxygen1, nbC+oxygen2 ) < cut_off_bond
                    bond_matrix[nbC+oxygen1,nbC+oxygen2] = 1
                    bond_matrix[nbC+oxygen2,nbC+oxygen1] = 1
                end
            end
        end

        # Cleaning weird 3-member rings
        for atom1=1:nb_atoms-1
            for atom2=atom1+1:nb_atoms
                if bond_matrix[atom1,atom2] == 1
                    # For each bonded pair of atoms, check whether they are
                    # both bonded to another atom, forming unatural 3-member
                    # ring
                    for atom3=1:nb_atoms
                        if atom3 == atom1 || atom3 == atom1
                            continue
                        end
                        if bond_matrix[atom3,atom1] == 1 && bond_matrix[atom3,atom2] == 1
                            if cell_mod.distance(traj[step_sim],cell,atom1,atom3) > cell_mod.distance(traj[step_sim],cell,atom1,atom2) && cell_mod.distance(traj[step_sim],cell,atom1,atom3) > cell_mod.distance(traj[step_sim],cell,atom2,atom3)
                                bond_matrix[atom1,atom3]=0
                                bond_matrix[atom3,atom1]=0
                            elseif cell_mod.distance(traj[step_sim],cell,atom1,atom2) > cell_mod.distance(traj[step_sim],cell,atom2,atom3)
                                bond_matrix[atom1,atom2]=0
                                bond_matrix[atom2,atom1]=0
                            else
                                bond_matrix[atom2,atom3]=0
                                bond_matrix[atom3,atom2]=0
                            end
                        end
                    end
                end
            end
        end


        # Compute coord matrix
        for carbon=1:nbC
            # compute coord
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
                    coord_matrix[step_sim,carbon,4+count_coord]=sum(bond_matrix[nbC+oxygen,:])
                    count_coord += 1
                end
            end
            # sort coord
            for i=1:4
                for j=i+1:4
                    if coord_matrix[step_sim,carbon,i] < coord_matrix[step_sim,carbon,j]
                        stock=coord_matrix[step_sim,carbon,i]
                        coord_matrix[step_sim,carbon,i]=coord_matrix[step_sim,carbon,j]
                        coord_matrix[step_sim,carbon,j]=stock
                    end
                end
            end
            for i=5:n_dim
                for j=i+1:n_dim
                    if coord_matrix[step_sim,carbon,i] < coord_matrix[step_sim,carbon,j]
                        stock=coord_matrix[step_sim,carbon,i]
                        coord_matrix[step_sim,carbon,i]=coord_matrix[step_sim,carbon,j]
                        coord_matrix[step_sim,carbon,j]=stock
                    end
                end
            end
        end
    end
    return coord_matrix
end
function defineStatesCoordinances()
    states=zeros(39,2)
    # CO2
    states[1,:]=[0,1]
    states[2,:]=[0,2]
    states[3,:]=[0,3]
    states[4,:]=[0,4]
    states[5,:]=[1,1]
    states[6,:]=[1,2]
    states[7,:]=[1,3]
    states[8,:]=[1,4]
    states[9,:]=[2,1]
    states[10,:]=[2,2]
    states[11,:]=[2,3]
    return states
end
function defineStatesExtendedCoordinances()
    states=zeros(86,9)
    states[1,:]=[-1,-1,-1,-1,1,-1,-1,-1,-1]
    states[2,:]=[-1,-1,-1,-1,2,-1,-1,-1,-1]
    states[3,:]=[-1,-1,-1,-1,3,-1,-1,-1,-1]
    # CO2
    states[4,:]=[-1,-1,-1,-1,1,1,-1,-1,-1]
    states[5,:]=[-1,-1,-1,-1,2,1,-1,-1,-1]
    states[6,:]=[-1,-1,-1,-1,2,2,-1,-1,-1]
    states[7,:]=[-1,-1,-1,-1,3,1,-1,-1,-1]
    states[8,:]=[-1,-1,-1,-1,3,2,-1,-1,-1]
    states[9,:]=[-1,-1,-1,-1,3,3,-1,-1,-1]
    # CO3
    states[10,:]=[-1,-1,-1,-1,2,2,2,-1,-1]
    states[11,:]=[-1,-1,-1,-1,2,2,1,-1,-1]
    states[12,:]=[-1,-1,-1,-1,2,1,1,-1,-1]
    states[13,:]=[-1,-1,-1,-1,1,1,1,-1,-1]
    #
    states[14,:]=[-1,-1,-1,-1,3,1,1,-1,-1]
    states[15,:]=[-1,-1,-1,-1,3,2,1,-1,-1]
    states[16,:]=[-1,-1,-1,-1,3,2,2,-1,-1]
    states[17,:]=[-1,-1,-1,-1,3,3,1,-1,-1]
    states[18,:]=[-1,-1,-1,-1,3,3,2,-1,-1]
    states[19,:]=[-1,-1,-1,-1,3,3,3,-1,-1]
    # CO4
    states[20,:]=[-1,-1,-1,-1,2,2,2,2,-1]
    states[21,:]=[-1,-1,-1,-1,2,2,2,1,-1]
    states[22,:]=[-1,-1,-1,-1,2,2,1,1,-1]
    states[23,:]=[-1,-1,-1,-1,2,1,1,1,-1]
    states[24,:]=[-1,-1,-1,-1,1,1,1,1,-1]
    states[25,:]=[-1,-1,-1,-1,3,1,1,1,-1]
    states[26,:]=[-1,-1,-1,-1,3,2,1,1,-1]
    states[27,:]=[-1,-1,-1,-1,3,2,2,1,-1]
    states[28,:]=[-1,-1,-1,-1,3,2,2,2,-1]
    states[29,:]=[-1,-1,-1,-1,3,3,1,1,-1]
    states[30,:]=[-1,-1,-1,-1,3,3,2,1,-1]
    states[31,:]=[-1,-1,-1,-1,3,3,2,2,-1]
    states[32,:]=[-1,-1,-1,-1,3,3,3,1,-1]
    states[33,:]=[-1,-1,-1,-1,3,3,3,2,-1]
    states[34,:]=[-1,-1,-1,-1,3,3,3,3,-1]
    # CO5
    states[35,:]=[-1,-1,-1,-1,2,2,2,2,2]
    states[36,:]=[-1,-1,-1,-1,2,2,2,2,1]
    states[37,:]=[-1,-1,-1,-1,2,2,2,1,1]
    states[38,:]=[-1,-1,-1,-1,2,2,1,1,1]
    states[39,:]=[-1,-1,-1,-1,2,1,1,1,1]
    # CCO1
    states[40,:]=[2,-1,-1,-1,2,-1,-1,-1,-1]
    states[41,:]=[2,-1,-1,-1,1,-1,-1,-1,-1]
    states[42,:]=[3,-1,-1,-1,2,-1,-1,-1,-1]
    states[43,:]=[3,-1,-1,-1,1,-1,-1,-1,-1]
    states[44,:]=[4,-1,-1,-1,2,-1,-1,-1,-1]
    states[45,:]=[4,-1,-1,-1,1,-1,-1,-1,-1]
    # CCO2
    states[46,:]=[2,-1,-1,-1,2,2,-1,-1,-1]
    states[47,:]=[2,-1,-1,-1,2,1,-1,-1,-1]
    states[48,:]=[2,-1,-1,-1,1,1,-1,-1,-1]
    states[49,:]=[3,-1,-1,-1,2,2,-1,-1,-1]
    states[50,:]=[3,-1,-1,-1,2,1,-1,-1,-1]
    states[51,:]=[3,-1,-1,-1,1,1,-1,-1,-1]
    states[52,:]=[4,-1,-1,-1,2,2,-1,-1,-1]
    states[53,:]=[4,-1,-1,-1,2,1,-1,-1,-1]
    states[54,:]=[4,-1,-1,-1,1,1,-1,-1,-1]
    # CCO3
    states[55,:]=[2,-1,-1,-1,2,2,2,-1,-1]
    states[56,:]=[2,-1,-1,-1,2,2,1,-1,-1]
    states[57,:]=[2,-1,-1,-1,2,1,1,-1,-1]
    states[58,:]=[2,-1,-1,-1,1,1,1,-1,-1]
    states[59,:]=[3,-1,-1,-1,2,2,2,-1,-1]
    states[60,:]=[3,-1,-1,-1,2,2,1,-1,-1]
    states[61,:]=[3,-1,-1,-1,2,1,1,-1,-1]
    states[62,:]=[3,-1,-1,-1,1,1,1,-1,-1]
    states[63,:]=[4,-1,-1,-1,2,2,2,-1,-1]
    states[64,:]=[4,-1,-1,-1,2,2,1,-1,-1]
    states[65,:]=[4,-1,-1,-1,2,1,1,-1,-1]
    states[66,:]=[4,-1,-1,-1,1,1,1,-1,-1]
    #
    states[67,:]=[2,-1,-1,-1,2,2,2,1,-1]
    states[68,:]=[2,-1,-1,-1,2,2,1,1,-1]
    states[69,:]=[2,-1,-1,-1,2,1,1,1,-1]
    states[70,:]=[2,-1,-1,-1,1,1,1,1,-1]
    states[71,:]=[3,-1,-1,-1,2,2,2,1,-1]
    states[72,:]=[3,-1,-1,-1,2,2,1,1,-1]
    states[73,:]=[3,-1,-1,-1,2,1,1,1,-1]
    states[74,:]=[3,-1,-1,-1,1,1,1,1,-1]
    states[75,:]=[4,-1,-1,-1,2,2,2,1,-1]
    states[76,:]=[4,-1,-1,-1,2,2,1,1,-1]
    states[77,:]=[4,-1,-1,-1,2,1,1,1,-1]
    states[78,:]=[4,-1,-1,-1,1,1,1,1,-1]
    #
    states[79,:]=[2,-1,-1,-1,2,2,2,2,-1]
    states[80,:]=[2,-1,-1,-1,2,2,2,1,-1]
    states[81,:]=[2,-1,-1,-1,2,2,1,1,-1]
    states[82,:]=[2,-1,-1,-1,2,1,1,1,-1]
    states[83,:]=[3,-1,-1,-1,2,2,2,2,-1]
    states[84,:]=[3,-1,-1,-1,2,2,2,1,-1]
    states[85,:]=[3,-1,-1,-1,2,2,1,1,-1]
    states[86,:]=[3,-1,-1,-1,2,1,1,1,-1]
    return states
end
function defineStatesExtendedCoordinancesO()
    states=zeros(351,6)
    # CO2
    states[1,:]=[1,-1,-1,-1,-1,-1]
    states[2,:]=[2,-1,-1,-1,-1,-1]
    states[3,:]=[3,-1,-1,-1,-1,-1]
    states[4,:]=[4,-1,-1,-1,-1,-1]
    #
    states[5,:]=[1,1,-1,-1,-1,-1]
    states[6,:]=[2,1,-1,-1,-1,-1]
    states[7,:]=[2,2,-1,-1,-1,-1]
    states[8,:]=[2,3,-1,-1,-1,-1]
    states[9,:]=[2,4,-1,-1,-1,-1]
    states[10,:]=[3,1,-1,-1,-1,-1]
    states[11,:]=[3,2,-1,-1,-1,-1]
    states[12,:]=[3,3,-1,-1,-1,-1]
    states[13,:]=[3,4,-1,-1,-1,-1]
    states[14,:]=[4,1,-1,-1,-1,-1]
    states[15,:]=[4,2,-1,-1,-1,-1]
    states[16,:]=[4,3,-1,-1,-1,-1]
    states[17,:]=[4,4,-1,-1,-1,-1]
    #
    states[18,:]=[1,1,1,-1,-1,-1]
    states[19,:]=[2,2,1,-1,-1,-1]
    states[20,:]=[2,2,2,-1,-1,-1]
    states[21,:]=[2,2,3,-1,-1,-1]
    states[22,:]=[2,3,1,-1,-1,-1]
    states[23,:]=[2,3,2,-1,-1,-1]
    states[24,:]=[2,3,3,-1,-1,-1]
    states[25,:]=[2,4,1,-1,-1,-1]
    states[26,:]=[2,4,2,-1,-1,-1]
    states[27,:]=[2,4,3,-1,-1,-1]
    states[28,:]=[3,2,1,-1,-1,-1]
    states[29,:]=[3,2,2,-1,-1,-1]
    states[30,:]=[3,3,2,-1,-1,-1]
    states[31,:]=[3,3,3,-1,-1,-1]
    states[32,:]=[4,1,1,-1,-1,-1]
    states[33,:]=[4,2,1,-1,-1,-1]
    states[34,:]=[4,2,2,-1,-1,-1]
    states[35,:]=[4,3,2,-1,-1,-1]
    states[36,:]=[4,4,1,-1,-1,-1]
    states[37,:]=[4,4,2,-1,-1,-1]
    states[38,:]=[4,4,3,-1,-1,-1]
    states[39,:]=[4,4,4,-1,-1,-1]
    #
    states[40,:]=[1,1,1,1,-1,-1]
    states[41,:]=[2,1,1,1,-1,-1]
    states[42,:]=[2,2,1,1,-1,-1]
    states[43,:]=[2,2,2,1,-1,-1]
    states[44,:]=[2,2,2,2,-1,-1]
    states[45,:]=[3,1,1,1,-1,-1]
    states[46,:]=[3,2,1,1,-1,-1]
    states[48,:]=[3,2,2,1,-1,-1]
    states[49,:]=[3,2,2,2,-1,-1]
    states[50,:]=[3,3,1,1,-1,-1]
    states[51,:]=[3,3,2,1,-1,-1]
    states[52,:]=[3,3,2,2,-1,-1]
    states[53,:]=[3,3,3,1,-1,-1]
    states[54,:]=[3,3,3,2,-1,-1]
    states[55,:]=[3,3,3,3,-1,-1]
    states[56,:]=[4,1,1,1,-1,-1]
    states[57,:]=[4,2,1,1,-1,-1]
    states[58,:]=[4,2,2,1,-1,-1]
    states[59,:]=[4,2,2,2,-1,-1]
    states[60,:]=[4,3,1,1,-1,-1]
    states[61,:]=[4,3,2,1,-1,-1]
    states[62,:]=[4,3,2,2,-1,-1]
    states[63,:]=[4,3,3,1,-1,-1]
    states[64,:]=[4,3,3,2,-1,-1]
    states[65,:]=[4,3,3,3,-1,-1]
    states[66,:]=[4,4,1,1,-1,-1]
    states[67,:]=[4,4,2,1,-1,-1]
    states[68,:]=[4,4,2,2,-1,-1]
    states[69,:]=[4,4,3,2,-1,-1]
    states[70,:]=[4,4,3,3,-1,-1]
    states[71,:]=[4,4,4,4,-1,-1]
    ###
    states[72,:]=[1,-1,-1,-1,1,-1]
    states[73,:]=[2,-1,-1,-1,1,-1]
    states[74,:]=[3,-1,-1,-1,1,-1]
    states[75,:]=[4,-1,-1,-1,1,-1]
    #
    states[76,:]=[1,1,-1,-1,1,-1]
    states[77,:]=[2,1,-1,-1,1,-1]
    states[78,:]=[2,2,-1,-1,1,-1]
    states[79,:]=[2,3,-1,-1,1,-1]
    states[80,:]=[2,4,-1,-1,1,-1]
    states[81,:]=[3,1,-1,-1,1,-1]
    states[82,:]=[3,2,-1,-1,1,-1]
    states[83,:]=[3,3,-1,-1,1,-1]
    states[84,:]=[3,4,-1,-1,1,-1]
    states[85,:]=[4,1,-1,-1,1,-1]
    states[86,:]=[4,2,-1,-1,1,-1]
    states[87,:]=[4,3,-1,-1,1,-1]
    states[88,:]=[4,4,-1,-1,1,-1]
    #
    states[89,:]=[1,1,1,-1,1,-1]
    states[90,:]=[2,2,1,-1,1,-1]
    states[91,:]=[2,2,2,-1,1,-1]
    states[92,:]=[2,2,3,-1,1,-1]
    states[93,:]=[2,3,1,-1,1,-1]
    states[94,:]=[2,3,2,-1,1,-1]
    states[95,:]=[2,3,3,-1,1,-1]
    states[96,:]=[2,4,1,-1,1,-1]
    states[97,:]=[2,4,2,-1,1,-1]
    states[98,:]=[2,4,3,-1,1,-1]
    states[99,:]=[3,2,1,-1,1,-1]
    states[100,:]=[3,2,2,-1,1,-1]
    states[101,:]=[3,3,2,-1,1,-1]
    states[102,:]=[3,3,3,-1,1,-1]
    states[103,:]=[4,1,1,-1,1,-1]
    states[104,:]=[4,2,1,-1,1,-1]
    states[105,:]=[4,2,2,-1,1,-1]
    states[106,:]=[4,3,2,-1,1,-1]
    states[107,:]=[4,4,1,-1,1,-1]
    states[108,:]=[4,4,2,-1,1,-1]
    states[109,:]=[4,4,3,-1,1,-1]
    states[110,:]=[4,4,4,-1,1,-1]
    #
    states[111,:]=[1,1,1,1,1,-1]
    states[112,:]=[2,1,1,1,1,-1]
    states[113,:]=[2,2,1,1,1,-1]
    states[114,:]=[2,2,2,1,1,-1]
    states[115,:]=[2,2,2,2,1,-1]
    states[116,:]=[3,1,1,1,1,-1]
    states[117,:]=[3,2,1,1,1,-1]
    states[118,:]=[3,2,2,1,1,-1]
    states[119,:]=[3,2,2,2,1,-1]
    states[120,:]=[3,3,1,1,1,-1]
    states[121,:]=[3,3,2,1,1,-1]
    states[122,:]=[3,3,2,2,1,-1]
    states[123,:]=[3,3,3,1,1,-1]
    states[124,:]=[3,3,3,2,1,-1]
    states[125,:]=[3,3,3,3,1,-1]
    states[126,:]=[4,1,1,1,1,-1]
    states[217,:]=[4,2,1,1,1,-1]
    states[128,:]=[4,2,2,1,1,-1]
    states[129,:]=[4,2,2,2,1,-1]
    states[130,:]=[4,3,1,1,1,-1]
    states[131,:]=[4,3,2,1,1,-1]
    states[132,:]=[4,3,2,2,1,-1]
    states[133,:]=[4,3,3,1,1,-1]
    states[134,:]=[4,3,3,2,1,-1]
    states[135,:]=[4,3,3,3,1,-1]
    states[136,:]=[4,4,1,1,1,-1]
    states[137,:]=[4,4,2,1,1,-1]
    states[138,:]=[4,4,2,2,1,-1]
    states[139,:]=[4,4,3,2,1,-1]
    states[140,:]=[4,4,3,3,1,-1]
    states[141,:]=[4,4,4,4,1,-1]
    ####
    states[142,:]=[1,-1,-1,-1,2,-1]
    states[143,:]=[2,-1,-1,-1,2,-1]
    states[144,:]=[3,-1,-1,-1,2,-1]
    states[145,:]=[4,-1,-1,-1,2,-1]
    #
    states[146,:]=[1,1,-1,-1,2,-1]
    states[147,:]=[2,1,-1,-1,2,-1]
    states[148,:]=[2,2,-1,-1,2,-1]
    states[149,:]=[2,3,-1,-1,2,-1]
    states[150,:]=[2,4,-1,-1,2,-1]
    states[151,:]=[3,1,-1,-1,2,-1]
    states[152,:]=[3,2,-1,-1,2,-1]
    states[153,:]=[3,3,-1,-1,2,-1]
    states[154,:]=[3,4,-1,-1,2,-1]
    states[155,:]=[4,1,-1,-1,2,-1]
    states[156,:]=[4,2,-1,-1,2,-1]
    states[157,:]=[4,3,-1,-1,2,-1]
    states[158,:]=[4,4,-1,-1,2,-1]
    #
    states[159,:]=[1,1,1,-1,2,-1]
    states[160,:]=[2,2,1,-1,2,-1]
    states[161,:]=[2,2,2,-1,2,-1]
    states[162,:]=[2,2,3,-1,2,-1]
    states[163,:]=[2,3,1,-1,2,-1]
    states[164,:]=[2,3,2,-1,2,-1]
    states[165,:]=[2,3,3,-1,2,-1]
    states[166,:]=[2,4,1,-1,2,-1]
    states[167,:]=[2,4,2,-1,2,-1]
    states[168,:]=[2,4,3,-1,2,-1]
    states[169,:]=[3,2,1,-1,2,-1]
    states[170,:]=[3,2,2,-1,2,-1]
    states[171,:]=[3,3,2,-1,2,-1]
    states[172,:]=[3,3,3,-1,2,-1]
    states[173,:]=[4,1,1,-1,2,-1]
    states[174,:]=[4,2,1,-1,2,-1]
    states[175,:]=[4,2,2,-1,2,-1]
    states[176,:]=[4,3,2,-1,2,-1]
    states[177,:]=[4,4,1,-1,2,-1]
    states[178,:]=[4,4,2,-1,2,-1]
    states[179,:]=[4,4,3,-1,2,-1]
    states[180,:]=[4,4,4,-1,2,-1]
    #
    states[181,:]=[1,1,1,1,2,-1]
    states[182,:]=[2,1,1,1,2,-1]
    states[183,:]=[2,2,1,1,2,-1]
    states[184,:]=[2,2,2,1,2,-1]
    states[185,:]=[2,2,2,2,2,-1]
    states[186,:]=[3,1,1,1,2,-1]
    states[187,:]=[3,2,1,1,2,-1]
    states[188,:]=[3,2,2,1,2,-1]
    states[189,:]=[3,2,2,2,2,-1]
    states[190,:]=[3,3,1,1,2,-1]
    states[191,:]=[3,3,2,1,2,-1]
    states[192,:]=[3,3,2,2,2,-1]
    states[193,:]=[3,3,3,1,2,-1]
    states[194,:]=[3,3,3,2,2,-1]
    states[195,:]=[3,3,3,3,2,-1]
    states[196,:]=[4,1,1,1,2,-1]
    states[197,:]=[4,2,1,1,2,-1]
    states[198,:]=[4,2,2,1,2,-1]
    states[199,:]=[4,2,2,2,2,-1]
    states[200,:]=[4,3,1,1,2,-1]
    states[201,:]=[4,3,2,1,2,-1]
    states[202,:]=[4,3,2,2,2,-1]
    states[203,:]=[4,3,3,1,2,-1]
    states[204,:]=[4,3,3,2,2,-1]
    states[205,:]=[4,3,3,3,2,-1]
    states[206,:]=[4,4,1,1,2,-1]
    states[207,:]=[4,4,2,1,2,-1]
    states[208,:]=[4,4,2,2,2,-1]
    states[209,:]=[4,4,3,2,2,-1]
    states[210,:]=[4,4,3,3,2,-1]
    states[211,:]=[4,4,4,4,2,-1]
    ####
    states[212,:]=[1,-1,-1,-1,2,1]
    states[213,:]=[2,-1,-1,-1,2,1]
    states[214,:]=[3,-1,-1,-1,2,1]
    states[215,:]=[4,-1,-1,-1,2,1]
    #
    states[216,:]=[1,1,-1,-1,2,1]
    states[217,:]=[2,1,-1,-1,2,1]
    states[218,:]=[2,2,-1,-1,2,1]
    states[219,:]=[2,3,-1,-1,2,1]
    states[220,:]=[2,4,-1,-1,2,1]
    states[221,:]=[3,1,-1,-1,2,1]
    states[222,:]=[3,2,-1,-1,2,1]
    states[223,:]=[3,3,-1,-1,2,1]
    states[224,:]=[3,4,-1,-1,2,1]
    states[225,:]=[4,1,-1,-1,2,1]
    states[226,:]=[4,2,-1,-1,2,1]
    states[227,:]=[4,3,-1,-1,2,1]
    states[228,:]=[4,4,-1,-1,2,1]
    #
    states[229,:]=[1,1,1,-1,2,1]
    states[230,:]=[2,2,1,-1,2,1]
    states[231,:]=[2,2,2,-1,2,1]
    states[232,:]=[2,2,3,-1,2,1]
    states[233,:]=[2,3,1,-1,2,1]
    states[234,:]=[2,3,2,-1,2,1]
    states[235,:]=[2,3,3,-1,2,1]
    states[236,:]=[2,4,1,-1,2,1]
    states[237,:]=[2,4,2,-1,2,1]
    states[238,:]=[2,4,3,-1,2,1]
    states[239,:]=[3,2,1,-1,2,1]
    states[240,:]=[3,2,2,-1,2,1]
    states[241,:]=[3,3,2,-1,2,1]
    states[242,:]=[3,3,3,-1,2,1]
    states[243,:]=[4,1,1,-1,2,1]
    states[244,:]=[4,2,1,-1,2,1]
    states[245,:]=[4,2,2,-1,2,1]
    states[246,:]=[4,3,2,-1,2,1]
    states[247,:]=[4,4,1,-1,2,1]
    states[248,:]=[4,4,2,-1,2,1]
    states[249,:]=[4,4,3,-1,2,1]
    states[250,:]=[4,4,4,-1,2,1]
    #
    states[251,:]=[1,1,1,1,2,1]
    states[252,:]=[2,1,1,1,2,1]
    states[253,:]=[2,2,1,1,2,1]
    states[254,:]=[2,2,2,1,2,1]
    states[255,:]=[2,2,2,2,2,1]
    states[256,:]=[3,1,1,1,2,1]
    states[257,:]=[3,2,1,1,2,1]
    states[258,:]=[3,2,2,1,2,1]
    states[259,:]=[3,2,2,2,2,1]
    states[260,:]=[3,3,1,1,2,1]
    states[261,:]=[3,3,2,1,2,1]
    states[262,:]=[3,3,2,2,2,1]
    states[263,:]=[3,3,3,1,2,1]
    states[264,:]=[3,3,3,2,2,1]
    states[265,:]=[3,3,3,3,2,1]
    states[266,:]=[4,1,1,1,2,1]
    states[267,:]=[4,2,1,1,2,1]
    states[268,:]=[4,2,2,1,2,1]
    states[269,:]=[4,2,2,2,2,1]
    states[270,:]=[4,3,1,1,2,1]
    states[271,:]=[4,3,2,1,2,1]
    states[272,:]=[4,3,2,2,2,1]
    states[273,:]=[4,3,3,1,2,1]
    states[274,:]=[4,3,3,2,2,1]
    states[275,:]=[4,3,3,3,2,1]
    states[276,:]=[4,4,1,1,2,1]
    states[277,:]=[4,4,2,1,2,1]
    states[278,:]=[4,4,2,2,2,1]
    states[279,:]=[4,4,3,2,2,1]
    states[280,:]=[4,4,3,3,2,1]
    states[281,:]=[4,4,4,4,2,1]
    ####
    states[292,:]=[1,-1,-1,-1,2,2]
    states[293,:]=[2,-1,-1,-1,2,2]
    states[294,:]=[3,-1,-1,-1,2,2]
    states[295,:]=[4,-1,-1,-1,2,2]
    #
    states[296,:]=[1,1,-1,-1,2,2]
    states[297,:]=[2,1,-1,-1,2,2]
    states[298,:]=[2,2,-1,-1,2,2]
    states[299,:]=[2,3,-1,-1,2,2]
    states[300,:]=[2,4,-1,-1,2,2]
    states[301,:]=[3,1,-1,-1,2,2]
    states[302,:]=[3,2,-1,-1,2,2]
    states[303,:]=[3,3,-1,-1,2,2]
    states[304,:]=[3,4,-1,-1,2,2]
    states[305,:]=[4,1,-1,-1,2,2]
    states[306,:]=[4,2,-1,-1,2,2]
    states[307,:]=[4,3,-1,-1,2,2]
    states[308,:]=[4,4,-1,-1,2,2]
    #
    states[309,:]=[1,1,1,-1,2,2]
    states[310,:]=[2,2,1,-1,2,2]
    states[311,:]=[2,2,2,-1,2,2]
    states[312,:]=[2,2,3,-1,2,2]
    states[313,:]=[2,3,1,-1,2,2]
    states[314,:]=[2,3,2,-1,2,2]
    states[315,:]=[2,3,3,-1,2,2]
    states[316,:]=[2,4,1,-1,2,2]
    states[317,:]=[2,4,2,-1,2,2]
    states[318,:]=[2,4,3,-1,2,2]
    states[319,:]=[3,2,1,-1,2,2]
    states[320,:]=[3,2,2,-1,2,2]
    states[321,:]=[3,3,2,-1,2,2]
    states[322,:]=[3,3,3,-1,2,2]
    states[323,:]=[4,1,1,-1,2,2]
    states[324,:]=[4,2,1,-1,2,2]
    states[325,:]=[4,2,2,-1,2,2]
    states[326,:]=[4,3,2,-1,2,2]
    states[327,:]=[4,4,1,-1,2,2]
    states[328,:]=[4,4,2,-1,2,2]
    states[329,:]=[4,4,3,-1,2,2]
    states[330,:]=[4,4,4,-1,2,2]
    #
    states[331,:]=[1,1,1,1,2,2]
    states[332,:]=[2,1,1,1,2,2]
    states[333,:]=[2,2,1,1,2,2]
    states[334,:]=[2,2,2,1,2,2]
    states[335,:]=[2,2,2,2,2,2]
    states[336,:]=[3,1,1,1,2,2]
    states[337,:]=[3,2,1,1,2,2]
    states[338,:]=[3,2,2,1,2,2]
    states[339,:]=[3,2,2,2,2,2]
    states[340,:]=[3,3,1,1,2,2]
    states[341,:]=[3,3,2,1,2,2]
    states[342,:]=[3,3,2,2,2,2]
    states[343,:]=[3,3,3,1,2,2]
    states[344,:]=[3,3,3,2,2,2]
    states[345,:]=[3,3,3,3,2,2]
    states[346,:]=[4,1,1,1,2,2]
    states[347,:]=[4,2,1,1,2,2]
    states[348,:]=[4,2,2,1,2,2]
    states[349,:]=[4,2,2,2,2,2]
    states[350,:]=[4,3,1,1,2,2]
    states[351,:]=[4,3,2,1,2,2]
    states[352,:]=[4,3,2,2,2,2]
    states[353,:]=[4,3,3,1,2,2]
    states[354,:]=[4,3,3,2,2,2]
    states[355,:]=[4,3,3,3,2,2]
    states[356,:]=[4,4,1,1,2,2]
    states[357,:]=[4,4,2,1,2,2]
    states[358,:]=[4,4,2,2,2,2]
    states[359,:]=[4,4,3,2,2,2]
    states[350,:]=[4,4,3,3,2,2]
    states[351,:]=[4,4,4,4,2,2]
    return states
end
function assignDataToStates( data::Array{T1,3}, states::Array{T2,2} , Err::T3 ) where { T1 <: Real, T2 <: Real , T3 <: Bool }
    nb_data_point=size(data)[1]
    nb_series = size(data)[2]
    dim_data = size(data)[3]
    nb_states = size(states)[1]
    state_matrix=ones(Int, nb_data_point, nb_series )*(-1)
    count_states=zeros( nb_states )
    unused=0
    for i=1:nb_data_point
        print("Assigning data to states - Progress: ",i/nb_data_point*100,"%\n")
        for j=1:nb_series
            for l=1:nb_states
                d=0
                for k=1:dim_data
                    d+= ( data[i,j,k] - states[l,k])*(data[i,j,k] - states[l,k])
                end
                if d == 0
                    count_states[l] += 1
                    state_matrix[i,j] = l
                    break
                end
            end
            if state_matrix[i,j] == -1
                if Err
                    print("Out ",i," ",j," ")
                    for k=1:dim_data
                        print(k," ",data[i,j,k]," ")
                    end
                end
                print("\n")
                unused += 1
            end
        end
    end
    return state_matrix, count_states/sum(count_states)*100, unused/(unused+sum(count_states))*100
end
function isolateSignificantStates( old_states::Array{T1,2}, percent_states::Vector{T2}, cut_off_states::T3 ) where { T1 <: Real, T2 <: Real , T3 <: Real }
    old_nb_states=size(old_states)[1]
    dim=size(old_states)[2]
    new_nb_states=0
    states_kept=zeros(0,dim)
    for i=1:old_nb_states
        if percent_states[i] > cut_off_states
            states_kept=[ states_kept ; transpose( old_states[i,:]) ]
        end
    end
    return states_kept
end
function transitionMatrix( states::Array{T1,2}, state_matrix::Array{T2,2}, min_lag::T3, max_lag::T4, d_lag::T5) where { T1 <: Real, T2 <: Real, T3 <: Real, T4<:Int, T5 <: Int }

    nb_states=size(states)[1]
    nb_data_point=size(state_matrix)[1]
    nb_series = size(state_matrix)[2]
    nb_lag_points=Int(trunc((max_lag-min_lag)/d_lag))

    states_transition_probability=zeros(Float64,nb_states,nb_states,nb_lag_points)

    # Chappman Kolmogorov test
    count_lag=1
    for lag=min_lag:d_lag:max_lag-1
        print("Chappman Kolmogorov Test - Progress: ",lag/max_lag*100,"%\n")
        for i=1:nb_series
            for j=lag+1:nb_data_point
                if  state_matrix[j-lag,i] == -1 ||  state_matrix[j,i] == -1
                    continue
                end
                states_transition_probability[ state_matrix[j-lag,i], state_matrix[j,i], count_lag ] += 1
            end
        end
        count_lag += 1
    end

    # Normalization
    for lag=1:nb_lag_points
        for i=1:nb_states
            states_transition_probability[i,:,lag] /= sum( states_transition_probability[i,:,lag] )
        end
    end

    return states_transition_probability
end
function chappmanKormologov( transition_matrix::Array{T1,3} ) where { T1 <: Real }
    nb_states   = size(transition_matrix)[1]
    nb_lag_time = size(transition_matrix)[3]
    nb_lag_compare = Int(trunc(nb_lag_time/2))
    transition_matrix_kolmo= zeros(nb_states,nb_states,nb_lag_compare)
    for lag=1:nb_lag_compare
        for i=1:nb_states
            for j=1:nb_states
                for k=1:nb_states
                    transition_matrix_kolmo[i,j,lag] += transition_matrix[i,k,lag]*transition_matrix[k,j,lag]
                end
            end
        end
    end
    return transition_matrix_kolmo
end
function writeStates( file::T1 , states::Array{T2,2}, percent::Vector{T3}) where { T1 <: AbstractString, T2 <: Real, T3 <: Real }
    file_out=open(file,"w")
    n_dim = size( states)[2]
    nb_states=size(states)[1]
    for i=1:nb_states
        for j=1:n_dim
            write(file_out,string(states[i,j]," "))
        end
        write(file_out,string(percent[i],"\n"))
    end
    close(file_out)
    return
end

# Folder for data
folder_base="/media/moogmt/Stock/CO2/AIMD/Liquid/PBE-MT/"
#folder_base="/home/moogmt/CO2/CO2_AIMD/"

# Thermo data
Volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.25,9.3,9.35,9.375,9.4,9.5,9.8,10.0]
Temperatures=[2000,2500,3000]
Cut_Off=[1.75]

# Number of atoms
nbC=32
nbO=nbC*2

cut_off_bond = 1.75

min_lag=1
max_lag=5001
d_lag=5
unit=0.005




T=3000
V=9.0

folder_in=string(folder_base,V,"/",T,"K/")
file=string(folder_in,"TRAJEC_wrapped.xyz")

folder_out=string(folder_in,"Data/")
#folder_out=string(folder_in)


print("Computing Data\n")
traj=filexyz.readFastFile(file)
cell=cell_mod.Cell_param(V,V,V)

states=defineStatesExtendedCoordinances()
data=buildCoordinationMatrix( traj , cell , cut_off_bond )
state_matrix, percent, unused_percent = assignDataToStates( data , states, true )

states = isolateSignificantStates( states, percent, 0.0 )
state_matrix, percent, unused_percent = assignDataToStates( data , states , false)
writeStates(string(folder_out,"markov_initial_states.dat"),states,percent)

cut_off_states = 0.1
states = isolateSignificantStates( states, percent, cut_off_states )
state_matrix, percent, unused_percent = assignDataToStates( data , states , false)
transition_matrix = transitionMatrix( states, state_matrix, min_lag, max_lag, d_lag )
transition_matrix_CK = chappmanKormologov( transition_matrix )
writeStates(string(folder_out,"markov_final_states-",percent,".dat"),states,percent)

nb_states=size(states)[1]

for j=1:nb_states
    file_out=open(string(folder_out,"C_markov_CK_test-",cut_off_bond,"-",j,"-part1.dat"),"w")
    for i=1:2:size(transition_matrix)[3]
        write(file_out,string(i*unit*d_lag," "))
        for k=1:nb_states
            write(file_out,string(transition_matrix[j,k,i]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
end
for j=1:nb_states
    file_out=open(string(folder_out,"C_markov_CK_test-",cut_off_bond,"-",j,"-part2.dat"),"w")
    for i=1:size(transition_matrix_CK)[3]
        write(file_out,string(2*i*unit*d_lag," "))
        for k=1:nb_states
            write(file_out,string(transition_matrix_CK[j,k,i]," "))
        end
        write(file_out,string("\n"))
    end
    close(file_out)
end


# nb_states=size(transition_matrix)[1]
# for state=1:nb_states
#     if percent[state] > 0.05
#         print("Progress: ",state/nb_states*100,"%\n")
#         distances=[]
#         for step_sim=1:nb_steps
#             for carbon=1:nbC
#                 if state_matrix[step_sim,carbon] == state
#                     min_dist=V
#                     for atom=1:96
#                         dist=cell_mod.distance(traj[step_sim],cell,carbon,atom)
#                         if atom == carbon
#                             continue
#                         end
#                         if min_dist > dist
#                             min_dist = dist
#                         end
#                     end
#                     push!(distances,min_dist)
#                 end
#             end
#         end
#         min_v=0.8
#         max_v=2
#         nb_points=200
#         delta=(max_v-min_v)/nb_points
#         hist1D=zeros(nb_points)
#         for j=1:size(distances)[1]
#             for i=1:nb_points
#                 if distances[j] > min_v + (i-1)*delta && distances[j] < min_v + i*delta
#                     hist1D[i] += 1
#                 end
#             end
#         end
#         hist1D /= sum(hist1D)
#         file_distance=open(string(folder_out,"distances-",coordinances[state],"-",state,".dat"),"w")
#         for i=1:nb_points
#             write(file_distance,string(min_v+i*delta," ",hist1D[i],"\n"))
#         end
#         close(file_distance)
#     end
# end
