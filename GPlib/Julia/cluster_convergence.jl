# Loading file
include("contactmatrix.jl")

functionnals=["PBE-MT","BLYP","PBE-Godecker"]
volumes=[8.82,9.0,9.05,9.1,9.15,9.2,9.3,9.325,9.35,9.375,9.4,9.5,9.8]
temperatures=[2000,2500,3000,3500]
cut_off=[1.6,1.7,1.8]

V=volumes[1]
T=temperatures[3]
func=functionnals[1]
stride=10
co=cut_off[1]

nbC=32
nbO=2*nbC

function kmenoid_clustering( n_clusters, n_structures, distance_matrix )
    cluster_sizes=zeros(n_clusters)
    cluster_centers=zeros(n_clusters)
    cluster_members=zeros(n_clusters,number_structures)
    v=zeros(n_structures)
    dl=zeros(n_structures)
    good=ones(n_structures)
    ass=zeros(n_structures)
    found=false

    # Find the cluster centers
    cluster_centers[1]=Int(trunc(rand()*n_structure))+1
    good[cluster_centers[1]] = 0
    for i=2:n_clusters
        for k=1:nb_structures
            dmin=1.0
            for j=1:i-1
                dtmp=distance_matrix[k,cluster_centers[j]]
                if dtmp < dmin
                    dtmp=dmin
                end
            end
            proba[k] = dmin^2
        end
        dtmp=sum(proba)
        proba/=dtmp
        for k=2:nb_steps
            proba[k] += proba[k-1]
        end
        found=false
        while ! found
            r=rand()
            if r < proba[1]
                if good[1] == 0
                    continue
                else
                    cluster_center[i]=1
                    good[1]=0
                    found = true
                end
            else
                for k=2:n_structures
                    if good[k] == 0
                        continue
                    else
                        if r > proba[k-1] && r < proba[k]
                            cluster_centers[i]=k
                            good[k]=0
                            found=true
                            break
                        end
                    end
                end
            end
        end
    end

end
