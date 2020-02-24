#INFORMATION
#The input of this script is a pandas dataframe containing the ab-initio data_path
#Under the PARAMETERS section the descriptor and training method can be defined
#The output is a dataframe that contains the predicted energies as well as the trained model
#(the trained model is accesible trough the metadata)


#Python Imports

#ASE imports
from ase import Atoms

#Other
import numpy as np
import pandas as pd
import os
import random
import dill as pickle
#import pickle
pd.set_option('max_colwidth', 20)
pd.set_option('max_rows', 9)
pd.set_option('max_columns',9)

from sklearn.model_selection import cross_val_score
from sklearn.gaussian_process.kernels import RBF
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import KMeans

############################CODE###############################################

def select_descriptor(data, descriptor_metadata)  :
    print('creating descriptor....')
    from dscribe.descriptors import CoulombMatrix
    if descriptor_metadata['descriptor_type'] == 'Coulomb_matrix':
        scaling = descriptor_metadata['scaling']
        descriptor = CoulombMatrix(
            n_atoms_max=data['configuration'][0].get_number_of_atoms(),
            flatten=True,
            permutation = descriptor_metadata['permutation'])

        rounding=5
        l = data['configuration'].tolist()
        features = [descriptor.create(l[i]) for i in range(len(l))]
        features = np.round(np.array(features),rounding)
        data = pd.concat([data,
                            pd.DataFrame(features)],
                            axis=1)

        descriptor_metadata['descriptor'] = descriptor

        features = [i for i in range(descriptor_metadata['descriptor'].get_number_of_features())]
        labels = ['energy']
        features_and_labels = labels + features

        descriptor_metadata['features'] = features
        descriptor_metadata['labels'] = labels
        descriptor_metadata['features_and_labels'] = features_and_labels

        data['metadata'].at[1] = descriptor_metadata
        return data

    elif descriptor_metadata['descriptor_type'] == 'PIV':
        def switching_OO(x):
            n = 8
            d0 = 4.5
            return x
        def switching_OH_plus(x):
            n = 8
            d0 = 2.3
            #return 1/(1+(x/d0)**n)
            return x
        def switching_HH(x):
            return x
        def switching_HH_plus(x):
            return switching_HH(x)
        def switching_OH(x):
            return switching_HH(x)

        def PIV(configuration):
            distances = configuration.get_all_distances()
            np.fill_diagonal(distances, 0)
            distances[np.tril_indices(distances.shape[0], -1)] = 0

            OO = distances[0,1].flatten()
            OH_plus = distances[0:2, 2].flatten()
            OH= distances[0:2, 3:7].flatten()
            HH_plus = distances[2, 3:7].flatten()
            HH = distances[3:6, 4:7].flatten()
            HH = np.delete(HH, (3,6,7))

            OO = switching_OO(OO)
            OH_plus = switching_OH_plus(OH_plus)
            OH = switching_OH(OH)
            HH_plus = switching_HH_plus(HH_plus)
            HH = switching_HH(HH)

            OO = np.sort(OO)
            OH_plus = np.sort(OH_plus)
            OH = np.sort(OH)
            HH_plus = np.sort(HH_plus)
            HH = np.sort(HH)

            PIV = np.concatenate((OO, OH_plus, OH, HH_plus, HH), axis=None)
            return PIV

        l = data['configuration'].tolist()
        features = [PIV(l[i]) for i in range(len(l))]
        data = pd.concat([data,
                        pd.DataFrame(features)],
                        axis=1)

        features = [i for i in range(21)]
        labels = ['energy']
        features_and_labels = labels + features
        descriptor_metadata['features'] = features
        descriptor_metadata['labels'] = labels
        descriptor_metadata['features_and_labels'] = features_and_labels
        descriptor_metadata['switching_OO'] = switching_OO
        descriptor_metadata['switching_OH_plus'] = switching_OH_plus
        descriptor_metadata['switching_OH'] = switching_OH
        descriptor_metadata['switching_HH_plus'] = switching_HH_plus
        descriptor_metadata['switching_HH'] = switching_HH
        descriptor_metadata['descriptor'] = PIV

        data['metadata'].at[1] = descriptor_metadata
        return data

    elif descriptor_metadata['descriptor_type'] == 'CM_with_PIV_sorting':
        def PIV(configuration):
            def w_diag(x):
                return 0.5*np.power(x,2.4)
            distances = configuration.get_all_distances()
            for i in range(len(distances)):
                distances[i,i] = 1
            def switching_function(x):
                return 1/x
            distances = switching_function(distances)
            distances[0,0] = w_diag(8)
            distances[1,1] = w_diag(8)
            for i in range(2,7):
                distances[i,i] = w_diag(1)
            distances[0,1] = distances[0,1]*64
            distances[1,0] = distances[1,0]*64
            distances[0:2, 2:7] = distances[0:2, 2:7]*8
            distances[2:7, 0:2] = distances[2:7, 0:2]*8
            OO = distances[0:2,0:2].flatten()
            HH1 = distances[0:7, 2:7].flatten()
            HH2 = distances[2:7, 0:2].flatten()
            HH = np.append(HH2,HH1)

            OO = np.sort(OO)
            HH = np.sort(HH)

            PIV = np.concatenate((OO, HH), axis=None)
            return PIV

        l = data['configuration'].tolist()
        features = [PIV(l[i]) for i in range(len(l))]
        data = pd.concat([data,
                        pd.DataFrame(features)],
                        axis=1)

        features = [i for i in range(49)]
        labels = ['energy']
        features_and_labels = labels + features

        descriptor_metadata['features'] = features
        descriptor_metadata['labels'] = labels
        descriptor_metadata['features_and_labels'] = features_and_labels

        data['metadata'].at[1] = descriptor_metadata
        return data

        data['metadata'].at[1] = descriptor_metadata
        return data

    elif descriptor_metadata['descriptor_type'] == 'PIV_without_H_plus':
        def PIV(configuration):
            distances = configuration.get_all_distances()
            distances[np.tril_indices(distances.shape[0], -1)] = 0
            distances
            OO = distances[0:2,0:2].flatten()
            HH1 = distances[0:7, 2:7].flatten()
            HH2 = distances[2:7, 0:2].flatten()
            HH = np.append(HH2,HH1)

            OO = np.sort(OO)
            HH = np.sort(HH)
            OO = OO[OO != 0]
            HH = HH[HH != 0]
            PIV = np.concatenate((OO, HH), axis=None)
            return PIV

        l = data['configuration'].tolist()
        features = [PIV(l[i]) for i in range(len(l))]
        data = pd.concat([data,
                        pd.DataFrame(features)],
                        axis=1)

        features = [i for i in range(21)]
        labels = ['energy']
        features_and_labels = labels + features

        descriptor_metadata['features'] = features
        descriptor_metadata['labels'] = labels
        descriptor_metadata['features_and_labels'] = features_and_labels

        data['metadata'].at[1] = descriptor_metadata
        return data

        data['metadata'].at[1] = descriptor_metadata
        return data

    elif descriptor_metadata['descriptor_type'] == 'PIV_with_CM_diagonal_and_weighting':
        def PIV(configuration):
            def w_diag(x):
                return 0.5*np.power(x,2.4)
            #This PIV has 28 elements
            distances = configuration.get_all_distances()
            for i in range(len(distances)):
                    distances[i,i] = 1
            def switching_function(x):
                    return 1/x
            distances = switching_function(distances)
            distances[0,0] = w_diag(8)
            distances[1,1] = w_diag(8)
            for i in range(2,7):
                distances[i,i] = w_diag(1)
            distances[np.tril_indices(distances.shape[0], -1)] = 0
            distances
            OO = distances[0:2,0:2].flatten()
            HH1 = distances[0:7, 2:7].flatten()
            HH2 = distances[2:7, 0:2].flatten()
            HH = np.append(HH2,HH1)

            OO = np.sort(OO)
            HH = np.sort(HH)

            OO = OO[OO != 0]
            HH = HH[HH != 0]

            PIV = np.concatenate((OO, HH), axis=None)
            return PIV

        l = data['configuration'].tolist()
        features = [PIV(l[i]) for i in range(len(l))]
        data = pd.concat([data,
                        pd.DataFrame(features)],
                        axis=1)

        features = [i for i in range(28)]
        labels = ['energy']
        features_and_labels = labels + features

        descriptor_metadata['features'] = features
        descriptor_metadata['labels'] = labels
        descriptor_metadata['features_and_labels'] = features_and_labels

        data['metadata'].at[1] = descriptor_metadata
        return data

        data['metadata'].at[1] = descriptor_metadata
        return data

    else:
        print('NO DESCRIPTOR SELECTED!')
        descriptor_metadata['descriptor_type'] = 'NO_DESCRIPTOR_SELECTED'
        data['metadata'].at[1] = descriptor_metadata
        return data

def select_training_points(data, sampling_metadata):
    print('selecting training points....')

    features = data['metadata'].at[1]['features']
    labels = data['metadata'].at[1]['labels']
    features_and_labels = data['metadata'].at[1]['features_and_labels']


    number_of_training_points = sampling_metadata['number_of_training_points']

    if sampling_metadata['sampling_method'] == 'random':
        index_train = data.sample(number_of_training_points, random_state=0).index
        train = data.loc[index_train]
        data['train'] = False
        data.at[index_train, 'train'] = True
    elif sampling_metadata['sampling_method'] == 'Kmeans':
        number_of_clusters = sampling_metadata['number_of_clusters']
        kmeans = KMeans(n_clusters=number_of_clusters, random_state=0).fit(data[features])
        data['cluster'] = kmeans.predict(data[features])
        data['train'] = False
        for c in range(number_of_clusters):
            points_to_sample_per_cluster = int((number_of_training_points/number_of_clusters))
            index = data.index[(data['cluster'] == c)]
            r = random.sample(list(index.values), points_to_sample_per_cluster)
            data.set_value(r, 'train', True)
    else:
        sampling_metadata['sampling_method'] = 'no sampling method selected!'
        print('no sampling method selected!')

    data['metadata'].at[2] = sampling_metadata

    return data

def select_ML_model(data, model_metadata):
    print('selecting ML model....')

    features = data['metadata'].at[1]['features']
    labels = data['metadata'].at[1]['labels']
    features_and_labels = data['metadata'].at[1]['features_and_labels']

    gs = GridSearchCV(cv=3,
                      param_grid=model_metadata['param_grid'],
                      estimator=model_metadata['estimator'],
                      scoring=model_metadata['scoring'],
                      n_jobs=-1,
                      verbose=4,
                      return_train_score=True)

    gs.fit(data[features][data['train']==True],
           data[labels][data['train']==True])

    data.at[3, 'metadata'] = model_metadata
    data.at[3, 'metadata']['scores'] = gs.cv_results_
    data.at[3, 'metadata']['best_model'] = gs.best_estimator_
    data.at[3, 'metadata']['best_model_fitted'] = gs.best_estimator_.fit(
                                                            data[features][data['train']==True],
                                                            data[labels][data['train']==True])
    return data

def make_prediction(data):
    print('predicting the energies for all data points...')

    features = data['metadata'].at[1]['features']
    labels = data['metadata'].at[1]['labels']
    features_and_labels = data['metadata'].at[1]['features_and_labels']
    model = data.at[3, 'metadata']['best_model_fitted']
    data['prediction'] = model.predict(data[features])
    return data

def create_model(data,
                 descriptor_metadata,
                 sampling_metadata,
                 model_metadata,
                 save_path):

    print('creating model for parameters:')
    print('-----------------------------------------------------------')
    print(data.at[0, 'metadata']['system_name'])
    print(data.at[0, 'metadata']['temperature_of_simulation'])
    print(data.at[0, 'metadata']['energy_calculation_method'])
    print(descriptor_metadata['descriptor_type'])
    print('#Training points: ' +str(sampling_metadata['number_of_training_points']))
    print('Sampling method: ' +sampling_metadata['sampling_method'])
    print('-----------------------------------------------------------')
    print('model will be saved to: '+save_path)
    print('-----------------------------------------------------------')

    data = select_descriptor(data, descriptor_metadata)
    data = select_training_points(data, sampling_metadata)
    data = select_ML_model(data, model_metadata)
    data = make_prediction(data)

    print('saving dataframe to '+ save_path)
    pickle.dump(data, open(save_path, 'wb'))
    return data

##################PARAMETERS###############################################

for n in [1000]:
    data_path = '../1_data_preparation/prepared_data/CCFF_zundel_50K_cleaned'
    #data_import = pickle.load(open(data_path, 'rb'))
    data = pickle.load(open(data_path, 'rb'))

    #data = data_import.copy()

    descriptor_metadata = {
    #options:
    #CM_with_PIV_sorting
    #PIV_with_CM_diagonal_and_weighting
    #PIV_without_H_plus
    #PIV
    #Coulomb_matrix
        'descriptor_type' : 'Coulomb_matrix',
        'permutation' : 'sorted_l2',
        'scaling' : False
    }

    sampling_metadata = {
        'number_of_training_points':n,
        'sampling_method': 'random',
        'number_of_clusters': 0
    }

    model_metadata = {
        'estimator' : KernelRidge(),
        'cv' : 3,
        'scoring' : 'neg_mean_squared_error',
        'param_grid':{
            'kernel':[RBF(length_scale =l) for l in np.arange(5,82,4)],
            'alpha':[10**x for x in range(-20,-6)]
        }
    }

    save_folder = 'trained_models'
    save_name = descriptor_metadata['descriptor_type']+"_"+ \
                     data['metadata'][0]['temperature_of_simulation'] +"_"+ \
                     str(sampling_metadata['number_of_training_points'])+"_"+ \
                     sampling_metadata['sampling_method']
    save_path = os.path.join(save_folder, save_name)
    if os.path.isfile(save_path) == False:
        create_model(data,
                     descriptor_metadata,
                     sampling_metadata,
                     model_metadata,
                     save_path)
    else:
        print('-----------------------------------------------------------')
        print('file: '+ save_path+' already exists!')
        print('-----------------------------------------------------------')

print('###################################################')
print('#                  JOB COMPLETED!                 #')
print('###################################################')
