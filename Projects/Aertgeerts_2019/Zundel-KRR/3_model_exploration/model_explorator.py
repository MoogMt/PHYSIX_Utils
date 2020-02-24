#PYTHON IMPORTS

#ASE imports
import ase
from ase.build import molecule
from ase import Atoms
from dscribe.descriptors import CoulombMatrix

#Other
import pandas as pd
import os
from dscribe.descriptors import SOAP
from sklearn.manifold import MDS
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import random
import matplotlib.cm as colormap
from IPython.display import clear_output
%matplotlib inline
from ipywidgets import IntProgress
from scipy.optimize import minimize
from IPython.display import display
from numpy import linalg as LA
import dill as pickle
#import pickle

#sklearn imports
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.gaussian_process.kernels import RBF
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.metrics import r2_score

pd.set_option('max_colwidth', 20)
pd.set_option('max_rows', 9)
pd.set_option('max_columns',9)

#Supress warnings
#import warnings; warnings.simplefilter('ignore')

#IMPORT & CALCULATION FUNCTIONS
def import_model(model_path):
    data = pickle.load(open(model_path, 'rb'))
    data['index'] = data.index
    return data

def model_error(data):
    rounding = 4
    train = data.loc[data['train'] == True]
    test = data.loc[data['train'] == False]
    MAE_train = mean_absolute_error(train['energy'],train['prediction'])
    MAE_test = mean_absolute_error(test['energy'],test['prediction'])
    MSE_train = mean_squared_error(train['energy'],train['prediction'])
    MSE_test = mean_squared_error(test['energy'],test['prediction'])
    R2_train = r2_score(train['energy'],train['prediction'])
    R2_test = r2_score(test['energy'],test['prediction'])

    errors = {'MAE_train' : np.round(MAE_train, rounding),
              'MAE_test' : np.round(MAE_test, rounding),
              'MSE_train' : np.round(MSE_train, rounding),
              'MSE_test' : np.round(MSE_test, rounding),
             'R2_train' : np.round(R2_train, rounding),
             'R2_test' : np.round(R2_test, rounding)}

    return errors

def add_distances(data):
    #Add dOH and dOO for plot
    data['dO1O2'] = [data['configuration'].iloc[n].get_all_distances()[0][1] for n in range(len(data))]
    data['dHO1'] = [data['configuration'].iloc[n].get_all_distances()[0][2] for n in range(len(data))]
    return data

def import_3D_downfold(path='../downfold/downfold_data/downfolded_3D'):
    # IMPORT 3D DOWNFOLD DATA
    down_3D = pd.read_csv(path)
    down_3D['config'] = [np.matrix(down_3D['position'].iloc[i]).reshape(7,3) for i in range(len(down_3D))]
    down_3D = down_3D.rename(columns={'dOO':'dO1O2', 'dOH':'dHO1'})
    down_3D = down_3D.drop(['position', 'success'], axis=1)
    return down_3D

def import_2D_downfold(path='../Downfold/downfold_data/downfolded_2D'):
    # IMPORT 2D DOWNFOLD DATA
    down_2D = pd.read_csv(path)
    down_2D['config'] = [np.matrix(down_2D['position'].iloc[i]).reshape(7,3) for i in range(len(down_2D))]
    down_2D = down_2D.rename(columns={'dOO':'dO1O2'})
    down_2D = down_2D.drop(['position', 'success'], axis=1)
    return down_2D

def ML_potential(config, data):
    model = data['metadata'][3]['best_model_fitted']
    if data['metadata'][1]['descriptor_type'] == 'Coulomb_matrix':
        descriptor = CoulombMatrix(
        n_atoms_max=7,
        flatten=True,
        permutation = 'sorted_l2')
        x = Atoms('O2H5',positions=config)
        X = descriptor.create(x)
        energy = model.predict(X)[0][0]
        return energy

    if data['metadata'][1]['descriptor_type'] == 'PIV':
        descriptor = data['metadata'][1]['descriptor']
        x = Atoms('O2H5', positions=config)
        X = descriptor(x)
        energy = model.predict(X)[0][0]
        return energy

def downfold_2D_prediction(data, down_2D):
    down_2D['prediction'] = [ML_potential(down_2D['config'][i], data) for i in range(len(down_2D))]
    return down_2D

def downfold_3D_prediction(data, down_3D):
    down_3D['prediction'] = [ML_potential(down_3D['config'][i], data) for i in range(len(down_3D))]
    down_3D['difference'] = down_3D['energy'] - down_3D['prediction']
    return down_3D

#PLOT FUNCTIONS

def plot_symmetry(data, plot_position=624):
    symmetry_train = fig.add_subplot(plot_position)
    data.loc[data['train']==False].plot(
            kind='scatter',
            x='dO1O2',
            y='dHO1',
            s=20,
            c=test_point_color,
            label='Test',
            ax=symmetry_train)

    data.loc[data['train']==True].plot(
            kind='scatter',
            x='dO1O2',
            y='dHO1',
            s=20,
            alpha=0.4,
            c=train_point_color,
            label='Train',
            ax=symmetry_train)

    symmetry_train.set_xlabel(r"$dO_1O_2 (Bohr)$", fontsize=fontsize_axes_labels)
    symmetry_train.set_ylabel(r"${dHO_1} (Bohr)$", fontsize=fontsize_axes_labels)
    symmetry_train.legend(fontsize=fontsize_legend)
    return symmetry_train

def plot_text(data, plot_position=621):
    #----------System Information-------
    text = fig.add_subplot(plot_position)
    error = model_error(data)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    system_text = (" \n"
            " SYSTEM INFORMATION \n"
            " "+data['metadata'][0]['system_name']+"\n"
            " "+data['metadata'][0]['energy_calculation_method']+"\n"
            " "+data['metadata'][0]['temperature_of_simulation']+"\n"
            ""
           )
    descriptor_text = (" \n"
            " DESCRIPTOR INFORMATION \n"
            " "+data['metadata'][1]['descriptor_type']+"\n"
            " "+data['metadata'][1]['permutation']+"\n"
            " scaling = "+str(data['metadata'][1]['scaling'])+"\n"
            " \n"
            ""
           )
    if data['metadata'][2]['sampling_method'] == 'Kmeans':
        sampling_text = (" \n"
                " SAMPLING INFORMATION \n"
                " #training points = "+str(data['metadata'][2]['number_of_training_points'])+"\n"
                " sampling method = "+data['metadata'][2]['sampling_method']+"\n"
                " number_of_clusters = "+str(data['metadata'][2]['number_of_clusters'])+"\n"
                ""
               )
    else:
        sampling_text = (" \n"
                " SAMPLING INFORMATION \n"
                " #training points = "+str(data['metadata'][2]['number_of_training_points'])+"\n"
                " sampling method = "+data['metadata'][2]['sampling_method']+"\n"
                ""
               )
    model_text = (" \n"
            " MODEL INFORMATION \n"
            " Alpha = "+str(data['metadata'][3]['best_model'].get_params()['alpha'])+"\n"
            " Kernel Width = "+str(data['metadata'][3]['best_model'].get_params()['kernel__length_scale'])+"\n"
            " MAE_train="+str(error['MAE_train'])+"\n"
            " MAE_test="+str(error['MAE_test'])+"\n"
            " MSE_train="+str(error['MSE_train'])+"\n"
            " MSE_test="+str(error['MSE_test'])+"\n"
            " R2_train="+str(error['R2_train'])+"\n"
            " R2_test="+str(error['R2_test'])+"\n"
            ""
           )
    text.text(0, 0,
            system_text,
            fontsize=fontsize_text_boxes,
            verticalalignment='bottom',
            horizontalalignment='left',
            bbox=props
    )
    text.text(0.5, 0,
            descriptor_text,
            fontsize=fontsize_text_boxes,
            verticalalignment='bottom',
            horizontalalignment='left',
            bbox=props
    )
    text.text(1, 0,
            sampling_text,
            fontsize=fontsize_text_boxes,
            verticalalignment='bottom',
            horizontalalignment='left',
            bbox=props
    )
    text.text(1.5,0,
            model_text,
            fontsize=fontsize_text_boxes,
            verticalalignment='bottom',
            horizontalalignment='left',
            bbox=props
    )
    text.set_yticklabels([])
    text.set_xticklabels([])
    text.axis('off')
    return text

def plot_symmetry_energy(data, plot_position = 623):
    symmetry_energy = fig.add_subplot(plot_position)

    data.plot(kind='scatter',
            x='dO1O2',
            y='dHO1',
            s=20,
            c='energy',
            alpha=0.5,
            cmap = plt.get_cmap("jet"),
            colorbar=True,
            ax=symmetry_energy)

    symmetry_energy.set_xlabel(r"$dO_1O_2 (Bohr)$", fontsize=fontsize_axes_labels)
    symmetry_energy.set_ylabel(r"${dHO_1} (Bohr)$", fontsize=fontsize_axes_labels)
    return symmetry_energy

def plot_symmetry_train(data, plot_position=624):
    symmetry_train = fig.add_subplot(plot_position)

    data.loc[data['train']==False].plot(
            kind='scatter',
            x='dO1O2',
            y='dHO1',
            s=20,
            c=test_point_color,
            label='Test',
            ax=symmetry_train)

    data.loc[data['train']==True].plot(
            kind='scatter',
            x='dO1O2',
            y='dHO1',
            s=20,
            alpha=0.4,
            c=train_point_color,
            label='Train',
            ax=symmetry_train)

    symmetry_train.set_xlabel(r"$dO_1O_2 (Bohr)$", fontsize=fontsize_axes_labels)
    symmetry_train.set_ylabel(r"${dHO_1} (Bohr)$", fontsize=fontsize_axes_labels)
    symmetry_train.legend(fontsize=fontsize_legend)
    return symmetry_train

def plot_ground_truth_vs_prediction(data, plot_position=625):
    ground_vs_energy = fig.add_subplot(plot_position)

    data.loc[data['train']==False].plot(
            kind='scatter',
            x='prediction',
            y='energy',
            color=test_point_color,
            label='Test',
            alpha=0.3,
            ax=ground_vs_energy)

    data.loc[data['train']==True].plot(
            kind='scatter',
            x='prediction',
            y='energy',
            color=train_point_color,
            label='Train',
            alpha = 0.2,
            ax=ground_vs_energy)

    ground_vs_energy.set_xlabel("Predicted Energy", fontsize=fontsize_axes_labels)
    ground_vs_energy.set_ylabel("Ground Truth Energy", fontsize=fontsize_axes_labels)
    ground_vs_energy.legend(fontsize=fontsize_legend)
    return ground_vs_energy

def plot_histogram(data, plot_position=626):
    hist = fig.add_subplot(plot_position)

    data.plot(
            kind='hist',
            x='index',
            bins=20,
            y='energy',
            color=test_point_color,
            label='All data',
            density=True,
            ax=hist)

    data.loc[data['train']==True].plot(
            kind='hist',
            x='index',
            bins=20,
            y='energy',
            color=train_point_color,
            label='Train',
            density=True,
            alpha= 1,
            histtype=u'step',
            ax=hist)

    hist.set_xlabel('Energy (Kcal/mol)', fontsize=fontsize_axes_labels)
    hist.set_ylabel('Counts', fontsize=fontsize_axes_labels)
    hist.legend(fontsize=fontsize_legend)
    return hist

def plot_trajectory(data, plot_position=627, begin=0, end=100):
    traj = fig.add_subplot(plot_position)

    data.loc[begin:end].plot(
            kind='line',
            x='index',
            y='energy',
            color='blue',
            label='Ground Truth',
            ax=traj)

    data.loc[begin:end].plot(
            x='index',
            y='prediction',
            color='red',
            linestyle='-',
            label='prediction',
            ax=traj)

    traj.set_xlabel("Time", fontsize=fontsize_axes_labels)
    traj.set_ylabel("Energy", fontsize=fontsize_axes_labels)
    traj.legend(fontsize=fontsize_legend)

    return traj

def plot_2D_downfold(data,down_2D, plot_position=628):
    down_plot_2D = fig.add_subplot(plot_position)

    down_2D.loc[down_2D['dO1O2'] < 5].plot(
            kind='line',
            x='dO1O2',
            y='energy',
            color='blue',
            label='Ground Truth',
            ax=down_plot_2D)

    down_2D.loc[down_2D['dO1O2'] < 5].plot(
            kind='line',
            x='dO1O2',
            y='prediction',
            color='red',
            label='Prediction',
            ax=down_plot_2D)

    down_plot_2D.set_xlabel("dOO", fontsize=fontsize_axes_labels)
    down_plot_2D.set_xlabel("dOO", fontsize=fontsize_axes_labels)
    down_plot_2D.set_ylabel("Energy", fontsize=fontsize_axes_labels)
    down_plot_2D.legend(fontsize=fontsize_legend)
    return down_plot_2D

def plot_3D_downfold_truth(data, down_3D, plot_position=629):
    down_plot_3D_truth = fig.add_subplot(plot_position)

    down_3D['index'] = down_3D.index

    down_3D.plot(kind='scatter',
            x='dO1O2',
            y='dHO1',
            c='energy',
            s=100,
            alpha=0.8,
            cmap = plt.get_cmap("jet"),
            colorbar=True,
            marker='s',
            clim=(0,5),
            ax=down_plot_3D_truth)


    down_plot_3D_truth.set_xlabel("dOO (Bohr)", fontsize=fontsize_axes_labels)
    down_plot_3D_truth.set_ylabel("dOH", fontsize=fontsize_axes_labels)
    return down_plot_3D_truth

def plot_3D_downfold(data, down_3D, plot_position=[6,2,10]):
    down_plot_3D = fig.add_subplot(plot_position[0], plot_position[1], plot_position[2])

    down_3D['index'] = down_3D.index

    max_truth = down_3D['energy'].max()
    min_truth = down_3D['energy'].min()

    down_3D.loc[(down_3D['prediction']<max_truth) & (down_3D['prediction']>min_truth)].plot(
            kind='scatter',
            x='dO1O2',
            y='dHO1',
            c='prediction',
            s=100,
            alpha=0.8,
            cmap = plt.get_cmap("jet"),
            colorbar=True,
            marker='s',
            ax=down_plot_3D)


    down_plot_3D.set_xlabel("dOO (Bohr)", fontsize=18)
    down_plot_3D.set_ylabel("dOH", fontsize=18)
    down_plot_3D.title.set_size(20)
    return down_plot_3D

def plot_SW(data, plot_position=[6,2,11]):
    d = np.linspace(0,8,100)
    OO = data['metadata'][1]['switching_OO'](d)
    OH = data['metadata'][1]['switching_OH'](d)
    HH = data['metadata'][1]['switching_HH'](d)
    OH_plus = data['metadata'][1]['switching_OH_plus'](d)
    HH_plus = data['metadata'][1]['switching_HH_plus'](d)
    SW_df = pd.DataFrame({'OO':OO,
                         'OH':OH,
                         'HH':HH,
                         'OH_plus':OH_plus,
                         'HH_plus':HH_plus,
                         'distance': d})
    SW_plot = fig.add_subplot(plot_position[0], plot_position[1], plot_position[2])

    SW_df.plot(
            kind='line',
            x='distance',
            y='OO',
            color='blue',
            label='OO',
            ax=SW_plot)

    SW_df.plot(
            kind='line',
            x='distance',
            y='OH',
            color='orange',
            label='OH',
            ax=SW_plot)

    SW_df.plot(
            kind='line',
            x='distance',
            y='HH',
            color='red',
            label='HH',
            ax=SW_plot)

    SW_df.plot(
            kind='line',
            x='distance',
            y='HH_plus',
            color='black',
            linestyle='--',
            label='HH+',
            ax=SW_plot)

    SW_df.plot(
            kind='line',
            x='distance',
            y='OH_plus',
            color='green',
            label='OH+',
            ax=SW_plot)

    SW_plot.set_xlabel("Distance (Bohr)", fontsize=fontsize_axes_labels)
    SW_plot.set_ylabel("Switching Function", fontsize=fontsize_axes_labels)
    SW_plot.legend(fontsize=fontsize_legend)
    return SW_plot

def save_plot(data, save_folder = '../results/'):
    plot_save_folder =  save_folder
    plot_save_name = data['metadata'][1]['descriptor_type']+"_"+ \
                     data['metadata'][0]['temperature_of_simulation'] +"_"+ \
                     str(data['metadata'][2]['number_of_training_points'])+"_"+ \
                     data['metadata'][2]['sampling_method']
    plot_save_path = os.path.join(plot_save_folder, plot_save_name) +'.png'

    if os.path.isfile(plot_save_path) == False:
        print('SAVING PLOT TO: '+plot_save_path)
        plt.savefig(plot_save_path, bbox_inches='tight')
    else:
        print(plot_save_path + ' ALREADY EXISTS!')

#EXECUTE CODE ONE TIME

#PATHS FOR IMPORT AND EXPORT
model_path = '../2_model_training/trained_models/Coulomb_matrix_50K_1000_random'

#CALSULATIONS
data = import_model(model_path)
data = add_distances(data)

#PLOTS
fontsize_axes_labels = 20
fontsize_legend = 20
fontsize_text_boxes = 20
train_point_color = 'lightsalmon'
test_point_color = 'cornflowerblue'
begin=0
end=200
fig = plt.figure(figsize=(30,45))
plot_text(data)
#plot_symmetry(data)
plot_symmetry_energy(data)
plot_symmetry_train(data)
plot_ground_truth_vs_prediction(data)
plot_histogram(data)
plot_trajectory(data, begin=begin, end=end)
if data['metadata'][1]['descriptor_type'] == 'PIV':
    plot_SW(data)
#save_plot(data, save_folder)
fig.show()
