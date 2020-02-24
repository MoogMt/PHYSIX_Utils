#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 18:30:00 2019

@author: julienh
"""
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
import matplotlib.pyplot as plt
import descriptors
from predictors import get_prediction
from ase import Atoms
import os


#PLOTS


def make_report(data,metadata):
    particles, model,scaler, folder_to_downfolded,path_to_output_images,time_launched,history=[metadata[x] for x in ['particles','model','scaler','folder_to_downfolded','path_to_output','datetime','history']]

#----------------------------IMPORT----------------------------------------------------------    
#--------------------------EXTRA CODE-----------------------------------------------------
    #Error of the model
    def model_error():
        rounding = 4
        train = data.loc[data['is_train'] == True]
        test = data.loc[data['is_train'] == False]
        range_energy = data['energy'].max() - data['energy'].min()
        MAE_train = mean_absolute_error(train['energy'],train['prediction'])
        MAE_train_percent = MAE_train/range_energy*100
        MAE_test = mean_absolute_error(test['energy'],test['prediction'])
        MAE_test_percent = MAE_test/range_energy*100
        MSE_train = mean_squared_error(train['energy'],train['prediction'])
        MSE_test = mean_squared_error(test['energy'],test['prediction'])
        R2_train = r2_score(train['energy'],train['prediction'])
        R2_test = r2_score(test['energy'],test['prediction'])

        errors = {'MAE_train' : np.round(MAE_train, rounding),
                  'MAE_train_percent' : np.round(MAE_train_percent,rounding),
                  'MAE_test' : np.round(MAE_test, rounding),
                  'MAE_test_percent' : np.round(MAE_test_percent,rounding),
                  'MSE_train' : np.round(MSE_train, rounding),
                  'MSE_test' : np.round(MSE_test, rounding),
                 'R2_train' : np.round(R2_train, rounding),
                 'R2_test' : np.round(R2_test, rounding)}

        return errors
    error = model_error()

    #Add dOH and dOO for plot
    data['dO1O2'] = [data['molec'].iloc[n].get_all_distances()[0][1] for n in range(len(data))]
    data['dHO1'] = [data['molec'].iloc[n].get_all_distances()[0][2] for n in range(len(data))] #Faux quand switch H
    
    ## IMPORT 3D DOWNFOLD DATA
    down_3D = pd.read_csv(folder_to_downfolded+'downfolded_3D_pickle')
    down_3D['pos'] = [np.matrix(down_3D['position'].iloc[i]).reshape(7,3) for i in range(len(down_3D))]
    down_3D = down_3D.rename(columns={'dOO':'dO1O2', 'dOH':'dHO1'})
    down_3D = down_3D.drop(['position', 'success'], axis=1)
    molec3D = np.empty(down_3D.index.max()+1,dtype=object)
    for i_config in range(down_3D.index.max()+1):
        molec3D[i_config] = Atoms(numbers=particles, positions=down_3D['pos'][i_config]) 
    down_3D['molec']=molec3D
    down_3D, scaler = getattr(descriptors, 'create_data_'+metadata['descriptor_type'])(down_3D,metadata)

    # IMPORT 2D DOWNFOLD DATA
    down_2D = pd.read_csv(folder_to_downfolded+'downfolded_2D_pickle')
    down_2D['pos'] = [np.matrix(down_2D['position'].iloc[i]).reshape(7,3) for i in range(len(down_2D))]
    down_2D = down_2D.rename(columns={'dOO':'dO1O2'})
    down_2D = down_2D.drop(['position', 'success'], axis=1)
    molec2D = np.empty(down_2D.index.max()+1,dtype=object)
    for i_config in range(down_2D.index.max()+1):
        molec2D[i_config] = Atoms(numbers=particles, positions=down_2D['pos'][i_config]) 
    down_2D['molec']=molec2D
    down_2D, scaler = getattr(descriptors, 'create_data_'+metadata['descriptor_type'])(down_2D,metadata)

    # CREATE 2D AND 3D PREDICTION
    down_2D['prediction'] = get_prediction(model, down_2D['descriptor'])
    down_3D['prediction'] = get_prediction(model, down_3D['descriptor'])
    down_3D['prediction'][0] = down_3D['energy'].max()
    down_3D['prediction'][-1] = down_3D['energy'].min()

##############################
    if metadata['scaler'] == False:
        isscaled = "False"
    else:
        isscaled = 'True'
        
    best_epoch = np.argmin(metadata['history'].history['val_loss'])

    try:
        metadata['old_N_feature']
    except KeyError:
        old_N_feature = metadata['N_feature']
    else:
        old_N_feature = metadata['old_N_feature']

    
    
    
    
#----------------------------PLOTS----------------------------------------------------------
    
    #Plot options
    fontsize_axes_labels = 20
    fontsize_legend = 20
    fontsize_text_boxes = 20
    train_point_color = 'lightsalmon'
    test_point_color = 'cornflowerblue'
    begin=0
    end=100
    fig = plt.figure(figsize=(30,45))

    data['index'] = data.index

    #----------System Information-------
    text = fig.add_subplot(521)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    
    system_text = (" \n"
            " SYSTEM INFORMATION \n"
            " "+metadata['system_name']+"\n"
            " "+metadata['energy_calculation_method']+"\n"
            " "+str(metadata['temperature'])+'K'+"\n"
            " test data : "+str(metadata['test_data'])+"\n"

            ""
           )
    descriptor_text = (" \n"
            " DESCRIPTOR INFORMATION \n"
            " "+metadata['descriptor_type']+"\n"
            " "+str(metadata["N_feature"])+"\n"
            " nodes H"+str(metadata["N_nodes_H"])+"\n"
            " nodes O"+str(metadata["N_nodes_O"])+"\n"
            #" "+data['metadata'][1]['permutation']+"\n"
            #" scaling = "+str(data['metadata'][1]['scaling'])+"\n"
            " \n"
            ""
           )
    sampling_text = (" \n"
            " SAMPLING INFORMATION \n"
            " #training points = "+str(int(metadata['total_time']*(1-metadata["test_size"])))+"\n"
            " scaling : "+isscaled+"\n"
            " N_PCA = "+str(metadata['N_PCA'])+"\n"
            " old nb feature = "+str(old_N_feature)+"\n"
            
            ""
           )
    
    model_text = (" \n"
            " MODEL INFORMATION \n"
            " Loss function = "+metadata["loss_function"]+"\n"
            " MAE_train="+str(error['MAE_train'])+"\n"
            " MAE_train_percent="+str(error['MAE_train_percent'])+"% \n"
            " MAE_test="+str(error['MAE_test'])+"\n"
            " MAE_test_percent="+str(error['MAE_test_percent'])+"% \n"
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

    #----------Symmetry Energy-------------
    symmetry_energy = fig.add_subplot(523)

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

    #---------Symmetry Train----------------
    history_loss = fig.add_subplot(524)
    history_loss.plot(history.history['loss'],label='train',)
    history_loss.plot(history.history['val_loss'],label='test')
    
    history_loss.plot(best_epoch, history.history['loss'][best_epoch],'*',markersize=14)
    history_loss.plot(best_epoch, history.history['val_loss'][best_epoch],'*',markersize=14)
    
    history_loss.set_yscale('log')
    history_loss.legend(fontsize=fontsize_legend)
    history_loss.set_xlabel('Epochs', fontsize=fontsize_axes_labels)
    history_loss.set_ylabel('Error', fontsize=fontsize_axes_labels)
    """
    data.loc[data['is_train']==False].plot(
            kind='scatter', 
            x='dO1O2', 
            y='dHO1', 
            s=20,
            c=test_point_color,
            label='Test',
            ax=symmetry_train)

    data.loc[data['is_train']==True].plot(
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
    """
    #-----Ground truth vs prediction---------
    ground_vs_energy = fig.add_subplot(525)

    data.loc[data['is_train']==False].plot(
            kind='scatter',
            x='prediction',
            y='energy', 
            color=test_point_color,
            label='Test',
            alpha=0.3,
            ax=ground_vs_energy)

    data.loc[data['is_train']==True].plot(
            kind='scatter',
            x='prediction',
            y='energy', 
            color=train_point_color,
            label='Train',
            alpha = 0.2,
            ax=ground_vs_energy)
    
    ground_vs_energy.plot([data['energy'].min(),data['energy'].max()],[data['energy'].min(),data['energy'].max()],
                           color = 'black', label='True')    

    ground_vs_energy.set_xlabel("Predicted Energy", fontsize=fontsize_axes_labels)
    ground_vs_energy.set_ylabel("Ground Truth Energy", fontsize=fontsize_axes_labels)
    ground_vs_energy.legend(fontsize=fontsize_legend)

    #-------------Histogram----------------
    hist = fig.add_subplot(526)

    data.plot(
            kind='hist', 
            x='index', 
            bins=20,
            y='energy', 
            color=test_point_color,
            label='All data',
            density=True,
            ax=hist)

    data.loc[data['is_train']==True].plot(
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

    #-------Trajectory Plot
    traj = fig.add_subplot(527)

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

    #-----------2D DOWNFOLD PLOT--------

    down_plot_2D = fig.add_subplot(528)

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
    
    down_plot_2D.axvline(data["dO1O2"].min(),linestyle='--')
    down_plot_2D.axvline(data["dO1O2"].max(),linestyle="--")

    down_plot_2D.set_xlabel("dOO", fontsize=fontsize_axes_labels)
    down_plot_2D.set_xlabel("dOO", fontsize=fontsize_axes_labels)
    down_plot_2D.set_ylabel("Energy", fontsize=fontsize_axes_labels)
    down_plot_2D.legend(fontsize=fontsize_legend)

    #----------3D DOWNFOLDED TRUTH
    down_plot_3D_truth = fig.add_subplot(529)

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
            ax=down_plot_3D_truth)


    down_plot_3D_truth.set_xlabel("dOO (Bohr)", fontsize=fontsize_axes_labels)
    down_plot_3D_truth.set_ylabel("dOH", fontsize=fontsize_axes_labels)

    #-----------3D Downfold
    down_plot_3D = fig.add_subplot(5,2,10)

    down_3D['index'] = down_3D.index

    down_3D.plot(kind='scatter',
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

    #--------Show and save Plot-----------------------------
    """
    plot_save_name = metadata['descriptor']+"_"+ \
                     metadata['temperature'] +"_"+ \
                     str(metadata['tot_time'])
    plot_save_path = os.path.join(path_to_output_images, plot_save_name) +'.png'
    """
    plt.savefig(path_to_output_images+time_launched+'.png', bbox_inches='tight')
    
    fig.show()