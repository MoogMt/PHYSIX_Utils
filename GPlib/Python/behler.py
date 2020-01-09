#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 13:32:52 2020

@author: julienh with modification from moogmt
"""

import numpy as np
import keras

# Neural Net default parameters
default_n_species=1
default_activation_fct = 'tanh'  # Activation function in the dense hidden layers
default_loss_fct = 'mean_squared_error' # Loss function in the NN
default_optimizer = 'adam'                    # Choice of optimizers for training of the NN weights 
default_n_epochs = 1000                  # Number of epoch for optimization?
default_patience = 100                  # Patience for convergence
default_n_nodes_per_layer= 80           # Number of nodes per hidden layer
default_n_hidden_layer=2                # Number of hidden layers
default_n_nodes_structure=np.ones((default_n_species,default_n_hidden_layer))*default_n_nodes_per_layer # Structure of the NNs (overrides the two precedent ones)
default_dropout_coef=np.zeros((default_n_hidden_layer+1,default_n_species)) # Dropout for faster convergence (can be desactivated) 
default_replace_inputs=False
default_plot_network=True
default_path_plot_network="./"

def handleNNOption( input_label, default_value, metadata, replace ):
    if not input_label in metadata or replace :
        metadata[input_label] = default_value
    return metadata

def buildNetwork( metadata,
                  activation_fct=default_activation_fct,
                  loss_fct=default_loss_fct,
                  optimizer=default_optimizer,
                  n_epochs=default_n_epochs,
                  patience=default_patience,
                  n_nodes_per_layer=default_n_nodes_per_layer,
                  n_hidden_layer=default_n_hidden_layer,
                  n_nodes_structure=default_n_nodes_structure,
                  dropout_coef=default_dropout_coef,
                  replace_inputs=default_replace_inputs
                 ):
    
    #Neural Net metadata
    #=========================================================================#
    metadata=handleNNOption("activation_fct", activation_fct, metadata, replace_inputs )
    metadata=handleNNOption("loss_fct", loss_fct, metadata, replace_inputs )
    metadata=handleNNOption("optimizer", optimizer, metadata, replace_inputs )
    metadata=handleNNOption("n_epochs", n_epochs, metadata, replace_inputs )
    metadata=handleNNOption("patience", patience, metadata, replace_inputs )
    metadata=handleNNOption("n_nodes_structure", n_hidden_layer, metadata, replace_inputs )
    metadata=handleNNOption("dropout_coef", dropout_coef, metadata, replace_inputs )
    #=========================================================================#
    
    # Construction subnetwork
    #=========================================================================#
    specie_subnets=[]
    for specie in range(metadata["n_species"]):
        specie_subnets=np.append(specie_subnets,keras.Sequential(name=str(metadata["species"][specie]+"_subnet")))
        specie_subnets[specie].add(keras.layers.Dropout(metadata["dropout_coef"][specie,0]))
        layer=1
        for n_node in metadata["n_nodes_structure"][specie,:] :
            if n_node > 0:
                specie_subnets[specie].add(keras.layers.Dense(n_node,activation=metadata["activation_fct"],kernel_constraint=keras.constraints.maxnorm(3) ) )
                specie_subnets[specie].add(keras.layers.Dropout(metadata["dropout_coef"][specie,layer]))
                layer += 1
            else:
                break
        specie_subnets[specie].add(keras.layers.Dense(1,activation="linear"))
    #=========================================================================#
    
    # Linking subnets
    #=========================================================================#
    all_input_layers=[]
    all_subnets=[]
    for specie in range(metadata["n_species"]):
        count_=1
        for atom in range(metadata["start_species"][specie],metadata["start_species"][specie]+metadata["nb_element_species"][specie]):
            all_input_layers.append( keras.layers.Input(shape=(metadata["n_features"],),name=str(metadata["species"][specie]+"_input"+str(count_) )) )            
            all_subnets.append( specie_subnets[specie](atom) )
            count_+=1
    added_layer = keras.layers.Add(name="Addition")
    #=========================================================================#

    return keras.models.Model(inputs=all_input_layers),added_layer  

def predict(model, input_data ):
        predictions = model.predict([np.stack(input_data.str[0].as_matrix()),
               np.stack(input_data.str[1].as_matrix()),
               np.stack(input_data.str[2].as_matrix()),
               np.stack(input_data.str[3].as_matrix()),
               np.stack(input_data.str[4].as_matrix()),
               np.stack(input_data.str[5].as_matrix()),
               np.stack(input_data.str[6].as_matrix()),
               ])
        return predictions

def train(model, input_train, output_train, input_test, output_test, metadata ):

    callback_ = keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=metadata["patience"],restore_best_weights=True)

    # fit model
    metadata['history'] = model.fit([np.stack(input_train.str[0].as_matrix()),
               np.stack(input_train.str[1].as_matrix()),
               np.stack(input_train.str[2].as_matrix()),
               np.stack(input_train.str[3].as_matrix()),
               np.stack(input_train.str[4].as_matrix()),
               np.stack(input_train.str[5].as_matrix()),
               np.stack(input_train.str[6].as_matrix()),
               ], 
                output_train,
                validation_data=([np.stack(input_test.str[0].as_matrix()),
                   np.stack(input_test.str[1].as_matrix()),
                   np.stack(input_test.str[2].as_matrix()),
                   np.stack(input_test.str[3].as_matrix()),
                   np.stack(input_test.str[4].as_matrix()),
                   np.stack(input_test.str[5].as_matrix()),
                   np.stack(input_test.str[6].as_matrix()),
                   ],output_test),
                epochs=metadata["epochs"],verbose=1,callbacks=[callback_])  # CHANGE EPOCHS
    
    mean_error = model.evaluate([np.stack(output_test.str[0].as_matrix()),
               np.stack(input_test.str[1].as_matrix()),
               np.stack(input_test.str[2].as_matrix()),
               np.stack(input_test.str[3].as_matrix()),
               np.stack(input_test.str[4].as_matrix()),
               np.stack(input_test.str[5].as_matrix()),
               np.stack(input_test.str[6].as_matrix()),
               ],output_test)[0]

    if metadata["save_model"] :
        model.save(metadata["path_saved_model"])
        
    return input_train, output_train, input_test, output_test, mean_error, metadata
 
 

