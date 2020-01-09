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
default_dropout_coef=np.zeros((default_n_species,default_n_hidden_layer+1)) # Dropout for faster convergence (can be desactivated) 
default_restore_weights=True
default_replace_inputs=False
default_plot_network=True
default_path_plot_network="./"
default_saved_model=False

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
            all_subnets.append( specie_subnets[specie](all_input_layers[atom]) )
            count_+=1
    added_layer = keras.layers.Add(name="Addition")
    #=========================================================================#

    return keras.models.Model(inputs=all_input_layers ,outputs=added_layer(all_subnets) )   

def predict(model, input_data ): 
        return model.predict(input_data)

def train(model, input_train, output_train, input_test, output_test, metadata,
                            n_epochs = default_n_epochs,
                            patience = default_patience,
                            restore_weights = default_restore_weights,
                            saved_model = default_saved_model
                            ):
    if not "n_epochs" in metadata:
        metadata["n_epochs"] = n_epochs
    if not "patience" in metadata :
        metadata["patience"] = patience
    if not "restore_weights" in metadata :
        metadata["restore_weights"] = restore_weights
    
    callback_ = keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=metadata["patience"] ,restore_best_weights=metadata["restore_weights"])

    # fit model
    metadata['history'] = model.fit( input_train, output_train, validation_data=(input_test,output_test), epochs=metadata["n_epochs"],verbose=1,callbacks=[callback_])  # CHANGE EPOCHS
    
    mean_error = model.evaluate(input_test,output_test)[0]

    if metadata["saved_model"] :
        model.save(metadata["path_saved_model"])
        
    return input_train, output_train, input_test, output_test, mean_error, metadata
 
 

