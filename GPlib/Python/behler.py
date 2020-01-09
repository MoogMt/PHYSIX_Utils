#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 13:32:52 2020

@author: julienh with modification from moogmt
"""

import numpy as np
import keras

# Neural Net Default Parameters
#==============================================================================
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
default_path_plot_network="./network_plot.png"
default_saved_model=False
#==============================================================================

# HANDLING OPTIONS OF NN
#==============================================================================
def handleNNOption( input_label, default_value, metadata, replace ):
    if not input_label in metadata or replace :
        metadata[input_label] = default_value
    return metadata
#==============================================================================

#===================
# BUILDING NETWORK
#==============================================================================
def buildNetwork( metadata,
                 # Optionnal Arguments
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
#==============================================================================


#============
# PREDICTION
#==============================================================================
default_save_prediction = False
default_save_prediction_path = "./prediction.dat"
def predict(model, metadata, input_, output_,
            # Optionnal
            save_prediction=default_save_prediction, 
            save_prediction_path=default_save_prediction_path 
            ): 
    # CHECK 
    if not "save_prediction" in metadata:
        metadata["save_prediction"] = save_prediction
    if not "save_prediction_path" in metadata:
        metadata["save_prediction_path"] = save_prediction_path
        
    # PREDICTION
    prediction=model.predict(input_)    
    
    # OPTIONNAL SAVE RESULTS
    if metadata["save_prediction"]:
        file_out=open(metadata["save_prediction_path"],"w")
        nb_data=np.shape(input_)[1]
        for i in range(nb_data):
            file_out.write(str(output_[i])+" "+str(prediction[i][0])+"\n")
        file_out.close()
    return prediction
#==============================================================================

#=================
# TRAINING MODEL
#==============================================================================
def train(model, input_train, output_train, input_test, output_test, metadata,
          # OPTIONNAL ARGS
          n_epochs = default_n_epochs,
          patience = default_patience,
          restore_weights = default_restore_weights,
          saved_model = default_saved_model
          ):
    
    # CHECK
    if not "n_epochs" in metadata:
        metadata["n_epochs"] = n_epochs
    if not "patience" in metadata :
        metadata["patience"] = patience
    if not "restore_weights" in metadata :
        metadata["restore_weights"] = restore_weights
    
    # Fit Parameters
    callback_ = keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=metadata["patience"] ,restore_best_weights=metadata["restore_weights"])

    # Actual Fitting
    metadata['history'] = model.fit( input_train, output_train, validation_data=(input_test,output_test), epochs=metadata["n_epochs"],verbose=1,callbacks=[callback_])  # CHANGE EPOCHS
    
    # Compute the mean error
    mean_error = model.evaluate(input_test,output_test)[0]

    # Optionnal save model
    if metadata["saved_model"] :
        model.save(metadata["path_saved_model"])
        
    return model, mean_error, metadata
#==============================================================================

# Plotting errors
#==============================================================================
    
#==============================================================================
