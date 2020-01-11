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
default_n_nodes_per_layer= 80           # Number of nodes per hidden layer
default_n_hidden_layer=2                # Number of hidden layers
default_n_nodes_structure=np.ones((default_n_species,default_n_hidden_layer))*default_n_nodes_per_layer # Structure of the NNs (overrides the two precedent ones)
default_dropout_coef=np.zeros((default_n_species,default_n_hidden_layer+1)) # Dropout for faster convergence (can be desactivated) 
default_restore_weights=True
default_replace_inputs=False
default_plot_network=True
default_path_plot_network="./network_plot.png"
#==============================================================================

# HANDLING OPTIONS OF NN
#==============================================================================
def handleNNOption( input_label, default_value, metadata):
    if not input_label in metadata :
        metadata[input_label] = default_value
    return metadata
#==============================================================================

#===================
# BUILDING NETWORK
#==============================================================================
def buildNetwork( metadata,
                 # Optionnal Arguments
                  activation_fct=default_activation_fct,
                  n_nodes_per_layer=default_n_nodes_per_layer,
                  n_hidden_layer=default_n_hidden_layer,
                  n_nodes_structure=default_n_nodes_structure,
                  dropout_coef=default_dropout_coef,
                  replace_inputs=default_replace_inputs
                 ):
    
    #Neural Net metadata
    #=========================================================================#
    metadata=handleNNOption("activation_fct", activation_fct, metadata, replace_inputs )
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
                specie_subnets[specie].add( keras.layers.Dense(n_node,activation=metadata["activation_fct"], kernel_constraint=keras.constraints.maxnorm(3) ) )
                specie_subnets[specie].add( keras.layers.Dropout(metadata["dropout_coef"][specie,layer]) )
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
def predict(model, input_, output_): 
    # PREDICTION
    prediction_raw=model.predict(input_)
    nb_data=np.shape(prediction_raw)[0]
    prediction_treat=np.zeros(nb_data)
    for i in range(nb_data):
        prediction_treat[i] = prediction_raw[i][0]
    return prediction_treat
#==============================================================================

#=================
# TRAINING MODEL
#==============================================================================
default_batch_size=10
default_verbose_train=0
default_saved_model=False
default_n_epochs = 1000                  # Number of epoch for optimization?
default_patience = 100                  # Patience for convergence
def train(model, input_train, output_train, input_test, output_test, metadata,
          # OPTIONNAL ARGS
          n_epochs = default_n_epochs,
          batch_size=default_batch_size,
          patience = default_patience,
          restore_weights = default_restore_weights,
          saved_model = default_saved_model,
          verbose_train=default_verbose_train
          ):
        
    metadata=handleNNOption("n_epochs", default_n_epochs, metadata )
    metadata=handleNNOption("verbose_train", verbose_train, metadata )
    metadata=handleNNOption("patience", default_patience, metadata )
    metadata=handleNNOption("batch_size", batch_size, metadata )
    metadata=handleNNOption("restore_weights", default_restore_weights, metadata )
    metadata=handleNNOption("saved_model", default_saved_model, metadata )
    
    # Fit Parameters
    callback_ = keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=metadata["verbose_train"], patience=metadata["patience"] ,restore_best_weights=metadata["restore_weights"])
    
    # Actual Fitting
    metadata['history'] = model.fit( input_train, output_train, validation_data=(input_test,output_test), epochs=metadata["n_epochs"],verbose=1,callbacks=[callback_])  # CHANGE EPOCHS
    
    # Compute the mean error
    mean_error = model.evaluate(input_test,output_test)[0]

    # Optionnal save model
    if metadata["saved_model"] :
        model.save(metadata["path_saved_model"])
        
    return model, mean_error, metadata
#==============================================================================


from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score

# Computing statistical errors
#==============================================================================
default_rounding=4 # Numerical precision on the error
def computeErrors( output_train, output_test, prediction_train, prediction_test, metadata,
                  # Optionnal args
                  rounding=default_rounding):
    
    metadata_stat_errors={}
    
    # Energie range
    metadata_stat_errors["range_energy"] =  output_train.max()-output_train.min()
    # Mean Absolute Errors
    metadata_stat_errors["MAE_train"] = np.round( mean_absolute_error( output_train, prediction_train ), rounding )
    metadata_stat_errors["MAE_test"]  = np.round( mean_absolute_error( output_test,  prediction_test  ), rounding )
    # Mean Squared Errors
    metadata_stat_errors["MSE_train"] = np.round( mean_squared_error( output_train, prediction_train ), rounding )
    metadata_stat_errors["MSE_test"]  = np.round( mean_squared_error( output_test,  prediction_test  ), rounding )
    # R2 Errors
    metadata_stat_errors["R2_train"] = np.round( r2_score( output_train, prediction_train ), rounding )
    metadata_stat_errors["R2_test"]  = np.round( r2_score( output_test,  prediction_test  ), rounding )
    # Rouding
    metadata_stat_errors["rounding"] = rounding
    
    return metadata_stat_errors
#==============================================================================

# Writting Statistical errors to disk
#==============================================================================
default_write_from_scratch = True
def writeStatErrors( metadata, path_out_stat, 
                    #Optionnal
                    write_from_scratch=default_write_from_scratch # Choose False if you want to append to file
                    ):
    # Choosing writing mode 
    write_="w"
    if not write_from_scratch:
        write_="w+"
    # Writting data
    file_out=open(path_out_stat,write_)
    file_out.write(str(metadata["MAE_train"])+" ",)
    file_out.write(str(metadata["MAE_test"])+" ")
    file_out.write(str(np.round(metadata["MAE_train"]/metadata["range_energy"]*100, metadata["rounding"] ))+ " ",)
    file_out.write(str(np.round(metadata["MAE_test"] /metadata["range_energy"]*100, metadata["rounding"] ))+ " ")
    file_out.write(str(metadata["MSE_train"])+" ",)
    file_out.write(str(metadata["MSE_test"])+" ")
    file_out.write(str(metadata["R2_train"])+" ",)
    file_out.write(str(metadata["R2_test"])+"\n")
    file_out.close()
    return
#==============================================================================

# Plotting prediction against reality           
#==============================================================================
def writeComparativePrediction( file_path, output_, prediction_ ):
    file_out=open(file_path,"w")
    nb_data=np.shape(output_)[0]
    for i in range(nb_data):
        file_out.write( str(str(output_[i])+" "+str(prediction_[i])+"\n") )
    file_out.close()
    return 
#==============================================================================
    
# ALL IN ONE
#==============================================================================
default_path_folder_save="./"
default_suffix_write=""
default_optimizer = 'adam'                    # Choice of optimizers for training of the NN weights 
default_loss_fct = 'mean_squared_error' # Loss function in the NN
def buildTrainPredict(metadata,input_train,input_test,output_train,output_test, 
                           loss_fct=default_loss_fct,
                           optimizer=default_optimizer,
                           path_folder_save=default_path_folder_save,
                           suffix_write=default_suffix_write):
    
    metadata=handleNNOption( "optimizer", optimizer, metadata )
    metadata=handleNNOption( "loss_fct", loss_fct, metadata )
    metadata=handleNNOption( "path_folder_save", path_folder_save, metadata )
    metadata=handleNNOption( "suffix_write", suffix_write, metadata )
    
    # BUILDING
    #=============================================================================#
    model=buildNetwork(metadata)
    # Compile the network
    model.compile(loss=metadata["loss_fct"], optimizer=metadata["optimizer"], metrics=['mse'])
    # Plot the network
    if metadata["plot_network"]:
        keras.utils.plot_model(model,to_file=metadata["path_plot_network"],show_shapes=True, show_layer_names=True)
    #=============================================================================#

    # TRAINING NETWORK
    #=============================================================================#
    model, mean_error, metadata = train(model,input_train,output_train,input_test,output_test,metadata)
    #=============================================================================#

    # PREDICTION
    #===============================================================================
    # Make the prediction
    predictions_train = predict( model, input_train, output_train )
    predictions_test  = predict( model, input_test,  output_test  )
    # Compute Statistical errors
    metadata_stat = computeErrors( output_train, output_test, predictions_train, predictions_test, metadata )
    path_stat_err=str(metadata["path_folder_save"]+"StatisticalErr_"+metadata["suffix_write"])
    writeStatErrors( metadata_stat, path_stat_err)    
    #===============================================================================    
    return model, metadata, metadata_stat, predictions_train, predictions_test
    