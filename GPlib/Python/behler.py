#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 13:32:52 2020

@author: julienh with modification from moogmt
"""

import numpy as np
import keras

#===================
# BUILDING NETWORK
default_n_species=1
default_activation_fct = 'tanh'  # Activation function in the dense hidden layers
default_n_nodes_per_layer= 80           # Number of nodes per hidden layer
default_n_hidden_layer=2                # Number of hidden layers
default_n_nodes_structure=np.ones((default_n_species,default_n_hidden_layer))*default_n_nodes_per_layer # Structure of the NNs (overrides the two precedent ones)
default_dropout_coef=np.zeros((default_n_species,default_n_hidden_layer+1)) # Dropout for faster convergence (can be desactivated) 
default_replace_inputs=False
default_plot_network=True
default_path_plot_network="./network_plot.png"
default_kernel_constraint=None
default_out_number=1
#=============================================================================
def buildSpecieSubNetwork( specie, n_nodes, drop_out_rate=default_dropout_coef, activation_fct=default_activation_fct, kernel_constraint=default_kernel_constraint, out_number=default_out_number ):
    subnet=keras.Sequential( name=str( specie+"_subnet" ) )
    subnet.add(keras.layers.Dropout( drop_out_rate[0] ))
    layer=1
    for node in n_nodes:
        if node > 0:
            subnet.add( keras.layers.Dense( node, activation=activation_fct, kernel_constraint=kernel_constraint ) )
            subnet.add( keras.layers.Dropout( drop_out_rate[ layer ] ) )
            layer += 1
        else:
            break
    subnet.add( keras.layers.Dense( out_number, activation="linear" ) )
    return subnet
#==============================================================================
def buildAllAtomsNetworks( species, nb_species, start_species, nb_element_species, n_features, specie_subnets ):
    all_layers   = []
    all_networks = []
    count_all=0
    for specie in range( nb_species ):
        count_ = 1
        for atom in range( start_species[specie], start_species[specie]+nb_element_species[specie] ):
            all_layers.append( keras.layers.Input( shape=(n_features,), name =str( species[specie] + "_input" + str(count_) ) ) )
            all_networks.append( specie_subnets[specie]( all_layers[ count_all ] ) )
            count_ += 1
            count_all += 1 
    return all_networks, all_layers
#==============================================================================
def buildNetwork( species, nb_species, n_features, start_species, nb_element_species, nodes_structure, drop_out_rate, activation_fct, kernel_constraint ):
    
    # Construction subnetwork
    #=========================================================================#
    specie_subnets=[]
    for specie in range( nb_species ):
        specie_subnets.append( buildSpecieSubNetwork( species[specie], nodes_structure[specie,:], drop_out_rate[specie,:], activation_fct, kernel_constraint, 1 ) )
    #=========================================================================#
    
    # Building All Networks
    #=========================================================================#
    all_networks, all_layers = buildAllAtomsNetworks( species, nb_species, start_species, nb_element_species, n_features, specie_subnets )
    #=========================================================================#

    # Final Addition Layer
    #==========================================================================
    final_layer = keras.layers.Add( name="Addition" )
    #==========================================================================

    return keras.models.Model( inputs=all_layers, outputs=final_layer( all_networks ) )   

#==============================================================================


#============
# PREDICTION
#==============================================================================
def predict(model, input_ ): 
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
default_verbose_train=1
default_n_epochs=1000
default_patience = 100                  # Patience for convergence
default_saved_model=False
default_path_saved_model="./saved_model"
default_batch_size=None
default_restore_weights=True
def train( model, input_train, output_train, input_test, output_test, 
          n_epochs=default_n_epochs, 
          batch_size=default_batch_size, 
          patience = default_patience, 
          restore_weights = default_restore_weights, 
          saved_model = default_saved_model, 
          path_saved_model = default_path_saved_model, 
          verbose_train=default_verbose_train ):
    
    # Fit Parameters
    callback_ = keras.callbacks.EarlyStopping( monitor='val_loss', mode='min', verbose=verbose_train, patience=patience ,restore_best_weights=restore_weights ) # Parameters for Early Stopping
    
    # Actual Fitting
    history = model.fit( input_train, output_train, validation_data=(input_test,output_test), epochs=n_epochs, verbose=verbose_train, callbacks=[callback_]) # Training pocedure
    
    # Optionnal save model
    if saved_model :
        model.save( path_saved_model )
        
    return model, history
#==============================================================================


from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score

# Computing statistical errors
#==============================================================================
default_rounding=4 # Numerical precision on the error
def computeErrors( output_train, output_test, prediction_train, prediction_test, 
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
def writeStatErrors( metadata, path_out_stat ):
    # Writting data
    file_out=open(path_out_stat,"w")
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
default_early_stop_metric=['mse']
default_plot_network = False
default_path_plot_network = "./plot_network.png"
default_suffix_write=""
def buildTrainPredictWrite(input_train,input_test,output_train,output_test, 
                           species, 
                           nb_species, 
                           n_features, 
                           start_species, 
                           nb_element_species, 
                           nodes_structure, 
                           drop_out_rate, 
                           n_epochs=default_n_epochs,
                           batch_size=default_batch_size,
                           kernel_constraint=default_kernel_constraint,
                           activation_fct = default_activation_fct,
                           loss_fct=default_loss_fct,
                           optimizer=default_optimizer,
                           path_folder_save=default_path_folder_save,
                           early_stop_metric=default_early_stop_metric,
                           plot_network=default_plot_network,
                           path_plot_network=default_path_plot_network,
                           suffix_write=default_suffix_write
                           ):
    
    # BUILDING
    #=============================================================================#
    model=buildNetwork( species, nb_species, n_features, start_species, nb_element_species, nodes_structure, drop_out_rate, activation_fct, kernel_constraint )
    # Compile the network
    model.compile(loss=loss_fct, optimizer=optimizer, metrics=early_stop_metric)
    # Plot the network
    if plot_network:
        keras.utils.plot_model(model,to_file=path_plot_network,show_shapes=True, show_layer_names=True)
    #=============================================================================#

    # TRAINING NETWORK
    #=============================================================================#
    model, history = train(model,input_train,output_train,input_test,output_test,n_epochs,batch_size,verbose_train=1)
    #=============================================================================#

    # PREDICTION
    #===============================================================================
    # Make the prediction
    predictions_train = predict( model, input_train )
    predictions_test  = predict( model, input_test  )
    # Compute Statistical errors
    metadata_stat = computeErrors( output_train, output_test, predictions_train, predictions_test )
    path_stat_err=str( path_folder_save + "StatisticalErr_" + suffix_write)
    writeStatErrors( metadata_stat, path_stat_err)    
    #===============================================================================    
    return model, metadata_stat, predictions_train, predictions_test

def getAtomicEnergy( specie, start_specie, nb_element_specie, n_nodes, drop_out_rate,  input_, model_general, activation_fct=default_activation_fct, loss_fct=default_loss_fct, optimizer=default_optimizer, kernel_constraint=default_kernel_constraint, early_stop_metric=default_early_stop_metric ):
    # Parameters
    nb_data=np.shape(input_)[1]
    n_features=np.shape(input_)[2]
    # Reshaping input
    input_reshape=np.array( input_[start_specie:start_specie+nb_element_specie]).reshape(nb_data*nb_element_specie,n_features)
    # Buildint structure
    network_structure = buildSpecieSubNetwork( specie, n_nodes, drop_out_rate, activation_fct, kernel_constraint, 1 )
    # Define Inputs
    input_layers = [ keras.layers.Input( shape=(n_features,), name =str( specie + "_input_energy" ) )  ]
    # Creating model
    model_specie=keras.models.Model( inputs=input_layers, outputs=network_structure( input_layers ) )   
    # Compiling model 
    model_specie.compile(loss=loss_fct, optimizer=optimizer, metrics=default_early_stop_metric)
    # Getting weights from previous model
    model_specie.set_weights ( model_general.get_layer( specie+"_subnet" ).get_weights() )
    # Prediction
    return predict( model_specie, input_reshape ) 


#==============================================================================
def buildNetwork_OLD( metadata,
                 # Optionnal Arguments
                  activation_fct=default_activation_fct,
                  n_nodes_per_layer=default_n_nodes_per_layer,
                  n_hidden_layer=default_n_hidden_layer,
                  n_nodes_structure=default_n_nodes_structure,
                  dropout_coef=default_dropout_coef,
                  replace_inputs=default_replace_inputs
                 ):
    
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