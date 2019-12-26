#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:55:21 2019

@author: julienh
"""

import numpy as np
import keras


def get_prediction(model,x_data):
        predictions = model.predict([np.stack(x_data.str[0].as_matrix()),
               np.stack(x_data.str[1].as_matrix()),
               np.stack(x_data.str[2].as_matrix()),
               np.stack(x_data.str[3].as_matrix()),
               np.stack(x_data.str[4].as_matrix()),
               np.stack(x_data.str[5].as_matrix()),
               np.stack(x_data.str[6].as_matrix()),
               ])
        return predictions

def energy_predictor_dropout(data,metadata):
    act_func,loss_func, optimizer, epochs, N_nodes_O, N_nodes_H, particles, path_to_output,date,import_from,patience = [metadata[x] for x in ['activation_function','loss_function','optimizer','epochs', 'N_nodes_O', 'N_nodes_H', 'particles', 'path_to_output','datetime','import_from','patience']]
        
    N_part, N_features = data["descriptor"][0].shape
    
    # Masks
    x_train = (data.loc[data['is_train']==True])['descriptor']
    y_train = (data.loc[data['is_train']==True])['energy']
    x_test = (data.loc[data['is_train']==False])['descriptor']
    y_test = (data.loc[data['is_train']==False])['energy']


# Creating layers

           
    if type(import_from) == type(None):     
                
        H_subnet = keras.Sequential(name='H_subnet')
        H_subnet.add(keras.layers.Dropout(0.2))
        for node in N_nodes_H:
            H_subnet.add(keras.layers.Dense(node,activation=act_func,kernel_constraint=keras.constraints.maxnorm(3)))
            H_subnet.add(keras.layers.Dropout(0.5))

        H_subnet.add(keras.layers.Dense(1,activation='linear'))

        O_subnet = keras.Sequential(name='O_subnet')
        O_subnet.add(keras.layers.Dropout(0.2))
        
        for node in N_nodes_O:
            O_subnet.add(keras.layers.Dense(node,activation=act_func,kernel_constraint=keras.constraints.maxnorm(3)))
            O_subnet.add(keras.layers.Dropout(0.5))
        O_subnet.add(keras.layers.Dense(1,activation='linear'))

        
        all_input_layers = []
        all_subnet = []
        count_H = 1
        count_O = 1
    
        for i_part in range(N_part):
            if particles[i_part] == 8:
                all_input_layers.append(keras.layers.Input(shape=(N_features,),name="O_input"+str(count_O)))
                count_O+=1
                all_subnet.append(O_subnet(all_input_layers[i_part]))
            elif particles[i_part] == 1:
                all_input_layers.append(keras.layers.Input(shape=(N_features,),name="H_input"+str(count_H)))
                all_subnet.append(H_subnet(all_input_layers[i_part]))

                count_H+=1
                
        #Add the subnets together
        added_layer = keras.layers.Add(name='Addition')
        added =added_layer(all_subnet)
        
        #Just a small modification, helps to get better fits
        out=added
    
        model = keras.models.Model(inputs=all_input_layers ,outputs=out)                                 #CHANGE ADDED
        model.compile(loss=loss_func, optimizer=optimizer, metrics=['accuracy'])
        
    else:
        
        model = keras.models.load_model(path_to_output+import_from+'/model')
        
        
        
    keras.utils.plot_model(model,to_file=path_to_output+date+'/network.png')


    es = keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience,restore_best_weights=True)


    # fit model
    history = model.fit([np.stack(x_train.str[0].as_matrix()),
               np.stack(x_train.str[1].as_matrix()),
               np.stack(x_train.str[2].as_matrix()),
               np.stack(x_train.str[3].as_matrix()),
               np.stack(x_train.str[4].as_matrix()),
               np.stack(x_train.str[5].as_matrix()),
               np.stack(x_train.str[6].as_matrix()),
               ], 
                y_train,
                validation_data=([np.stack(x_test.str[0].as_matrix()),
                   np.stack(x_test.str[1].as_matrix()),
                   np.stack(x_test.str[2].as_matrix()),
                   np.stack(x_test.str[3].as_matrix()),
                   np.stack(x_test.str[4].as_matrix()),
                   np.stack(x_test.str[5].as_matrix()),
                   np.stack(x_test.str[6].as_matrix()),
                   ],y_test),
                epochs=epochs,verbose=1,callbacks=[es])                                    # CHANGE EPOCHS
    
    mean_error = model.evaluate([np.stack(x_test.str[0].as_matrix()),
               np.stack(x_test.str[1].as_matrix()),
               np.stack(x_test.str[2].as_matrix()),
               np.stack(x_test.str[3].as_matrix()),
               np.stack(x_test.str[4].as_matrix()),
               np.stack(x_test.str[5].as_matrix()),
               np.stack(x_test.str[6].as_matrix()),
               ],y_test)[0]
    
    data['prediction'] = get_prediction(model, data['descriptor'])
    metadata['history'] = history
    model.save(path_to_output+date+'/model')        
    return data, model, mean_error



def energy_predictor(data,metadata):
    act_func,loss_func, optimizer, epochs, N_nodes_O, N_nodes_H, particles, path_to_output,date,import_from,patience = [metadata[x] for x in ['activation_function','loss_function','optimizer','epochs', 'N_nodes_O', 'N_nodes_H', 'particles', 'path_to_output','datetime','import_from','patience']]
        
    N_part, N_features = data["descriptor"][0].shape
    
    x_train = (data.loc[data['is_train']==True])['descriptor']
    y_train = (data.loc[data['is_train']==True])['energy']
    x_test = (data.loc[data['is_train']==False])['descriptor']
    y_test = (data.loc[data['is_train']==False])['energy']


# Creating layers

           
    if type(import_from) == type(None):     
                
        H_subnet = keras.Sequential(name='H_subnet')
        for node in N_nodes_H:
            H_subnet.add(keras.layers.Dense(node,activation=act_func))
        H_subnet.add(keras.layers.Dense(1,activation='linear'))

        O_subnet = keras.Sequential(name='O_subnet')
        for node in N_nodes_O:
            O_subnet.add(keras.layers.Dense(node,activation=act_func))
        O_subnet.add(keras.layers.Dense(1,activation='linear'))

        
        all_input_layers = []
        all_subnet = []
        count_H = 1
        count_O = 1
    
        for i_part in range(N_part):
            if particles[i_part] == 8:
                all_input_layers.append(keras.layers.Input(shape=(N_features,),name="O_input"+str(count_O)))
                count_O+=1
                all_subnet.append(O_subnet(all_input_layers[i_part]))
            elif particles[i_part] == 1:
                all_input_layers.append(keras.layers.Input(shape=(N_features,),name="H_input"+str(count_H)))
                all_subnet.append(H_subnet(all_input_layers[i_part]))

                count_H+=1
                
        #Add the subnets together
        added_layer = keras.layers.Add(name='Addition')
        added =added_layer(all_subnet)
        
        #Just a small modification, helps to get better fits
        out=added
    
        model = keras.models.Model(inputs=all_input_layers ,outputs=out)                                 #CHANGE ADDED
        model.compile(loss=loss_func, optimizer=optimizer, metrics=['accuracy'])
        
    else:
        
        model = keras.models.load_model(path_to_output+import_from+'/model')
        
        
        
    keras.utils.plot_model(model,to_file=path_to_output+date+'/network.png')


    es = keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience,restore_best_weights=True)


    # fit model
    history = model.fit([np.stack(x_train.str[0].as_matrix()),
               np.stack(x_train.str[1].as_matrix()),
               np.stack(x_train.str[2].as_matrix()),
               np.stack(x_train.str[3].as_matrix()),
               np.stack(x_train.str[4].as_matrix()),
               np.stack(x_train.str[5].as_matrix()),
               np.stack(x_train.str[6].as_matrix()),
               ], 
                y_train,
                validation_data=([np.stack(x_test.str[0].as_matrix()),
                   np.stack(x_test.str[1].as_matrix()),
                   np.stack(x_test.str[2].as_matrix()),
                   np.stack(x_test.str[3].as_matrix()),
                   np.stack(x_test.str[4].as_matrix()),
                   np.stack(x_test.str[5].as_matrix()),
                   np.stack(x_test.str[6].as_matrix()),
                   ],y_test),
                epochs=epochs,verbose=1,callbacks=[es])                                    # CHANGE EPOCHS
    
    mean_error = model.evaluate([np.stack(x_test.str[0].as_matrix()),
               np.stack(x_test.str[1].as_matrix()),
               np.stack(x_test.str[2].as_matrix()),
               np.stack(x_test.str[3].as_matrix()),
               np.stack(x_test.str[4].as_matrix()),
               np.stack(x_test.str[5].as_matrix()),
               np.stack(x_test.str[6].as_matrix()),
               ],y_test)[0]
    
    data['prediction'] = get_prediction(model, data['descriptor'])
    metadata['history'] = history
    model.save(path_to_output+date+'/model')        
    return data, model, mean_error





def energy_predictor_outdatedd(data,metadata):
    act_func,loss_func, optimizer, epochs, N_nodes_O, N_nodes_H, N_nodes_H_middle, particles, path_to_output,date,import_from = [metadata[x] for x in ['activation_function','loss_function','optimizer','epochs', 'N_nodes_O', 'N_nodes_H', 'N_nodes_H_middle', 'particles', 'path_to_output','datetime','import_from']]
        
    N_part, N_features = data["descriptor"][0].shape
    
    x_train = (data.loc[data['is_train']==True])['descriptor']
    y_train = (data.loc[data['is_train']==True])['energy']
    x_test = (data.loc[data['is_train']==False])['descriptor']
    y_test = (data.loc[data['is_train']==False])['energy']
    
    
    if type(import_from) == type(None):
        
        #Should H_middle have a different subnet
        if type(N_nodes_H_middle) == type(None):
            different_H_middle = False
        else:
            different_H_middle = True
        #Creating the subnets
        all_input_layers = []
        all_subnet = []
        count_H = 1
        count_O = 1
        
        for i_part in range(N_part):
            if particles[i_part] == 8:
                all_input_layers.append(keras.layers.Input(shape=(N_features,),name="O_input"+str(count_O)))
                subnet = keras.Sequential(name='O_subnet_'+str(count_O)+'__'+'-'.join(str(e) for e in N_nodes_O))
                count_O+=1
                for node in N_nodes_O:
                    subnet.add(keras.layers.Dense(node,activation=act_func))
                    
            if particles[i_part] == 1:
                if different_H_middle:
                    all_input_layers.append(keras.layers.Input(shape=(N_features,),name="H_middle_input"))
                    subnet = keras.Sequential(name='H_middle_subnet__'+'-'.join(str(e) for e in N_nodes_H_middle))
    
                    for node in N_nodes_H_middle:
                        subnet.add(keras.layers.Dense(node,activation=act_func))
                    different_H_middle = False
                else:
                    all_input_layers.append(keras.layers.Input(shape=(N_features,),name="H_input"+str(count_H)))
                    subnet = keras.Sequential(name='H_subnet_'+str(count_H)+"__"+'-'.join(str(e) for e in N_nodes_H))
                    count_H+=1
    
                    for node in N_nodes_H:
                        subnet.add(keras.layers.Dense(node,activation=act_func))
                        
            subnet.add(keras.layers.Dense(1,activation='linear'))
    
            all_subnet.append(subnet(all_input_layers[i_part]))
                
        #Add the subnets together
        added_layer = keras.layers.Add(name='Addition')
        added =added_layer(all_subnet)
        

        out=added
    
        model = keras.models.Model(inputs=all_input_layers ,outputs=out)                                 #CHANGE ADDED
        model.compile(loss=loss_func, optimizer=optimizer, metrics=['accuracy'])
        
    else:
        
        model = keras.models.load_model(path_to_output+import_from+'/model')
        
        
        
    keras.utils.plot_model(model,to_file=path_to_output+date+'/network.png')

    model.fit([np.stack(x_train.str[0].as_matrix()),
               np.stack(x_train.str[1].as_matrix()),
               np.stack(x_train.str[2].as_matrix()),
               np.stack(x_train.str[3].as_matrix()),
               np.stack(x_train.str[4].as_matrix()),
               np.stack(x_train.str[5].as_matrix()),
               np.stack(x_train.str[6].as_matrix()),
               ], y_train,epochs=epochs,verbose=1)                                    # CHANGE EPOCHS
    
    mean_error = model.evaluate([np.stack(x_test.str[0].as_matrix()),
               np.stack(x_test.str[1].as_matrix()),
               np.stack(x_test.str[2].as_matrix()),
               np.stack(x_test.str[3].as_matrix()),
               np.stack(x_test.str[4].as_matrix()),
               np.stack(x_test.str[5].as_matrix()),
               np.stack(x_test.str[6].as_matrix()),
               ],y_test)[0]
    
    data['prediction'] = get_prediction(model, data['descriptor'])
    model.save(path_to_output+date+'/model')        
    return data, model, mean_error







def energy_predictor_outdated(data,
                     N_nodes = 60, act_func= 'tanh',linear_adjust_end=True, 
                     loss_func = 'mean_absolute_error', optimizer='adam', epochs=100):
    
    N_part, nb_features = data["descriptor"][0].shape

    x_train = (data.loc[data['is_train']==True])['descriptor']
    y_train = (data.loc[data['is_train']==True])['energy']
    x_test = (data.loc[data['is_train']==False])['descriptor']
    y_test = (data.loc[data['is_train']==False])['energy']
    
        
    input_layers = np.empty(N_part,dtype=object)
    for i_layer in range(N_part):
        input_layers[i_layer] = keras.layers.Input(shape=(nb_features,))
        
    H_layers = np.empty(3,dtype=object)
    H_layers[0] = keras.layers.Dense(N_nodes,activation=act_func)                              # CHANGE LAYERS 
    H_layers[1] = keras.layers.Dense(N_nodes,activation=act_func)
    H_layers[2] = keras.layers.Dense(1,activation='linear')
    """
    H__middle_layers = np.empty(3,dtype=object)
    H__middle_layers[0] = keras.layers.Dense(N_nodes,activation='tanh')                              # CHANGE LAYERS 
    H__middle_layers[1] = keras.layers.Dense(N_nodes,activation='tanh')
    H__middle_layers[2] = keras.layers.Dense(1,activation='linear')
    """
    O_layers = np.empty(3,dtype=object)
    O_layers[0] = keras.layers.Dense(N_nodes,activation=act_func)   
    O_layers[1] = keras.layers.Dense(N_nodes,activation=act_func)
    O_layers[2] = keras.layers.Dense(1,activation='linear')
    
    added_layer = keras.layers.Add()
    added =added_layer([    O_layers[2](O_layers[1](O_layers[0](input_layers[0]))),
                            O_layers[2](O_layers[1](O_layers[0](input_layers[1]))),  
                            H_layers[2](H_layers[1](H_layers[0](input_layers[2]))), #       CHANGE MIDDLE LAYER
                            H_layers[2](H_layers[1](H_layers[0](input_layers[3]))),
                            H_layers[2](H_layers[1](H_layers[0](input_layers[4]))),
                            H_layers[2](H_layers[1](H_layers[0](input_layers[5]))),
                            H_layers[2](H_layers[1](H_layers[0](input_layers[6])))])
    
    if linear_adjust_end:
        out_layer = keras.layers.Dense(1,activation="linear")
        out = out_layer(added)
    else:
        out=added

    model = keras.models.Model(inputs=[input_layers[0],
                                       input_layers[1],
                                       input_layers[2],
                                       input_layers[3],
                                       input_layers[4],
                                       input_layers[5],
                                       input_layers[6]]
                                    ,outputs=out)                                 #CHANGE ADDED
    model.compile(loss=loss_func, optimizer=optimizer, metrics=['accuracy'])
    
    
    model.fit([np.stack(x_train.str[0].as_matrix()),
               np.stack(x_train.str[1].as_matrix()),
               np.stack(x_train.str[2].as_matrix()),
               np.stack(x_train.str[3].as_matrix()),
               np.stack(x_train.str[4].as_matrix()),
               np.stack(x_train.str[5].as_matrix()),
               np.stack(x_train.str[6].as_matrix()),
               ], y_train,epochs=epochs,verbose=1)                                    # CHANGE EPOCHS
    
    mean_error = model.evaluate([np.stack(x_test.str[0].as_matrix()),
               np.stack(x_test.str[1].as_matrix()),
               np.stack(x_test.str[2].as_matrix()),
               np.stack(x_test.str[3].as_matrix()),
               np.stack(x_test.str[4].as_matrix()),
               np.stack(x_test.str[5].as_matrix()),
               np.stack(x_test.str[6].as_matrix()),
               ],y_test)[0]
    
    data['prediction'] = get_prediction(model, data['descriptor'])
        
    
    
    return data, model, mean_error



