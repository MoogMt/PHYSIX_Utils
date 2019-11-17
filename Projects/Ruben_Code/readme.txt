Instructions on how to use the Ruben code.

Step 0: get some ab-initio data
       
save the trajectory in a xyz file
save the energies in a text file where the potential energy column is labeled with 'potential energy'.

Step 1: convert the ab-initio data to a pandas dataframe

use the data_preperation.ipynb file to convert the raw text data to a pandas object
(jupyter notebook is needed to open .ipynb files)

Step 2: train the model on the pandas object

adjust the parameters on the model_trainer.py file under the #parameters# section to select the descriptor, amount of training points, training model (only Kernel ridge is included), parameters for grid search, ... 

The output of this python script is a new dataframe that contains the trained model (accessible in the metadata of the dataframe) and all the other metadata, as well as the predicted energies for each point in the dataframe.

Step 3: open the dataframe containing the predicted energies in the data_explorator.ipynb file to visualise the results

Note that many functions are specifically coded for the Zundel ion.

The Ruben_code folder contains a data_set (in the 0_raw_data folder) that can be used to go trough the different steps



