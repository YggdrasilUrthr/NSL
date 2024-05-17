import numpy as np
import tensorflow as tf                                            

from tensorflow import keras

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras import backend as K
from tensorflow.keras.utils import get_custom_objects

class XYTrigModel:

    def __init__(self, N_Epochs, N_Train, N_Valid, Sigma, N_Layers, Shape: list, optimizer, loss, activation, verbose = True):
        self.N_Epochs = N_Epochs
        self.N_Train = N_Train
        self.N_Valid = N_Valid
        self.Sigma = Sigma
        self.N_Layers = N_Layers
        self.Shape = Shape
        self.optimizer = optimizer
        self.loss = loss
        self.activation = activation
        self.verbose = verbose     

        if(len(Shape) < N_Layers):

            while len(Shape) != N_Layers:
                Shape.append(1)

        self.GenerateTrainData()

    def F(self, x, y):                                            # Target function
        return np.sin(x ** 2 + y ** 2)

    def GenerateTrainData(self):

        # generate training inputs
        np.random.seed(0)  
        x_train = np.random.uniform(-3/2, 3/2, self.N_Train)
        y_train = np.random.uniform(-3/2, 3/2, self.N_Train)
        x_valid = np.random.uniform(-3/2, 3/2, 50)
        y_valid = np.random.uniform(-3/2, 3/2, 50)
        x_valid.sort()
        y_valid.sort

        self.train_data = np.asarray([x_train, y_train]).transpose()
        self.valid_data = np.asarray([x_valid, y_valid]).transpose()
        self.z_target = self.F(self.valid_data[:, 0], self.valid_data[:, 1]) # ideal (target) polynomial function

        self.z_train = np.random.normal(self.F(self.train_data[:, 0], self.train_data[:, 1]), self.Sigma)
        self.z_valid = np.random.normal(self.F(self.valid_data[:, 0], self.valid_data[:, 1]), self.Sigma)

    def TrainModel(self):

        self.model = tf.keras.Sequential()
        self.model.add(Dense(self.Shape[0], input_shape=(2,), activation = self.activation))              # Input layer

        for i in range(1, self.N_Layers):
            self.model.add(Dense(self.Shape[i], activation = self.activation))

        self.model.add(Dense(1, activation = self.activation))                                            # Output layer

        # compile the model choosing optimizer, loss and metrics objects
        self.model.compile(optimizer = self.optimizer, loss = self.loss, metrics=['mse'])

        # fit the model using training dataset
        # over N_Epochs epochs of 32 batch size each
        # report training progress against validation data
        self.history = self.model.fit(x=self.train_data, y=self.z_train, 
            batch_size=32, epochs=self.N_Epochs,
            shuffle=True,
            validation_data=(self.valid_data, self.z_valid), verbose = self.verbose)
        
    def GetInfo(self, full = False):

        # get a summary of our composed model
        if full:
            self.model.summary()
        # return weights and biases
        print(self.model.get_weights())

    def EvaluateAgainstValid(self):             # Get evaluation score against validation data
        
        # evaluate model
        score = self.model.evaluate(self.valid_data, self.z_valid, batch_size=32, verbose = self.verbose)
        return score
    
    def EvaluateAgainstTarget(self):            # Get evaluation score against target data
        
        # evaluate model with the exact curve
        score = self.model.evaluate(self.valid_data, self.z_target, batch_size=32, verbose = self.verbose)
        return score