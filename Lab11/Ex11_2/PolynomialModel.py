import numpy as np
import tensorflow as tf                                            

from tensorflow import keras

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras import backend as K
from tensorflow.keras.utils import get_custom_objects

class PolynomialModel:

    def __init__(self, N_Epochs, N_Train, Sigma, N_Layers, Shape: list, optimizer, loss, activation, verbose = True):
        self.N_Epochs = N_Epochs
        self.N_Train = N_Train
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

    def TP(self, x):                                            # Target Polynomial
        return 4 - 3 * x - 2 * (x ** 2) + 3 * (x ** 3) 

    def GenerateTrainData(self):

        # generate training inputs
        np.random.seed(0)  
        self.x_train = np.random.uniform(-1, 1, self.N_Train)
        self.x_valid = np.random.uniform(-1, 1, 50)
        self.x_valid.sort()
        self.y_target = self.TP(self.x_valid) # ideal (target) polynomial function

        self.y_train = np.random.normal(self.TP(self.x_train), self.Sigma) # actual measures from which we want to guess regression parameters
        self.y_valid = np.random.normal(self.TP(self.x_valid), self.Sigma)

    def TrainModel(self):

        self.model = tf.keras.Sequential()
        self.model.add(Dense(self.Shape[0], input_shape=(1,), activation = self.activation))              # Input layer

        for i in range(1, self.N_Layers):
            self.model.add(Dense(self.Shape[i], activation = self.activation))

        self.model.add(Dense(1, activation = self.activation))                                            # Output layer

        # compile the model choosing optimizer, loss and metrics objects
        self.model.compile(optimizer = self.optimizer, loss = self.loss, metrics=['mse'])

        # fit the model using training dataset
        # over N_Epochs epochs of 32 batch size each
        # report training progress against validation data
        self.history = self.model.fit(x=self.x_train, y=self.y_train, 
            batch_size=32, epochs=self.N_Epochs,
            shuffle=True,
            validation_data=(self.x_valid, self.y_valid), verbose = self.verbose)
        
    def GetInfo(self, full = False):

        # get a summary of our composed model
        if full:
            self.model.summary()
        # return weights and biases
        print(self.model.get_weights())

    def EvaluateAgainstValid(self):             # Get evaluation score against validation data
        
        # evaluate model
        score = self.model.evaluate(self.x_valid, self.y_valid, batch_size=32, verbose = self.verbose)
        return score
    
    def EvaluateAgainstTarget(self):            # Get evaluation score against target data
        
        # evaluate model with the exact curve
        score = self.model.evaluate(self.x_valid, self.y_target, batch_size=32, verbose = self.verbose)
        return score