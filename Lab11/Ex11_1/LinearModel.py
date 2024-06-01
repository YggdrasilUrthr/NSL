import numpy as np
import tensorflow as tf
from tensorflow import keras

seed=0
np.random.seed(seed) # fix random seed
tf.random.set_seed(seed)

from keras.models import Sequential
from keras.layers import Dense, Activation
from keras import backend as K
from keras.utils import get_custom_objects

class LinearModel:

    def __init__(self, N_Epochs, N_Train, Sigma, m, b, verbose = True):
        self.N_Epochs = N_Epochs
        self.N_Train = N_Train
        self.Sigma = Sigma
        self.m = m
        self.b = b
        self.verbose = verbose
        self.history = None

        self.GenerateTrainData()

    def GenerateTrainData(self):

        # generate training inputs
        np.random.seed(0)
        self.x_train = np.random.uniform(-1, 1, self.N_Train)
        self.x_valid = np.random.uniform(-1, 1, 50)
        self.x_valid.sort()
        self.y_target = self.m * self.x_valid + self.b # ideal (target) linear function

        self.y_train = np.random.normal(self.m * self.x_train + self.b, self.Sigma) # actual measures from which we want to guess regression parameters
        self.y_valid = np.random.normal(self.m * self.x_valid + self.b, self.Sigma)

    def TrainModel(self):

        self.model = tf.keras.Sequential()
        self.model.add(Dense(1, input_shape=(1,)))

        # compile the model choosing optimizer, loss and metrics objects
        self.model.compile(optimizer='sgd', loss='mse', metrics=['mse'])

        # fit the model using training dataset
        # over N_Epochs epochs of 32 batch size each
        # report training progress against validation data
        self.history = self.model.fit(x=self.x_train, y=self.y_train, 
            batch_size=32, epochs=self.N_Epochs,
            shuffle=True, # a good idea is to shuffle input before at each epoch
            validation_data=(self.x_valid, self.y_valid), verbose = self.verbose)

    def GetInfo(self, full = False):
        
        # get a summary of our composed model
        if full:
            self.model.summary()
        # return weights and biases
        print("m: {}, b: {}".format(self.model.get_weights()[0], self.model.get_weights()[1]))

    def EvaluateAgainstValid(self):             # Get evaluation score against validation data
        
        # evaluate model
        score = self.model.evaluate(self.x_valid, self.y_valid, batch_size=32, verbose = self.verbose)
        return score
    
    def EvaluateAgainstTarget(self):            # Get evaluation score against target data
        
        # evaluate model with the exact curve
        score = self.model.evaluate(self.x_valid, self.y_target, batch_size=32, verbose = self.verbose)
        return score
    
    def SaveModel(self, filename_model = './LM.keras', filename_history = './hist.npy'):
        self.model.save(filename_model)
        np.save(filename_history, self.history)

    def LoadModel(self, filename_model = './LM.keras', filename_history = './hist.npy'):
        print("Loading pre-trained model...")
        self.model = keras.models.load_model(filename_model)
        self.history = np.load(filename_history, allow_pickle = True).item()        