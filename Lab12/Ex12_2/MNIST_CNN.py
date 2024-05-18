#import os
import tensorflow as tf
from tensorflow import keras
import os
import numpy as np
import matplotlib.pyplot as plt
seed=0
np.random.seed(seed) # fix random seed
tf.random.set_seed(seed)

from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPooling2D
from keras.optimizers import SGD, Adam, RMSprop, Adagrad, Adadelta, Adam, Adamax, Nadam
from keras import saving

class CNN:

    def __init__(self, batch_size, epochs, optimizer = SGD(), verbose = True):
        self.batch_size = batch_size
        self.epochs = epochs
        self.optimizer = optimizer
        self.verbose = verbose

        self.load_dataset()

    def load_dataset(self):

        # input image dimensions
        self.img_rows, self.img_cols = 28, 28 # number of pixels 
        # output
        self.num_classes = 10 # 10 digits

        # the data, split between train and test sets
        (X_train, Y_train), (X_test, Y_test) = mnist.load_data()
        print('X_train shape:', X_train.shape)
        print('Y_train shape:', Y_train.shape)

        # reshape data, depending on Keras backend
        if keras.backend.image_data_format() == 'channels_first':
            X_train = X_train.reshape(X_train.shape[0], 1, self.img_rows, self.img_cols)
            X_test = X_test.reshape(X_test.shape[0], 1, self.img_rows, self.img_cols)
            self.input_shape = (1, self.img_rows, self.img_cols)
        else:
            X_train = X_train.reshape(X_train.shape[0], self.img_rows, self.img_cols, 1)
            X_test = X_test.reshape(X_test.shape[0], self.img_rows, self.img_cols, 1)
            self.input_shape = (self.img_rows, self.img_cols, 1)
            
        print('X_train shape:', X_train.shape)
        print('Y_train shape:', Y_train.shape)
        print()
        print(X_train.shape[0], 'train samples')
        print(X_test.shape[0], 'test samples')

        # cast to floats
        X_train = X_train.astype('float32')
        X_test = X_test.astype('float32')

        # rescale data in interval [0,1]
        X_train /= 255
        X_test /= 255

        # convert class vectors to binary class matrices, e.g. for use with categorical_crossentropy
        Y_train = keras.utils.to_categorical(Y_train, self.num_classes)
        Y_test = keras.utils.to_categorical(Y_test, self.num_classes)

        self.X_train = X_train
        self.Y_train = Y_train
        self.X_test = X_test
        self.Y_test = Y_test

    def create_CNN(self):
        model = Sequential()                                                                            # instantiate model
        model.add(Conv2D(6, kernel_size=(5, 5),
                    activation='relu',
                    input_shape = self.input_shape))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        model.add(Conv2D(16, kernel_size = (5, 5), activation = 'relu'))
        model.add(MaxPooling2D(pool_size = (2, 2)))
        model.add(Flatten())
        model.add(Dense(256, activation = 'relu'))                                                      # add a dense all-to-all relu layer                                                        # add a dense all-to-all relu layer
        model.add(Dropout(0.5))                                                                         # apply dropout with rate 0.5
        model.add(Dense(self.num_classes, activation='softmax'))                                        # soft-max layer
        self.model = model
        print('Model architecture created successfully!')

    def compile_model(self):
        self.create_CNN()                                                                   # create the model
        self.model.compile(loss=keras.losses.categorical_crossentropy,                      # compile the model
                    optimizer=self.optimizer,
                    metrics=['acc'])
        
    def evaluate_CNN(self):

        score = self.model.evaluate(self.X_test, self.Y_test, verbose=self.verbose)         # evaluate model
        self.score = score
        # print performance
        print()
        print('Test loss:', score[0])
        print('Test accuracy:', score[1])

        return (self.score, self.history)

    def train_CNN(self):
        self.compile_model()                                                                # create the deep neural net
        self.history = self.model.fit(self.X_train, self.Y_train,                           # train CNN and store training info in history
                batch_size=self.batch_size,
                epochs=self.epochs,
                verbose=self.verbose,
                validation_data=(self.X_test, self.Y_test))
        
        return self.evaluate_CNN()

    def CNN_predict(self):
        self.predictions = self.model.predict(self.X_test)
        return self.predictions
    
    def save_model(self, filename_model = './CNN.keras', filename_history = './hist.npy'):
        self.model.save(filename_model)
        np.save(filename_history, self.history)

    def load_model(self, filename_model = './CNN.keras', filename_history = './hist.npy'):
        print("Loading pre-trained model...")
        self.model = keras.models.load_model(filename_model)
        self.history = np.load(filename_history, allow_pickle = True).item()