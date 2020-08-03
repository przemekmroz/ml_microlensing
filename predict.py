import tensorflow as tf
from tensorflow import keras

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('model')
parser.add_argument('features')
args = parser.parse_args()

#-----Load model from HDF5 file----#

model = keras.models.load_model(args.model)

#-----Load features and labels-----#

test_features = np.loadtxt(args.features)

#-----Calculate classifier's predictions and save to a text file */

preds = model.predict(test_features)

np.savetxt('predictions.txt',preds,fmt='%7.5f')

