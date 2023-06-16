import numpy as np

# Define the RNA alphabet
rna_alphabet = ['A', 'C', 'G', 'T']


def rna_to_one_hot(rna_sequence):
    """
    Converts an RNA sequence to a 1-hot encoded matrix.

    Parameters:
        rna_sequence (str): RNA sequence to be converted.

    Returns:
        numpy.ndarray: 1-hot encoded matrix of the RNA sequence.
    """
    # Initialize an empty matrix with the correct dimensions
    one_hot = np.zeros((len(rna_sequence), len(rna_alphabet)), dtype=int)

    # Iterate through each nucleotide in the RNA sequence
    for i, nucleotide in enumerate(rna_sequence):
        if nucleotide == "N":
            continue
        # Find the index of the nucleotide in the RNA alphabet
        j = rna_alphabet.index(nucleotide)
        # Set the corresponding element in the 1-hot matrix to 1
        one_hot[i, j] = 1

    return one_hot

# Read File
# change fasta file into a list of sequences
# fasta reader
def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    return lines

# fasta parser


def parse_fasta(lines):
    seq_list = []
    seq_name = [int(line[-2]) for line in lines if line.startswith('>')]
    seq = ''
    for line, index in zip(lines, range(len(lines))):
        if index == len(lines) - 1:
            seq += line.strip()
            seq_list.append(seq)
        if line.startswith('>'):
            seq_list.append(seq)
            seq = ''
            continue
        else:
            seq += line.strip()
    for i in seq_list:
        if i == '':
            seq_list.remove(i)
    return seq_list, seq_name

# Sim 1
x_train, y_train = parse_fasta(read_fasta('./sim1/train.fasta'))
x_train = [rna_to_one_hot(i) for i in x_train]
x_train = np.array(x_train)
y_train = np.array(y_train)

x_validation, y_validation = parse_fasta(read_fasta('./sim1/validation.fasta'))
x_validation = [rna_to_one_hot(i) for i in x_validation]
x_validation = np.array(x_validation)
y_validation = np.array(y_validation)

x_test, y_test = parse_fasta(read_fasta('./sim1/test.fasta'))
x_test = [rna_to_one_hot(i) for i in x_test]
x_test = np.array(x_test)
y_test = np.array(y_test)

from tensorflow import keras as keras
from keras.models import Sequential
from keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from keras.callbacks import EarlyStopping, History

# Define the model
sequence_length = x_train.shape[1]
model = Sequential()
model.add(Conv1D(filters=32, kernel_size=10, activation='relu',
          input_shape=(sequence_length, 4)))
model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(units=64, activation='LeakyReLU'))
model.add(Dense(units=1, activation='sigmoid'))

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy',
              metrics=['accuracy'])

# Train the model
model.fit(x_train, y_train, verbose=1,  
 validation_data=(x_validation, y_validation),batch_size=128, epochs=100,  
       callbacks=[EarlyStopping(patience=10, monitor="val_loss", 
 restore_best_weights=True), History()])  

# Evaluate the model
model.evaluate(x_test, y_test)

from sklearn.metrics import roc_curve
from sklearn.metrics import auc
y_pred = model.predict(x_test).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred)
auc_keras = auc(fpr_keras, tpr_keras)
auc_keras
import matplotlib.pyplot as plt
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
plt.show()

#Sim 2
from keras.layers import Conv1D, MaxPooling1D, Flatten, Dense
from keras.models import Sequential
from tensorflow import keras as keras
x_train, y_train = parse_fasta(read_fasta('./sim2/train.fasta'))
x_train = [rna_to_one_hot(i) for i in x_train]
x_train = np.array(x_train)
y_train = np.array(y_train)

x_validation, y_validation = parse_fasta(read_fasta('./sim2/validation.fasta'))
x_validation = [rna_to_one_hot(i) for i in x_validation]
x_validation = np.array(x_validation)
y_validation = np.array(y_validation)

x_test, y_test = parse_fasta(read_fasta('./sim2/test.fasta'))
x_test = [rna_to_one_hot(i) for i in x_test]
x_test = np.array(x_test)
y_test = np.array(y_test)


# Define the model
sequence_length = x_train.shape[1]
model = Sequential()
model.add(Conv1D(filters=32, kernel_size=20, activation='relu',
          input_shape=(sequence_length, 4)))
model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(units=64, activation='LeakyReLU'))
model.add(Dense(units=1, activation='sigmoid'))

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy',
              metrics=['accuracy'])

# Train the model
model.fit(x_train, y_train, verbose=1,  
 validation_data=(x_validation, y_validation),batch_size=128, epochs=100,  
       callbacks=[EarlyStopping(patience=10, monitor="val_loss", 
 restore_best_weights=True), History()])  

# Evaluate the model
model.evaluate(x_test, y_test)

from sklearn.metrics import roc_curve
from sklearn.metrics import auc
y_pred = model.predict(x_test).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred)
auc_keras = auc(fpr_keras, tpr_keras)
auc_keras

import matplotlib.pyplot as plt
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
plt.show()

# Sim 6
from keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout
from keras.models import Sequential
from tensorflow import keras as keras
x_train, y_train = parse_fasta(read_fasta('./sim6/train.fasta'))
x_train = [rna_to_one_hot(i) for i in x_train]
x_train = np.array(x_train)
y_train = np.array(y_train)

x_validation, y_validation = parse_fasta(read_fasta('./sim6/validation.fasta'))
x_validation = [rna_to_one_hot(i) for i in x_validation]
x_validation = np.array(x_validation)
y_validation = np.array(y_validation)

x_test, y_test = parse_fasta(read_fasta('./sim6/test.fasta'))
x_test = [rna_to_one_hot(i) for i in x_test]
x_test = np.array(x_test)
y_test = np.array(y_test)


# Define the model
sequence_length = x_train.shape[1]
model = Sequential()
model.add(Conv1D(filters=32, kernel_size=20, activation='relu',
          input_shape=(sequence_length, 4)))


model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(units=64, activation='LeakyReLU'))
model.add(Dropout(rate=0.5))
model.add(Dense(units=1, activation='sigmoid'))

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy',
              metrics=['accuracy'])

# Train the model
model.fit(x_train, y_train, verbose=1,  
 validation_data=(x_validation, y_validation),batch_size=128, epochs=100,  
       callbacks=[EarlyStopping(patience=10, monitor="val_loss", 
 restore_best_weights=True), History()])  

# Evaluate the model
model.evaluate(x_test, y_test)

from sklearn.metrics import roc_curve
from sklearn.metrics import auc
y_pred = model.predict(x_test).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test, y_pred)
auc_keras = auc(fpr_keras, tpr_keras)
auc_keras

import matplotlib.pyplot as plt
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
plt.show()

# Sim 7
from keras import regularizers
from keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout
from keras.models import Sequential
from tensorflow import keras as keras
x_train, y_train = parse_fasta(read_fasta('./sim7/train.fasta'))
x_train = [rna_to_one_hot(i) for i in x_train]
x_train = np.array(x_train)
y_train = np.array(y_train)

x_validation, y_validation = parse_fasta(read_fasta('./sim7/validation.fasta'))
x_validation = [rna_to_one_hot(i) for i in x_validation]
x_validation = np.array(x_validation)
y_validation = np.array(y_validation)

x_test, y_test = parse_fasta(read_fasta('./sim7/test.fasta'))
x_test = [rna_to_one_hot(i) for i in x_test]
x_test = np.array(x_test)
y_test = np.array(y_test)


# Define the model
sequence_length = x_train.shape[1]
model = Sequential()
model.add(Conv1D(filters=32, kernel_size=20, activation='relu',
          input_shape=(sequence_length, 4)))

model.add(MaxPooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(units=512, activation='LeakyReLU'))
model.add(Dropout(rate=0.5))
model.add(Dense(units=128, activation='LeakyReLU'))
model.add(Dropout(rate=0.5))
model.add(Dense(units=32, activation='LeakyReLU'))
model.add(Dropout(rate=0.5))
model.add(Dense(units=1, activation='sigmoid'))

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy',
              metrics=['accuracy'])

# Train the model
model.fit(x_train, y_train, verbose=1,  
 validation_data=(x_validation, y_validation),batch_size=128, epochs=100,  
       callbacks=[EarlyStopping(patience=10, monitor="val_loss", 
 restore_best_weights=True), History()])  

# Evaluate the model
model.evaluate(x_test, y_test)
