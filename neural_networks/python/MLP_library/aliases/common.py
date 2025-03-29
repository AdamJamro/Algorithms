from numpy import (
    ndarray
)

type Matrix = list[list[float]] | ndarray
type Tensor3D = list[list[list[float]]] | ndarray
type Vector = list[float] | ndarray


class ActivationFunction:
    SIGMOID = 'sigmoid'
    RELU = 'relu'
    SOFTMAX = 'softmax'
    TANH = 'tanh'
    DEFAULT = RELU

    @staticmethod
    def normalize(function_type: str, neurons: Vector) -> Vector:
        if function_type == ActivationFunction.SIGMOID:
            return ActivationFunction.sigmoid(neurons)
        elif function_type == ActivationFunction.RELU:
            return ActivationFunction.relu(neurons)
        elif function_type == ActivationFunction.SOFTMAX:
            return ActivationFunction.softmax(neurons)
        elif function_type == ActivationFunction.TANH:
            return ActivationFunction.tanh(neurons)
        else:
            raise ValueError(f"Unknown activation function: {function_type}")

    @staticmethod
    def sigmoid(neurons: Vector) -> Vector:
        pass

    @staticmethod
    def relu(neurons: Vector) -> Vector:
        pass

    @staticmethod
    def softmax(neurons: Vector) -> Vector:
        pass

    @staticmethod
    def tanh(neurons: Vector) -> Vector:
        pass


class LossFunction:
    BINARY_CROSSENTROPY = 'binary_crossentropy'
    CATEGORICAL_CROSSENTROPY = 'categorical_crossentropy'
    MEAN_SQUARED_ERROR = 'mean_squared_error'
    DEFAULT = BINARY_CROSSENTROPY


class Optimizer:
    SGD = 'sgd'
    ADAM = 'adam'
    RMSPROP = 'rmsprop'
    ADAGRAD = 'adagrad'
    RAW = 'raw'
    DEFAULT = ADAM


class Metrics:
    ACCURACY = 'accuracy'
    PRECISION = 'precision'
    RECALL = 'recall'
    F1_SCORE = 'f1_score'
    AUC = 'auc'
    DEFAULT = ACCURACY



