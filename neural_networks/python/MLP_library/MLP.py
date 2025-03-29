import cupy
import numpy as np

from aliases.common import (
    Tensor3D, Matrix, Vector,
    ActivationFunction, LossFunction, Optimizer, Metrics
)
from random_weights.generate import random_biases, random_weights

type WeightsGradient = Tensor3D
type BiasesGradient = Matrix


class Layer:
    def __init__(self, neurons: Vector =None, weights: Matrix =None, biases: Vector =None, activation_function : str = ActivationFunction.DEFAULT):
        if neurons is None:
            neurons = []
        if biases is None:
            biases = []
        if weights is None:
            weights = []
        self.neurons: Vector = neurons
        self.weights: Matrix = weights
        self.biases: Vector = biases
        self.activation: str = activation_function

    def normalize(self):
        ActivationFunction.normalize(self.activation, self.neurons)


class MLP:
    def __init__(self, input_size: int):
        self.biases_gradient: BiasesGradient = []
        self.weights_gradient: WeightsGradient = []

        self.metrics = Metrics.DEFAULT
        self.optimizer = Optimizer.DEFAULT
        self.loss = LossFunction.DEFAULT

        self.input_layer: Vector = np.zeros(input_size)
        self.inner_layers: list[Layer] = []

        self.dimensions: list[int] = [input_size]

    def add_layers(self, sizes: list[int]):
        for size in sizes:
            self.add_layer(size)

    def add_layer(self, size: int, activation_function: ActivationFunction = ActivationFunction.DEFAULT):
        new_layer = Layer(
            neurons=np.zeros(size),
            weights=random_weights(size, self.dimensions[-1]),
            biases=random_biases(size),
            activation_function=activation_function
        )
        self.inner_layers.append(new_layer)
        self.dimensions.append(size)

    def compile(self, loss=LossFunction.DEFAULT, optimizer=Optimizer.DEFAULT, metrics=Metrics.DEFAULT):
        self.loss = loss
        self.optimizer = optimizer
        self.metrics = metrics
        for layer in self.inner_layers:
            self.weights_gradient.append(np.zeros_like(layer.weights))
            self.biases_gradient.append(np.zeros_like(layer.biases))

    def fit(self, x: list[Vector], y: list[Vector], epochs: int = 1, batch_size: int = 32):
        if len(x) != len(y):
            print("training data and labels must have same length!")
            return

        for epoch in range(epochs):
            self.grad_zero()
            for i in range(0, len(x), batch_size):
                x_batch = x[i:min(i + batch_size, len(x))]
                y_batch = y[i:min(i + batch_size, len(y))]
                self.train_on_batch(x_batch, y_batch)
            self.step()

    def train_on_batch(self, x: list[Vector], y: list[Vector]) -> None:
        for i in range(len(x)):
            self.forward(x[i])
            self.backward(y[i]) # backward updates the gradient for us

    def forward(self, x: Vector) -> None:
        if type(x) != np.ndarray:
            x = np.array(x)
        if len(self.inner_layers) == 0:
            raise Exception("No layers were added!")

        # here we would need to optionally copy our weight matrices onto GPU
        # probably by several workers run in parallel to different segments of the GPU

        self.input_layer = x
        self.inner_layers[0].neurons = (
                self.inner_layers[0].weights @ self.input_layer
                + self.inner_layers[0].biases
        )

        for (index, layer) in enumerate(self.inner_layers[1:]):
            previous_layer = self.inner_layers[index]
            layer.neurons = (
                layer.weights @ previous_layer.neurons
                + layer.biases
            )
            layer.normalize()


    # performs one backward pass updating the gradient
    def backward(self, y: list[float]) -> None:
        pass

    def step(self):
        pass

    def grad_zero(self):
        pass

    def output_layer(self) -> Vector:
        return self.inner_layers[-1].neurons


# from joblib import Parallel, delayed
# import cupy as cp
# import time


# def gpu_task(matrix, vector, gpu_id):
#     with cp.cuda.Device(gpu_id):
#         square = cp.asarray(matrix)
#         line = cp.array(vector)
#         return cp.asnumpy(square @ line)


if __name__ == "__main__":
    model = MLP(700)
    model.add_layers([16,16,10])
    model.compile()
    model.forward([1] * 700)
    print(model.output_layer())

    # v_col = np.array([[1],[2],[3]])
    # v_row = np.array([1,2,3])
    # mat = np.array([np.arange(3), np.arange(3) + 2])
    # print (mat @ v_row)
    # print (mat[1][0])
    # # v_col.reshape(1, 3)
    # # v_row.reshape(1, 3)
    # print(mat)
    # print(v_col)
    # print(v_row)
    # dimension_size = 10000
    # sizes = [10, 100, 1000]
    # num_gpus = cp.cuda.runtime.getDeviceCount()
    # print(num_gpus)
    # print()
    # for size in sizes:
    #     matrices = [np.random.rand(dimension_size, dimension_size) for _ in range(size)]
    #     vectors = [np.random.rand(dimension_size) for _ in range(size)]
    #     start = time.time()
    #     for i, (m, v) in enumerate(zip(matrices, vectors)):
    #         gpu_task(m, v, i % num_gpus)
    #
    #     stop = time.time()
    #     print(stop - start)
    #
    # print()
    # for size in sizes:
    #     matrices = [np.random.rand(dimension_size, dimension_size) for _ in range(size)]
    #     vectors = [np.random.rand(dimension_size) for _ in range(size)]
    #     start = time.time()
    #     for i, (m, v) in enumerate(zip(matrices, vectors)):
    #         dummy = m @ v
    #     stop = time.time()
    #     print (stop - start)



