import os.path
import numpy as np
import matplotlib.pyplot as plt
import combivep.settings as cbv_const
from combivep.engine.mlp import Mlp

class Trainer(Mlp):
    """

    This class is to produce parameters which will be used later by
    the Predictor class

    """


    def __init__(self,
                 training_dataset,
                 validation_dataset,
                 seed=cbv_const.DEFAULT_SEED, 
                 n_hidden_nodes=cbv_const.DEFAULT_HIDDEN_NODES, 
                 figure_dir=cbv_const.DEFAULT_FIGURE_DIR):
        Mlp.__init__(self, training_dataset.n_features,
                                        seed=seed,
                                        n_hidden_nodes=n_hidden_nodes)

        self.__training_dataset   = training_dataset
        self.__validation_dataset = validation_dataset
        self.__n_hidden_nodes     = n_hidden_nodes
        self.__figure_dir         = figure_dir

    def train(self, iterations=cbv_const.DEFAULT_ITERATIONS):
        training_error   = []
        validation_error = []
        running_round    = 0
        best_validation_error = 0.99
        while True:
            #tune up model parameters
            out = self.forward_propagation(self.__training_dataset)
            self.backward_propagation(self.__training_dataset)
            weights1, weights2 = self.weight_update(self.__training_dataset)
            errors = self.calc_error(out, self.__training_dataset.targets)
            training_error.append(np.sum(np.absolute(errors),
                                         axis=1
                                         ).item(0)
                                  )

            #evaluate model using validation dataset
            out = self.forward_propagation(self.__validation_dataset)
            errors = self.calc_error(out, self.__validation_dataset.targets)
            validation_error.append(np.sum(np.absolute(errors),
                                           axis=1
                                           ).item(0)
                                    )

            #check ending condition
            #(acceptable error rate and not much improvement in each iteration)
            current_validation_error = validation_error[len(validation_error)-1]
            if (current_validation_error < cbv_const.MAX_ALLOWED_ERROR):
                improvement = best_validation_error - current_validation_error
                if (improvement < cbv_const.MIN_IMPROVEMENT):
                    break

            #otherwise save parameters and record last error
            best_validation_error = validation_error[len(validation_error)-1]
            self.best_weights1 = weights1
            self.best_weights2 = weights2

            #check if it reach maximum iteration
            running_round += 1
            if running_round >= iterations:
                break

        if self.__figure_dir:
            self.__save_figure(training_error, validation_error)

        return best_validation_error

    def __save_figure(self, training_error, validation_error):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(training_error, label='training')
        ax.plot(validation_error, label='validation')
        ax.set_ylabel('average error')
        ax.set_xlabel('iterations')
        ax.legend(bbox_to_anchor=(0, 0, 0.98, 0.98), loc=1, borderaxespad=0.1)
        file_name = "%02d.eps" % (self.__n_hidden_nodes)
        fig.savefig(os.path.join(self.__figure_dir, file_name))

class Predictor(Mlp):
    """

    This class is to predict a probability how a variant likely to be deleterious.

    """


    def __init__(self):
        pass

    def predict(self, dataset):
        return self.forward_propagation(dataset)


