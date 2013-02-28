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
                 seed=cbv_const.DFLT_SEED,
                 n_hidden_nodes=cbv_const.DFLT_HIDDEN_NODES,
                 figure_dir=cbv_const.DFLT_FIGURE_DIR):
        Mlp.__init__(self,
                     training_dataset.n_features,
                     seed=seed,
                     n_hidden_nodes=n_hidden_nodes)

        self.training_data   = training_dataset
        self.validation_data = validation_dataset
        self.n_hidden_nodes  = n_hidden_nodes
        self.figure_dir      = figure_dir

    def train(self, iterations=cbv_const.DFLT_ITERATIONS):
        training_error   = []
        validation_error = []
        running_round    = 0
        best_validation_error = 0.99
        while True:
            #tune up model parameters
            training_out = self.forward_propagation(self.training_data)
            self.backward_propagation(self.training_data)
            weights1, weights2 = self.weight_update(self.training_data)
            errors = self.calc_error(training_out,
                                     self.training_data.targets)
            training_error.append(np.sum(np.absolute(errors),
                                         axis=1
                                         ).item(0)
                                  )

            #evaluate model using validation dataset
            validation_out = self.forward_propagation(self.validation_data)
            errors = self.calc_error(validation_out,
                                     self.validation_data.targets)
            validation_error.append(np.sum(np.absolute(errors),
                                           axis=1
                                           ).item(0)
                                    )

            #check ending condition
            #(acceptable error rate and not much improvement in each iteration)
            current_validation_error = validation_error[-1]
            if (current_validation_error < cbv_const.MAX_ALLOWED_ERROR):
                improvement = best_validation_error - current_validation_error
                if (improvement < cbv_const.MIN_IMPROVEMENT):
                    self.best_weights1 = weights1
                    self.best_weights2 = weights2
                    self.min_training_out = np.amin(training_out)
                    self.max_training_out = np.amax(training_out)
                    break

            #otherwise save parameters and record last error
            best_validation_error = validation_error[-1]

            #check if it reach maximum iteration
            running_round += 1
            if running_round >= iterations:
                self.best_weights1 = weights1
                self.best_weights2 = weights2
                self.min_training_out = np.amin(training_out)
                self.max_training_out = np.amax(training_out)
                break

        if self.figure_dir:
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
        file_name = "%02d.eps" % (self.n_hidden_nodes)
        fig.savefig(os.path.join(self.figure_dir, file_name))


class Predictor(Mlp):
    """

    This class is to predict a probability
    how a variant likely to be deleterious.

    """

    def __init__(self):
        pass

    def predict(self, dataset):
        out = self.forward_propagation(dataset)

        #scale output
        out = self.__scale(out)
#        self.max_training_out = 9
#        self.min_training_out = 3
#        self.__scale(np.matrix([0.3, 0.4, 0.6, 0.9]))

        #bring back those that go over boundaries
        out = np.where(out > 1, 1, out)
        out = np.where(out < 0, 0, out)

        return out

    def __scale(self, vals):
        """scale vals according to min and max training output"""

#        print vals
#        left_scale = 0.5 / (0.5-self.min_trainin_out)
#        print left_scale
#        low_vals = vals[vals < 0.5]


        scale = self.max_training_out - self.min_training_out
        mid_training_out = (self.max_training_out+self.min_training_out) / 2
        vals = np.add(np.divide(np.subtract(vals,
                                            mid_training_out
                                            ),
                                scale
                                ),
                      0.5
                      )
        return vals
