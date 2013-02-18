import numpy as np
import combivep.settings as cbv_const


class Mlp(object):
    """MultiLayer Perceptron class"""

    def __init__(self,
                 n_features,
                 seed=cbv_const.DFLT_SEED,
                 n_hidden_nodes=cbv_const.DFLT_HIDDEN_NODES):
        object.__init__(self)
        #set initial configuration values and memorize input
        self.n_features = n_features
        self.n_hidden_nodes = n_hidden_nodes
        self.best_weights1 = []
        self.best_weights2 = []

        #set initial values of weight matrixs to random small values
        np.random.seed(seed)
        self.weights1 = 0.01 * np.random.rand(n_hidden_nodes,
                                              n_features+1)
        self.weights2 = 0.01 * np.random.rand(1,
                                              n_hidden_nodes+1)

        #set initial values of momentum matrixs to zeros
        self.momentums1 = np.zeros((n_hidden_nodes, n_features+1))
        self.momentums2 = np.zeros((1, n_hidden_nodes+1))

    def forward_propagation(self, dataset):
        #calculate sum of product in the hidden layer
        in1 = np.dot(self.weights1,
                     np.concatenate((dataset.feature_vectors,
                                     np.ones((1, dataset.n_data))
                                     ),
                                    axis=0
                                    )
                     )

        #calculate outputs of hidden layer using non-linear function
        self.out1 = np.concatenate((2/(1+np.exp(-in1))-1,
                                    np.ones((1, dataset.n_data))
                                    ),
                                   axis=0
                                   )

        #calculate sum of product in the output node
        in2 = np.dot(self.weights2, self.out1)

        #calculate output of mlp using non-linear function
        self.out2 = 1/(1+np.exp(-in2))

        #return prediction result
        return self.out2

    def backward_propagation(self, training_dataset):
        model_error = np.multiply(self.calc_error(self.out2,
                                                  training_dataset.targets),
                                  training_dataset.n_data/2)
        self.error_sig_output = np.multiply(model_error,
                                            np.multiply((1-self.out2),
                                                        self.out2
                                                        )
                                            )
        self.error_sig_hidden = np.multiply(np.dot(self.weights2.T,
                                                   self.error_sig_output
                                                   ),
                                            np.multiply((1+self.out1),
                                                        (1-self.out1)
                                                        )
                                            ) * 0.5
        self.error_sig_hidden = self.error_sig_hidden[0:self.n_hidden_nodes]

        return np.sum(np.absolute(model_error), axis=1).item(0)

    def weight_update(self,
                      training_dataset,
                      coefficient=cbv_const.MLP_COEFFICIENT,
                      step_size=cbv_const.STEP_SIZE):
        features = training_dataset.feature_vectors
        n_data = training_dataset.n_data
        transformed_input_features = np.concatenate((features,
                                                     np.ones((1, n_data))
                                                     ),
                                                    axis=0
                                                    ).T
        self.momentums1 = np.subtract((self.momentums1*coefficient),
                                      (np.dot(self.error_sig_hidden,
                                              transformed_input_features
                                              )
                                       )*(1-coefficient)
                                      )
        self.momentums2 = np.subtract((self.momentums2*coefficient),
                                      (np.dot(self.error_sig_output,
                                              self.out1.T
                                              )
                                       )*(1-coefficient)
                                      )
        self.weights1 = np.add(self.weights1,
                               np.multiply(self.momentums1,
                                           step_size)
                               )
        self.weights2 = np.add(self.weights2,
                               np.multiply(self.momentums2,
                                           step_size)
                               )
        return self.weights1, self.weights2

    def calc_error(self, actual_output, expected_output):
        error = np.subtract(actual_output, expected_output)
        n_positive_samples = expected_output[expected_output == 1].shape[0]
        n_negative_samples = expected_output[expected_output == 0].shape[0]
        error = np.where(expected_output == 1,
                         np.divide(error, n_positive_samples*2),
                         error)
        error = np.where(expected_output == 0,
                         np.divide(error, n_negative_samples*2),
                         error)
        return error

    def export_best_parameters(self, params_file=cbv_const.USER_PARAMS_FILE):
        np.savez(params_file,
                 best_weights1=self.best_weights1,
                 best_weights2=self.best_weights2)

    def import_parameters(self, params_file=cbv_const.USER_PARAMS_FILE):
        params = np.load(params_file)
        self.weights1 = params['best_weights1']
        self.weights2 = params['best_weights2']
