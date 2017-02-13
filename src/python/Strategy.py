'''@cond
   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2016 MoDeNa Consortium, All rights reserved.

License
    This file is part of Modena.

    The Modena interface library is free software; you can redistribute it
    and/or modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.
@endcond'''

"""
@namespace python.Strategy
@brief     Module providing strategies
@details

@author    Henrik Rusche
@author    Sigve Karolius
@author    Mandar Thombre
@copyright 2014-2016, MoDeNa Project. GNU Public License.
"""

import six
import abc
import sys
import copy
import modena
from fireworks import Firework, Workflow, FWAction, FireTaskBase, ScriptTask
from fireworks.utilities.fw_serializers import FWSerializable, \
    recursive_serialize, recursive_deserialize, serialize_fw
import fireworks.utilities.fw_serializers as fw_serializers
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.utilities.fw_serializers import load_object
from collections import defaultdict
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects.vectors import FloatVector
from numpy import array
from numpy.random import choice, seed
from blessings import Terminal

# Import R libraries
rinterface.initr()
nlmrt = importr('nlmrt')
lhs = importr('lhs')

# Create terminal for colour output
term = Terminal()

##
# @addtogroup python_interface_library
# @{

class StrategyBaseClass(defaultdict, FWSerializable):
    """
    @brief   Base class for all strategies
    @details
             The purpose of the base class is to ensure that all strategy types
             used in MoDeNa are embedded correctly into the FireWorks workflow.
    """


    @abc.abstractmethod
    def newPoints(self, model):
        """Method which adds new points to the database."""
        raise NotImplementedError('newPoints not implemented!')


    @abc.abstractmethod
    def workflow(self, model):
        """
        @brief    Method which adds new points to the database.
        @details
                  f
        @param    model modena SurrogateModel object
        """
        raise NotImplementedError('workflow not implemented!')


    @serialize_fw
    @recursive_serialize
    def to_dict(self):
        """
        @brief Required by FireWorks to deserialise objects
        """
        return dict(self)


    @classmethod
    @recursive_deserialize
    def from_dict(cls, m_dict):
        """
        @brief Required by FireWorks to serialise objects
        """
        return cls(m_dict)


    def __repr__(self):
        return '<{}>:{}'.format(self.fw_name, dict(self))


class InitialisationStrategy(StrategyBaseClass):
    """
    @brief    Parent class for the initialisation strategies.
    @details
              The purpose of the initialisation strategy is to initialise the
              surrogate model, i.e. compile the source code and obtain a set of
              validated parameters.
    """

    def __init__(self, *args, **kwargs):
        """
        @brief Constructor
        """
        dict.__init__(self, *args, **kwargs)


    def workflow(self, model):
        """
        @brief    Create a FireWorks Workflow object performing initialisation.
        @details
                  The workflow

        @param model surrogate model object.

        @return Workflow object
        """
        ## Call the newPoints method to receive a list of dictionaries each
        #  dictionary representing one data point.
        p = self.newPoints(model)
        if len(p):
            wf = model.exactTasks(p)
            wf.append_wf(
                model.parameterFittingStrategy().workflow(model),
                wf.leaf_fw_ids
            )
            return wf

        elif not len(p) and len(model.substituteModels):
            wf = Workflow([])
            for sm in model.substituteModels:
                wf.append_wf(
                    sm.initialisationStrategy().workflow(sm),
                    []
                )
            return wf

        else:
            return Workflow([])


class OutOfBoundsStrategy(StrategyBaseClass):
    """
    @brief    Base class for the out of bounds strategies.
    @details
              Classes inheriting this class must implement the newPoints
    """
    def __init__(self, *args, **kwargs):
        """Constructor"""
        dict.__init__(self, *args, **kwargs)



    def workflow(self, model, **kwargs):
        """
        @brief    Generating a workflow
        @details
                  The workflow generated
                  1. Extend and sample domain
                  2. Perform detailed simulations
                  3. Perform parameter fitting
        @returns wf Workflow object.
        """
        wf = model.exactTasks(self.newPoints(model, **kwargs))
        wf.append_wf(
            model.parameterFittingStrategy().workflow(model),
            wf.leaf_fw_ids
        )
        return wf


class ImproveErrorStrategy(StrategyBaseClass):
    """
    @brief    Base class for strategies 'fixing' the error of a surrogate model
    @details
              Im
    """

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    def workflow(self, model, **kwargs):
        wf = model.exactTasks(self.newPoints(model))
        wf.append_wf(
            model.parameterFittingStrategy().workflow(model),
            wf.leaf_fw_ids
        )
        return wf


class ParameterFittingStrategy(StrategyBaseClass):
    """
    @brief   Base Class for creating parameter fitting strategies.
    @details
    """
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    def workflow(self, model):
        return Workflow(
            [
                Firework(
                    ParameterFitting(surrogateModelId=model._id),
                    name='parameter fitting'
                )
            ]
        )


    def errorTest(model, parameters, testIndices):
        """
        @param model MoDeNa surrogate model
        @param parameters list
        @param testIndicies
        """
        def fitData(testIndices):
            for i in testPoint:
                yield i

        # Instantiate the surrogate model
        cModel = modena.libmodena.modena_model_t(
             model,
             parameters=list(parameters)
        )

        return max(
            abs(i) for i in model.error(
                cModel,
                idxGenerator=fitData(testPoint),
                checkBounds=False
            )
        )


    # errorFit function can only take a single arguemnt (parameters) when it
    # is called from R. Using wrapper class instead!
    class errorFit:
        """
        """
        def __init__(self, *args, **kwargs):
            self.model = args[0]
            self.testPoint = args[1]

        def function(parameters):

            def fitData(n, testPoint):
                for i in xrange(n):
                    if i not in testPoint:
                         yield i

            # Instantiate the surrogate model
            cModel = modena.libmodena.modena_model_t(
                model=model,
                parameters=list(parameters)
            )

            return FloatVector(
                list(
                    model.error(
                        cModel,
                        idxGenerator=fitData(model.nSamples, self.testPoint),
                        checkBounds=False
                    )
                )
            )





class SamplingStrategy(StrategyBaseClass):
    """
    @brief    Base class for Sampling strategies (DoE).
    @details
              Sampling
    """


    def samplePoints(self, model, sr, nPoints):
        """
        @brief    Generate "n" sample points in a domain
        @details
                  The sample points are used as inputs to detailed simulations
        @param    model -- SurrogateModel -- Required | Surrogate Model
        @param    sr -- dictionary -- Required |
                  sample range: { 'key1': {'min': float, 'max': float}, ... }
        @param    nPoints -- int -- Required | Number of sample points
        @returns  dictionary
                  {'key1': [ * , ^ , < ] , 'key2': [ * , ^ , < ] , ... }
        """

        sampleRange = {
            k: {
                'min': min(model.fitData[k]),
                'max': max(model.fitData[k])
            } for k in model.inputs.keys()
        }

        points = array(lhs.randomLHS(nPoints, len(sampleRange))).tolist()

        return {
            key: [
                sr[key]['min'] +
                (sr[key]['max'] - sr[key]['min']) * points[i][j]
                for i in xrange(nPoints)
            ] for j, key in enumerate(sr)
        }






@explicit_serialize
class InitialPoints(InitialisationStrategy):
    """
    @brief    Initialise by performing detailed simulations at a set of points.
    @details
              The sample data points are specified by the user as a dictionary,
    """

    def __init__(self, *args, **kwargs):
        InitialisationStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model):
        return self['initialPoints']


@explicit_serialize
class InitialRange(InitialisationStrategy, SamplingStrategy):
    """
    @brief    Initialise by performing detailed simulations inside a range
    @details
              The range specified by the user is sampled using LHS sampling.
    """

    def __init__(self, *args, **kwargs):
        InitialisationStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model):
        """
        @brief Range
        @param model MoDeNa surrogate model
        """
        sampleRange = self['initialRange']
        return self.samplePoints(model, sampleRange, 10)#self['nNewPoints']


@explicit_serialize
class InitialData(InitialisationStrategy):
    """
    @brief    Initialise a SurrogateModel given a dataset of input-output data
    @details
              The purpose of this strategy is to initialise a surrogate model
              by providing a set of initial data-points, this can be results
              from validated simulations data or experimental data.

              The idea is to instantiate a surrogate model in a domain where
              the input-output behaviour is known and let the framework handle
              the expansion of the surrogate model beyond the initial domain.
    """

    def __init__(self, *args, **kwargs):
        InitialisationStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model):
        """
        @brief    Function providing the user-specified initial points.
        @details
                  This strategy requires that the data points are given in a
                  dictionary structure.
        @returns  List
        """
        return self['initialData']


    def workflow(self, model):
        """
        """
        # Get initial data
        points = self.newPoints(model)

        # Save initial data in database
        model.updateFitDataFromFwSpec(points)
        model.updateMinMax()
        model.save()

        wf = Workflow( [], name='initialising to dataset')
        wf.append_wf( model.parameterFittingStrategy().workflow(model),
                      wf.leaf_fw_ids)

        return wf


@explicit_serialize
class EmptyInitialisationStrategy(InitialisationStrategy):
    """
    @brief    Empty initialisation strategy, used by Forward Mapping Models.
    @details
              The strategy is used in SurrogateModel.ForwardMappingModel as the
              default initialisation strategy. The reason is that the parent
              class requires a surrogate model to implement an initialisation
              strategy.
    """
    def newPoints(self, model):
        return []


@explicit_serialize
class ExtendSpaceStochasticSampling(OutOfBoundsStrategy, SamplingStrategy):
    """
    @brief    Class for extending the design space using stochastic sampling.
    @details
              Strategy used to extend the domain of a surrogate model.
    """
    def __init__(self, *args, **kwargs):
        OutOfBoundsStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model, **kwargs):
        # Get a sampling range around the outside point and create points
        sampleRange, limitPoint = model.extendedRange(kwargs['outsidePoint'])
        sp = self.samplePoints(model, sampleRange, self['nNewPoints']-1)
        return  { k: v + [limitPoint[k]] for k, v in sp.iteritems() }


@explicit_serialize
class StochasticSampling(ImproveErrorStrategy, SamplingStrategy):
    """
    @brief    Design of experiments class, Monte Carlo sampling.
    @details
              Th
    """
    def __init__(self, *args, **kwargs):
        ImproveErrorStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model):
        """
        @brief    Add samples to the current range if needed by ParFit.
        """
        # Get a sampling range from fitData. Note: Cannot use MinMax must not
        # be updated, yet
        sampleRange = {
            k: {
                'min': min(model.fitData[k]),
                'max': max(model.fitData[k])
            } for k in model.inputs.keys()
        }

        return self.samplePoints(model, sampleRange, self['nNewPoints'])



#from scipy.optimize import leastsq

#from sklearn.model_selection import LeaveOneOut

# leastsq(func, x0, args=(), Dfun=None, full_output=0, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=None, factor=100, diag=None)

@explicit_serialize
class LeaveOneOut(ParameterFittingStrategy):

    def __init__(self, *args, **kwargs):
        """
        """
        ParameterFittingStrategy.__init__(self, *args, **kwargs)


    def split(self, nSamples):
        """
        """
        vs = range(nSamples) # validation samples
        ts = ( tuple(vs[0:i] + vs[i+1:]) for (i, si) in enumerate(vs) )
        return vs


    def validation(sel):
        pass


    def fit(self, model, valIndices):
        pass


    def newPointsFWAction(self, model, **kwargs):
        """
        * Get training and validation sets
        * Regression
        * Validation
        * if valid, return empty workflow, else perform DoE.
        """

        #loo = LeaveOneOut()
        points = [ [ 0.0 for j in xrange(model.nSamples) ] for i in \
                                 xrange(model.surrogateFunction.inputs_size())]
        for (k, v) in model.inputs.iteritems():
            points[model.inputs_argPos(k)] = inputs[k]

        # Transpose
        points = zip(*points)
        testIndicies = self.split()


@explicit_serialize
class Test(ParameterFittingStrategy):


    def __init__(self, *args, **kwargs):
        """
        """
        ParameterFittingStrategy.__init__(self, *args, **kwargs)


    def validationSets(self, n, testIndices):
        for i in xrange(n):
            if i not in testIndices:
                 yield i



    def split(self, nSamples):
        """
        s = [ ( ( training_set ), validation_point ), ... ]
        """
        s = range(nSamples) # validation samples
        training_samples = ( tuple(s[0:i]+s[i+1:]) for (i, si) in enumerate(s) )
        return ( training_samples, s )


    def fit(self, model, testIndices):
        """
        """
        # ------------------------------ Function --------------------------- #
        def errorFit(parameters):

            def fitData(n, testIndices):
                for i in xrange(n):
                    if i not in testIndices:
                         yield i

            # Instantiate the surrogate model
            cModel = modena.libmodena.modena_model_t(model=model,parameters=list(parameters))

            return FloatVector(list(model.error(cModel,idxGenerator=fitData(model.nSamples, testIndices),checkBounds=False)))
        # ------------------------------------------------------------------- #

        new_parameters = model.parameters
        if not len(new_parameters):
            new_parameters = [None] * len(model.surrogateFunction.parameters)
            for k, v in model.surrogateFunction.parameters.iteritems():
                new_parameters[v.argPos] = (v.min + v.max)/2

        # make objects usable in R
        R_par = FloatVector(new_parameters)
        R_res = rinterface.rternalize(errorFit)

        max_parameters = [None]*len(new_parameters)
        min_parameters = [None]*len(new_parameters)
        for k, v in model.surrogateFunction.parameters.iteritems():
            min_parameters[v.argPos] = v.min
            max_parameters[v.argPos] = v.max

        # perform fitting (nonlinear MSSQ)
        nlfb = nlmrt.nlfb(start=R_par,resfn=R_res,jacfn=rinterface.NULL,trace=rinterface.FALSE,lower=FloatVector(min_parameters),upper=FloatVector(max_parameters),maskidx=rinterface.NULL)

        # optimised coefficients and sum of squares
        nlfb_coeffs = nlfb[nlfb.names.index('coefficients')]
        nlfb_ssqres = nlfb[nlfb.names.index('ssquares')]
        new_parameters = list(nlfb_coeffs)

        return new_parameters


    def validate(self, model, parameters, testIndices):
        """
        """

        # ------------------------------ Function --------------------------- #
        def errorTest(parameters):

            def fitData(testIndices):
                for i in testIndices:
                     yield i

            # Instantiate the surrogate model
            cModel = modena.libmodena.modena_model_t(model,parameters=list(parameters))

            return max(abs(i) for i in model.error(cModel,idxGenerator=fitData(testIndices),checkBounds=False))
        # ------------------------------------------------------------------- #
        return errorTest(parameters)



    def newPointsFWAction(self, model, **kwargs):
        """
        * Get training and validation sets
        * Regression
        * Validation
        * if valid, return empty workflow, else perform DoE.
        """

        training_sets, validation_sets = self.split(model.nSamples)

        parameters = [ self.fit(model, [v_set] ) for v_set in validation_sets ]
        errors = [ self.validate(model, pi, [vi] ) for (pi, vi) in zip(parameters, validation_sets) ]

        maxError = min(errors)
        new_parameters = parameters[errors.index(maxError)]

        print 'Maximum Error = %s' % maxError
        if maxError > self['maxError']:
            print('Parameters ' + term.red + 'not' + term.normal + ' valid, adding samples.')
            print('current parameters = [%s]' % ', '.join('%g' % k for k in new_parameters))

            # Update database
            model.save()

            return FWAction(detours=self['improveErrorStrategy'].workflow(model))

        else:
            print('old parameters = [%s]' % ', '.join('%g' % k for k in model.parameters))
            print('new parameters = [%s]' % ', '.join('%g' % k for k in new_parameters))

            model["parameters"] = new_parameters
            model.updateMinMax()
            model.save()

            # return nothing to restart normal operation
            return FWAction()



@explicit_serialize
class NonLinFitWithErrorContol(ParameterFittingStrategy):
    """
    @brief    Parameter fitting class, non-linear least squares regression.
    @details
              The Strategy
    """

    def __init__(self, *args, **kwargs):
        """
        @todo access tuple correctly
        """
        #if '_fw_name' in args[0]:
        #    ParameterFittingStrategy.__init__(self, *args, **kwargs)

        #if not kwargs.has_key('improveErrorStrategy'):
        #    raise Exception('Need improveErrorStrategy')
        #if not isinstance(
        #    kwargs['improveErrorStrategy'], ImproveErrorStrategy
        #):
        #    raise TypeError('Need improveErrorStrategy')

        ParameterFittingStrategy.__init__(self, *args, **kwargs)


    def newPointsFWAction(self, model, **kwargs):
        # Make sure we get new samples in deterministic manner
        seed(model.nSamples)

        # TODO: The rest of this function should become a method in model (or /
        #                                                       strategy class)

        # Create indices a subset (~20% of samples) for testing
        testIndices = set(
            choice(
                model.nSamples,
                size=max(1, self['testDataPercentage']*model.nSamples),
                replace=False
            )
        )

        # ------------------------------ Function --------------------------- #
        def errorFit(parameters):

            def fitData(n, testIndices):
                for i in xrange(n):
                    if i not in testIndices:
                         yield i

            # Instantiate the surrogate model
            cModel = modena.libmodena.modena_model_t(
                model=model,
                parameters=list(parameters)
            )

            return FloatVector(
                list(
                    model.error(
                        cModel,
                        idxGenerator=fitData(model.nSamples, testIndices),
                        checkBounds=False
                    )
                )
            )
        # ------------------------------------------------------------------- #

        # ------------------------------ Function --------------------------- #
        def errorTest(parameters):

            def fitData(testIndices):
                for i in testIndices:
                     yield i

            # Instantiate the surrogate model
            cModel = modena.libmodena.modena_model_t(
                model,
                parameters=list(parameters)
            )

            return max(
                abs(i) for i in model.error(
                    cModel,
                    idxGenerator=fitData(testIndices),
                    checkBounds=False
                )
            )
        # ------------------------------------------------------------------- #

        new_parameters = model.parameters
        if not len(new_parameters):
            new_parameters = [None] * len(model.surrogateFunction.parameters)
            for k, v in model.surrogateFunction.parameters.iteritems():
                new_parameters[v.argPos] = (v.min + v.max)/2

        # make objects usable in R
        R_par = FloatVector(new_parameters)
        R_res = rinterface.rternalize(errorFit)

        max_parameters = [None]*len(new_parameters)
        min_parameters = [None]*len(new_parameters)
        for k, v in model.surrogateFunction.parameters.iteritems():
            min_parameters[v.argPos] = v.min
            max_parameters[v.argPos] = v.max

        # perform fitting (nonlinear MSSQ)
        nlfb = nlmrt.nlfb(
            start=R_par,
            resfn=R_res,
            jacfn=rinterface.NULL,
            trace=rinterface.FALSE,
            lower=FloatVector(min_parameters),
            upper=FloatVector(max_parameters),
            maskidx=rinterface.NULL
        )
        del max_parameters
        del min_parameters

        # optimised coefficients and sum of squares
        nlfb_coeffs = nlfb[nlfb.names.index('coefficients')]
        nlfb_ssqres = nlfb[nlfb.names.index('ssquares')]
        new_parameters = list(nlfb_coeffs)

        # The following code block will check the error and call the ImproveEr-
        # rorStrategy to add more points to the design of experiments if the m-
        # odel is not validated

        maxError = errorTest(new_parameters)

        print 'Maximum Error = %s' % maxError
        if maxError > self['maxError']:
            print(
                'Parameters ' + term.red + 'not' + term.normal
              + ' valid, adding samples.'
            )
            print(
                'current parameters = [%s]' % ', '.join(
                    '%g' % k for k in new_parameters
                )
            )

            # Update database
            model.save()

            return FWAction(
                detours=self['improveErrorStrategy'].workflow(model)
            )

        else:
            print(
                'old parameters = [%s]' % ', '.join(
                    '%g' % k for k in model.parameters
                )
            )
            print(
                'new parameters = [%s]' % ', '.join(
                    '%g' % k for k in new_parameters
                )
            )

            # Update database
            # TODO: This is not save if the model is updated by another fitting
            #                                         task running concurrently
            model.parameters = new_parameters
            model.updateMinMax()
            model.save()

            # return nothing to restart normal operation
            return FWAction()


@explicit_serialize
class NonLinFitToPointWithSmallestError(ParameterFittingStrategy):
    """
    Performs parameter fitting of a set of samples and returns the parameters
    that yield the smallest error.

    @todo The strategy does **not** allow for error to be improved, but this
          can be changed in the future.
    """

    def __init__(self, *args, **kwargs):

        # TODO: access tuple correctly
        #if '_fw_name' in args[0]:
        #    ParameterFittingStrategy.__init__(self, *args, **kwargs)

        #if not kwargs.has_key('improveErrorStrategy'):
        #    raise Exception('Need improveErrorStrategy')
        #if not isinstance(
        #    kwargs['improveErrorStrategy'], ImproveErrorStrategy
        #):
        #    raise TypeError('Need improveErrorStrategy')

        ParameterFittingStrategy.__init__(self, *args, **kwargs)


    def newPointsFWAction(self, model, **kwargs):

        new_parameters = model.parameters
        if not len(new_parameters):
            new_parameters = [None] * len(model.surrogateFunction.parameters)
            for k, v in model.surrogateFunction.parameters.iteritems():
                new_parameters[v.argPos] = (v.min + v.max)/2

        max_parameters = [None]*len(new_parameters)
        min_parameters = [None]*len(new_parameters)
        for k, v in model.surrogateFunction.parameters.iteritems():
            min_parameters[v.argPos] = v.min
            max_parameters[v.argPos] = v.max

        maxError = 1000
        coeffs = None
        for i in xrange(model.nSamples):
            testPoint = [i]

            # -------------------------- Function --------------------------- #
            def errorTest(parameters):

                def fitData(testIndices):
                    for i in testPoint:
                         yield i

                # Instantiate the surrogate model
                cModel = modena.libmodena.modena_model_t(
                    model,
                    parameters=list(parameters)
                )

                return max(
                    abs(i) for i in model.error(
                        cModel,
                        idxGenerator=fitData(testPoint),
                        checkBounds=False
                    )
                )

            # -------------------------- Function --------------------------- #
            def errorFit(parameters):

                def fitData(n, testPoint):
                    for i in xrange(n):
                        if i not in testPoint:
                             yield i

                # Instantiate the surrogate model
                cModel = modena.libmodena.modena_model_t(
                    model=model,
                    parameters=list(parameters)
                )

                return FloatVector(
                    list(
                        model.error(
                            cModel,
                            idxGenerator=fitData(model.nSamples, testPoint),
                            checkBounds=False
                        )
                    )
                )
            # --------------------------------------------------------------- #

            # make objects usable in R
            R_res = rinterface.rternalize(errorFit)
            R_par = FloatVector(new_parameters)

            # perform fitting (nonlinear MSSQ)
            nlfb = nlmrt.nlfb(
                start=R_par,
                resfn=R_res,
                jacfn=rinterface.NULL,
                trace=rinterface.FALSE,
                lower=FloatVector(min_parameters),
                upper=FloatVector(max_parameters),
                maskidx=rinterface.NULL
            )

            parameterError = errorTest(nlfb[nlfb.names.index('coefficients')])
            if parameterError < maxError:
                maxError = parameterError
                new_parameters = list(nlfb[nlfb.names.index('coefficients')])
                coeffs = nlfb


        # optimised coefficients and sum of squares
        nlfb_coeffs = coeffs[nlfb.names.index('coefficients')]
        nlfb_ssqres = coeffs[nlfb.names.index('ssquares')]

        print 'Maximum Error = %s' % maxError
        print(
            'old parameters = [%s]' % ', '.join(
                '%g' % k for k in model.parameters
            )
        )
        print(
            'new parameters = [%s]' % ', '.join(
                '%g' % k for k in new_parameters
            )
        )

        # Update database
        # TODO: This is not save if the model is updated by another fitting
        #                                         task running concurrently
        model.parameters = new_parameters
        model.updateMinMax()
        model.save()

        del parameterError
        del max_parameters
        del min_parameters
        del coeffs

        # return nothing to restart normal operation
        return FWAction()


@explicit_serialize
class Initialisation(FireTaskBase):
    """
    @brief    Defines a computational, i.e. Firetask, performing initialisation
    @details
              The firework loads the surrogate model and serialises the
              initialisation strategy.

    @author Henrik Rusche
    """

    def __init__(self, *args, **kwargs):
        FireTaskBase.__init__(self, *args, **kwargs)

        if kwargs.has_key('surrogateModel'):
            if isinstance(kwargs['surrogateModel'], modena.SurrogateModel):
                self['surrogateModelId'] = kwargs['surrogateModel']['_id']
                del self['surrogateModel']


    def run_task(self, fw_spec):
        """
        @brief    Method called by Fireworks in order to run the Firetask
        @params   fw_spec (dict) parameters passed to the Firetask
        """
        try:
            print term.cyan + 'Performing initialisation' + term.normal
            model = modena.SurrogateModel.load(self['surrogateModelId'])
            return FWAction(
                detours=model.initialisationStrategy().workflow(model)
            )
        except:
            import traceback
            traceback.print_exc()
            return FWAction(defuse_workflow=True)


@explicit_serialize
class ParameterFitting(FireTaskBase):
    """
    @brief    Defines the computational task performing parameter estimation
    @details
              The purpose of this class is to load the "parameter estimation
              strategy" from a surrogate model.

    @author Henrik Rusche
    """

    def __init__(self, *args, **kwargs):
        FireTaskBase.__init__(self, *args, **kwargs)

        if kwargs.has_key('surrogateModel'):
            if isinstance(kwargs['surrogateModel'], modena.SurrogateModel):
                self['surrogateModelId'] = kwargs['surrogateModel']['_id']
                del self['surrogateModel']


    def run_task(self, fw_spec):
        """
        @brief    Method called by Fireworks in order to run the Firetask
        @params   fw_spec (dict) parameters passed to the Firetask
        """
        try:
            model = modena.SurrogateModel.load(self['surrogateModelId'])
            print(
                term.cyan
              + 'Performing parameter fitting for model %s' % model._id
              + term.normal
            );
            model.updateFitDataFromFwSpec(fw_spec)
            return model.parameterFittingStrategy().newPointsFWAction(model)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return FWAction(defuse_workflow=True)


class MoDeNaLaunchPad(object):


    def __init__(self):
        pass


    def resolve(self):
        pass


    def run_task(self):
        pass


class OutOfBounds(Exception):
    def __init__(self, *args, **kwargs):
        print args
        print kwargs
        print(
            term.cyan
          + '%s out-of-bounds, executing outOfBoundsStrategy for model %s'
            % (args[0], args[1]._id)
          + term.normal
        )
        super(OutOfBounds, self).__init__(*args, **kwargs)


class ParametersNotValid(Exception):
    pass


class TerminateWorkflow(Exception):
    pass


class ModifyWorkflow(Exception):
    pass


@explicit_serialize
class ModenaFireTask(FireTaskBase):
    """
    @brief    Defines a computational task for detailed simulations.
    @details
    """

    def executeAndCatchExceptions(self, op, text):
        """
        @brief    Method executing tasks and catching callbacks
        @details
                  The
        """
        try:
            op()

        except OutOfBounds as e:
            model = e.args[1]
            print(
                term.cyan
              + '%s out-of-bounds, executing outOfBoundsStrategy for model %s'
                % (text, model._id)
              + term.normal
            )

            # Continue with exact tasks, parameter estimation and (finally) this
            # task in order to resume normal operation
            wf = model.outOfBoundsStrategy().workflow(
                model,
                outsidePoint=model.outsidePoint
            )
            wf.append_wf(
                Workflow([Firework(self)], name='original task'),
                wf.leaf_fw_ids
            )
            raise ModifyWorkflow(FWAction(detours=wf))

        except ParametersNotValid as e:
            model = e.args[1]
            print(
                term.cyan
              + '%s is not initialised, ' % text
              + 'executing initialisationStrategy for model %s' % model._id
              + term.normal
            )

            # Continue with exact tasks, parameter estimation and (finally) this
            # task in order to resume normal operation
            wf = model.initialisationStrategy().workflow(model)
            wf.append_wf(
                Workflow([Firework(self)], name='original task'),
                wf.leaf_fw_ids
            )
            raise ModifyWorkflow(FWAction(detours=wf))


    def run_task(self, fw_spec):
        """
        @brief    Method called by Fireworks in order to run the Firetask
        @params   fw_spec (dict) parameters passed to the Firetask
        @details
                  The Firetask performs **one** detailed simulation.

                  The Firetask reads the input, i.e. "point" from the "fw_spec"
                  input.

                  The execution of the simulation is wrapped inside a lambda
                  function and sent to the method executeAndCatchExceptions
                  which captures and handles error callbacks.

        @returns  FWAction object telling Fireworks how to proceed.
        """
        try:
            print(
                term.yellow
              + 'Performing exact simulation (microscopic code recipe) for model '
              + self['modelId']
              + term.normal
            )

            p = self['point']

            print(
                term.yellow
              + 'point = {%s}' % ', '.join(
                    '%s: %g' % (k, v) for (k, v) in p.iteritems()
                )
              + term.normal
            )

            model = modena.SurrogateModel.load(self['modelId'])
            oldP = copy.copy(p)
            for m in model.substituteModels:
                self.executeAndCatchExceptions(
                    lambda: p.update(m.callModel(p)),
                    'Substituted Model'
                )

            if not len(p) == len(oldP):
                print(
                    term.yellow
                  + 'values added by substitution = {%s}' % ', '.join(
                      '%s: %g' % (k, v) for (k, v) in p.iteritems()
                      if k not in oldP
                    )
                  + term.normal
                )

            self.executeAndCatchExceptions(
                lambda: self.task(fw_spec),
                'Model'
            )
            return FWAction(mod_spec=[{'_push': self['point']}])

        except ModifyWorkflow as e:
            return e.args[0]

        except Exception as e:
            print(term.red + e.args[0] + term.normal)
            import traceback
            traceback.print_exc()
            return FWAction(defuse_workflow=True)

        return FWAction()


    def handleReturnCode(self, returnCode):
        """
        @brief    Handle return code caught in executeAndCatchExceptions
        @params   returnCode (integer) error code from simulation
        @details
                  The method is called with an integer, i.e. the error code, as
                  the argument and raises the appropriate Python exception if
                  the error code is MoDeNa-related.

                  | Error code | Exception                 |
                  | ---------- | ------------------------- |
                  | 200        | Model is Out of bounds    |
                  | 201        | Model is not in database  |
                  | 202        | Parameters not validated  |

        """
        if returnCode > 0:
            print(term.red + 'return code = %i' % returnCode + term.normal)

        if returnCode == 200:
            try:
                # TODO
                # Finding the 'failing' model using the outsidePoint will fail
                # eventually fail when running in parallel. Need to pass id of
                # calling FireTask. However, this requires additional code in the
                # library as well as cooperation of the recipie

                model = modena.SurrogateModel.loadFailing()
            except:
                raise TerminateWorkflow(
                    'Exact task raised OutOfBounds signal, '
                  + 'but failing model could not be determined',
                    returnCode
                )

            raise OutOfBounds(
                'Exact task raised OutOfBounds signal',
                model,
                returnCode
            )

        elif returnCode == 201:
            try:
                model = modena.SurrogateModel.loadFromModule()
            except:
                raise TerminateWorkflow(
                    'Exact task raised LoadFromModule signal, '
                  + 'but failing model could not be determined',
                    returnCode
                )
            raise ParametersNotValid(
                'Exact task raised LoadFromModule signal',
                model,
                returnCode
            )

        elif returnCode == 202:
            try:
                model = modena.SurrogateModel.loadParametersNotValid()
            except:
                raise TerminateWorkflow(
                    'Exact task raised ParametersNotValid, '
                  + 'but failing model could not be determined',
                    returnCode
                )
            raise ParametersNotValid(
                'Exact task raised ParametersNotValid',
                model,
                returnCode
            )

        elif returnCode > 0:
            raise TerminateWorkflow(
                'An unknow error occurred calling exact simulation',
                returnCode
            )


@explicit_serialize
class BackwardMappingScriptTask(ModenaFireTask, ScriptTask):
    """
    @brief  FireTask that starts a macroscopic code and catches its return code
    @author Henrik Rusche
    """
    required_params = ['script']

    def run_task(self, fw_spec):
        """
        """
        try:
            print(
                term.yellow
                + 'Performing backward mapping simulation '
                + '(macroscopic code recipe)'
                + term.normal
            )
            self.executeAndCatchExceptions(
                lambda: self.task(fw_spec),
                'Model'
            )
            print(term.green + 'Success - We are done' + term.normal)

        except ModifyWorkflow as e:
            return e.args[0]

        except Exception as e:
            print(term.red + e.args[0] + term.normal)
            import traceback
            traceback.print_exc()
            return FWAction(defuse_workflow=True)

        return FWAction()


    def task(self, fw_spec):

        self['defuse_bad_rc'] = True

        # Execute the macroscopic code by calling function in base class
        ret = ScriptTask.run_task(self, fw_spec)
        self.handleReturnCode(ret.stored_data['returncode'])

##
# @} # end of python_interface_library

