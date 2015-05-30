'''

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014 MoDeNa Consortium, All rights reserved.

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

Description
    Module providing strategies

Authors
    Henrik Rusche

Contributors
    Sigve Karolius
    Mandar Thombre
'''

import six
import abc
import modena
from fireworks.core.firework import FireTaskMeta
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_serializers import FWSerializable, \
    recursive_serialize, recursive_deserialize, serialize_fw
from fireworks.utilities.fw_utilities import explicit_serialize
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


class InitialisationStrategy(defaultdict, FWSerializable):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPoints(self):
        raise NotImplementedError("newPoints not implemented!")


    @abc.abstractmethod
    def newPointsWorkflow(self):
        raise NotImplementedError("newPointsWorkflow not implemented!")


    @serialize_fw
    @recursive_serialize
    def to_dict(self):
        return dict(self)


    @classmethod
    @recursive_deserialize
    def from_dict(cls, m_dict):
        return cls(m_dict)


    def __repr__(self):
        return '<{}>:{}'.format(self.fw_name, dict(self))


class OutOfBoundsStrategy(defaultdict, FWSerializable):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPointsFWAction(self, model, caller, **kwargs):
        raise NotImplementedError("newPointsFWAction not implemented!")


    @serialize_fw
    @recursive_serialize
    def to_dict(self):
        return dict(self)


    @classmethod
    @recursive_deserialize
    def from_dict(cls, m_dict):
        return cls(m_dict)


    def __repr__(self):
        return '<{}>:{}'.format(self.fw_name, dict(self))


class ImproveErrorStrategy(defaultdict, FWSerializable):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPointsFWAction(self, model):
        raise NotImplementedError("newPointsFWAction not implemented!")


    @serialize_fw
    @recursive_serialize
    def to_dict(self):
        return dict(self)


    @classmethod
    @recursive_deserialize
    def from_dict(cls, m_dict):
        return cls(m_dict)


    def __repr__(self):
        return '<{}>:{}'.format(self.fw_name, dict(self))


class ParameterFittingStrategy(dict, FWSerializable):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPointsFWAction(self, model, **kwargs):
        raise NotImplementedError("newPointsFWAction not implemented!")


    @serialize_fw
    @recursive_serialize
    def to_dict(self):
        return dict(self)


    @classmethod
    @recursive_deserialize
    def from_dict(cls, m_dict):
        return cls(m_dict)


    def __repr__(self):
        return '<{}>:{}'.format(self.fw_name, dict(self))


class SamplingStrategy():

    def newPoints(self):
        raise NotImplementedError("newPoints not implemented!")


    def newPointsWorkflow(self):
        raise NotImplementedError("newPointsWorkflow not implemented!")


    def samplePoints(self, model, sampleRange, nPoints):

        points = array(lhs.randomLHS(nPoints, len(sampleRange))).tolist()

        sr = sampleRange
        return {
            key: [
                sr[key]['min'] +
                (sr[key]['max'] - sr[key]['min']) * points[i][j]
                for i in xrange(nPoints)
            ] for j, key in enumerate(sr)
        }


@explicit_serialize
class InitialPoints(InitialisationStrategy):

    def __init__(self, *args, **kwargs):
        InitialisationStrategy.__init__(self, *args, **kwargs)


    def newPoints(self):
        return self['initialPoints']


    def newPointsWorkflow(self):
        return self['initialPoints']


@explicit_serialize
class ExtendSpaceStochasticSampling(OutOfBoundsStrategy, SamplingStrategy):

    def __init__(self, *args, **kwargs):
        OutOfBoundsStrategy.__init__(self, *args, **kwargs)


    def newPointsFWAction(self, model, caller, **kwargs):
        # Get a sampling range around the outside point and create points
        sampleRange = model.extendedRange(kwargs['outsidePoint'])
        points = self.samplePoints(model, sampleRange, self['nNewPoints'])

        # Continue with exact tasks, parameter estimation and (finally) this
        # task in order to resume normal operation
        (tl, dep, last) = model.exactTasks(points)
        tl.append(Firework(caller))
        dep[last] = tl[-1]

        return FWAction(
            detours=Workflow(tl, dep, name="test workflow")
        )


@explicit_serialize
class StochasticSampling(ImproveErrorStrategy, SamplingStrategy):

    def __init__(self, *args, **kwargs):
        ImproveErrorStrategy.__init__(self, *args, **kwargs)


    def newPointsFWAction(self, model):
        '''
        The function serves the following purposes:
        It will add samples to the current range if needed by ParFit.
        '''
        # Get a sampling range from fitData. Note: Cannot use MinMax must not
        # be updated, yet
        sampleRange = {
            k: {
                'min': min(model.fitData[k]),
                'max': max(model.fitData[k])
            } for k in model.inputs
        }

        points = self.samplePoints(model, sampleRange, self['nNewPoints'])

        # Continue with exact tasks, parameter estimation and (finally) this
        # task in order to resume normal operation
        (tl, dep, last) = model.exactTasks(points)

        return FWAction(
            detours=Workflow(tl, dep, name="test workflow")
        )


@explicit_serialize
class NonLinFitWithErrorContol(ParameterFittingStrategy):

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

        # ------------------------------ Function ----------------------------#
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
                        idxGenerator=fitData(model.nSamples, testIndices)
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
                    idxGenerator=fitData(testIndices)
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

        # perform fitting (nonlinear MSSQ)
        nlfb = nlmrt.nlfb(
            R_par,
            R_res,
            jacfn=rinterface.NULL,
            trace=rinterface.FALSE,
            maskidx=rinterface.NULL
        )


        # optimised coefficients and sum of squares
        nlfb_coeffs = nlfb[nlfb.names.index('coefficients')]
        nlfb_ssqres = nlfb[nlfb.names.index('ssquares')]
        new_parameters = list(nlfb_coeffs)

        # The following code block will check the error and call the ImproveEr-
        # rorStrategy to add more points to the design of experiments if the m-
        # odel is not validated

        maxError = errorTest(new_parameters)

        print "Maximum Error = %s" % maxError
        if maxError > self['maxError']:
            print 'Parameters ' + term.red + 'not' + term.normal + \
                ' valid, adding samples.'
            print 'current parameters = [%s]' % ' '.join(
                '%g' % k for k in new_parameters
            )

            # Update database
            model.save()

            return self['improveErrorStrategy'].newPointsFWAction(model)

        else:
            print('old parameters = [%s]' % ' '.join(
                '%g' % k for k in model.parameters)
            )
            print('new parameters = [%s]' % ' '.join(
                '%g' % k for k in new_parameters)
            )

            # Update database
            # TODO: This is not save if the model is updated by another fitting
            #                                         task running concurrently
            model.parameters = new_parameters
            model.updateMinMax()
            model.save()

            # return nothing to restart normal operation
            return FWAction()

