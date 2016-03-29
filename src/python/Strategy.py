'''@cond

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2015 MoDeNa Consortium, All rights reserved.

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
@file
Module providing strategies

@author    Henrik Rusche
@author    Sigve Karolius
@author    Mandar Thombre
@copyright 2014-2015, MoDeNa Project. GNU Public License.
"""

import six
import abc
import sys
import copy
import modena
from fireworks.core.firework import FireTaskMeta
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

class Workflow2(Workflow):
    """Workflow2, expanding on "Workflow" from "FireWorks"

    @brief The class specifically allows the user to add dependencies
    """
    def __init__(self, *args, **kwargs):
        Workflow.__init__(self, *args, **kwargs)

    def addAfterAll(self, wf):
        """Method which adds a the workflow "wf" as a dependency to the
        fireworks.

        @param wf "Workflow2" object.

        @var updated_ids 
        @var root_ids: Root "FireWorks" of this workflow
        @var leaf_ids: Leaf FireWorks, i.e. those with no children
        @var my_leaf_ids
        """
        updated_ids = []

        root_ids = wf.root_fw_ids
        leaf_ids = wf.leaf_fw_ids

        my_leaf_ids = self.leaf_fw_ids

        for new_fw in wf.fws:
            if new_fw.fw_id > 0:
                raise ValueError(
                    'FireWorks to add must use a negative fw_id! Got fw_id: '
                    '{}'.format(
                        new_fw.fw_id))

            self.id_fw[new_fw.fw_id] = new_fw  # add new_fw to id_fw

            if new_fw.fw_id in leaf_ids:
                self.links[new_fw.fw_id] = []
            else:
                self.links[new_fw.fw_id] = wf.links[new_fw.fw_id]
            updated_ids.append(new_fw.fw_id)

        for fw_id in my_leaf_ids:
            self.links[fw_id].extend(root_ids)  # add the root id as my child

        for new_fw in wf.fws:
            updated_ids = self.refresh(new_fw.fw_id, set(updated_ids))

    def addNoLink(self, wf):
        """Method adding "Workflow" without dependency."""
        updated_ids = []

        leaf_ids = wf.leaf_fw_ids

        for new_fw in wf.fws:
            if new_fw.fw_id > 0:
                raise ValueError(
                    'FireWorks to add must use a negative fw_id! Got fw_id: '
                    '{}'.format(
                        new_fw.fw_id))

            self.id_fw[new_fw.fw_id] = new_fw  # add new_fw to id_fw

            if new_fw.fw_id in leaf_ids:
                self.links[new_fw.fw_id] = []
            else:
                self.links[new_fw.fw_id] = wf.links[new_fw.fw_id]
            updated_ids.append(new_fw.fw_id)

        for new_fw in wf.fws:
            updated_ids = self.refresh(new_fw.fw_id, set(updated_ids))

    def printWf(self):
        for k, v in self.id_fw.iteritems():
            print k, v.spec
        print sorted(self.root_fw_ids)
        print sorted(self.leaf_fw_ids)


class InitialisationStrategy(defaultdict, FWSerializable):
    """Parent class for the initialisation strategies.

    defaultdict: subclass of type(dict), overrides method __missing__ in dict
                 and uses a method "defaultfactory" in order to automatically
                 return a dictionary "key" instead of "KeyError".
    FWSerializable: Creates a serializable object within "FireWorks"
    """
    def __init__(self, *args, **kwargs):
        """Constructor"""
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPoints(self):
        """Method which adds new points to the database."""
        raise NotImplementedError('newPoints not implemented!')


    def workflow(self, model):
        """Method creating a Workflow2 object for the initialisation.

        @param model surrogate model object.

        @var p list of dicts each representing inputs for a computation.
        @var wf Workflow2 object containing FireTasks for every point in "p".

        @return Workflow2 object
        """
        p = self.newPoints()
        if len(p):
            wf = model.exactTasks(p)
            wf.addAfterAll(model.parameterFittingStrategy().workflow(model))
            return wf

        elif not len(p) and len(model.substituteModels):
            wf = Workflow2([])
            for sm in model.substituteModels:
                wf.addNoLink(sm.initialisationStrategy().workflow(sm))
            return wf

        else:
            return Workflow2([])


    @serialize_fw
    @recursive_serialize
    def to_dict(self):
        """Method used by FireWorks to deserialise the object instance."""
        return dict(self)


    @classmethod
    @recursive_deserialize
    def from_dict(cls, m_dict):
        """Method used by FireWorks to deserialise all insatnces."""
        return cls(m_dict)


    def __repr__(self):
        return '<{}>:{}'.format(self.fw_name, dict(self))


class OutOfBoundsStrategy(defaultdict, FWSerializable):
    """Parent class for the out of bounds strategies.

    defaultdict: subclass of type(dict), overrides method __missing__ in dict
                 and uses a method "defaultfactory" in order to automatically
                 return a dictionary "key" instead of "KeyError".
    FWSerializable: Creates a serializable object within "FireWorks"
    """
    def __init__(self, *args, **kwargs):
        """Constructor"""
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPoints(self):
        raise NotImplementedError('newPoints not implemented!')


    def workflow(self, model, **kwargs):
        """Method generating the workflow for the 'out of bounds strategy'.

        @returns wf Workflow2 object.
        """
        wf = model.exactTasks(self.newPoints(model, **kwargs))
        wf.addAfterAll(model.parameterFittingStrategy().workflow(model))
        return wf


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


def recursive_serialize2(func):
    """
    a decorator to add FW serializations keys
    see documentation of FWSerializable for more details
    <https://pythonhosted.org/FireWorks/fireworks.utilities.html#fireworks.utilities.fw_serializers.FWSerializable>
    """

    def _decorator(self, *args, **kwargs):
        m_dict = func(self, *args, **kwargs)
        m_dict = fw_serializers.recursive_dict(m_dict)
        return m_dict

    return _decorator


def recursive_deserialize2(func):
    """
    a decorator to add FW serializations keys
    see documentation of FWSerializable for more details:
    <https://pythonhosted.org/FireWorks/fireworks.utilities.html#fireworks.utilities.fw_serializers.FWSerializable>
    """

    def _decorator(self, *args, **kwargs):
        new_args = [a for a in args]
        new_args[0] = {k: fw_serializers._recursive_load(v) for k, v in args[0].items()}
        m_dict = func(self, *new_args, **kwargs)
        print 'recursive_deserialize2'
        print m_dict
        return m_dict

    return _decorator


class ImproveErrorStrategy(defaultdict, FWSerializable):
    """Base class for strategies 'fixing' the error of a surrogate model."""
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    def workflow(self, model, **kwargs):
        wf = model.exactTasks(self.newPoints(model))
        wf.addAfterAll(model.parameterFittingStrategy().workflow(model))
        return wf


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
    """Base Class for creating parameter fitting strategies."""
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPointsFWAction(self, model, **kwargs):
        raise NotImplementedError('newPointsFWAction not implemented!')

    def workflow(self, model):
        return Workflow2(
            [
                Firework(
                    ParameterFitting(surrogateModelId=model._id),
                    name='parameter fitting'
                )
            ]
        )


    def errorTest(model, parameters, testIndices):

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
    """Base class for Sampling strategies (DoE)."""

    def newPoints(self):
        raise NotImplementedError('newPoints not implemented!')


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
    """Class for initialisation of a surrogate model by fitting it to
    user-specified points.
    """
    def __init__(self, *args, **kwargs):
        InitialisationStrategy.__init__(self, *args, **kwargs)


    def newPoints(self):
        return self['initialPoints']


@explicit_serialize
class InitialData(InitialisationStrategy):
    """Class initialising a SurrogateModel given a dataset of input-output
    relations.
    """
    def __init__(self, *args, **kwargs):
        InitialisationStrategy.__init__(self, *args, **kwargs)


    def newPoints(self):
        return self['initialData']


    def workflow(self, model):

        et = load_object({'_fw_name': '{{modena.Strategy.InitialDataPoints}}'})

        points = self.newPoints()

        e = six.next(six.itervalues(points))
        p = { k:[0]*len(points[k]) for k in points }
        for i in xrange(len(e)):
             for k in points:
                p[k][i] = points[k][i]

        t = et
        t['point'] = p
        t['indices'] = indices
        t['modelId'] = self._id
        fw = Firework(t)
        wf = Workflow2( [fw], name='initialising to dataset')

        wf.addAfterAll(model.parameterFittingStrategy().workflow(model))

        return wf


@explicit_serialize
class EmptyInitialisationStrategy(InitialisationStrategy):
    """Empty initialisation strategy, used by Forward Mapping Models."""
    def newPoints(self):
        return []


@explicit_serialize
class ExtendSpaceStochasticSampling(OutOfBoundsStrategy, SamplingStrategy):
    """Class for extending the design space using stochastic sampling."""
    def __init__(self, *args, **kwargs):
        OutOfBoundsStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model, **kwargs):
        # Get a sampling range around the outside point and create points
        sampleRange, limitPoint = model.extendedRange(kwargs['outsidePoint'])
        sp = self.samplePoints(model, sampleRange, self['nNewPoints']-1)
        return  { k: v + [limitPoint[k]] for k, v in sp.iteritems() }


@explicit_serialize
class StochasticSampling(ImproveErrorStrategy, SamplingStrategy):
    """Design of experiments class, Monte Carlo sampling."""
    def __init__(self, *args, **kwargs):
        ImproveErrorStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model):
        """
        The function serves the following purposes:
        It will add samples to the current range if needed by ParFit.
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


@explicit_serialize
class NonLinFitWithErrorContol(ParameterFittingStrategy):
    """Parameter fitting class, non-linear least squares regression."""
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

            # ------------------------------ Function --------------------------- #
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
            # ------------------------------------------------------------------- #

            # ------------------------------ Function ----------------------------#
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
            # ------------------------------------------------------------------- #

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
    A FireTask that performs the initialisation

    @author Henrik Rusche
    """

    def __init__(self, *args, **kwargs):
        FireTaskBase.__init__(self, *args, **kwargs)

        if kwargs.has_key('surrogateModel'):
            if isinstance(kwargs['surrogateModel'], modena.SurrogateModel):
                self['surrogateModelId'] = kwargs['surrogateModel']['_id']
                del self['surrogateModel']


    def run_task(self, fw_spec):
        print term.cyan + 'Performing initialisation' + term.normal

        model=modena.SurrogateModel.load(self['surrogateModelId'])

        return FWAction(
            detours=model.initialisationStrategy().workflow(model)
        )


@explicit_serialize
class ParameterFitting(FireTaskBase):
    """
    A FireTask that performs parameter fitting

    @author Henrik Rusche
    """

    def __init__(self, *args, **kwargs):
        FireTaskBase.__init__(self, *args, **kwargs)

        if kwargs.has_key('surrogateModel'):
            if isinstance(kwargs['surrogateModel'], modena.SurrogateModel):
                self['surrogateModelId'] = kwargs['surrogateModel']['_id']
                del self['surrogateModel']


    def run_task(self, fw_spec):
        print term.cyan + 'Performing parameter fitting' + term.normal

        model=modena.SurrogateModel.load(self['surrogateModelId'])

        model.updateFitDataFromFwSpec(fw_spec)

        return model.parameterFittingStrategy().newPointsFWAction(model)


@explicit_serialize
class InitialDataPoints(FireTaskBase):
    """
    Pushes **all** "points" to the next firework.

    @todo See the method `updateFitDataFromFwSpec` in SurrogateModel.py
          This class returns:
                     {key: array}
          this results in the following representation in fw_spec:
                     {key: [array]}
          **How to fix this?**
          The temporary solution checks the __class__ in `updateFitData...`
          however this is neither beautiful or fool proof.
    """

    def run_task(self, fw_spec):
        return FWAction(mod_spec=[{'_push': self['point']}])


class OutOfBounds(Exception):
    pass

class ParametersNotValid(Exception):
    pass
'''
    def __init__(self, *args):
        super(ParametersNotValid, self).__init__(args)
        print self.args[0]
        print 'In exception'
    def __str__(self):
        return repr(self.value)
'''

@explicit_serialize
class ModenaFireTask(FireTaskBase):
    """
    """

    def outOfBounds(self):
        """
        """

        try:
            # TODO
            # Finding the 'failing' model using the outsidePoint will fail
            # eventually fail when running in parallel. Need to pass id of
            # calling FireTask. However, this requires additional code in the
            # library as well as cooperation of the recipie
            model = modena.SurrogateModel.loadFailing()

            # Continue with exact tasks, parameter estimation and (finally) this
            # task in order to resume normal operation
            wf = model.outOfBoundsStrategy().workflow(
                model,
                outsidePoint= model.outsidePoint
            )
            wf.addAfterAll(
                Workflow2([Firework(self)], name='original task')
            )
            return FWAction(detours=wf)

        except:
            import traceback
            traceback.print_exc()
            return FWAction(defuse_children=True)


    def parametersNotValid(self, model):
        """
        """
        try:
            # Continue with exact tasks, parameter estimation and (finally) this
            # task in order to resume normal operation
            wf = model.initialisationStrategy().workflow(model)
            wf.addAfterAll(
                Workflow2([Firework(self)], name='original task')
            )
            return FWAction(detours=wf)

        except:
            import traceback
            traceback.print_exc()
            return FWAction(defuse_children=True)


    def run_task(self, fw_spec):
        """
        """
        if 'point' in self:
            print(
                term.yellow
              + 'Performing exact simulation (microscopic code recipe)'
              + term.normal
            )

            p = self['point']

            print(
                term.yellow + 'point = {%s}' % ', '.join(
                    '%s: %g' % (k, v) for (k, v) in p.iteritems()
                )
              + term.normal
            )

            model = modena.SurrogateModel.load(self['modelId'])
            oldP = copy.copy(p)
            for m in model.substituteModels:
                try:
                    p.update(m.callModel(p))

                except OutOfBounds:
                    print(
                        term.red
                      + 'Substituted model out-of-bounds, executing outOfBoundsStrategy.'
                      + term.normal
                    )
                    return self.outOfBounds()

                except ParametersNotValid, e:
                    print(
                        term.red
                      + 'Substituted model is not initialised, executing initialisationStrategy.'
                      + term.normal
                    )
                    return self.parametersNotValid(e.args[1])

            if not len(p) == len(oldP):
                print(
                    term.yellow
                  + 'values added by substitution = {%s}' % ', '.join(
                        '%s: %g' % (k, v) for (k, v) in p.iteritems()
                        if k not in oldP
                    )
                  + term.normal
                )

            try:
                self.task(fw_spec)
                return FWAction(mod_spec=[{'_push': self['point']}])

            except OutOfBounds:
                print(
                    term.red
                  + 'Model out-of-bounds, executing outOfBoundsStrategy.'
                  + term.normal
                )
                return self.outOfBounds()

            except ParametersNotValid, e:
                print term.cyan + 'Performing Initialisation' + term.normal
                return self.parametersNotValid(e.args[1])

        else:
            try:
                print(
                    term.yellow
                  + 'Performing backward mapping simulation (macroscopic code recipe)'
                  + term.normal
                )
                self.task(fw_spec)

            except OutOfBounds:
                print(
                    term.red
                  + 'Model out-of-bounds, executing outOfBoundsStrategy.'
                  + term.normal
                )
                return self.outOfBounds()

            except ParametersNotValid, e:
                print(
                    term.cyan
                  + 'Model not initialised, executing initialisationStrategy'
                  + term.normal
                )
                return self.parametersNotValid(e.args[1])

            print('Success - We are done')
            return FWAction()

        return FWAction()


    def handleReturnCode(self, returnCode):

        # Analyse return code and raise appropriate exception
        if returnCode > 0:
            print(term.red + 'return code = %i' % returnCode + term.normal)

        if returnCode == 200:
            raise OutOfBounds('Exact task of model returned 200')

        elif returnCode == 201:
            model = modena.SurrogateModel.loadFromModule()
            raise ParametersNotValid('Exact task of model returned 201', model)

        elif returnCode > 0:
            print('An unknow error occurred')
            sys.exit(returnCode)


@explicit_serialize
class BackwardMappingScriptTask(ModenaFireTask, ScriptTask):
    """
    A FireTask that starts a macroscopic code and catches its return code.
    @author Henrik Rusche
    """
    required_params = ['script']

    def task(self, fw_spec):

        self['defuse_bad_rc'] = True

        # Execute the macroscopic code by calling function in base class
        ret = ScriptTask.run_task(self, fw_spec)
        self.handleReturnCode(ret.stored_data['returncode'])

##
# @} # end of python_interface_library

