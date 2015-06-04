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
import sys
import modena
from fireworks.core.firework import FireTaskMeta
from fireworks import Firework, Workflow, FWAction, FireTaskBase, ScriptTask
from fireworks.utilities.fw_serializers import FWSerializable, \
    recursive_serialize, recursive_deserialize, serialize_fw
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


class Workflow2(Workflow):

    def __init__(self, *args, **kwargs):
        Workflow.__init__(self, *args, **kwargs)

    def addAfterAll(self, wf):
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

        for fw_id in my_leaf_ids:
            self.links[fw_id].extend(root_ids)  # add the root id as my child

    def addNoLink(self, wf):
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

        print self.links


class InitialisationStrategy(defaultdict, FWSerializable):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPoints(self):
        raise NotImplementedError("newPoints not implemented!")


    def workflow(self, model):
        p = self.newPoints()
        if len(p):
            wf = model.exactTasks(p)
            wf.addAfterAll(model.parameterFittingStrategy().workflow(model))
            return wf
        else:
            return Workflow2([])


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
    def newPoints(self):
        raise NotImplementedError("newPoints not implemented!")


    def workflow(self, model, **kwargs):
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


class ImproveErrorStrategy(defaultdict, FWSerializable):

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

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)


    @abc.abstractmethod
    def newPointsFWAction(self, model, **kwargs):
        raise NotImplementedError("newPointsFWAction not implemented!")

    def workflow(self, model):
        return Workflow2(
            [
                Firework(
                    ParameterFitting(surrogateModelId=model._id),
                    name='parameter fitting'
                )
            ]
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

    def newPoints(self):
        raise NotImplementedError("newPoints not implemented!")


    def samplePoints(self, model, sampleRange, nPoints):

        points = array(lhs.randomLHS(nPoints, len(sampleRange))).tolist()

        sr = sampleRange
        return {
            key: [
                sr[key]['min'] +
                (sr[key]['max'] - sr[key]['min']) * points[i][j]
                for i in xrange(nPoints)
            ] for (j, key) in enumerate(sr)
        }


@explicit_serialize
class InitialPoints(InitialisationStrategy):

    def __init__(self, *args, **kwargs):
        InitialisationStrategy.__init__(self, *args, **kwargs)


    def newPoints(self):
        return self['initialPoints']


@explicit_serialize
class InitialData(InitialisationStrategy):
    """
    Class initialising a SurrogateModel given a dataset of input-output
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

             for m in model.substituteModels:
                p.update(m.callModel(p))

        t = et
        t['point'] = p
        fw = Firework(t)
        wf = Workflow2( [fw], name='initialising to dataset')

        # pf = load_object({'_fw_name': '{{modena.Strategy.NonLinFitToPointWithSmallestError}}'})
        # NonLinFitToPointWithSmallestError
        wf.addAfterAll(model.parameterFittingStrategy().workflow(model))

        return wf


@explicit_serialize
class EmptyInitialisationStrategy(InitialisationStrategy):

    def newPoints(self):
        return []


@explicit_serialize
class ExtendSpaceStochasticSampling(OutOfBoundsStrategy, SamplingStrategy):

    def __init__(self, *args, **kwargs):
        OutOfBoundsStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model, **kwargs):
        # Get a sampling range around the outside point and create points
        sampleRange = model.extendedRange(kwargs['outsidePoint'])
        return self.samplePoints(model, sampleRange, self['nNewPoints'])


@explicit_serialize
class StochasticSampling(ImproveErrorStrategy, SamplingStrategy):

    def __init__(self, *args, **kwargs):
        ImproveErrorStrategy.__init__(self, *args, **kwargs)


    def newPoints(self, model):
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

        return self.samplePoints(model, sampleRange, self['nNewPoints'])


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

        print "Maximum Error = %s" % maxError
        if maxError > self['maxError']:
            print 'Parameters ' + term.red + 'not' + term.normal + \
                ' valid, adding samples.'
            print 'current parameters = [%s]' % ' '.join(
                '%g' % k for k in new_parameters
            )

            # Update database
            model.save()

            return FWAction(
                detours=self['improveErrorStrategy'].workflow(model)
            )

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


@explicit_serialize
class NonLinFitToPointWithSmallestError(ParameterFittingStrategy):
    """
    Performs parameter fitting of a set of samples and returns the parameters
    that yield the smallest error.

    TODO: The strategy does **not** allow for error to be improved, but this
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
                        idxGenerator=fitData(testPoint)
                    ))
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
                            idxGenerator=fitData(model.nSamples, testPoint)
                        )))
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

        print "Maximum Error = %s" % maxError
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

        del parameterError
        del max_parameters
        del min_parameters
        del coeffs

        # return nothing to restart normal operation
        return FWAction()


@explicit_serialize
class Initialisation(FireTaskBase):
    '''
    @summary: A FireTask that performs the initialisation

    @author: Henrik Rusche
    '''

    def __init__(self, *args, **kwargs):
        FireTaskBase.__init__(self, *args, **kwargs)

        if kwargs.has_key('surrogateModel'):
            if isinstance(kwargs['surrogateModel'], modena.SurrogateModel):
                self['surrogateModelId'] = kwargs['surrogateModel']['_id']
                del self['surrogateModel']


    def run_task(self, fw_spec):
        print term.cyan + "Performing initialisation" + term.normal

        model=modena.SurrogateModel.load(self['surrogateModelId'])

        return FWAction(
            detours=model.initialisationStrategy().workflow(model)
        )


@explicit_serialize
class ParameterFitting(FireTaskBase):
    '''
    @summary: A FireTask that performs parameter fitting

    @author: Henrik Rusche
    '''

    def __init__(self, *args, **kwargs):
        FireTaskBase.__init__(self, *args, **kwargs)

        if kwargs.has_key('surrogateModel'):
            if isinstance(kwargs['surrogateModel'], modena.SurrogateModel):
                self['surrogateModelId'] = kwargs['surrogateModel']['_id']
                del self['surrogateModel']


    def run_task(self, fw_spec):
        print term.cyan + "Performing parameter fitting" + term.normal

        model=modena.SurrogateModel.load(self['surrogateModelId'])

        model.updateFitDataFromFwSpec(fw_spec)

        return model.parameterFittingStrategy().newPointsFWAction(model)


class BackwardMappingTask:

    def handleReturnCode(self, returnCode):

        # Analyse return code
        print('return code = %i' % returnCode)
        if returnCode == 201:
            print term.cyan + "Performing Initialisation" + term.normal
            # TODO
            # Finding the 'failing' model using the outsidePoint will fail
            # eventually fail when running in parallel. Need to pass id of
            # calling FireTask. However, this requires additional code in the
            # library as well as cooperation of the recipie
            model = modena.SurrogateModel.loadFromModule()

            # Continue with exact tasks, parameter estimation and (finally) this
            # task in order to resume normal operation

            print model.initialisationStrategy()
            wf = model.initialisationStrategy().workflow(model)
            wf.addAfterAll(
                Workflow2([Firework(self)], name='original task')
            )

            return FWAction(detours=wf)

        elif returnCode == 200:
            print term.cyan + "Performing Design of Experiments" + term.normal
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

        elif returnCode > 0:
            print('An error occurred')
            sys.exit(returnCode)

        else:
            print('Success - We are done')
            return FWAction()


@explicit_serialize
class InitialDataPoints(FireTaskBase):
    """
    Pushes **all** "points" to the next firework.

    TODO: See the method `updateFitDataFromFwSpec` in SurrogateModel.py
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


@explicit_serialize
class BackwardMappingScriptTask(ScriptTask, BackwardMappingTask):
    """
    A FireTask that starts a macroscopic code and catches its return code.
    @author: Henrik Rusche
    """
    required_params = ['script']

    def run_task(self, fw_spec):
        print(
            term.yellow
          + "Performing backward mapping simulation (macroscopic code recipe)"
          + term.normal
        )

        self['defuse_bad_rc'] = True

        # Execute the macroscopic code by calling function in base class
        ret = super(BackwardMappingScriptTask, self).run_task(fw_spec)

        return self.handleReturnCode(ret.stored_data['returncode'])

