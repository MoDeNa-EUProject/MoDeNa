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
    Module providing functions and models

Authors
    Henrik Rusche

Contributors
    Sigve Karolius
    Mandar Thombre
'''

import os
import six
import abc
import hashlib
import modena
from modena.Strategy import *
import weakref
import re
import sys
from mongoengine import *
from mongoengine.document import TopLevelDocumentMetaclass
from fireworks import Firework, FWAction, FireTaskBase
from collections import defaultdict
from blessings import Terminal


# Create connection to database
MODENA_URI = os.environ.get('MODENA_URI', 'mongodb://localhost:27017/test')
(uri, database) = MODENA_URI.rsplit('/', 1);
connect(database, host=MODENA_URI)

# Create terminal for colour output
term = Terminal()


def existsAndHasArgPos(i, name):
    return name in i and hasattr(i[name], 'argPos')


def checkAndConvertType(kwargs, name, cls):
    if not name in kwargs:
        raise Exception('Need ' + name)
    elif not isinstance(kwargs[name], cls):
        raise TypeError('Need ' + name)
    else:
        kwargs['meth_' + name] = kwargs[name].to_dict()
        del kwargs[name]


def loadType(obj, name, cls):
    n = '___' + name
    if(hasattr(obj, n)):
        return getattr(obj, n)
    else:
        var = getattr(obj, 'meth_' + name)
        var = load_object(var)
        setattr(obj, n, var)
        return var


class EmbDoc(DynamicEmbeddedDocument):
    meta = {'allow_inheritance': False}


class GrowingList(list):
    def __setitem__(self, index, value):
        if index >= len(self):
            self.extend([None]*(index + 1 - len(self)))
        list.__setitem__(self, index, value)


# Fitting data is not stored here to allow excluding it in load since it
# is not possible to exclude inputs.*.fitData
class MinMax(EmbeddedDocument):
    min = FloatField(required=True)
    max = FloatField(required=True)
    meta = {'allow_inheritance': False}


class MinMaxOpt(EmbeddedDocument):
    min = FloatField()
    max = FloatField()
    meta = {'allow_inheritance': False}


class MinMaxArgPos(EmbeddedDocument):
    min = FloatField(required=True, default=None)
    max = FloatField(required=True, default=None)
    argPos = IntField(required=True)
    meta = {'allow_inheritance': False}


class MinMaxArgPosOpt(EmbeddedDocument):
    min = FloatField()
    max = FloatField()
    argPos = IntField()
    meta = {'allow_inheritance': False}


class combinedMeta(abc.ABCMeta, TopLevelDocumentMetaclass):
   pass


class SurrogateFunction(DynamicDocument):

    # Database definition
    name = StringField(primary_key=True)
    parameters = MapField(EmbeddedDocumentField(MinMaxArgPos))
    functionName = StringField(required=True)
    libraryName = StringField(required=True)
    meta = {'allow_inheritance': True}

    @abc.abstractmethod
    def __init__(self, *args, **kwargs):

        for k, v in kwargs['inputs'].iteritems():
            if not isinstance(v, MinMaxArgPos):
                kwargs['inputs'][k] = MinMaxArgPos(**v)

        for k, v in kwargs['outputs'].iteritems():
            if not isinstance(v, MinMaxArgPos):
                kwargs['outputs'][k] = MinMaxArgPos(**v)

        for k, v in kwargs['parameters'].iteritems():
            if not isinstance(v, MinMaxArgPos):
                kwargs['parameters'][k] = MinMaxArgPos(**v)

        DynamicDocument.__init__(self, **kwargs)

    def modena_function(self):
        return pymodena.modena_function_new(
            str(self.functionName),
            str(self.libraryName)
        )


    @classmethod
    def load(self, surrogateFunctionId):
        return self.objects.get(_id=surrogateFunctionId)


class CFunction(SurrogateFunction):

    def __init__(self, *args, **kwargs):
        if kwargs.has_key('_cls'):
            SurrogateFunction.__init__(self, *args, **kwargs)
        if kwargs.has_key('libraryName'):
            SurrogateFunction.__init__(self, *args, **kwargs)
        else:
            if not kwargs.has_key('Ccode'):
                raise Exception('Need Ccode')

            ln = self.compileCcode(kwargs['Ccode'])
            fn = re.search(
                'void\s*(.*)\s*\('
                '\s*const\s*double\s*\*\s*parameters\s*,'
                '\s*const\s*double\s*\*\s*inherited_inputs\s*,'
                '\s*const\s*double\s*\*\s*inputs\s*,'
                '\s*double\s*\*\s*outputs\s*\)',
                kwargs['Ccode']
            ).group(1)
            fn = fn.strip(' \t\n\r')

            kwargs['name'] = fn
            kwargs['libraryName'] = ln
            kwargs['functionName'] = fn

            SurrogateFunction.__init__(self, *args, **kwargs)
            self.save()

    def compileCcode(self, Ccode):
        """
        Helper function to compile a model into local library
        """

        m = hashlib.md5()
        m.update(Ccode)
        h = m.hexdigest()
        d = "func_" + h
        ln = "%s/%s/lib%s.so" % (os.getcwd(), d, h)

        if(True or not os.path.exists(ln)):
            if(not os.path.isdir(d)): os.mkdir(d)
            os.chdir(d)

            f = open('%s.c' % h, 'w')
            f.write(Ccode)
            f.close()

            f = open('CMakeLists.txt', 'w')
            f.write("""
cmake_minimum_required (VERSION 2.6)
project (%(h)s C)

find_package(MODENA REQUIRED)
find_package(LTDL REQUIRED)

add_library(%(h)s MODULE %(h)s.c)
target_link_libraries(%(h)s MODENA::modena ${LTDL_LIBRARIES})

install(TARGETS %(h)s DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
""" % {"h": h})
            f.close()

            from subprocess import call
            call(["cmake", "."])
            call(["make"])
            os.chdir('..')

            return ln


class SurrogateModel(DynamicDocument):

    # List of all instances (for initialisation)
    ___refs___ = []

    # Database definition
    _id = StringField(primary_key=True)
    surrogateFunction = ReferenceField(SurrogateFunction, required=True)
    parameters = ListField(FloatField())
    meta = {'allow_inheritance': True}

    def __init__(self, *args, **kwargs):
        self.___refs___.append(weakref.ref(self))


    def outputsToModels(self):
        o = { k: self for k in self.outputs }
        for m in self.substituteModels:
            o.update(m.outputsToModels())
        return o


    def inputs_argPos(self, name):
        if existsAndHasArgPos(self.inputs, name):
            return self.inputs[name].argPos
        elif existsAndHasArgPos(self.surrogateFunction.inputs, name):
            return self.surrogateFunction.inputs[name].argPos
        raise Exception(name + ' not found')


    def outputs_argPos(self, name):
        if existsAndHasArgPos(self.outputs, name):
            return self.outputs[name].argPos
        elif existsAndHasArgPos(self.surrogateFunction.outputs, name):
            return self.surrogateFunction.outputs[name].argPos
        raise Exception(name + ' not found')


    def inputs_max_argPos(self):
        return max(self.inputs_argPos(k) for k in self.inputs)


    def calculate_maps(self, sm):
        map_outputs = []
        map_inputs = []

        for k in self.inputs:
            try:
                map_inputs.extend(
                    [self.inputs_argPos(k), sm.inputs_argPos(k)]
                )
            except:
                pass

        for k, v in self.surrogateFunction.inputs.iteritems():
            try:
                map_outputs.extend([sm.outputs_argPos(k), v.argPos])
            except:
                pass

        return map_outputs, map_inputs


    def minMax(self):
        len = 1 + self.inputs_max_argPos()
        minValues = [-9e99] * len
        maxValues = [9e99] * len
        for k, v in self.inputs.iteritems():
            minValues[self.inputs_argPos(k)] = v.min
            maxValues[self.inputs_argPos(k)] = v.max

        return minValues, maxValues


    def __getattribute__(self, name):
        if name.startswith( '___' ):
            return object.__getattribute__(self, name)
        else:
            return super(DynamicDocument, self).__getattribute__(name)


    def __setattribute__(self, name, value):
        if name.startswith( '___' ):
            object.__setattribute__(self, name, value)
        else:
            super(DynamicDocument, self).__setattribute__(name, value)


    @classmethod
    def exceptionModelLoad(self, surrogateModelId):
        # TODO
        # Finding the 'unitialised' models using this method will fail
        # eventually fail when running in parallel. Need to pass id of
        # calling FireTask. However, this requires additional code in the
        # library as well as cooperation of the recipie
        collection = self._get_collection()
        collection.update(
            { '_id': surrogateModelId },
            { '_id': surrogateModelId },
            upsert=True
        )
        return 201


    def exceptionOutOfBounds(self, oPoint):
        oPointDict = {
            k: oPoint[v.argPos]
            for k, v in self.inputs.iteritems()
        }
        self.outsidePoint = EmbDoc(**oPointDict)
        self.save()
        return 200


    def callModel(self, inputs):
        # Instantiate the surrogate model
        cModel = modena.libmodena.modena_model_t(model=self)

        in_i = list()
        i = [0] * (1 + self.inputs_max_argPos())

        # Set inputs
        for k, v in self.inputs.iteritems():
            i[self.inputs_argPos(k)] = inputs[k]

        # Call the surrogate model
        out = cModel.call(in_i, i)

        outputs = {
            k: out[v.argPos]
            for k, v in self.surrogateFunction.outputs.iteritems()
        };

        return outputs


    @classmethod
    def load(self, surrogateModelId):
        # Removed temporarily, probably bug in mongo engine
        #return self.objects.exclude('fitData').get(_id=surrogateModelId)
        return self.objects.get(_id=surrogateModelId)


    @classmethod
    def loadFailing(self):
        # Removed temporarily, probably bug in mongo engine
        #return self.objects(
        #    __raw__={'outsidePoint': { '$exists': True}}
        #).exclude('fitData').first()
        return self.objects(
            __raw__={'outsidePoint': { '$exists': True}}
        ).first()


    @classmethod
    def loadFromModule(self):
        collection = self._get_collection()
        doc = collection.find_one({ '_cls': { '$exists': False}})
        modelId = doc['_id']
        mod = __import__(modelId)
        # TODO:
        # Give a better name to the variable a model is imported from
        return mod.m


    @classmethod
    def get_instances(self):
        for inst_ref in self.___refs___:
            inst = inst_ref()
            if inst is not None:
                yield inst


class ForwardMappingModel(SurrogateModel):

    # Database definition
    inputs = MapField(EmbeddedDocumentField(MinMaxArgPosOpt))
    outputs = MapField(EmbeddedDocumentField(MinMaxArgPosOpt))
    substituteModels = ListField(ReferenceField(SurrogateModel))
    meta = {'allow_inheritance': True}

    def __init__(self, *args, **kwargs):
        SurrogateModel.__init__(self, *args, **kwargs)

        if kwargs.has_key('_cls'):
            DynamicDocument.__init__(self, *args, **kwargs)

        else:
            if not kwargs.has_key('_id'):
                raise Exception('Need _id')
            if not kwargs.has_key('surrogateFunction'):
                raise Exception('Need surrogateFunction')
            if not isinstance(kwargs['surrogateFunction'], SurrogateFunction):
                raise TypeError('Need surrogateFunction')

            kwargs['inherited_inputs'] = 0

            kwargs['fitData'] = {}
            kwargs['inputs'] = {}
            for k, v in kwargs['surrogateFunction'].inputs.iteritems():
                kwargs['inputs'][k] = v.to_mongo()

            kwargs['outputs'] = {}
            for k, v in kwargs['surrogateFunction'].outputs.iteritems():
                kwargs['fitData'][k] = []
                kwargs['outputs'][k] = MinMaxArgPosOpt(**{})

            subOutputs = {}
            for m in kwargs['substituteModels']:
                if not isinstance(m, SurrogateModel):
                    raise TypeError(
                        'Elements of substituteModels '
                        'must be derived from SurrogateModel'
                    )
                subOutputs.update(m.outputsToModels())

            nInp = len(kwargs['inputs'])
            replaced = {}
            while True:
                found = None
                for o in subOutputs:
                    if o in kwargs['inputs']:
                        found = o
                        break

                if found == None:
                    break

                del kwargs['inputs'][o]
                for k, v in subOutputs[o].inputs.iteritems():
                    if not k in kwargs['inputs']:
                        kwargs['inputs'][k] = { 'argPos': nInp }
                        nInp += 1

            nInputs = 0
            for k, v in kwargs['inputs'].iteritems():
                kwargs['fitData'][k] = []
                kwargs['inputs'][k] = MinMaxArgPosOpt(**v)

            if 'initialisationStrategy' not in kwargs:
                kwargs['initialisationStrategy'] = EmptyInitialisationStrategy()
            
            checkAndConvertType(
                kwargs,
                'initialisationStrategy',
                InitialisationStrategy
            );

            DynamicDocument.__init__(self, *args, **kwargs)
            self.save()

    def updateMinMax(self):
        if not self.nSamples:
            for v in self.inputs.values():
                v.min = 9e99
                v.max = -9e99

            for v in self.outputs.values():
                v.min = 9e99
                v.max = -9e99

        for k, v in self.inputs.iteritems():
            v.min = min(self.fitData[k])
            v.max = max(self.fitData[k])

        for k, v in self.outputs.iteritems():
            v.min = min(self.fitData[k])
            v.max = max(self.fitData[k])

    def updateFitDataFromFwSpec(self, fw_spec):
        # Load the fitting data
        # Removed temporarily, probably bug in mongo engine
        #self.reload('fitData')

        for k in self.inputs:
            if fw_spec[k][0].__class__ == list:
                self.fitData[k].extend(fw_spec[k][0])
            else:
                self.fitData[k].extend(fw_spec[k])
                
        for k in self.outputs:
            if fw_spec[k][0].__class__ == list:
                self.fitData[k].extend(fw_spec[k][0])
            else:
                self.fitData[k].extend(fw_spec[k])

        # Get first set
        firstSet = six.next(six.itervalues(self.fitData))
        self.nSamples = len(firstSet)

    def error(self, cModel, **kwargs):
        idxGenerator = kwargs.pop('idxGenerator', xrange(self.nSamples))

        in_i = list()
        i = [0] * (1 + self.inputs_max_argPos())

        # TODO: Deal with multivalued functions
        output = self.fitData[six.next(six.iterkeys(self.outputs))]

        for j in idxGenerator:
            # Load inputs
            for k, v in self.inputs.iteritems():
                i[v.argPos] = self.fitData[k][j]

            # Call the surrogate model
            # out = cModel.call(in_i, i, checkBounds=False)
            out = cModel.call(in_i, i)

            #print "%i %f - %f = %f" % (j, out[0], output[j], out[0] - output[j])
            yield out[0] - output[j]

    def exactTasks(self, points):
        '''
        Return an empty workflow
        '''
        return Workflow2([])

    def initialisationStrategy(self):
        return loadType(
            self,
            'initialisationStrategy',
            InitialisationStrategy
        )


    def parameterFittingStrategy(self):
        return load_object({'_fw_name': '{{modena.Strategy.NonLinFitToPointWithSmallestError}}'})

#     def initialisationStrategy(self):
#         '''
#         Return an empty workflow
#         '''
# 
#         return EmptyInitialisationStrategy()


class BackwardMappingModel(SurrogateModel):

    # Database definition
    inputs = MapField(EmbeddedDocumentField(MinMaxArgPosOpt))
    outputs = MapField(EmbeddedDocumentField(MinMaxArgPosOpt))
    fitData = MapField(ListField(FloatField(required=True)))
    substituteModels = ListField(ReferenceField(SurrogateModel))
    outsidePoint = EmbeddedDocumentField(EmbDoc)
    meta = {'allow_inheritance': True}


    # TODO: Should be able to have
    # __init__(self, _id, surrogateFunction, *args, **kwargs):
    # But this is not working with the mongoengine 0.8.7
    # Try again when > 0.9.0 comes out
    def __init__(self, *args, **kwargs):
        SurrogateModel.__init__(self, *args, **kwargs)

        if kwargs.has_key('_cls'):
            DynamicDocument.__init__(self, *args, **kwargs)

        else:
            if not kwargs.has_key('_id'):
                raise Exception('Need _id')
            if not kwargs.has_key('surrogateFunction'):
                raise Exception('Need surrogateFunction')
            if not isinstance(kwargs['surrogateFunction'], SurrogateFunction):
                raise TypeError('Need surrogateFunction')

            kwargs['inherited_inputs'] = 0

            kwargs['fitData'] = {}
            kwargs['inputs'] = {}
            for k, v in kwargs['surrogateFunction'].inputs.iteritems():
                kwargs['inputs'][k] = v.to_mongo()

            kwargs['outputs'] = {}
            for k, v in kwargs['surrogateFunction'].outputs.iteritems():
                kwargs['fitData'][k] = []
                kwargs['outputs'][k] = MinMaxArgPosOpt(**{})

            subOutputs = {}
            for m in kwargs['substituteModels']:
                if not isinstance(m, SurrogateModel):
                    raise TypeError(
                        'Elements of substituteModels '
                        'must be derived from SurrogateModel'
                    )
                subOutputs.update(m.outputsToModels())

            nInp = len(kwargs['inputs'])
            replaced = {}
            while True:
                found = None
                for o in subOutputs:
                    if o in kwargs['inputs']:
                        found = o
                        break

                if found == None:
                    break

                del kwargs['inputs'][o]
                for k, v in subOutputs[o].inputs.iteritems():
                    if not k in kwargs['inputs']:
                        kwargs['inputs'][k] = { 'argPos': nInp }
                        nInp += 1

            nInputs = 0
            for k, v in kwargs['inputs'].iteritems():
                kwargs['fitData'][k] = []
                kwargs['inputs'][k] = MinMaxArgPosOpt(**v)

            checkAndConvertType(kwargs, 'exactTask', FireTaskBase);

            checkAndConvertType(
                kwargs,
                'initialisationStrategy',
                InitialisationStrategy
            );

            checkAndConvertType(
                kwargs,
                'outOfBoundsStrategy',
                OutOfBoundsStrategy
            );

            checkAndConvertType(
                kwargs,
                'parameterFittingStrategy',
                ParameterFittingStrategy
            );

            DynamicDocument.__init__(self, *args, **kwargs)
            self.save()


    def exactTasks(self, points):
        '''
        Build a workflow to excute an exactTask for each point
        '''

        # De-serialise the exact task from dict
        et = load_object(self.meth_exactTask)

        tl = []
        e = six.next(six.itervalues(points))
        for i in xrange(len(e)):
            p = {}
            for k in points:
                p[k] = points[k][i]

            for m in self.substituteModels:
                p.update(m.callModel(p))

            t = et
            t['point'] = p
            fw = Firework(t)

            tl.append(fw)

        return Workflow2(tl, name='exact tasks for new points')


    def initialisationStrategy(self):
        return loadType(
            self,
            'initialisationStrategy',
            InitialisationStrategy
        )


    def parameterFittingStrategy(self):
        return loadType(
            self,
            'parameterFittingStrategy',
            ParameterFittingStrategy
        )


    def outOfBoundsStrategy(self):
        return loadType(
            self,
            'outOfBoundsStrategy',
            OutOfBoundsStrategy
        )


    def updateFitDataFromFwSpec(self, fw_spec):
        # Load the fitting data
        # Removed temporarily, probably bug in mongo engine
        #self.reload('fitData')

        for k in self.inputs:
            if fw_spec[k][0].__class__ == list:
                self.fitData[k].extend(fw_spec[k][0])
            else:
                self.fitData[k].extend(fw_spec[k])
                
        for k in self.outputs:
            if fw_spec[k][0].__class__ == list:
                self.fitData[k].extend(fw_spec[k][0])
            else:
                self.fitData[k].extend(fw_spec[k])

        # Get first set
        firstSet = six.next(six.itervalues(self.fitData))
        self.nSamples = len(firstSet)

    def updateMinMax(self):
        if not self.nSamples:
            for v in self.inputs.values():
                v.min = 9e99
                v.max = -9e99

            for v in self.outputs.values():
                v.min = 9e99
                v.max = -9e99

        for k, v in self.inputs.iteritems():
            v.min = min(self.fitData[k])
            v.max = max(self.fitData[k])

        for k, v in self.outputs.iteritems():
            v.min = min(self.fitData[k])
            v.max = max(self.fitData[k])
    
    def error(self, cModel, **kwargs):
        idxGenerator = kwargs.pop('idxGenerator', xrange(self.nSamples))

        in_i = list()
        i = [0] * (1 + self.inputs_max_argPos())

        # TODO: Deal with multivalued functions
        output = self.fitData[six.next(six.iterkeys(self.outputs))]

        for j in idxGenerator:
            # Load inputs
            for k, v in self.inputs.iteritems():
                 i[v.argPos] = self.fitData[k][j]

            # Call the surrogate model
            # out = cModel.call(in_i, i, checkBounds=False)
            out = cModel.call(in_i, i)

            #print "%i %f - %f = %f" % (j, out[0], output[j], out[0] - output[j])
            yield out[0] - output[j]

    def extendedRange(self, outsidePoint, expansion_factor=1.2):
        """
        @summary: Method expanding the design space. The method ONLY operates
                  on 'self.dict', this means that the database is NOT updated.
                  This is performed afterwards by 'run_task'.

                  The method will update the 'inputRanges' key in the
                  'self.dict'. Moreover, it will ensure that the min/max values
                  in'sampleRange' are consistent, meaning that the sampling is
                  performed in the correct region.

                              +------------+....+
                              |  *    *    |    .
                              |         *  |    .
                              |     *      |    .
                              |  *     *   |  X .
                              +------------+....+ <- new global max
                              ^            ^
                        global min      new min (temporary, only for sampling)

        @author:   Sigve Karolius
        @coauthor: Mandar Thombre
        @TODO:     Document...
        """

        sampleRange = {}

        for k, v in self.inputs.iteritems():
            sampleRange[k] = {}
            outsideValue = outsidePoint[k]

            # If the value outside point is outside the range, set the
            # "localdict" max to the outside point value

            if outsideValue > v['max']:
                sampleRange[k]['min'] = v['max']
                sampleRange[k]['max'] = min(
                    outsideValue*expansion_factor,
                    self.surrogateFunction.inputs[k].max
                )

            elif outsideValue < v['min']:
                sampleRange[k]['min'] = max(
                    outsideValue/expansion_factor,
                    self.surrogateFunction.inputs[k].min
                )
                sampleRange[k]['max'] = v['min']

            else:
                sampleRange[k]['min'] = v['min']
                sampleRange[k]['max'] = v['max']

        return sampleRange

