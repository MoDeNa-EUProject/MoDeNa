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
Module providing functions and models

@author    Henrik Rusche
@author    Sigve Karolius
@author    Mandar Thombre
@copyright 2014-2015, MoDeNa Project. GNU Public License.
"""

import os
import six
import abc
import hashlib
import modena
from modena.Strategy import *
import weakref
import re
import random
from mongoengine import *
from mongoengine.document import TopLevelDocumentMetaclass
from fireworks import Firework, FireTaskBase
from collections import defaultdict


# Create connection to database
MODENA_URI = os.environ.get('MODENA_URI', 'mongodb://localhost:27017/test')
(uri, database) = MODENA_URI.rsplit('/', 1)
connect(database, host=MODENA_URI)

##
# @addtogroup python_interface_library
# @{

def existsAndHasArgPos(i, name):
    if name not in i or 'argPos' not in i[name]:
        raise Exception('[%s][\'argPos\'] not found' % name)
    return i[name]['argPos']


def checkAndConvertType(kwargs, name, cls):
    try:
        if not isinstance(kwargs[name], cls):
            raise TypeError('%s must be of type %s' % (name, cls))
        kwargs['meth_' + name] = kwargs[name].to_dict()
        del kwargs[name]
    except:
        raise Exception('%s not found' % name)


def loadType(obj, name, cls):
    #print 'In loadType', name
    n = '___' + name
    if hasattr(obj, n):
        return getattr(obj, n)
    else:
        var = getattr(obj, 'meth_' + name)
        #print obj._get_changed_fields()
        var = load_object(var)
        #print obj._get_changed_fields()
        setattr(obj, n, var)
        return var


class EmbDoc(DynamicEmbeddedDocument):
    meta = {'allow_inheritance': False}


class GrowingList(list):
    def __setitem__(self, index, value):
        if index >= len(self):
            self.extend([None]*(index + 1 - len(self)))
        list.__setitem__(self, index, value)


class IndexSet(Document):

    # Database definition
    name = StringField(primary_key=True)
    names = ListField(StringField(required=True))
    meta = {'allow_inheritance': True}

    @abc.abstractmethod
    def __init__(self, *args, **kwargs):
        self.___index___ = {j: i for i, j in enumerate(kwargs['names'])}
        Document.__init__(self, *args, **kwargs)
        self.save()


    def get_name(self, index):
        try:
            return self.names[index]
        except:
            raise Exception('%i is not in index set %s' % (index, self.name))


    def get_index(self, name):
        try:
            return self.___index___[name]
        except:
            raise Exception('%s is not in index set %s' % (name, self.name))


    def iterator_end(self):
        return len(self.names)


    @classmethod
    def exceptionLoad(self, indexSetId):
        return 401


    @classmethod
    def load(self, indexSetId):
        return self.objects.get(name=indexSetId)


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
    #index = IndexSet()
    #index = ReferenceField(IndexSet)
    argPos = IntField(required=True)
    meta = {'allow_inheritance': False}

    def __init__(self, *args, **kwargs):
        EmbeddedDocument.__init__(self, **kwargs)


class MinMaxArgPosOpt(EmbeddedDocument):
    min = FloatField()
    max = FloatField()
    argPos = IntField()
    meta = {'allow_inheritance': False}


class SurrogateFunction(DynamicDocument):

    # Database definition
    name = StringField(primary_key=True)
    parameters = MapField(EmbeddedDocumentField(MinMaxArgPos))
    functionName = StringField(required=True)
    libraryName = StringField(required=True)
    indices = MapField(ReferenceField(IndexSet))
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

        for k in self.inputs.keys():
            self.checkVariableName(k)

        for k in self.outputs.keys():
            self.checkVariableName(k)

        for k in self.parameters.keys():
            self.checkVariableName(k)


    def indexSet(self, name):
        return self.indices[name]


    def checkVariableName(self, name):
        m = re.search('\[(.*)\]', name)
        if m and not m.group(1) in self.indices:
            raise Exception('Index %s not defined' % m.group(1))


    @classmethod
    def exceptionLoad(self, surrogateFunctionId):
        return 201


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
        d = 'func_' + h
        ln = '%s/%s/lib%s.so' % (os.getcwd(), d, h)

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

set(CMAKE_BUILD_TYPE Debug)

find_package(MODENA REQUIRED)
find_package(LTDL REQUIRED)

add_library(%(h)s MODULE %(h)s.c)
target_link_libraries(%(h)s MODENA::modena ${LTDL_LIBRARIES})

install(TARGETS %(h)s DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
""" % {'h': h})
            f.close()

            from subprocess import call
            call(['cmake', '.'])
            call(['make'])
            os.chdir('..')

            return ln

class Function(CFunction):
    """
    @todo This is a draft for a class that can parse simple functions and
          write the Ccode that is compiled by CFunction.
    """
    def __init__(self, *args, **kwargs):
        if kwargs.has_key('_cls'):
            SurrogateFunction.__init__(self, *args, **kwargs)
        if kwargs.has_key('libraryName'):
            SurrogateFunction.__init__(self, *args, **kwargs)
        else:
            # This is a bad check, make a better one...
            if not kwargs.has_key('function'):
                raise Exception('Algebraic representation not found')

        # lambda function writing inputs, parameters
        cDouble = lambda VAR: '\n'.join(['const double %s = %s[%s];' \
                                         %(V, VAR, kwargs[VAR][V]['argPos'] )\
                                                        for V in kwargs[VAR]])

        # lambda function parsing 'function' and writing outputs
        outPut = lambda OUT: '\n'.join(['outputs[%s] = %s;' \
                                 %(kwargs['outputs'][O]['argPos'],\
                                           self.Parse(kwargs['function'][O]))\
                                                     for O in kwargs[OUT] ] )

        # Main body of the Ccode
        Ccode='''
#include "modena.h"
#include "math.h"

void {name}
(
    const double* parameters,
    const double* inherited_inputs,
    const double* inputs,
    double *outputs
)
{{
{inputs}
{parameters}
{outputs}
}}
'''
        kwargs['Ccode'] = Ccode.format(\
                              name=kwargs['function']['name'],\
                              inputs=cDouble('inputs'),\
                              parameters=cDouble('parameters'),\
                              outputs=outPut('outputs')\
                          )

        CFunction.__init__(self, *args, **kwargs)

    def Parse(self, formula, debug=False, model='', stack={}, delim=0, \
              var=r'[A-Za-z]+\d*',add=r'\+',sub=r'-',mul=r'\*',\
              div=r'/',pow=r'\^',dig=r'\d+\.?\d*'\
        ):
        operators=r'%s|%s|%s|%s|%s' %(add,sub,mul,div,pow)
        ldel=r'\('
        rdel=r'\)'

        #Test explicitly for empty string. Returning error.
        # ----------------------- function ------------------------------ #
        empty = re.match('\s',formula)
        if empty:
            print 'Error: The string is empty'
            return

        formula = re.sub(r'\s+','',formula)

        # Initialise a dictionary stack.
        stack = stack or {}

        # Python has no  switch - case construct.  Match all possibilities first and
        # test afterwards:
        re_var = re.match(var,formula)
        re_dig = re.match(dig,formula)
        re_ldel = re.match(ldel,formula)
        re_rdel = re.match(rdel,formula)
        re_oper = re.match(operators,formula)

        # Parameter followed by an optional number. Allow 'p' or 'p0' as variable names
        if re_var:
            tail = formula[len(re_var.group(0)):]
            head = re_var.group(0)

        elif re_dig:
            tail = formula[len(re_dig.group(0)):]
            head = re_dig.group(0)

        elif re_oper:
            head = re_oper.group(0)
            tail = formula[1:]

        # Left delimiter.
        elif re_ldel:
            head = re_ldel.group(0)
            tail  = formula[1:]
            delim += 1

        # Right delimiter followed by an optional number (default is 1).
        elif re_rdel:
            head = re_rdel.group(0)
            tail  = formula[len(re_rdel.group(0)):]
            delim -= 1

            # Testing if there is a parenthesis imbalance.
            if delim < 0:
                raise Exception('Unmatched parenthesis.')

        # Wrong syntax. Returning an error message.
        else:
            raise Exception('The expression syntax is not suported.')

        model += head

        # The formula has not been consumed yet. Continue recursive parsing.
        if len(tail) > 0:
            return self.Parse(tail,debug,model,stack,delim)

        # Nothing left to parse. Stop recursion.
        else:
            return model


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


    def parseIndices(self):
        indices = {}
        m = re.search('.*\[(.*)\]', self._id)
        if m:
            for exp in m.group(1).split(','):
                m = re.search('(.*)=(.*)', exp)
                if m:
                    indices[m.group(1)] = m.group(2)
                else:
                    raise Exception('Unable to parse %s' % exp)

        return indices


    def outputsToModels(self):
        o = { k: self for k in self.outputs }
        for m in self.substituteModels:
            o.update(m.outputsToModels())
        return o


    def inputs_argPos(self, name):
        try:
            return existsAndHasArgPos(self.inputs, name)
        except:
            try:
                return existsAndHasArgPos(
                   self.surrogateFunction.inputs,
                   name
                )
            except:
                raise Exception(name + ' not found in inputs')


    def outputs_argPos(self, name):
        try:
            return existsAndHasArgPos(self.outputs, name)
        except:
            try:
                return existsAndHasArgPos(
                    self.surrogateFunction.outputs,
                    name
                )
            except:
                raise Exception(name + ' not found in outputs')


    def parameters_argPos(self, name):
        try:
            return existsAndHasArgPos(self.parameters, name)
        except:
            try:
                return existsAndHasArgPos(
                    self.surrogateFunction.parameters,
                    name
                )
            except:
                raise Exception(name + ' not found in parameters')


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

        #print 'maps: output =', map_outputs, 'input =', map_inputs
        return map_outputs, map_inputs


    def minMax(self):
        minValues = [-9e99] * len(self.surrogateFunction.inputs)
        maxValues = [9e99] * len(self.surrogateFunction.inputs)
        for k, v in self.inputs.iteritems():
            minValues[self.inputs_argPos(k)] = v.min
            maxValues[self.inputs_argPos(k)] = v.max

        #print 'min =', minValues, 'max =', maxValues
        return minValues, maxValues


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
            v.max = max(max(self.fitData[k]), v.min*1.000001)

        for k, v in self.outputs.iteritems():
            v.min = min(self.fitData[k])
            v.max = max(self.fitData[k])


    def error(self, cModel, **kwargs):
        idxGenerator = kwargs.pop('idxGenerator', xrange(self.nSamples))
        checkBounds = kwargs.pop('checkBounds', True)

        in_i = list()
        i = [0] * len(self.surrogateFunction.inputs)

        # TODO: Deal with multivalued functions
        output = self.fitData[six.next(six.iterkeys(self.outputs))]

        for idx in idxGenerator:
            # Load inputs
            for k, v in self.inputs.iteritems():
                i[v.argPos] = self.fitData[k][idx]

            #print 'i = {', ', '.join('%s: %g' % (
            #    k, self.fitData[k][idx]
            #) for k in self.inputs.keys()), '}'
            #print 'i =', str(i)

            # Call the surrogate model
            out = cModel.call(in_i, i, checkBounds= checkBounds)

            #print '%i %g - %g = %g' % (
            #    idx, out[0], output[idx], out[0] - output[idx]
            #)
            yield out[0] - output[idx]


    def __getattribute__(self, name):
        if name.startswith( '___' ):
            return object.__getattribute__(self, name)
        else:
            return super(SurrogateModel, self).__getattribute__(name)


    def __setattribute__(self, name, value):
        if name.startswith( '___' ):
            object.__setattribute__(self, name, value)
        else:
            super(SurrogateModel, self).__setattribute__(name, value)


    @classmethod
    def exceptionLoad(self, surrogateModelId):
        """
        @todo Finding the 'unintialised' models using this method will fail
              eventually fail when running in parallel. Need to pass id of
              calling FireTask. However, this requires additional code in the
              library as well as cooperation of the recipie
        """
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
        i = [0] * len(self.surrogateFunction.inputs)

        # Set inputs
        for k, v in self.inputs.iteritems():
            i[self.inputs_argPos(k)] = inputs[k]

        # Call the surrogate model
        out = cModel.call(in_i, i)

        outputs = {
            k: out[v.argPos]
            for k, v in self.surrogateFunction.outputs.iteritems()
        }

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

            m = re.search('.*\[(.*)\]', kwargs['_id'])
            if m:
                sf = kwargs['surrogateFunction']
                for exp in m.group(1).split(','):
                    m = re.search('(.*)=(.*)', exp)
                    if m:
                        sf.indices[m.group(1)].get_index(m.group(2))
                    else:
                        raise Exception('Unable to parse %s' % exp)

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
                kwargs['initialisationStrategy'] = \
                    EmptyInitialisationStrategy()

            checkAndConvertType(
                kwargs,
                'initialisationStrategy',
                InitialisationStrategy
            )

            DynamicDocument.__init__(self, *args, **kwargs)
            self.save()


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


    def exactTasks(self, points):
        """
        Return an empty workflow
        """
        return Workflow2([])


    def initialisationStrategy(self):
        return loadType(
            self,
            'initialisationStrategy',
            InitialisationStrategy
        )


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

            #print ' '.join(kwargs['inputs'].keys())

            subOutputs = {}
            for m in kwargs['substituteModels']:
                if not isinstance(m, SurrogateModel):
                    raise TypeError(
                        'Elements of substituteModels '
                        'must be derived from SurrogateModel'
                    )
                subOutputs.update(m.outputsToModels())
            #print subOutputs

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
            #print replaced
            #print ' '.join(kwargs['inputs'].keys())

            nInputs = 0
            for k, v in kwargs['inputs'].iteritems():
                kwargs['fitData'][k] = []
                kwargs['inputs'][k] = MinMaxArgPosOpt(**v)

            checkAndConvertType(kwargs, 'exactTask', FireTaskBase)

            checkAndConvertType(
                kwargs,
                'initialisationStrategy',
                InitialisationStrategy
            )

            checkAndConvertType(
                kwargs,
                'outOfBoundsStrategy',
                OutOfBoundsStrategy
            )

            checkAndConvertType(
                kwargs,
                'parameterFittingStrategy',
                ParameterFittingStrategy
            )

            DynamicDocument.__init__(self, *args, **kwargs)

            indices = self.parseIndices()
            for k, v in indices.iteritems():
                kwargs['surrogateFunction'].indices[k].get_index(v)

            self.save()


    def exactTasks(self, points):
        """
        Build a workflow to excute an exactTask for each point
        """

        indices = self.parseIndices()

        # De-serialise the exact task from dict
        et = load_object(self.meth_exactTask)

        tl = []
        e = six.next(six.itervalues(points))
        for i in xrange(len(e)):
            p = { k: points[k][i] for k in points }

            t = et
            t['point'] = p
            t['indices'] = indices
            t['modelId'] = self._id
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
        pfs = loadType(
            self,
            'parameterFittingStrategy',
            ParameterFittingStrategy
        )

        # Nasty hack to work around a bug somewhere in mongoengine or fireworks
        self._changed_fields = filter(
            lambda a: a != u'improveErrorStrategy._fw_name', self._changed_fields
        )

        return pfs


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


    def extendedRange(self, outsidePoint, expansion_factor=1.2):
        """
                  Method expanding the design space. The method ONLY operates
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
 
        @param    outsidePoint     The point that where found to be outside (X)
        @retval   expansion_factor The ratio that is used to expand the space beyond X
    
        @author   Sigve Karolius
        @author   Mandar Thombre
        @todo     Document...
        """

        sampleRange = {}
        limitPoint = {}

        for k, v in self.inputs.iteritems():
            sampleRange[k] = {}
            outsideValue = outsidePoint[k]

            # If the value outside point is outside the range, set the
            # "localdict" max to the outside point value

            if outsideValue > v['max']:
                if outsideValue > self.surrogateFunction.inputs[k].max:
                    raise OutOfBounds('outside point is larger than function min for %s' % k)
                sampleRange[k]['min'] = v['max']
                sampleRange[k]['max'] = min(
                    outsideValue*expansion_factor,
                    self.surrogateFunction.inputs[k].max
                )
                limitPoint[k] = sampleRange[k]['max']

            elif outsideValue < v['min']:
                if outsideValue < self.surrogateFunction.inputs[k].min:
                    raise OutOfBounds('outside point is smaller than function max for %s' % k)
                sampleRange[k]['min'] = max(
                    outsideValue/expansion_factor,
                    self.surrogateFunction.inputs[k].min
                )
                sampleRange[k]['max'] = v['min']
                limitPoint[k] = sampleRange[k]['min']

            else:
                sampleRange[k]['min'] = v['min']
                sampleRange[k]['max'] = v['max']
                limitPoint[k] = random.uniform(v['min'], v['max'])

        return sampleRange, limitPoint

##
# @} # end of python_interface_library

