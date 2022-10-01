import numpy as np
import sympy as sy
import sys
from phypylib.general_relativity import *
from collections import Counter
import ast


class Field():
    def __init__(self, field, x, y=None, z=None, t=None):
        self.field = field                                                              # field function
        self.t = t                                                                      # time coordinate
        self.x = x                                                                      # first space coordinate
        self.y = y                                                                      # second space coordinate
        self.z = z                                                                      # third space coordinate
        coords = [self.t, self.x, self.y, self.z]
        self.coords = [coord for coord in coords if coord is not None]                  # list of spacetime coordinates
        self.dim = len(self.coords)


class RealScalarField4D(Field):
    """Real scalar field in 4 spacetime dimensions.

    parameter
    ---------
    name : string
        name of the field

    m : sympy symbol
        symbol for the mass term. Can also be set m=0

    metric : spacetime metric from general_relativity.Metric()
        specified metric of spacetime

    # TODO: 3- and 2-dim realscalarfield using this class -> edit MinkowskiMetric()
    """

    def __init__(self, name, m=sy.symbols("m", real=True), metric=MinkowskiMetric()):
        try:
            isinstance(name, str)
        except:
            sys.exit("RealScalarField4D: name variable is a string")

        self.m = m                                                                      # mass term
        self.fieldtype = "4-dim real scalar field"                                      # field type
        self.metric = metric                                                            # metric
        self.name = name                                                                # name as string
        t = self.metric.t
        x = self.metric.x
        y = self.metric.y
        z = self.metric.z
        field = sy.Function(name, real=True)(t, x, y, z)
        super(RealScalarField4D, self).__init__(field=field, x=x, y=y, z=z, t=t)

    def gr_dalembert_operator(self, retG=True, simplify=True, latex=False):
        """Calculates the d'Alembert operator acting on the field Nabla_mu*Nabla^mu*field for a given metric.

        parameter
        ---------
        retG : bool
            if True the result is returned

        simplify : bool
            if True the sympy.simplify function is used

        latex : bool
            if True the result will be printed in latex format

        return
        ------
        dalembert
        """

        covariantpartial = self.metric.covariant_partial(retC=True)
        contravariantpartial = self.metric.contravariant_partial(retC=True)

        # dalembert acting on phi
        dalembert = 0
        for i in range(self.dim):
            for j in range(self.dim):
                dalembert += sy.diff(contravariantpartial[i, j]*sy.diff(self.field, self.coords[j]), self.coords[i])
        if simplify:
            dalembert = sy.simplify(dalembert)

        # correction term from Christoffel symbols
        self.metric.christoffel_symbols(retC=False)
        for i in range(self.dim):
            for j in range(self.dim):
                dalembert += self.metric.cs[i][i, j] * contravariantpartial[i, j] * sy.diff(self.field, self.coords[j])

        if simplify:
            dalembert = sy.simplify(dalembert)
        self.dalembert = dalembert

        if latex:
            sy.print_latex(dalembert)

        if retG:
            return dalembert

    def klein_gordon(self, retK=True, simplify=True, latex=False):
        """Calculate the Klein-Gordon equation for a massive real scalar field in 4D using the given metric by initializing
        the field.

        parameter
        ---------
        retK : bool
            if True the result is returned

        simplify : bool
            if True the sympy.simplify function is used

        latex : bool
            if True the result will be printed in latex format

        return
        ------
        kleingordon
            how to read: 0 = kleingordon
        """
        self.gr_dalembert_operator(retG=False, simplify=simplify, latex=False)

        kleingordon = self.dalembert + self.m**2*self.field
        self.kleingordon = kleingordon

        if latex:
            sy.print_latex(kleingordon)

        if retK:
            return kleingordon


class WickContraction():
    def __init__(self, fields, mode="console", ignore=None):
        # Check for even number of fields
        if len(fields) % 2 != 0:
            sys.exit("Wick contractions yields zero since an off number of fields were given.")
        # Replace fields with integers and collect them in a list. E.g. ["phi_1", "phi_1", "phi_2", "phi_3"] yields
        # [1, 1, 2, 3]
        name_str_list = []                                                              # Create list of names as strings
        for field in fields:
            name_str_list.append(field.name)
        field_indices = []                                                              # Create dict matching every entry
        dict = {name_str_list[0]:1}                                                     # of name_str_list to the
        index = 2                                                                       # corresponding number
        for i in range(1, len(name_str_list)):
            if name_str_list[i] not in dict.keys():
                dict.update({name_str_list[i]:index})
                index += 1
        for name in name_str_list:                                                      # throw field-indices in a list
            field_indices.append(dict[name])
        self.field_indices = field_indices

        self.contractions(field_indices=field_indices, mode=mode, ignore=ignore)

    def contractions(self, field_indices, mode="console", ignore=None):
        index_list = [i for i in range(len(field_indices))]
        res = self.pairgroup(index_list=index_list)
        for i in range(len(res)):
            for j in range(len(res[i])):
                res[i][j] = field_indices[res[i][j]]

        self.count_all_multiples(res=res)

        if ignore is None:
            self.output(mode=mode)

        elif ignore == "vac":
            self.is_graph_connected()

        else:
            sys.exit("WickContraction: ignore does not support the value: " + ignore)

    def pairgroup(self, index_list):
        if len(index_list) == 2:
            return [[index_list[0], index_list[1]]]
        pairedList = []
        startTupleList = self.get_start_tuples(index_list=index_list)
        for i in range(len(startTupleList)):
            res = self.pairgroup(startTupleList[i].friends)
            for j in range(len(res)):
                res2 = res[j]
                comboList = []
                comboList.extend(startTupleList[i].tuple)
                comboList.extend(res2)
                pairedList.append(comboList)
        return pairedList

    def get_start_tuples(self, index_list):
        startTupleList = []
        for tupleIndex in range(1, len(index_list)):
            startTuple = [index_list[0], index_list[tupleIndex]]
            missingNumsList = []
            for missingTupleIndex in range(1, len(index_list)):
                if not missingTupleIndex == tupleIndex:
                    missingNumsList.append(index_list[missingTupleIndex])
            tuple = TupleAndMissingFriends(tuple=startTuple, friends=missingNumsList)
            startTupleList.append(tuple)
        return startTupleList

    def count_all_multiples(self, res):
        pairedList = []
        for pseudo in range(len(res)):
            pairedList.append([])
        for i in range(len(res)):
            for j in range(0, len(res[i]), 2):
                pairedList[i].append(res[i][j:j+2])

        for item in pairedList:
            item.sort(key = lambda x: ((x[0], x[1])))

        stringList = []
        for item in pairedList:
            stringList.append("{}".format(item))
        dic = Counter(stringList)
        uniqueResSet = set(stringList)
        uniqueResList = list(uniqueResSet)
        multiplierList = []
        for item in uniqueResList:
            multiplier = dic[item]
            multiplierList.append(multiplier)

        multiplierList.reverse()
        uniqueResList.reverse()

        self.multiplierList = multiplierList
        self.uniqueResList = uniqueResList

    def output(self, mode="console"):
        if mode == "console":
            print("<0|T{}|0> =".format(self.field_indices))
            print("")
            for i in range(len(self.multiplierList)):
                if i != len(self.multiplierList)-1:
                    print("{} x {} +".format(self.multiplierList[i], self.uniqueResList[i]))
                else:
                    print("{} x {}".format(self.multiplierList[i], self.uniqueResList[i]))

        else:
            print("Else")

    def is_graph_connected(self):
        for graph_str in self.uniqueResList:
            graph = ast.literal_eval(graph_str)                                         # string to list .. thats kinda silly
            for tuple in graph:
                pass # do stuff

class TupleAndMissingFriends():
    def __init__(self, tuple, friends):
        self.tuple = tuple
        self.friends = friends
