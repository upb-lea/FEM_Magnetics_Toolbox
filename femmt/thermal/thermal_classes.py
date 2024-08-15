"""Classes for thermal simulation."""

# Python standard libraries
from typing import List, Tuple, Dict


class ConstraintPro:
    """Boundary constraints."""

    # For boundary constraints, the tuple contains (key, region, value)
    boundary_constraints: List[Tuple[str, str, str]]

    def __init__(self):
        self.boundary_constraints = []

    @staticmethod
    def create_boundary_if_string(flag: int, region: int, value: float):
        """
        Create the boundary region.

        :param flag: flag as an integer
        :type flag: int
        :param region:
        :type region: int
        :param value: value to write in the file
        :type value: float
        """
        return f"\t\tIf({flag}==1)\n\t\t\t{{ Region {region} ; Type Assign; Value {value} ; }}\n\t\tEndIf\n"

    def add_boundary_constraint(self, more_constraints: List):
        """
        Add boundary constraint.

        :param more_constraints: constraints
        :type more_constraints: List
        """
        for constraint in more_constraints:
            self.boundary_constraints.append(constraint)

    def create_file(self, file_path: str):
        """
        Create the constraint.pro file.

        :param file_path: file path
        :type file_path: str
        """
        with open(file_path, "w") as fd:
            fd.write("Constraint {\n  { Name Temperature ;\n\tCase {\n")
            for bc in self.boundary_constraints:
                fd.write(ConstraintPro.create_boundary_if_string(bc[0], bc[1], bc[2]))
            fd.write("\t\t}\n\t}\n}")


class GroupPro:
    """Define the group.pro file for the thermal simulation."""

    regions: dict

    def __init__(self):
        self.regions = {}

    def add_regions(self, more_regions: Dict):
        """
        Add given regions to the group.

        :param more_regions:
        :type more_regions: Dict
        """
        self.regions.update(more_regions)

    def create_file(self, file_path: str, air_gaps_enabled: bool = False, insulation_enabled: bool = False):
        """
        Create some onelab specific result files.

        :param file_path: file path
        :type file_path: str
        :param air_gaps_enabled: "True" to enable air gaps
        :type air_gaps_enabled: bool
        :param insulation_enabled: "True" to enable insulations
        :type insulation_enabled: bool
        """
        with open(file_path, "w") as fd:
            fd.write("Group {\n")
            for key, value in self.regions.items():
                fd.write(f"\t{key} = Region[{value}];\n")
            air_gaps = ""
            if air_gaps_enabled:
                air_gaps = ", air_gaps"
            insulation = ""
            if insulation_enabled:
                insulation = ", insulation"

            fd.write(f"\tCold = Region[{{air, case_top, case_top_right, case_right, case_bot_right, case_bot{air_gaps}{insulation}}}];\n")
            fd.write("\tWarm = Region[{core_total, windings_total}];\n")
            fd.write("\tTotal = Region[{Warm, Cold}];\n")
            fd.write("}")


class ParametersPro:
    """Define the parameter.pro file for the thermal simulation."""

    parameters: dict

    def __init__(self):
        self.parameters = {}

    def add_to_parameters(self, more_parameters: Dict):
        """
        Add more_parameters to existing parameters.

        :param more_parameters:
        :type more_parameters: Dict
        """
        self.parameters.update(more_parameters)

    def create_file(self, file_path: str):
        """
        Create Parameter.pro file.

        :param file_path: Path to the *.pro file
        :type file_path: str
        """
        with open(file_path, "w") as fd:
            for key, value in self.parameters.items():
                if isinstance(value, str):
                    line = f"{key} = \"{value}\";"
                else:
                    line = f"{key} = {value};"

                fd.write(line + "\n")


class FunctionPro:
    """For creating a function.pro."""

    k: dict
    q_vol: dict

    def __init__(self):
        self.k = {}
        self.q_vol = {}

    @staticmethod
    def dict_as_function_str(name: str, dct: Dict):
        """
        Write dictionary as a string.

        :param name: name
        :type name: str
        :param dct: Dictionary
        :type dct: Dict
        """
        dict_as_str = ""
        for key, value in dct.items():
            dict_as_str += f"\t{name}[{key}] = {value};\n"

        return dict_as_str

    def add_dicts(self, k, q_vol):
        """
        Update a key  to the given conductivity and q_vol attributes.

        Order is important: k, qVol.
        :param k: conductivity
        :type k: float
        :param q_vol:
        :type q_vol: float

        """
        if k is not None:
            self.k.update(k)
        if q_vol is not None:
            self.q_vol.update(q_vol)

    def create_file(self, file_path: str):
        """
        Create .pro file.

        :param file_path: File path to generate the .pro file
        :type file_path: str
        """
        with open(file_path, 'w') as fd:
            fd.write("Function {\n")
            fd.write(FunctionPro.dict_as_function_str("k", self.k))
            fd.write(FunctionPro.dict_as_function_str("qVol", self.q_vol))
            fd.write("}")


class PostOperationPro:
    """For creating a post_operation.pro."""

    statements: List[str]    

    def __init__(self):
        self.statements = []

    def add_on_point_statement(self, field, x: float, y: float, formatting, file_name: str, point_name: str, append: bool = False):
        """
        Add x,y point to statement, what is written to the .pro file.

        :param field:
        :type field:
        :param x: x coordinate
        :type x: float
        :param y: y coordinate
        :type y: float
        :param formatting:
        :type formatting:
        :param file_name: file name
        :type file_name: str
        :param point_name: point name
        :type point_name: str
        :param append:
        :type append: bool
        """
        format_str = ""
        if formatting is not None:
            format_str = f"Format {formatting}, "

        append_str = "> " if append else ""

        self.statements.append(
            f"Print [ {field}, OnPoint {{{x}, {y}, 0}}, {format_str}Name \"{point_name}\", File {append_str}\"{file_name}\" ];")

    def add_on_elements_of_statement(self, field, region, file_name: str, formatting=None, depth=None,
                                     name=None, append=False):
        """
        Add elements to a statement for post-operations.

        :param field:
        :type field:
        :param region:
        :type region:
        :param file_name:
        :type file_name: str
        :param formatting:
        :type formatting:
        :param depth:
        :type depth:
        :param name:
        :type name:
        :param append:
        :type append:
        """
        format_str = ""
        depth_str = ""
        name_str = ""

        if formatting is not None:
            format_str = f"Format {formatting}, "

        if depth is not None:
            depth_str = f"Depth {depth}, "

        if name is not None:
            name_str = f"Name \"{name}\", "

        append_str = "> " if append else ""

        self.statements.append(
            f"Print [ {field}, OnElementsOf {region}, {format_str}{depth_str}{name_str}File {append_str}\"{file_name}\"];")

    def create_file(self, file_path: str):
        """
        Create file.

        :param file_path: file path
        :type file_path: str
        """
        with open(file_path, 'w') as fd:
            fd.write("PostOperation map UsingPost The {\n")
            for statement in self.statements:
                fd.write(f"\t{statement}\n")
            fd.write("}")
