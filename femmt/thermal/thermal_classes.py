from typing import List, Tuple

class ConstraintPro:
    # For boundary contstraints, the tuple contains (key, region, value)
    boundary_constraints: List[Tuple[str, str, str]]

    def __init__(self):
        self.boundary_constraints = []

    @staticmethod
    def create_boundary_if_string(flag, region, value):
        return f"\t\tIf({flag}==1)\n\t\t\t{{ Region {region} ; Type Assign; Value {value} ; }}\n\t\tEndIf\n"

    def add_boundary_constraint(self, more_constraints):
        for constraint in more_constraints:
            self.boundary_constraints.append(constraint)

    def create_file(self, file_path):
        with open(file_path, "w") as fd:
            fd.write("Constraint {\n  { Name Temperature ;\n\tCase {\n")
            for bc in self.boundary_constraints:
                fd.write(ConstraintPro.create_boundary_if_string(bc[0], bc[1], bc[2]))
            fd.write("\t\t}\n\t}\n}")

class GroupPro:
    regions: dict

    def __init__(self):
        self.regions = {}

    def add_regions(self, more_regions):
        self.regions.update(more_regions)

    def create_file(self, file_path, air_gaps_enabled = False):
        with open(file_path, "w") as fd:
            fd.write("Group {\n")
            for key, value in self.regions.items():
                fd.write(f"\t{key} = Region[{value}];\n")
            if air_gaps_enabled:
                fd.write("\tCold = Region[{air, air_gaps, isolation, case_top, case_top_right, case_right, case_bot_right, case_bot}];\n")
            else:
                fd.write("\tCold = Region[{air, case}];\n")
            fd.write("\tWarm = Region[{core, windings_total}];\n")
            fd.write("\tTotal = Region[{Warm, Cold}];\n")
            fd.write("}")

class ParametersPro:
    """
    For creating a parameters.pro
    """
    parameters: dict

    def __init__(self):
        self.parameters = {}

    def add_to_parameters(self, more_parameters):
        self.parameters.update(more_parameters)

    def create_file(self, file_path):
        """
        file_path: Path to the *.pro file
        parameters: Dict of parameters
        """

        with open(file_path, "w") as fd:
            for key, value in self.parameters.items():
                line = ""''""
                if type(value) is str:
                    line = f"{key} = \"{value}\";"
                else:
                    line = f"{key} = {value};"

                fd.write(line + "\n")

class FunctionPro:
    """
    For creating a function.pro
    """
    k: dict
    q_vol: dict

    def __init__(self):
        self.k = {}
        self.q_vol = {}

    @staticmethod
    def dict_as_fuction_str(name, dct):
        str = ""
        for key, value in dct.items():
            str += f"\t{name}[{key}] = {value};\n"

        return str

    def add_dicts(self, k, q_vol):
        """
        Order is important: k, qVol
        """
        if k is not None:
            self.k.update(k)
        if q_vol is not None:
            self.q_vol.update(q_vol)

    def create_file(self, file_path):
        with open(file_path, 'w') as fd:
            fd.write("Function {\n")
            fd.write(FunctionPro.dict_as_fuction_str("k", self.k))
            fd.write(FunctionPro.dict_as_fuction_str("qVol", self.q_vol))
            fd.write("}")

class PostOperationPro:
    """
    For creating a post_operation.pro
    """
    statements: List[str]    

    def __init__(self):
        self.statements = []

    def add_on_point_statement(self, field, x, y, format, file_name, point_name, append = False):
        format_str = ""
        if format is not None:
            format_str = f"Format {format}, "

        append_str = "> " if append else ""

        self.statements.append(f"Print [ {field}, OnPoint {{{x}, {y}, 0}}, {format_str}Name \"{point_name}\", File {append_str}\"{file_name}\" ];")

    def add_on_elements_of_statement(self, field, region, file_name, format = None, depth = None, name = None, append = False):
        format_str = ""
        depth_str = ""
        name_str = ""

        if format is not None:
            format_str = f"Format {format}, "

        if depth is not None:
            depth_str = f"Depth {depth}, "

        if name is not None:
            name_str = f"Name \"{name}\", "

        append_str = "> " if append else ""

        self.statements.append(f"Print [ {field}, OnElementsOf {region}, {format_str}{depth_str}{name_str}File {append_str}\"{file_name}\"];")

    def create_file(self, file_path):
        with open(file_path, 'w') as fd:
            fd.write("PostOperation map UsingPost The {\n")
            for statement in self.statements:
                fd.write(f"\t{statement}\n")
            fd.write("}")