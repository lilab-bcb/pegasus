from docopt import docopt
from typing import List, Union


class Base:
    """Base class for the commands"""

    def __init__(self, command_args: List[str]):
        self.args = docopt(self.__doc__, argv=command_args)

    def split_string(self, astring: str, sep: str = ",") -> List[str]:
        return astring.split(sep) if astring is not None else []

    def convert_to_int(self, value: str) -> int:
        return int(value) if value is not None else None

    def convert_to_float(self, value: str) -> float:
    	return float(value) if value is not None else None

    def convert_to_list(self, astring: str, sep: str = ",", converter: object = str) -> Union[List[object], None]:
        if astring is None:
            return None
        fields = [converter(x) for x in astring.split(sep)]
        if len(fields) == 1:
            return fields[0]
        return fields

    def execute(self):
        raise NotImplementedError
