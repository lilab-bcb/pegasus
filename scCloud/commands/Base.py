from docopt import docopt
from typing import List
import logging

logger = logging.getLogger('sccloud')


class Base:
    """Base class for the commands"""

    def __init__(self, command_args: List[str]):
        self.args = docopt(self.__doc__, argv=command_args)
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.StreamHandler())


    def split_string(self, astring: str, sep: str = ",") -> List[str]:
        return astring.split(sep) if astring is not None else []

    def convert_to_int(self, value: str) -> int:
        return int(value) if value is not None else None

    def execute(self):
        raise NotImplementedError
