import os
from copy import deepcopy

import yaml


class BaseConfig(yaml.YAMLObject):
    hidden_fields = ['_filepath']
    yaml_tag = u"!BaseConfig"
    yaml_loader = yaml.SafeLoader
    yaml_dumper = yaml.SafeDumper
    yaml_flow_style = False

    def __init__(self, filepath, **kwargs):
        super(BaseConfig, self).__init__()
        self._filepath = filepath

    def __str__(self):
        return "Config:\n{}".format(self.__dict__)

    @property
    def filepath(self):
        return self._filepath

    def save(self):
        with open(self._filepath, "wb+") as f:
            yaml.safe_dump(self, f)
            # cfg = yaml.safe_dump(self, f)
        return self

    def load(self):
        # Load variables from yaml into the config class
        if not os.path.exists(self._filepath):
            self.save()
            print("Created a new empty config file!")
        with open(self._filepath, "r") as f:
            obj = yaml.load(f)
        return obj

    @classmethod
    def to_yaml(cls, dumper, data):
        new_data = deepcopy(data)
        for item in cls.hidden_fields:
            del new_data.__dict__[item]
        return dumper.represent_yaml_object(cls.yaml_tag, new_data, cls,
                                            flow_style=cls.yaml_flow_style)

    # Static constructor method in class for clarity
    @classmethod
    def cfg_constructor(cls, loader, node):
        # values = loader.construct_mapping(node, deep=True)
        # return cls(values)
        instance = cls.__new__(cls)
        yield instance
        state = loader.construct_mapping(node, deep=True)
        instance.__init__(loader.name, **state)
