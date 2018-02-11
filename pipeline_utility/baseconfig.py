import os
from copy import deepcopy

import yaml


class BaseConfig(yaml.YAMLObject):

    def __init__(self, filepath, **kwargs):
        self.yaml_tag = unicode("!" + self.__class__.__name__)
        super(BaseConfig, self).__init__()
        self._hidden_fields = ['_hidden_fields', '_filepath', 'yaml_tag']
        self._filepath = filepath

        yaml.add_constructor(self.yaml_tag, self.__class__.from_yaml)
        yaml.add_representer(self.__class__, self.to_yaml)

    def __str__(self):
        return "Config:\n{}".format(self.__dict__)

    def __iter__(self):
        return iter(self.__dict__.items())

    @property
    def filepath(self):
        return self._filepath

    def save(self):
        with open(self._filepath, "wb+") as f:
            # This flow style makes the document human-readable
            yaml.dump(self, f, default_flow_style=False)
            # cfg = yaml.safe_dump(self, f)
        return self

    def load(self):
        # Load variables from yaml into the config class
        # If the file does not exist, create an empty file with default values
        if not os.path.exists(self._filepath):
            self.save()
            print("Created a new empty config file: {0}".format(self._filepath))
            print("Fill out the config and rerun your program!")
            exit(0)
        with open(self._filepath, "r") as f:
            obj = yaml.load(f)
        return obj

    @classmethod
    def to_yaml(cls, dumper, data):
        """
        This overriding function is used to convert the class into a YAML representation.
        The override creates a deep copy of the class and removes any fields marked hidden
        so they won't appear but will continue to exist.
        :param dumper: The default dumper
        :param data: The data represented in the class
        :return:
        """
        new_data = deepcopy(data)
        for item in new_data._hidden_fields:
            if item in new_data.__dict__:
                del new_data.__dict__[item]
        return dumper.represent_yaml_object(data.yaml_tag, new_data, cls,
                                            flow_style=cls.yaml_flow_style)

    # Similar to to_yaml but this will construct the class
    @classmethod
    def from_yaml(cls, loader, node):
        instance = cls.__new__(cls)
        yield instance
        state = loader.construct_mapping(node, deep=True)
        instance.__init__(loader.name, **state)
