"""Attribute descriptors for ReCo."""

import logging
import os


class NameDesc:
    """Descriptor for class name attribute."""
    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        return instance.__dict__[self._name]

    def __set__(self, instance, value):
        try:
            instance.__dict__[self._name] = value
        except ValueError as exc:
            raise ValueError("Attribute 'name' must be a string!") from exc


class LoggerDesc:
    """
    Descriptor for logger attribute.
    """
    def __get__(self, obj, object_type=None):
        return self.value

    def __set__(self, obj, value):
        if not isinstance(value, logging.Logger):
            raise ValueError("'logger' must be a logging.Logger object!")
        self.value = value


class FileDesc:
    """
    Descriptor for files.
    """
    def __set_name__(self, owner, name):
        self.public_name = name
        self.private_name = f"_{name}"

    def __get__(self, obj, object_type=None):
        return getattr(obj, self.private_name)

    def __set__(self, obj, value):
        if value is not None:
            if not os.path.exists(value):
                raise FileNotFoundError(f"File '{value}' does not exist!")
            if not isinstance(value, str):
                raise ValueError(f"Attribute '{value}' must be a string!")
        setattr(obj, self.private_name, value)
