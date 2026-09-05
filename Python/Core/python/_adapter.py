import inspect
import functools
from typing import Optional, Callable, Dict, Any
from pathlib import Path


def _applyLogging(cls, cfg, args, kwargs):
    """Resolve the logging arguments onto `cfg`, and return the remaining args.

    Logging is configured either by a level, which stays an unnamed logger so
    the component names it after itself, or by a ready-made logger the caller
    has named. Everything else is an ordinary config field.
    """
    if not hasattr(cfg, "logger"):
        # not a framework component; its logging stays a constructor argument
        return args

    import acts

    # a level may also arrive positionally, right after the config
    rest = []
    for a in args:
        if isinstance(a, acts.logging.Level) and "level" not in kwargs:
            kwargs["level"] = a
        elif isinstance(a, acts.Logger) and "logger" not in kwargs:
            kwargs["logger"] = a
        else:
            rest.append(a)

    if "level" in kwargs and "logger" in kwargs:
        raise TypeError(f"{cls.__name__}(): pass either level= or logger=, not both")

    # both are ordinary config fields, so set them however the config arrived
    for key in ("level", "logger"):
        if key in kwargs:
            setattr(cfg, key, kwargs.pop(key))
    return rest


def _make_config_adapter(fn):
    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        cls = type(self)

        if len(args) > 0 and isinstance(args[0], inspect.unwrap(cls.Config)):
            cfg = args[0]
            rest = _applyLogging(cls, cfg, args[1:], kwargs)
            fn(self, cfg, *rest, **kwargs)
            return

        if "config" in kwargs:
            cfg = kwargs.pop("config")
            rest = _applyLogging(cls, cfg, args, kwargs)
            fn(self, cfg, *rest, **kwargs)
            return

        cfg = cls.Config()
        args = _applyLogging(cls, cfg, args, kwargs)

        _kwargs = {}
        for k, v in kwargs.items():
            if isinstance(v, Path):
                v = str(v)

            if not hasattr(cfg, k):
                if not hasattr(cfg, "logger"):
                    # not a framework component: its remaining kwargs are
                    # constructor arguments, not config fields
                    _kwargs[k] = v
                    continue
                members = inspect.getmembers(
                    type(cfg), lambda a: not inspect.isroutine(a)
                )
                members = [m for m, _ in members if not m.startswith("_")]
                raise TypeError(
                    "{}(): unexpected keyword argument {!r}. {} has the "
                    "following properties: {}".format(
                        cls.__name__, k, type(cfg).__name__, ", ".join(members)
                    )
                )

            try:
                setattr(cfg, k, v)
            except TypeError as e:
                raise RuntimeError(
                    "{}: Failed to set {}={}".format(type(cfg), k, v)
                ) from e

        fn(self, cfg, *args, **_kwargs)

    return wrapped


def _make_config_constructor(
    cls, proc: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None
):
    fn = cls.__init__

    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        _kwargs = {}
        for k in list(kwargs.keys()):
            if hasattr(cls, k):
                _kwargs[k] = kwargs.pop(k)

        fn(self, *args, **kwargs)

        if proc is not None:
            _kwargs = proc(_kwargs)
        for k, v in _kwargs.items():
            setattr(self, k, v)

    return wrapped


def _patchKwargsConstructor(
    cls, proc: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None
):
    cls.__init__ = _make_config_constructor(cls, proc)


def _patch_config(m):
    for name, cls in inspect.getmembers(m, inspect.isclass):
        if name == "Config":
            _patchKwargsConstructor(cls)

        if hasattr(cls, "Config"):
            cls.__init__ = _make_config_adapter(cls.__init__)
            _patchKwargsConstructor(cls.Config)
