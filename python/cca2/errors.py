"""Named refusals.

MATLAB raises ``error('PLQ:conjQ:nonConvexPiece', ...)``; the port raises ``CCA2Error`` carrying
the SAME identifier string, so a MATLAB ``verifyError(..., id)`` ports to
``pytest.raises(CCA2Error)`` plus a check of ``.identifier``.
"""


class CCA2Error(Exception):
    def __init__(self, identifier, message=""):
        super().__init__(f"{identifier}: {message}" if message else identifier)
        self.identifier = identifier
        self.message = message


def error(identifier, message="", *args):
    """Raise a named refusal (MATLAB ``error(id, fmt, args...)``)."""
    if args:
        try:
            message = message % args
        except (TypeError, ValueError):
            message = message + " " + " ".join(str(a) for a in args)
    raise CCA2Error(identifier, message)
